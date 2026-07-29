#include <mpi.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <vector>

#include "Snapshot_Heartbeat.hpp"
#include "Snapshot_IO.hpp"
#include "Snapshot_Shared_State.hpp"

using namespace DaMaSCUS_SUN;

namespace
{
bool PathExists(const std::string& path)
{
	struct stat info;
	return lstat(path.c_str(), &info) == 0;
}

void RemoveTree(const std::string& path)
{
	struct stat info;
	if(lstat(path.c_str(), &info) != 0)
		return;
	if(!S_ISDIR(info.st_mode))
	{
		std::remove(path.c_str());
		return;
	}

	DIR* directory = opendir(path.c_str());
	if(directory != nullptr)
	{
		while(dirent* entry = readdir(directory))
		{
			const std::string name(entry->d_name);
			if(name == "." || name == "..")
				continue;
			RemoveTree(path + "/" + name);
		}
		closedir(directory);
	}
	rmdir(path.c_str());
}

std::string ReadAll(const std::string& path)
{
	std::ifstream file(path, std::ios::binary);
	std::ostringstream contents;
	contents << file.rdbuf();
	return contents.str();
}

bool DirectoryContainsTemporaryFile(const std::string& path)
{
	DIR* directory = opendir(path.c_str());
	if(directory == nullptr)
		return false;
	bool found = false;
	while(dirent* entry = readdir(directory))
	{
		if(std::string(entry->d_name).find(".tmp.") != std::string::npos)
		{
			found = true;
			break;
		}
	}
	closedir(directory);
	return found;
}

void Check(bool condition, const std::string& message, int rank, int& failures)
{
	if(condition)
		return;
	std::cerr << "[snapshot-mpi rank " << rank << "] " << message << std::endl;
	failures++;
}

bool BroadcastString(std::string& value, int rank)
{
	int length = (rank == 0) ? static_cast<int>(value.size()) : 0;
	if(MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS || length <= 0)
		return false;
	std::vector<char> bytes(static_cast<size_t>(length));
	if(rank == 0)
		std::memcpy(bytes.data(), value.data(), static_cast<size_t>(length));
	if(MPI_Bcast(bytes.data(), length, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		return false;
	if(rank != 0)
		value.assign(bytes.begin(), bytes.end());
	return true;
}

void SeedSyntheticState(SnapshotSharedState& shared_state, uint64_t run_id, int rank)
{
	shared_state.Initialize(run_id, rank);
	TrajectoryBincount bincount;
	bincount.is_captured = (rank % 2 == 0);
	const double value = static_cast<double>(rank + 1);
	bincount.dt_hist[0] = value;
	bincount.v2dt_hist[0] = 10.0 * value;

	SnapshotEvaporationProgressEntry event;
	event.trajectory_id = static_cast<uint64_t>(1000 + rank);
	event.completion_wall_time_sec = value;
	event.lifetime_unbinding_sec = static_cast<double>(4 - rank);
	shared_state.RecordCompletedTrajectory(
		bincount,
		bincount.is_captured,
		!bincount.is_captured,
		std::vector<SnapshotEvaporationProgressEntry>{event});
}
}

int main(int argc, char* argv[])
{
	int provided = MPI_THREAD_SINGLE;
	if(MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) != MPI_SUCCESS)
		return 1;

	int rank = 0;
	int processes = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &processes);
	int failures = 0;
	Check(processes == 4, "test requires exactly four MPI ranks", rank, failures);
	Check(provided >= MPI_THREAD_FUNNELED, "MPI_THREAD_FUNNELED was not provided", rank, failures);

	int preflight_failure = (failures == 0) ? 0 : 1;
	MPI_Allreduce(MPI_IN_PLACE, &preflight_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if(preflight_failure != 0)
	{
		MPI_Finalize();
		return 1;
	}

	std::string root;
	int directory_ok = 1;
	if(rank == 0)
	{
		char path_template[] = "/tmp/damascus_snapshot_mpi_XXXXXX";
		char* created = mkdtemp(path_template);
		if(created == nullptr)
			directory_ok = 0;
		else
		{
			root = created;
			const std::string snapshot_root = root + "/snapshot";
			const std::string rank_snapshot_dir = snapshot_root + "/rank_snapshot";
			if(mkdir(snapshot_root.c_str(), 0700) != 0 || mkdir(rank_snapshot_dir.c_str(), 0700) != 0)
				directory_ok = 0;
		}
	}
	MPI_Bcast(&directory_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(directory_ok == 0 || !BroadcastString(root, rank))
	{
		if(rank == 0 && !root.empty())
			RemoveTree(root);
		MPI_Finalize();
		return 1;
	}

	const std::string snapshot_root = root + "/snapshot/";
	const std::string rank_snapshot_dir = snapshot_root + "rank_snapshot/";
	const uint64_t run_id = 2026071001ULL;
	const double interval_seconds = 10.0;
	SnapshotSharedState shared_state;
	SeedSyntheticState(shared_state, run_id, rank);

	std::unique_ptr<SnapshotHeartbeat> heartbeat;
	int construction_failed = 0;
	try
	{
		heartbeat.reset(new SnapshotHeartbeat(
			shared_state,
			rank,
			processes,
			run_id,
			snapshot_root,
			rank_snapshot_dir,
			interval_seconds,
			0.5,
			1.0e-40));
	}
	catch(const std::exception& error)
	{
		std::cerr << "[snapshot-mpi rank " << rank << "] heartbeat construction failed: " << error.what() << std::endl;
		construction_failed = 1;
	}
	MPI_Allreduce(MPI_IN_PLACE, &construction_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if(construction_failed != 0)
	{
		if(rank == 0)
			RemoveTree(root);
		MPI_Finalize();
		return 1;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	using Clock = std::chrono::steady_clock;
	const Clock::time_point epoch = Clock::now();
	heartbeat->Start(epoch);

	const Clock::time_point finish_deadline = epoch + std::chrono::milliseconds(10900);
	std::this_thread::sleep_until(finish_deadline);
	const double elapsed_seconds = 1.0e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - epoch).count();
	heartbeat->MarkDoneAndWriteFinal(elapsed_seconds);

	MPI_Barrier(MPI_COMM_WORLD);
	double final_elapsed_seconds = elapsed_seconds;
	MPI_Allreduce(MPI_IN_PLACE, &final_elapsed_seconds, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	bool finalized = true;
	if(rank == 0)
		finalized = heartbeat->FinalizeAfterAllRanksDone(final_elapsed_seconds);
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0)
	{
		const std::string report_path = SnapshotTextFilePath(snapshot_root, 1, interval_seconds);
		const std::string evaporation_path = SnapshotEvaporationTimeFilePath(snapshot_root, 1, interval_seconds);
		const std::string report = ReadAll(report_path);
		const std::string evaporation = ReadAll(evaporation_path);
		Check(finalized, "finalization did not merge every due snapshot", rank, failures);
		Check(PathExists(report_path), "merged snapshot_10s.txt is missing", rank, failures);
		Check(PathExists(evaporation_path), "merged snapshot evaporation file is missing", rank, failures);
		Check(report.find("# snapshot_status = merged") != std::string::npos, "snapshot report is not merged", rank, failures);
		Check(report.find("# snapshot_run_id = 2026071001") != std::string::npos, "snapshot report has the wrong run id", rank, failures);
		Check(report.find("# ready_ranks = 0,1,2,3") != std::string::npos, "snapshot report does not list all four ready ranks", rank, failures);
		Check(report.find("# missing_ranks = \n") != std::string::npos, "snapshot report still lists a missing rank", rank, failures);
		Check(report.find("# total_trajectories = 4") != std::string::npos, "snapshot report has the wrong trajectory total", rank, failures);
		Check(report.find("# captured_particles = 2") != std::string::npos, "snapshot report has the wrong capture total", rank, failures);
		Check(report.find("# valid_trajectories = 4") != std::string::npos, "snapshot report has the wrong classified total", rank, failures);
		Check(report.find("# capture_rate_raw = 0.50000000") != std::string::npos, "snapshot report has the wrong raw capture rate", rank, failures);
		Check(report.find("# capture_rate_valid = 0.50000000") != std::string::npos, "snapshot report has the wrong valid capture rate", rank, failures);
		Check(evaporation.find("# snapshot_status = merged") != std::string::npos, "evaporation report is not merged", rank, failures);

		for(int check_rank = 0; check_rank < 4; check_rank++)
		{
			Check(!PathExists(SnapshotRankCheckpointPath(rank_snapshot_dir, check_rank, 1, interval_seconds)), "rank checkpoint survived final cleanup", rank, failures);
			Check(!PathExists(SnapshotRankFinalPath(rank_snapshot_dir, check_rank)), "rank final state survived final cleanup", rank, failures);
		}
		Check(!PathExists(rank_snapshot_dir), "rank_snapshot directory survived final cleanup", rank, failures);
		Check(!PathExists(SnapshotTextFilePath(snapshot_root, 2, interval_seconds)), "unexpected snapshot_20s.txt was written", rank, failures);
		Check(!DirectoryContainsTemporaryFile(snapshot_root), "temporary snapshot file survived atomic publication", rank, failures);
	}

	int global_failures = 0;
	MPI_Allreduce(&failures, &global_failures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0)
		RemoveTree(root);
	MPI_Finalize();
	return (global_failures == 0) ? 0 : 1;
}

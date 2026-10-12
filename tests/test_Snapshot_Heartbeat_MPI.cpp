#include <mpi.h>

#include <array>
#include <cmath>
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
	const double value = static_cast<double>(rank + 1);
	TrajectoryBincount captured;
	captured.is_captured = true;
	captured.termination_reason = TrajectoryTerminationReason::OutwardEscape;
	captured.dt_hist[0] = value;
	captured.v2dt_hist[0] = 10.0 * value;
	shared_state.RecordCompletedTrajectory(captured, true, false);
	TrajectoryBincount never;
	never.termination_reason = TrajectoryTerminationReason::OutwardEscape;
	never.transit_dt_hist[0] = 2.0 * value;
	never.transit_v2dt_hist[0] = 20.0 * value;
	// A unique shell makes each MPI rank's contribution observable.
	never.transit_dt_hist[rank + 1] = value;
	never.transit_v2dt_hist[rank + 1] = 10.0 * value;
	shared_state.RecordCompletedTrajectory(never, false, true);
}

bool ReadHistogramBin(const std::string& path, int requested_bin, std::array<double, 8>& values)
{
	std::ifstream file(path);
	std::string line;
	while(std::getline(file, line))
	{
		if(line.empty() || line[0] == '#') continue;
		std::istringstream row(line);
		int bin = -1;
		double lower = 0.0, upper = 0.0;
		row >> bin >> lower >> upper;
		for(double& value : values) row >> value;
		if(bin == requested_bin)
		{
			std::string extra;
			return static_cast<bool>(row) && !(row >> extra);
		}
	}
	return false;
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
		const std::string report = ReadAll(report_path);
		Check(finalized, "finalization did not merge every due snapshot", rank, failures);
		Check(PathExists(report_path), "merged snapshot_10s.txt is missing", rank, failures);
		Check(!PathExists(snapshot_root + "snapshot_10s_evaporation_times.txt"), "obsolete snapshot evaporation file was written", rank, failures);
		Check(report.find("# snapshot_status = merged") != std::string::npos, "snapshot report is not merged", rank, failures);
		Check(report.find("# snapshot_run_id = 2026071001") != std::string::npos, "snapshot report has the wrong run id", rank, failures);
		Check(report.find("# ready_ranks = 0,1,2,3") != std::string::npos, "snapshot report does not list all four ready ranks", rank, failures);
		Check(report.find("# missing_ranks = \n") != std::string::npos, "snapshot report still lists a missing rank", rank, failures);
		Check(report.find("# total_trajectories = 8") != std::string::npos, "snapshot report has the wrong trajectory total", rank, failures);
		Check(report.find("# captured_particles = 4") != std::string::npos, "snapshot report has the wrong capture total", rank, failures);
		Check(report.find("# valid_trajectories = 8") != std::string::npos, "snapshot report has the wrong classified total", rank, failures);
		Check(report.find("# capture_rate_raw = 0.50000000") != std::string::npos, "snapshot report has the wrong raw capture rate", rank, failures);
		Check(report.find("# capture_rate_valid = 0.50000000") != std::string::npos, "snapshot report has the wrong valid capture rate", rank, failures);
		Check(report.find("# never_captured_bincount_samples = 4") != std::string::npos, "never-captured denominator omitted an MPI rank", rank, failures);
		Check(report.find("# residence_bincount_samples = 4") != std::string::npos, "captured denominator omitted an MPI rank", rank, failures);
		std::array<double, 8> moments{};
		Check(ReadHistogramBin(report_path, 0, moments), "snapshot histogram must contain exactly 11 columns", rank, failures);
		Check(moments[0] == 10.0 && moments[1] == 100.0, "captured moments were merged incorrectly", rank, failures);
		Check(moments[4] == 20.0 && moments[5] == 200.0, "never-captured moments were merged incorrectly", rank, failures);
		Check(std::abs(moments[6] - std::sqrt(20.0)) < 1e-8 && std::abs(moments[7] - std::sqrt(2000.0)) < 1e-8,
		      "never-captured standard errors use the wrong population", rank, failures);
		for(int check_rank = 0; check_rank < 4; ++check_rank)
		{
			Check(ReadHistogramBin(report_path, check_rank + 1, moments), "rank-specific histogram shell is missing", rank, failures);
			Check(moments[4] == check_rank + 1 && moments[5] == 10.0 * (check_rank + 1),
			      "rank-specific never-captured contribution is missing or duplicated", rank, failures);
		}

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

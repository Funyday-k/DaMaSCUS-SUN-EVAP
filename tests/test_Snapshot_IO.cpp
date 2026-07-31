#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
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

size_t CountOccurrences(const std::string& contents, const std::string& needle)
{
	size_t count = 0;
	size_t position = 0;
	while((position = contents.find(needle, position)) != std::string::npos)
	{
		count++;
		position += needle.size();
	}
	return count;
}

bool ReadHistogramBin(const std::string& path, int requested_bin, std::array<double, 4>& values)
{
	std::ifstream file(path);
	std::string line;
	while(std::getline(file, line))
	{
		if(line.empty() || line[0] == '#')
			continue;
		std::istringstream row(line);
		int bin = -1;
		double radius_lower_rsun = 0.0;
		double radius_upper_rsun = 0.0;
		row >> bin;
		row >> radius_lower_rsun >> radius_upper_rsun;
		for(double& value : values)
			row >> value;
		if(bin == requested_bin)
			return static_cast<bool>(row);
	}
	return false;
}

SnapshotRankState MakeMergeState(uint64_t run_id, int rank)
{
	SnapshotRankState state;
	state.run_id = run_id;
	state.snapshot_index = 1;
	state.rank = rank;
	state.rank_elapsed_wall_sec = 10.0;
	state.local_total = 1;
	state.local_classified = 1;

	const double value = static_cast<double>(rank + 1);
	if(rank % 2 == 0)
	{
		state.local_captured = 1;
		state.bincount_captured_samples = 1;
		state.captured_dt_hist[0] = value;
		state.captured_v2dt_hist[0] = 10.0 * value;
		state.captured_dt_sq_hist[0] = value * value;
		state.captured_v2dt_sq_hist[0] = 100.0 * value * value;
	}
	SnapshotEvaporationProgressEntry event;
	event.trajectory_id = static_cast<uint64_t>(100 + rank);
	event.completion_wall_time_sec = value;
	event.lifetime_unbinding_sec = static_cast<double>(4 - rank);
	state.new_evaporation_events.push_back(event);
	return state;
}

SnapshotRankState MakeRoundTripState(uint64_t run_id)
{
	SnapshotRankState state;
	state.run_id = run_id;
	state.snapshot_index = 3;
	state.rank = 2;
	state.trajectory_in_progress = 1;
	state.local_captured = 1;
	state.local_total = 2;
	state.local_classified = 2;
	state.local_numerical_failures = 1;
	state.bincount_captured_samples = 1;
	state.current_trajectory_id = 901;
	state.rank_elapsed_wall_sec = 30.125;
	state.current_trajectory_wall_sec = 12.5;
	state.current_trajectory_simulated_elapsed_sec = 44.0;
	state.current_trajectory_scatterings = 17;
	state.current_trajectory_captured = 1;
	state.current_trajectory_dt_hist[0] = 0.5;
	state.current_trajectory_v2dt_hist[0] = 5.0;
	state.captured_dt_hist[1] = 2.0;
	state.captured_v2dt_hist[1] = 20.0;
	state.captured_dt_sq_hist[1] = 4.0;
	state.captured_v2dt_sq_hist[1] = 400.0;

	SnapshotEvaporationProgressEntry first;
	first.trajectory_id = 41;
	first.completion_wall_time_sec = 5.0;
	first.lifetime_unbinding_sec = 9.0;
	state.new_evaporation_events.push_back(first);
	SnapshotEvaporationProgressEntry second;
	second.trajectory_id = 42;
	second.completion_wall_time_sec = 7.5;
	second.lifetime_unbinding_sec = 8.0;
	state.new_evaporation_events.push_back(second);
	return state;
}

void ExpectStatesEqual(const SnapshotRankState& expected, const SnapshotRankState& actual)
{
	EXPECT_EQ(expected.run_id, actual.run_id);
	EXPECT_EQ(expected.snapshot_index, actual.snapshot_index);
	EXPECT_EQ(expected.rank, actual.rank);
	EXPECT_EQ(expected.done, actual.done);
	EXPECT_EQ(expected.trajectory_in_progress, actual.trajectory_in_progress);
	EXPECT_EQ(expected.local_captured, actual.local_captured);
	EXPECT_EQ(expected.local_total, actual.local_total);
	EXPECT_EQ(expected.local_classified, actual.local_classified);
	EXPECT_EQ(expected.local_numerical_failures, actual.local_numerical_failures);
	EXPECT_EQ(expected.bincount_captured_samples, actual.bincount_captured_samples);
	EXPECT_EQ(expected.current_trajectory_id, actual.current_trajectory_id);
	EXPECT_DOUBLE_EQ(expected.rank_elapsed_wall_sec, actual.rank_elapsed_wall_sec);
	EXPECT_DOUBLE_EQ(expected.current_trajectory_wall_sec, actual.current_trajectory_wall_sec);
	EXPECT_DOUBLE_EQ(expected.current_trajectory_simulated_elapsed_sec, actual.current_trajectory_simulated_elapsed_sec);
	EXPECT_EQ(expected.current_trajectory_scatterings, actual.current_trajectory_scatterings);
	EXPECT_EQ(expected.current_trajectory_captured, actual.current_trajectory_captured);
	EXPECT_EQ(expected.current_trajectory_dt_hist, actual.current_trajectory_dt_hist);
	EXPECT_EQ(expected.current_trajectory_v2dt_hist, actual.current_trajectory_v2dt_hist);
	EXPECT_EQ(expected.captured_dt_hist, actual.captured_dt_hist);
	EXPECT_EQ(expected.captured_v2dt_hist, actual.captured_v2dt_hist);
	EXPECT_EQ(expected.captured_dt_sq_hist, actual.captured_dt_sq_hist);
	EXPECT_EQ(expected.captured_v2dt_sq_hist, actual.captured_v2dt_sq_hist);
	ASSERT_EQ(expected.new_evaporation_events.size(), actual.new_evaporation_events.size());
	for(size_t index = 0; index < expected.new_evaporation_events.size(); index++)
	{
		EXPECT_EQ(expected.new_evaporation_events[index].trajectory_id, actual.new_evaporation_events[index].trajectory_id);
		EXPECT_DOUBLE_EQ(expected.new_evaporation_events[index].completion_wall_time_sec, actual.new_evaporation_events[index].completion_wall_time_sec);
		EXPECT_DOUBLE_EQ(expected.new_evaporation_events[index].lifetime_unbinding_sec, actual.new_evaporation_events[index].lifetime_unbinding_sec);
		EXPECT_DOUBLE_EQ(expected.new_evaporation_events[index].r_capture_rsun, actual.new_evaporation_events[index].r_capture_rsun);
		EXPECT_DOUBLE_EQ(expected.new_evaporation_events[index].E_capture_eV, actual.new_evaporation_events[index].E_capture_eV);
		EXPECT_DOUBLE_EQ(expected.new_evaporation_events[index].dE_capture_eV, actual.new_evaporation_events[index].dE_capture_eV);
	}
}

class SnapshotIOTest : public ::testing::Test
{
  protected:
	void SetUp() override
	{
		char path_template[] = "/tmp/damascus_snapshot_io_XXXXXX";
		char* created = mkdtemp(path_template);
		ASSERT_NE(nullptr, created);
		root = created;
		snapshot_root = root + "/snapshot/";
		rank_snapshot_dir = snapshot_root + "rank_snapshot/";
		ASSERT_EQ(0, mkdir(snapshot_root.c_str(), 0700));
		ASSERT_EQ(0, mkdir(rank_snapshot_dir.c_str(), 0700));
	}

	void TearDown() override
	{
		RemoveTree(root);
	}

	std::string root;
	std::string snapshot_root;
	std::string rank_snapshot_dir;
};
}

TEST_F(SnapshotIOTest, AcceptsOnlyPositiveIntegerSnapshotIntervals)
{
	EXPECT_TRUE(IsValidSnapshotIntervalSeconds(1.0));
	EXPECT_TRUE(IsValidSnapshotIntervalSeconds(10.0));
	EXPECT_TRUE(IsValidSnapshotIntervalSeconds(60.0));
	EXPECT_FALSE(IsValidSnapshotIntervalSeconds(0.0));
	EXPECT_FALSE(IsValidSnapshotIntervalSeconds(-1.0));
	EXPECT_FALSE(IsValidSnapshotIntervalSeconds(0.5));
	EXPECT_FALSE(IsValidSnapshotIntervalSeconds(10.25));
	EXPECT_FALSE(IsValidSnapshotIntervalSeconds(std::numeric_limits<double>::quiet_NaN()));
	EXPECT_FALSE(IsValidSnapshotIntervalSeconds(std::numeric_limits<double>::infinity()));

	EXPECT_EQ(20, SnapshotTimeLabelSeconds(2, 10.0));
	EXPECT_EQ(snapshot_root + "snapshot_10s.txt", SnapshotTextFilePath(snapshot_root, 1, 10.0));
	EXPECT_EQ(snapshot_root + "snapshot_10s_evaporation_times.txt", SnapshotEvaporationTimeFilePath(snapshot_root, 1, 10.0));
	EXPECT_EQ(rank_snapshot_dir + "snapshot_10s_rank3.bin", SnapshotRankCheckpointPath(rank_snapshot_dir, 3, 1, 10.0));
	EXPECT_THROW(SnapshotTimeLabelSeconds(1, 10.25), std::invalid_argument);

	SnapshotSharedState shared_state;
	shared_state.Initialize(77, 0);
	EXPECT_THROW(
		SnapshotHeartbeat(shared_state, 0, 1, 77, snapshot_root, rank_snapshot_dir, 10.25, 1.0, 1.0e-40),
		std::invalid_argument);
}

TEST_F(SnapshotIOTest, BinaryRankStateRoundTripsAndRejectsCorruption)
{
	const uint64_t run_id = 0x123456789ULL;
	const SnapshotRankState expected = MakeRoundTripState(run_id);
	const std::string path = rank_snapshot_dir + "roundtrip.bin";
	ASSERT_TRUE(WriteSnapshotRankState(path, expected));

	SnapshotRankState actual;
	ASSERT_TRUE(ReadSnapshotRankState(path, run_id, actual));
	ExpectStatesEqual(expected, actual);
	EXPECT_FALSE(ReadSnapshotRankState(path, run_id + 1, actual));

	ASSERT_TRUE(WriteSnapshotRankState(path, expected));
	{
		std::fstream file(path, std::ios::binary | std::ios::in | std::ios::out);
		ASSERT_TRUE(file.is_open());
		const char corrupt_magic = 0;
		file.write(&corrupt_magic, 1);
	}
	EXPECT_FALSE(ReadSnapshotRankState(path, run_id, actual));

	ASSERT_TRUE(WriteSnapshotRankState(path, expected));
	{
		std::ofstream file(path, std::ios::binary | std::ios::app);
		ASSERT_TRUE(file.is_open());
		const char trailing_byte = 1;
		file.write(&trailing_byte, 1);
	}
	EXPECT_FALSE(ReadSnapshotRankState(path, run_id, actual));

	ASSERT_TRUE(WriteSnapshotRankState(path, expected));
	struct stat info;
	ASSERT_EQ(0, stat(path.c_str(), &info));
	ASSERT_GT(info.st_size, 1);
	ASSERT_EQ(0, truncate(path.c_str(), info.st_size - 1));
	EXPECT_FALSE(ReadSnapshotRankState(path, run_id, actual));
}

TEST_F(SnapshotIOTest, SharedStatePublishesCurrentTrajectoryProgress)
{
	SnapshotSharedState shared_state;
	shared_state.Initialize(501, 4);
	shared_state.BeginTrajectory(77, 100.0);
	std::array<double, TOTAL_BINS> dt_hist{};
	std::array<double, TOTAL_BINS> v2dt_hist{};
	dt_hist[3] = 1.0;
	v2dt_hist[3] = 4.0;
	dt_hist[4] = 1.5;
	v2dt_hist[4] = 6.0;
	shared_state.PublishCurrentTrajectoryProgress(dt_hist, v2dt_hist, 102.5);
	shared_state.UpdateCurrentScatterings(9);
	shared_state.MarkCurrentCaptured(true);

	size_t evaporation_end = 0;
	const SnapshotRankState running = shared_state.CopyForSnapshot(1, 10.0, 10.0, 0, evaporation_end);
	EXPECT_EQ(1, running.trajectory_in_progress);
	EXPECT_EQ(77U, running.current_trajectory_id);
	EXPECT_EQ(1, running.current_trajectory_captured);
	EXPECT_GE(running.current_trajectory_wall_sec, 0.0);
	EXPECT_DOUBLE_EQ(2.5, running.current_trajectory_simulated_elapsed_sec);
	EXPECT_EQ(9U, running.current_trajectory_scatterings);
	EXPECT_DOUBLE_EQ(1.0, running.current_trajectory_dt_hist[3]);
	EXPECT_DOUBLE_EQ(4.0, running.current_trajectory_v2dt_hist[3]);
	EXPECT_DOUBLE_EQ(1.5, running.current_trajectory_dt_hist[4]);
	EXPECT_DOUBLE_EQ(6.0, running.current_trajectory_v2dt_hist[4]);

	TrajectoryBincount completed;
	completed.is_captured = true;
	completed.termination_reason = TrajectoryTerminationReason::OutwardEscape;
	completed.dt_hist = dt_hist;
	completed.v2dt_hist = v2dt_hist;
	shared_state.RecordCompletedTrajectory(completed, true, false, {});
	TrajectoryBincount failed;
	failed.termination_reason = TrajectoryTerminationReason::NumericalFailure;
	shared_state.BeginTrajectory(78, 103.0);
	shared_state.RecordCompletedTrajectory(failed, false, false, {});
	const SnapshotRankState idle = shared_state.CopyForSnapshot(2, 20.0, 20.0, 0, evaporation_end);
	EXPECT_EQ(0, idle.trajectory_in_progress);
	EXPECT_EQ(0U, idle.current_trajectory_id);
	EXPECT_DOUBLE_EQ(0.0, idle.current_trajectory_wall_sec);
	EXPECT_DOUBLE_EQ(0.0, idle.current_trajectory_simulated_elapsed_sec);
	EXPECT_EQ(0U, idle.current_trajectory_scatterings);
	EXPECT_EQ(1U, idle.local_numerical_failures);
	EXPECT_EQ(1U, idle.bincount_captured_samples);
	EXPECT_DOUBLE_EQ(1.0, idle.captured_dt_hist[3]);
	EXPECT_DOUBLE_EQ(4.0, idle.captured_v2dt_hist[3]);
	EXPECT_DOUBLE_EQ(0.0, idle.current_trajectory_dt_hist[3]);
	EXPECT_DOUBLE_EQ(0.0, idle.current_trajectory_v2dt_hist[3]);
}

TEST_F(SnapshotIOTest, TextReportListsEachMpiRankActivity)
{
	const uint64_t run_id = 606;
	const double interval = 10.0;

	SnapshotRankState running_uncaptured;
	running_uncaptured.run_id = run_id;
	running_uncaptured.snapshot_index = 1;
	running_uncaptured.rank = 0;
	running_uncaptured.trajectory_in_progress = 1;
	running_uncaptured.current_trajectory_id = 1001;
	running_uncaptured.rank_elapsed_wall_sec = 10.0;
	running_uncaptured.current_trajectory_wall_sec = 3.5;
	running_uncaptured.current_trajectory_simulated_elapsed_sec = 125.0;
	running_uncaptured.current_trajectory_scatterings = 7;

	SnapshotRankState running_captured = running_uncaptured;
	running_captured.rank = 1;
	running_captured.current_trajectory_id = 2001;
	running_captured.current_trajectory_wall_sec = 8.0;
	running_captured.current_trajectory_simulated_elapsed_sec = 2500.0;
	running_captured.current_trajectory_scatterings = 19;
	running_captured.current_trajectory_captured = 1;
	running_captured.current_trajectory_dt_hist[TOTAL_BINS - 1] = 42.0;
	running_captured.current_trajectory_v2dt_hist[TOTAL_BINS - 1] = 84.0;

	SnapshotRankState idle;
	idle.run_id = run_id;
	idle.snapshot_index = 1;
	idle.rank = 2;
	idle.rank_elapsed_wall_sec = 10.0;

	SnapshotRankState done;
	done.run_id = run_id;
	done.snapshot_index = 0;
	done.rank = 3;
	done.done = 1;
	done.rank_elapsed_wall_sec = 5.0;

	ASSERT_TRUE(WriteSnapshotRankState(
		SnapshotRankCheckpointPath(rank_snapshot_dir, 0, 1, interval), running_uncaptured));
	ASSERT_TRUE(WriteSnapshotRankState(
		SnapshotRankCheckpointPath(rank_snapshot_dir, 1, 1, interval), running_captured));
	ASSERT_TRUE(WriteSnapshotRankState(
		SnapshotRankCheckpointPath(rank_snapshot_dir, 2, 1, interval), idle));
	ASSERT_TRUE(WriteSnapshotRankState(SnapshotRankFinalPath(rank_snapshot_dir, 3), done));

	const SnapshotMergeResult merged = TryWriteSnapshot(
		snapshot_root, rank_snapshot_dir, 1, interval, 4, run_id, 0.5, 1.0e-40, false);
	ASSERT_EQ(SnapshotMergeStatus::Merged, merged.status);
	const std::string report = ReadAll(SnapshotTextFilePath(snapshot_root, 1, interval));
	// The completed counters remain unchanged, while the already captured
	// in-progress trajectory contributes its provisional residence prefix.
	EXPECT_NE(std::string::npos, report.find("# total_trajectories = 0"));
	EXPECT_NE(std::string::npos, report.find("# captured_particles = 0"));
	EXPECT_NE(std::string::npos, report.find("# valid_trajectories = 0"));
	EXPECT_NE(std::string::npos, report.find("# in_progress_bincount_included = 1"));
	EXPECT_NE(std::string::npos, report.find("# in_progress_bincount_captured_samples = 1"));
	EXPECT_NE(std::string::npos, report.find("# in_progress_bincount_is_provisional = 1"));
	EXPECT_NE(std::string::npos, report.find("# residence_bincount_samples = 1"));
	EXPECT_NE(std::string::npos, report.find("\t4.2000000000e+01\t8.4000000000e+01\t"));
	EXPECT_NE(std::string::npos, report.find("# [MPI rank status]"));
	EXPECT_NE(std::string::npos, report.find(
		"# rank  state  trajectory_id  trajectory_wall_s  simulated_elapsed_s  scatterings  observed_at_wall_s"));
	EXPECT_NE(std::string::npos, report.find(
		"# 0\trunning_uncaptured\t1001\t3.5000000000e+00\t1.2500000000e+02\t7\t1.0000000000e+01"));
	EXPECT_NE(std::string::npos, report.find(
		"# 1\trunning_captured\t2001\t8.0000000000e+00\t2.5000000000e+03\t19\t1.0000000000e+01"));
	EXPECT_NE(std::string::npos, report.find(
		"# 2\tidle_or_waiting\t0\t0.0000000000e+00\t0.0000000000e+00\t0\t1.0000000000e+01"));
	EXPECT_NE(std::string::npos, report.find(
		"# 3\tdone\t0\t0.0000000000e+00\t0.0000000000e+00\t0\t5.0000000000e+00"));
}

TEST_F(SnapshotIOTest, CompletedOuterDomainRemovalContributesThroughFivePointTwoAU)
{
	const uint64_t run_id = 607;
	const double interval = 10.0;
	SnapshotSharedState shared_state;
	shared_state.Initialize(run_id, 0);
	shared_state.BeginTrajectory(9001, 0.0);

	TrajectoryBincount removed;
	removed.is_captured = true;
	removed.outer_domain_removed = true;
	removed.termination_reason =
	    TrajectoryTerminationReason::OuterDomainRemoval;
	removed.dt_hist[TOTAL_BINS - 1] = 12.5;
	removed.v2dt_hist[TOTAL_BINS - 1] = 25.0;
	shared_state.RecordCompletedTrajectory(
	    removed, true, false, {});

	size_t evaporation_end = 0;
	const SnapshotRankState checkpoint =
	    shared_state.CopyForSnapshot(
	        1, interval, interval, 0, evaporation_end);
	ASSERT_TRUE(WriteSnapshotRankState(
	    SnapshotRankCheckpointPath(
	        rank_snapshot_dir, 0, 1, interval),
	    checkpoint));

	const SnapshotMergeResult merged = TryWriteSnapshot(
	    snapshot_root,
	    rank_snapshot_dir,
	    1,
	    interval,
	    1,
	    run_id,
	    0.5,
	    1.0e-40,
	    false);
	ASSERT_EQ(SnapshotMergeStatus::Merged, merged.status);
	const std::string report =
	    ReadAll(SnapshotTextFilePath(snapshot_root, 1, interval));
	EXPECT_NE(
	    std::string::npos,
	    report.find("# residence_bincount_samples = 1"));
	EXPECT_NE(
	    std::string::npos,
	    report.find("# in_progress_bincount_included = 0"));
	EXPECT_NE(
	    std::string::npos,
	    report.find("\t1.2500000000e+01\t2.5000000000e+01\t"));
}

TEST_F(SnapshotIOTest, FourLogicalRanksProgressFromPartialToMergedWithoutDowngrade)
{
	const uint64_t run_id = 20260710;
	const double interval = 10.0;
	for(int rank = 0; rank < 3; rank++)
	{
		ASSERT_TRUE(WriteSnapshotRankState(
			SnapshotRankCheckpointPath(rank_snapshot_dir, rank, 1, interval),
			MakeMergeState(run_id, rank)));
	}

	const SnapshotMergeResult partial = TryWriteSnapshot(
		snapshot_root, rank_snapshot_dir, 1, interval, 4, run_id, 0.5, 1.0e-40, true);
	EXPECT_EQ(SnapshotMergeStatus::Partial, partial.status);
	EXPECT_TRUE(partial.cleanup_succeeded);
	EXPECT_EQ((std::vector<int>{0, 1, 2}), partial.ready_ranks);
	EXPECT_EQ((std::vector<int>{3}), partial.missing_ranks);

	const std::string report_path = SnapshotTextFilePath(snapshot_root, 1, interval);
	const std::string evaporation_path = SnapshotEvaporationTimeFilePath(snapshot_root, 1, interval);
	const std::string partial_report = ReadAll(report_path);
	EXPECT_NE(std::string::npos, partial_report.find("# snapshot_status = partial"));
	EXPECT_NE(std::string::npos, partial_report.find("# ready_ranks = 0,1,2"));
	EXPECT_NE(std::string::npos, partial_report.find("# missing_ranks = 3"));
	EXPECT_NE(std::string::npos, partial_report.find("# total_trajectories = 3"));
	EXPECT_NE(std::string::npos, partial_report.find("# captured_particles = 2"));
	EXPECT_NE(std::string::npos, partial_report.find("# capture_rate_raw = 0.66666667"));
	for(int rank = 0; rank < 3; rank++)
		EXPECT_TRUE(PathExists(SnapshotRankCheckpointPath(rank_snapshot_dir, rank, 1, interval)));

	ASSERT_TRUE(WriteSnapshotRankState(
		SnapshotRankCheckpointPath(rank_snapshot_dir, 3, 1, interval),
		MakeMergeState(run_id, 3)));
	const SnapshotMergeResult merged = TryWriteSnapshot(
		snapshot_root, rank_snapshot_dir, 1, interval, 4, run_id, 0.5, 1.0e-40, false);
	EXPECT_EQ(SnapshotMergeStatus::Merged, merged.status);
	EXPECT_TRUE(merged.cleanup_succeeded);
	EXPECT_EQ((std::vector<int>{0, 1, 2, 3}), merged.ready_ranks);
	EXPECT_TRUE(merged.missing_ranks.empty());

	const std::string merged_report = ReadAll(report_path);
	const std::string merged_evaporation = ReadAll(evaporation_path);
	EXPECT_NE(std::string::npos, merged_report.find("# snapshot_status = merged"));
	EXPECT_NE(std::string::npos, merged_report.find("# ready_ranks = 0,1,2,3"));
	EXPECT_NE(std::string::npos, merged_report.find("# missing_ranks = \n"));
	EXPECT_NE(std::string::npos, merged_report.find("# total_trajectories = 4"));
	EXPECT_NE(std::string::npos, merged_report.find("# captured_particles = 2"));
	EXPECT_NE(std::string::npos, merged_report.find("# capture_rate_raw = 0.50000000"));

	std::array<double, 4> histogram{};
	ASSERT_TRUE(ReadHistogramBin(report_path, 0, histogram));
	EXPECT_DOUBLE_EQ(4.0, histogram[0]);
	EXPECT_DOUBLE_EQ(40.0, histogram[1]);
	EXPECT_NEAR(std::sqrt(2.0), histogram[2], 1.0e-9);
	EXPECT_NEAR(std::sqrt(200.0), histogram[3], 1.0e-8);

	const size_t rank3_event = merged_evaporation.find("3\t103\t1.0000000000e+00");
	const size_t rank2_event = merged_evaporation.find("2\t102\t2.0000000000e+00");
	const size_t rank1_event = merged_evaporation.find("1\t101\t3.0000000000e+00");
	const size_t rank0_event = merged_evaporation.find("0\t100\t4.0000000000e+00");
	ASSERT_NE(std::string::npos, rank3_event);
	ASSERT_NE(std::string::npos, rank2_event);
	ASSERT_NE(std::string::npos, rank1_event);
	ASSERT_NE(std::string::npos, rank0_event);
	EXPECT_LT(rank3_event, rank2_event);
	EXPECT_LT(rank2_event, rank1_event);
	EXPECT_LT(rank1_event, rank0_event);

	for(int rank = 0; rank < 4; rank++)
		EXPECT_FALSE(PathExists(SnapshotRankCheckpointPath(rank_snapshot_dir, rank, 1, interval)));

	ASSERT_TRUE(WriteSnapshotRankState(
		SnapshotRankCheckpointPath(rank_snapshot_dir, 0, 1, interval),
		MakeMergeState(run_id, 0)));
	const SnapshotMergeResult repeat = TryWriteSnapshot(
		snapshot_root, rank_snapshot_dir, 1, interval, 4, run_id, 0.5, 1.0e-40, true);
	EXPECT_EQ(SnapshotMergeStatus::Merged, repeat.status);
	EXPECT_TRUE(repeat.cleanup_succeeded);
	EXPECT_EQ(merged_report, ReadAll(report_path));
	EXPECT_EQ(merged_evaporation, ReadAll(evaporation_path));
	EXPECT_FALSE(PathExists(SnapshotRankCheckpointPath(rank_snapshot_dir, 0, 1, interval)));
}

TEST_F(SnapshotIOTest, LateBoundaryEventIsPublishedOnceByNextCheckpoint)
{
	const uint64_t run_id = 88001;
	const double interval = 10.0;
	SnapshotSharedState shared_state;
	shared_state.Initialize(run_id, 0);

	size_t first_event_end = 0;
	const SnapshotRankState first_checkpoint = shared_state.CopyForSnapshot(
		1, 10.0, 10.0, 0, first_event_end);
	EXPECT_TRUE(first_checkpoint.new_evaporation_events.empty());
	EXPECT_EQ(0U, first_event_end);

	TrajectoryBincount bincount;
	bincount.is_captured = true;
	bincount.dt_hist[0] = 1.0;
	bincount.v2dt_hist[0] = 2.0;
	SnapshotEvaporationProgressEntry late_event;
	late_event.trajectory_id = 777;
	late_event.completion_wall_time_sec = 9.9;
	late_event.lifetime_unbinding_sec = 3.0;
	shared_state.RecordCompletedTrajectory(
		bincount, true, false, std::vector<SnapshotEvaporationProgressEntry>{late_event});

	size_t second_event_end = first_event_end;
	const SnapshotRankState second_checkpoint = shared_state.CopyForSnapshot(
		2, 20.0, 20.0, first_event_end, second_event_end);
	ASSERT_EQ(1U, second_checkpoint.new_evaporation_events.size());
	EXPECT_EQ(1U, second_event_end);
	ASSERT_TRUE(WriteSnapshotRankState(
		SnapshotRankCheckpointPath(rank_snapshot_dir, 0, 2, interval), second_checkpoint));

	const SnapshotMergeResult merged = TryWriteSnapshot(
		snapshot_root, rank_snapshot_dir, 2, interval, 1, run_id, 0.5, 1.0e-40, false);
	EXPECT_EQ(SnapshotMergeStatus::Merged, merged.status);
	const std::string evaporation = ReadAll(
		SnapshotEvaporationTimeFilePath(snapshot_root, 2, interval));
	EXPECT_EQ(1U, CountOccurrences(evaporation, "0\t777\t"));
	EXPECT_NE(std::string::npos, evaporation.find("# completed_evaporation_events_assigned_to_snapshot = 1"));
}

TEST_F(SnapshotIOTest, ReportsRawAndClassifiedCaptureRatesSeparately)
{
	const uint64_t run_id = 88003;
	const double interval = 10.0;
	SnapshotRankState state;
	state.run_id = run_id;
	state.snapshot_index = 1;
	state.rank = 0;
	state.rank_elapsed_wall_sec = 10.0;
	state.local_total = 2;
	state.local_classified = 1;
	state.local_numerical_failures = 1;
	state.local_captured = 1;
	state.bincount_captured_samples = 1;
	state.captured_dt_hist[0] = 1.0;
	state.captured_v2dt_hist[0] = 2.0;
	state.captured_dt_sq_hist[0] = 1.0;
	state.captured_v2dt_sq_hist[0] = 4.0;
	ASSERT_TRUE(WriteSnapshotRankState(
		SnapshotRankCheckpointPath(rank_snapshot_dir, 0, 1, interval), state));

	const SnapshotMergeResult merged = TryWriteSnapshot(
		snapshot_root, rank_snapshot_dir, 1, interval, 1, run_id, 0.5, 1.0e-40, false);
	ASSERT_EQ(SnapshotMergeStatus::Merged, merged.status);
	const std::string report = ReadAll(SnapshotTextFilePath(snapshot_root, 1, interval));
	EXPECT_NE(std::string::npos, report.find("# total_trajectories = 2"));
	EXPECT_NE(std::string::npos, report.find("# valid_trajectories = 1"));
	EXPECT_NE(std::string::npos, report.find("# numerical_failures = 1"));
	EXPECT_NE(std::string::npos, report.find("# unresolved_not_captured_trajectories = 1"));
	EXPECT_NE(std::string::npos, report.find("# capture_rate_raw = 0.50000000"));
	EXPECT_NE(std::string::npos, report.find("# capture_rate_valid = 1.00000000"));
}

TEST_F(SnapshotIOTest, FinalStateAssignsEventsToFirstEligibleSnapshotWithoutRepeating)
{
	const uint64_t run_id = 88002;
	const double interval = 10.0;
	SnapshotRankState final_state;
	final_state.run_id = run_id;
	final_state.rank = 0;
	final_state.done = 1;
	final_state.snapshot_index = 0;
	final_state.rank_elapsed_wall_sec = 15.0;
	final_state.local_total = 2;
	final_state.local_classified = 2;

	SnapshotEvaporationProgressEntry first_event;
	first_event.trajectory_id = 901;
	first_event.completion_wall_time_sec = 5.0;
	first_event.lifetime_unbinding_sec = 2.0;
	SnapshotEvaporationProgressEntry second_event;
	second_event.trajectory_id = 902;
	second_event.completion_wall_time_sec = 15.0;
	second_event.lifetime_unbinding_sec = 1.0;
	final_state.new_evaporation_events = {first_event, second_event};
	ASSERT_TRUE(WriteSnapshotRankState(
		SnapshotRankFinalPath(rank_snapshot_dir, 0), final_state));

	const SnapshotMergeResult too_early = TryWriteSnapshot(
		snapshot_root, rank_snapshot_dir, 1, interval, 1, run_id, 0.5, 1.0e-40, false);
	EXPECT_EQ(SnapshotMergeStatus::NoRanksReady, too_early.status);
	EXPECT_FALSE(PathExists(SnapshotEvaporationTimeFilePath(snapshot_root, 1, interval)));

	const SnapshotMergeResult first_eligible = TryWriteSnapshot(
		snapshot_root, rank_snapshot_dir, 2, interval, 1, run_id, 0.5, 1.0e-40, false);
	EXPECT_EQ(SnapshotMergeStatus::Merged, first_eligible.status);
	const std::string second_report = ReadAll(
		SnapshotEvaporationTimeFilePath(snapshot_root, 2, interval));
	EXPECT_EQ(1U, CountOccurrences(second_report, "0\t901\t"));
	EXPECT_EQ(1U, CountOccurrences(second_report, "0\t902\t"));

	const SnapshotMergeResult later = TryWriteSnapshot(
		snapshot_root, rank_snapshot_dir, 3, interval, 1, run_id, 0.5, 1.0e-40, false);
	EXPECT_EQ(SnapshotMergeStatus::Merged, later.status);
	const std::string third_report = ReadAll(
		SnapshotEvaporationTimeFilePath(snapshot_root, 3, interval));
	EXPECT_EQ(0U, CountOccurrences(third_report, "0\t901\t"));
	EXPECT_EQ(0U, CountOccurrences(third_report, "0\t902\t"));

	SnapshotRankState invalid_final = final_state;
	invalid_final.rank_elapsed_wall_sec = 14.0;
	EXPECT_FALSE(WriteSnapshotRankState(rank_snapshot_dir + "invalid_final.bin", invalid_final));
}

TEST_F(SnapshotIOTest, IncompleteFinalizationPreservesRecoverableFinalStates)
{
	const uint64_t run_id = 99001;
	SnapshotRankState final_state = MakeMergeState(run_id, 0);
	final_state.snapshot_index = 0;
	final_state.done = 1;
	final_state.rank_elapsed_wall_sec = 120.0;
	const std::string final_path = SnapshotRankFinalPath(rank_snapshot_dir, 0);
	ASSERT_TRUE(WriteSnapshotRankState(final_path, final_state));

	SnapshotSharedState shared_state;
	shared_state.Initialize(run_id, 0);
	SnapshotHeartbeat heartbeat(
		shared_state, 0, 1, run_id, snapshot_root, rank_snapshot_dir, 60.0, 1.0, 1.0e-40);

	EXPECT_FALSE(heartbeat.FinalizeAfterAllRanksDone(60.0));
	EXPECT_TRUE(PathExists(final_path));
	EXPECT_TRUE(PathExists(rank_snapshot_dir));
}

TEST_F(SnapshotIOTest, FinalizationReportsCleanupFailure)
{
	const uint64_t run_id = 99002;
	const std::string blocked_final_path = SnapshotRankFinalPath(rank_snapshot_dir, 0);
	ASSERT_EQ(0, mkdir(blocked_final_path.c_str(), 0700));
	{
		std::ofstream blocker(blocked_final_path + "/keep");
		ASSERT_TRUE(blocker.is_open());
		blocker << "keep";
	}

	SnapshotSharedState shared_state;
	shared_state.Initialize(run_id, 0);
	SnapshotHeartbeat heartbeat(
		shared_state, 0, 1, run_id, snapshot_root, rank_snapshot_dir, 10.0, 1.0, 1.0e-40);

	EXPECT_FALSE(heartbeat.FinalizeAfterAllRanksDone(0.0));
	EXPECT_TRUE(PathExists(blocked_final_path));
}

TEST_F(SnapshotIOTest, FinalizationRetriesMergedCheckpointCleanup)
{
	const uint64_t run_id = 99003;
	const double interval = 10.0;
	const std::string checkpoint_path = SnapshotRankCheckpointPath(
		rank_snapshot_dir, 0, 1, interval);
	ASSERT_TRUE(WriteSnapshotRankState(checkpoint_path, MakeMergeState(run_id, 0)));

	const SnapshotMergeResult initial_merge = TryWriteSnapshot(
		snapshot_root, rank_snapshot_dir, 1, interval, 1, run_id, 0.5, 1.0e-40, false);
	ASSERT_EQ(SnapshotMergeStatus::Merged, initial_merge.status);
	ASSERT_TRUE(initial_merge.cleanup_succeeded);
	ASSERT_FALSE(PathExists(checkpoint_path));

	ASSERT_EQ(0, mkdir(checkpoint_path.c_str(), 0700));
	{
		std::ofstream blocker(checkpoint_path + "/keep");
		ASSERT_TRUE(blocker.is_open());
		blocker << "keep";
	}

	SnapshotSharedState shared_state;
	shared_state.Initialize(run_id, 0);
	SnapshotHeartbeat heartbeat(
		shared_state, 0, 1, run_id, snapshot_root, rank_snapshot_dir, interval, 1.0, 1.0e-40);

	EXPECT_FALSE(heartbeat.FinalizeAfterAllRanksDone(interval));
	EXPECT_TRUE(PathExists(checkpoint_path));

	RemoveTree(checkpoint_path);
	ASSERT_FALSE(PathExists(checkpoint_path));
	EXPECT_TRUE(heartbeat.FinalizeAfterAllRanksDone(interval));
	EXPECT_FALSE(PathExists(rank_snapshot_dir));
}

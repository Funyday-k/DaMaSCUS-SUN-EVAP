#include <gtest/gtest.h>

#include <chrono>
#include <cstdio>
#include <dirent.h>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>

#include "Snapshot_Heartbeat.hpp"
#include "Snapshot_IO.hpp"

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
			if(name != "." && name != "..")
				RemoveTree(path + "/" + name);
		}
		closedir(directory);
	}
	rmdir(path.c_str());
}

std::string ReadAll(const std::string& path)
{
	std::ifstream stream(path);
	std::ostringstream contents;
	contents << stream.rdbuf();
	return contents.str();
}

class SnapshotHeartbeatRetry : public ::testing::Test
{
  protected:
	void SetUp() override
	{
		char root_template[] = "/tmp/damascus_snapshot_retry_XXXXXX";
		char* created = mkdtemp(root_template);
		ASSERT_NE(created, nullptr);
		root_ = std::string(created) + "/";
		ranks_ = root_ + "ranks/";
		ASSERT_EQ(mkdir(ranks_.c_str(), 0700), 0);
		shared_.Initialize(777, 0);
		heartbeat_.reset(new SnapshotHeartbeat(
			shared_, 0, 2, 777, root_, ranks_, 1.0, 1.0, 1.0e-40));
	}

	void TearDown() override
	{
		if(heartbeat_)
			heartbeat_->Stop();
		if(!root_.empty())
			RemoveTree(root_);
	}

	SnapshotRankState RankState(int rank) const
	{
		SnapshotRankState state;
		state.run_id = 777;
		state.rank = rank;
		state.snapshot_index = 1;
		state.rank_elapsed_wall_sec = 1.0;
		state.local_total = 1;
		state.local_classified = 1;
		return state;
	}

	void Start()
	{
		epoch_ = std::chrono::steady_clock::now();
		ASSERT_TRUE(heartbeat_->Start(epoch_));
	}

	void SleepUntil(int milliseconds) const
	{
		std::this_thread::sleep_until(epoch_ + std::chrono::milliseconds(milliseconds));
	}

	template <typename Predicate>
	bool WaitUntil(Predicate predicate, int milliseconds) const
	{
		const auto deadline = epoch_ + std::chrono::milliseconds(milliseconds);
		while(std::chrono::steady_clock::now() < deadline)
		{
			if(predicate())
				return true;
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
		return predicate();
	}

	std::string root_;
	std::string ranks_;
	SnapshotSharedState shared_;
	std::unique_ptr<SnapshotHeartbeat> heartbeat_;
	std::chrono::steady_clock::time_point epoch_;
};

TEST_F(SnapshotHeartbeatRetry, DelayedShardMergesDuringRunAfterInitialRetry)
{
	Start();
	SleepUntil(2250);
	const std::string report_path = SnapshotTextFilePath(root_, 1, 1.0);
	ASSERT_NE(ReadAll(report_path).find("# snapshot_status = partial"), std::string::npos);
	ASSERT_TRUE(WriteSnapshotRankState(
		SnapshotRankCheckpointPath(ranks_, 1, 1, 1.0), RankState(1)));

	// Inspect while the worker is still running. Neither Stop nor finalization
	// may repair the report on behalf of the periodic retry path being tested.
	ASSERT_TRUE(WaitUntil([&] {
		return ReadAll(report_path).find("# snapshot_status = merged") != std::string::npos;
	}, 4250));
	const std::string report = ReadAll(report_path);
	EXPECT_NE(report.find("# ready_ranks = 0,1"), std::string::npos);
	EXPECT_NE(report.find("# total_trajectories = 1\n"), std::string::npos);
	EXPECT_TRUE(WaitUntil([&] {
		return !PathExists(SnapshotRankCheckpointPath(ranks_, 0, 1, 1.0))
		       && !PathExists(SnapshotRankCheckpointPath(ranks_, 1, 1, 1.0));
	}, 4250));
}

TEST_F(SnapshotHeartbeatRetry, FailedCleanupIsRetriedWithoutLosingMergedReport)
{
	for(int rank = 0; rank < 2; ++rank)
		ASSERT_TRUE(WriteSnapshotRankState(
			SnapshotRankCheckpointPath(ranks_, rank, 1, 1.0), RankState(rank)));
	const SnapshotMergeResult initial = TryWriteSnapshot(
		root_, ranks_, 1, 1.0, 2, 777, 1.0, 1.0e-40, false);
	ASSERT_EQ(initial.status, SnapshotMergeStatus::Merged);
	ASSERT_TRUE(initial.cleanup_succeeded);
	const std::string report_path = SnapshotTextFilePath(root_, 1, 1.0);
	const std::string merged_report = ReadAll(report_path);

	// A nonempty directory at a checkpoint path fails remove() even when the
	// test runs with elevated privileges, unlike permission-based fixtures.
	const std::string obstacle = SnapshotRankCheckpointPath(ranks_, 1, 1, 1.0);
	ASSERT_EQ(mkdir(obstacle.c_str(), 0700), 0);
	{
		std::ofstream child(obstacle + "/keep");
		child << "temporarily prevents checkpoint cleanup";
		ASSERT_TRUE(child.good());
	}
	Start();
	SleepUntil(2250);
	ASSERT_TRUE(PathExists(obstacle));
	EXPECT_EQ(ReadAll(report_path), merged_report);
	RemoveTree(obstacle);
	{
		std::ofstream checkpoint(obstacle);
		checkpoint << "cleanup can now remove this leftover";
		ASSERT_TRUE(checkpoint.good());
	}
	ASSERT_TRUE(WaitUntil([&] { return !PathExists(obstacle); }, 4250));
	EXPECT_EQ(ReadAll(report_path), merged_report);
}
}

#ifndef __Snapshot_Heartbeat_hpp_
#define __Snapshot_Heartbeat_hpp_

#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <thread>

#include "Snapshot_Shared_State.hpp"

namespace DaMaSCUS_SUN
{

class SnapshotHeartbeat
{
  public:
	SnapshotHeartbeat(
		SnapshotSharedState& shared_state,
		int mpi_rank,
		int mpi_processes,
		uint64_t run_id,
		const std::string& snapshot_root,
		const std::string& rank_snapshot_dir,
		double snapshot_interval,
		double mass_gev,
		double sigma_cm2);
	~SnapshotHeartbeat();

	bool Start(const std::chrono::steady_clock::time_point& epoch);
	void MarkDoneAndWriteFinal(double computing_time_sec);
	bool FinalizeAfterAllRanksDone(double final_elapsed_sec);
	void Stop();

  private:
	void ThreadMain() noexcept;
	void RunThread();
	bool WriteRankCheckpoint(
		int snapshot_index,
		double target_wall_sec,
		double actual_elapsed_sec,
		bool& deadline_missed);
	SnapshotMergeCache& AcquireMergeCache(int snapshot_index);
	SnapshotMergeResult TryMergeIndex(int snapshot_index, bool allow_partial);
	void ProcessMergeTick(int snapshot_index, bool local_checkpoint_missed);
	double ElapsedSinceStart() const;

	SnapshotSharedState& shared_state_;
	int mpi_rank_;
	int mpi_processes_;
	uint64_t run_id_;
	std::string snapshot_root_;
	std::string rank_snapshot_dir_;
	double snapshot_interval_;
	double mass_gev_;
	double sigma_cm2_;
	double missed_tolerance_sec_;

	std::thread worker_;
	std::mutex control_mutex_;
	std::condition_variable cv_;
	bool started_ = false;
	bool stopping_ = false;
	bool stop_requested_ = false;
	bool worker_failed_ = false;
	std::chrono::steady_clock::time_point start_time_;

	size_t first_uncommitted_evaporation_entry_ = 0;
	int last_rank_snapshot_written_ = 0;
	int highest_snapshot_index_seen_ = 0;
	std::set<int> retry_merge_indices_;
	std::set<int> unresolved_merge_indices_;
	std::map<int, SnapshotMergeCache> merge_caches_;
};

}	// namespace DaMaSCUS_SUN

#endif

#include "Snapshot_Heartbeat.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <vector>

namespace DaMaSCUS_SUN
{
namespace
{
constexpr size_t MAX_SNAPSHOT_MERGE_CACHE_ENTRIES = 32;
}

SnapshotHeartbeat::SnapshotHeartbeat(
	SnapshotSharedState& shared_state,
	int mpi_rank,
	int mpi_processes,
	uint64_t run_id,
	const std::string& snapshot_root,
	const std::string& rank_snapshot_dir,
	double snapshot_interval,
	double mass_gev,
	double sigma_cm2)
: shared_state_(shared_state),
  mpi_rank_(mpi_rank),
  mpi_processes_(mpi_processes),
  run_id_(run_id),
  snapshot_root_(snapshot_root),
  rank_snapshot_dir_(rank_snapshot_dir),
  snapshot_interval_(snapshot_interval),
  mass_gev_(mass_gev),
  sigma_cm2_(sigma_cm2),
  missed_tolerance_sec_(std::max(1.0, 0.1 * snapshot_interval))
{
	if(!IsValidSnapshotIntervalSeconds(snapshot_interval_))
		throw std::invalid_argument("SnapshotHeartbeat: snapshot interval must be a positive integer number of seconds.");
}

SnapshotHeartbeat::~SnapshotHeartbeat()
{
	Stop();
}

bool SnapshotHeartbeat::Start(const std::chrono::steady_clock::time_point& epoch)
{
	std::lock_guard<std::mutex> lock(control_mutex_);
	if(started_)
		return true;
	start_time_ = epoch;
	stopping_ = false;
	stop_requested_ = false;
	worker_failed_ = false;
	retry_merge_indices_.clear();
	unresolved_merge_indices_.clear();
	merge_caches_.clear();
	highest_snapshot_index_seen_ = 0;
	first_uncommitted_evaporation_entry_ = 0;
	last_rank_snapshot_written_ = 0;
	try
	{
		worker_ = std::thread(&SnapshotHeartbeat::ThreadMain, this);
		started_ = true;
		return true;
	}
	catch(const std::exception& error)
	{
		worker_failed_ = true;
		stop_requested_ = true;
		std::cerr << "Warning in SnapshotHeartbeat::Start(): rank " << mpi_rank_
		          << " failed to start heartbeat thread: " << error.what() << std::endl;
		return false;
	}
	catch(...)
	{
		worker_failed_ = true;
		stop_requested_ = true;
		std::cerr << "Warning in SnapshotHeartbeat::Start(): rank " << mpi_rank_
		          << " failed to start heartbeat thread with an unknown exception." << std::endl;
		return false;
	}
}

void SnapshotHeartbeat::MarkDoneAndWriteFinal(double computing_time_sec)
{
	Stop();
	try
	{
		SnapshotRankState final_state = shared_state_.CopyFinal(computing_time_sec, first_uncommitted_evaporation_entry_);
		final_state.snapshot_index = last_rank_snapshot_written_;
		if(!WriteSnapshotRankState(SnapshotRankFinalPath(rank_snapshot_dir_, mpi_rank_), final_state))
		{
			std::cerr << "Warning in SnapshotHeartbeat::MarkDoneAndWriteFinal(): rank "
			          << mpi_rank_ << " failed to write final snapshot state." << std::endl;
		}
	}
	catch(const std::exception& error)
	{
		std::cerr << "Warning in SnapshotHeartbeat::MarkDoneAndWriteFinal(): rank "
		          << mpi_rank_ << " failed with exception: " << error.what() << std::endl;
	}
	catch(...)
	{
		std::cerr << "Warning in SnapshotHeartbeat::MarkDoneAndWriteFinal(): rank "
		          << mpi_rank_ << " failed with an unknown exception." << std::endl;
	}
}

bool SnapshotHeartbeat::FinalizeAfterAllRanksDone(double final_elapsed_sec)
{
	try
	{
		if(!std::isfinite(final_elapsed_sec) || final_elapsed_sec < 0.0)
			throw std::invalid_argument("final elapsed time is not finite and non-negative");

		const double final_index_value = std::floor(final_elapsed_sec / snapshot_interval_);
		if(final_index_value > static_cast<double>(std::numeric_limits<int>::max()))
			throw std::overflow_error("final snapshot index exceeds int range");
		const int final_snapshot_index = std::max(0, static_cast<int>(final_index_value));
		const int highest_target_index = std::max(final_snapshot_index, highest_snapshot_index_seen_);

		std::set<int> final_targets = retry_merge_indices_;
		final_targets.insert(unresolved_merge_indices_.begin(), unresolved_merge_indices_.end());
		if(highest_snapshot_index_seen_ < highest_target_index)
		{
			for(int index = highest_snapshot_index_seen_ + 1; ; index++)
			{
				final_targets.insert(index);
				if(index == highest_target_index)
					break;
			}
		}

		bool all_targets_merged = true;
		for(const int index : final_targets)
		{
			if(index <= 0 || index > highest_target_index)
				continue;
			const SnapshotMergeResult result = TryMergeIndex(index, false);
			if(result.status == SnapshotMergeStatus::Merged && result.cleanup_succeeded)
			{
				retry_merge_indices_.erase(index);
				unresolved_merge_indices_.erase(index);
			}
			else
			{
				unresolved_merge_indices_.insert(index);
				all_targets_merged = false;
			}
		}

		if(all_targets_merged)
		{
			if(CleanupFinalSnapshotStates(rank_snapshot_dir_, mpi_processes_))
				return true;
			std::cerr << "Warning in SnapshotHeartbeat::FinalizeAfterAllRanksDone(): "
			             "all reports merged, but final rank state cleanup failed."
			          << std::endl;
			return false;
		}

		std::cerr << "Warning in SnapshotHeartbeat::FinalizeAfterAllRanksDone(): "
		          << unresolved_merge_indices_.size()
		          << " snapshot index(es) remain incomplete; preserving final rank states for retry and diagnosis."
		          << std::endl;
		return false;
	}
	catch(const std::exception& error)
	{
		std::cerr << "Warning in SnapshotHeartbeat::FinalizeAfterAllRanksDone(): "
		          << error.what() << "; preserving final rank states for retry and diagnosis." << std::endl;
		return false;
	}
	catch(...)
	{
		std::cerr << "Warning in SnapshotHeartbeat::FinalizeAfterAllRanksDone(): unknown exception; "
		             "preserving final rank states for retry and diagnosis." << std::endl;
		return false;
	}
}

void SnapshotHeartbeat::Stop()
{
	std::thread worker_to_join;
	{
		std::unique_lock<std::mutex> lock(control_mutex_);
		if(!started_)
			return;
		if(stopping_)
		{
			cv_.wait(lock, [&] { return !stopping_; });
			return;
		}
		stopping_ = true;
		stop_requested_ = true;
		worker_to_join = std::move(worker_);
	}
	cv_.notify_all();
	if(worker_to_join.joinable())
		worker_to_join.join();
	{
		std::lock_guard<std::mutex> lock(control_mutex_);
		started_ = false;
		stopping_ = false;
	}
	cv_.notify_all();
}

void SnapshotHeartbeat::ThreadMain() noexcept
{
	try
	{
		RunThread();
	}
	catch(const std::exception& error)
	{
		{
			std::lock_guard<std::mutex> lock(control_mutex_);
			worker_failed_ = true;
			stop_requested_ = true;
		}
		std::cerr << "Warning in SnapshotHeartbeat::ThreadMain(): rank " << mpi_rank_
		          << " heartbeat stopped after exception: " << error.what() << std::endl;
		cv_.notify_all();
	}
	catch(...)
	{
		{
			std::lock_guard<std::mutex> lock(control_mutex_);
			worker_failed_ = true;
			stop_requested_ = true;
		}
		std::cerr << "Warning in SnapshotHeartbeat::ThreadMain(): rank " << mpi_rank_
		          << " heartbeat stopped after an unknown exception." << std::endl;
		cv_.notify_all();
	}
}

void SnapshotHeartbeat::RunThread()
{
	int next_index = 1;
	while(true)
	{
		const auto deadline = start_time_ + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
			std::chrono::duration<double>(next_index * snapshot_interval_));

		{
			std::unique_lock<std::mutex> lock(control_mutex_);
			if(cv_.wait_until(lock, deadline, [&] { return stop_requested_; }))
				break;
		}

		const double target_elapsed = next_index * snapshot_interval_;
		const double actual_elapsed = ElapsedSinceStart();
		bool checkpoint_missed = actual_elapsed > target_elapsed + missed_tolerance_sec_;
		if(checkpoint_missed)
		{
			if(mpi_rank_ == 0)
				WriteMissedSnapshotMarker(snapshot_root_, next_index, snapshot_interval_, run_id_, actual_elapsed);
		}
		else
		{
			bool missed_while_copying = false;
			WriteRankCheckpoint(next_index, target_elapsed, actual_elapsed, missed_while_copying);
			checkpoint_missed = missed_while_copying;
			if(checkpoint_missed && mpi_rank_ == 0)
				WriteMissedSnapshotMarker(snapshot_root_, next_index, snapshot_interval_, run_id_, ElapsedSinceStart());
		}

		if(mpi_rank_ == 0)
			ProcessMergeTick(next_index, checkpoint_missed);
		if(next_index == std::numeric_limits<int>::max())
			throw std::overflow_error("snapshot index exceeds int range");
		next_index++;
	}
}

bool SnapshotHeartbeat::WriteRankCheckpoint(
	int snapshot_index,
	double target_wall_sec,
	double actual_elapsed_sec,
	bool& deadline_missed)
{
	deadline_missed = false;
	size_t evaporation_entry_end = first_uncommitted_evaporation_entry_;
	// Publish at most one budget worth of completed events per checkpoint. The
	// remainder is published by a later checkpoint, whose window filter admits
	// earlier completion times, so a backlog drains instead of ratcheting up
	// until every write exceeds the size the reader accepts.
	SnapshotRankState state = shared_state_.CopyForSnapshot(
		snapshot_index,
		actual_elapsed_sec,
		target_wall_sec,
		first_uncommitted_evaporation_entry_,
		evaporation_entry_end,
		SnapshotEvaporationEventsPerCheckpoint());
	state.snapshot_index = snapshot_index;
	state.rank_elapsed_wall_sec = ElapsedSinceStart();
	if(state.rank_elapsed_wall_sec > target_wall_sec + missed_tolerance_sec_)
	{
		deadline_missed = true;
		return false;
	}

	const std::string checkpoint_path = SnapshotRankCheckpointPath(rank_snapshot_dir_, mpi_rank_, snapshot_index, snapshot_interval_);
	if(WriteSnapshotRankState(checkpoint_path, state))
	{
		first_uncommitted_evaporation_entry_ = evaporation_entry_end;
		last_rank_snapshot_written_ = snapshot_index;
		return true;
	}

	std::cerr << "Warning in SnapshotHeartbeat::WriteRankCheckpoint(): rank "
	          << mpi_rank_ << " failed to write " << checkpoint_path << std::endl;
	return false;
}

SnapshotMergeCache& SnapshotHeartbeat::AcquireMergeCache(int snapshot_index)
{
	SnapshotMergeCache& cache = merge_caches_[snapshot_index];
	// Dropping an accumulation only costs a full reload of that index later, so
	// the map may be bounded without affecting the merged result.
	while(merge_caches_.size() > MAX_SNAPSHOT_MERGE_CACHE_ENTRIES)
	{
		auto lowest = merge_caches_.begin();
		if(lowest->first == snapshot_index)
			lowest = std::next(lowest);
		if(lowest == merge_caches_.end())
			break;
		merge_caches_.erase(lowest);
	}
	return cache;
}

SnapshotMergeResult SnapshotHeartbeat::TryMergeIndex(int snapshot_index, bool allow_partial)
{
	SnapshotMergeResult result = TryWriteSnapshotCached(
		snapshot_root_,
		rank_snapshot_dir_,
		snapshot_index,
		snapshot_interval_,
		mpi_processes_,
		run_id_,
		mass_gev_,
		sigma_cm2_,
		AcquireMergeCache(snapshot_index),
		allow_partial);
	if(result.status == SnapshotMergeStatus::Merged)
		merge_caches_.erase(snapshot_index);
	return result;
}

void SnapshotHeartbeat::ProcessMergeTick(int snapshot_index, bool local_checkpoint_missed)
{
	const std::vector<int> previous_retries(retry_merge_indices_.begin(), retry_merge_indices_.end());
	if(local_checkpoint_missed)
		unresolved_merge_indices_.insert(snapshot_index);
	else
		retry_merge_indices_.insert(snapshot_index);
	highest_snapshot_index_seen_ = std::max(highest_snapshot_index_seen_, snapshot_index);

	for(const int index : previous_retries)
	{
		const SnapshotMergeResult result = TryMergeIndex(index, true);
		if(result.status == SnapshotMergeStatus::Merged && result.cleanup_succeeded)
		{
			retry_merge_indices_.erase(index);
			unresolved_merge_indices_.erase(index);
		}
		else if(result.status != SnapshotMergeStatus::Merged)
		{
			unresolved_merge_indices_.insert(index);
			retry_merge_indices_.erase(index);
		}
	}

	if(local_checkpoint_missed)
		return;

	const SnapshotMergeResult result = TryMergeIndex(snapshot_index, true);
	if(result.status == SnapshotMergeStatus::Merged && result.cleanup_succeeded)
		retry_merge_indices_.erase(snapshot_index);
}

double SnapshotHeartbeat::ElapsedSinceStart() const
{
	return 1.0e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(
		std::chrono::steady_clock::now() - start_time_).count();
}

}	// namespace DaMaSCUS_SUN

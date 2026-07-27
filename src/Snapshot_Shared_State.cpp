#include "Snapshot_Shared_State.hpp"

#include <algorithm>
#include <cmath>

namespace DaMaSCUS_SUN
{
namespace
{
bool IsNumericalTermination(TrajectoryTerminationReason reason)
{
	return reason == TrajectoryTerminationReason::NumericalFailure
	    || reason == TrajectoryTerminationReason::NonFiniteState
	    || reason == TrajectoryTerminationReason::SpeedLimit
	    || reason == TrajectoryTerminationReason::EnergyDriftEscape
	    || reason == TrajectoryTerminationReason::Unknown;
}
}

void SnapshotSharedState::Initialize(uint64_t run_id, int rank)
{
	std::lock_guard<std::mutex> lock(mutex_);
	run_id_ = run_id;
	rank_ = rank;
	trajectory_in_progress_ = false;
	current_trajectory_captured_ = false;
	current_trajectory_id_ = 0;
	current_trajectory_wall_start_ = std::chrono::steady_clock::time_point();
	current_trajectory_simulation_start_sec_ = 0.0;
	current_trajectory_simulated_elapsed_sec_ = 0.0;
	current_trajectory_scatterings_ = 0;
	current_dt_hist_.fill(0.0);
	current_v2dt_hist_.fill(0.0);
	completed_trajectories_ = 0;
	captured_particles_ = 0;
	classified_trajectories_ = 0;
	numerical_failures_ = 0;
	bincount_captured_samples_ = 0;
	bincount_not_captured_samples_ = 0;
	captured_dt_hist_.fill(0.0);
	captured_v2dt_hist_.fill(0.0);
	captured_dt_sq_hist_.fill(0.0);
	captured_v2dt_sq_hist_.fill(0.0);
	not_captured_dt_hist_.fill(0.0);
	not_captured_v2dt_hist_.fill(0.0);
	not_captured_dt_sq_hist_.fill(0.0);
	not_captured_v2dt_sq_hist_.fill(0.0);
	evaporation_events_.clear();
}

void SnapshotSharedState::BeginTrajectory(uint64_t trajectory_id, double initial_simulated_time_sec)
{
	std::lock_guard<std::mutex> lock(mutex_);
	trajectory_in_progress_ = true;
	current_trajectory_captured_ = false;
	current_trajectory_id_ = trajectory_id;
	current_trajectory_wall_start_ = std::chrono::steady_clock::now();
	current_trajectory_simulation_start_sec_ = std::isfinite(initial_simulated_time_sec) ? initial_simulated_time_sec : 0.0;
	current_trajectory_simulated_elapsed_sec_ = 0.0;
	current_trajectory_scatterings_ = 0;
	current_dt_hist_.fill(0.0);
	current_v2dt_hist_.fill(0.0);
}

void SnapshotSharedState::AddCurrentBincountInterval(
	const std::vector<BincountContribution>& contributions,
	double simulated_time_sec)
{
	std::lock_guard<std::mutex> lock(mutex_);
	if(!trajectory_in_progress_)
		return;
	if(std::isfinite(simulated_time_sec))
		current_trajectory_simulated_elapsed_sec_ =
			std::max(0.0, simulated_time_sec - current_trajectory_simulation_start_sec_);
	for(const BincountContribution& contribution : contributions)
	{
		if(contribution.bin < 0
		   || !std::isfinite(contribution.dt_sec) || contribution.dt_sec <= 0.0
		   || !std::isfinite(contribution.v2dt_km2_per_sec)
		   || contribution.v2dt_km2_per_sec < 0.0)
			continue;
		if(static_cast<std::size_t>(contribution.bin) >= TOTAL_BINS)
			continue;
		current_dt_hist_[contribution.bin] += contribution.dt_sec;
		current_v2dt_hist_[contribution.bin] += contribution.v2dt_km2_per_sec;
	}
}

void SnapshotSharedState::UpdateCurrentSimulationTime(double simulated_time_sec)
{
	std::lock_guard<std::mutex> lock(mutex_);
	if(!trajectory_in_progress_ || !std::isfinite(simulated_time_sec))
		return;
	current_trajectory_simulated_elapsed_sec_ =
		std::max(0.0, simulated_time_sec - current_trajectory_simulation_start_sec_);
}

void SnapshotSharedState::UpdateCurrentScatterings(uint64_t scatterings)
{
	std::lock_guard<std::mutex> lock(mutex_);
	if(trajectory_in_progress_)
		current_trajectory_scatterings_ = scatterings;
}

void SnapshotSharedState::MarkCurrentCaptured(bool captured)
{
	std::lock_guard<std::mutex> lock(mutex_);
	if(trajectory_in_progress_)
		current_trajectory_captured_ = captured;
}

void SnapshotSharedState::PublishCurrentTrajectoryProgress(
	const std::array<double, TOTAL_BINS>& dt_hist,
	const std::array<double, TOTAL_BINS>& v2dt_hist,
	double simulated_time_sec)
{
	std::lock_guard<std::mutex> lock(mutex_);
	if(!trajectory_in_progress_)
		return;
	if(std::isfinite(simulated_time_sec))
		current_trajectory_simulated_elapsed_sec_ =
			std::max(0.0, simulated_time_sec - current_trajectory_simulation_start_sec_);
	// Published snapshots must satisfy the rank-state validity invariants, so a
	// non-finite or negative bin is dropped exactly as the per-step path did.
	for(std::size_t bin = 0; bin < TOTAL_BINS; bin++)
	{
		const double dt = dt_hist[bin];
		const double v2dt = v2dt_hist[bin];
		current_dt_hist_[bin] = (std::isfinite(dt) && dt > 0.0) ? dt : 0.0;
		current_v2dt_hist_[bin] = (std::isfinite(v2dt) && v2dt >= 0.0) ? v2dt : 0.0;
	}
}

void SnapshotSharedState::RecordCompletedTrajectory(
	const TrajectoryBincount& bincount,
	bool count_as_captured_bincount_sample,
	bool count_as_not_captured_bincount_sample,
	const std::vector<SnapshotEvaporationProgressEntry>& new_evaporation_events)
{
	std::lock_guard<std::mutex> lock(mutex_);
	completed_trajectories_++;
	if(bincount.is_captured)
		captured_particles_++;
	if(bincount.is_captured || count_as_not_captured_bincount_sample)
		classified_trajectories_++;
	if(IsNumericalTermination(bincount.termination_reason))
		numerical_failures_++;

	if(count_as_captured_bincount_sample)
	{
		bincount_captured_samples_++;
		for(std::size_t bin = 0; bin < TOTAL_BINS; bin++)
		{
			captured_dt_hist_[bin] += bincount.dt_hist[bin];
			captured_v2dt_hist_[bin] += bincount.v2dt_hist[bin];
			captured_dt_sq_hist_[bin] += bincount.dt_hist[bin] * bincount.dt_hist[bin];
			captured_v2dt_sq_hist_[bin] += bincount.v2dt_hist[bin] * bincount.v2dt_hist[bin];
		}
	}

	if(count_as_not_captured_bincount_sample)
	{
		bincount_not_captured_samples_++;
		for(std::size_t bin = 0; bin < TOTAL_BINS; bin++)
		{
			not_captured_dt_hist_[bin] += bincount.dt_hist[bin];
			not_captured_v2dt_hist_[bin] += bincount.v2dt_hist[bin];
			not_captured_dt_sq_hist_[bin] += bincount.dt_hist[bin] * bincount.dt_hist[bin];
			not_captured_v2dt_sq_hist_[bin] += bincount.v2dt_hist[bin] * bincount.v2dt_hist[bin];
		}
	}

	evaporation_events_.insert(evaporation_events_.end(), new_evaporation_events.begin(), new_evaporation_events.end());
	ClearCurrentTrajectoryLocked();
}

SnapshotRankState SnapshotSharedState::CopyForSnapshot(
	int snapshot_index,
	double rank_elapsed_wall_sec,
	double target_wall_sec,
	size_t evaporation_begin,
	size_t& evaporation_end,
	size_t max_new_evaporation_events) const
{
	std::lock_guard<std::mutex> lock(mutex_);
	evaporation_begin = std::min(evaporation_begin, evaporation_events_.size());
	evaporation_end = evaporation_begin;
	while(evaporation_end < evaporation_events_.size()
	      && evaporation_events_[evaporation_end].completion_wall_time_sec <= target_wall_sec)
		evaporation_end++;
	// Events are appended in completion order, so a prefix is also a
	// wall-time-ordered prefix. Whatever is left over is picked up by a later
	// checkpoint, whose window filter admits earlier completion times as well.
	if(evaporation_end - evaporation_begin > max_new_evaporation_events)
		evaporation_end = evaporation_begin + max_new_evaporation_events;
	return CopyLocked(snapshot_index, false, rank_elapsed_wall_sec, evaporation_begin, evaporation_end);
}

SnapshotRankState SnapshotSharedState::CopyFinal(double computing_time_sec, size_t evaporation_begin) const
{
	std::lock_guard<std::mutex> lock(mutex_);
	return CopyLocked(-1, true, computing_time_sec, evaporation_begin, evaporation_events_.size());
}

void SnapshotSharedState::ClearCurrentTrajectoryLocked()
{
	trajectory_in_progress_ = false;
	current_trajectory_captured_ = false;
	current_trajectory_id_ = 0;
	current_trajectory_wall_start_ = std::chrono::steady_clock::time_point();
	current_trajectory_simulation_start_sec_ = 0.0;
	current_trajectory_simulated_elapsed_sec_ = 0.0;
	current_trajectory_scatterings_ = 0;
	current_dt_hist_.fill(0.0);
	current_v2dt_hist_.fill(0.0);
}

SnapshotRankState SnapshotSharedState::CopyLocked(
	int snapshot_index,
	bool done,
	double rank_elapsed_wall_sec,
	size_t evaporation_begin,
	size_t evaporation_end) const
{
	SnapshotRankState state;
	state.run_id = run_id_;
	state.rank = rank_;
	state.snapshot_index = snapshot_index;
	state.done = done ? 1 : 0;
	state.trajectory_in_progress = (!done && trajectory_in_progress_) ? 1 : 0;
	state.local_captured = captured_particles_;
	state.local_total = completed_trajectories_;
	state.local_classified = classified_trajectories_;
	state.local_numerical_failures = numerical_failures_;
	state.bincount_captured_samples = bincount_captured_samples_;
	state.bincount_not_captured_samples = bincount_not_captured_samples_;
	state.current_trajectory_id = state.trajectory_in_progress ? current_trajectory_id_ : 0;
	state.rank_elapsed_wall_sec = rank_elapsed_wall_sec;
	state.current_trajectory_captured = (state.trajectory_in_progress && current_trajectory_captured_) ? 1 : 0;
	if(state.trajectory_in_progress)
	{
		state.current_trajectory_wall_sec = 1.0e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(
			std::chrono::steady_clock::now() - current_trajectory_wall_start_).count();
		state.current_trajectory_simulated_elapsed_sec = current_trajectory_simulated_elapsed_sec_;
		state.current_trajectory_scatterings = current_trajectory_scatterings_;
		state.current_trajectory_dt_hist = current_dt_hist_;
		state.current_trajectory_v2dt_hist = current_v2dt_hist_;
	}
	state.captured_dt_hist = captured_dt_hist_;
	state.captured_v2dt_hist = captured_v2dt_hist_;
	state.captured_dt_sq_hist = captured_dt_sq_hist_;
	state.captured_v2dt_sq_hist = captured_v2dt_sq_hist_;
	state.not_captured_dt_hist = not_captured_dt_hist_;
	state.not_captured_v2dt_hist = not_captured_v2dt_hist_;
	state.not_captured_dt_sq_hist = not_captured_dt_sq_hist_;
	state.not_captured_v2dt_sq_hist = not_captured_v2dt_sq_hist_;

	evaporation_begin = std::min(evaporation_begin, evaporation_events_.size());
	evaporation_end = std::min(evaporation_end, evaporation_events_.size());
	state.new_evaporation_events.reserve(evaporation_end - evaporation_begin);
	for(size_t index = evaporation_begin; index < evaporation_end; index++)
		state.new_evaporation_events.push_back(evaporation_events_[index]);
	return state;
}

SnapshotRecorder::SnapshotRecorder(SnapshotSharedState& state)
: state_(state)
{
}

void SnapshotRecorder::BeginTrajectory(uint64_t trajectory_id, double initial_simulated_time_sec)
{
	state_.BeginTrajectory(trajectory_id, initial_simulated_time_sec);
}

void SnapshotRecorder::AddCurrentBincountInterval(
	const std::vector<BincountContribution>& contributions,
	double simulated_time_sec)
{
	state_.AddCurrentBincountInterval(contributions, simulated_time_sec);
}

void SnapshotRecorder::UpdateCurrentSimulationTime(double simulated_time_sec)
{
	state_.UpdateCurrentSimulationTime(simulated_time_sec);
}

void SnapshotRecorder::PublishCurrentTrajectoryProgress(
	const std::array<double, TOTAL_BINS>& dt_hist,
	const std::array<double, TOTAL_BINS>& v2dt_hist,
	double simulated_time_sec)
{
	state_.PublishCurrentTrajectoryProgress(dt_hist, v2dt_hist, simulated_time_sec);
}

void SnapshotRecorder::UpdateCurrentScatterings(uint64_t scatterings)
{
	state_.UpdateCurrentScatterings(scatterings);
}

void SnapshotRecorder::MarkCurrentCaptured(bool captured)
{
	state_.MarkCurrentCaptured(captured);
}

}	// namespace DaMaSCUS_SUN

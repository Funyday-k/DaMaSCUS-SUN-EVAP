#ifndef __Snapshot_Shared_State_hpp_
#define __Snapshot_Shared_State_hpp_

#include <chrono>
#include <limits>
#include <mutex>
#include <vector>

#include "Snapshot_IO.hpp"
#include "Simulation_Trajectory.hpp"

namespace DaMaSCUS_SUN
{

class SnapshotSharedState
{
  public:
	void Initialize(uint64_t run_id, int rank);
	void BeginTrajectory(uint64_t trajectory_id, double initial_simulated_time_sec = 0.0);
	void AddCurrentBincountInterval(
		const std::vector<BincountContribution>& contributions,
		double simulated_time_sec);
	void UpdateCurrentSimulationTime(double simulated_time_sec);
	void UpdateCurrentScatterings(uint64_t scatterings);
	void MarkCurrentCaptured(bool captured);

	// Replace the in-progress trajectory histogram wholesale instead of
	// streaming every integrator step through the mutex. The simulator already
	// keeps the identical per-trajectory accumulation, so publishing it at a
	// coarse cadence costs one lock per publish rather than one lock per step.
	void PublishCurrentTrajectoryProgress(
		const std::array<double, TOTAL_BINS>& dt_hist,
		const std::array<double, TOTAL_BINS>& v2dt_hist,
		double simulated_time_sec);

	void RecordCompletedTrajectory(
		const TrajectoryBincount& bincount,
		bool count_as_residence_sample,
		bool physically_classified_uncaptured,
		const std::vector<SnapshotEvaporationProgressEntry>& new_evaporation_events);

	// max_new_evaporation_events bounds how many completed events a single
	// checkpoint publishes. Committing a bounded prefix keeps an oversized
	// backlog from making every later checkpoint fail permanently.
	SnapshotRankState CopyForSnapshot(
		int snapshot_index,
		double rank_elapsed_wall_sec,
		double target_wall_sec,
		size_t evaporation_begin,
		size_t& evaporation_end,
		size_t max_new_evaporation_events = std::numeric_limits<size_t>::max()) const;

	SnapshotRankState CopyFinal(double computing_time_sec, size_t evaporation_begin) const;

  private:
	void ClearCurrentTrajectoryLocked();
	SnapshotRankState CopyLocked(int snapshot_index, bool done, double rank_elapsed_wall_sec, size_t evaporation_begin, size_t evaporation_end) const;

	mutable std::mutex mutex_;
	uint64_t run_id_ = 0;
	int rank_ = 0;

	bool trajectory_in_progress_ = false;
	bool current_trajectory_captured_ = false;
	uint64_t current_trajectory_id_ = 0;
	std::chrono::steady_clock::time_point current_trajectory_wall_start_{};
	double current_trajectory_simulation_start_sec_ = 0.0;
	double current_trajectory_simulated_elapsed_sec_ = 0.0;
	uint64_t current_trajectory_scatterings_ = 0;
	std::array<double, TOTAL_BINS> current_dt_hist_{};
	std::array<double, TOTAL_BINS> current_v2dt_hist_{};

	uint64_t completed_trajectories_ = 0;
	uint64_t captured_particles_ = 0;
	uint64_t classified_trajectories_ = 0;
	uint64_t numerical_failures_ = 0;
	uint64_t bincount_captured_samples_ = 0;

	std::array<double, TOTAL_BINS> captured_dt_hist_{};
	std::array<double, TOTAL_BINS> captured_v2dt_hist_{};
	std::array<double, TOTAL_BINS> captured_dt_sq_hist_{};
	std::array<double, TOTAL_BINS> captured_v2dt_sq_hist_{};

	std::vector<SnapshotEvaporationProgressEntry> evaporation_events_;
};

class SnapshotRecorder
{
  public:
	explicit SnapshotRecorder(SnapshotSharedState& state);

	void BeginTrajectory(uint64_t trajectory_id, double initial_simulated_time_sec = 0.0);
	void AddCurrentBincountInterval(
		const std::vector<BincountContribution>& contributions,
		double simulated_time_sec);
	void UpdateCurrentSimulationTime(double simulated_time_sec);
	void PublishCurrentTrajectoryProgress(
		const std::array<double, TOTAL_BINS>& dt_hist,
		const std::array<double, TOTAL_BINS>& v2dt_hist,
		double simulated_time_sec);
	void UpdateCurrentScatterings(uint64_t scatterings);
	void MarkCurrentCaptured(bool captured);

  private:
	SnapshotSharedState& state_;
};

}	// namespace DaMaSCUS_SUN

#endif

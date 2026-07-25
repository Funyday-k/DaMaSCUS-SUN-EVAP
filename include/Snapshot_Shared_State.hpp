#ifndef __Snapshot_Shared_State_hpp_
#define __Snapshot_Shared_State_hpp_

#include <chrono>
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

	void RecordCompletedTrajectory(
		const TrajectoryBincount& bincount,
		bool count_as_captured_bincount_sample,
		bool count_as_not_captured_bincount_sample,
		const std::vector<SnapshotEvaporationProgressEntry>& new_evaporation_events);

	SnapshotRankState CopyForSnapshot(
		int snapshot_index,
		double rank_elapsed_wall_sec,
		double target_wall_sec,
		size_t evaporation_begin,
		size_t& evaporation_end) const;

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
	std::array<double, NUM_BINS> current_dt_hist_{};
	std::array<double, NUM_BINS> current_v2dt_hist_{};

	uint64_t completed_trajectories_ = 0;
	uint64_t captured_particles_ = 0;
	uint64_t classified_trajectories_ = 0;
	uint64_t numerical_failures_ = 0;
	uint64_t bincount_captured_samples_ = 0;
	uint64_t bincount_not_captured_samples_ = 0;

	std::array<double, NUM_BINS> captured_dt_hist_{};
	std::array<double, NUM_BINS> captured_v2dt_hist_{};
	std::array<double, NUM_BINS> captured_dt_sq_hist_{};
	std::array<double, NUM_BINS> captured_v2dt_sq_hist_{};
	std::array<double, NUM_BINS> not_captured_dt_hist_{};
	std::array<double, NUM_BINS> not_captured_v2dt_hist_{};
	std::array<double, NUM_BINS> not_captured_dt_sq_hist_{};
	std::array<double, NUM_BINS> not_captured_v2dt_sq_hist_{};

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
	void UpdateCurrentScatterings(uint64_t scatterings);
	void MarkCurrentCaptured(bool captured);

  private:
	SnapshotSharedState& state_;
};

}	// namespace DaMaSCUS_SUN

#endif

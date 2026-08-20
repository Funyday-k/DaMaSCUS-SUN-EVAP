#ifndef __Snapshot_IO_hpp_
#define __Snapshot_IO_hpp_

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "Data_Generation.hpp"

namespace DaMaSCUS_SUN
{

struct SnapshotEvaporationProgressEntry
{
	uint64_t trajectory_id = 0;
	double completion_wall_time_sec = 0.0;
	double lifetime_unbinding_sec = -1.0;
	double r_capture_rsun = -1.0;
	double E_capture_eV = 0.0;
	double dE_capture_eV = 0.0;
};

struct SnapshotRankedEvaporationEntry
{
	int rank = -1;
	SnapshotEvaporationProgressEntry entry;
};

struct SnapshotRankState
{
	uint64_t run_id = 0;
	int32_t snapshot_index = 0;
	int32_t rank = 0;
	int32_t done = 0;
	int32_t trajectory_in_progress = 0;
	uint64_t local_captured = 0;
	uint64_t local_total = 0;
	uint64_t local_classified = 0;
	uint64_t local_numerical_failures = 0;
	uint64_t bincount_captured_samples = 0;
	uint64_t current_trajectory_id = 0;
	double rank_elapsed_wall_sec = 0.0;
	double current_trajectory_wall_sec = 0.0;
	double current_trajectory_simulated_elapsed_sec = 0.0;
	uint64_t current_trajectory_scatterings = 0;
	int32_t current_trajectory_captured = 0;
	std::array<double, TOTAL_BINS> current_trajectory_dt_hist{};
	std::array<double, TOTAL_BINS> current_trajectory_v2dt_hist{};
	std::array<double, TOTAL_BINS> captured_dt_hist{};
	std::array<double, TOTAL_BINS> captured_v2dt_hist{};
	std::array<double, TOTAL_BINS> captured_dt_sq_hist{};
	std::array<double, TOTAL_BINS> captured_v2dt_sq_hist{};
	std::vector<SnapshotEvaporationProgressEntry> new_evaporation_events;
};

enum class SnapshotMergeStatus
{
	NoRanksReady,
	Partial,
	Merged
};

struct SnapshotMergeResult
{
	SnapshotMergeStatus status = SnapshotMergeStatus::NoRanksReady;
	bool cleanup_succeeded = true;
	std::vector<int> ready_ranks;
	std::vector<int> missing_ranks;
};

struct SnapshotReportState
{
	struct RankProgress
	{
		int rank = -1;
		int done = 0;
		int trajectory_in_progress = 0;
		int current_trajectory_captured = 0;
		uint64_t current_trajectory_id = 0;
		uint64_t current_trajectory_scatterings = 0;
		double rank_elapsed_wall_sec = 0.0;
		double current_trajectory_wall_sec = 0.0;
		double current_trajectory_simulated_elapsed_sec = 0.0;
	};

	int snapshot_index = 0;
	long long snapshot_time_label = 0;
	double snapshot_interval_seconds = 0.0;
	uint64_t total_trajectories = 0;
	uint64_t captured_particles = 0;
	uint64_t classified_trajectories = 0;
	uint64_t numerical_failures = 0;
	uint64_t snapshot_bincount_captured_samples = 0;
	uint64_t in_progress_bincount_captured_samples = 0;
	std::array<double, TOTAL_BINS> captured_dt_hist{};
	std::array<double, TOTAL_BINS> captured_v2dt_hist{};
	std::array<double, TOTAL_BINS> captured_dt_sq_hist{};
	std::array<double, TOTAL_BINS> captured_v2dt_sq_hist{};
	std::vector<SnapshotRankedEvaporationEntry> new_evaporation_events;
	std::vector<RankProgress> rank_progress;
};

// Rank checkpoint and final files are immutable once atomically renamed, so a
// rank already folded into the report never has to be read again across the
// retries of one snapshot index.
struct SnapshotMergeCache
{
	SnapshotReportState report;
	std::vector<char> rank_accumulated;
	bool initialized = false;
};

long long SnapshotTimeLabelSeconds(int snapshot_index, double interval_seconds);
std::string SnapshotTextFilePath(const std::string& snapshot_root, int snapshot_index, double interval_seconds);
std::string SnapshotEvaporationTimeFilePath(const std::string& snapshot_root, int snapshot_index, double interval_seconds);
std::string SnapshotRankCheckpointPath(const std::string& rank_snapshot_dir, int rank, int snapshot_index, double interval_seconds);
std::string SnapshotRankFinalPath(const std::string& rank_snapshot_dir, int rank);

SnapshotEvaporationProgressEntry MakeSnapshotEvaporationProgressEntry(const CompactEvaporationEvent& event);

// Upper bound on the completed evaporation events a single rank checkpoint may
// publish. A larger backlog is drained over consecutive checkpoints.
size_t SnapshotEvaporationEventsPerCheckpoint();

bool WriteSnapshotRankState(const std::string& path, const SnapshotRankState& state);
bool ReadSnapshotRankState(const std::string& path, uint64_t expected_run_id, SnapshotRankState& state);
bool CleanupSnapshotCheckpoints(const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes);
bool CleanupFinalSnapshotStates(const std::string& rank_snapshot_dir, int mpi_processes);

SnapshotMergeResult TryWriteSnapshot(
	const std::string& snapshot_root,
	const std::string& rank_snapshot_dir,
	int snapshot_index,
	double interval_seconds,
	int mpi_processes,
	uint64_t run_id,
	double mass_gev,
	double sigma_cm2,
	bool allow_partial);

SnapshotMergeResult TryWriteSnapshotCached(
	const std::string& snapshot_root,
	const std::string& rank_snapshot_dir,
	int snapshot_index,
	double interval_seconds,
	int mpi_processes,
	uint64_t run_id,
	double mass_gev,
	double sigma_cm2,
	SnapshotMergeCache& cache,
	bool allow_partial);

bool WriteMissedSnapshotMarker(
	const std::string& snapshot_root,
	int snapshot_index,
	double interval_seconds,
	uint64_t run_id,
	double actual_write_wall_time_sec);

}	// namespace DaMaSCUS_SUN

#endif

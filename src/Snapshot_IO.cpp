#include "Snapshot_IO.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

namespace DaMaSCUS_SUN
{
namespace
{
constexpr uint64_t SNAPSHOT_RANK_STATE_MAGIC = 0x4453534e41503031ULL;
constexpr uint32_t SNAPSHOT_RANK_STATE_VERSION = 6;
constexpr uint32_t SNAPSHOT_RANK_STATE_HEADER_BYTES = sizeof(uint64_t) + 2 * sizeof(uint32_t);
constexpr uint64_t SNAPSHOT_EVAPORATION_EVENTS_PER_CHECKPOINT = 10000000ULL;
// Only an absurd declared event count is rejected here: ReadSnapshotRankState
// cross-checks the count against the real file size and guards the resize with
// try/catch, so a small cap would be redundant defense that instead turns an
// oversized backlog into a permanent write failure.
constexpr uint64_t MAX_SNAPSHOT_EVAPORATION_EVENTS = 1000000000ULL;

uint64_t SnapshotRankStateFixedBytes()
{
	return SNAPSHOT_RANK_STATE_HEADER_BYTES
	     + sizeof(uint64_t)
	     + 4 * sizeof(int32_t)
	     + 8 * sizeof(uint64_t)
	     + 3 * sizeof(double)
	     + sizeof(int32_t)
	     + sizeof(uint64_t)  // fixed radial-bin count
	     + sizeof(uint64_t);
}

uint64_t SnapshotEvaporationEntryBytes()
{
	return sizeof(uint64_t) + 5 * sizeof(double);
}

template<typename T>
void WriteBinaryValue(std::ofstream& file, const T& value)
{
	file.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<typename T, std::size_t N>
void WriteBinaryArray(std::ofstream& file, const std::array<T, N>& values)
{
	file.write(
	    reinterpret_cast<const char*>(values.data()),
	    static_cast<std::streamsize>(N * sizeof(T)));
}

template<typename T>
void ReadBinaryValue(std::ifstream& file, T& value)
{
	file.read(reinterpret_cast<char*>(&value), sizeof(T));
}

template<typename T, std::size_t N>
void ReadBinaryArray(std::ifstream& file, std::array<T, N>& values)
{
	file.read(
	    reinterpret_cast<char*>(values.data()),
	    static_cast<std::streamsize>(N * sizeof(T)));
}

void WriteSnapshotEvaporationEntryBinary(std::ofstream& file, const SnapshotEvaporationProgressEntry& entry)
{
	WriteBinaryValue(file, entry.trajectory_id);
	WriteBinaryValue(file, entry.completion_wall_time_sec);
	WriteBinaryValue(file, entry.lifetime_unbinding_sec);
	WriteBinaryValue(file, entry.r_capture_rsun);
	WriteBinaryValue(file, entry.E_capture_eV);
	WriteBinaryValue(file, entry.dE_capture_eV);
}

void ReadSnapshotEvaporationEntryBinary(std::ifstream& file, SnapshotEvaporationProgressEntry& entry)
{
	ReadBinaryValue(file, entry.trajectory_id);
	ReadBinaryValue(file, entry.completion_wall_time_sec);
	ReadBinaryValue(file, entry.lifetime_unbinding_sec);
	ReadBinaryValue(file, entry.r_capture_rsun);
	ReadBinaryValue(file, entry.E_capture_eV);
	ReadBinaryValue(file, entry.dE_capture_eV);
}

template<std::size_t N>
bool IsValidHistogram(const std::array<double, N>& values)
{
	return std::all_of(values.begin(), values.end(), [](double value)
	{
		return std::isfinite(value) && value >= 0.0;
	});
}

bool IsValidRankState(const SnapshotRankState& state)
{
	if(state.snapshot_index < 0 || state.rank < 0)
		return false;
	if((state.done != 0 && state.done != 1)
	   || (state.trajectory_in_progress != 0 && state.trajectory_in_progress != 1)
	   || (state.current_trajectory_captured != 0 && state.current_trajectory_captured != 1))
		return false;
	if(state.local_captured > state.local_classified
	   || state.local_classified > state.local_total
	   || state.local_numerical_failures > state.local_total
	   || state.bincount_captured_samples > state.local_total
	   || state.bincount_not_captured_samples > state.local_total
	   || state.bincount_captured_samples > state.local_total - state.bincount_not_captured_samples)
		return false;
	if(!std::isfinite(state.rank_elapsed_wall_sec) || state.rank_elapsed_wall_sec < 0.0
	   || !std::isfinite(state.current_trajectory_wall_sec) || state.current_trajectory_wall_sec < 0.0
	   || !std::isfinite(state.current_trajectory_simulated_elapsed_sec)
	   || state.current_trajectory_simulated_elapsed_sec < 0.0)
		return false;
	if(state.done && state.trajectory_in_progress)
		return false;
	for(const auto& entry : state.new_evaporation_events)
	{
		if(entry.completion_wall_time_sec > state.rank_elapsed_wall_sec)
			return false;
	}
	if(state.trajectory_in_progress && state.current_trajectory_id == 0)
		return false;
	if(!state.trajectory_in_progress
	   && (state.current_trajectory_id != 0
	       || state.current_trajectory_captured != 0
	       || state.current_trajectory_wall_sec != 0.0
	       || state.current_trajectory_simulated_elapsed_sec != 0.0
	       || state.current_trajectory_scatterings != 0))
		return false;

	return IsValidHistogram(state.current_trajectory_dt_hist)
	    && IsValidHistogram(state.current_trajectory_v2dt_hist)
	    && IsValidHistogram(state.captured_dt_hist)
	    && IsValidHistogram(state.captured_v2dt_hist)
	    && IsValidHistogram(state.captured_dt_sq_hist)
	    && IsValidHistogram(state.captured_v2dt_sq_hist)
	    && IsValidHistogram(state.not_captured_dt_hist)
	    && IsValidHistogram(state.not_captured_v2dt_hist)
	    && IsValidHistogram(state.not_captured_dt_sq_hist)
	    && IsValidHistogram(state.not_captured_v2dt_sq_hist);
}

bool IsValidEvaporationEvents(const std::vector<SnapshotEvaporationProgressEntry>& entries)
{
	double previous_completion_time = -1.0;
	for(const auto& entry : entries)
	{
		if(entry.trajectory_id == 0
		   || !std::isfinite(entry.completion_wall_time_sec)
		   || entry.completion_wall_time_sec < 0.0
		   || !std::isfinite(entry.lifetime_unbinding_sec)
		   || entry.lifetime_unbinding_sec < 0.0
		   || !std::isfinite(entry.r_capture_rsun)
		   || entry.r_capture_rsun < -1.0
		   || !std::isfinite(entry.E_capture_eV)
		   || !std::isfinite(entry.dE_capture_eV)
		   || entry.completion_wall_time_sec < previous_completion_time)
			return false;
		previous_completion_time = entry.completion_wall_time_sec;
	}
	return true;
}

CompactEvaporationEvent MakeLogEvent(int rank, const SnapshotEvaporationProgressEntry& entry)
{
	CompactEvaporationEvent event;
	event.rank = rank;
	event.trajectory_id = entry.trajectory_id;
	event.completion_wall_time_sec = entry.completion_wall_time_sec;
	event.lifetime_unbinding = entry.lifetime_unbinding_sec;
	event.r_capture_rsun = entry.r_capture_rsun;
	event.E_capture_eV = entry.E_capture_eV;
	event.dE_capture_eV = entry.dE_capture_eV;
	return event;
}

double SnapshotBinError(double sum, double sum_sq, double count)
{
	if(count <= 1.0)
		return 0.0;

	const double mean = sum / count;
	double variance = sum_sq / count - mean * mean;
	if(variance < 0.0)
		variance = 0.0;
	return std::sqrt(count * variance);
}

struct SnapshotRateEstimate
{
	double rate = 0.0;
	double standard_error = 0.0;
	double ci_lower = 0.0;
	double ci_upper = 0.0;
};

SnapshotRateEstimate EstimateSnapshotRate(uint64_t trials, uint64_t successes)
{
	SnapshotRateEstimate estimate;
	if(trials == 0)
		return estimate;

	const double N = static_cast<double>(trials);
	estimate.rate = static_cast<double>(std::min(trials, successes)) / N;
	estimate.standard_error = std::sqrt(estimate.rate * (1.0 - estimate.rate) / N);
	const double z = 1.96;
	const double denominator = 1.0 + z * z / N;
	const double center = estimate.rate + z * z / (2.0 * N);
	const double spread = z * std::sqrt(estimate.rate * (1.0 - estimate.rate) / N + z * z / (4.0 * N * N));
	estimate.ci_lower = std::max(0.0, (center - spread) / denominator);
	estimate.ci_upper = std::min(1.0, (center + spread) / denominator);
	return estimate;
}

bool IsCompletedEvaporationEvent(const CompactEvaporationEvent& event)
{
	return std::isfinite(event.lifetime_unbinding) && event.lifetime_unbinding >= 0.0;
}

bool EvaporationEventOrder(const CompactEvaporationEvent& lhs, const CompactEvaporationEvent& rhs)
{
	if(lhs.lifetime_unbinding != rhs.lifetime_unbinding)
		return lhs.lifetime_unbinding < rhs.lifetime_unbinding;
	if(lhs.rank != rhs.rank)
		return lhs.rank < rhs.rank;
	return lhs.trajectory_id < rhs.trajectory_id;
}

void WriteEvaporationLogEvent(std::ostream& file, const CompactEvaporationEvent& event)
{
	file << event.rank << "\t" << event.trajectory_id << "\t" << std::scientific << std::setprecision(10)
	     << event.lifetime_unbinding << "\t" << event.r_capture_rsun
	     << "\t" << event.E_capture_eV << "\t" << event.dE_capture_eV << "\n";
}

void WriteEvaporationLogFileHeader(std::ofstream& file, double mass_gev, double sigma_cm2)
{
	file << "# DaMaSCUS-SUN snapshot evaporation times\n";
	file << "# format_version = 4\n";
	file << "# DIAGNOSTIC_ONLY = 1\n";
	file << "# NOT_FOR_FINAL_SURVIVAL_ANALYSIS = 1\n";
	file << "# completion_time_selected = 1\n";
	file << "# DM_mass_GeV = " << std::scientific << std::setprecision(6) << mass_gev << "\n";
	file << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(6) << sigma_cm2 << "\n";
	file << "# sorted_by = lifetime_unbinding_sec rank trajectory_id\n";
	file << "# rank trajectory_id lifetime_unbinding_sec r_capture_Rsun E_capture_eV dE_capture_eV\n";
}

void WriteEvaporationLogEvents(std::ostream& file, const std::vector<CompactEvaporationEvent>& events)
{
	std::vector<CompactEvaporationEvent> sorted_events = events;
	std::sort(sorted_events.begin(), sorted_events.end(), EvaporationEventOrder);
	for(const auto& event : sorted_events)
	{
		if(IsCompletedEvaporationEvent(event))
			WriteEvaporationLogEvent(file, event);
	}
}

bool WriteTextFileAtomically(const std::string& path, int unique_tag, const std::function<void(std::ofstream&)>& writer)
{
	const std::string tmp_path = path + ".tmp." + std::to_string(getpid()) + "." + std::to_string(unique_tag);
	std::ofstream file(tmp_path, std::ios::trunc);
	if(!file.is_open())
		return false;

	writer(file);
	file.close();

	if(!file.good())
	{
		std::remove(tmp_path.c_str());
		return false;
	}

	if(std::rename(tmp_path.c_str(), path.c_str()) != 0)
	{
		std::remove(tmp_path.c_str());
		return false;
	}

	return true;
}

std::string JoinInts(const std::vector<int>& values)
{
	if(values.empty())
		return "";

	std::ostringstream stream;
	for(size_t i = 0; i < values.size(); i++)
	{
		if(i != 0)
			stream << ",";
		stream << values[i];
	}
	return stream.str();
}

bool SnapshotTextFileIsMerged(const std::string& path, uint64_t expected_run_id)
{
	std::ifstream file(path);
	if(!file.is_open())
		return false;

	bool merged = false;
	bool run_id_matches = false;
	std::string line;
	while(std::getline(file, line))
	{
		if(line == "# snapshot_status = merged")
			merged = true;
		const std::string run_id_prefix = "# snapshot_run_id = ";
		if(line.compare(0, run_id_prefix.size(), run_id_prefix) == 0)
		{
			std::istringstream value_stream(line.substr(run_id_prefix.size()));
			uint64_t run_id = 0;
			if((value_stream >> run_id) && value_stream.eof())
				run_id_matches = run_id == expected_run_id;
		}
		if(!line.empty() && line[0] != '#')
			break;
	}
	return merged && run_id_matches;
}

const char* SnapshotRankActivityLabel(const SnapshotReportState::RankProgress& progress)
{
	if(progress.done)
		return "done";
	if(!progress.trajectory_in_progress)
		return "idle_or_waiting";
	return progress.current_trajectory_captured ? "running_captured" : "running_uncaptured";
}

void AccumulateSnapshotReportState(SnapshotReportState& report, const SnapshotRankState& state)
{
	SnapshotReportState::RankProgress progress;
	progress.rank = state.rank;
	progress.done = state.done;
	progress.trajectory_in_progress = state.trajectory_in_progress;
	progress.current_trajectory_captured = state.current_trajectory_captured;
	progress.current_trajectory_id = state.current_trajectory_id;
	progress.current_trajectory_scatterings = state.current_trajectory_scatterings;
	progress.rank_elapsed_wall_sec = state.rank_elapsed_wall_sec;
	progress.current_trajectory_wall_sec = state.current_trajectory_wall_sec;
	progress.current_trajectory_simulated_elapsed_sec = state.current_trajectory_simulated_elapsed_sec;
	report.rank_progress.push_back(progress);

	report.total_trajectories += state.local_total;
	report.captured_particles += state.local_captured;
	report.classified_trajectories += state.local_classified;
	report.numerical_failures += state.local_numerical_failures;
	report.snapshot_bincount_captured_samples += state.bincount_captured_samples;
	report.snapshot_bincount_not_captured_samples += state.bincount_not_captured_samples;

	const double lower_wall_time = (report.snapshot_time_label > report.snapshot_interval_seconds) ? (report.snapshot_time_label - report.snapshot_interval_seconds) : 0.0;
	const double upper_wall_time = static_cast<double>(report.snapshot_time_label);
	const bool first_final_snapshot = state.done
	                               && state.snapshot_index < report.snapshot_index
	                               && (state.snapshot_index + 1 == report.snapshot_index
	                                   || state.rank_elapsed_wall_sec > lower_wall_time);
	for(const auto& entry : state.new_evaporation_events)
	{
		// Checkpoint entries are already assigned by the recorder cursor. A late
		// boundary event must remain in that publication instead of being dropped.
		const bool belongs_to_snapshot = entry.completion_wall_time_sec <= upper_wall_time
		                               && (!state.done
		                                   || first_final_snapshot
		                                   || entry.completion_wall_time_sec > lower_wall_time);
		if(belongs_to_snapshot)
		{
			SnapshotRankedEvaporationEntry ranked_entry;
			ranked_entry.rank = state.rank;
			ranked_entry.entry = entry;
			report.new_evaporation_events.push_back(ranked_entry);
		}
	}

	for(std::size_t bin = 0; bin < TOTAL_BINS; bin++)
	{
		report.captured_dt_hist[bin] += state.captured_dt_hist[bin];
		report.captured_v2dt_hist[bin] += state.captured_v2dt_hist[bin];
		report.captured_dt_sq_hist[bin] += state.captured_dt_sq_hist[bin];
		report.captured_v2dt_sq_hist[bin] += state.captured_v2dt_sq_hist[bin];
		report.not_captured_dt_hist[bin] += state.not_captured_dt_hist[bin];
		report.not_captured_v2dt_hist[bin] += state.not_captured_v2dt_hist[bin];
		report.not_captured_dt_sq_hist[bin] += state.not_captured_dt_sq_hist[bin];
		report.not_captured_v2dt_sq_hist[bin] += state.not_captured_v2dt_sq_hist[bin];
	}

	// In-progress trajectories have not passed the final validity/domain
	// classification yet. Keep their progress in the rank-status section, but
	// do not publish provisional residence time into the physical histogram.
}

SnapshotMergeResult LoadSnapshotReportState(
	const std::string& rank_snapshot_dir,
	int snapshot_index,
	double interval_seconds,
	int mpi_processes,
	uint64_t run_id,
	SnapshotMergeCache& cache)
{
	const size_t rank_count = static_cast<size_t>(std::max(0, mpi_processes));
	if(!cache.initialized
	   || cache.report.snapshot_index != snapshot_index
	   || cache.rank_accumulated.size() != rank_count)
	{
		cache.report = SnapshotReportState();
		cache.report.snapshot_index = snapshot_index;
		cache.report.snapshot_time_label = SnapshotTimeLabelSeconds(snapshot_index, interval_seconds);
		cache.report.snapshot_interval_seconds = interval_seconds;
		cache.rank_accumulated.assign(rank_count, 0);
		cache.initialized = true;
	}

	SnapshotReportState& report = cache.report;
	SnapshotMergeResult result;
	bool accumulated_any = false;
	for(int rank = 0; rank < mpi_processes; rank++)
	{
		if(cache.rank_accumulated[static_cast<size_t>(rank)] != 0)
		{
			result.ready_ranks.push_back(rank);
			continue;
		}

		SnapshotRankState state;
		const std::string checkpoint_path = SnapshotRankCheckpointPath(rank_snapshot_dir, rank, snapshot_index, interval_seconds);
		if(ReadSnapshotRankState(checkpoint_path, run_id, state)
		   && state.rank == rank
		   && state.snapshot_index == snapshot_index
		   && !state.done)
		{
			AccumulateSnapshotReportState(report, state);
			cache.rank_accumulated[static_cast<size_t>(rank)] = 1;
			accumulated_any = true;
			result.ready_ranks.push_back(rank);
			continue;
		}

		const std::string final_path = SnapshotRankFinalPath(rank_snapshot_dir, rank);
		if(ReadSnapshotRankState(final_path, run_id, state)
		   && state.rank == rank
		   && state.done
		   && state.snapshot_index < snapshot_index
		   && state.rank_elapsed_wall_sec <= report.snapshot_time_label)
		{
			AccumulateSnapshotReportState(report, state);
			cache.rank_accumulated[static_cast<size_t>(rank)] = 1;
			accumulated_any = true;
			result.ready_ranks.push_back(rank);
			continue;
		}

		// A rank missing now can still show up later through its final state, so
		// only already accumulated ranks may skip the file reads.
		result.missing_ranks.push_back(rank);
	}

	// Ranks are folded in whenever they first become readable, so restore the
	// ascending rank order the report file has always been written in.
	if(accumulated_any)
	{
		std::sort(report.rank_progress.begin(), report.rank_progress.end(),
			[](const SnapshotReportState::RankProgress& lhs, const SnapshotReportState::RankProgress& rhs)
			{
				return lhs.rank < rhs.rank;
			});
	}

	if(result.ready_ranks.empty())
		result.status = SnapshotMergeStatus::NoRanksReady;
	else if(result.missing_ranks.empty())
		result.status = SnapshotMergeStatus::Merged;
	else
		result.status = SnapshotMergeStatus::Partial;

	return result;
}

void WriteReportHeader(
	std::ofstream& file,
	double mass_gev,
	double sigma_cm2,
	uint64_t total_trajectories,
	uint64_t captured_particles,
	uint64_t classified_trajectories,
	uint64_t numerical_failures,
	long long snapshot_time_label,
	double snapshot_interval_seconds)
{
	file << "# snapshot_target_wall_time_s = " << snapshot_time_label << "\n";
	file << "# snapshot_interval_s = " << std::fixed << std::setprecision(3) << snapshot_interval_seconds << "\n";
	file << "# DM_mass_GeV = " << std::scientific << std::setprecision(6) << mass_gev << "\n";
	file << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(6) << sigma_cm2 << "\n";
	file << "# total_trajectories = " << total_trajectories << "\n";
	file << "# captured_particles = " << captured_particles << "\n";
	file << "# bincount_integration = " << BincountIntegrationScheme() << "\n";
	file << "# bincount_dense_position_tolerance_km = "
	     << std::scientific << std::setprecision(6)
	     << BincountDensePositionToleranceKm() << "\n";
	file << "# valid_trajectories = " << classified_trajectories << "\n";
	file << "# numerical_failures = " << numerical_failures << "\n";
	file << "# unresolved_not_captured_trajectories = "
	     << (total_trajectories - classified_trajectories) << "\n";

	const SnapshotRateEstimate raw = EstimateSnapshotRate(total_trajectories, captured_particles);
	const SnapshotRateEstimate valid = EstimateSnapshotRate(classified_trajectories, captured_particles);
	file << "# capture_rate = " << std::fixed << std::setprecision(8) << raw.rate << "\n";
	file << "# capture_rate_raw = " << std::fixed << std::setprecision(8) << raw.rate << "\n";
	file << "# capture_rate_err = " << std::fixed << std::setprecision(8) << raw.standard_error << "\n";
	file << "# capture_rate_raw_err = " << std::fixed << std::setprecision(8) << raw.standard_error << "\n";
	file << "# capture_rate_CI_95_lower = " << std::fixed << std::setprecision(8) << raw.ci_lower << "\n";
	file << "# capture_rate_CI_95_upper = " << std::fixed << std::setprecision(8) << raw.ci_upper << "\n";
	file << "# capture_rate_raw_CI_95_lower = " << std::fixed << std::setprecision(8) << raw.ci_lower << "\n";
	file << "# capture_rate_raw_CI_95_upper = " << std::fixed << std::setprecision(8) << raw.ci_upper << "\n";
	file << "# capture_rate_valid = " << std::fixed << std::setprecision(8) << valid.rate << "\n";
	file << "# capture_rate_valid_err = " << std::fixed << std::setprecision(8) << valid.standard_error << "\n";
	file << "# capture_rate_valid_CI_95_lower = " << std::fixed << std::setprecision(8) << valid.ci_lower << "\n";
	file << "# capture_rate_valid_CI_95_upper = " << std::fixed << std::setprecision(8) << valid.ci_upper << "\n";
}

bool WriteSnapshotEvaporationTimeFile(
	const std::string& snapshot_root,
	int snapshot_index,
	double interval_seconds,
	uint64_t run_id,
	double mass_gev,
	double sigma_cm2,
	const std::vector<SnapshotRankedEvaporationEntry>& entries,
	const SnapshotMergeResult& merge_result)
{
	std::vector<CompactEvaporationEvent> events;
	events.reserve(entries.size());
	for(const auto& entry : entries)
		events.push_back(MakeLogEvent(entry.rank, entry.entry));

	const std::string status = (merge_result.status == SnapshotMergeStatus::Merged) ? "merged" : "partial";
	const std::string path = SnapshotEvaporationTimeFilePath(snapshot_root, snapshot_index, interval_seconds);
	return WriteTextFileAtomically(path, snapshot_index, [&](std::ofstream& file)
	{
		file << "# snapshot_status = " << status << "\n";
		file << "# snapshot_run_id = " << run_id << "\n";
		file << "# snapshot_target_wall_time_s = " << SnapshotTimeLabelSeconds(snapshot_index, interval_seconds) << "\n";
		file << "# snapshot_interval_s = " << std::fixed << std::setprecision(3) << interval_seconds << "\n";
		file << "# ready_ranks = " << JoinInts(merge_result.ready_ranks) << "\n";
		file << "# missing_ranks = " << JoinInts(merge_result.missing_ranks) << "\n";
		file << "# completed_evaporation_events_assigned_to_snapshot = 1\n";
		file << "# event_assignment = first published checkpoint after completion\n";
		WriteEvaporationLogFileHeader(file, mass_gev, sigma_cm2);
		WriteEvaporationLogEvents(file, events);
	});
}

bool WriteSnapshotReportFile(
	const std::string& snapshot_root,
	int snapshot_index,
	double interval_seconds,
	uint64_t run_id,
	double mass_gev,
	double sigma_cm2,
	const SnapshotReportState& report,
	const SnapshotMergeResult& merge_result)
{
	if(merge_result.status != SnapshotMergeStatus::Merged
	   && SnapshotTextFileIsMerged(SnapshotTextFilePath(snapshot_root, snapshot_index, interval_seconds), run_id))
		return true;

	if(!WriteSnapshotEvaporationTimeFile(snapshot_root, snapshot_index, interval_seconds, run_id, mass_gev, sigma_cm2, report.new_evaporation_events, merge_result))
		return false;

	const std::string status = (merge_result.status == SnapshotMergeStatus::Merged) ? "merged" : "partial";
	return WriteTextFileAtomically(SnapshotTextFilePath(snapshot_root, snapshot_index, interval_seconds), snapshot_index, [&](std::ofstream& file)
	{
		file << "# snapshot_status = " << status << "\n";
		file << "# snapshot_run_id = " << run_id << "\n";
		file << "# Cumulative snapshot report\n";
		file << "# ready_ranks = " << JoinInts(merge_result.ready_ranks) << "\n";
		file << "# missing_ranks = " << JoinInts(merge_result.missing_ranks) << "\n";
		if(merge_result.status == SnapshotMergeStatus::Partial)
			file << "# WARNING = partial snapshot; not all ranks reached this checkpoint\n";
		WriteReportHeader(
			file,
			mass_gev,
			sigma_cm2,
			report.total_trajectories,
			report.captured_particles,
			report.classified_trajectories,
			report.numerical_failures,
			report.snapshot_time_label,
			report.snapshot_interval_seconds);

		file << "#\n";
		file << "# [MPI rank status]\n";
		file << "# rank  state  trajectory_id  trajectory_wall_s  simulated_elapsed_s  scatterings  observed_at_wall_s\n";
		for(const auto& progress : report.rank_progress)
		{
			file << "# " << progress.rank << "\t" << SnapshotRankActivityLabel(progress)
			     << "\t" << progress.current_trajectory_id
			     << "\t" << std::scientific << std::setprecision(10) << progress.current_trajectory_wall_sec
			     << "\t" << progress.current_trajectory_simulated_elapsed_sec
			     << "\t" << progress.current_trajectory_scatterings
			     << "\t" << progress.rank_elapsed_wall_sec << "\n";
		}

		const double snapshot_captured_samples = static_cast<double>(report.snapshot_bincount_captured_samples);
		const double snapshot_not_captured_samples = static_cast<double>(report.snapshot_bincount_not_captured_samples);
		file << "#\n";
		file << "# [Bincount histogram]\n";
		file << "# base_grid_bins = " << NUM_BINS << "\n";
		file << "# exterior_log_bins = " << EXTERIOR_LOG_BINS << "\n";
		file << "# total_radial_bins = " << TOTAL_BINS << "\n";
		file << "# radial_bin_width_Rsun = " << std::scientific << std::setprecision(10)
		     << BIN_WIDTH_KM / R_SUN_KM << "\n";
		file << "# exterior_grid = logarithmic\n";
		file << "# radial_domain_max_AU = " << RADIAL_DOMAIN_MAX_AU << "\n";
		file << "# radial_extent_Rsun = " << RADIAL_DOMAIN_MAX_RSUN << "\n";
		file << "# in_progress_bincount_included = 0\n";
		file << "# bin_index  r_lower_Rsun  r_upper_Rsun  cap_dt[s]  cap_v2dt[km2/s]  cap_err_dt[s]  cap_err_v2dt[km2/s]  not_cap_dt[s]  not_cap_v2dt[km2/s]  not_cap_err_dt[s]  not_cap_err_v2dt[km2/s]\n";
		for(std::size_t bin = 0; bin < report.captured_dt_hist.size(); bin++)
		{
			const double cap_err_dt = SnapshotBinError(report.captured_dt_hist[bin], report.captured_dt_sq_hist[bin], snapshot_captured_samples);
			const double cap_err_v2dt = SnapshotBinError(report.captured_v2dt_hist[bin], report.captured_v2dt_sq_hist[bin], snapshot_captured_samples);
			const double not_cap_err_dt = SnapshotBinError(report.not_captured_dt_hist[bin], report.not_captured_dt_sq_hist[bin], snapshot_not_captured_samples);
			const double not_cap_err_v2dt = SnapshotBinError(report.not_captured_v2dt_hist[bin], report.not_captured_v2dt_sq_hist[bin], snapshot_not_captured_samples);

			file << bin << "\t" << std::scientific << std::setprecision(10)
			     << BincountBinLowerKm(bin) / R_SUN_KM << "\t"
			     << BincountBinUpperKm(bin) / R_SUN_KM << "\t"
			     << report.captured_dt_hist[bin] << "\t" << report.captured_v2dt_hist[bin]
			     << "\t" << cap_err_dt << "\t" << cap_err_v2dt
			     << "\t" << report.not_captured_dt_hist[bin] << "\t" << report.not_captured_v2dt_hist[bin]
			     << "\t" << not_cap_err_dt << "\t" << not_cap_err_v2dt << "\n";
		}
	});
}
}

bool IsValidSnapshotIntervalSeconds(double interval_seconds)
{
	return std::isfinite(interval_seconds)
	    && interval_seconds >= 1.0
	    && interval_seconds <= static_cast<double>(std::numeric_limits<int>::max())
	    && std::floor(interval_seconds) == interval_seconds;
}

long long SnapshotTimeLabelSeconds(int snapshot_index, double interval_seconds)
{
	if(snapshot_index <= 0)
		throw std::invalid_argument("snapshot index must be positive");
	if(!IsValidSnapshotIntervalSeconds(interval_seconds))
		throw std::invalid_argument("snapshot interval must be a positive integer number of seconds");

	const long long interval = static_cast<long long>(interval_seconds);
	if(static_cast<long long>(snapshot_index) > std::numeric_limits<long long>::max() / interval)
		throw std::overflow_error("snapshot time label exceeds long long range");
	return static_cast<long long>(snapshot_index) * interval;
}

std::string SnapshotTextFilePath(const std::string& snapshot_root, int snapshot_index, double interval_seconds)
{
	return snapshot_root + "snapshot_" + std::to_string(SnapshotTimeLabelSeconds(snapshot_index, interval_seconds)) + "s.txt";
}

std::string SnapshotEvaporationTimeFilePath(const std::string& snapshot_root, int snapshot_index, double interval_seconds)
{
	return snapshot_root + "snapshot_" + std::to_string(SnapshotTimeLabelSeconds(snapshot_index, interval_seconds)) + "s_evaporation_times.txt";
}

std::string SnapshotRankCheckpointPath(const std::string& rank_snapshot_dir, int rank, int snapshot_index, double interval_seconds)
{
	return rank_snapshot_dir + "snapshot_" + std::to_string(SnapshotTimeLabelSeconds(snapshot_index, interval_seconds)) + "s_rank" + std::to_string(rank) + ".bin";
}

std::string SnapshotRankFinalPath(const std::string& rank_snapshot_dir, int rank)
{
	return rank_snapshot_dir + "rank" + std::to_string(rank) + "_final.bin";
}

SnapshotEvaporationProgressEntry MakeSnapshotEvaporationProgressEntry(const CompactEvaporationEvent& event)
{
	SnapshotEvaporationProgressEntry entry;
	entry.trajectory_id = static_cast<uint64_t>(event.trajectory_id);
	entry.completion_wall_time_sec = event.completion_wall_time_sec;
	entry.lifetime_unbinding_sec = event.lifetime_unbinding;
	entry.r_capture_rsun = event.r_capture_rsun;
	entry.E_capture_eV = event.E_capture_eV;
	entry.dE_capture_eV = event.dE_capture_eV;
	return entry;
}

size_t SnapshotEvaporationEventsPerCheckpoint()
{
	return static_cast<size_t>(SNAPSHOT_EVAPORATION_EVENTS_PER_CHECKPOINT);
}

bool WriteSnapshotRankState(const std::string& path, const SnapshotRankState& state)
{
	if(!IsValidRankState(state)
	   || state.new_evaporation_events.size() > MAX_SNAPSHOT_EVAPORATION_EVENTS
	   || !IsValidEvaporationEvents(state.new_evaporation_events))
		return false;

	const std::string tmp_path = path + ".tmp." + std::to_string(getpid()) + "." + std::to_string(state.rank);
	std::ofstream file(tmp_path, std::ios::binary | std::ios::trunc);
	if(!file.is_open())
		return false;

	WriteBinaryValue(file, SNAPSHOT_RANK_STATE_MAGIC);
	WriteBinaryValue(file, SNAPSHOT_RANK_STATE_VERSION);
	WriteBinaryValue(file, SNAPSHOT_RANK_STATE_HEADER_BYTES);
	WriteBinaryValue(file, state.run_id);
	WriteBinaryValue(file, state.snapshot_index);
	WriteBinaryValue(file, state.rank);
	WriteBinaryValue(file, state.done);
	WriteBinaryValue(file, state.trajectory_in_progress);
	WriteBinaryValue(file, state.local_captured);
	WriteBinaryValue(file, state.local_total);
	WriteBinaryValue(file, state.local_classified);
	WriteBinaryValue(file, state.local_numerical_failures);
	WriteBinaryValue(file, state.bincount_captured_samples);
	WriteBinaryValue(file, state.bincount_not_captured_samples);
	WriteBinaryValue(file, state.current_trajectory_id);
	WriteBinaryValue(file, state.rank_elapsed_wall_sec);
	WriteBinaryValue(file, state.current_trajectory_wall_sec);
	WriteBinaryValue(file, state.current_trajectory_simulated_elapsed_sec);
	WriteBinaryValue(file, state.current_trajectory_scatterings);
	WriteBinaryValue(file, state.current_trajectory_captured);
	const uint64_t radial_bin_count = TOTAL_BINS;
	WriteBinaryValue(file, radial_bin_count);
	WriteBinaryArray(file, state.current_trajectory_dt_hist);
	WriteBinaryArray(file, state.current_trajectory_v2dt_hist);
	WriteBinaryArray(file, state.captured_dt_hist);
	WriteBinaryArray(file, state.captured_v2dt_hist);
	WriteBinaryArray(file, state.captured_dt_sq_hist);
	WriteBinaryArray(file, state.captured_v2dt_sq_hist);
	WriteBinaryArray(file, state.not_captured_dt_hist);
	WriteBinaryArray(file, state.not_captured_v2dt_hist);
	WriteBinaryArray(file, state.not_captured_dt_sq_hist);
	WriteBinaryArray(file, state.not_captured_v2dt_sq_hist);
	const uint64_t event_count = static_cast<uint64_t>(state.new_evaporation_events.size());
	WriteBinaryValue(file, event_count);
	for(const auto& entry : state.new_evaporation_events)
		WriteSnapshotEvaporationEntryBinary(file, entry);
	file.close();

	if(!file.good())
	{
		std::remove(tmp_path.c_str());
		return false;
	}
	if(std::rename(tmp_path.c_str(), path.c_str()) != 0)
	{
		std::remove(tmp_path.c_str());
		return false;
	}
	return true;
}

bool ReadSnapshotRankState(const std::string& path, uint64_t expected_run_id, SnapshotRankState& state)
{
	std::ifstream file(path, std::ios::binary);
	if(!file.is_open())
		return false;

	file.seekg(0, std::ios::end);
	const std::streamoff file_size_value = file.tellg();
	const uint64_t minimum_file_size =
	    SnapshotRankStateFixedBytes() + 10ULL * TOTAL_BINS * sizeof(double);
	if(file_size_value < 0 || static_cast<uint64_t>(file_size_value) < minimum_file_size)
		return false;
	const uint64_t file_size = static_cast<uint64_t>(file_size_value);
	file.seekg(0, std::ios::beg);

	uint64_t magic = 0;
	uint32_t version = 0;
	uint32_t header_bytes = 0;
	ReadBinaryValue(file, magic);
	ReadBinaryValue(file, version);
	ReadBinaryValue(file, header_bytes);
	if(!file
	   || magic != SNAPSHOT_RANK_STATE_MAGIC
	   || version != SNAPSHOT_RANK_STATE_VERSION
	   || header_bytes != SNAPSHOT_RANK_STATE_HEADER_BYTES)
		return false;

	state = SnapshotRankState();
	ReadBinaryValue(file, state.run_id);
	if(!file || state.run_id != expected_run_id)
		return false;
	ReadBinaryValue(file, state.snapshot_index);
	ReadBinaryValue(file, state.rank);
	ReadBinaryValue(file, state.done);
	ReadBinaryValue(file, state.trajectory_in_progress);
	ReadBinaryValue(file, state.local_captured);
	ReadBinaryValue(file, state.local_total);
	ReadBinaryValue(file, state.local_classified);
	ReadBinaryValue(file, state.local_numerical_failures);
	ReadBinaryValue(file, state.bincount_captured_samples);
	ReadBinaryValue(file, state.bincount_not_captured_samples);
	ReadBinaryValue(file, state.current_trajectory_id);
	ReadBinaryValue(file, state.rank_elapsed_wall_sec);
	ReadBinaryValue(file, state.current_trajectory_wall_sec);
	ReadBinaryValue(file, state.current_trajectory_simulated_elapsed_sec);
	ReadBinaryValue(file, state.current_trajectory_scatterings);
	ReadBinaryValue(file, state.current_trajectory_captured);
	uint64_t radial_bin_count = 0;
	ReadBinaryValue(file, radial_bin_count);
	if(!file || radial_bin_count != TOTAL_BINS)
		return false;
	ReadBinaryArray(file, state.current_trajectory_dt_hist);
	ReadBinaryArray(file, state.current_trajectory_v2dt_hist);
	ReadBinaryArray(file, state.captured_dt_hist);
	ReadBinaryArray(file, state.captured_v2dt_hist);
	ReadBinaryArray(file, state.captured_dt_sq_hist);
	ReadBinaryArray(file, state.captured_v2dt_sq_hist);
	ReadBinaryArray(file, state.not_captured_dt_hist);
	ReadBinaryArray(file, state.not_captured_v2dt_hist);
	ReadBinaryArray(file, state.not_captured_dt_sq_hist);
	ReadBinaryArray(file, state.not_captured_v2dt_sq_hist);

	uint64_t event_count = 0;
	ReadBinaryValue(file, event_count);
	if(!file || event_count > MAX_SNAPSHOT_EVAPORATION_EVENTS)
		return false;
	const uint64_t entry_bytes = SnapshotEvaporationEntryBytes();
	const uint64_t histogram_bytes =
	    10ULL * radial_bin_count * sizeof(double);
	if(event_count > (std::numeric_limits<uint64_t>::max()
	                  - SnapshotRankStateFixedBytes() - histogram_bytes) / entry_bytes)
		return false;
	const uint64_t expected_file_size =
	    SnapshotRankStateFixedBytes() + histogram_bytes + event_count * entry_bytes;
	if(file_size != expected_file_size)
		return false;

	state.new_evaporation_events.clear();
	try
	{
		state.new_evaporation_events.resize(static_cast<size_t>(event_count));
	}
	catch(const std::exception&)
	{
		state.new_evaporation_events.clear();
		return false;
	}
	for(size_t i = 0; i < state.new_evaporation_events.size(); i++)
		ReadSnapshotEvaporationEntryBinary(file, state.new_evaporation_events[i]);

	if(!file)
		return false;
	return IsValidRankState(state) && IsValidEvaporationEvents(state.new_evaporation_events);
}

bool CleanupSnapshotCheckpoints(const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes)
{
	bool cleanup_succeeded = true;
	for(int rank = 0; rank < mpi_processes; rank++)
	{
		const std::string path = SnapshotRankCheckpointPath(rank_snapshot_dir, rank, snapshot_index, interval_seconds);
		errno = 0;
		if(std::remove(path.c_str()) != 0 && errno != ENOENT)
		{
			cleanup_succeeded = false;
			std::cerr << "Warning in CleanupSnapshotCheckpoints(): failed to remove " << path << std::endl;
		}
	}
	return cleanup_succeeded;
}

bool CleanupFinalSnapshotStates(const std::string& rank_snapshot_dir, int mpi_processes)
{
	bool cleanup_succeeded = true;
	for(int rank = 0; rank < mpi_processes; rank++)
	{
		const std::string path = SnapshotRankFinalPath(rank_snapshot_dir, rank);
		errno = 0;
		if(std::remove(path.c_str()) != 0 && errno != ENOENT)
		{
			cleanup_succeeded = false;
			std::cerr << "Warning in CleanupFinalSnapshotStates(): failed to remove " << path << std::endl;
		}
	}

	errno = 0;
	if(rmdir(rank_snapshot_dir.c_str()) != 0 && errno != ENOENT)
	{
		cleanup_succeeded = false;
		std::cerr << "Warning in CleanupFinalSnapshotStates(): failed to remove directory "
		          << rank_snapshot_dir << std::endl;
	}
	return cleanup_succeeded;
}

SnapshotMergeResult TryWriteSnapshot(
	const std::string& snapshot_root,
	const std::string& rank_snapshot_dir,
	int snapshot_index,
	double interval_seconds,
	int mpi_processes,
	uint64_t run_id,
	double mass_gev,
	double sigma_cm2,
	bool allow_partial)
{
	SnapshotMergeCache cache;
	return TryWriteSnapshotCached(
		snapshot_root,
		rank_snapshot_dir,
		snapshot_index,
		interval_seconds,
		mpi_processes,
		run_id,
		mass_gev,
		sigma_cm2,
		cache,
		allow_partial);
}

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
	bool allow_partial)
{
	SnapshotMergeResult merged_result;
	merged_result.status = SnapshotMergeStatus::Merged;
	if(SnapshotTextFileIsMerged(SnapshotTextFilePath(snapshot_root, snapshot_index, interval_seconds), run_id)
	   && SnapshotTextFileIsMerged(SnapshotEvaporationTimeFilePath(snapshot_root, snapshot_index, interval_seconds), run_id))
	{
		merged_result.cleanup_succeeded = CleanupSnapshotCheckpoints(
			rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes);
		return merged_result;
	}

	SnapshotMergeResult result = LoadSnapshotReportState(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes, run_id, cache);
	if(result.status == SnapshotMergeStatus::NoRanksReady)
		return result;
	if(result.status == SnapshotMergeStatus::Partial && !allow_partial)
		return result;

	if(!WriteSnapshotReportFile(snapshot_root, snapshot_index, interval_seconds, run_id, mass_gev, sigma_cm2, cache.report, result))
	{
		result.status = SnapshotMergeStatus::NoRanksReady;
		return result;
	}

	if(result.status == SnapshotMergeStatus::Merged)
		result.cleanup_succeeded = CleanupSnapshotCheckpoints(
			rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes);
	return result;
}

bool WriteMissedSnapshotMarker(
	const std::string& snapshot_root,
	int snapshot_index,
	double interval_seconds,
	uint64_t run_id,
	double actual_write_wall_time_sec)
{
	const std::string path = SnapshotTextFilePath(snapshot_root, snapshot_index, interval_seconds);
	if(SnapshotTextFileIsMerged(path, run_id))
		return true;

	return WriteTextFileAtomically(path, snapshot_index, [&](std::ofstream& file)
	{
		file << "# snapshot_status = missed_writer_delay\n";
		file << "# snapshot_run_id = " << run_id << "\n";
		file << "# snapshot_target_wall_time_s = " << SnapshotTimeLabelSeconds(snapshot_index, interval_seconds) << "\n";
		file << "# snapshot_interval_s = " << std::fixed << std::setprecision(3) << interval_seconds << "\n";
		file << "# actual_write_wall_time_s = " << std::fixed << std::setprecision(3) << actual_write_wall_time_sec << "\n";
		file << "# note = heartbeat woke too late; no state rewind is available\n";
	});
}

}	// namespace DaMaSCUS_SUN

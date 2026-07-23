#include "Data_Generation.hpp"

#include <dirent.h>
#include <errno.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <mpi.h>
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <type_traits>
#include <unistd.h>
#include <utility>
#include <vector>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/Astronomy.hpp"
#include "MPI_Tail_Work.hpp"
#include "Snapshot_Heartbeat.hpp"
#include "Snapshot_Shared_State.hpp"
#include "version.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

namespace
{
bool Has_Positive_Evaporation_Time(double t_evap)
{
	return std::isfinite(t_evap) && t_evap > 0.0;
}

void MPI_Trace_Point(int mpi_rank, const std::string& label)
{
	if(std::getenv("DAMASCUS_MPI_TRACE") == nullptr)
		return;
	std::cerr << "[mpi-trace rank " << mpi_rank << "] " << label << std::endl;
}

constexpr unsigned long int CAPTURE_MODE_MPI_SYNC_INTERVAL = 32UL;
constexpr unsigned long int NORMAL_MODE_MPI_SYNC_INTERVAL_FALLBACK = 128UL;
constexpr double NUMERICAL_FAILURE_WARNING_FRACTION = 1.0e-4;
constexpr double NUMERICAL_FAILURE_ABORT_FRACTION = 1.0e-2;

unsigned long int Normal_Mode_MPI_Sync_Interval(double sigma_cm2)
{
	if(!std::isfinite(sigma_cm2) || sigma_cm2 <= 0.0)
		return NORMAL_MODE_MPI_SYNC_INTERVAL_FALLBACK;
	if(sigma_cm2 >= 1.0e-35)
		return 64UL;
	if(sigma_cm2 >= 1.0e-36)
		return 128UL;
	if(sigma_cm2 >= 1.0e-37)
		return 1024UL;
	if(sigma_cm2 >= 1.0e-38)
		return 8192UL;
	if(sigma_cm2 >= 1.0e-39)
		return 65536UL;
	return 1048576UL;
}

bool Is_Completed_Evaporation_Record(const EvaporationRecord& rec)
{
	return rec.survival_valid
	    && rec.event_observed
	    && std::isfinite(rec.lifetime_unbinding)
	    && rec.lifetime_unbinding >= 0.0;
}

bool Is_Completed_Evaporation_Event(const CompactEvaporationEvent& event)
{
	return std::isfinite(event.lifetime_unbinding) && event.lifetime_unbinding >= 0.0;
}

int TerminationReason_Index(TrajectoryTerminationReason reason)
{
	int index = static_cast<int>(reason);
	if(index < 0 || index >= TRAJECTORY_TERMINATION_REASON_COUNT)
		return static_cast<int>(TrajectoryTerminationReason::Unknown);
	return index;
}

bool Completed_Outward_Escape(TrajectoryTerminationReason reason)
{
	return reason == TrajectoryTerminationReason::OutwardEscape;
}

bool Is_Computational_Truncation(TrajectoryTerminationReason reason)
{
	return reason == TrajectoryTerminationReason::WallTimeLimit
	    || reason == TrajectoryTerminationReason::MaxFreeSteps
	    || reason == TrajectoryTerminationReason::MaxScatterings;
}

bool Is_Numerical_Termination(TrajectoryTerminationReason reason)
{
	return reason == TrajectoryTerminationReason::NumericalFailure
	    || reason == TrajectoryTerminationReason::NonFiniteState
	    || reason == TrajectoryTerminationReason::SpeedLimit
	    || reason == TrajectoryTerminationReason::EnergyDriftEscape
	    || reason == TrajectoryTerminationReason::Unknown;
}

bool Build_Evaporation_Record(const TrajectoryBincount& bincount, int mpi_rank, unsigned long int trajectory_id, double completion_wall_time_sec, EvaporationRecord& rec)
{
	if(!bincount.is_captured || !std::isfinite(bincount.t_capture))
		return false;

	const bool numerically_invalid_escape =
	    bincount.numerically_invalid_escape || bincount.termination_reason == TrajectoryTerminationReason::EnergyDriftEscape;
	const bool survival_valid = bincount.survival_valid && !numerically_invalid_escape;
	const bool event_observed = survival_valid && bincount.event_observed;
	const double t_termination = std::isfinite(bincount.t_termination) ? bincount.t_termination : bincount.t_capture;
	const double lifetime_unbinding = event_observed ? (bincount.t_final_unbinding_scatter - bincount.t_capture) : -1.0;
	const double lifetime_boundary = bincount.boundary_escape_observed ? (bincount.t_boundary_escape - bincount.t_capture) : -1.0;
	double observed_lifetime = event_observed ? lifetime_unbinding : (t_termination - bincount.t_capture);
	if(!survival_valid)
		observed_lifetime = std::numeric_limits<double>::quiet_NaN();
	else if(!std::isfinite(observed_lifetime) || observed_lifetime < 0.0)
		observed_lifetime = 0.0;

	rec.rank = mpi_rank;
	rec.trajectory_id = trajectory_id;
	rec.completion_wall_time_sec = completion_wall_time_sec;
	rec.t_evap = event_observed ? lifetime_unbinding : std::numeric_limits<double>::quiet_NaN();
	rec.t_capture = bincount.t_capture;
	rec.t_first_unbinding_scatter = std::isfinite(bincount.t_first_unbinding_scatter) ? bincount.t_first_unbinding_scatter : -1.0;
	rec.t_final_unbinding_scatter = std::isfinite(bincount.t_final_unbinding_scatter) ? bincount.t_final_unbinding_scatter : -1.0;
	rec.t_boundary_escape = std::isfinite(bincount.t_boundary_escape) ? bincount.t_boundary_escape : -1.0;
	rec.t_termination = t_termination;
	rec.observed_lifetime = observed_lifetime;
	rec.lifetime_unbinding = (std::isfinite(lifetime_unbinding) && lifetime_unbinding >= 0.0) ? lifetime_unbinding : -1.0;
	rec.lifetime_boundary = (std::isfinite(lifetime_boundary) && lifetime_boundary >= 0.0) ? lifetime_boundary : -1.0;
	rec.r_first_negative_km = bincount.r_first_negative_km;
	rec.E_first_negative_eV = bincount.E_first_negative_eV;
	rec.dE_first_negative_from_prev_eV = bincount.dE_first_negative_from_prev_eV;
	rec.r_first_unbinding_km = bincount.r_first_unbinding_km;
	rec.E_first_unbinding_eV = bincount.E_first_unbinding_eV;
	rec.r_final_unbinding_km = bincount.r_final_unbinding_km;
	rec.E_final_unbinding_eV = bincount.E_final_unbinding_eV;
	rec.r_boundary_escape_km = bincount.r_boundary_escape_km;
	rec.vr_boundary_escape_km_s = bincount.vr_boundary_escape_km_s;
	rec.E_boundary_escape_eV = bincount.E_boundary_escape_eV;
	rec.event_observed = event_observed;
	rec.boundary_escape_observed = survival_valid && bincount.boundary_escape_observed;
	rec.survival_valid = survival_valid;
	rec.numerically_invalid_escape = numerically_invalid_escape;
	// Normal evaporation output only accepts completed unbinding events. A
	// captured trajectory that has not completed is invalid for this estimator,
	// not a right-censored observation.
	rec.censored = false;
	rec.truncated = false;
	rec.termination_reason = bincount.termination_reason;
	rec.max_free_energy_drift_eV = bincount.max_free_energy_drift_eV;
	rec.max_free_energy_drift_rel = bincount.max_free_energy_drift_rel;
	rec.number_of_scatterings = bincount.number_of_scatterings;
	rec.number_of_bound_to_unbound = bincount.number_of_bound_to_unbound;
	rec.number_of_recaptures = bincount.number_of_recaptures;
	rec.number_of_integrator_steps_after_capture = bincount.number_of_integrator_steps_after_capture;
	rec.min_energy_after_capture_eV = bincount.min_energy_after_capture_eV;
	rec.max_radius_after_capture_km = bincount.max_radius_after_capture_km;
	rec.time_inside_sun_after_capture_sec = bincount.time_inside_sun_after_capture_sec;
	rec.time_outside_sun_after_capture_sec = bincount.time_outside_sun_after_capture_sec;
	return true;
}

const char* Termination_Reason_Key(TrajectoryTerminationReason reason)
{
	switch(reason)
	{
		case TrajectoryTerminationReason::OutwardEscape: return "outward_escape";
		case TrajectoryTerminationReason::Scatter: return "scatter";
		case TrajectoryTerminationReason::WallTimeLimit: return "wall_time_limit";
		case TrajectoryTerminationReason::MaxFreeSteps: return "max_free_steps";
		case TrajectoryTerminationReason::MaxScatterings: return "max_scatterings";
		case TrajectoryTerminationReason::NonFiniteState: return "non_finite_state";
		case TrajectoryTerminationReason::SpeedLimit: return "speed_limit";
		case TrajectoryTerminationReason::NumericalFailure: return "numerical_failure";
		case TrajectoryTerminationReason::CaptureMode: return "capture_mode";
		case TrajectoryTerminationReason::EnergyDriftEscape: return "energy_drift_escape";
		case TrajectoryTerminationReason::Unknown:
		default: return "unknown";
	}
}

bool Diagnostic_Trace_Selected(uint64_t trace_seed, int rank, uint64_t trajectory_id, double rate)
{
	if(rate <= 0.0)
		return false;
	if(rate >= 1.0)
		return true;
	uint64_t value = trace_seed ^ (static_cast<uint64_t>(static_cast<uint32_t>(rank)) << 32) ^ trajectory_id;
	value += 0x9e3779b97f4a7c15ULL;
	value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
	value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
	value ^= value >> 31;
	const double unit = static_cast<double>(value >> 11) * (1.0 / 9007199254740992.0);
	return unit < rate;
}

const char* Diagnostic_Status(const EvaporationRecord& rec)
{
	if(rec.event_observed)
		return "escaped";
	if(Is_Numerical_Termination(rec.termination_reason))
		return "numerical_failure";
	return "censored";
}

bool Make_Compact_Evaporation_Event(const EvaporationRecord& rec, CompactEvaporationEvent& event)
{
	if(!Is_Completed_Evaporation_Record(rec))
		return false;

	event.rank = rec.rank;
	event.trajectory_id = rec.trajectory_id;
	event.completion_wall_time_sec = rec.completion_wall_time_sec;
	event.lifetime_unbinding = rec.lifetime_unbinding;
	return true;
}

struct BinomialRateEstimate
{
	double rate = 0.0;
	double standard_error = 0.0;
	double ci_lower = 0.0;
	double ci_upper = 0.0;
};

BinomialRateEstimate Estimate_Binomial_Rate(uint64_t trials, uint64_t successes)
{
	BinomialRateEstimate estimate;
	if(trials == 0)
		return estimate;

	const double N = static_cast<double>(trials);
	estimate.rate = static_cast<double>(std::min(successes, trials)) / N;
	estimate.standard_error = sqrt(estimate.rate * (1.0 - estimate.rate) / N);
	const double z = 1.96;
	const double denominator = 1.0 + z * z / N;
	const double center = estimate.rate + z * z / (2.0 * N);
	const double spread = z * sqrt(estimate.rate * (1.0 - estimate.rate) / N + z * z / (4.0 * N * N));
	estimate.ci_lower = std::max(0.0, (center - spread) / denominator);
	estimate.ci_upper = std::min(1.0, (center + spread) / denominator);
	return estimate;
}

const char* Stop_Reason_Key(SimulationStopReason reason)
{
	switch(reason)
	{
		case SimulationStopReason::MaxTrajectoriesReached:
			return "max_trajectories_reached";
		case SimulationStopReason::CaptureTargetNotReached:
			return "capture_target_not_reached";
		case SimulationStopReason::InitialShiftFailureFractionExceeded:
			return "initial_shift_failure_fraction_exceeded";
		case SimulationStopReason::InvalidTrajectoryFractionExceeded:
			return "invalid_trajectory_fraction_exceeded";
		case SimulationStopReason::None:
		default:
			return "none";
	}
}

const char* Stop_Reason_Display(SimulationStopReason reason)
{
	switch(reason)
	{
		case SimulationStopReason::MaxTrajectoriesReached:
			return "max_trajectories reached";
		case SimulationStopReason::CaptureTargetNotReached:
			return "capture target not reached";
		case SimulationStopReason::InitialShiftFailureFractionExceeded:
			return "initial shift failure fraction exceeded";
		case SimulationStopReason::InvalidTrajectoryFractionExceeded:
			return "numerical/computational invalid-trajectory fraction exceeded";
		case SimulationStopReason::None:
		default:
			return "none";
	}
}

double Snapshot_Bin_Error(double sum, double sum_sq, double count)
{
	if(count <= 1.0)
		return 0.0;

	double mean = sum / count;
	double variance = sum_sq / count - mean * mean;
	if(variance < 0.0)
		variance = 0.0;

	return sqrt(count * variance);
}

void Write_Report_Header(
	std::ofstream& file,
	double mass_gev,
	double sigma_cm2,
	uint64_t total_trajectories,
	uint64_t captured_particles,
	SimulationStopReason stop_reason)
{
	file << "# DM_mass_GeV = " << std::scientific << std::setprecision(6) << mass_gev << "\n";
	file << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(6) << sigma_cm2 << "\n";
	file << "# total_trajectories = " << total_trajectories << "\n";
	file << "# captured_particles = " << captured_particles << "\n";
	file << "# bincount_integration = " << BincountIntegrationScheme() << "\n";
	file << "# bincount_dense_position_tolerance_km = "
	     << std::scientific << std::setprecision(6)
	     << BincountDensePositionToleranceKm() << "\n";

	const BinomialRateEstimate raw = Estimate_Binomial_Rate(total_trajectories, captured_particles);
	file << "# capture_rate = " << std::fixed << std::setprecision(8) << raw.rate << "\n";
	file << "# capture_rate_raw = " << std::fixed << std::setprecision(8) << raw.rate << "\n";
	file << "# capture_rate_err = " << std::fixed << std::setprecision(8) << raw.standard_error << "\n";
	file << "# capture_rate_raw_err = " << std::fixed << std::setprecision(8) << raw.standard_error << "\n";
	file << "# capture_rate_CI_95_lower = " << std::fixed << std::setprecision(8) << raw.ci_lower << "\n";
	file << "# capture_rate_CI_95_upper = " << std::fixed << std::setprecision(8) << raw.ci_upper << "\n";
	file << "# capture_rate_raw_CI_95_lower = " << std::fixed << std::setprecision(8) << raw.ci_lower << "\n";
	file << "# capture_rate_raw_CI_95_upper = " << std::fixed << std::setprecision(8) << raw.ci_upper << "\n";

	if(stop_reason != SimulationStopReason::None)
		file << "# EARLY_STOP: " << Stop_Reason_Key(stop_reason) << "\n";
}

std::string Join_Path(const std::string& directory, const std::string& name)
{
	if(directory.empty() || directory.back() == '/')
		return directory + name;
	return directory + "/" + name;
}

std::string Evaporation_Log_Path_From_Output_Dir(const std::string& output_dir)
{
	return Join_Path(output_dir, "evaporation_times.txt");
}

void Write_Evaporation_Log_File_Header(std::ofstream& file, double mass_gev, double sigma_cm2)
{
	file << "# DaMaSCUS-SUN evaporation times\n";
	file << "# format_version = 3\n";
	file << "# DM_mass_GeV = " << std::scientific << std::setprecision(6) << mass_gev << "\n";
	file << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(6) << sigma_cm2 << "\n";
	file << "# sorted_by = lifetime_unbinding_sec rank trajectory_id\n";
	file << "# rank trajectory_id lifetime_unbinding_sec\n";
}

void Write_Evaporation_Log_Event(std::ostream& file, const CompactEvaporationEvent& event)
{
	file << event.rank << "\t" << event.trajectory_id << "\t" << std::scientific << std::setprecision(10)
	     << event.lifetime_unbinding << "\n";
}

bool Evaporation_Event_Order(const CompactEvaporationEvent& lhs, const CompactEvaporationEvent& rhs)
{
	if(lhs.lifetime_unbinding != rhs.lifetime_unbinding)
		return lhs.lifetime_unbinding < rhs.lifetime_unbinding;
	if(lhs.rank != rhs.rank)
		return lhs.rank < rhs.rank;
	return lhs.trajectory_id < rhs.trajectory_id;
}

void Write_Evaporation_Log_Events(std::ostream& file, const std::vector<CompactEvaporationEvent>& events)
{
	std::vector<CompactEvaporationEvent> sorted_events = events;
	std::sort(sorted_events.begin(), sorted_events.end(), Evaporation_Event_Order);
	for(const auto& event : sorted_events)
	{
		if(Is_Completed_Evaporation_Event(event))
			Write_Evaporation_Log_Event(file, event);
	}
}

bool Write_Final_Evaporation_Time_File(const std::string& path, double mass_gev, double sigma_cm2, const std::vector<CompactEvaporationEvent>& events)
{
	const std::string tmp_path = path + ".final.tmp." + std::to_string(getpid());
	std::ofstream file(tmp_path, std::ios::out | std::ios::trunc);
	if(!file.is_open())
		return false;
	Write_Evaporation_Log_File_Header(file, mass_gev, sigma_cm2);
	Write_Evaporation_Log_Events(file, events);
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

bool Remove_Path_Recursive(const std::string& path)
{
	struct stat info;
	if(lstat(path.c_str(), &info) != 0)
		return errno == ENOENT;

	if(S_ISDIR(info.st_mode))
	{
		DIR* dir = opendir(path.c_str());
		if(dir == NULL)
			return false;

		bool success = true;
		struct dirent* entry = NULL;
		while((entry = readdir(dir)) != NULL)
		{
			const std::string name = entry->d_name;
			if(name == "." || name == "..")
				continue;
			if(!Remove_Path_Recursive(Join_Path(path, name)))
				success = false;
		}
		closedir(dir);

		if(rmdir(path.c_str()) != 0 && errno != ENOENT)
			success = false;

		return success;
	}

	return std::remove(path.c_str()) == 0 || errno == ENOENT;
}

bool Clear_Directory_Contents(const std::string& directory)
{
	DIR* dir = opendir(directory.c_str());
	if(dir == NULL)
		return errno == ENOENT;

	bool success = true;
	struct dirent* entry = NULL;
	while((entry = readdir(dir)) != NULL)
	{
		const std::string name = entry->d_name;
		if(name == "." || name == "..")
			continue;
		if(!Remove_Path_Recursive(Join_Path(directory, name)))
			success = false;
	}
	closedir(dir);
	return success;
}

bool Ensure_Directory_Exists(const std::string& directory)
{
	if(directory.empty())
		return false;

	std::string normalized = directory;
	while(normalized.size() > 1 && normalized.back() == '/')
		normalized.pop_back();

	struct stat info;
	if(stat(normalized.c_str(), &info) == 0)
		return S_ISDIR(info.st_mode);

	const std::size_t separator = normalized.find_last_of('/');
	if(separator != std::string::npos)
	{
		std::string parent = normalized.substr(0, separator);
		if(parent.empty())
			parent = "/";
		if(parent != normalized && !Ensure_Directory_Exists(parent))
			return false;
	}

	if(mkdir(normalized.c_str(), 0755) != 0 && errno != EEXIST)
		return false;
	if(stat(normalized.c_str(), &info) != 0)
		return false;
	return S_ISDIR(info.st_mode);
}

void Build_MPI_Gatherv_Layout(const std::vector<int>& item_counts, int fields_per_item, std::vector<int>& recv_counts, std::vector<int>& displacements, int& total_items)
{
	recv_counts.resize(item_counts.size());
	displacements.resize(item_counts.size());
	for(size_t i = 0; i < item_counts.size(); i++)
	{
		recv_counts[i] = item_counts[i] * fields_per_item;
		displacements[i] = (i == 0) ? 0 : displacements[i-1] + recv_counts[i-1];
	}
	total_items = std::accumulate(item_counts.begin(), item_counts.end(), 0);
}
}

bool TrajectoryTraceSelected(uint64_t trace_seed, int rank, uint64_t trajectory_id, double rate)
{
	return Diagnostic_Trace_Selected(trace_seed, rank, trajectory_id, rate);
}

Simulation_Data::Simulation_Data(unsigned int sample_size, unsigned int max_trajectories, double u_min, unsigned int iso_rings)
: requested_captured_particles(sample_size),
  normal_mode_mpi_sync_interval(NORMAL_MODE_MPI_SYNC_INTERVAL_FALLBACK),
  number_of_trajectories(0), number_of_free_particles(0), number_of_reflected_particles(0), number_of_captured_particles(0),
  number_of_completed_outward_escapes(0),
  number_of_complete_evaporation_particles(0), number_of_censored_captured_particles(0),
  number_of_invalid_survival_captured_particles(0),
  number_of_initial_shift_failures(0), number_of_final_reflection_shift_failures(0), number_of_numerical_failures(0),
  number_of_computational_truncations(0),
  total_number_of_scatterings(0), average_number_of_scatterings(0.0),
  mpi_sync_rounds(0), final_mpi_round_trajectories(0), mpi_tail_trajectories(0), capture_target_overshoot(0),
  computing_time(0.0), early_stopped(false), early_stop_reason(SimulationStopReason::None),
  mpi_rank(0), mpi_processes(1), isoreflection_rings(iso_rings), minimum_speed_threshold(u_min),
  number_of_data_points(std::vector<unsigned long int>(iso_rings, 0)),
  data(iso_rings, std::vector<libphysica::DataPoint>())
{
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if(max_trajectories == 0)
        max_trajectories_per_rank = std::numeric_limits<unsigned long int>::max();  // no limit
    else
        max_trajectories_per_rank = (max_trajectories + mpi_processes - 1) / mpi_processes;

    captured_dt_hist.fill(0.0);
    captured_v2dt_hist.fill(0.0);
    not_captured_dt_hist.fill(0.0);
    not_captured_v2dt_hist.fill(0.0);
	captured_dt_sq_hist.fill(0.0);
	captured_v2dt_sq_hist.fill(0.0);
	not_captured_dt_sq_hist.fill(0.0);
	not_captured_v2dt_sq_hist.fill(0.0);
}

void Simulation_Data::Configure(double initial_radius, unsigned int min_scattering, unsigned long int max_scattering, unsigned long int max_free_steps)
{
	initial_and_final_radius      = initial_radius;
	minimum_number_of_scatterings = min_scattering;
	maximum_number_of_scatterings = max_scattering;
	maximum_free_time_steps       = max_free_steps;
}

void Simulation_Data::Configure_Trajectory_Diagnostics(const TrajectoryDiagnosticConfig& config)
{
	if(!std::isfinite(config.trace_rate) || config.trace_rate < 0.0 || config.trace_rate > 1.0)
		throw std::invalid_argument("trajectory diagnostic trace rate must lie in [0, 1]");
	trajectory_diagnostic_config = config;
	if(!trajectory_diagnostic_config.summary_enabled)
	{
		trajectory_diagnostic_config.events_enabled = false;
		trajectory_diagnostic_config.trace_rate = 0.0;
	}
	evaporation_diagnostics_enabled = trajectory_diagnostic_config.summary_enabled;
}

void Simulation_Data::Generate_Data(obscura::DM_Particle& DM, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, SnapshotConfig snapshot_cfg, unsigned int fixed_seed, bool capture_mode)
{
	if(capture_mode)
		snapshot_cfg.enabled = false;
	if(snapshot_cfg.enabled && !IsValidSnapshotIntervalSeconds(snapshot_cfg.interval_seconds))
		throw std::invalid_argument("snapshot interval must be a positive integer number of seconds");
	if(snapshot_cfg.enabled)
	{
		int mpi_thread_level = MPI_THREAD_SINGLE;
		MPI_Query_thread(&mpi_thread_level);
		if(mpi_thread_level < MPI_THREAD_FUNNELED)
			throw std::runtime_error("snapshot heartbeat requires MPI_THREAD_FUNNELED or stronger thread support");
	}
	normal_mode_mpi_sync_interval = capture_mode ? 0UL : Normal_Mode_MPI_Sync_Interval(In_Units(DM.Sigma_Proton(), cm * cm));
	diagnostic_base_seed = fixed_seed;
	diagnostic_run_id = 0;
	if(evaporation_diagnostics_enabled && mpi_rank == 0)
	{
		diagnostic_run_id = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::microseconds>(
			std::chrono::system_clock::now().time_since_epoch()).count());
	}
	if(evaporation_diagnostics_enabled)
		MPI_Bcast(&diagnostic_run_id, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

	auto time_start = std::chrono::steady_clock::now();
	unsigned long int local_captured = 0;
	unsigned long int local_total = 0;
	bool initial_shift_failure_warning_emitted = false;

	// Configure the simulator
	Trajectory_Simulator simulator(solar_model, maximum_free_time_steps, maximum_number_of_scatterings, initial_and_final_radius);
	simulator.max_trajectory_wall_time_sec = snapshot_cfg.max_trajectory_wall_time_sec;
	simulator.Enable_Capture_Mode(capture_mode);
	if(fixed_seed != 0)
	{
		const uint64_t rank_seed = static_cast<uint64_t>(fixed_seed) + 1000003ULL * static_cast<uint64_t>(mpi_rank);
		simulator.Fix_PRNG_Seed(static_cast<unsigned int>(rank_seed));
	}

	// Snapshot configuration
	const double snapshot_interval = (snapshot_cfg.interval_seconds > 0.0) ? snapshot_cfg.interval_seconds : 60.0;
	const double snapshot_mass_gev = In_Units(DM.mass, GeV);
	const double snapshot_sigma_cm2 = In_Units(DM.Sigma_Proton(), cm * cm);
	const std::string output_root = g_top_level_dir + "results_" + std::to_string(log10(snapshot_mass_gev)) + "_" + std::to_string(log10(snapshot_sigma_cm2)) + "/";
	std::string snapshot_root;
	std::string rank_snapshot_dir;
	uint64_t snapshot_run_id = 0;
	if(snapshot_cfg.enabled)
	{
		snapshot_root = output_root + "snapshot/";
		rank_snapshot_dir = snapshot_root + "rank_snapshot/";
		int snapshot_init_ok = 1;
		if(mpi_rank == 0)
		{
			if(!Ensure_Directory_Exists(output_root))
				snapshot_init_ok = 0;

			if(snapshot_init_ok
			   && (!Ensure_Directory_Exists(snapshot_root)
			       || !Clear_Directory_Contents(snapshot_root)
			       || !Ensure_Directory_Exists(snapshot_root)
			       || !Ensure_Directory_Exists(rank_snapshot_dir)))
			{
				std::cerr << "Warning in Generate_Data(): failed to initialize snapshot directory; snapshots disabled for this run." << std::endl;
				snapshot_init_ok = 0;
			}
			if(snapshot_init_ok)
			{
				// Do not leave a previous run's final reports beside the new snapshots.
				std::remove(Join_Path(output_root, "bincount.txt").c_str());
				std::remove(Evaporation_Log_Path_From_Output_Dir(output_root).c_str());
				std::remove(Join_Path(output_root, "evaporation_diagnostics.txt").c_str());
			}

			if(snapshot_init_ok)
				snapshot_run_id = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
		}
		MPI_Bcast(&snapshot_init_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if(!snapshot_init_ok)
		{
			snapshot_cfg.enabled = false;
		}
		else
		{
			MPI_Barrier(MPI_COMM_WORLD);
			int rank_snapshot_dirs_ok =
			    (Ensure_Directory_Exists(snapshot_root)
			     && Ensure_Directory_Exists(rank_snapshot_dir)) ? 1 : 0;
			MPI_Allreduce(MPI_IN_PLACE, &rank_snapshot_dirs_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
			if(!rank_snapshot_dirs_ok)
			{
				if(mpi_rank == 0)
					std::cerr << "Warning in Generate_Data(): failed to initialize rank snapshot directories; snapshots disabled for this run." << std::endl;
				snapshot_cfg.enabled = false;
			}
			else
				MPI_Bcast(&snapshot_run_id, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
		}
	}
	std::unique_ptr<SnapshotSharedState> snapshot_state;
	std::unique_ptr<SnapshotRecorder> snapshot_recorder;
	std::unique_ptr<SnapshotHeartbeat> snapshot_heartbeat;
	if(snapshot_cfg.enabled)
	{
		int local_objects_ready = 1;
		try
		{
			snapshot_state.reset(new SnapshotSharedState());
			snapshot_state->Initialize(snapshot_run_id, mpi_rank);
			snapshot_recorder.reset(new SnapshotRecorder(*snapshot_state));
			snapshot_heartbeat.reset(new SnapshotHeartbeat(
				*snapshot_state,
				mpi_rank,
				mpi_processes,
				snapshot_run_id,
				snapshot_root,
				rank_snapshot_dir,
				snapshot_interval,
				snapshot_mass_gev,
				snapshot_sigma_cm2));
		}
		catch(const std::exception& error)
		{
			local_objects_ready = 0;
			std::cerr << "Warning in Generate_Data(): rank " << mpi_rank
			          << " failed to initialize snapshot heartbeat: " << error.what() << std::endl;
		}
		catch(...)
		{
			local_objects_ready = 0;
			std::cerr << "Warning in Generate_Data(): rank " << mpi_rank
			          << " failed to initialize snapshot heartbeat with an unknown exception." << std::endl;
		}

		int all_objects_ready = local_objects_ready;
		MPI_Allreduce(MPI_IN_PLACE, &all_objects_ready, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
		if(all_objects_ready)
		{
			simulator.Set_Snapshot_Recorder(snapshot_recorder.get());
			MPI_Barrier(MPI_COMM_WORLD);
			time_start = std::chrono::steady_clock::now();
			int all_threads_started = snapshot_heartbeat->Start(time_start) ? 1 : 0;
			MPI_Allreduce(MPI_IN_PLACE, &all_threads_started, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
			if(!all_threads_started)
			{
				snapshot_heartbeat->Stop();
				simulator.Set_Snapshot_Recorder(nullptr);
				MPI_Barrier(MPI_COMM_WORLD);
				if(mpi_rank == 0 && !Clear_Directory_Contents(snapshot_root))
					std::cerr << "Warning in Generate_Data(): failed to clear snapshot files after heartbeat startup rollback." << std::endl;
				MPI_Barrier(MPI_COMM_WORLD);
				snapshot_cfg.enabled = false;
			}
		}

		if(!all_objects_ready || !snapshot_cfg.enabled)
		{
			if(snapshot_heartbeat)
				snapshot_heartbeat->Stop();
			simulator.Set_Snapshot_Recorder(nullptr);
			snapshot_heartbeat.reset();
			snapshot_recorder.reset();
			snapshot_state.reset();
			snapshot_cfg.enabled = false;
			if(mpi_rank == 0)
				std::cerr << "Warning in Generate_Data(): snapshot heartbeat setup failed on at least one rank; "
				             "snapshots are disabled for this run." << std::endl;
		}
	}
	if(!snapshot_cfg.enabled)
		time_start = std::chrono::steady_clock::now();

	auto elapsed_since_start = [&]()
	{
		return 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - time_start).count();
	};
	early_stopped = false;
	early_stop_reason = SimulationStopReason::None;

	unsigned long int global_captured = 0;
	const bool unlimited_trajectory_budget =
	    max_trajectories_per_rank == std::numeric_limits<unsigned long int>::max();
	const std::vector<double> progress_milestones = {0.0, 0.01, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 1.0};
	size_t next_progress_milestone = 0;
	bool progress_line_printed = false;
	unsigned long int last_progress_line_captured = 0;
	auto print_progress_update = [&](unsigned long int captured_particles, bool force)
	{
		if(mpi_rank != 0)
			return;

		const double denominator = static_cast<double>(std::max(1u, requested_captured_particles));
		const double progress = std::min(1.0, static_cast<double>(captured_particles) / denominator);

		bool should_print = force;
		while(next_progress_milestone < progress_milestones.size()
		      && progress + 1.0e-12 >= progress_milestones[next_progress_milestone])
		{
			should_print = true;
			next_progress_milestone++;
		}
		if(!should_print)
			return;
		if(force && progress_line_printed && captured_particles == last_progress_line_captured)
			return;

		const double time_elapsed = elapsed_since_start();
		const double captured_particle_rate = (time_elapsed > 0.0) ? static_cast<double>(captured_particles) / time_elapsed : 0.0;
		libphysica::Print_Progress_Bar(progress, 0, 44, time_elapsed);
		std::cout << " captured_particles=" << captured_particles << "/" << requested_captured_particles
		          << " captured_particle_rate[1/s]=" << libphysica::Round(captured_particle_rate)
		          << std::endl;
		progress_line_printed = true;
		last_progress_line_captured = captured_particles;
	};
	print_progress_update(global_captured, true);

	const unsigned long int mpi_sync_interval = capture_mode ? CAPTURE_MODE_MPI_SYNC_INTERVAL : normal_mode_mpi_sync_interval;
	auto select_trajectory_batch = [&](unsigned long int remaining_captures, unsigned long int& total_assigned)
	{
		total_assigned = 0;
		unsigned long int local_round_capacity = 0;
		if(local_total < max_trajectories_per_rank)
		{
			const unsigned long int local_capacity = max_trajectories_per_rank - local_total;
			local_round_capacity = std::min(local_capacity, mpi_sync_interval);
		}

		std::vector<unsigned long int> round_capacities(mpi_processes, 0);
		MPI_Allgather(&local_round_capacity, 1, MPI_UNSIGNED_LONG, round_capacities.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

		unsigned long int total_round_capacity = 0;
		for(const auto capacity : round_capacities)
			total_round_capacity += capacity;
		if(total_round_capacity == 0)
			return 0UL;

		const unsigned long int attempt_budget = capture_mode ? std::min(remaining_captures, total_round_capacity) : total_round_capacity;
		unsigned long int local_batch = 0;
		for(int rank = 0; rank < mpi_processes; rank++)
		{
			if(total_assigned >= attempt_budget)
				break;
			const unsigned long int assigned = std::min(round_capacities[rank], attempt_budget - total_assigned);
			if(rank == mpi_rank)
				local_batch = assigned;
			total_assigned += assigned;
		}
		return local_batch;
	};

	while(global_captured < requested_captured_particles)
	{
		const unsigned long int remaining_captures = requested_captured_particles - global_captured;
		unsigned long int assigned_trajectories = 0;
		const unsigned long int local_trajectories_this_round =
		    select_trajectory_batch(remaining_captures, assigned_trajectories);
		if(assigned_trajectories == 0)
		{
			early_stopped = true;
			early_stop_reason = SimulationStopReason::MaxTrajectoriesReached;
			break;
		}
		mpi_sync_rounds++;
		final_mpi_round_trajectories = assigned_trajectories;

		auto simulate_one_trajectory = [&]()
		{
			const unsigned long int trajectory_id = local_total + 1;
			const bool trace_selected = trajectory_diagnostic_config.events_enabled
			                         && Diagnostic_Trace_Selected(trajectory_diagnostic_config.trace_seed,
			                                                      mpi_rank, trajectory_id,
			                                                      trajectory_diagnostic_config.trace_rate);
			const std::string rng_state_before_initial_conditions = trace_selected
			                                                       ? simulator.Serialize_PRNG_State()
			                                                       : std::string();
			simulator.Enable_Diagnostic_Trace(trace_selected);
			Event IC = Initial_Conditions(halo_model, solar_model, simulator.PRNG);
			const bool initial_shift_ok = Hyperbolic_Kepler_Shift(IC, initial_and_final_radius);
			const std::string rng_state_before_simulation = trace_selected
			                                                ? simulator.Serialize_PRNG_State()
			                                                : std::string();
			if(!initial_shift_ok)
			{
				number_of_initial_shift_failures++;
				number_of_numerical_failures++;
			}
			if(initial_shift_ok)
				simulator.current_trajectory_id = trajectory_id - 1;
			Trajectory_Result trajectory = initial_shift_ok
			    ? simulator.Simulate(IC, DM, mpi_rank)
			    : [&]() {
				TrajectoryBincount failed_bincount;
				failed_bincount.termination_reason = TrajectoryTerminationReason::NumericalFailure;
				failed_bincount.survival_valid = false;
				return Trajectory_Result(IC, IC, 0, failed_bincount);
			}();
			const double trajectory_completion_wall_time_sec =
			    (!capture_mode && trajectory.bincount.is_captured) ? elapsed_since_start() : 0.0;

			local_total++;
			number_of_trajectories++;
			total_number_of_scatterings += static_cast<uint64_t>(trajectory.number_of_scatterings);
			const bool completed_outward_escape = Completed_Outward_Escape(trajectory.bincount.termination_reason);
			if(initial_shift_ok && Is_Numerical_Termination(trajectory.bincount.termination_reason))
				number_of_numerical_failures++;
			if(Is_Computational_Truncation(trajectory.bincount.termination_reason))
				number_of_computational_truncations++;
			std::vector<SnapshotEvaporationProgressEntry> trajectory_snapshot_evaporation_events;

			if(trajectory.bincount.is_captured)
			{
				number_of_captured_particles++;
				local_captured++;
				if(!capture_mode)
				{
					if(trajectory.bincount.event_observed)
						number_of_complete_evaporation_particles++;
					else
						number_of_invalid_survival_captured_particles++;
				}

				if(!capture_mode)
				{
					if(trace_selected)
					{
						trajectory_diagnostic_events.insert(trajectory_diagnostic_events.end(),
						                                    trajectory.diagnostic_events.begin(),
						                                    trajectory.diagnostic_events.end());
						TrajectoryReplayRecord replay;
						replay.rank = mpi_rank;
						replay.trajectory_id = trajectory_id;
						replay.initial_time_s = In_Units(IC.time, sec);
						for(size_t component = 0; component < 3; component++)
						{
							replay.initial_position_km[component] = In_Units(IC.position[component], km);
							replay.initial_velocity_km_s[component] = In_Units(IC.velocity[component], km / sec);
						}
						replay.rng_state_before_initial_conditions = rng_state_before_initial_conditions;
						replay.rng_state_before_simulation = rng_state_before_simulation;
						trajectory_replay_records.push_back(std::move(replay));
					}
					for(int b = 0; b < NUM_BINS; b++)
					{
						captured_dt_hist[b]   += trajectory.bincount.dt_hist[b];
						captured_v2dt_hist[b] += trajectory.bincount.v2dt_hist[b];
						captured_dt_sq_hist[b]   += trajectory.bincount.dt_hist[b] * trajectory.bincount.dt_hist[b];
						captured_v2dt_sq_hist[b] += trajectory.bincount.v2dt_hist[b] * trajectory.bincount.v2dt_hist[b];
					}

					EvaporationRecord rec;
					if(Build_Evaporation_Record(trajectory.bincount, mpi_rank, number_of_trajectories, trajectory_completion_wall_time_sec, rec))
					{
						if(evaporation_diagnostics_enabled)
							evaporation_records.push_back(rec);
						CompactEvaporationEvent event;
						if(Make_Compact_Evaporation_Event(rec, event))
						{
							compact_evaporation_events.push_back(event);
							trajectory_snapshot_evaporation_events.push_back(MakeSnapshotEvaporationProgressEntry(event));
						}
					}
				}
			}
			else if(completed_outward_escape)
			{
				number_of_completed_outward_escapes++;
				if(!capture_mode)
				{
					for(int b = 0; b < NUM_BINS; b++)
					{
						not_captured_dt_hist[b]   += trajectory.bincount.dt_hist[b];
						not_captured_v2dt_hist[b] += trajectory.bincount.v2dt_hist[b];
						not_captured_dt_sq_hist[b]   += trajectory.bincount.dt_hist[b] * trajectory.bincount.dt_hist[b];
						not_captured_v2dt_sq_hist[b] += trajectory.bincount.v2dt_hist[b] * trajectory.bincount.v2dt_hist[b];
					}

					if(trajectory.Particle_Free())
						number_of_free_particles++;
					else if(trajectory.Particle_Reflected())
					{
						number_of_reflected_particles++;
						if(Hyperbolic_Kepler_Shift(trajectory.final_event, 1.0 * AU))
						{
							const double v_final = trajectory.final_event.Speed();
							if(trajectory.number_of_scatterings >= minimum_number_of_scatterings
							   && v_final > KDE_boundary_correction_factor * minimum_speed_threshold)
							{
								const unsigned int isoreflection_ring =
								    (isoreflection_rings == 1)
								    ? 0
								    : trajectory.final_event.Isoreflection_Ring(obscura::Sun_Velocity(), isoreflection_rings);
								data[isoreflection_ring].push_back(libphysica::DataPoint(v_final));
							}
						}
						else
						{
							number_of_final_reflection_shift_failures++;
							number_of_numerical_failures++;
						}
					}
				}
			}

			if(snapshot_state)
			{
				const bool count_as_captured_bincount_sample = (!capture_mode && trajectory.bincount.is_captured);
				const bool count_as_not_captured_bincount_sample = (!capture_mode && completed_outward_escape && !trajectory.bincount.is_captured);
				snapshot_state->RecordCompletedTrajectory(
					trajectory.bincount,
					count_as_captured_bincount_sample,
					count_as_not_captured_bincount_sample,
					trajectory_snapshot_evaporation_events);
			}
		};

		for(unsigned long int trajectory_in_round = 0; trajectory_in_round < local_trajectories_this_round; trajectory_in_round++)
			simulate_one_trajectory();

		// A blocking reduction here makes every early rank idle behind the rank
		// whose final scheduled trajectory has the longest tail. In normal mode,
		// rendezvous nonblockingly and use that tail time for more independent
		// trajectories. All extra results are included in the reduction below.
		if(!capture_mode)
		{
			const unsigned long int extra_capacity =
			    (local_total < max_trajectories_per_rank)
			    ? (max_trajectories_per_rank - local_total)
			    : 0UL;
			mpi_tail_trajectories += Perform_MPI_Tail_Work(
				MPI_COMM_WORLD,
				extra_capacity,
				simulate_one_trajectory);
		}

		const std::array<unsigned long int, 5> local_progress = {
		    local_captured, number_of_trajectories, number_of_initial_shift_failures,
		    number_of_numerical_failures, number_of_computational_truncations};
		std::array<unsigned long int, 5> global_progress = {{0, 0, 0, 0, 0}};
		MPI_Allreduce(local_progress.data(), global_progress.data(), static_cast<int>(global_progress.size()), MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		global_captured = global_progress[0];
		const unsigned long int global_attempted_trajectories = global_progress[1];
		const unsigned long int global_initial_shift_failures = global_progress[2];
		const unsigned long int global_numerical_failures = global_progress[3];
		const unsigned long int global_computational_truncations = global_progress[4];

		if(global_attempted_trajectories > 0)
		{
			const double initial_shift_failure_fraction =
			    static_cast<double>(global_initial_shift_failures) / static_cast<double>(global_attempted_trajectories);
			if(initial_shift_failure_fraction > NUMERICAL_FAILURE_ABORT_FRACTION)
			{
				if(mpi_rank == 0)
					std::cerr << "Error in Generate_Data(): initial Kepler shift failure fraction "
					          << initial_shift_failure_fraction
					          << " exceeds abort threshold "
					          << NUMERICAL_FAILURE_ABORT_FRACTION
					          << ". Stopping this run to avoid biased capture statistics." << std::endl;
				early_stopped = true;
				early_stop_reason = SimulationStopReason::InitialShiftFailureFractionExceeded;
				break;
			}
			if(!initial_shift_failure_warning_emitted
			   && initial_shift_failure_fraction > NUMERICAL_FAILURE_WARNING_FRACTION)
			{
				if(mpi_rank == 0)
					std::cerr << "Warning in Generate_Data(): initial Kepler shift failure fraction "
					          << initial_shift_failure_fraction
					          << " exceeds warning threshold "
					          << NUMERICAL_FAILURE_WARNING_FRACTION
					          << "." << std::endl;
				initial_shift_failure_warning_emitted = true;
			}

			const double invalid_trajectory_fraction =
			    static_cast<double>(global_numerical_failures + global_computational_truncations)
			    / static_cast<double>(global_attempted_trajectories);
			if(unlimited_trajectory_budget
			   && invalid_trajectory_fraction > NUMERICAL_FAILURE_ABORT_FRACTION)
			{
				if(mpi_rank == 0)
					std::cerr << "Error in Generate_Data(): invalid trajectory fraction "
					          << invalid_trajectory_fraction
					          << " exceeds abort threshold "
					          << NUMERICAL_FAILURE_ABORT_FRACTION
					          << " with an unlimited trajectory budget. Stopping to avoid a non-progressing run."
					          << std::endl;
				early_stopped = true;
				early_stop_reason = SimulationStopReason::InvalidTrajectoryFractionExceeded;
				break;
			}
		}

		print_progress_update(global_captured, false);
	}
	capture_target_overshoot = (global_captured > requested_captured_particles)
	                         ? global_captured - requested_captured_particles
	                         : 0UL;

	if(global_captured < requested_captured_particles)
	{
		early_stopped = true;
		if(early_stop_reason == SimulationStopReason::None)
			early_stop_reason = SimulationStopReason::CaptureTargetNotReached;
	}

	auto time_end  = std::chrono::steady_clock::now();
	computing_time = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();

	if(snapshot_heartbeat)
	{
		snapshot_heartbeat->MarkDoneAndWriteFinal(computing_time);

		MPI_Barrier(MPI_COMM_WORLD);
		double final_snapshot_elapsed = computing_time;
		MPI_Allreduce(MPI_IN_PLACE, &final_snapshot_elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if(mpi_rank == 0)
			snapshot_heartbeat->FinalizeAfterAllRanksDone(final_snapshot_elapsed);
	}

	print_progress_update(global_captured, global_captured > 0 || !early_stopped);
	if(mpi_rank == 0)
		std::cout << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	Perform_MPI_Reductions(capture_mode);
}

void Simulation_Data::Perform_MPI_Reductions(bool capture_mode)
{
	MPI_Trace_Point(mpi_rank, "enter Perform_MPI_Reductions");
	std::array<unsigned long int, 8> primary_counters = {{
	    number_of_trajectories,
	    number_of_captured_particles,
	    number_of_completed_outward_escapes,
	    number_of_initial_shift_failures,
	    number_of_final_reflection_shift_failures,
	    number_of_numerical_failures,
	    number_of_computational_truncations,
	    mpi_tail_trajectories}};
	MPI_Trace_Point(mpi_rank, "before allreduce primary counters");
	MPI_Allreduce(MPI_IN_PLACE, primary_counters.data(), static_cast<int>(primary_counters.size()), MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	number_of_trajectories = primary_counters[0];
	number_of_captured_particles = primary_counters[1];
	number_of_completed_outward_escapes = primary_counters[2];
	number_of_initial_shift_failures = primary_counters[3];
	number_of_final_reflection_shift_failures = primary_counters[4];
	number_of_numerical_failures = primary_counters[5];
	number_of_computational_truncations = primary_counters[6];
	mpi_tail_trajectories = primary_counters[7];
	MPI_Trace_Point(mpi_rank, "before allreduce total scatterings");
	MPI_Allreduce(MPI_IN_PLACE, &total_number_of_scatterings, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
	average_number_of_scatterings = (number_of_trajectories > 0)
	                               ? static_cast<double>(total_number_of_scatterings) / number_of_trajectories
	                               : 0.0;

	int global_stop_reason = static_cast<int>(early_stop_reason);
	MPI_Trace_Point(mpi_rank, "before allreduce early_stop_reason");
	MPI_Allreduce(MPI_IN_PLACE, &global_stop_reason, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	early_stop_reason = static_cast<SimulationStopReason>(global_stop_reason);
	early_stopped = early_stop_reason != SimulationStopReason::None;

	if(capture_mode)
	{
		MPI_Trace_Point(mpi_rank, "before allreduce computing_time capture");
		MPI_Allreduce(MPI_IN_PLACE, &computing_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		MPI_Trace_Point(mpi_rank, "leave Perform_MPI_Reductions capture");
		return;
	}

	std::array<unsigned long int, 5> result_counters = {{
	    number_of_free_particles,
	    number_of_reflected_particles,
	    number_of_complete_evaporation_particles,
	    number_of_censored_captured_particles,
	    number_of_invalid_survival_captured_particles}};
	MPI_Trace_Point(mpi_rank, "before allreduce result counters");
	MPI_Allreduce(MPI_IN_PLACE, result_counters.data(), static_cast<int>(result_counters.size()), MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	number_of_free_particles = result_counters[0];
	number_of_reflected_particles = result_counters[1];
	number_of_complete_evaporation_particles = result_counters[2];
	number_of_censored_captured_particles = result_counters[3];
	number_of_invalid_survival_captured_particles = result_counters[4];

	// Reduce bincount histograms across all ranks
	MPI_Trace_Point(mpi_rank, "before allreduce captured_dt_hist");
	MPI_Allreduce(MPI_IN_PLACE, captured_dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce captured_v2dt_hist");
	MPI_Allreduce(MPI_IN_PLACE, captured_v2dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce not_captured_dt_hist");
	MPI_Allreduce(MPI_IN_PLACE, not_captured_dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce not_captured_v2dt_hist");
	MPI_Allreduce(MPI_IN_PLACE, not_captured_v2dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce captured_dt_sq_hist");
	MPI_Allreduce(MPI_IN_PLACE, captured_dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce captured_v2dt_sq_hist");
	MPI_Allreduce(MPI_IN_PLACE, captured_v2dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce not_captured_dt_sq_hist");
	MPI_Allreduce(MPI_IN_PLACE, not_captured_dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce not_captured_v2dt_sq_hist");
	MPI_Allreduce(MPI_IN_PLACE, not_captured_v2dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "after histogram allreduces");
	if(evaporation_diagnostics_enabled)
	{
		const int local_evap_count = static_cast<int>(evaporation_records.size());
		std::vector<int> evap_counts(mpi_processes, 0);
		MPI_Trace_Point(mpi_rank, "before gather diagnostic evap counts");
		MPI_Gather(&local_evap_count, 1, MPI_INT, mpi_rank == 0 ? evap_counts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

		constexpr int EVAPORATION_MPI_INT_FIELDS = 8;
		constexpr int EVAPORATION_MPI_UINT_FIELDS = 5;
		constexpr int EVAPORATION_MPI_DOUBLE_FIELDS = 26;
		std::vector<int> local_evap_ints(local_evap_count * EVAPORATION_MPI_INT_FIELDS);
		std::vector<unsigned long long> local_evap_uints(local_evap_count * EVAPORATION_MPI_UINT_FIELDS);
		std::vector<double> local_evap_doubles(local_evap_count * EVAPORATION_MPI_DOUBLE_FIELDS);
		for(int i = 0; i < local_evap_count; i++)
		{
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i] = evaporation_records[i].rank;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 1] = evaporation_records[i].event_observed ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 2] = evaporation_records[i].boundary_escape_observed ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 3] = evaporation_records[i].survival_valid ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 4] = evaporation_records[i].numerically_invalid_escape ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 5] = evaporation_records[i].censored ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 6] = evaporation_records[i].truncated ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 7] = TerminationReason_Index(evaporation_records[i].termination_reason);

			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i] = static_cast<unsigned long long>(evaporation_records[i].trajectory_id);
			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 1] = static_cast<unsigned long long>(evaporation_records[i].number_of_scatterings);
			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 2] = static_cast<unsigned long long>(evaporation_records[i].number_of_bound_to_unbound);
			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 3] = static_cast<unsigned long long>(evaporation_records[i].number_of_recaptures);
			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 4] = static_cast<unsigned long long>(evaporation_records[i].number_of_integrator_steps_after_capture);

			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i] = evaporation_records[i].completion_wall_time_sec;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 1] = evaporation_records[i].t_evap;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 2] = evaporation_records[i].t_capture;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 3] = evaporation_records[i].t_final_unbinding_scatter;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 4] = evaporation_records[i].t_boundary_escape;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 5] = evaporation_records[i].t_termination;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 6] = evaporation_records[i].observed_lifetime;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 7] = evaporation_records[i].lifetime_unbinding;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 8] = evaporation_records[i].lifetime_boundary;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 9] = evaporation_records[i].r_first_negative_km;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 10] = evaporation_records[i].E_first_negative_eV;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 11] = evaporation_records[i].dE_first_negative_from_prev_eV;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 12] = evaporation_records[i].max_free_energy_drift_eV;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 13] = evaporation_records[i].max_free_energy_drift_rel;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 14] = evaporation_records[i].t_first_unbinding_scatter;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 15] = evaporation_records[i].r_first_unbinding_km;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 16] = evaporation_records[i].E_first_unbinding_eV;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 17] = evaporation_records[i].r_final_unbinding_km;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 18] = evaporation_records[i].E_final_unbinding_eV;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 19] = evaporation_records[i].r_boundary_escape_km;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 20] = evaporation_records[i].vr_boundary_escape_km_s;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 21] = evaporation_records[i].E_boundary_escape_eV;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 22] = evaporation_records[i].min_energy_after_capture_eV;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 23] = evaporation_records[i].max_radius_after_capture_km;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 24] = evaporation_records[i].time_inside_sun_after_capture_sec;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 25] = evaporation_records[i].time_outside_sun_after_capture_sec;
		}

		std::vector<int> recv_counts, displacements;
		std::vector<int> global_evap_ints;
		std::vector<unsigned long long> global_evap_uints;
		std::vector<double> global_evap_doubles;
		int total_evap = 0;
		if(mpi_rank == 0)
		{
			Build_MPI_Gatherv_Layout(evap_counts, EVAPORATION_MPI_INT_FIELDS, recv_counts, displacements, total_evap);
			global_evap_ints.resize(total_evap * EVAPORATION_MPI_INT_FIELDS);
		}
		MPI_Trace_Point(mpi_rank, "before gatherv diagnostic evap ints");
		MPI_Gatherv(local_evap_ints.data(), local_evap_count * EVAPORATION_MPI_INT_FIELDS, MPI_INT,
		            mpi_rank == 0 ? global_evap_ints.data() : nullptr,
		            mpi_rank == 0 ? recv_counts.data() : nullptr,
		            mpi_rank == 0 ? displacements.data() : nullptr,
		            MPI_INT, 0, MPI_COMM_WORLD);

		if(mpi_rank == 0)
		{
			Build_MPI_Gatherv_Layout(evap_counts, EVAPORATION_MPI_UINT_FIELDS, recv_counts, displacements, total_evap);
			global_evap_uints.resize(total_evap * EVAPORATION_MPI_UINT_FIELDS);
		}
		MPI_Trace_Point(mpi_rank, "before gatherv diagnostic evap uints");
		MPI_Gatherv(local_evap_uints.data(), local_evap_count * EVAPORATION_MPI_UINT_FIELDS, MPI_UNSIGNED_LONG_LONG,
		            mpi_rank == 0 ? global_evap_uints.data() : nullptr,
		            mpi_rank == 0 ? recv_counts.data() : nullptr,
		            mpi_rank == 0 ? displacements.data() : nullptr,
		            MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

		if(mpi_rank == 0)
		{
			Build_MPI_Gatherv_Layout(evap_counts, EVAPORATION_MPI_DOUBLE_FIELDS, recv_counts, displacements, total_evap);
			global_evap_doubles.resize(total_evap * EVAPORATION_MPI_DOUBLE_FIELDS);
		}
		MPI_Trace_Point(mpi_rank, "before gatherv diagnostic evap doubles");
		MPI_Gatherv(local_evap_doubles.data(), local_evap_count * EVAPORATION_MPI_DOUBLE_FIELDS, MPI_DOUBLE,
		            mpi_rank == 0 ? global_evap_doubles.data() : nullptr,
		            mpi_rank == 0 ? recv_counts.data() : nullptr,
		            mpi_rank == 0 ? displacements.data() : nullptr,
		            MPI_DOUBLE, 0, MPI_COMM_WORLD);

		evaporation_records.clear();
		compact_evaporation_events.clear();
		if(mpi_rank == 0)
		{
			evaporation_records.resize(total_evap);
			for(int i = 0; i < total_evap; i++)
			{
				evaporation_records[i].rank = global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i];
				evaporation_records[i].trajectory_id = static_cast<unsigned long int>(global_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i]);
				evaporation_records[i].completion_wall_time_sec = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i];
				evaporation_records[i].t_evap = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 1];
				evaporation_records[i].t_capture = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 2];
				evaporation_records[i].t_final_unbinding_scatter = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 3];
				evaporation_records[i].t_boundary_escape = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 4];
				evaporation_records[i].t_termination = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 5];
				evaporation_records[i].observed_lifetime = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 6];
				evaporation_records[i].lifetime_unbinding = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 7];
				evaporation_records[i].lifetime_boundary = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 8];
				evaporation_records[i].r_first_negative_km = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 9];
				evaporation_records[i].E_first_negative_eV = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 10];
				evaporation_records[i].dE_first_negative_from_prev_eV = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 11];
				evaporation_records[i].event_observed = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 1] != 0);
				evaporation_records[i].boundary_escape_observed = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 2] != 0);
				evaporation_records[i].survival_valid = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 3] != 0);
				evaporation_records[i].numerically_invalid_escape = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 4] != 0);
				evaporation_records[i].censored = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 5] != 0);
				evaporation_records[i].truncated = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 6] != 0);
				evaporation_records[i].termination_reason = static_cast<TrajectoryTerminationReason>(TerminationReason_Index(static_cast<TrajectoryTerminationReason>(global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 7])));
				evaporation_records[i].max_free_energy_drift_eV = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 12];
				evaporation_records[i].max_free_energy_drift_rel = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 13];
				evaporation_records[i].number_of_scatterings = static_cast<unsigned long int>(global_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 1]);
				evaporation_records[i].number_of_bound_to_unbound = static_cast<unsigned long int>(global_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 2]);
				evaporation_records[i].number_of_recaptures = static_cast<unsigned long int>(global_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 3]);
				evaporation_records[i].number_of_integrator_steps_after_capture = static_cast<unsigned long int>(global_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 4]);
				evaporation_records[i].t_first_unbinding_scatter = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 14];
				evaporation_records[i].r_first_unbinding_km = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 15];
				evaporation_records[i].E_first_unbinding_eV = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 16];
				evaporation_records[i].r_final_unbinding_km = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 17];
				evaporation_records[i].E_final_unbinding_eV = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 18];
				evaporation_records[i].r_boundary_escape_km = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 19];
				evaporation_records[i].vr_boundary_escape_km_s = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 20];
				evaporation_records[i].E_boundary_escape_eV = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 21];
				evaporation_records[i].min_energy_after_capture_eV = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 22];
				evaporation_records[i].max_radius_after_capture_km = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 23];
				evaporation_records[i].time_inside_sun_after_capture_sec = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 24];
				evaporation_records[i].time_outside_sun_after_capture_sec = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 25];
			}
			compact_evaporation_events.reserve(evaporation_records.size());
			for(const auto& record : evaporation_records)
			{
				CompactEvaporationEvent event;
				if(Make_Compact_Evaporation_Event(record, event))
					compact_evaporation_events.push_back(event);
			}
		}
	}
	else
	{
		const int local_evap_count = static_cast<int>(compact_evaporation_events.size());
		std::vector<int> evap_counts(mpi_processes, 0);
		MPI_Trace_Point(mpi_rank, "before gather compact evap counts");
		MPI_Gather(&local_evap_count, 1, MPI_INT, mpi_rank == 0 ? evap_counts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

		int max_evap_count = 0;
		MPI_Trace_Point(mpi_rank, "before allreduce compact evap max count");
		MPI_Allreduce(&local_evap_count, &max_evap_count, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		std::vector<CompactEvaporationEvent> local_evap_padded(max_evap_count);
		for(int i = 0; i < local_evap_count; i++)
			local_evap_padded[i] = compact_evaporation_events[i];

		std::vector<CompactEvaporationEvent> global_evap_padded;
		if(mpi_rank == 0)
			global_evap_padded.resize(static_cast<size_t>(mpi_processes) * static_cast<size_t>(max_evap_count));

		const int padded_evap_bytes = max_evap_count * static_cast<int>(sizeof(CompactEvaporationEvent));
		MPI_Trace_Point(mpi_rank, "before gather compact evap padded bytes");
		MPI_Gather(max_evap_count == 0 ? nullptr : static_cast<void*>(local_evap_padded.data()),
		            padded_evap_bytes,
		            MPI_BYTE,
		            mpi_rank == 0 && max_evap_count > 0 ? static_cast<void*>(global_evap_padded.data()) : nullptr,
		            padded_evap_bytes,
		            MPI_BYTE, 0, MPI_COMM_WORLD);

		compact_evaporation_events.clear();
		if(mpi_rank == 0)
		{
			int total_evap = std::accumulate(evap_counts.begin(), evap_counts.end(), 0);
			compact_evaporation_events.reserve(total_evap);
			for(int rank = 0; rank < mpi_processes; rank++)
			{
				const size_t rank_offset = static_cast<size_t>(rank) * static_cast<size_t>(max_evap_count);
				for(int i = 0; i < evap_counts[rank]; i++)
					compact_evaporation_events.push_back(global_evap_padded[rank_offset + static_cast<size_t>(i)]);
			}
		}
	}

	if(evaporation_diagnostics_enabled)
	{
		static_assert(std::is_trivially_copyable<TrajectoryDiagnosticEvent>::value,
		              "trajectory diagnostic events must remain MPI byte-copyable");
		const int local_event_count = static_cast<int>(trajectory_diagnostic_events.size());
		std::vector<int> event_counts(mpi_processes, 0);
		MPI_Gather(&local_event_count, 1, MPI_INT,
		           mpi_rank == 0 ? event_counts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
		std::vector<int> event_recv_counts;
		std::vector<int> event_displacements;
		int total_event_count = 0;
		std::vector<TrajectoryDiagnosticEvent> global_events;
		if(mpi_rank == 0)
		{
			Build_MPI_Gatherv_Layout(event_counts, static_cast<int>(sizeof(TrajectoryDiagnosticEvent)),
			                         event_recv_counts, event_displacements, total_event_count);
			global_events.resize(static_cast<size_t>(total_event_count));
		}
		MPI_Gatherv(trajectory_diagnostic_events.empty() ? nullptr : trajectory_diagnostic_events.data(),
		            local_event_count * static_cast<int>(sizeof(TrajectoryDiagnosticEvent)), MPI_BYTE,
		            mpi_rank == 0 && !global_events.empty() ? global_events.data() : nullptr,
		            mpi_rank == 0 ? event_recv_counts.data() : nullptr,
		            mpi_rank == 0 ? event_displacements.data() : nullptr,
		            MPI_BYTE, 0, MPI_COMM_WORLD);
		trajectory_diagnostic_events.clear();
		if(mpi_rank == 0)
			trajectory_diagnostic_events.swap(global_events);

		auto encode_rng_state = [](std::string state) {
			std::replace(state.begin(), state.end(), ' ', ',');
			return state;
		};
		std::ostringstream local_replay_stream;
		local_replay_stream << std::scientific << std::setprecision(17);
		for(const auto& replay : trajectory_replay_records)
		{
			local_replay_stream << replay.rank << '\t' << replay.trajectory_id << '\t' << replay.initial_time_s;
			for(double value : replay.initial_position_km)
				local_replay_stream << '\t' << value;
			for(double value : replay.initial_velocity_km_s)
				local_replay_stream << '\t' << value;
			local_replay_stream << '\t' << encode_rng_state(replay.rng_state_before_initial_conditions)
			                    << '\t' << encode_rng_state(replay.rng_state_before_simulation) << '\n';
		}
		const std::string local_replay_text = local_replay_stream.str();
		const int local_replay_bytes = static_cast<int>(local_replay_text.size());
		std::vector<int> replay_byte_counts(mpi_processes, 0);
		MPI_Gather(&local_replay_bytes, 1, MPI_INT,
		           mpi_rank == 0 ? replay_byte_counts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
		std::vector<int> replay_recv_counts;
		std::vector<int> replay_displacements;
		int total_replay_bytes = 0;
		std::vector<char> global_replay_text;
		if(mpi_rank == 0)
		{
			Build_MPI_Gatherv_Layout(replay_byte_counts, 1, replay_recv_counts, replay_displacements, total_replay_bytes);
			global_replay_text.resize(static_cast<size_t>(total_replay_bytes));
		}
		MPI_Gatherv(local_replay_text.empty() ? nullptr : local_replay_text.data(), local_replay_bytes, MPI_CHAR,
		            mpi_rank == 0 && !global_replay_text.empty() ? global_replay_text.data() : nullptr,
		            mpi_rank == 0 ? replay_recv_counts.data() : nullptr,
		            mpi_rank == 0 ? replay_displacements.data() : nullptr,
		            MPI_CHAR, 0, MPI_COMM_WORLD);
		trajectory_replay_records.clear();
		if(mpi_rank == 0)
		{
			std::string combined(global_replay_text.begin(), global_replay_text.end());
			std::istringstream lines(combined);
			std::string line;
			while(std::getline(lines, line))
			{
				if(line.empty())
					continue;
				TrajectoryReplayRecord replay;
				std::string initial_rng_state;
				std::string simulation_rng_state;
				std::istringstream fields(line);
				fields >> replay.rank >> replay.trajectory_id >> replay.initial_time_s;
				for(double& value : replay.initial_position_km)
					fields >> value;
				for(double& value : replay.initial_velocity_km_s)
					fields >> value;
				fields >> initial_rng_state >> simulation_rng_state;
				if(!fields)
					throw std::runtime_error("Perform_MPI_Reductions(): failed to parse trajectory replay record.");
				std::replace(initial_rng_state.begin(), initial_rng_state.end(), ',', ' ');
				std::replace(simulation_rng_state.begin(), simulation_rng_state.end(), ',', ' ');
				replay.rng_state_before_initial_conditions = std::move(initial_rng_state);
				replay.rng_state_before_simulation = std::move(simulation_rng_state);
				trajectory_replay_records.push_back(std::move(replay));
			}
		}
	}

	MPI_Trace_Point(mpi_rank, "before allreduce computing_time final");
	MPI_Allreduce(MPI_IN_PLACE, &computing_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "leave Perform_MPI_Reductions");
}

void Simulation_Data::Write_Output_Files(const std::string& output_dir, obscura::DM_Particle& DM)
{
	if(mpi_rank != 0)
		return;

	if(!Ensure_Directory_Exists(output_dir))
	{
		std::cerr << "Warning in Write_Output_Files(): failed to create output directory " << output_dir << std::endl;
		return;
	}

	double mass_gev = In_Units(DM.mass, GeV);
	double sigma_cm2 = In_Units(DM.Sigma_Proton(), cm * cm);
	auto unresolved_not_captured_trajectories = [&]() {
		const unsigned long int classified_trajectories = Valid_Trajectories();
		return (number_of_trajectories > classified_trajectories)
		       ? (number_of_trajectories - classified_trajectories)
		       : 0UL;
	};

	auto write_header = [&](std::ofstream& f) {
		Write_Report_Header(f, mass_gev, sigma_cm2, number_of_trajectories, number_of_captured_particles, early_stop_reason);
		f << "# valid_trajectories = " << Valid_Trajectories() << "\n";
		f << "# completed_outward_escapes = " << number_of_completed_outward_escapes << "\n";
		f << "# unresolved_not_captured_trajectories = " << unresolved_not_captured_trajectories() << "\n";
		f << "# numerical_failures = " << number_of_numerical_failures << "\n";
		f << "# computational_truncations = " << number_of_computational_truncations << "\n";
		f << "# initial_shift_failures = " << number_of_initial_shift_failures << "\n";
		f << "# final_reflection_shift_failures = " << number_of_final_reflection_shift_failures << "\n";
		f << "# normal_mode_mpi_sync_interval = " << normal_mode_mpi_sync_interval << "\n";
		f << "# mpi_sync_rounds = " << mpi_sync_rounds << "\n";
		f << "# final_mpi_round_trajectories = " << final_mpi_round_trajectories << "\n";
		f << "# mpi_tail_trajectories = " << mpi_tail_trajectories << "\n";
		f << "# capture_target_overshoot = " << capture_target_overshoot << "\n";
		f << "# total_scatterings = " << total_number_of_scatterings << "\n";
		f << "# average_scatterings = " << std::setprecision(12) << average_number_of_scatterings << "\n";
		f << "# simulation_time_seconds = " << std::setprecision(12) << computing_time << "\n";
		const BinomialRateEstimate valid = Estimate_Binomial_Rate(Valid_Trajectories(), number_of_captured_particles);
		f << "# capture_rate_valid = " << std::fixed << std::setprecision(8) << valid.rate << "\n";
		f << "# capture_rate_valid_err = " << std::fixed << std::setprecision(8) << valid.standard_error << "\n";
		f << "# capture_rate_valid_CI_95_lower = " << std::fixed << std::setprecision(8) << valid.ci_lower << "\n";
		f << "# capture_rate_valid_CI_95_upper = " << std::fixed << std::setprecision(8) << valid.ci_upper << "\n";
		f << "# numerical_failure_rate = " << std::fixed << std::setprecision(8) << Numerical_Failure_Ratio() << "\n";
	};

	std::remove((output_dir + "/captured_bincount.txt").c_str());
	std::remove((output_dir + "/not_captured_bincount.txt").c_str());
	std::remove((output_dir + "/evaporation_diagnostics.txt").c_str());
	std::remove((output_dir + "/run_metadata.json").c_str());
	std::remove((output_dir + "/trajectory_summary.tsv").c_str());
	std::remove((output_dir + "/trajectory_events.tsv").c_str());
	std::remove((output_dir + "/evaporation_" + "summary.txt").c_str());
	std::remove((output_dir + "/evaporation_" + "mode_summary.txt").c_str());
	std::remove((output_dir + "/evaporation_" + "mode_" + "bincount.txt").c_str());
	std::remove((output_dir + "/computation_" + "time_summary.txt").c_str());

	// 1. Merged bincount
	{
		std::ofstream f(output_dir + "/bincount.txt");
		write_header(f);
		const unsigned long int physical_not_captured_particles = number_of_completed_outward_escapes;
		const unsigned long int excluded_not_captured_particles = unresolved_not_captured_trajectories();
		f << "# not_captured_bincount_samples = " << physical_not_captured_particles << "\n";
		if(excluded_not_captured_particles > 0)
			f << "# excluded_incomplete_not_captured = " << excluded_not_captured_particles << "\n";
		f << "# bin_index  cap_dt[s]  cap_v2dt[km2/s]  cap_err_dt[s]  cap_err_v2dt[km2/s]  not_cap_dt[s]  not_cap_v2dt[km2/s]  not_cap_err_dt[s]  not_cap_err_v2dt[km2/s]\n";
		double N_cap = static_cast<double>(number_of_captured_particles);
		double N_nc = static_cast<double>(physical_not_captured_particles);
		for(int b = 0; b < NUM_BINS; b++)
		{
			double cap_err_dt = Snapshot_Bin_Error(captured_dt_hist[b], captured_dt_sq_hist[b], N_cap);
			double cap_err_v2dt = Snapshot_Bin_Error(captured_v2dt_hist[b], captured_v2dt_sq_hist[b], N_cap);
			double not_cap_err_dt = Snapshot_Bin_Error(not_captured_dt_hist[b], not_captured_dt_sq_hist[b], N_nc);
			double not_cap_err_v2dt = Snapshot_Bin_Error(not_captured_v2dt_hist[b], not_captured_v2dt_sq_hist[b], N_nc);

			f << b << "\t" << std::scientific << std::setprecision(10)
			  << captured_dt_hist[b] << "\t" << captured_v2dt_hist[b]
			  << "\t" << cap_err_dt << "\t" << cap_err_v2dt
			  << "\t" << not_captured_dt_hist[b] << "\t" << not_captured_v2dt_hist[b]
			  << "\t" << not_cap_err_dt << "\t" << not_cap_err_v2dt << "\n";
		}
		f.close();
	}

	// 2. Final evaporation-time list.  This is intentionally the only final
	// evaporation report; snapshot files are intermediate progress reports.
	bool evaporation_times_ok = Write_Final_Evaporation_Time_File(Evaporation_Log_Path_From_Output_Dir(output_dir), mass_gev, sigma_cm2, compact_evaporation_events);
	if(!evaporation_times_ok)
		std::cerr << "Warning in Write_Output_Files(): failed to write evaporation_times.txt" << std::endl;

	if(evaporation_diagnostics_enabled)
	{
		std::vector<EvaporationRecord> records = evaporation_records;
		std::sort(records.begin(), records.end(), [](const EvaporationRecord& lhs, const EvaporationRecord& rhs) {
			if(lhs.rank != rhs.rank)
				return lhs.rank < rhs.rank;
			return lhs.trajectory_id < rhs.trajectory_id;
		});

		bool unique_keys = true;
		bool time_invariants = true;
		bool legacy_reconciliation = true;
		bool escape_radius_invariant = true;
		bool residence_time_invariant = true;
		bool event_sequence_invariant = true;
		bool event_count_invariant = true;
		bool trace_selection_invariant = true;
		bool replay_state_invariant = true;
		std::sort(trajectory_diagnostic_events.begin(), trajectory_diagnostic_events.end(),
		          [](const TrajectoryDiagnosticEvent& lhs, const TrajectoryDiagnosticEvent& rhs) {
			          if(lhs.rank != rhs.rank) return lhs.rank < rhs.rank;
			          if(lhs.trajectory_id != rhs.trajectory_id) return lhs.trajectory_id < rhs.trajectory_id;
			          return lhs.event_index < rhs.event_index;
		          });
		std::unordered_map<std::string, std::vector<const TrajectoryDiagnosticEvent*> > events_by_key;
		for(const auto& event : trajectory_diagnostic_events)
			events_by_key[std::to_string(event.rank) + ":" + std::to_string(event.trajectory_id)].push_back(&event);
		std::unordered_map<std::string, const TrajectoryReplayRecord*> replay_by_key;
		for(const auto& replay : trajectory_replay_records)
			replay_by_key[std::to_string(replay.rank) + ":" + std::to_string(replay.trajectory_id)] = &replay;
		std::unordered_set<std::string> completed_keys;
		for(const auto& event : compact_evaporation_events)
		{
			completed_keys.insert(std::to_string(event.rank) + ":" + std::to_string(event.trajectory_id));
		}
		for(size_t index = 0; index < records.size(); index++)
		{
			const auto& rec = records[index];
			if(index > 0 && rec.rank == records[index - 1].rank
			   && rec.trajectory_id == records[index - 1].trajectory_id)
				unique_keys = false;
			if(rec.t_first_unbinding_scatter >= 0.0 && rec.t_first_unbinding_scatter < rec.t_capture)
				time_invariants = false;
			if(rec.event_observed
			   && (!(rec.t_final_unbinding_scatter >= rec.t_capture)
			       || !(rec.t_boundary_escape >= rec.t_final_unbinding_scatter)))
				time_invariants = false;
			const std::string key = std::to_string(rec.rank) + ":" + std::to_string(rec.trajectory_id);
			const bool expected_trace = trajectory_diagnostic_config.events_enabled
			                         && Diagnostic_Trace_Selected(trajectory_diagnostic_config.trace_seed,
			                                                      rec.rank, rec.trajectory_id,
			                                                      trajectory_diagnostic_config.trace_rate);
			const bool actual_trace = events_by_key.count(key) != 0;
			if(expected_trace != actual_trace)
				trace_selection_invariant = false;
			if(expected_trace)
			{
				const auto replay_it = replay_by_key.find(key);
				if(replay_it == replay_by_key.end()
				   || replay_it->second->rng_state_before_initial_conditions.empty()
				   || replay_it->second->rng_state_before_simulation.empty())
					replay_state_invariant = false;
			}
			if(rec.event_observed)
			{
				const double escape_radius_rsun = rec.r_boundary_escape_km / R_SUN_KM;
				if(!std::isfinite(escape_radius_rsun)
				   || std::fabs(escape_radius_rsun - In_Units(initial_and_final_radius, rSun)) > 1.0e-10)
					escape_radius_invariant = false;
			}
			const double observation_end = rec.event_observed ? rec.t_boundary_escape : rec.t_termination;
			const double observed_duration = observation_end - rec.t_capture;
			const double residence_duration = rec.time_inside_sun_after_capture_sec
			                                + rec.time_outside_sun_after_capture_sec;
			const double residence_tolerance = std::max(1.0e-6, 1.0e-10 * std::fabs(observed_duration));
			if(!std::isfinite(observed_duration) || observed_duration < 0.0
			   || !std::isfinite(residence_duration)
			   || std::fabs(residence_duration - observed_duration) > residence_tolerance)
				residence_time_invariant = false;
			if(expected_trace)
			{
				const auto event_it = events_by_key.find(key);
				if(event_it == events_by_key.end())
					event_count_invariant = false;
				else
				{
					uint64_t candidate_count = 0;
					uint64_t recapture_count = 0;
					const auto& trajectory_events = event_it->second;
					for(size_t event_index = 0; event_index < trajectory_events.size(); event_index++)
					{
						const auto& event = *trajectory_events[event_index];
						candidate_count += event.event_type == TrajectoryDiagnosticEventType::CandidateUnbinding ? 1 : 0;
						recapture_count += event.event_type == TrajectoryDiagnosticEventType::Recapture ? 1 : 0;
						if(event.event_index != event_index)
							event_sequence_invariant = false;
						if(event_index > 0)
						{
							const auto& previous = *trajectory_events[event_index - 1];
							if(event.t_s < previous.t_s
							   || event.scatter_index < previous.scatter_index
							   || event.step_index < previous.step_index)
								event_sequence_invariant = false;
						}
					}
					if(candidate_count != rec.number_of_bound_to_unbound
					   || recapture_count != rec.number_of_recaptures)
						event_count_invariant = false;
				}
			}
			if(rec.event_observed != (completed_keys.count(key) != 0))
				legacy_reconciliation = false;
		}
		if(completed_keys.size() != number_of_complete_evaporation_particles)
			legacy_reconciliation = false;

		const std::string git_commit = GIT_COMMIT_HASH;
		const bool git_dirty = git_commit.find("-dirty") != std::string::npos;
#if defined(__clang__)
		const std::string compiler = std::string("Clang ") + __clang_version__;
#elif defined(__GNUC__)
		const std::string compiler = std::string("GCC ") + __VERSION__;
#else
		const std::string compiler = "unknown";
#endif
#ifdef NDEBUG
		const char* build_type = "Release";
#else
		const char* build_type = "Debug";
#endif
		{
			std::ofstream metadata(output_dir + "/run_metadata.json");
			metadata << "{\n"
			         << "  \"schema_version\": \"trajectory-diagnostic-v2\",\n"
			         << "  \"run_id\": \"" << diagnostic_run_id << "\",\n"
			         << "  \"git_branch\": \"" << GIT_BRANCH << "\",\n"
			         << "  \"git_commit\": \"" << git_commit << "\",\n"
			         << "  \"git_dirty\": " << (git_dirty ? "true" : "false") << ",\n"
			         << "  \"build_type\": \"" << build_type << "\",\n"
			         << "  \"compiler\": \"" << compiler << "\",\n"
			         << "  \"mass_GeV\": " << std::scientific << std::setprecision(17) << mass_gev << ",\n"
			         << "  \"sigma_cm2\": " << sigma_cm2 << ",\n"
			         << "  \"R_escape_Rsun\": " << In_Units(initial_and_final_radius, rSun) << ",\n"
			         << "  \"mpi_size\": " << mpi_processes << ",\n"
			         << "  \"rng_algorithm\": \"std::mt19937\",\n"
			         << "  \"base_seed\": " << diagnostic_base_seed << ",\n"
			         << "  \"rng_stream_definition\": \"rank seed = base_seed + 1000003*rank; rng_counter is the rank-local trajectory id\",\n"
			         << "  \"integrator\": \"adaptive RK45\",\n"
			         << "  \"rk_position_tolerance_km\": " << RK45PositionToleranceKm() << ",\n"
			         << "  \"rk_velocity_tolerance_km_s\": " << RK45VelocityToleranceKmPerSec() << ",\n"
			         << "  \"rk_phase_tolerance\": " << RK45PhaseTolerance() << ",\n"
			         << "  \"rk_absolute_max_step_s\": " << RK45AbsoluteMaxStepSec() << ",\n"
			         << "  \"bincount_integration\": \"" << BincountIntegrationScheme() << "\",\n"
			         << "  \"bincount_dense_position_tolerance_km\": " << BincountDensePositionToleranceKm() << ",\n"
			         << "  \"interpolation_points\": " << trajectory_diagnostic_config.interpolation_points << ",\n"
			         << "  \"max_optical_depth_step\": " << NormalModeMaxOpticalDepthStep() << ",\n"
			         << "  \"optical_depth_relative_tolerance\": " << OpticalDepthRelativeTolerance() << ",\n"
			         << "  \"energy_definition\": \"0.5*m_chi*(v^2-v_escape(r)^2); bound iff energy < 0\",\n"
			         << "  \"energy_unit\": \"eV\",\n"
			         << "  \"length_unit\": \"Rsun\",\n"
			         << "  \"velocity_unit\": \"km/s\",\n"
			         << "  \"angular_momentum_unit\": \"km^2/s\",\n"
			         << "  \"n_scatter_total_definition\": \"all trajectory scatters, including scatters before first capture\",\n"
			         << "  \"stop_conditions\": {\"max_free_steps\": " << maximum_free_time_steps
			         << ", \"max_scatterings\": " << maximum_number_of_scatterings << "},\n"
			         << "  \"trace_selection_rule\": \"stable hash(trace_seed, rank, trajectory_id) < trace_rate; no physics RNG calls\",\n"
			         << "  \"trace_seed\": " << trajectory_diagnostic_config.trace_seed << ",\n"
			         << "  \"trace_rate\": " << std::fixed << std::setprecision(8) << trajectory_diagnostic_config.trace_rate << ",\n"
			         << "  \"events_scope\": \"real-time scatter, state transition, solar crossing, escape and termination events\",\n"
			         << "  \"event_order_definition\": \"event_index is strictly increasing; physical time, scatter_index and step_index are non-decreasing because paired pre/post events can share an integration instant\",\n"
			         << "  \"summary_unique_keys\": " << (unique_keys ? "true" : "false") << ",\n"
			         << "  \"summary_time_invariants\": " << (time_invariants ? "true" : "false") << ",\n"
			         << "  \"escape_radius_invariant\": " << (escape_radius_invariant ? "true" : "false") << ",\n"
			         << "  \"residence_time_invariant\": " << (residence_time_invariant ? "true" : "false") << ",\n"
			         << "  \"event_sequence_invariant\": " << (event_sequence_invariant ? "true" : "false") << ",\n"
			         << "  \"event_count_invariant\": " << (event_count_invariant ? "true" : "false") << ",\n"
			         << "  \"trace_selection_invariant\": " << (trace_selection_invariant ? "true" : "false") << ",\n"
			         << "  \"replay_state_invariant\": " << (replay_state_invariant ? "true" : "false") << ",\n"
			         << "  \"legacy_evaporation_reconciliation\": " << (legacy_reconciliation ? "true" : "false") << "\n"
			         << "}\n";
		}

		auto nan_if_missing = [](double value) {
			return (std::isfinite(value) && value >= 0.0)
			       ? value : std::numeric_limits<double>::quiet_NaN();
		};
		{
			std::ofstream summary(output_dir + "/trajectory_summary.tsv");
			summary << "run_id\trank\ttrajectory_id\trng_stream\trng_counter\tstatus\ttermination_reason\tevent_observed"
			        << "\tt_capture_s\tt_first_unbinding_s\tt_final_unbinding_s\tt_escape_s\tt_censor_s"
			        << "\tlifetime_first_unbinding_s\tlifetime_final_unbinding_s\tlifetime_validated_escape_s"
			        << "\tr_capture_Rsun\tE_capture_eV\tr_first_unbinding_Rsun\tE_first_unbinding_eV"
			        << "\tr_final_unbinding_Rsun\tE_final_unbinding_eV\tr_escape_Rsun\tvr_escape_km_s\tE_escape_eV"
			        << "\tn_scatter_total\tn_bound_to_unbound\tn_recapture\tmin_energy_after_capture_eV\tmax_r_Rsun"
			        << "\ttime_inside_sun_s\ttime_outside_sun_s\tn_integrator_steps\tmax_abs_ballistic_energy_drift_eV"
			        << "\tmax_scaled_ballistic_energy_drift\ttrace_written"
			        << "\treplay_initial_time_s\treplay_initial_x_km\treplay_initial_y_km\treplay_initial_z_km"
			        << "\treplay_initial_vx_km_s\treplay_initial_vy_km_s\treplay_initial_vz_km_s"
			        << "\trng_state_before_initial_conditions\trng_state_before_simulation\n";
			summary << std::scientific << std::setprecision(17);
			for(const auto& rec : records)
			{
				const std::string key = std::to_string(rec.rank) + ":" + std::to_string(rec.trajectory_id);
				const auto replay_it = replay_by_key.find(key);
				const bool traced = replay_it != replay_by_key.end();
				const TrajectoryReplayRecord* replay = traced ? replay_it->second : nullptr;
				const double t_first = nan_if_missing(rec.t_first_unbinding_scatter);
				const double t_final = rec.event_observed ? nan_if_missing(rec.t_final_unbinding_scatter)
				                                          : std::numeric_limits<double>::quiet_NaN();
				const double t_escape = rec.event_observed ? nan_if_missing(rec.t_boundary_escape)
				                                           : std::numeric_limits<double>::quiet_NaN();
				const double t_censor = rec.event_observed ? std::numeric_limits<double>::quiet_NaN()
				                                           : nan_if_missing(rec.t_termination);
				const double lifetime_first = std::isfinite(t_first) ? t_first - rec.t_capture
				                                                        : std::numeric_limits<double>::quiet_NaN();
				summary << diagnostic_run_id << '\t' << rec.rank << '\t' << rec.trajectory_id << '\t'
				        << rec.rank << '\t' << rec.trajectory_id << '\t' << Diagnostic_Status(rec) << '\t'
				        << Termination_Reason_Key(rec.termination_reason) << '\t' << (rec.event_observed ? 1 : 0) << '\t'
				        << rec.t_capture << '\t' << t_first << '\t' << t_final << '\t' << t_escape << '\t' << t_censor << '\t'
				        << lifetime_first << '\t'
				        << (rec.event_observed ? rec.lifetime_unbinding : std::numeric_limits<double>::quiet_NaN()) << '\t'
				        << (rec.event_observed ? rec.lifetime_boundary : std::numeric_limits<double>::quiet_NaN()) << '\t'
				        << rec.r_first_negative_km / R_SUN_KM << '\t' << rec.E_first_negative_eV << '\t'
				        << nan_if_missing(rec.r_first_unbinding_km) / R_SUN_KM << '\t' << rec.E_first_unbinding_eV << '\t'
				        << (rec.event_observed ? rec.r_final_unbinding_km / R_SUN_KM : std::numeric_limits<double>::quiet_NaN()) << '\t'
				        << (rec.event_observed ? rec.E_final_unbinding_eV : std::numeric_limits<double>::quiet_NaN()) << '\t'
				        << (rec.event_observed ? rec.r_boundary_escape_km / R_SUN_KM : std::numeric_limits<double>::quiet_NaN()) << '\t'
				        << (rec.event_observed ? rec.vr_boundary_escape_km_s : std::numeric_limits<double>::quiet_NaN()) << '\t'
				        << (rec.event_observed ? rec.E_boundary_escape_eV : std::numeric_limits<double>::quiet_NaN()) << '\t'
				        << rec.number_of_scatterings << '\t' << rec.number_of_bound_to_unbound << '\t'
				        << rec.number_of_recaptures << '\t' << rec.min_energy_after_capture_eV << '\t'
				        << rec.max_radius_after_capture_km / R_SUN_KM << '\t'
				        << rec.time_inside_sun_after_capture_sec << '\t' << rec.time_outside_sun_after_capture_sec << '\t'
				        << rec.number_of_integrator_steps_after_capture << '\t' << rec.max_free_energy_drift_eV << '\t'
				        << rec.max_free_energy_drift_rel << '\t' << (traced ? 1 : 0) << '\t'
				        << (traced ? replay->initial_time_s : std::numeric_limits<double>::quiet_NaN());
				for(size_t component = 0; component < 3; component++)
					summary << '\t' << (traced ? replay->initial_position_km[component] : std::numeric_limits<double>::quiet_NaN());
				for(size_t component = 0; component < 3; component++)
					summary << '\t' << (traced ? replay->initial_velocity_km_s[component] : std::numeric_limits<double>::quiet_NaN());
				summary << '\t' << (traced ? replay->rng_state_before_initial_conditions : std::string())
				        << '\t' << (traced ? replay->rng_state_before_simulation : std::string()) << '\n';
			}
		}

		{
			std::ofstream events(output_dir + "/trajectory_events.tsv");
			events << "run_id\trank\ttrajectory_id\tevent_index\tscatter_index\tstep_index\tevent_type\tt_s"
			       << "\tr_Rsun\tvr_km_s\tspeed_km_s\tx\ty\tz\tvx\tvy\tvz\tenergy_eV\tangular_momentum"
			       << "\tis_bound\tinside_sun\tcandidate_active\ttarget_species\tballistic_energy_drift_eV\n";
			events << std::scientific << std::setprecision(17);
			for(const auto& event : trajectory_diagnostic_events)
			{
				events << diagnostic_run_id << '\t' << event.rank << '\t' << event.trajectory_id << '\t'
				       << event.event_index << '\t' << event.scatter_index << '\t' << event.step_index << '\t'
				       << TrajectoryDiagnosticEventTypeKey(event.event_type) << '\t' << event.t_s << '\t'
				       << event.r_km / R_SUN_KM << '\t' << event.vr_km_s << '\t' << event.speed_km_s;
				for(double value : event.position_km)
					events << '\t' << value;
				for(double value : event.velocity_km_s)
					events << '\t' << value;
				events << '\t' << event.energy_eV << '\t' << event.angular_momentum_km2_s
				       << '\t' << event.is_bound << '\t' << event.inside_sun << '\t' << event.candidate_active
				       << '\t' << event.target_species << '\t' << event.ballistic_energy_drift_eV << '\n';
			}
		}

		if(!unique_keys || !time_invariants || !legacy_reconciliation || !escape_radius_invariant
		   || !residence_time_invariant || !event_sequence_invariant || !event_count_invariant
		   || !trace_selection_invariant || !replay_state_invariant)
			std::cerr << "Warning in Write_Output_Files(): trajectory diagnostic invariant check failed; inspect run_metadata.json" << std::endl;
	}


}

double Simulation_Data::Free_Ratio() const
{
	return (number_of_trajectories > 0) ? 1.0 * number_of_free_particles / number_of_trajectories : 0.0;
}
double Simulation_Data::Capture_Ratio() const
{
	return (number_of_trajectories > 0) ? 1.0 * number_of_captured_particles / number_of_trajectories : 0.0;
}
double Simulation_Data::Reflection_Ratio(int isoreflection_ring) const
{
	if(isoreflection_ring < 0)
		return (number_of_trajectories > 0) ? 1.0 * number_of_reflected_particles / number_of_trajectories : 0.0;
	else
		return (number_of_trajectories > 0) ? 1.0 * data[isoreflection_ring].size() / number_of_trajectories : 0.0;
}

unsigned long int Simulation_Data::Valid_Trajectories() const
{
	return number_of_captured_particles + number_of_completed_outward_escapes;
}

double Simulation_Data::Free_Ratio_Valid() const
{
	const unsigned long int valid_trajectories = Valid_Trajectories();
	return (valid_trajectories > 0) ? 1.0 * number_of_free_particles / valid_trajectories : 0.0;
}

double Simulation_Data::Capture_Ratio_Valid() const
{
	const unsigned long int valid_trajectories = Valid_Trajectories();
	return (valid_trajectories > 0) ? 1.0 * number_of_captured_particles / valid_trajectories : 0.0;
}

double Simulation_Data::Reflection_Ratio_Valid(int isoreflection_ring) const
{
	const unsigned long int valid_trajectories = Valid_Trajectories();
	if(valid_trajectories == 0)
		return 0.0;
	if(isoreflection_ring < 0)
		return 1.0 * number_of_reflected_particles / valid_trajectories;
	else
		return 1.0 * data[isoreflection_ring].size() / valid_trajectories;
}

double Simulation_Data::Numerical_Failure_Ratio() const
{
	return (number_of_trajectories > 0) ? 1.0 * number_of_numerical_failures / number_of_trajectories : 0.0;
}

double Simulation_Data::Minimum_Speed() const
{
	return KDE_boundary_correction_factor * minimum_speed_threshold;
}

double Simulation_Data::Lowest_Speed(unsigned int iso_ring) const
{
	return (*std::min_element(data[iso_ring].begin(), data[iso_ring].end())).value;
}

double Simulation_Data::Highest_Speed(unsigned int iso_ring) const
{
	return (*std::max_element(data[iso_ring].begin(), data[iso_ring].end())).value;
}

void Simulation_Data::Print_Capture_Mode_Summary(unsigned int mpi_rank)
{
	if(mpi_rank == 0)
	{
		const BinomialRateEstimate raw = Estimate_Binomial_Rate(number_of_trajectories, number_of_captured_particles);
		const BinomialRateEstimate valid = Estimate_Binomial_Rate(Valid_Trajectories(), number_of_captured_particles);

		std::cout << SEPARATOR
		          << "CAPTURE MODE summary" << std::endl
		          << std::endl
		          << "Termination condition:\t\tpost-scatter E < 0" << std::endl
		          << "File output:\t\t\tdisabled" << std::endl
		          << "Simulated trajectories:\t\t" << number_of_trajectories << std::endl
		          << "Capture-classified trajectories:\t" << Valid_Trajectories() << std::endl
		          << "Unresolved non-captures:\t\t" << (number_of_trajectories - Valid_Trajectories()) << std::endl
		          << "Captured count:\t\t\t" << number_of_captured_particles << std::endl
		          << "Capture rate raw:\t\t" << std::fixed << std::setprecision(8) << raw.rate << std::endl
		          << "Capture rate valid:\t\t" << std::fixed << std::setprecision(8) << valid.rate << std::endl
		          << "Capture rate raw 95% CI:\t[" << std::fixed << std::setprecision(8) << raw.ci_lower << ", " << raw.ci_upper << "]" << std::endl
		          << "Capture rate valid 95% CI:\t[" << std::fixed << std::setprecision(8) << valid.ci_lower << ", " << valid.ci_upper << "]" << std::endl
		          << "Numerical failure count:\t" << number_of_numerical_failures << std::endl
		          << "Computational truncations:\t" << number_of_computational_truncations << std::endl
		          << "Initial shift failures:\t\t" << number_of_initial_shift_failures << std::endl
		          << "Final reflection shift failures:\t" << number_of_final_reflection_shift_failures << std::endl
		          << "Numerical failure rate:\t\t" << std::fixed << std::setprecision(8) << Numerical_Failure_Ratio() << std::endl;

		if(early_stopped)
			std::cout << "*** EARLY STOP: " << Stop_Reason_Display(early_stop_reason) << " ***" << std::endl;
		if(Numerical_Failure_Ratio() > NUMERICAL_FAILURE_WARNING_FRACTION)
			std::cout << "*** WARNING: numerical failure rate exceeded "
			          << NUMERICAL_FAILURE_WARNING_FRACTION << " ***" << std::endl;

		const double captured_particle_rate = (computing_time > 0.0) ? number_of_captured_particles / computing_time : 0.0;
		std::cout << "Captured particle rate [1/s]:\t" << libphysica::Round(captured_particle_rate) << std::endl
		          << "Simulation time [s]:\t\t" << std::fixed << std::setprecision(6) << computing_time << std::endl
		          << "Simulation time:\t\t" << libphysica::Time_Display(computing_time) << std::endl
		          << SEPARATOR << std::endl;
	}
}

void Simulation_Data::Print_Summary(unsigned int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR
				  << "Simulation data summary" << std::endl
				  << std::endl
				  << "Results:" << std::endl
				  << "Simulated trajectories:\t\t" << number_of_trajectories << std::endl
				  << "Capture-classified trajectories:\t" << Valid_Trajectories() << std::endl
				  << "Unresolved non-captures:\t\t" << (number_of_trajectories - Valid_Trajectories()) << std::endl
				  << "Average # of scatterings:\t" << libphysica::Round(average_number_of_scatterings) << std::endl
				  << "Free particles raw [%]:\t\t" << libphysica::Round(100.0 * Free_Ratio()) << std::endl
				  << "Reflected particles raw [%]:\t" << libphysica::Round(100.0 * Reflection_Ratio()) << std::endl
				  << "Captured particles raw [%]:\t" << libphysica::Round(100.0 * Capture_Ratio()) << std::endl
				  << "Free particles valid [%]:\t" << libphysica::Round(100.0 * Free_Ratio_Valid()) << std::endl
				  << "Reflected particles valid [%]:\t" << libphysica::Round(100.0 * Reflection_Ratio_Valid()) << std::endl
				  << "Captured particles valid [%]:\t" << libphysica::Round(100.0 * Capture_Ratio_Valid()) << std::endl
				  << "Captured count:\t\t\t" << number_of_captured_particles << std::endl
				  << "Normal-mode MPI sync interval:\t" << normal_mode_mpi_sync_interval << std::endl
				  << "MPI sync rounds:\t\t" << mpi_sync_rounds << std::endl
				  << "Final MPI round trajectories:\t" << final_mpi_round_trajectories << std::endl
				  << "MPI tail trajectories:\t\t" << mpi_tail_trajectories << std::endl
				  << "Capture target overshoot:\t" << capture_target_overshoot << std::endl
				  << "Total # of scatterings:\t" << total_number_of_scatterings << std::endl
				  << "Numerical failure count:\t" << number_of_numerical_failures << std::endl
				  << "Computational truncations:\t" << number_of_computational_truncations << std::endl
				  << "Initial shift failures:\t\t" << number_of_initial_shift_failures << std::endl
				  << "Final reflection shift failures:\t" << number_of_final_reflection_shift_failures << std::endl
				  << "Numerical failure rate:\t\t" << std::fixed << std::setprecision(6) << Numerical_Failure_Ratio() << std::endl
				  << "Complete evaporation count:\t" << number_of_complete_evaporation_particles << std::endl
				  << "Censored captured count:\t" << number_of_censored_captured_particles << std::endl
				  << "Invalid survival count:\t\t" << number_of_invalid_survival_captured_particles << std::endl;

		// Raw and classified-sample capture-rate intervals.
		{
			const BinomialRateEstimate raw = Estimate_Binomial_Rate(number_of_trajectories, number_of_captured_particles);
			const BinomialRateEstimate valid = Estimate_Binomial_Rate(Valid_Trajectories(), number_of_captured_particles);
			std::cout << "Capture rate raw error (1σ):\t" << std::fixed << std::setprecision(6) << raw.standard_error << std::endl
			          << "Capture rate raw 95% CI:\t[" << raw.ci_lower << ", " << raw.ci_upper << "]" << std::endl
			          << "Capture rate valid error (1σ):\t" << valid.standard_error << std::endl
			          << "Capture rate valid 95% CI:\t[" << valid.ci_lower << ", " << valid.ci_upper << "]" << std::endl;
		}

		if(early_stopped)
			std::cout << "*** EARLY STOP: " << Stop_Reason_Display(early_stop_reason) << " ***" << std::endl;
		if(Numerical_Failure_Ratio() > NUMERICAL_FAILURE_WARNING_FRACTION)
			std::cout << "*** WARNING: numerical failure rate exceeded "
			          << NUMERICAL_FAILURE_WARNING_FRACTION << " ***" << std::endl;

		// Median for observed unbinding events only; censored records belong in survival analysis.
		std::vector<double> observed_unbinding_lifetimes;
		if(evaporation_diagnostics_enabled)
		{
			for(const auto& rec : evaporation_records)
			{
				if(rec.survival_valid && rec.event_observed && Has_Positive_Evaporation_Time(rec.lifetime_unbinding))
					observed_unbinding_lifetimes.push_back(rec.lifetime_unbinding);
			}
		}
		else
		{
			for(const auto& event : compact_evaporation_events)
			{
				if(Has_Positive_Evaporation_Time(event.lifetime_unbinding))
					observed_unbinding_lifetimes.push_back(event.lifetime_unbinding);
			}
		}
		if(!observed_unbinding_lifetimes.empty())
		{
			std::sort(observed_unbinding_lifetimes.begin(), observed_unbinding_lifetimes.end());
			double median;
			size_t n = observed_unbinding_lifetimes.size();
			if(n % 2 == 0)
				median = 0.5 * (observed_unbinding_lifetimes[n/2 - 1] + observed_unbinding_lifetimes[n/2]);
			else
				median = observed_unbinding_lifetimes[n/2];
			std::cout << "Observed unbinding lifetime median [s]:\t" << std::scientific << std::setprecision(4) << median << " (" << observed_unbinding_lifetimes.size() << " events)" << std::endl;
		}
		else
		{
			std::cout << "Observed unbinding lifetime median:\tN/A (no positive observed unbinding events)" << std::endl;
		}

		std::cout << std::endl
				  << "Trajectory rate [1/s]:\t\t" << libphysica::Round(1.0 * number_of_trajectories / computing_time) << std::endl
				  << "Captured particle rate [1/s]:\t" << libphysica::Round(1.0 * number_of_captured_particles / computing_time) << std::endl
				  << "Simulation time [s]:\t\t" << std::fixed << std::setprecision(6) << computing_time << std::endl
				  << "Simulation time:\t\t" << libphysica::Time_Display(computing_time) << std::endl;


		std::cout << SEPARATOR << std::endl;
	}
}

}	// namespace DaMaSCUS_SUN

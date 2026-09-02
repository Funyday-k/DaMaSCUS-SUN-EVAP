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
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <thread>
#include <type_traits>
#include <unistd.h>
#include <utility>
#include <vector>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/Astronomy.hpp"
#include "MPI_Utilities.hpp"
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

constexpr double NUMERICAL_FAILURE_WARNING_FRACTION = 1.0e-4;
constexpr double INITIAL_SHIFT_FAILURE_ABORT_FRACTION = 1.0e-2;

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

int NumericalFailureDetail_Index(
	TrajectoryNumericalFailureDetail detail)
{
	int index = static_cast<int>(detail);
	if(index < 0
	   || index >= TRAJECTORY_NUMERICAL_FAILURE_DETAIL_COUNT)
		return static_cast<int>(
		    TrajectoryNumericalFailureDetail::None);
	return index;
}

bool Completed_Outward_Escape(TrajectoryTerminationReason reason)
{
	return reason == TrajectoryTerminationReason::OutwardEscape;
}

bool Is_Computational_Truncation(TrajectoryTerminationReason reason)
{
	// A wall-time stop is an intentional observation cutoff: its accepted
	// captured path remains usable for the residence bincount, and the scheduler
	// must replace it with a new trajectory rather than count it as invalid work.
	return reason == TrajectoryTerminationReason::MaxFreeSteps
	    || reason == TrajectoryTerminationReason::MaxScatterings;
}

bool Is_Numerical_Termination(TrajectoryTerminationReason reason)
{
	return TrajectoryTerminationInvalidatesResidenceBincount(reason);
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
	// Wall-time termination is a valid computational censor: its accepted
	// residence prefix is retained, but no compact evaporation-time event is
	// emitted. A radial-domain removal remains a separate physical exclusion.
	rec.censored = survival_valid
	            && !event_observed
	            && bincount.termination_reason
	               == TrajectoryTerminationReason::WallTimeLimit;
	rec.outer_domain_removed =
	    bincount.termination_reason == TrajectoryTerminationReason::OuterDomainRemoval;
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
	rec.number_of_bound_exterior_arcs = bincount.number_of_bound_exterior_arcs;
	rec.first_bound_exit_kepler_period_sec = bincount.first_bound_exit_kepler_period_sec;
	rec.last_bound_exit_kepler_period_sec = bincount.last_bound_exit_kepler_period_sec;
	rec.max_bound_exit_kepler_period_sec = bincount.max_bound_exit_kepler_period_sec;
	rec.first_bound_exit_exterior_time_sec = bincount.first_bound_exit_exterior_time_sec;
	rec.last_bound_exit_exterior_time_sec = bincount.last_bound_exit_exterior_time_sec;
	rec.max_bound_exit_exterior_time_sec = bincount.max_bound_exit_exterior_time_sec;
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
		case TrajectoryTerminationReason::OuterDomainRemoval: return "outer_domain_removal";
		case TrajectoryTerminationReason::Unknown:
		default: return "unknown";
	}
}

const char* Invalid_Trajectory_Stage_Key(InvalidTrajectoryStage stage)
{
	switch(stage)
	{
		case InvalidTrajectoryStage::InitialShift: return "initial_shift";
		case InvalidTrajectoryStage::Propagation: return "propagation";
		case InvalidTrajectoryStage::InvalidCapturedSurvival: return "invalid_captured_survival";
		default: return "unknown";
	}
}

std::string Serialize_PRNG_State(const std::mt19937& state)
{
	std::ostringstream stream;
	stream << state;
	return stream.str();
}

std::string Encode_PRNG_State(std::string state)
{
	std::replace(state.begin(), state.end(), ' ', ',');
	return state;
}

std::string Decode_PRNG_State(std::string state)
{
	std::replace(state.begin(), state.end(), ',', ' ');
	return state;
}

void Print_Termination_Reason_Summary(
	const std::array<unsigned long int, TRAJECTORY_TERMINATION_REASON_COUNT>& captured,
	const std::array<unsigned long int, TRAJECTORY_TERMINATION_REASON_COUNT>& uncaptured)
{
	std::cout << "Termination reasons (captured / uncaptured):" << std::endl;
	for(int reason_index = 0;
	    reason_index < TRAJECTORY_TERMINATION_REASON_COUNT;
	    reason_index++)
	{
		const size_t index = static_cast<size_t>(reason_index);
		if(captured[index] == 0 && uncaptured[index] == 0)
			continue;
		const auto reason = static_cast<TrajectoryTerminationReason>(reason_index);
		std::cout << "\t" << Termination_Reason_Key(reason) << ":\t"
		          << captured[index] << " / " << uncaptured[index] << std::endl;
	}
}

void Print_Numerical_Failure_Detail_Summary(
	const std::array<unsigned long int, TRAJECTORY_NUMERICAL_FAILURE_DETAIL_COUNT>& captured,
	const std::array<unsigned long int, TRAJECTORY_NUMERICAL_FAILURE_DETAIL_COUNT>& uncaptured)
{
	bool any_detail = false;
	for(int detail_index = 1;
	    detail_index < TRAJECTORY_NUMERICAL_FAILURE_DETAIL_COUNT;
	    detail_index++)
	{
		const size_t index = static_cast<size_t>(detail_index);
		any_detail =
		    any_detail
		    || captured[index] != 0
		    || uncaptured[index] != 0;
	}
	if(!any_detail)
		return;
	std::cout
	    << "Numerical failure details (captured / uncaptured):"
	    << std::endl;
	for(int detail_index = 1;
	    detail_index < TRAJECTORY_NUMERICAL_FAILURE_DETAIL_COUNT;
	    detail_index++)
	{
		const size_t index = static_cast<size_t>(detail_index);
		if(captured[index] == 0 && uncaptured[index] == 0)
			continue;
		const auto detail =
		    static_cast<TrajectoryNumericalFailureDetail>(
		        detail_index);
		std::cout
		    << "\t"
		    << TrajectoryNumericalFailureDetailKey(detail)
		    << ":\t"
		    << captured[index]
		    << " / "
		    << uncaptured[index]
		    << std::endl;
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

std::size_t Residence_Jackknife_Block(
	uint64_t base_seed,
	int rank,
	uint64_t trajectory_id)
{
	uint64_t value =
	    base_seed
	    ^ (static_cast<uint64_t>(
	           static_cast<uint32_t>(rank))
	       << 32)
	    ^ trajectory_id;
	value += 0x9e3779b97f4a7c15ULL;
	value =
	    (value ^ (value >> 30))
	    * 0xbf58476d1ce4e5b9ULL;
	value =
	    (value ^ (value >> 27))
	    * 0x94d049bb133111ebULL;
	value ^= value >> 31;
	return static_cast<std::size_t>(
	    value % RESIDENCE_JACKKNIFE_BLOCKS);
}

const char* Diagnostic_Status(const EvaporationRecord& rec)
{
	if(rec.event_observed)
		return "escaped";
	if(rec.termination_reason == TrajectoryTerminationReason::OuterDomainRemoval)
		return "outer_domain_removed";
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
	event.r_capture_rsun = rec.r_first_negative_km / R_SUN_KM;
	event.E_capture_eV = rec.E_first_negative_eV;
	event.dE_capture_eV = rec.dE_first_negative_from_prev_eV;
	event.number_of_bound_exterior_arcs = rec.number_of_bound_exterior_arcs;
	event.first_bound_exit_kepler_period_sec = rec.first_bound_exit_kepler_period_sec;
	event.last_bound_exit_kepler_period_sec = rec.last_bound_exit_kepler_period_sec;
	event.max_bound_exit_kepler_period_sec = rec.max_bound_exit_kepler_period_sec;
	event.first_bound_exit_exterior_time_sec = rec.first_bound_exit_exterior_time_sec;
	event.last_bound_exit_exterior_time_sec = rec.last_bound_exit_exterior_time_sec;
	event.max_bound_exit_exterior_time_sec = rec.max_bound_exit_exterior_time_sec;
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
	file << "# capture_rate_raw = " << std::fixed << std::setprecision(8) << raw.rate << "\n";
	file << "# capture_rate_raw_err = " << std::fixed << std::setprecision(8) << raw.standard_error << "\n";
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
	file << "# format_version = 5\n";
	file << "# DM_mass_GeV = " << std::scientific << std::setprecision(6) << mass_gev << "\n";
	file << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(6) << sigma_cm2 << "\n";
	file << "# sorted_by = lifetime_unbinding_sec rank trajectory_id\n";
	file << "# rank trajectory_id lifetime_unbinding_sec r_capture_Rsun E_capture_eV dE_capture_eV"
	     << " n_bound_exterior_arcs P_kepler_first_bound_exit_sec"
	     << " P_kepler_last_bound_exit_sec P_kepler_max_bound_exit_sec"
	     << " t_exterior_first_bound_exit_sec t_exterior_last_bound_exit_sec"
	     << " t_exterior_max_bound_exit_sec\n";
}

void Write_Evaporation_Log_Event(std::ostream& file, const CompactEvaporationEvent& event)
{
	file << event.rank << "\t" << event.trajectory_id << "\t" << std::scientific << std::setprecision(10)
	     << event.lifetime_unbinding << "\t" << event.r_capture_rsun
	     << "\t" << event.E_capture_eV << "\t" << event.dE_capture_eV
	     << "\t" << event.number_of_bound_exterior_arcs
	     << "\t" << event.first_bound_exit_kepler_period_sec
	     << "\t" << event.last_bound_exit_kepler_period_sec
	     << "\t" << event.max_bound_exit_kepler_period_sec
	     << "\t" << event.first_bound_exit_exterior_time_sec
	     << "\t" << event.last_bound_exit_exterior_time_sec
	     << "\t" << event.max_bound_exit_exterior_time_sec << "\n";
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
  maximum_trajectories(
      max_trajectories == 0
      ? std::numeric_limits<uint64_t>::max()
      : static_cast<uint64_t>(max_trajectories)),
  number_of_trajectories(0), number_of_free_particles(0), number_of_reflected_particles(0), number_of_captured_particles(0),
  number_of_completed_outward_escapes(0),
  number_of_complete_evaporation_particles(0), number_of_residence_samples(0),
  number_of_censored_captured_particles(0),
  number_of_outer_domain_removed_particles(0),
  number_of_invalid_survival_captured_particles(0),
  number_of_initial_shift_failures(0), number_of_final_reflection_shift_failures(0), number_of_numerical_failures(0),
  number_of_computational_truncations(0),
  total_number_of_scatterings(0), average_number_of_scatterings(0.0),
  mpi_scheduler_work_claims(0), mpi_scheduler_peak_in_flight(0),
  capture_target_overshoot(0),
  computing_time(0.0), early_stopped(false),
  early_stop_reason(SimulationStopReason::None),
  mpi_rank(0), mpi_processes(1), isoreflection_rings(iso_rings), minimum_speed_threshold(u_min),
  number_of_data_points(std::vector<unsigned long int>(iso_rings, 0)),
  data(iso_rings, std::vector<libphysica::DataPoint>())
{
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
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

void Simulation_Data::Generate_Data(obscura::DM_Particle& DM, Celestial_Model& solar_model, obscura::DM_Distribution& halo_model, SnapshotConfig snapshot_cfg, unsigned int fixed_seed, bool capture_mode)
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
	unsigned long int local_total = 0;
	bool initial_shift_failure_warning_emitted = false;

    // One grid is shared by simulator, snapshots, aggregation,
    // jackknife storage, and final output.
    bincount_radial_grid = Bincount_Radial_Grid(
        In_Units(solar_model.Radius(), km));

    const std::size_t radial_bin_count =
        bincount_radial_grid.Bin_Count();

    captured_dt_hist.assign(radial_bin_count, 0.0);
    captured_v2dt_hist.assign(radial_bin_count, 0.0);
    captured_dt_sq_hist.assign(radial_bin_count, 0.0);
    captured_v2dt_sq_hist.assign(radial_bin_count, 0.0);

    residence_jackknife_block_dt_hist.assign(
        RESIDENCE_JACKKNIFE_BLOCKS * radial_bin_count,
        0.0);
    residence_jackknife_block_v2dt_hist.assign(
        RESIDENCE_JACKKNIFE_BLOCKS * radial_bin_count,
        0.0);

    // Configure the simulator with that same runtime grid.
    Trajectory_Simulator simulator(
        solar_model,
        bincount_radial_grid,
        maximum_free_time_steps,
        maximum_number_of_scatterings,
        initial_and_final_radius);
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
                    snapshot_state.reset(new SnapshotSharedState(
                        bincount_radial_grid.Bin_Count()));
                    snapshot_state->Initialize(snapshot_run_id, mpi_rank);
                    snapshot_recorder.reset(
                        new SnapshotRecorder(*snapshot_state));
                    snapshot_heartbeat.reset(new SnapshotHeartbeat(
                            *snapshot_state,
                            mpi_rank,
                            mpi_processes,
                            snapshot_run_id,
                            snapshot_root,
                            rank_snapshot_dir,
                            snapshot_interval,
                            snapshot_mass_gev,
                            snapshot_sigma_cm2,
                            bincount_radial_grid));
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

	unsigned long int global_target_samples = 0;
	const std::vector<double> progress_milestones = {0.0, 0.01, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 1.0};
	size_t next_progress_milestone = 0;
	bool progress_line_printed = false;
	unsigned long int last_progress_line_samples = 0;
	auto print_progress_update = [&](unsigned long int accepted_samples, bool force)
	{
		if(mpi_rank != 0)
			return;

		const double denominator = static_cast<double>(std::max(1u, requested_captured_particles));
		const double progress = std::min(1.0, static_cast<double>(accepted_samples) / denominator);

		bool should_print = force;
		while(next_progress_milestone < progress_milestones.size()
		      && progress + 1.0e-12 >= progress_milestones[next_progress_milestone])
		{
			should_print = true;
			next_progress_milestone++;
		}
		if(!should_print)
			return;
		if(force && progress_line_printed && accepted_samples == last_progress_line_samples)
			return;

		const double time_elapsed = elapsed_since_start();
		const double sample_rate = (time_elapsed > 0.0) ? static_cast<double>(accepted_samples) / time_elapsed : 0.0;
		libphysica::Print_Progress_Bar(progress, 0, 44, time_elapsed);
		std::cout << (capture_mode ? " captured_particles=" : " valid_evaporation_samples=")
		          << accepted_samples << "/" << requested_captured_particles
		          << (capture_mode ? " captured_particle_rate[1/s]=" : " valid_evaporation_sample_rate[1/s]=")
		          << libphysica::Round(sample_rate)
		          << std::endl;
		progress_line_printed = true;
		last_progress_line_samples = accepted_samples;
	};
	print_progress_update(global_target_samples, true);

	auto simulate_one_trajectory = [&]() -> MPIWorkOutcome
	{
			const unsigned long int trajectory_id = local_total + 1;
			const bool trace_selected = trajectory_diagnostic_config.events_enabled
			                         && Diagnostic_Trace_Selected(trajectory_diagnostic_config.trace_seed,
			                                                      mpi_rank, trajectory_id,
			                                                      trajectory_diagnostic_config.trace_rate);
			// Copying the engine is cheap compared with serializing it. Keep
			// exact pre-trajectory states in memory, and serialize only the rare
			// invalid trajectories (or explicitly selected diagnostic traces).
			const std::mt19937 replay_rng_before_initial_conditions = simulator.PRNG;
			const std::string rng_state_before_initial_conditions = trace_selected
			                                                       ? simulator.Serialize_PRNG_State()
			                                                       : std::string();
			simulator.Enable_Diagnostic_Trace(trace_selected);
			Event IC = Initial_Conditions(halo_model, solar_model, simulator.PRNG);
			const bool initial_shift_ok = Hyperbolic_Kepler_Shift(IC, initial_and_final_radius);
			const std::mt19937 replay_rng_before_simulation = simulator.PRNG;
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
			const size_t jackknife_block =
			    Residence_Jackknife_Block(
			        diagnostic_base_seed,
			        mpi_rank,
			        trajectory_id);
			bool invalid_recorded = false;
			auto record_invalid_trajectory = [&](InvalidTrajectoryStage stage)
			{
				if(capture_mode || invalid_recorded)
					return;
				InvalidTrajectoryRecord record;
				record.rank = mpi_rank;
				record.trajectory_id = trajectory_id;
				record.failure_stage = stage;
				record.termination_reason = trajectory.bincount.termination_reason;
				record.numerical_failure_detail =
				    trajectory.bincount.numerical_failure_detail;
				record.initial_shift_ok = initial_shift_ok;
				record.is_captured = trajectory.bincount.is_captured;
				record.survival_valid = trajectory.bincount.survival_valid;
				record.event_observed = trajectory.bincount.event_observed;
				record.number_of_scatterings = trajectory.number_of_scatterings;
				record.number_of_bound_to_unbound = trajectory.bincount.number_of_bound_to_unbound;
				record.number_of_recaptures = trajectory.bincount.number_of_recaptures;
				record.t_capture_s = trajectory.bincount.t_capture;
				record.t_termination_s = trajectory.bincount.t_termination;
				const double final_radius = trajectory.final_event.Radius();
				const double final_speed = trajectory.final_event.Speed();
				record.final_r_rsun = std::isfinite(final_radius)
				                    ? In_Units(final_radius, g_body_radius)
				                    : std::numeric_limits<double>::quiet_NaN();
				record.final_speed_km_s = std::isfinite(final_speed)
				                        ? In_Units(final_speed, km / sec)
				                        : std::numeric_limits<double>::quiet_NaN();
				if(std::isfinite(final_radius) && final_radius > 0.0)
				{
					const double final_radial_velocity =
					    trajectory.final_event.position.Dot(trajectory.final_event.velocity) / final_radius;
					record.final_vr_km_s = In_Units(final_radial_velocity, km / sec);
					const double final_escape_speed = solar_model.Local_Escape_Speed(final_radius);
					const double final_energy =
					    0.5 * DM.mass * (final_speed * final_speed - final_escape_speed * final_escape_speed);
					record.final_energy_eV = In_Units(final_energy, eV);
				}
				record.max_r_after_capture_rsun =
				    std::isfinite(trajectory.bincount.max_radius_after_capture_km)
				    ? trajectory.bincount.max_radius_after_capture_km / R_SUN_KM
				    : std::numeric_limits<double>::quiet_NaN();
				record.max_free_energy_drift_eV = trajectory.bincount.max_free_energy_drift_eV;
				record.max_free_energy_drift_rel = trajectory.bincount.max_free_energy_drift_rel;
				record.failure_energy_before_step_eV =
				    trajectory.bincount.failure_energy_before_step_eV;
				record.failure_energy_after_step_eV =
				    trajectory.bincount.failure_energy_after_step_eV;
				record.failure_energy_at_boundary_eV =
				    trajectory.bincount.failure_energy_at_boundary_eV;
				record.failure_reference_energy_eV =
				    trajectory.bincount.failure_reference_energy_eV;
				record.failure_boundary_vr_km_s =
				    trajectory.bincount.failure_boundary_vr_km_s;
				record.failure_attempted_step_s =
				    trajectory.bincount.failure_attempted_step_s;
				record.failure_accepted_step_s =
				    trajectory.bincount.failure_accepted_step_s;
				record.initial_time_s = In_Units(IC.time, sec);
				for(size_t component = 0; component < 3; component++)
				{
					record.initial_position_km[component] = In_Units(IC.position[component], km);
					record.initial_velocity_km_s[component] = In_Units(IC.velocity[component], km / sec);
				}
				record.rng_state_before_initial_conditions =
				    trace_selected ? rng_state_before_initial_conditions
				                   : Serialize_PRNG_State(replay_rng_before_initial_conditions);
				record.rng_state_before_simulation =
				    trace_selected ? rng_state_before_simulation
				                   : Serialize_PRNG_State(replay_rng_before_simulation);
				invalid_trajectory_records.push_back(std::move(record));
				invalid_recorded = true;
			};

			local_total++;
			number_of_trajectories++;
			jackknife_attempted_counts[jackknife_block]++;
			total_number_of_scatterings += static_cast<uint64_t>(trajectory.number_of_scatterings);
			const int termination_reason_index =
			    TerminationReason_Index(trajectory.bincount.termination_reason);
			if(trajectory.bincount.is_captured)
				captured_termination_reason_counts[static_cast<size_t>(termination_reason_index)]++;
			else
				uncaptured_termination_reason_counts[static_cast<size_t>(termination_reason_index)]++;
			const int failure_detail_index =
			    NumericalFailureDetail_Index(
			        trajectory.bincount.numerical_failure_detail);
			if(trajectory.bincount.is_captured)
				captured_numerical_failure_detail_counts[
				    static_cast<size_t>(failure_detail_index)]++;
			else
				uncaptured_numerical_failure_detail_counts[
				    static_cast<size_t>(failure_detail_index)]++;
			const bool completed_outward_escape = Completed_Outward_Escape(trajectory.bincount.termination_reason);
			if(initial_shift_ok && Is_Numerical_Termination(trajectory.bincount.termination_reason))
				number_of_numerical_failures++;
			if(Is_Computational_Truncation(trajectory.bincount.termination_reason))
				number_of_computational_truncations++;
			if(!initial_shift_ok
			   || Is_Numerical_Termination(
			       trajectory.bincount.termination_reason)
			   || Is_Computational_Truncation(
			       trajectory.bincount.termination_reason))
				jackknife_invalid_counts[jackknife_block]++;
			if(!initial_shift_ok)
				record_invalid_trajectory(InvalidTrajectoryStage::InitialShift);
			else if(Is_Numerical_Termination(trajectory.bincount.termination_reason)
			        || Is_Computational_Truncation(trajectory.bincount.termination_reason))
				record_invalid_trajectory(InvalidTrajectoryStage::Propagation);
			std::vector<SnapshotEvaporationProgressEntry> trajectory_snapshot_evaporation_events;
			const bool accepted_evaporation_sample =
			    !capture_mode
			    && trajectory.bincount.is_captured
			    && trajectory.bincount.survival_valid
			    && trajectory.bincount.event_observed
			    && trajectory.bincount.termination_reason
			       == TrajectoryTerminationReason::OutwardEscape;
			const bool accepted_residence_sample =
			    !capture_mode
			    && trajectory.bincount.is_captured
			    && !TrajectoryTerminationInvalidatesResidenceBincount(
			        trajectory.bincount.termination_reason);

			if(trajectory.bincount.is_captured)
			{
				number_of_captured_particles++;
				jackknife_captured_counts[jackknife_block]++;
				if(capture_mode)
				{
					// Capture itself is the accepted target in this mode.
				}
				else if(accepted_evaporation_sample)
				{
					number_of_complete_evaporation_particles++;
				}
				else if(trajectory.bincount.termination_reason
				        == TrajectoryTerminationReason::OuterDomainRemoval)
				{
					number_of_outer_domain_removed_particles++;
					jackknife_outer_domain_removed_counts[
					    jackknife_block]++;
				}
				else if(trajectory.bincount.termination_reason
				        == TrajectoryTerminationReason::WallTimeLimit)
				{
					// This trajectory contributes its accepted residence prefix but
					// has no observed evaporation time. It is intentionally replaced
					// by another work-queue claim.
					number_of_censored_captured_particles++;
				}
				else
				{
					number_of_invalid_survival_captured_particles++;
					record_invalid_trajectory(InvalidTrajectoryStage::InvalidCapturedSurvival);
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
                                          if(accepted_residence_sample)
                                          {
                                                  if(trajectory.bincount.dt_hist.size()
                                                         != captured_dt_hist.size()
                                                     || trajectory.bincount.v2dt_hist.size()
                                                         != captured_v2dt_hist.size())
                                                  {
                                                          throw std::runtime_error(
                                                              "trajectory bincount histogram size "
                                                              "does not match Simulation_Data radial grid");
                                                  }

                                                  number_of_residence_samples++;
                                                  jackknife_residence_sample_counts[
                                                      jackknife_block]++;

                                                  const std::size_t radial_bin_count =
                                                      captured_dt_hist.size();

                                                  for(std::size_t b = 0;
                                                      b < radial_bin_count;
                                                      b++)
                                                  {
                                                          captured_dt_hist[b] +=
                                                              trajectory.bincount.dt_hist[b];
                                                          captured_v2dt_hist[b] +=
                                                              trajectory.bincount.v2dt_hist[b];

                                                          captured_dt_sq_hist[b] +=
                                                              trajectory.bincount.dt_hist[b]
                                                              * trajectory.bincount.dt_hist[b];
                                                          captured_v2dt_sq_hist[b] +=
                                                              trajectory.bincount.v2dt_hist[b]
                                                              * trajectory.bincount.v2dt_hist[b];

                                                          const std::size_t block_offset =
                                                              jackknife_block
                                                              * radial_bin_count
                                                              + b;

                                                          residence_jackknife_block_dt_hist[
                                                              block_offset] +=
                                                              trajectory.bincount.dt_hist[b];
                                                          residence_jackknife_block_v2dt_hist[
                                                              block_offset] +=
                                                              trajectory.bincount.v2dt_hist[b];
                                                  }
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
				jackknife_completed_escape_counts[
				    jackknife_block]++;
				if(!capture_mode)
				{
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
							}
					}
				}
			}

			if(snapshot_state)
			{
				const bool count_as_residence_sample =
				    accepted_residence_sample;
				const bool physically_classified_uncaptured =
				    !capture_mode
				    && completed_outward_escape
				    && !trajectory.bincount.is_captured;
				snapshot_state->RecordCompletedTrajectory(
					trajectory.bincount,
					count_as_residence_sample,
					physically_classified_uncaptured,
					trajectory_snapshot_evaporation_events);
			}

			MPIWorkOutcome outcome;
			outcome.accepted_sample =
			    capture_mode
			    ? trajectory.bincount.is_captured
			    : accepted_evaporation_sample;
			outcome.initial_shift_failure = !initial_shift_ok;
			outcome.numerical_failure =
			    !initial_shift_ok
			    || (initial_shift_ok
			        && Is_Numerical_Termination(
			           trajectory.bincount.termination_reason));
			outcome.computational_truncation =
			    Is_Computational_Truncation(
			        trajectory.bincount.termination_reason);
			return outcome;
		};

	MPIWorkQueue work_queue(
	    requested_captured_particles,
	    maximum_trajectories,
	    INITIAL_SHIFT_FAILURE_ABORT_FRACTION,
	    MPI_COMM_WORLD);
	while(true)
	{
		MPIWorkQueueState observed;
		const MPIWorkClaimResult claim =
		    work_queue.TryClaim(&observed);
		global_target_samples =
		    static_cast<unsigned long int>(
		        observed.accepted_samples);
		if(claim == MPIWorkClaimResult::Stop)
			break;
		if(claim == MPIWorkClaimResult::Wait)
		{
			// Only exact-target tail pressure or an exhausted global trajectory
			// budget can leave a rank without a claim. A short local sleep keeps
			// polling from consuming a full CPU while active ranks finish.
			std::this_thread::sleep_for(
			    std::chrono::milliseconds(1));
			continue;
		}

		const MPIWorkQueueState completed =
		    work_queue.Complete(simulate_one_trajectory());
		global_target_samples =
		    static_cast<unsigned long int>(
		        completed.accepted_samples);
		print_progress_update(global_target_samples, false);
	}

	const MPIWorkQueueState final_work_state =
	    work_queue.Finalize();
	global_target_samples =
	    static_cast<unsigned long int>(
	        final_work_state.accepted_samples);
	mpi_scheduler_work_claims = final_work_state.work_claims;
	mpi_scheduler_peak_in_flight =
	    final_work_state.peak_in_flight;
	switch(final_work_state.stop_reason)
	{
		case MPIWorkStopReason::MaxTrajectoriesReached:
			early_stop_reason =
			    SimulationStopReason::MaxTrajectoriesReached;
			break;
		case MPIWorkStopReason::
		    InitialShiftFailureFractionExceeded:
			early_stop_reason =
			    SimulationStopReason::
			        InitialShiftFailureFractionExceeded;
			break;
		case MPIWorkStopReason::None:
		default:
			break;
	}

	if(final_work_state.completed_trajectories > 0)
	{
		const double attempted = static_cast<double>(
		    final_work_state.completed_trajectories);
		const double initial_shift_failure_fraction =
		    static_cast<double>(
		        final_work_state.initial_shift_failures)
		    / attempted;
		if(mpi_rank == 0
		   && final_work_state.stop_reason
		      == MPIWorkStopReason::
		          InitialShiftFailureFractionExceeded)
		{
			std::cerr
			    << "Error in Generate_Data(): initial Kepler shift "
			    << "failure fraction "
			    << initial_shift_failure_fraction
			    << " exceeds abort threshold "
			    << INITIAL_SHIFT_FAILURE_ABORT_FRACTION
			    << ". Stopping this run to avoid biased capture "
			    << "statistics." << std::endl;
		}
		else if(mpi_rank == 0
		        && !initial_shift_failure_warning_emitted
		        && initial_shift_failure_fraction
		           > NUMERICAL_FAILURE_WARNING_FRACTION)
		{
			std::cerr
			    << "Warning in Generate_Data(): initial Kepler shift "
			    << "failure fraction "
			    << initial_shift_failure_fraction
			    << " exceeds warning threshold "
			    << NUMERICAL_FAILURE_WARNING_FRACTION
			    << "." << std::endl;
			initial_shift_failure_warning_emitted = true;
		}
	}
	capture_target_overshoot = (global_target_samples > requested_captured_particles)
	                         ? global_target_samples - requested_captured_particles
	                         : 0UL;

	if(global_target_samples < requested_captured_particles)
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

	print_progress_update(global_target_samples, global_target_samples > 0 || !early_stopped);
	if(mpi_rank == 0)
		std::cout << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	Perform_MPI_Reductions(capture_mode);
}

void Simulation_Data::Perform_MPI_Reductions(bool capture_mode)
{
	MPI_Trace_Point(mpi_rank, "enter Perform_MPI_Reductions");
	std::array<unsigned long int, 7> primary_counters = {{
	    number_of_trajectories,
	    number_of_captured_particles,
	    number_of_completed_outward_escapes,
	    number_of_initial_shift_failures,
	    number_of_final_reflection_shift_failures,
	    number_of_numerical_failures,
	    number_of_computational_truncations}};
	MPI_Trace_Point(mpi_rank, "before allreduce primary counters");
	MPI_Allreduce(MPI_IN_PLACE, primary_counters.data(), static_cast<int>(primary_counters.size()), MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	number_of_trajectories = primary_counters[0];
	number_of_captured_particles = primary_counters[1];
	number_of_completed_outward_escapes = primary_counters[2];
	number_of_initial_shift_failures = primary_counters[3];
	number_of_final_reflection_shift_failures = primary_counters[4];
	number_of_numerical_failures = primary_counters[5];
	number_of_computational_truncations = primary_counters[6];
	MPI_Trace_Point(mpi_rank, "before allreduce captured termination reasons");
	MPI_Allreduce(
	    MPI_IN_PLACE,
	    captured_termination_reason_counts.data(),
	    TRAJECTORY_TERMINATION_REASON_COUNT,
	    MPI_UNSIGNED_LONG,
	    MPI_SUM,
	    MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce uncaptured termination reasons");
	MPI_Allreduce(
	    MPI_IN_PLACE,
	    uncaptured_termination_reason_counts.data(),
	    TRAJECTORY_TERMINATION_REASON_COUNT,
	    MPI_UNSIGNED_LONG,
	    MPI_SUM,
	    MPI_COMM_WORLD);
	MPI_Trace_Point(
	    mpi_rank,
	    "before allreduce captured numerical failure details");
	MPI_Allreduce(
	    MPI_IN_PLACE,
	    captured_numerical_failure_detail_counts.data(),
	    TRAJECTORY_NUMERICAL_FAILURE_DETAIL_COUNT,
	    MPI_UNSIGNED_LONG,
	    MPI_SUM,
	    MPI_COMM_WORLD);
	MPI_Trace_Point(
	    mpi_rank,
	    "before allreduce uncaptured numerical failure details");
	MPI_Allreduce(
	    MPI_IN_PLACE,
	    uncaptured_numerical_failure_detail_counts.data(),
	    TRAJECTORY_NUMERICAL_FAILURE_DETAIL_COUNT,
	    MPI_UNSIGNED_LONG,
	    MPI_SUM,
	    MPI_COMM_WORLD);
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

	std::array<unsigned long int, 7> result_counters = {{
	    number_of_free_particles,
	    number_of_reflected_particles,
	    number_of_complete_evaporation_particles,
	    number_of_residence_samples,
	    number_of_censored_captured_particles,
	    number_of_outer_domain_removed_particles,
	    number_of_invalid_survival_captured_particles}};
	MPI_Trace_Point(mpi_rank, "before allreduce result counters");
	MPI_Allreduce(MPI_IN_PLACE, result_counters.data(), static_cast<int>(result_counters.size()), MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	number_of_free_particles = result_counters[0];
	number_of_reflected_particles = result_counters[1];
	number_of_complete_evaporation_particles = result_counters[2];
	number_of_residence_samples = result_counters[3];
	number_of_censored_captured_particles = result_counters[4];
	number_of_outer_domain_removed_particles = result_counters[5];
	number_of_invalid_survival_captured_particles = result_counters[6];

	auto reduce_jackknife_counts =
	    [&](std::array<unsigned long int, RESIDENCE_JACKKNIFE_BLOCKS>& counts,
	        const char* trace_point)
	{
		MPI_Trace_Point(mpi_rank, trace_point);
		MPI_Allreduce(
		    MPI_IN_PLACE,
		    counts.data(),
		    static_cast<int>(counts.size()),
		    MPI_UNSIGNED_LONG,
		    MPI_SUM,
		    MPI_COMM_WORLD);
	};
	reduce_jackknife_counts(
	    jackknife_attempted_counts,
	    "before allreduce jackknife attempted counts");
	reduce_jackknife_counts(
	    jackknife_captured_counts,
	    "before allreduce jackknife captured counts");
	reduce_jackknife_counts(
	    jackknife_completed_escape_counts,
	    "before allreduce jackknife completed escape counts");
	reduce_jackknife_counts(
	    jackknife_residence_sample_counts,
	    "before allreduce jackknife residence sample counts");
	reduce_jackknife_counts(
	    jackknife_invalid_counts,
	    "before allreduce jackknife invalid counts");
	reduce_jackknife_counts(
	    jackknife_outer_domain_removed_counts,
	    "before allreduce jackknife outer-domain counts");

   const int histogram_count =
       static_cast<int>(captured_dt_hist.size());
	MPI_Trace_Point(mpi_rank, "before allreduce captured_dt_hist");
	MPI_Allreduce(MPI_IN_PLACE, captured_dt_hist.data(), histogram_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce captured_v2dt_hist");
	MPI_Allreduce(MPI_IN_PLACE, captured_v2dt_hist.data(), histogram_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce captured_dt_sq_hist");
	MPI_Allreduce(MPI_IN_PLACE, captured_dt_sq_hist.data(), histogram_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce captured_v2dt_sq_hist");
	MPI_Allreduce(MPI_IN_PLACE, captured_v2dt_sq_hist.data(), histogram_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	const int jackknife_histogram_count =
	    static_cast<int>(
           residence_jackknife_block_dt_hist.size());
	MPI_Trace_Point(
	    mpi_rank,
	    "before allreduce jackknife dt histograms");
	MPI_Allreduce(
	    MPI_IN_PLACE,
	    residence_jackknife_block_dt_hist.data(),
	    jackknife_histogram_count,
	    MPI_DOUBLE,
	    MPI_SUM,
	    MPI_COMM_WORLD);
	MPI_Trace_Point(
	    mpi_rank,
	    "before allreduce jackknife v2dt histograms");
	MPI_Allreduce(
	    MPI_IN_PLACE,
	    residence_jackknife_block_v2dt_hist.data(),
	    jackknife_histogram_count,
	    MPI_DOUBLE,
	    MPI_SUM,
	    MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "after histogram allreduces");
	if(evaporation_diagnostics_enabled)
	{
		const int local_evap_count = static_cast<int>(evaporation_records.size());
		std::vector<int> evap_counts(mpi_processes, 0);
		MPI_Trace_Point(mpi_rank, "before gather diagnostic evap counts");
		MPI_Gather(&local_evap_count, 1, MPI_INT, mpi_rank == 0 ? evap_counts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

		constexpr int EVAPORATION_MPI_INT_FIELDS = 8;
		constexpr int EVAPORATION_MPI_UINT_FIELDS = 6;
		constexpr int EVAPORATION_MPI_DOUBLE_FIELDS = 32;
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
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 6] = evaporation_records[i].outer_domain_removed ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 7] = TerminationReason_Index(evaporation_records[i].termination_reason);

			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i] = static_cast<unsigned long long>(evaporation_records[i].trajectory_id);
			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 1] = static_cast<unsigned long long>(evaporation_records[i].number_of_scatterings);
			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 2] = static_cast<unsigned long long>(evaporation_records[i].number_of_bound_to_unbound);
			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 3] = static_cast<unsigned long long>(evaporation_records[i].number_of_recaptures);
			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 4] = static_cast<unsigned long long>(evaporation_records[i].number_of_integrator_steps_after_capture);
			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 5] = static_cast<unsigned long long>(evaporation_records[i].number_of_bound_exterior_arcs);

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
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 26] = evaporation_records[i].first_bound_exit_kepler_period_sec;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 27] = evaporation_records[i].last_bound_exit_kepler_period_sec;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 28] = evaporation_records[i].max_bound_exit_kepler_period_sec;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 29] = evaporation_records[i].first_bound_exit_exterior_time_sec;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 30] = evaporation_records[i].last_bound_exit_exterior_time_sec;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 31] = evaporation_records[i].max_bound_exit_exterior_time_sec;
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
				evaporation_records[i].outer_domain_removed = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 6] != 0);
				evaporation_records[i].termination_reason = static_cast<TrajectoryTerminationReason>(TerminationReason_Index(static_cast<TrajectoryTerminationReason>(global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 7])));
				evaporation_records[i].max_free_energy_drift_eV = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 12];
				evaporation_records[i].max_free_energy_drift_rel = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 13];
				evaporation_records[i].number_of_scatterings = static_cast<unsigned long int>(global_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 1]);
				evaporation_records[i].number_of_bound_to_unbound = static_cast<unsigned long int>(global_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 2]);
				evaporation_records[i].number_of_recaptures = static_cast<unsigned long int>(global_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 3]);
				evaporation_records[i].number_of_integrator_steps_after_capture = static_cast<unsigned long int>(global_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 4]);
				evaporation_records[i].number_of_bound_exterior_arcs = static_cast<unsigned long int>(global_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 5]);
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
				evaporation_records[i].first_bound_exit_kepler_period_sec = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 26];
				evaporation_records[i].last_bound_exit_kepler_period_sec = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 27];
				evaporation_records[i].max_bound_exit_kepler_period_sec = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 28];
				evaporation_records[i].first_bound_exit_exterior_time_sec = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 29];
				evaporation_records[i].last_bound_exit_exterior_time_sec = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 30];
				evaporation_records[i].max_bound_exit_exterior_time_sec = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 31];
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

	// Invalid trajectories are rare but scientifically high leverage because
	// excluding a long-lived capture can bias the residence measure. Gather a
	// compact, replayable ledger even when full trajectory diagnostics are off.
	std::ostringstream local_invalid_stream;
	local_invalid_stream << std::scientific << std::setprecision(17);
	for(const auto& record : invalid_trajectory_records)
	{
		local_invalid_stream
		    << record.rank << '\t'
		    << record.trajectory_id << '\t'
		    << static_cast<int>(record.failure_stage) << '\t'
		    << TerminationReason_Index(record.termination_reason) << '\t'
		    << NumericalFailureDetail_Index(
		           record.numerical_failure_detail)
		    << '\t'
		    << (record.initial_shift_ok ? 1 : 0) << '\t'
		    << (record.is_captured ? 1 : 0) << '\t'
		    << (record.survival_valid ? 1 : 0) << '\t'
		    << (record.event_observed ? 1 : 0) << '\t'
		    << record.number_of_scatterings << '\t'
		    << record.number_of_bound_to_unbound << '\t'
		    << record.number_of_recaptures << '\t'
		    << record.t_capture_s << '\t'
		    << record.t_termination_s << '\t'
		    << record.final_r_rsun << '\t'
		    << record.final_vr_km_s << '\t'
		    << record.final_speed_km_s << '\t'
		    << record.final_energy_eV << '\t'
		    << record.max_r_after_capture_rsun << '\t'
		    << record.max_free_energy_drift_eV << '\t'
		    << record.max_free_energy_drift_rel << '\t'
		    << record.failure_energy_before_step_eV << '\t'
		    << record.failure_energy_after_step_eV << '\t'
		    << record.failure_energy_at_boundary_eV << '\t'
		    << record.failure_reference_energy_eV << '\t'
		    << record.failure_boundary_vr_km_s << '\t'
		    << record.failure_attempted_step_s << '\t'
		    << record.failure_accepted_step_s << '\t'
		    << record.initial_time_s;
		for(double value : record.initial_position_km)
			local_invalid_stream << '\t' << value;
		for(double value : record.initial_velocity_km_s)
			local_invalid_stream << '\t' << value;
		local_invalid_stream
		    << '\t' << Encode_PRNG_State(record.rng_state_before_initial_conditions)
		    << '\t' << Encode_PRNG_State(record.rng_state_before_simulation)
		    << '\n';
	}
	const std::string local_invalid_text = local_invalid_stream.str();
	MPI_Trace_Point(mpi_rank, "before gather invalid ledger");
	const std::string combined_invalid_text =
	    Gather_MPI_Text_To_Root(local_invalid_text);
	MPI_Trace_Point(mpi_rank, "after gather invalid ledger");
	invalid_trajectory_records.clear();
	if(mpi_rank == 0)
	{
		std::istringstream lines(combined_invalid_text);
		std::string line;
		while(std::getline(lines, line))
		{
			if(line.empty())
				continue;
			InvalidTrajectoryRecord record;
			int failure_stage = 0;
			int termination_reason = 0;
			int numerical_failure_detail = 0;
			int initial_shift_ok = 0;
			int is_captured = 0;
			int survival_valid = 0;
			int event_observed = 0;
			std::string initial_rng_state;
			std::string simulation_rng_state;
			std::istringstream fields(line);
			fields
			    >> record.rank
			    >> record.trajectory_id
			    >> failure_stage
			    >> termination_reason
			    >> numerical_failure_detail
			    >> initial_shift_ok
			    >> is_captured
			    >> survival_valid
			    >> event_observed
			    >> record.number_of_scatterings
			    >> record.number_of_bound_to_unbound
			    >> record.number_of_recaptures;
			auto read_double = [&](double& value)
			{
				std::string token;
				fields >> token;
				if(!fields)
					return;
				char* end = nullptr;
				errno = 0;
				value = std::strtod(token.c_str(), &end);
				if(end == token.c_str() || *end != '\0' || errno == ERANGE)
					fields.setstate(std::ios::failbit);
			};
			read_double(record.t_capture_s);
			read_double(record.t_termination_s);
			read_double(record.final_r_rsun);
			read_double(record.final_vr_km_s);
			read_double(record.final_speed_km_s);
			read_double(record.final_energy_eV);
			read_double(record.max_r_after_capture_rsun);
			read_double(record.max_free_energy_drift_eV);
			read_double(record.max_free_energy_drift_rel);
			read_double(record.failure_energy_before_step_eV);
			read_double(record.failure_energy_after_step_eV);
			read_double(record.failure_energy_at_boundary_eV);
			read_double(record.failure_reference_energy_eV);
			read_double(record.failure_boundary_vr_km_s);
			read_double(record.failure_attempted_step_s);
			read_double(record.failure_accepted_step_s);
			read_double(record.initial_time_s);
			for(double& value : record.initial_position_km)
				read_double(value);
			for(double& value : record.initial_velocity_km_s)
				read_double(value);
			fields >> initial_rng_state >> simulation_rng_state;
			if(!fields)
				throw std::runtime_error(
				    "Perform_MPI_Reductions(): failed to parse invalid trajectory record.");
			record.failure_stage = static_cast<InvalidTrajectoryStage>(failure_stage);
			record.termination_reason = static_cast<TrajectoryTerminationReason>(
			    TerminationReason_Index(
			        static_cast<TrajectoryTerminationReason>(termination_reason)));
			record.numerical_failure_detail =
			    static_cast<TrajectoryNumericalFailureDetail>(
			        NumericalFailureDetail_Index(
			            static_cast<TrajectoryNumericalFailureDetail>(
			                numerical_failure_detail)));
			record.initial_shift_ok = initial_shift_ok != 0;
			record.is_captured = is_captured != 0;
			record.survival_valid = survival_valid != 0;
			record.event_observed = event_observed != 0;
			record.rng_state_before_initial_conditions =
			    Decode_PRNG_State(std::move(initial_rng_state));
			record.rng_state_before_simulation =
			    Decode_PRNG_State(std::move(simulation_rng_state));
			invalid_trajectory_records.push_back(std::move(record));
		}
		std::sort(
		    invalid_trajectory_records.begin(),
		    invalid_trajectory_records.end(),
		    [](const InvalidTrajectoryRecord& lhs, const InvalidTrajectoryRecord& rhs) {
			    if(lhs.rank != rhs.rank)
				    return lhs.rank < rhs.rank;
			    return lhs.trajectory_id < rhs.trajectory_id;
		    });
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

		std::ostringstream local_replay_stream;
		local_replay_stream << std::scientific << std::setprecision(17);
		for(const auto& replay : trajectory_replay_records)
		{
			local_replay_stream << replay.rank << '\t' << replay.trajectory_id << '\t' << replay.initial_time_s;
			for(double value : replay.initial_position_km)
				local_replay_stream << '\t' << value;
			for(double value : replay.initial_velocity_km_s)
				local_replay_stream << '\t' << value;
			local_replay_stream << '\t' << Encode_PRNG_State(replay.rng_state_before_initial_conditions)
			                    << '\t' << Encode_PRNG_State(replay.rng_state_before_simulation) << '\n';
		}
		const std::string local_replay_text = local_replay_stream.str();
		MPI_Trace_Point(mpi_rank, "before gather replay ledger");
		const std::string combined_replay_text =
		    Gather_MPI_Text_To_Root(local_replay_text);
		MPI_Trace_Point(mpi_rank, "after gather replay ledger");
		trajectory_replay_records.clear();
		if(mpi_rank == 0)
		{
			std::istringstream lines(combined_replay_text);
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
				replay.rng_state_before_initial_conditions =
				    Decode_PRNG_State(std::move(initial_rng_state));
				replay.rng_state_before_simulation =
				    Decode_PRNG_State(std::move(simulation_rng_state));
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
		f << "# accepted_evaporation_samples = " << number_of_complete_evaporation_particles << "\n";
		f << "# residence_samples = " << number_of_residence_samples << "\n";
		f << "# censored_captured = " << number_of_censored_captured_particles << "\n";
		f << "# outer_domain_removed_captured = " << number_of_outer_domain_removed_particles << "\n";
		f << "# invalid_survival_captured = " << number_of_invalid_survival_captured_particles << "\n";
		f << "# initial_shift_failures = " << number_of_initial_shift_failures << "\n";
		f << "# final_reflection_shift_failures = " << number_of_final_reflection_shift_failures << "\n";
		f << "# mpi_scheduler = dynamic_rma_work_queue_v1\n";
		f << "# mpi_scheduler_work_claims = " << mpi_scheduler_work_claims << "\n";
		f << "# mpi_scheduler_peak_in_flight = " << mpi_scheduler_peak_in_flight << "\n";
		f << "# capture_target_overshoot = " << capture_target_overshoot << "\n";
		f << "# sample_target_type = valid_complete_evaporation_within_radial_domain\n";
		f << "# invalid_trajectory_records = " << invalid_trajectory_records.size() << "\n";
		f << "# residence_jackknife_blocks = "
		  << RESIDENCE_JACKKNIFE_BLOCKS << "\n";
		f << "# residence_jackknife_assignment = "
		  << "splitmix64(base_seed,rank,trajectory_id)%"
		  << RESIDENCE_JACKKNIFE_BLOCKS << "\n";
		for(int reason_index = 0;
		    reason_index < TRAJECTORY_TERMINATION_REASON_COUNT;
		    reason_index++)
		{
			const auto reason = static_cast<TrajectoryTerminationReason>(reason_index);
			f << "# termination_" << Termination_Reason_Key(reason) << "_captured = "
			  << captured_termination_reason_counts[static_cast<size_t>(reason_index)] << "\n";
			f << "# termination_" << Termination_Reason_Key(reason) << "_uncaptured = "
			  << uncaptured_termination_reason_counts[static_cast<size_t>(reason_index)] << "\n";
		}
		for(int detail_index = 1;
		    detail_index
		        < TRAJECTORY_NUMERICAL_FAILURE_DETAIL_COUNT;
		    detail_index++)
		{
			const auto detail =
			    static_cast<TrajectoryNumericalFailureDetail>(
			        detail_index);
			f << "# numerical_failure_detail_"
			  << TrajectoryNumericalFailureDetailKey(detail)
			  << "_captured = "
			  << captured_numerical_failure_detail_counts[
			         static_cast<size_t>(detail_index)]
			  << "\n";
			f << "# numerical_failure_detail_"
			  << TrajectoryNumericalFailureDetailKey(detail)
			  << "_uncaptured = "
			  << uncaptured_numerical_failure_detail_counts[
			         static_cast<size_t>(detail_index)]
			  << "\n";
		}
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

	auto sum_counts = [](const auto& values) {
		return std::accumulate(values.begin(), values.end(), 0UL);
	};
	auto histogram_close = [](double block_sum, double total) {
		const double scale =
		    std::max({1.0, std::fabs(block_sum), std::fabs(total)});
		return std::fabs(block_sum - total) <= 1.0e-10 * scale;
	};

	const bool scalar_closure =
	    sum_counts(jackknife_attempted_counts) == number_of_trajectories
	    && sum_counts(jackknife_captured_counts)
	           == number_of_captured_particles
	    && sum_counts(jackknife_completed_escape_counts)
	           == number_of_completed_outward_escapes
	    && sum_counts(jackknife_residence_sample_counts)
	           == number_of_residence_samples
	    && sum_counts(jackknife_invalid_counts)
	           == number_of_numerical_failures
	                  + number_of_computational_truncations
	    && sum_counts(jackknife_outer_domain_removed_counts)
	           == number_of_outer_domain_removed_particles;
	if(!scalar_closure)
		throw std::runtime_error(
		    "residence jackknife scalar counts do not close");

   const std::size_t radial_bin_count = captured_dt_hist.size();

   if(captured_v2dt_hist.size() != radial_bin_count
      || captured_dt_sq_hist.size() != radial_bin_count
      || captured_v2dt_sq_hist.size() != radial_bin_count
      || residence_jackknife_block_dt_hist.size()
             != RESIDENCE_JACKKNIFE_BLOCKS * radial_bin_count
      || residence_jackknife_block_v2dt_hist.size()
             != RESIDENCE_JACKKNIFE_BLOCKS * radial_bin_count)
   {
           throw std::runtime_error(
               "runtime radial bincount storage dimensions are inconsistent");
   }

   for(std::size_t bin = 0; bin < radial_bin_count; bin++)
   {
           double dt_sum = 0.0;
           double v2dt_sum = 0.0;

           for(std::size_t block = 0;
               block < RESIDENCE_JACKKNIFE_BLOCKS;
               block++)
           {
                   const std::size_t index =
                       block * radial_bin_count + bin;

                   dt_sum += residence_jackknife_block_dt_hist[index];
                   v2dt_sum += residence_jackknife_block_v2dt_hist[index];
           }

           if(!histogram_close(dt_sum, captured_dt_hist[bin])
              || !histogram_close(v2dt_sum, captured_v2dt_hist[bin]))
           {
                   throw std::runtime_error(
                       "residence jackknife histogram does not close");
           }
   }
	const std::string bincount_path = output_dir + "/bincount.txt";
	const std::string bincount_temporary_path = bincount_path + ".tmp";
	const std::string jackknife_path =
	    output_dir + "/residence_jackknife_blocks.tsv";
	const std::string jackknife_temporary_path = jackknife_path + ".tmp";
	std::remove(bincount_temporary_path.c_str());
	std::remove(jackknife_temporary_path.c_str());

	// 1. Merged bincount
	{
		std::ofstream f(bincount_temporary_path);
		if(!f)
			throw std::runtime_error("failed to open temporary bincount.txt");
		write_header(f);
		const unsigned long int excluded_not_captured_particles = unresolved_not_captured_trajectories();
		f << "# residence_bincount_samples = " << number_of_residence_samples << "\n";
		if(excluded_not_captured_particles > 0)
			f << "# excluded_incomplete_not_captured = " << excluded_not_captured_particles << "\n";
           f << "# body_radius_km = "
             << std::scientific << std::setprecision(10)
             << bincount_radial_grid.Body_Radius_Km() << "\n";
           f << "# inner_grid_bins = "
             << bincount_radial_grid.Inner_Bin_Count() << "\n";
           f << "# exterior_bins = "
             << bincount_radial_grid.Exterior_Bin_Count() << "\n";
           f << "# total_radial_bins = "
             << radial_bin_count << "\n";
           f << "# radial_grid = uniform_inner_geometric_width_capped_exterior_v2\n";
           f << "# radial_inner_bin_width_Rbody = "
             << bincount_radial_grid.Inner_Bin_Width_Km()
                    / bincount_radial_grid.Body_Radius_Km() << "\n";
           f << "# exterior_initial_bin_width_Rbody = "
             << bincount_radial_grid.Inner_Bin_Width_Km()
                    / bincount_radial_grid.Body_Radius_Km() << "\n";
           f << "# exterior_bin_growth_factor = "
             << bincount_radial_grid.Exterior_Growth_Factor() << "\n";
           f << "# exterior_max_bin_width_Rbody = "
             << bincount_radial_grid.Exterior_Max_Bin_Width_Km()
                    / bincount_radial_grid.Body_Radius_Km() << "\n";
           f << "# radial_inner_extent_Rbody = "
             << bincount_radial_grid.Inner_Extent_Km()
                    / bincount_radial_grid.Body_Radius_Km() << "\n";
           f << "# radial_domain_max_AU = "
             << bincount_radial_grid.Domain_Max_Km()
                    / BINCOUNT_RADIAL_GRID_AU_KM << "\n";
           f << "# radial_extent_Rbody = "
             << bincount_radial_grid.Domain_Max_Km()
                    / bincount_radial_grid.Body_Radius_Km() << "\n";
           f << "# bin_index  r_lower_Rbody  r_upper_Rbody  "
                "residence_dt[s]  residence_v2dt[km2/s]  "
                "residence_err_dt[s]  residence_err_v2dt[km2/s]\n";

           const double residence_samples =
               static_cast<double>(number_of_residence_samples);

           for(std::size_t b = 0; b < radial_bin_count; b++)
           {
                   const double residence_err_dt =
                       Snapshot_Bin_Error(
                           captured_dt_hist[b],
                           captured_dt_sq_hist[b],
                           residence_samples);
                   const double residence_err_v2dt =
                       Snapshot_Bin_Error(
                           captured_v2dt_hist[b],
                           captured_v2dt_sq_hist[b],
                           residence_samples);

                   f << b << "\t"
                     << std::scientific << std::setprecision(10)
                     << bincount_radial_grid.Bin_Lower_Km(b)
                            / bincount_radial_grid.Body_Radius_Km()
                     << "\t"
                     << bincount_radial_grid.Bin_Upper_Km(b)
                            / bincount_radial_grid.Body_Radius_Km()
                     << "\t"
                     << captured_dt_hist[b]
                     << "\t"
                     << captured_v2dt_hist[b]
                     << "\t"
                     << residence_err_dt
                     << "\t"
                     << residence_err_v2dt
                     << "\n";
           }
		f.close();
		if(!f)
		{
			std::remove(bincount_temporary_path.c_str());
			throw std::runtime_error("failed to write temporary bincount.txt");
		}
	}

	// 2. Deterministic block totals for end-to-end delete-one-block
	// jackknife propagation.  These totals include both the capture
	// normalization counts and the accepted residence histograms.
	{
		std::ofstream blocks(jackknife_temporary_path);
		if(!blocks)
		{
			std::remove(bincount_temporary_path.c_str());
			throw std::runtime_error(
			    "failed to open temporary residence_jackknife_blocks.tsv");
		}
		blocks << "# DaMaSCUS-SUN residence jackknife blocks\n";
		blocks << "# format_version = 1\n";
		blocks << "# block_count = "
		       << RESIDENCE_JACKKNIFE_BLOCKS << "\n";
           blocks << "# radial_bins = "
                  << radial_bin_count << "\n";
		blocks << "# base_seed = " << diagnostic_base_seed << "\n";
		blocks << "# assignment = "
		          "splitmix64(base_seed,rank,trajectory_id)%"
		       << RESIDENCE_JACKKNIFE_BLOCKS << "\n";
		blocks << "# units = residence_dt:s "
		          "residence_v2dt:km2/s\n";
		for(std::size_t block = 0;
		    block < RESIDENCE_JACKKNIFE_BLOCKS;
		    block++)
		{
			blocks << "# block_" << block << "_attempted = "
			       << jackknife_attempted_counts[block] << "\n";
			blocks << "# block_" << block << "_captured = "
			       << jackknife_captured_counts[block] << "\n";
			blocks << "# block_" << block
			       << "_completed_uncaptured_escape = "
			       << jackknife_completed_escape_counts[block] << "\n";
			blocks << "# block_" << block << "_residence_samples = "
			       << jackknife_residence_sample_counts[block] << "\n";
			blocks << "# block_" << block << "_invalid = "
			       << jackknife_invalid_counts[block] << "\n";
			blocks << "# block_" << block
			       << "_outer_domain_removed = "
			       << jackknife_outer_domain_removed_counts[block] << "\n";
		}
		blocks << "block_id\tbin_index\tresidence_dt_s"
		          "\tresidence_v2dt_km2_s\n";
		blocks << std::scientific << std::setprecision(17);
		for(std::size_t block = 0;
		    block < RESIDENCE_JACKKNIFE_BLOCKS;
		    block++)
		{
                   for(std::size_t bin = 0;
                       bin < radial_bin_count;
                       bin++)
                   {
                           const std::size_t index =
                               block * radial_bin_count + bin;
				blocks << block << '\t' << bin << '\t'
				       << residence_jackknife_block_dt_hist[index]
				       << '\t'
				       << residence_jackknife_block_v2dt_hist[index]
				       << '\n';
			}
		}
		blocks.close();
		if(!blocks)
		{
			std::remove(bincount_temporary_path.c_str());
			std::remove(jackknife_temporary_path.c_str());
			throw std::runtime_error(
			    "failed to write temporary residence_jackknife_blocks.tsv");
		}
	}

	// Publish the companion jackknife file first and the bincount last.  Both
	// files are fully written and closed before either latest-output path is
	// replaced, so readers never observe a new bincount with a missing block
	// file.
	if(std::rename(jackknife_temporary_path.c_str(), jackknife_path.c_str()) != 0)
	{
		std::remove(bincount_temporary_path.c_str());
		std::remove(jackknife_temporary_path.c_str());
		throw std::runtime_error(
		    "failed to publish residence_jackknife_blocks.tsv");
	}
	if(std::rename(bincount_temporary_path.c_str(), bincount_path.c_str()) != 0)
	{
		std::remove(bincount_temporary_path.c_str());
		throw std::runtime_error("failed to publish bincount.txt");
	}
	std::cout << "Residence jackknife blocks: " << jackknife_path << std::endl;

	std::remove((output_dir + "/captured_bincount.txt").c_str());
	std::remove((output_dir + "/not_captured_bincount.txt").c_str());
	std::remove((output_dir + "/evaporation_diagnostics.txt").c_str());
	std::remove((output_dir + "/run_metadata.json").c_str());
	std::remove((output_dir + "/trajectory_summary.tsv").c_str());
	std::remove((output_dir + "/trajectory_events.tsv").c_str());
	std::remove((output_dir + "/invalid_trajectories.tsv").c_str());
	std::remove((output_dir + "/evaporation_" + "summary.txt").c_str());
	std::remove((output_dir + "/evaporation_" + "mode_summary.txt").c_str());
	std::remove((output_dir + "/evaporation_" + "mode_" + "bincount.txt").c_str());
	std::remove((output_dir + "/computation_" + "time_summary.txt").c_str());

	// 3. Final evaporation-time list.  This is intentionally the only final
	// evaporation report; snapshot files are intermediate progress reports.
	bool evaporation_times_ok = Write_Final_Evaporation_Time_File(Evaporation_Log_Path_From_Output_Dir(output_dir), mass_gev, sigma_cm2, compact_evaporation_events);
	if(!evaporation_times_ok)
		std::cerr << "Warning in Write_Output_Files(): failed to write evaporation_times.txt" << std::endl;

	// 4. Always-on replay ledger for trajectories excluded by numerical or
	// computational validity rules. A header-only file is written when no
	// invalid trajectory occurred so downstream checks never need to guess.
	{
		const std::string invalid_path = output_dir + "/invalid_trajectories.tsv";
		std::ofstream invalid(invalid_path);
		invalid << "# DaMaSCUS-SUN invalid trajectory replay ledger\n";
		invalid << "# format_version = 1\n";
		invalid << "# base_seed = " << diagnostic_base_seed << "\n";
		invalid << "# rank_seed_definition = base_seed + 1000003*rank\n";
		invalid << "# record_count = " << invalid_trajectory_records.size() << "\n";
		invalid << "# replay_definition = restore rng_state_before_simulation and simulate from the listed shifted initial state\n";
		invalid << "# rng_state_encoding = comma-separated std::mt19937 words\n";
		invalid << "# units = time:s position:km velocity:km/s radius:Rsun energy:eV\n";
			invalid << "rank\ttrajectory_id\tfailure_stage\ttermination_reason"
			        << "\tnumerical_failure_detail\tinitial_shift_ok"
			        << "\tis_captured\tsurvival_valid\tevent_observed\tn_scatter"
		        << "\tn_bound_to_unbound\tn_recapture\tt_capture_s\tt_termination_s"
		        << "\tfinal_r_Rsun\tfinal_vr_km_s\tfinal_speed_km_s\tfinal_energy_eV"
			        << "\tmax_r_after_capture_Rsun\tmax_free_energy_drift_eV"
			        << "\tmax_free_energy_drift_rel"
			        << "\tfailure_energy_before_step_eV"
			        << "\tfailure_energy_after_step_eV"
			        << "\tfailure_energy_at_boundary_eV"
			        << "\tfailure_reference_energy_eV"
			        << "\tfailure_boundary_vr_km_s"
			        << "\tfailure_attempted_step_s"
			        << "\tfailure_accepted_step_s"
			        << "\tinitial_time_s"
		        << "\tinitial_x_km\tinitial_y_km\tinitial_z_km"
		        << "\tinitial_vx_km_s\tinitial_vy_km_s\tinitial_vz_km_s"
		        << "\trng_state_before_initial_conditions\trng_state_before_simulation\n";
		invalid << std::scientific << std::setprecision(17);
		for(const auto& record : invalid_trajectory_records)
		{
			invalid
			    << record.rank << '\t'
				    << record.trajectory_id << '\t'
				    << Invalid_Trajectory_Stage_Key(record.failure_stage) << '\t'
				    << Termination_Reason_Key(record.termination_reason) << '\t'
				    << TrajectoryNumericalFailureDetailKey(
				           record.numerical_failure_detail)
				    << '\t'
				    << (record.initial_shift_ok ? 1 : 0) << '\t'
			    << (record.is_captured ? 1 : 0) << '\t'
			    << (record.survival_valid ? 1 : 0) << '\t'
			    << (record.event_observed ? 1 : 0) << '\t'
			    << record.number_of_scatterings << '\t'
			    << record.number_of_bound_to_unbound << '\t'
			    << record.number_of_recaptures << '\t'
			    << record.t_capture_s << '\t'
			    << record.t_termination_s << '\t'
			    << record.final_r_rsun << '\t'
			    << record.final_vr_km_s << '\t'
			    << record.final_speed_km_s << '\t'
			    << record.final_energy_eV << '\t'
				    << record.max_r_after_capture_rsun << '\t'
				    << record.max_free_energy_drift_eV << '\t'
				    << record.max_free_energy_drift_rel << '\t'
				    << record.failure_energy_before_step_eV << '\t'
				    << record.failure_energy_after_step_eV << '\t'
				    << record.failure_energy_at_boundary_eV << '\t'
				    << record.failure_reference_energy_eV << '\t'
				    << record.failure_boundary_vr_km_s << '\t'
				    << record.failure_attempted_step_s << '\t'
				    << record.failure_accepted_step_s << '\t'
				    << record.initial_time_s;
			for(double value : record.initial_position_km)
				invalid << '\t' << value;
			for(double value : record.initial_velocity_km_s)
				invalid << '\t' << value;
			invalid
			    << '\t' << Encode_PRNG_State(record.rng_state_before_initial_conditions)
			    << '\t' << Encode_PRNG_State(record.rng_state_before_simulation)
			    << '\n';
		}
		invalid.close();
		if(!invalid)
			std::cerr << "Warning in Write_Output_Files(): failed to write "
			          << invalid_path << std::endl;
		else
			std::cout << "Invalid trajectory ledger:\t" << invalid_path
			          << " (" << invalid_trajectory_records.size() << " records)"
			          << std::endl;
	}

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
		bool evaporation_event_reconciliation = true;
		bool escape_radius_invariant = true;
		bool residence_time_invariant = true;
		bool event_sequence_invariant = true;
		bool event_count_invariant = true;
		bool trace_selection_invariant = true;
		bool replay_state_invariant = true;
		bool bound_exit_orbit_invariant = true;
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
				   || std::fabs(escape_radius_rsun - In_Units(initial_and_final_radius, g_body_radius)) > 1.0e-10)
					escape_radius_invariant = false;
			}
			// Excluded/numerically invalid trajectories intentionally stop
			// before a complete observation window exists. The residence-time
			// invariant describes only records eligible for survival analysis.
			if(rec.survival_valid)
			{
				const double observation_end = rec.event_observed ? rec.t_boundary_escape : rec.t_termination;
				const double observed_duration = observation_end - rec.t_capture;
				const double residence_duration = rec.time_inside_sun_after_capture_sec
				                                + rec.time_outside_sun_after_capture_sec;
				const double residence_tolerance = std::max(1.0e-6, 1.0e-10 * std::fabs(observed_duration));
				if(!std::isfinite(observed_duration) || observed_duration < 0.0
				   || !std::isfinite(residence_duration)
				   || std::fabs(residence_duration - observed_duration) > residence_tolerance)
					residence_time_invariant = false;
			}
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
				evaporation_event_reconciliation = false;
			if(rec.number_of_bound_exterior_arcs == 0)
			{
				if(std::isfinite(rec.first_bound_exit_kepler_period_sec)
				   || std::isfinite(rec.last_bound_exit_kepler_period_sec)
				   || std::isfinite(rec.max_bound_exit_kepler_period_sec)
				   || std::isfinite(rec.first_bound_exit_exterior_time_sec)
				   || std::isfinite(rec.last_bound_exit_exterior_time_sec)
				   || std::isfinite(rec.max_bound_exit_exterior_time_sec))
					bound_exit_orbit_invariant = false;
			}
			else
			{
				const double periods[] = {
				    rec.first_bound_exit_kepler_period_sec,
				    rec.last_bound_exit_kepler_period_sec,
				    rec.max_bound_exit_kepler_period_sec,
				    rec.first_bound_exit_exterior_time_sec,
				    rec.last_bound_exit_exterior_time_sec,
				    rec.max_bound_exit_exterior_time_sec};
				for(double value : periods)
				if(!std::isfinite(value) || value <= 0.0)
					bound_exit_orbit_invariant = false;
				if(rec.first_bound_exit_kepler_period_sec > rec.max_bound_exit_kepler_period_sec
				   || rec.last_bound_exit_kepler_period_sec > rec.max_bound_exit_kepler_period_sec
				   || rec.first_bound_exit_exterior_time_sec > rec.max_bound_exit_exterior_time_sec
				   || rec.last_bound_exit_exterior_time_sec > rec.max_bound_exit_exterior_time_sec
				   || rec.first_bound_exit_exterior_time_sec > rec.first_bound_exit_kepler_period_sec
				   || rec.last_bound_exit_exterior_time_sec > rec.last_bound_exit_kepler_period_sec
				   || rec.max_bound_exit_exterior_time_sec > rec.max_bound_exit_kepler_period_sec)
					bound_exit_orbit_invariant = false;
			}
		}
		if(completed_keys.size() != number_of_complete_evaporation_particles)
			evaporation_event_reconciliation = false;

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
			         << "  \"schema_version\": \"trajectory-diagnostic-v5\",\n"
			         << "  \"run_id\": \"" << diagnostic_run_id << "\",\n"
			         << "  \"git_branch\": \"" << GIT_BRANCH << "\",\n"
			         << "  \"git_commit\": \"" << git_commit << "\",\n"
			         << "  \"git_dirty\": " << (git_dirty ? "true" : "false") << ",\n"
			         << "  \"build_type\": \"" << build_type << "\",\n"
			         << "  \"compiler\": \"" << compiler << "\",\n"
			         << "  \"mass_GeV\": " << std::scientific << std::setprecision(17) << mass_gev << ",\n"
			         << "  \"sigma_cm2\": " << sigma_cm2 << ",\n"
			         << "  \"trajectory_boundary_Rbody\": " << In_Units(initial_and_final_radius, g_body_radius) << ",\n"
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
                                   << "  \"radial_grid\": \"uniform_inner_geometric_width_capped_exterior_v2\",\n"
                                   << "  \"radial_body_radius_km\": "
                                   << bincount_radial_grid.Body_Radius_Km()
                                   << ",\n"
                                   << "  \"radial_inner_bins\": "
                                   << bincount_radial_grid.Inner_Bin_Count()
                                   << ",\n"
                                   << "  \"radial_inner_extent_Rbody\": "
                                   << bincount_radial_grid.Inner_Extent_Km()
                                          / bincount_radial_grid.Body_Radius_Km()
                                   << ",\n"
                                   << "  \"radial_total_bins\": "
                                   << bincount_radial_grid.Bin_Count()
                                   << ",\n"
                                   << "  \"radial_exterior_bins\": "
                                   << bincount_radial_grid.Exterior_Bin_Count()
                                   << ",\n"
                                   << "  \"radial_exterior_initial_bin_width_Rbody\": "
                                   << bincount_radial_grid.Inner_Bin_Width_Km()
                                          / bincount_radial_grid.Body_Radius_Km()
                                   << ",\n"
                                   << "  \"radial_exterior_bin_growth_factor\": "
                                   << bincount_radial_grid.Exterior_Growth_Factor()
                                   << ",\n"
                                   << "  \"radial_exterior_max_bin_width_Rbody\": "
                                   << bincount_radial_grid.Exterior_Max_Bin_Width_Km()
                                          / bincount_radial_grid.Body_Radius_Km()
                                   << ",\n"
                                   << "  \"outer_domain_removal_AU\": "
                                   << bincount_radial_grid.Domain_Max_Km()
                                          / BINCOUNT_RADIAL_GRID_AU_KM
                                   << ",\n"
			         << "  \"interpolation_points\": " << trajectory_diagnostic_config.interpolation_points << ",\n"
			         << "  \"max_optical_depth_step\": " << NormalModeMaxOpticalDepthStep() << ",\n"
			         << "  \"optical_depth_relative_tolerance\": " << OpticalDepthRelativeTolerance() << ",\n"
			         << "  \"energy_definition\": \"0.5*m_chi*(v^2-v_escape(r)^2); bound iff energy < 0\",\n"
			         << "  \"energy_unit\": \"eV\",\n"
			         << "  \"length_unit\": \"Rsun\",\n"
			         << "  \"velocity_unit\": \"km/s\",\n"
			         << "  \"angular_momentum_unit\": \"km^2/s\",\n"
			         << "  \"bound_exit_period_definition\": \"point-mass osculating Kepler period at a negative-energy outward crossing of 1.1 Rsun\",\n"
			         << "  \"bound_exit_exterior_time_definition\": \"analytic elapsed time: round trip through apoapsis for contained arcs; one-way to outward radial-domain removal for outer-domain arcs\",\n"
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
			         << "  \"bound_exit_orbit_invariant\": " << (bound_exit_orbit_invariant ? "true" : "false") << ",\n"
			         << "  \"evaporation_event_reconciliation\": " << (evaporation_event_reconciliation ? "true" : "false") << "\n"
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
			        << "\ttime_inside_sun_s\ttime_outside_sun_s\tn_bound_exterior_arcs"
			        << "\tP_kepler_first_bound_exit_s\tP_kepler_last_bound_exit_s\tP_kepler_max_bound_exit_s"
			        << "\tt_exterior_first_bound_exit_s\tt_exterior_last_bound_exit_s\tt_exterior_max_bound_exit_s"
			        << "\tn_integrator_steps\tmax_abs_ballistic_energy_drift_eV"
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
				        << rec.number_of_bound_exterior_arcs << '\t'
				        << rec.first_bound_exit_kepler_period_sec << '\t'
				        << rec.last_bound_exit_kepler_period_sec << '\t'
				        << rec.max_bound_exit_kepler_period_sec << '\t'
				        << rec.first_bound_exit_exterior_time_sec << '\t'
				        << rec.last_bound_exit_exterior_time_sec << '\t'
				        << rec.max_bound_exit_exterior_time_sec << '\t'
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

		if(!unique_keys || !time_invariants || !evaporation_event_reconciliation || !escape_radius_invariant
		   || !residence_time_invariant || !event_sequence_invariant || !event_count_invariant
		   || !trace_selection_invariant || !replay_state_invariant || !bound_exit_orbit_invariant)
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

unsigned long int Simulation_Data::Numerical_Failure_Detail_Count(
	TrajectoryNumericalFailureDetail detail,
	bool captured) const
{
	const size_t index = static_cast<size_t>(
	    NumericalFailureDetail_Index(detail));
	return captured
	     ? captured_numerical_failure_detail_counts[index]
	     : uncaptured_numerical_failure_detail_counts[index];
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
		Print_Termination_Reason_Summary(
		    captured_termination_reason_counts,
		    uncaptured_termination_reason_counts);
		Print_Numerical_Failure_Detail_Summary(
		    captured_numerical_failure_detail_counts,
		    uncaptured_numerical_failure_detail_counts);

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
				  << "MPI scheduler:\t\t\tdynamic RMA work queue" << std::endl
				  << "MPI work claims:\t\t" << mpi_scheduler_work_claims << std::endl
				  << "MPI peak in-flight tasks:\t" << mpi_scheduler_peak_in_flight << std::endl
				  << "Capture target overshoot:\t" << capture_target_overshoot << std::endl
				  << "Total # of scatterings:\t" << total_number_of_scatterings << std::endl
				  << "Numerical failure count:\t" << number_of_numerical_failures << std::endl
				  << "Computational truncations:\t" << number_of_computational_truncations << std::endl
				  << "Initial shift failures:\t\t" << number_of_initial_shift_failures << std::endl
				  << "Final reflection shift failures:\t" << number_of_final_reflection_shift_failures << std::endl
				  << "Numerical failure rate:\t\t" << std::fixed << std::setprecision(6) << Numerical_Failure_Ratio() << std::endl
				  << "Complete evaporation count:\t" << number_of_complete_evaporation_particles << std::endl
				  << "Residence sample count:\t\t" << number_of_residence_samples << std::endl
				  << "Censored captured count:\t" << number_of_censored_captured_particles << std::endl
				  << "Outer-domain removed captures:\t" << number_of_outer_domain_removed_particles << std::endl
				  << "Invalid survival count:\t\t" << number_of_invalid_survival_captured_particles << std::endl
				  << "Invalid ledger records:\t\t" << invalid_trajectory_records.size() << std::endl;
		const double invalid_trajectory_rate =
		    (number_of_trajectories > 0)
		    ? static_cast<double>(number_of_numerical_failures + number_of_computational_truncations)
		          / static_cast<double>(number_of_trajectories)
		    : 0.0;
		std::cout << "Numerical + truncation rate:\t" << std::fixed << std::setprecision(6)
		          << invalid_trajectory_rate << std::endl;
		Print_Termination_Reason_Summary(
		    captured_termination_reason_counts,
		    uncaptured_termination_reason_counts);
		Print_Numerical_Failure_Detail_Summary(
		    captured_numerical_failure_detail_counts,
		    uncaptured_numerical_failure_detail_counts);

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

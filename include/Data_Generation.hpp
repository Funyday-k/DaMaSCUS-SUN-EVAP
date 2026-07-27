#ifndef __Data_Generation_hpp_
#define __Data_Generation_hpp_

#include <vector>
#include <string>
#include <array>
#include <cstdint>
#include <limits>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/DM_Particle.hpp"

#include "Simulation_Trajectory.hpp"

namespace DaMaSCUS_SUN
{

bool TrajectoryTraceSelected(uint64_t trace_seed, int rank, uint64_t trajectory_id, double rate);

// Survival-analysis record for every captured trajectory.
struct EvaporationRecord
{
	int rank = -1;
	unsigned long int trajectory_id = 0;
	double completion_wall_time_sec = 0.0; // rank-local wall time when the trajectory finished [seconds]
	double t_evap = 0.0;    // finite only for observed unbinding events [seconds]
	double t_capture = -1.0;
	double t_first_unbinding_scatter = -1.0;
	double t_final_unbinding_scatter = -1.0;
	double t_boundary_escape = -1.0;
	double t_termination = -1.0;
	double observed_lifetime = 0.0;
	double lifetime_unbinding = -1.0;
	double lifetime_boundary = -1.0;
	double r_first_negative_km = -1.0;
	double E_first_negative_eV = 0.0;
	double dE_first_negative_from_prev_eV = 0.0;
	double r_first_unbinding_km = -1.0;
	double E_first_unbinding_eV = 0.0;
	double r_final_unbinding_km = -1.0;
	double E_final_unbinding_eV = 0.0;
	double r_boundary_escape_km = -1.0;
	double vr_boundary_escape_km_s = 0.0;
	double E_boundary_escape_eV = 0.0;
	bool event_observed = false;
	bool boundary_escape_observed = false;
	bool survival_valid = true;
	bool numerically_invalid_escape = false;
	bool censored = false;
	bool truncated = false;
	TrajectoryTerminationReason termination_reason = TrajectoryTerminationReason::Unknown;
	double max_free_energy_drift_eV = 0.0;
	double max_free_energy_drift_rel = 0.0;
	unsigned long int number_of_scatterings = 0;
	unsigned long int number_of_bound_to_unbound = 0;
	unsigned long int number_of_recaptures = 0;
	unsigned long int number_of_integrator_steps_after_capture = 0;
	double min_energy_after_capture_eV = 0.0;
	double max_radius_after_capture_km = 0.0;
	double time_inside_sun_after_capture_sec = 0.0;
	double time_outside_sun_after_capture_sec = 0.0;
	unsigned long int number_of_bound_exterior_arcs = 0;
	double first_bound_exit_kepler_period_sec = std::numeric_limits<double>::quiet_NaN();
	double last_bound_exit_kepler_period_sec = std::numeric_limits<double>::quiet_NaN();
	double max_bound_exit_kepler_period_sec = std::numeric_limits<double>::quiet_NaN();
	double first_bound_exit_exterior_time_sec = std::numeric_limits<double>::quiet_NaN();
	double last_bound_exit_exterior_time_sec = std::numeric_limits<double>::quiet_NaN();
	double max_bound_exit_exterior_time_sec = std::numeric_limits<double>::quiet_NaN();
};

struct CompactEvaporationEvent
{
	int rank = -1;
	uint64_t trajectory_id = 0;
	double completion_wall_time_sec = 0.0;
	double lifetime_unbinding = -1.0;
	double r_capture_rsun = -1.0;
	double E_capture_eV = 0.0;
	double dE_capture_eV = 0.0;
	uint64_t number_of_bound_exterior_arcs = 0;
	double first_bound_exit_kepler_period_sec = std::numeric_limits<double>::quiet_NaN();
	double last_bound_exit_kepler_period_sec = std::numeric_limits<double>::quiet_NaN();
	double max_bound_exit_kepler_period_sec = std::numeric_limits<double>::quiet_NaN();
	double first_bound_exit_exterior_time_sec = std::numeric_limits<double>::quiet_NaN();
	double last_bound_exit_exterior_time_sec = std::numeric_limits<double>::quiet_NaN();
	double max_bound_exit_exterior_time_sec = std::numeric_limits<double>::quiet_NaN();
};

struct TrajectoryReplayRecord
{
	int rank = -1;
	uint64_t trajectory_id = 0;
	double initial_time_s = 0.0;
	std::array<double, 3> initial_position_km{{0.0, 0.0, 0.0}};
	std::array<double, 3> initial_velocity_km_s{{0.0, 0.0, 0.0}};
	std::string rng_state_before_initial_conditions;
	std::string rng_state_before_simulation;
};

enum class SimulationStopReason
{
	None = 0,
	MaxTrajectoriesReached = 1,
	CaptureTargetNotReached = 2,
	InitialShiftFailureFractionExceeded = 3,
	InvalidTrajectoryFractionExceeded = 4
};

class Simulation_Data
{
  private:
	// Configuration
	unsigned int requested_captured_particles;
	unsigned long int max_trajectories_per_rank;
	unsigned long int normal_mode_mpi_sync_interval;
	double initial_and_final_radius = TRAJECTORY_BOUNDARY_RSUN * libphysica::natural_units::rSun;
	unsigned int minimum_number_of_scatterings = 1;
	unsigned long int maximum_number_of_scatterings = DEFAULT_MAXIMUM_SCATTERINGS;
	unsigned long int maximum_free_time_steps = DEFAULT_MAXIMUM_FREE_TIME_STEPS;

	// Results
	unsigned long int number_of_trajectories;
	unsigned long int number_of_free_particles;
	unsigned long int number_of_reflected_particles;
	unsigned long int number_of_captured_particles;
	unsigned long int number_of_completed_outward_escapes;
	unsigned long int number_of_complete_evaporation_particles;
	unsigned long int number_of_censored_captured_particles;
	unsigned long int number_of_radial_domain_excluded_particles;
	unsigned long int number_of_invalid_survival_captured_particles;
	unsigned long int number_of_initial_shift_failures;
	unsigned long int number_of_final_reflection_shift_failures;
	unsigned long int number_of_numerical_failures;
	unsigned long int number_of_computational_truncations;
	uint64_t total_number_of_scatterings;
	double average_number_of_scatterings;
	unsigned long int mpi_sync_rounds;
	unsigned long int final_mpi_round_trajectories;
	unsigned long int capture_target_overshoot;
	double computing_time;
	bool early_stopped;
	SimulationStopReason early_stop_reason;

	// Aggregated bincount histograms
	std::array<double, TOTAL_BINS> captured_dt_hist{};
	std::array<double, TOTAL_BINS> captured_v2dt_hist{};
	std::array<double, TOTAL_BINS> not_captured_dt_hist{};
	std::array<double, TOTAL_BINS> not_captured_v2dt_hist{};

	// Per-bin sum of squares for error estimation
	std::array<double, TOTAL_BINS> captured_dt_sq_hist{};      // Σ (per-traj dt)²
	std::array<double, TOTAL_BINS> captured_v2dt_sq_hist{};    // Σ (per-traj v²dt)²
	std::array<double, TOTAL_BINS> not_captured_dt_sq_hist{};
	std::array<double, TOTAL_BINS> not_captured_v2dt_sq_hist{};

	// Evaporation records
	std::vector<EvaporationRecord> evaporation_records;
	std::vector<CompactEvaporationEvent> compact_evaporation_events;
	std::vector<TrajectoryDiagnosticEvent> trajectory_diagnostic_events;
	std::vector<TrajectoryReplayRecord> trajectory_replay_records;
	bool evaporation_diagnostics_enabled = false;
	TrajectoryDiagnosticConfig trajectory_diagnostic_config;
	unsigned int diagnostic_base_seed = 0;
	uint64_t diagnostic_run_id = 0;

	// MPI
	int mpi_rank, mpi_processes;
	void Perform_MPI_Reductions(bool capture_mode);

	// For reflection spectrum compatibility (kept but not actively used in new logic)
	unsigned int isoreflection_rings;
	double minimum_speed_threshold;
	std::vector<unsigned long int> number_of_data_points;
	double KDE_boundary_correction_factor = 0.75;

  public:
	std::vector<std::vector<libphysica::DataPoint>> data;

	Simulation_Data(unsigned int sample_size, unsigned int max_trajectories, double u_min = 0.0, unsigned int iso_rings = 1);

	void Configure(double initial_radius, unsigned int min_scattering, unsigned long int max_scattering, unsigned long int max_free_steps = DEFAULT_MAXIMUM_FREE_TIME_STEPS);
	void Configure_Trajectory_Diagnostics(const TrajectoryDiagnosticConfig& config);

	void Generate_Data(obscura::DM_Particle& DM, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, SnapshotConfig snapshot_cfg = SnapshotConfig(), unsigned int fixed_seed = 0, bool capture_mode = false);

	// Output files
	void Write_Output_Files(const std::string& output_dir, obscura::DM_Particle& DM);

	double Free_Ratio() const;
	double Capture_Ratio() const;
	double Reflection_Ratio(int isoreflection_ring = -1) const;
	// Only physically classified captures and completed outward escapes
	// belong in the capture-rate normalization.
	unsigned long int Valid_Trajectories() const;
	double Free_Ratio_Valid() const;
	double Capture_Ratio_Valid() const;
	double Reflection_Ratio_Valid(int isoreflection_ring = -1) const;
	double Numerical_Failure_Ratio() const;

	double Minimum_Speed() const;
	double Lowest_Speed(unsigned int iso_ring = 0) const;
	double Highest_Speed(unsigned int iso_ring = 0) const;

	void Print_Summary(unsigned int mpi_rank = 0);
	void Print_Capture_Mode_Summary(unsigned int mpi_rank = 0);
};
}	// namespace DaMaSCUS_SUN
#endif

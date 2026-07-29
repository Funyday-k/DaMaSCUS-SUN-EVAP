#ifndef __Simulation_Trajectory_hpp_
#define __Simulation_Trajectory_hpp_

#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <array>
#include <chrono>
#include <cstdint>
#include <limits>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Particle.hpp"

#include "Simulation_Utilities.hpp"
#include "Solar_Model.hpp"

extern std::string g_top_level_dir;  // 从config文件读取的输出目录

namespace DaMaSCUS_SUN
{

class SnapshotRecorder;

// The numerical trajectory, injection surface, and exterior Kepler matching
// surface are deliberately identical.
constexpr double TRAJECTORY_BOUNDARY_RSUN = 1.1;
constexpr double R_SUN_KM = 6.957e5;  // km
constexpr double AU_KM = 1.495978707e8;  // IAU 2012 exact astronomical unit [km]
constexpr double BIN_WIDTH_KM = R_SUN_KM / 1000.0;  // 0.001 R_sun
constexpr int NUM_BINS = 1100;  // base grid through 1.1 R_sun
constexpr double BIN_MAX_KM = NUM_BINS * BIN_WIDTH_KM;
constexpr double RADIAL_DOMAIN_MAX_AU = 5.2;  // Jupiter-orbit cutoff
constexpr double RADIAL_DOMAIN_MAX_KM = RADIAL_DOMAIN_MAX_AU * AU_KM;
constexpr double RADIAL_DOMAIN_MAX_RSUN = RADIAL_DOMAIN_MAX_KM / R_SUN_KM;
constexpr std::size_t EXTERIOR_LOG_BINS = 512;
constexpr std::size_t TOTAL_BINS =
    static_cast<std::size_t>(NUM_BINS) + EXTERIOR_LOG_BINS;
constexpr unsigned long int DEFAULT_MAXIMUM_FREE_TIME_STEPS = 1000000000000UL;
constexpr unsigned long int DEFAULT_MAXIMUM_SCATTERINGS = 100000000000000UL;
constexpr int TRAJECTORY_TERMINATION_REASON_COUNT = 12;

struct BincountContribution
{
	int bin = -1;
	double dt_sec = 0.0;
	double v2dt_km2_per_sec = 0.0;
};

// Conservatively project one accepted trajectory interval onto the radial
// histogram. The trajectory dynamics are not re-integrated or otherwise
// changed; this only supplies dense-output quadrature for bincount.
void Compute_Bincount_Interval_Contributions(
	const Event& before,
	const Event& after,
	std::vector<BincountContribution>& contributions);

// The fixed histogram is uniform through 1.1 R_sun and logarithmic from there
// to the 5.2-AU Jupiter cutoff. These helpers are the single source of truth
// for writers, snapshots, and exact exterior Kepler shell integration.
double BincountBinLowerKm(std::size_t bin);
double BincountBinUpperKm(std::size_t bin);
int BincountBinIndexKm(double radius_km);

struct BoundKeplerExteriorArc
{
	Event terminal_event;
	double elapsed_time_sec = 0.0;
	double kepler_period_sec = 0.0;
	double apoapsis_km = 0.0;
	bool outer_domain_removed = false;
	std::array<double, TOTAL_BINS> dt_hist{};
	std::array<double, TOTAL_BINS> v2dt_hist{};
};

// Analytically propagate a negative-specific-energy outward crossing from the
// 1.1 R_sun matching surface. Orbits contained inside 5.2 AU return to the
// matching surface; larger orbits terminate at the outward 5.2-AU crossing.
// The returned histogram is a round trip in the first case and a one-way
// residence contribution in the second.
bool Compute_Bound_Kepler_Exterior_Arc(const Event& outward_event, BoundKeplerExteriorArc& arc);

enum class TrajectoryTerminationReason
{
	Unknown = 0,
	OutwardEscape = 1,
	Scatter = 2,
	WallTimeLimit = 3,
	MaxFreeSteps = 4,
	MaxScatterings = 5,
	NonFiniteState = 6,
	SpeedLimit = 7,
	NumericalFailure = 8,
	CaptureMode = 9,
	EnergyDriftEscape = 10,
	OuterDomainRemoval = 11
};

enum class TrajectoryNumericalFailureDetail
{
	None = 0,
	RK45ToleranceFailure = 1,
	NonFinitePropagationState = 2,
	InvalidScatteringRate = 3,
	OpticalDepthRetryExhausted = 4,
	UncapturedBoundMismatch = 5,
	BoundaryCrossingLocalizationFailure = 6,
	BoundaryDirectionMismatch = 7,
	BoundaryEnergyMismatch = 8,
	KeplerArcConstructionFailure = 9
};

constexpr int TRAJECTORY_NUMERICAL_FAILURE_DETAIL_COUNT = 10;

enum class TrajectoryDiagnosticEventType
{
	Capture = 0,
	ScatterPre = 1,
	ScatterPost = 2,
	CandidateUnbinding = 3,
	Recapture = 4,
	SunExit = 5,
	SunEntry = 6,
	EscapeValidated = 7,
	Censored = 8,
	NumericalFailure = 9
};

struct TrajectoryDiagnosticEvent
{
	int rank = -1;
	uint64_t trajectory_id = 0;
	uint64_t event_index = 0;
	uint64_t scatter_index = 0;
	uint64_t step_index = 0;
	TrajectoryDiagnosticEventType event_type = TrajectoryDiagnosticEventType::Capture;
	double t_s = 0.0;
	double r_km = 0.0;
	double vr_km_s = 0.0;
	double speed_km_s = 0.0;
	std::array<double, 3> position_km{{0.0, 0.0, 0.0}};
	std::array<double, 3> velocity_km_s{{0.0, 0.0, 0.0}};
	double energy_eV = 0.0;
	double angular_momentum_km2_s = 0.0;
	int is_bound = 0;
	int inside_sun = 0;
	int candidate_active = 0;
	char target_species[32] = {0};
	double ballistic_energy_drift_eV = 0.0;
};

const char* TrajectoryDiagnosticEventTypeKey(TrajectoryDiagnosticEventType type);
const char* TrajectoryNumericalFailureDetailKey(TrajectoryNumericalFailureDetail detail);
double RK45PositionToleranceKm();
double RK45VelocityToleranceKmPerSec();
double RK45PhaseTolerance();
double RK45AbsoluteMaxStepSec();
double NormalModeMaxOpticalDepthStep();
double OpticalDepthRelativeTolerance();
const char* BincountIntegrationScheme();
double BincountDensePositionToleranceKm();
double SnapshotProgressPublishWallIntervalSeconds();
bool SnapshotProgressPublishDue(
    unsigned long int accepted_steps_since_publish, double wall_seconds_since_publish, bool force);

bool TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason reason);

// Locate the first outward target-radius crossing on the cubic Hermite curve
// defined by two accepted trajectory states. Positions and velocities use
// libphysica natural units. Returns false if no outward crossing is bracketed.
bool Find_First_Outward_Hermite_Radius_Crossing(
	const Event& before,
	const Event& after,
	double target_radius,
	Event& crossing);

// Per-trajectory bincount result
struct TrajectoryBincount
{
	std::array<double, TOTAL_BINS> dt_hist{};     // Σ Δt per radial bin
	std::array<double, TOTAL_BINS> v2dt_hist{};   // Σ v²·Δt per radial bin

	// Capture/evaporation info
	bool is_captured = false;
	double t_capture = -1.0;         // first post-scatter transition to E < 0 [seconds]
	double t_last_bound = -1.0;      // latest time physically bound during free flight [seconds]
	double t_first_unbinding_scatter = std::numeric_limits<double>::quiet_NaN();  // first bound-to-unbound scatter [seconds]
	double t_final_unbinding_scatter = std::numeric_limits<double>::quiet_NaN();  // scatter that unbound the particle before a successful escape [seconds]
	double t_boundary_escape = std::numeric_limits<double>::quiet_NaN();          // outward escape boundary crossing [seconds]
	double t_termination = -1.0;     // final simulated time, event or censoring time [seconds]
	double r_first_negative_km = -1.0;  // radius when E first becomes negative [km]
	double E_first_negative_eV = std::numeric_limits<double>::quiet_NaN();  // first negative energy [eV]
	double dE_first_negative_from_prev_eV = std::numeric_limits<double>::quiet_NaN();  // E_now - E_previous_step [eV]
	double r_first_unbinding_km = std::numeric_limits<double>::quiet_NaN();
	double E_first_unbinding_eV = std::numeric_limits<double>::quiet_NaN();
	double r_final_unbinding_km = std::numeric_limits<double>::quiet_NaN();
	double E_final_unbinding_eV = std::numeric_limits<double>::quiet_NaN();
	double r_boundary_escape_km = std::numeric_limits<double>::quiet_NaN();
	double vr_boundary_escape_km_s = std::numeric_limits<double>::quiet_NaN();
	double E_boundary_escape_eV = std::numeric_limits<double>::quiet_NaN();
	bool event_observed = false;     // true if final unbinding was followed by outward escape
	bool boundary_escape_observed = false;  // true if the captured particle escaped through R_max with E >= 0
	bool survival_valid = true;      // false for captured trajectories invalid as survival samples
	bool numerically_invalid_escape = false;
	double max_free_energy_drift_eV = 0.0;
	double max_free_energy_drift_rel = 0.0;
	unsigned long int number_of_scatterings = 0;
	unsigned long int number_of_bound_to_unbound = 0;
	unsigned long int number_of_recaptures = 0;
	unsigned long int number_of_integrator_steps_after_capture = 0;
	double min_energy_after_capture_eV = std::numeric_limits<double>::quiet_NaN();
	double max_radius_after_capture_km = std::numeric_limits<double>::quiet_NaN();
	double time_inside_sun_after_capture_sec = 0.0;
	double time_outside_sun_after_capture_sec = 0.0;
	unsigned long int number_of_bound_exterior_arcs = 0;
	double first_bound_exit_kepler_period_sec = std::numeric_limits<double>::quiet_NaN();
	double last_bound_exit_kepler_period_sec = std::numeric_limits<double>::quiet_NaN();
	double max_bound_exit_kepler_period_sec = std::numeric_limits<double>::quiet_NaN();
	double first_bound_exit_exterior_time_sec = std::numeric_limits<double>::quiet_NaN();
	double last_bound_exit_exterior_time_sec = std::numeric_limits<double>::quiet_NaN();
	double max_bound_exit_exterior_time_sec = std::numeric_limits<double>::quiet_NaN();
	bool outer_domain_removed = false;  // true after physical removal at 5.2 AU
	TrajectoryTerminationReason termination_reason = TrajectoryTerminationReason::Unknown;
	TrajectoryNumericalFailureDetail numerical_failure_detail =
	    TrajectoryNumericalFailureDetail::None;
	double failure_energy_before_step_eV = std::numeric_limits<double>::quiet_NaN();
	double failure_energy_after_step_eV = std::numeric_limits<double>::quiet_NaN();
	double failure_energy_at_boundary_eV = std::numeric_limits<double>::quiet_NaN();
	double failure_reference_energy_eV = std::numeric_limits<double>::quiet_NaN();
	double failure_boundary_vr_km_s = std::numeric_limits<double>::quiet_NaN();
	double failure_attempted_step_s = std::numeric_limits<double>::quiet_NaN();
	double failure_accepted_step_s = std::numeric_limits<double>::quiet_NaN();

};

// Snapshot configuration: periodic cumulative bincount output
struct SnapshotConfig
{
	bool enabled = false;
	double interval_seconds = 60.0;  // wall-clock seconds between snapshots
	// 单条轨迹最长允许的 wall-clock 时间（秒）；超过即中止该轨迹。
	// 用于防止单条病态轨迹把整个 rank 卡死，从而导致 snapshot/MPI_Barrier 死锁。
	// 0 表示不限制。正式统计默认不限制。
	double max_trajectory_wall_time_sec = 0.0;
};

struct TrajectoryDiagnosticConfig
{
	bool summary_enabled = false;
	bool events_enabled = false;
	double trace_rate = 0.0;
	uint64_t trace_seed = 0;
	unsigned int interpolation_points = 0;
};

bool IsValidSnapshotIntervalSeconds(double interval_seconds);

// 1. Result of one trajectory
struct Trajectory_Result
{
	Event initial_event, final_event;
	unsigned long int number_of_scatterings;
	TrajectoryBincount bincount;
	std::vector<TrajectoryDiagnosticEvent> diagnostic_events;

	Trajectory_Result(const Event& event_ini, const Event& event_final, unsigned long int nScat, TrajectoryBincount bc,
	                  std::vector<TrajectoryDiagnosticEvent> events = std::vector<TrajectoryDiagnosticEvent>());

	bool Particle_Reflected() const;
	bool Particle_Free() const;
	bool Particle_Captured(Solar_Model& solar_model) const;

	void Print_Summary(Solar_Model& solar_model, unsigned int mpi_rank = 0);
};

// 2. Simulator
class Trajectory_Simulator
{
  private:
	Solar_Model solar_model;
	double v_max = 0.75;

	// Per-trajectory bincount accumulation
	TrajectoryBincount current_bincount;
	Event previous_bincount_event;
	bool has_previous_bincount_event;
	std::vector<BincountContribution> bincount_contribution_cache;
	double previous_capture_energy_eV;
	double free_flight_reference_energy_eV;
	double current_ballistic_energy_drift_eV;
	bool current_physical_bound_state;
	bool terminate_on_capture;

	void Accumulate_Bincount_Interval(const Event& before, const Event& after, double simulated_time_sec);
	void Reset_Bincount_Anchor(const Event& event);
	double Capture_Energy_eV(double radius, double speed, obscura::DM_Particle& DM);
	bool Update_Capture_State(double radius, double speed, double time, obscura::DM_Particle& DM, bool allow_new_capture);
	void Record_Diagnostic_Event(TrajectoryDiagnosticEventType type, const Event& event, obscura::DM_Particle& DM,
	                             const char* target_species = "");
	void Record_Surface_Crossing_Events(const Event& before, const Event& after, obscura::DM_Particle& DM);
	bool diagnostic_trace_enabled;
	uint64_t diagnostic_event_index;
	uint64_t diagnostic_scatter_index;
	uint64_t diagnostic_step_index;
	int last_scatter_target_index;
	std::vector<TrajectoryDiagnosticEvent> current_diagnostic_events;

	SnapshotRecorder* snapshot_recorder;
	bool trajectory_in_progress;
	bool track_trajectory_wall_time;
	std::chrono::steady_clock::time_point current_trajectory_wall_start;
	std::chrono::steady_clock::time_point last_snapshot_publish_wall_time;
	double accumulated_snapshot_overhead_sec;
	void Accumulate_Snapshot_Overhead(const std::chrono::steady_clock::time_point& operation_start);

	// Publish in batches, but bound staleness independently of the accepted-step
	// rate so a slow trajectory remains visible to the snapshot thread.
	unsigned long int steps_since_snapshot_publish;
	void Maybe_Publish_Snapshot_Progress(double simulated_time_sec, bool force);

	TrajectoryTerminationReason Propagate_Freely(Event& current_event, obscura::DM_Particle& DM);

	int Sample_Target(obscura::DM_Particle& DM, double r, double DM_speed);
	libphysica::Vector Sample_Target_Velocity(double temperature, double target_mass, const libphysica::Vector& vel_DM);
	libphysica::Vector New_DM_Velocity(double cos_scattering_angle, double DM_mass, double target_mass, libphysica::Vector& vel_DM, libphysica::Vector& vel_target);
	std::vector<double> rate_nuclei_cache;

  public:
	std::mt19937 PRNG;
	unsigned long int maximum_time_steps;
	unsigned long int maximum_scatterings;
	double maximum_distance;

	// 单条轨迹的 wall-clock 时间上限（秒）。超过后 Propagate_Freely 会返回 WallTimeLimit，
	// 防止任一 rank 被单条病态轨迹卡死从而阻塞 snapshot / MPI_Barrier。
	// 设为 0 表示不限制。正式统计默认不限制。
	double max_trajectory_wall_time_sec = 0.0;

	unsigned int current_mpi_rank;
	unsigned long int current_trajectory_id;

	Trajectory_Simulator(const Solar_Model& model, unsigned long int max_time_steps = DEFAULT_MAXIMUM_FREE_TIME_STEPS, unsigned long int max_scatterings = DEFAULT_MAXIMUM_SCATTERINGS, double max_distance = TRAJECTORY_BOUNDARY_RSUN * libphysica::natural_units::rSun);

	void Fix_PRNG_Seed(unsigned int fixed_seed);
	void Restore_PRNG_State(const std::string& serialized_state);
	std::string Serialize_PRNG_State() const;
	void Set_Snapshot_Recorder(SnapshotRecorder* recorder);
	void Enable_Capture_Mode(bool enabled);
	void Enable_Diagnostic_Trace(bool enabled);

	void Scatter(Event& current_event, obscura::DM_Particle& DM);
	Trajectory_Result Simulate(const Event& initial_condition, obscura::DM_Particle& DM, unsigned int mpi_rank);

	bool Trajectory_In_Progress() const;
	unsigned long int Current_Trajectory_ID() const;
	double Current_Trajectory_Wall_Time_Seconds() const;
	const TrajectoryBincount& Current_Trajectory_Bincount() const;
};

// 3. Equation of motion solution with Runge-Kutta-Fehlberg
class Free_Particle_Propagator
{
  private:
	double time, radius, phi, v_radial;
	double angular_momentum;
	libphysica::Vector axis_x, axis_y, axis_z;
	bool radial_mode;

	double dr_dt(double v);
	double dv_dt(double r, double mass);
	double dphi_dt(double r);
	double error_tolerances[3];

  public:
	double time_step = 0.1 * libphysica::natural_units::sec;

	// A Runge-Kutta step mutates only these scalars; the orbital basis, angular
	// momentum, and error tolerances are fixed at construction. Rolling a
	// rejected step back therefore needs no copy of the whole propagator.
	struct Scalar_State
	{
		double time;
		double radius;
		double phi;
		double v_radial;
		double time_step;
	};

	explicit Free_Particle_Propagator(const Event& event);

	bool Runge_Kutta_45_Step(Solar_Model& solar_model);
	bool Runge_Kutta_45_Step(double constant_mass);

	Scalar_State Save_Scalar_State() const;
	void Restore_Scalar_State(const Scalar_State& state);

	double Current_Time();
	double Current_Radius();
	double Current_Speed();

	Event Event_In_3D();
	// Overwrite an existing Event in place. libphysica::Vector allocates on every
	// construction and assignment, so the per-step propagation loop reuses one
	// Event buffer instead of building a fresh one.
	void Fill_Event_In_3D(Event& event) const;
};
}	// namespace DaMaSCUS_SUN

#endif

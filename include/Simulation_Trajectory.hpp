#ifndef __Simulation_Trajectory_hpp_
#define __Simulation_Trajectory_hpp_

#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <vector>
#include <array>
#include <chrono>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Particle.hpp"

#include "Simulation_Utilities.hpp"
#include "Solar_Model.hpp"

extern std::string g_top_level_dir;  // 从config文件读取的输出目录

namespace DaMaSCUS_SUN
{

// Bincount histogram constants
constexpr int NUM_BINS = 2000;
constexpr double R_SUN_KM = 6.957e5;  // km
constexpr double BIN_MAX_KM = 2.0 * R_SUN_KM;  // 2 R_sun in km
constexpr double BIN_WIDTH_KM = BIN_MAX_KM / NUM_BINS;  // ~695.7 km

// Per-trajectory bincount result
struct TrajectoryBincount
{
	std::array<double, NUM_BINS> dt_hist;     // Σ Δt per radial bin
	std::array<double, NUM_BINS> v2dt_hist;   // Σ v²·Δt per radial bin

	// Capture/evaporation info
	bool is_captured = false;
	double t_first_negative = -1.0;  // first time E <= 0 [seconds]
	double t_last_negative  = -1.0;  // last time E <= 0 [seconds]
	bool truncated = false;          // true if last step has E <= 0

	TrajectoryBincount()
	{
		dt_hist.fill(0.0);
		v2dt_hist.fill(0.0);
	}
};

// Snapshot configuration: periodic cumulative bincount output
struct SnapshotConfig
{
	bool enabled = false;
	double interval_seconds = 60.0;  // wall-clock seconds between snapshots
};

// 1. Result of one trajectory
struct Trajectory_Result
{
	Event initial_event, final_event;
	unsigned long int number_of_scatterings;
	unsigned long int total_rk45_steps;
	TrajectoryBincount bincount;

	Trajectory_Result(const Event& event_ini, const Event& event_final, unsigned long int nScat, unsigned long int rk45_steps, TrajectoryBincount bc);

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
	double prev_time_sec;       // previous step time in seconds (for dt calculation)
	double prev_r_km;           // previous step radius in km
	double prev_v2_km2s2;       // previous step v² in (km/s)²
	double prev_dt_sec;         // previous step dt in seconds (for last-step accumulation)

	void Accumulate_Bincount_Step(double r_km, double v2_km2s2, double dt_sec);

	// RK45 step counter for current trajectory
	unsigned long int total_rk45_steps_current_traj;
	std::function<void(const Trajectory_Simulator&)> snapshot_progress_callback;
	bool trajectory_in_progress;
	double current_trajectory_physical_time_sec;
	std::chrono::steady_clock::time_point current_trajectory_wall_start;

	void Publish_Snapshot_Progress() const;

	bool Propagate_Freely(Event& current_event, obscura::DM_Particle& DM);

	int Sample_Target(obscura::DM_Particle& DM, double r, double DM_speed);
	libphysica::Vector Sample_Target_Velocity(double temperature, double target_mass, const libphysica::Vector& vel_DM);
	libphysica::Vector New_DM_Velocity(double cos_scattering_angle, double DM_mass, double target_mass, libphysica::Vector& vel_DM, libphysica::Vector& vel_target);
	std::vector<double> rate_nuclei_cache;

  public:
	std::mt19937 PRNG;
	unsigned long int maximum_time_steps;
	unsigned long int maximum_scatterings;
	double maximum_distance;

	unsigned int current_mpi_rank;
	unsigned long int current_trajectory_id;

	Trajectory_Simulator(const Solar_Model& model, unsigned long int max_time_steps = 1e12, long int max_scatterings = 1e11, double max_distance = 2.0 * libphysica::natural_units::rSun);

	void Fix_PRNG_Seed(int fixed_seed);
	void Set_Snapshot_Progress_Callback(std::function<void(const Trajectory_Simulator&)> callback);

	void Scatter(Event& current_event, obscura::DM_Particle& DM);
	Trajectory_Result Simulate(const Event& initial_condition, obscura::DM_Particle& DM, unsigned int mpi_rank);

	bool Trajectory_In_Progress() const;
	unsigned long int Current_Trajectory_ID() const;
	double Current_Trajectory_Wall_Time_Seconds() const;
	double Current_Trajectory_Physical_Time_Seconds() const;
};

// 3. Equation of motion solution with Runge-Kutta-Fehlberg
class Free_Particle_Propagator
{
  private:
	double time, radius, phi, v_radial;
	double angular_momentum;
	libphysica::Vector axis_x, axis_y, axis_z;

	double dr_dt(double v);
	double dv_dt(double r, double mass);
	double dphi_dt(double r);
	double error_tolerances[3];

  public:
	double time_step = 0.1 * libphysica::natural_units::sec;

	explicit Free_Particle_Propagator(const Event& event);

	void Runge_Kutta_45_Step(Solar_Model& solar_model);
	void Runge_Kutta_45_Step(double constant_mass);

	double Current_Time();
	double Current_Radius();
	double Current_Speed();

	Event Event_In_3D();
};
}	// namespace DaMaSCUS_SUN

#endif

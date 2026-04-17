#ifndef __Simulation_Trajectory_hpp_
#define __Simulation_Trajectory_hpp_

#include <fstream>
#include <random>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Particle.hpp"

#include "Simulation_Utilities.hpp"
#include "Solar_Model.hpp"

#include <string>
extern std::string g_top_level_dir;  // 从config文件读取的输出目录
extern int g_recording_step_override;  // 步数记录覆盖值, 0=自动校准, >0=固定每N步记录

static long int* DMinBin= new long int[2000]; /////////////
static unsigned long int num_data_points;   ////////////

// 根据 DM 质量和截面返回推荐的轨迹记录间隔 (自然单位)
double Get_Recording_Interval(double mass_GeV, double sigma_cm2);

namespace DaMaSCUS_SUN
{

// 1. Result of one trajectory
struct Trajectory_Result
{
	Event initial_event, final_event;
	unsigned long int number_of_scatterings;

	Trajectory_Result(const Event& event_ini, const Event& event_final, unsigned long int nScat);

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

	unsigned int saved_trajectories, saved_trajectories_max;
	bool save_trajectories = false;
	double v_max		   = 0.75;

	bool Propagate_Freely(Event& current_event, obscura::DM_Particle& DM, std::ofstream& f);

	int Sample_Target(obscura::DM_Particle& DM, double r, double DM_speed);
	libphysica::Vector Sample_Target_Velocity(double temperature, double target_mass, const libphysica::Vector& vel_DM);
	libphysica::Vector New_DM_Velocity(double cos_scattering_angle, double DM_mass, double target_mass, libphysica::Vector& vel_DM, libphysica::Vector& vel_target);
	std::vector<double> rate_nuclei_cache;  // C: 预分配缓存，避免每次散射堆分配
	// --- 步数记录(替代时间记录，避免频闪混叠) ---
	int recording_step_interval;           // 每 N 步记录一次
	bool step_interval_calibrated;         // 是否已通过径向振荡完成校准
	int calibration_peri_count;            // 近日点检测计数
	unsigned long int peri_step_1;         // 第1个近日点的步数
	double calib_prev_r;                   // 自校准: 上一步半径
	double calib_prev_prev_r;              // 自校准: 上上步半径

	// --- 捕获检测：追踪轨迹中是否出现负能量 ---
	bool trajectory_has_negative_energy;
	unsigned int saved_trajectories_captured;
	unsigned int saved_trajectories_not_captured;

  public:
	std::mt19937 PRNG;
	unsigned long int maximum_time_steps;
	unsigned long int maximum_scatterings;
	double maximum_distance;

	Trajectory_Simulator(const Solar_Model& model, unsigned long int max_time_steps = 1e12, long int max_scatterings = 1e11, double max_distance = 2.0 * libphysica::natural_units::rSun);

	void Toggle_Trajectory_Saving(unsigned int max_trajectories = 500);
	void Fix_PRNG_Seed(int fixed_seed);

	unsigned int saved_trajectories_captured_count() const { return saved_trajectories_captured; }
	unsigned int saved_trajectories_not_captured_count() const { return saved_trajectories_not_captured; }

	void Scatter(Event& current_event, obscura::DM_Particle& DM);
	Trajectory_Result Simulate(const Event& initial_condition, obscura::DM_Particle& DM, unsigned int mpi_rank);////////
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
	double error_tolerances[3];  // B: 栈数组，避免每步堆分配

  public:
	double time_step = 0.1 * libphysica::natural_units::sec; //0.1 //statics is added by me

	explicit Free_Particle_Propagator(const Event& event);

	void Runge_Kutta_45_Step(double mass);

	double Current_Time();
	double Current_Radius();
	double Current_Speed();

	Event Event_In_3D();
};
}	// namespace DaMaSCUS_SUN

#endif

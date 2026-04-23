#include "Simulation_Trajectory.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <iostream>

#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/Astronomy.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

namespace
{
constexpr double RK45_MIN_STEP_FACTOR = 0.1;
constexpr double RK45_MAX_STEP_FACTOR = 4.0;

bool RK45_Errors_Within_Tolerance(const double errors[3], const double tolerances[3])
{
	for(int i = 0; i < 3; i++)
	{
		if(!std::isfinite(errors[i]) || errors[i] > tolerances[i])
			return false;
	}

	return true;
}

double RK45_Next_Step_Size(double current_step, const double errors[3], const double tolerances[3])
{
	double factor = RK45_MAX_STEP_FACTOR;

	for(int i = 0; i < 3; i++)
	{
		if(!std::isfinite(errors[i]) || errors[i] <= 0.0)
			continue;

		const double candidate = 0.84 * pow(tolerances[i] / errors[i], 0.25);
		if(std::isfinite(candidate))
			factor = std::min(factor, candidate);
	}

	factor = std::max(RK45_MIN_STEP_FACTOR, std::min(factor, RK45_MAX_STEP_FACTOR));
	return factor * current_step;
}
}

// 1. Result of one trajectory
Trajectory_Result::Trajectory_Result(const Event& event_ini, const Event& event_final, unsigned long int nScat, unsigned long int rk45_steps, TrajectoryBincount bc)
: initial_event(event_ini), final_event(event_final), number_of_scatterings(nScat), total_rk45_steps(rk45_steps), bincount(std::move(bc))
{
}

bool Trajectory_Result::Particle_Reflected() const
{
	double r	= final_event.Radius();
	double vesc = sqrt(2 * G_Newton * mSun / r);
	return r > rSun && final_event.Speed() > vesc && number_of_scatterings > 0;
}

bool Trajectory_Result::Particle_Free() const
{
	return number_of_scatterings == 0;
}

bool Trajectory_Result::Particle_Captured(Solar_Model& solar_model) const
{
	double r	= final_event.Radius();
	double vesc = solar_model.Local_Escape_Speed(r);
	return final_event.Speed() < vesc;
}

void Trajectory_Result::Print_Summary(Solar_Model& solar_model, unsigned int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR
				  << "Trajectory result summary" << std::endl
				  << std::endl
				  << "Number of scatterings:\t" << number_of_scatterings << std::endl
				  << "Simulation time [days]:\t" << libphysica::Round(In_Units(final_event.time, day)) << std::endl
				  << "Final radius [rSun]:\t" << libphysica::Round(In_Units(final_event.Radius(), rSun)) << std::endl
				  << "Final speed [km/sec]:\t" << libphysica::Round(In_Units(final_event.Speed(), km / sec)) << std::endl
				  << "Free particle:\t\t[" << (Particle_Free() ? "x" : " ") << "]" << std::endl
				  << "Captured:\t\t[" << (Particle_Captured(solar_model) ? "x" : " ") << "]" << std::endl
				  << "Reflection:\t\t[" << (Particle_Reflected() ? "x" : " ") << "]";

		if(Particle_Reflected())
		{
			double u_i_sqr = initial_event.Asymptotic_Speed_Sqr(solar_model);
			double u_f_sqr = final_event.Asymptotic_Speed_Sqr(solar_model);
			// 检查渐近速度平方是否为正
			if(u_i_sqr > 0.0 && u_f_sqr > 0.0)
			{
				double u_i = sqrt(u_i_sqr);
				double u_f = sqrt(u_f_sqr);
				std::cout << "\t(ratio u_f/u_i = " << libphysica::Round(u_f / u_i) << ")" << std::endl;
			}
			else
				std::cout << "\t(asymptotic speeds unavailable)" << std::endl;
		}
		else
			std::cout << std::endl;
		std::cout << SEPARATOR << std::endl;
	}
}

// 2. Simulator
Trajectory_Simulator::Trajectory_Simulator(const Solar_Model& model, unsigned long int max_time_steps, long int max_scatterings, double max_distance)
: solar_model(model), maximum_time_steps(max_time_steps), maximum_scatterings(max_scatterings), maximum_distance(max_distance), total_rk45_steps_current_traj(0), trajectory_in_progress(false), current_trajectory_physical_time_sec(0.0), current_mpi_rank(0), current_trajectory_id(0)
{
	std::random_device rd;
	PRNG.seed(rd());
	rate_nuclei_cache.resize(solar_model.target_isotopes.size());
}

void Trajectory_Simulator::Publish_Snapshot_Progress() const
{
	if(snapshot_progress_callback)
		snapshot_progress_callback(*this);
}

void Trajectory_Simulator::Set_Snapshot_Progress_Callback(std::function<void(const Trajectory_Simulator&)> callback)
{
	snapshot_progress_callback = std::move(callback);
}

bool Trajectory_Simulator::Trajectory_In_Progress() const
{
	return trajectory_in_progress;
}

unsigned long int Trajectory_Simulator::Current_Trajectory_ID() const
{
	return current_trajectory_id;
}

double Trajectory_Simulator::Current_Trajectory_Wall_Time_Seconds() const
{
	if(!trajectory_in_progress)
		return 0.0;

	return 1.0e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - current_trajectory_wall_start).count();
}

double Trajectory_Simulator::Current_Trajectory_Physical_Time_Seconds() const
{
	return current_trajectory_physical_time_sec;
}

// Accumulate one step into the current bincount
void Trajectory_Simulator::Accumulate_Bincount_Step(double r_km, double v2_km2s2, double dt_sec)
{
	if(dt_sec <= 0.0 || r_km < 0.0 || r_km >= BIN_MAX_KM)
		return;
	int bin_idx = static_cast<int>(r_km / BIN_WIDTH_KM);
	if(bin_idx < 0) bin_idx = 0;
	if(bin_idx >= NUM_BINS) return;
	current_bincount.dt_hist[bin_idx] += dt_sec;
	current_bincount.v2dt_hist[bin_idx] += v2_km2s2 * dt_sec;
}

bool Trajectory_Simulator::Propagate_Freely(Event& current_event, obscura::DM_Particle& DM)
{
	Free_Particle_Propagator particle_propagator(current_event);

	double minus_log_xi = -log(libphysica::Sample_Uniform(PRNG));
	bool success = false;
	unsigned long int time_steps = 0;

	while(time_steps < maximum_time_steps && !success)
	{
		time_steps++;
		double r_before = particle_propagator.Current_Radius();
		double t_before = particle_propagator.Current_Time();
		particle_propagator.Runge_Kutta_45_Step(solar_model);
		double actual_dt = particle_propagator.Current_Time() - t_before;
		double r_after = particle_propagator.Current_Radius();
		double v_after = particle_propagator.Current_Speed();

		if(v_after > v_max)
		{
			std::cerr << "\nWarning in Propagate_Freely(): DM speed exceeds the maximum of v_max = " << v_max << std::endl
					  << "\tAbort simulation." << std::endl;
			return false;
		}

		// --- Online bincount accumulation (every RK45 step) ---
		double t_now_sec = In_Units(particle_propagator.Current_Time(), sec);
		double r_now_km  = In_Units(r_after, km);
		double v_now_kms = In_Units(v_after, km / sec);
		double v2_now    = v_now_kms * v_now_kms;

		if(prev_time_sec >= 0.0)
		{
			double dt_sec = t_now_sec - prev_time_sec;
			// Accumulate previous step with forward difference dt
			Accumulate_Bincount_Step(prev_r_km, prev_v2_km2s2, dt_sec);
			prev_dt_sec = dt_sec;
		}

		// Energy check for capture detection
		double vesc = solar_model.Local_Escape_Speed(r_after);
		double E = 0.5 * DM.mass * (v_after * v_after - vesc * vesc);
		double E_eV = In_Units(E, eV);

		if(E_eV <= 0.0)
		{
			if(!current_bincount.is_captured)
			{
				current_bincount.is_captured = true;
				current_bincount.t_first_negative = t_now_sec;
			}
			current_bincount.t_last_negative = t_now_sec;
		}

		prev_time_sec = t_now_sec;
		prev_r_km     = r_now_km;
		prev_v2_km2s2 = v2_now;
		current_trajectory_physical_time_sec = t_now_sec;
		Publish_Snapshot_Progress();

		// Check for scatterings and reflection
		bool scattering = false;
		bool reflection = false;
		if(r_after < rSun)
		{
			if(v_after < 0.0)
			{
				std::cerr << "Warning: Negative velocity detected (v = " << v_after << ") at r = " << r_after << ", skipping scattering calculation." << std::endl;
				break;
			}
			double total_rate = solar_model.Total_DM_Scattering_Rate(DM, r_after, v_after);
			double time_step_max = (total_rate > 0.0) ? (0.1 / total_rate) : (1e30);
			if(particle_propagator.time_step > time_step_max)
				particle_propagator.time_step = time_step_max;
			minus_log_xi -= actual_dt * total_rate;
			if(minus_log_xi < 0.0)
				scattering = true;
		}
		else
		{
			if(r_before < maximum_distance && r_after > maximum_distance && v_after > solar_model.Local_Escape_Speed(r_after))
				reflection = true;
		}

		if(reflection || scattering)
			success = true;
	}

	// Accumulate the last step using the last known dt
	if(prev_dt_sec > 0.0)
	{
		Accumulate_Bincount_Step(prev_r_km, prev_v2_km2s2, prev_dt_sec);
	}

	// Accumulate RK45 steps for this Propagate_Freely call
	total_rk45_steps_current_traj += time_steps;

	current_event = particle_propagator.Event_In_3D();
	return success;
}


int Trajectory_Simulator::Sample_Target(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(r > rSun)
	{
		std::cerr << "Error in Trajectory_Simulator::Sample_Target(): r > rSun." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		// C: 复用预分配的 rate_nuclei_cache，避免每次散射事件堆分配
		for(unsigned int i = 0; i < solar_model.target_isotopes.size(); i++)
			rate_nuclei_cache[i] = solar_model.DM_Scattering_Rate_Nucleus(DM, r, DM_speed, i);
		double rate_electron = solar_model.DM_Scattering_Rate_Electron(DM, r, DM_speed);
		double total_rate	 = std::accumulate(rate_nuclei_cache.begin(), rate_nuclei_cache.end(), rate_electron);

		double xi = libphysica::Sample_Uniform(PRNG);
		// Electron
		double sum = rate_electron / total_rate;
		if(sum > xi)
			return -1;
		// Nuclei
		for(unsigned int i = 0; i < solar_model.target_isotopes.size(); i++)
		{
			sum += rate_nuclei_cache[i] / total_rate;
			if(sum > xi || i == solar_model.target_isotopes.size() - 1)
				return i;
		}
		std::cerr << "Error in Trajectory_Simulator::Sample_Target(): No target could be sampled." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

libphysica::Vector Trajectory_Simulator::Sample_Target_Velocity(double temperature, double target_mass, const libphysica::Vector& vel_DM)
{
	// Sampling algorithm taken from Romano & Walsh, "An improved target velocity sampling algorithm for free gas elastic scattering"
	double kappa = sqrt(target_mass / 2.0 / temperature);
	double vDM	 = vel_DM.Norm();
	// 1. Sample target speed vT and mu = cos alpha
	double y  = kappa * vDM;
	double x  = y;
	double mu = 1.0;
	while(sqrt(x * x + y * y - 2.0 * x * y * mu) / (x + y) < libphysica::Sample_Uniform(PRNG, 0.0, 1.0))
	{
		mu			= libphysica::Sample_Uniform(PRNG, -1.0, 1.0);
		double xi_1 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
		if(xi_1 < 2.0 / (sqrt(M_PI) * y + 2.0))
		{
			double xi_2 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
			double xi_3 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
			double z	= -log(xi_2 * xi_3);
			x			= sqrt(z);
		}
		else
		{
			double xi_2 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
			double xi_3 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
			double xi_4 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
			double z	= -log(xi_2) - pow(cos(M_PI / 2.0 * xi_3), 2.0) * log(xi_4);
			x			= sqrt(z);
		}
	}
	// 2. Construct the target velocity vel_T
	double vT		 = x / kappa;
	double cos_theta = mu;
	double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
	double phi		 = libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	double cos_phi	 = cos(phi);
	double sin_phi	 = sin(phi);

	libphysica::Vector unit_vector_DM = vel_DM.Normalized();
	double aux						  = sqrt(1.0 - pow(unit_vector_DM[2], 2.0));
	libphysica::Vector unit_vector_T({cos_theta * unit_vector_DM[0] + sin_theta / aux * (unit_vector_DM[0] * unit_vector_DM[2] * cos_phi - unit_vector_DM[1] * sin_phi),
									  cos_theta * unit_vector_DM[1] + sin_theta / aux * (unit_vector_DM[1] * unit_vector_DM[2] * cos_phi + unit_vector_DM[0] * sin_phi),
									  cos_theta * unit_vector_DM[2] - aux * cos_phi * sin_theta});

	return vT * unit_vector_T;
}

libphysica::Vector Trajectory_Simulator::New_DM_Velocity(double cos_scattering_angle, double DM_mass, double target_mass, libphysica::Vector& vel_DM, libphysica::Vector& vel_target)
{
	// Construction of n, the unit vector pointing into the direction of vfinal.
	double phi			 = libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	libphysica::Vector n = libphysica::Spherical_Coordinates(1.0, acos(cos_scattering_angle), phi, vel_DM);

	double relative_speed = (vel_target - vel_DM).Norm();

	return target_mass * relative_speed / (target_mass + DM_mass) * n + (DM_mass * vel_DM + target_mass * vel_target) / (target_mass + DM_mass);
}

void Trajectory_Simulator::Scatter(Event& current_event, obscura::DM_Particle& DM)
{
	double r = current_event.Radius();
	double v = current_event.Speed();
	// 1. Find target properties.
	int target_index = Sample_Target(DM, r, v);

	double target_mass;
	if(target_index == -1)
		target_mass = mElectron;
	else
		target_mass = solar_model.target_isotopes[target_index].mass;

	libphysica::Vector vel_target = Sample_Target_Velocity(solar_model.Temperature(r), target_mass, current_event.velocity);

	// 2. Sample the scattering angle
	double cos_alpha = (target_index == -1) ? DM.Sample_Scattering_Angle_Electron(PRNG, v, r) : DM.Sample_Scattering_Angle_Nucleus(PRNG, solar_model.target_isotopes[target_index], v, r);

	// 3. Construct the final DM velocity
	current_event.velocity = New_DM_Velocity(cos_alpha, DM.mass, target_mass, current_event.velocity, vel_target);
}

void Trajectory_Simulator::Fix_PRNG_Seed(int fixed_seed)
{
	PRNG.seed(fixed_seed);
}

Trajectory_Result Trajectory_Simulator::Simulate(const Event& initial_condition, obscura::DM_Particle& DM, unsigned int mpi_rank)
{
	Event current_event = initial_condition;
	long unsigned int number_of_scatterings = 0;

	// Initialize per-trajectory bincount
	current_bincount = TrajectoryBincount();
	prev_time_sec = -1.0;
	prev_r_km = 0.0;
	prev_v2_km2s2 = 0.0;
	prev_dt_sec = 0.0;
	current_mpi_rank = mpi_rank;
	current_trajectory_id++;
	total_rk45_steps_current_traj = 0;
	trajectory_in_progress = true;
	current_trajectory_physical_time_sec = In_Units(current_event.time, sec);
	current_trajectory_wall_start = std::chrono::steady_clock::now();
	Publish_Snapshot_Progress();

	while(Propagate_Freely(current_event, DM) && number_of_scatterings < maximum_scatterings)
	{
		if(current_event.Radius() < rSun)
		{
			Scatter(current_event, DM);
			number_of_scatterings++;
		}
		else
			break;
	}
	trajectory_in_progress = false;

	// Check truncation: if the last step still had negative energy
	if(current_bincount.is_captured)
	{
		// Compute energy at the final event
		double r_final = current_event.Radius();
		double v_final = current_event.Speed();
		double vesc_final = solar_model.Local_Escape_Speed(r_final);
		double E_final = 0.5 * DM.mass * (v_final * v_final - vesc_final * vesc_final);
		if(In_Units(E_final, eV) <= 0.0)
			current_bincount.truncated = true;
	}

	return Trajectory_Result(initial_condition, current_event, number_of_scatterings, total_rk45_steps_current_traj, current_bincount);
}  

// 3. Equation of motion solution with Runge-Kutta-Fehlberg
Free_Particle_Propagator::Free_Particle_Propagator(const Event& event)
{
	// 1. Coordinate system
	axis_x = event.position.Normalized();
	axis_z = event.position.Cross(event.velocity).Normalized();
	axis_y = axis_z.Cross(axis_x);

	// 2. Coordinates
	time			 = event.time;
	radius			 = event.Radius();
	phi				 = 0.0;
	v_radial		 = (radius == 0) ? event.Speed() : event.position.Dot(event.velocity) / radius;
	angular_momentum = (event.position.Cross(event.velocity)).Dot(axis_z);

	// 3. Error tolerances (fixed-size array, no heap allocation)
	error_tolerances[0] = 1.0 * km;
	error_tolerances[1] = 1.0e-3 * km / sec;
	error_tolerances[2] = 1.0e-7;
}

double Free_Particle_Propagator::dr_dt(double v)
{
	return v;
}

double Free_Particle_Propagator::dv_dt(double r, double mass)
{
	return angular_momentum * angular_momentum / r / r / r - G_Newton * mass / r / r;
}

double Free_Particle_Propagator::dphi_dt(double r)
{
	return angular_momentum / r / r;
}

void Free_Particle_Propagator::Runge_Kutta_45_Step(Solar_Model& solar_model)
{
	bool accepted = false;
	while(!accepted)
	{
	// RK coefficients:
	double k_r[6];
	double k_v[6];
	double k_p[6];
	double r_i;  // intermediate radius for each RK stage

	// Stage 0
	k_r[0] = time_step * dr_dt(v_radial);
	k_v[0] = time_step * dv_dt(radius, solar_model.Mass(radius));
	k_p[0] = time_step * dphi_dt(radius);

	// Stage 1
	r_i = radius + k_r[0] / 4.0;
	k_r[1] = time_step * dr_dt(v_radial + k_v[0] / 4.0);
	k_v[1] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	// k_p[1]=	dt*dphi_dt(radius+k_r[0]/4.0,J);

	// Stage 2
	r_i = radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1];
	k_r[2] = time_step * dr_dt(v_radial + 3.0 / 32.0 * k_v[0] + 9.0 / 32.0 * k_v[1]);
	k_v[2] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	k_p[2] = time_step * dphi_dt(r_i);

	// Stage 3
	r_i = radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2];
	k_r[3] = time_step * dr_dt(v_radial + 1932.0 / 2197.0 * k_v[0] - 7200.0 / 2197.0 * k_v[1] + 7296.0 / 2197.0 * k_v[2]);
	k_v[3] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	k_p[3] = time_step * dphi_dt(r_i);

	// Stage 4
	r_i = radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3];
	k_r[4] = time_step * dr_dt(v_radial + 439.0 / 216.0 * k_v[0] - 8.0 * k_v[1] + 3680.0 / 513.0 * k_v[2] - 845.0 / 4104.0 * k_v[3]);
	k_v[4] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	k_p[4] = time_step * dphi_dt(r_i);

	// Stage 5
	r_i = radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4];
	k_r[5] = time_step * dr_dt(v_radial - 8.0 / 27.0 * k_v[0] + 2.0 * k_v[1] - 3544.0 / 2565.0 * k_v[2] + 1859.0 / 4104.0 * k_v[3] - 11.0 / 40.0 * k_v[4]);
	k_v[5] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	k_p[5] = time_step * dphi_dt(r_i);

	// New values with Runge Kutta 4 and Runge Kutta 5
	double radius_4	  = radius + 25.0 / 216.0 * k_r[0] + 1408.0 / 2565.0 * k_r[2] + 2197.0 / 4104.0 * k_r[3] - 1.0 / 5.0 * k_r[4];
	double v_radial_4 = v_radial + 25.0 / 216.0 * k_v[0] + 1408.0 / 2565.0 * k_v[2] + 2197.0 / 4104.0 * k_v[3] - 1.0 / 5.0 * k_v[4];
	double phi_4	  = phi + 25.0 / 216.0 * k_p[0] + 1408.0 / 2565.0 * k_p[2] + 2197.0 / 4104.0 * k_p[3] - 1.0 / 5.0 * k_p[4];

	double radius_5	  = radius + 16.0 / 135.0 * k_r[0] + 6656.0 / 12825.0 * k_r[2] + 28561.0 / 56430.0 * k_r[3] - 9.0 / 50.0 * k_r[4] + 2.0 / 55.0 * k_r[5];
	double v_radial_5 = v_radial + 16.0 / 135.0 * k_v[0] + 6656.0 / 12825.0 * k_v[2] + 28561.0 / 56430.0 * k_v[3] - 9.0 / 50.0 * k_v[4] + 2.0 / 55.0 * k_v[5];
	double phi_5	  = phi + 16.0 / 135.0 * k_p[0] + 6656.0 / 12825.0 * k_p[2] + 28561.0 / 56430.0 * k_p[3] - 9.0 / 50.0 * k_p[4] + 2.0 / 55.0 * k_p[5];

	// Error and adapting the time step (B: 栈数组替代堆分配vector)
	double errors[3] = {fabs(radius_5 - radius_4), fabs(v_radial_5 - v_radial_4), fabs(phi_5 - phi_4)};
	double time_step_new = RK45_Next_Step_Size(time_step, errors, error_tolerances);

	// Check if errors fall below the tolerance
	// 自适应步长方法：根据误差调整 time_step
	if(RK45_Errors_Within_Tolerance(errors, error_tolerances))
	{
		time	  = time + time_step;
		radius	  = radius_4;
		// 边界检查：防止半径变为负数（数值误差）
		if(radius < 0.0)
			radius = 0.0;
		v_radial  = v_radial_4;
		phi		  = phi_4;
		time_step = time_step_new;
		accepted  = true;
	}
	else
	{
		time_step = time_step_new;
		// Loop continues with smaller time_step (was recursive call before Q1 fix)
	}
  }   // end while(!accepted)
}

// Overload with constant mass (for Kepler orbits outside the sun / testing)
void Free_Particle_Propagator::Runge_Kutta_45_Step(double constant_mass)
{
	bool accepted = false;
	while(!accepted)
	{
	double k_r[6];
	double k_v[6];
	double k_p[6];

	k_r[0] = time_step * dr_dt(v_radial);
	k_v[0] = time_step * dv_dt(radius, constant_mass);
	k_p[0] = time_step * dphi_dt(radius);

	k_r[1] = time_step * dr_dt(v_radial + k_v[0] / 4.0);
	k_v[1] = time_step * dv_dt(radius + k_r[0] / 4.0, constant_mass);

	k_r[2] = time_step * dr_dt(v_radial + 3.0 / 32.0 * k_v[0] + 9.0 / 32.0 * k_v[1]);
	k_v[2] = time_step * dv_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1], constant_mass);
	k_p[2] = time_step * dphi_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1]);

	k_r[3] = time_step * dr_dt(v_radial + 1932.0 / 2197.0 * k_v[0] - 7200.0 / 2197.0 * k_v[1] + 7296.0 / 2197.0 * k_v[2]);
	k_v[3] = time_step * dv_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2], constant_mass);
	k_p[3] = time_step * dphi_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2]);

	k_r[4] = time_step * dr_dt(v_radial + 439.0 / 216.0 * k_v[0] - 8.0 * k_v[1] + 3680.0 / 513.0 * k_v[2] - 845.0 / 4104.0 * k_v[3]);
	k_v[4] = time_step * dv_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3], constant_mass);
	k_p[4] = time_step * dphi_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3]);

	k_r[5] = time_step * dr_dt(v_radial - 8.0 / 27.0 * k_v[0] + 2.0 * k_v[1] - 3544.0 / 2565.0 * k_v[2] + 1859.0 / 4104.0 * k_v[3] - 11.0 / 40.0 * k_v[4]);
	k_v[5] = time_step * dv_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4], constant_mass);
	k_p[5] = time_step * dphi_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4]);

	double radius_4	  = radius + 25.0 / 216.0 * k_r[0] + 1408.0 / 2565.0 * k_r[2] + 2197.0 / 4104.0 * k_r[3] - 1.0 / 5.0 * k_r[4];
	double v_radial_4 = v_radial + 25.0 / 216.0 * k_v[0] + 1408.0 / 2565.0 * k_v[2] + 2197.0 / 4104.0 * k_v[3] - 1.0 / 5.0 * k_v[4];
	double phi_4	  = phi + 25.0 / 216.0 * k_p[0] + 1408.0 / 2565.0 * k_p[2] + 2197.0 / 4104.0 * k_p[3] - 1.0 / 5.0 * k_p[4];

	double radius_5	  = radius + 16.0 / 135.0 * k_r[0] + 6656.0 / 12825.0 * k_r[2] + 28561.0 / 56430.0 * k_r[3] - 9.0 / 50.0 * k_r[4] + 2.0 / 55.0 * k_r[5];
	double v_radial_5 = v_radial + 16.0 / 135.0 * k_v[0] + 6656.0 / 12825.0 * k_v[2] + 28561.0 / 56430.0 * k_v[3] - 9.0 / 50.0 * k_v[4] + 2.0 / 55.0 * k_v[5];
	double phi_5	  = phi + 16.0 / 135.0 * k_p[0] + 6656.0 / 12825.0 * k_p[2] + 28561.0 / 56430.0 * k_p[3] - 9.0 / 50.0 * k_p[4] + 2.0 / 55.0 * k_p[5];

	double errors[3] = {fabs(radius_5 - radius_4), fabs(v_radial_5 - v_radial_4), fabs(phi_5 - phi_4)};
	double time_step_new = RK45_Next_Step_Size(time_step, errors, error_tolerances);

	if(RK45_Errors_Within_Tolerance(errors, error_tolerances))
	{
		time	  = time + time_step;
		radius	  = radius_4;
		if(radius < 0.0)
			radius = 0.0;
		v_radial  = v_radial_4;
		phi		  = phi_4;
		time_step = time_step_new;
		accepted  = true;
	}
	else
	{
		time_step = time_step_new;
	}
  }
}

double Free_Particle_Propagator::Current_Time()
{
	return time;
}

double Free_Particle_Propagator::Current_Radius()
{
	// 边界检查：确保返回非负半径
	return (radius < 0.0) ? 0.0 : radius;
}

double Free_Particle_Propagator::Current_Speed()
{
	// 边界检查：确保 radius 非负
	double r = (radius < 0.0) ? 0.0 : radius;
	if(r == 0 || angular_momentum == 0)
		return fabs(v_radial);  // 返回速度绝对值
	else
		return sqrt(v_radial * v_radial + angular_momentum * angular_momentum / r / r);
}

Event Free_Particle_Propagator::Event_In_3D()
{
	double v_phi			= angular_momentum / pow(radius, 2);
	libphysica::Vector xNew = radius * (cos(phi) * axis_x + sin(phi) * axis_y);
	libphysica::Vector vNew = (v_radial * cos(phi) - v_phi * radius * sin(phi)) * axis_x + (v_radial * sin(phi) + radius * v_phi * cos(phi)) * axis_y;

	return Event(time, xNew, vNew);
}

}	// namespace DaMaSCUS_SUN

#include "Simulation_Trajectory.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
/////////
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
/////////

#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/Astronomy.hpp"

/////////
bool directoryExists(const std::string &path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0)
        return false;
    else if (info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}
/////////
namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

// 根据 DM 质量 (GeV) 和截面 (cm²) 返回推荐的记录间隔 (秒数, 非自然单位)
// 由 evaporation_theory.py 计算理论蒸发率 E_⊙(m_χ, σ) 生成
// 目标每条轨迹 ~10000 个数据点, Δt = τ_evap / N_target, 范围 [1, 1800] s
double Get_Recording_Interval(double mass_GeV, double sigma_cm2) {
    double log10_m = log10(mass_GeV);
    double log10_s = log10(sigma_cm2);

    // 查表: {log10(m/GeV), log10(σ/cm²), Δt(秒)}
    // 由 generate_recording_table.py 调用 evaporation_theory.py 自动生成
    struct Entry { double log10_m; double log10_s; double dt; };
    static const Entry table[] = {
        // m = 4.0 GeV (log10=0.602060)
        {  0.602060, -40.0, 1800.0}, {  0.602060, -39.0, 1800.0},
        {  0.602060, -38.0, 1800.0}, {  0.602060, -37.0, 1800.0},
        {  0.602060, -36.0, 1800.0}, {  0.602060, -35.0, 1800.0},
        {  0.602060, -34.0, 1800.0}, {  0.602060, -33.0, 1800.0},
        {  0.602060, -32.0, 1800.0}, {  0.602060, -31.0, 1800.0},
        {  0.602060, -30.0, 1800.0}, {  0.602060, -29.0, 1800.0},
        {  0.602060, -28.0, 1800.0}, {  0.602060, -27.0, 1800.0},
        {  0.602060, -26.0, 1800.0}, {  0.602060, -25.0, 1800.0},
        // m = 3.0 GeV (log10=0.477121)
        {  0.477121, -40.0, 1800.0}, {  0.477121, -39.0, 1800.0},
        {  0.477121, -38.0, 1800.0}, {  0.477121, -37.0, 1800.0},
        {  0.477121, -36.0, 1800.0}, {  0.477121, -35.0, 1800.0},
        {  0.477121, -34.0, 1800.0}, {  0.477121, -33.0, 1800.0},
        {  0.477121, -32.0, 1800.0}, {  0.477121, -31.0, 1800.0},
        {  0.477121, -30.0, 1800.0}, {  0.477121, -29.0, 1800.0},
        {  0.477121, -28.0, 1800.0}, {  0.477121, -27.0, 1800.0},
        {  0.477121, -26.0, 1800.0}, {  0.477121, -25.0, 1800.0},
        // m = 2.0 GeV (log10=0.301030)
        {  0.301030, -40.0, 1800.0}, {  0.301030, -39.0, 1800.0},
        {  0.301030, -38.0, 1800.0}, {  0.301030, -37.0, 1800.0},
        {  0.301030, -36.0, 1800.0}, {  0.301030, -35.0, 1800.0},
        {  0.301030, -34.0, 1800.0}, {  0.301030, -33.0, 1800.0},
        {  0.301030, -32.0, 1800.0}, {  0.301030, -31.0, 1800.0},
        {  0.301030, -30.0, 1800.0}, {  0.301030, -29.0, 1800.0},
        {  0.301030, -28.0, 1800.0}, {  0.301030, -27.0, 1800.0},
        {  0.301030, -26.0, 1800.0}, {  0.301030, -25.0, 1800.0},
        // m = 1.0 GeV (log10=0.000000)
        {  0.000000, -40.0, 1800.0}, {  0.000000, -39.0, 1800.0},
        {  0.000000, -38.0, 1800.0}, {  0.000000, -37.0, 1101.7},
        {  0.000000, -36.0,  146.7}, {  0.000000, -35.0,  114.4},
        {  0.000000, -34.0, 1800.0}, {  0.000000, -33.0, 1800.0},
        {  0.000000, -32.0, 1800.0}, {  0.000000, -31.0, 1800.0},
        {  0.000000, -30.0, 1800.0}, {  0.000000, -29.0, 1800.0},
        {  0.000000, -28.0, 1800.0}, {  0.000000, -27.0, 1800.0},
        {  0.000000, -26.0, 1800.0}, {  0.000000, -25.0, 1800.0},
        // m = 0.5 GeV (log10=-0.301030)
        { -0.301030, -40.0, 1800.0}, { -0.301030, -39.0, 1800.0},
        { -0.301030, -38.0,  835.9}, { -0.301030, -37.0,   85.5},
        { -0.301030, -36.0,   10.4}, { -0.301030, -35.0,    3.0},
        { -0.301030, -34.0,   16.2}, { -0.301030, -33.0,  251.4},
        { -0.301030, -32.0, 1800.0}, { -0.301030, -31.0, 1800.0},
        { -0.301030, -30.0, 1800.0}, { -0.301030, -29.0, 1800.0},
        { -0.301030, -28.0, 1800.0}, { -0.301030, -27.0, 1800.0},
        { -0.301030, -26.0, 1800.0}, { -0.301030, -25.0, 1800.0},
        // m = 0.1 GeV (log10=-1.000000)
        { -1.000000, -40.0, 1800.0}, { -1.000000, -39.0, 1800.0},
        { -1.000000, -38.0,  379.4}, { -1.000000, -37.0,   38.3},
        { -1.000000, -36.0,    3.7}, { -1.000000, -35.0,    1.0},
        { -1.000000, -34.0,    1.0}, { -1.000000, -33.0,    1.0},
        { -1.000000, -32.0,    3.0}, { -1.000000, -31.0,  100.3},
        { -1.000000, -30.0, 1800.0}, { -1.000000, -29.0, 1800.0},
        { -1.000000, -28.0, 1800.0}, { -1.000000, -27.0, 1800.0},
        { -1.000000, -26.0, 1800.0}, { -1.000000, -25.0, 1800.0},
    };
    static const int N = sizeof(table) / sizeof(table[0]);

    // 最近邻查找 (按 log10 空间的欧氏距离)
    double best_dist = 1e30;
    double best_dt = 600.0;  // 默认回退到原始值
    for (int i = 0; i < N; i++) {
        double dm = table[i].log10_m - log10_m;
        double ds = table[i].log10_s - log10_s;
        double d = dm * dm + ds * ds;
        if (d < best_dist) {
            best_dist = d;
            best_dt = table[i].dt;
        }
    }
    // 硬约束: [1, 1800] 秒
    if (best_dt < 1.0) best_dt = 1.0;
    if (best_dt > 1800.0) best_dt = 1800.0;
    return best_dt;
}

// 1. Result of one trajectory
Trajectory_Result::Trajectory_Result(const Event& event_ini, const Event& event_final, unsigned long int nScat)
: initial_event(event_ini), final_event(event_final), number_of_scatterings(nScat)
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
: solar_model(model), maximum_time_steps(max_time_steps), maximum_scatterings(max_scatterings), maximum_distance(max_distance)
{
	// Pseudo-random number generator
	std::random_device rd;
	PRNG.seed(rd());
	// C: 预分配核素散射率缓存
	rate_nuclei_cache.resize(solar_model.target_isotopes.size());
}

bool Trajectory_Simulator::Propagate_Freely(Event& current_event, obscura::DM_Particle& DM, std::ofstream& f)
{
	// 1. Define a equation-of-motion-solver in the orbital plane
	Free_Particle_Propagator particle_propagator(current_event);

	// 2. Simulate a free orbit
	double minus_log_xi			 = -log(libphysica::Sample_Uniform(PRNG));
	bool success				 = false; 
	unsigned long int time_steps = 0;
	//unsigned long int num_data_points = 0;   ////////////
	
	while(time_steps < maximum_time_steps && !success)
	{
		time_steps++;
		double r_before = particle_propagator.Current_Radius();
		particle_propagator.Runge_Kutta_45_Step(solar_model.Mass(r_before));
		double r_after = particle_propagator.Current_Radius();
		double v_after = particle_propagator.Current_Speed();
		double T_after = particle_propagator.Current_Time(); ////////////

		if(v_after > v_max)
		{
			std::cerr << "\nWarning in Propagate_Freely(): DM speed exceeds the maximum of v_max = " << v_max << std::endl
					  << "\tAbort simulation." << std::endl;
			return false;
		}

		double Time_point = 0; // (legacy, 不再使用 — 已改为基于步数的记录)

		// --- 运行时自校准: 检测径向振荡周期(近日点到近日点的步数) ---
		if(!step_interval_calibrated && time_steps > 3)
		{
			// 检测近日点: calib_prev_r 是局部最小值
			if(calib_prev_r < calib_prev_prev_r && calib_prev_r < r_after
			   && calib_prev_r < 0.9 * rSun)   // 确保在太阳内部
			{
				calibration_peri_count++;
				if(calibration_peri_count == 1)
				{
					peri_step_1 = time_steps - 1;
				}
				else if(calibration_peri_count == 2)
				{
					unsigned long int steps_per_orbit = (time_steps - 1) - peri_step_1;
					recording_step_interval = std::max(1, (int)(steps_per_orbit / 50));
					step_interval_calibrated = true;
				}
			}
			calib_prev_prev_r = calib_prev_r;
			calib_prev_r = r_after;
		}

		//if(save_trajectories && time_steps % 20 == 0)
		if(save_trajectories && time_steps % recording_step_interval == 0)   //////////////
		{
			num_data_points++;  ////////////
			Event event	 = particle_propagator.Event_In_3D();
			double r	 = event.Radius();
			double v	 = event.Speed();
			double vesc2 = solar_model.Local_Escape_Speed(r) * solar_model.Local_Escape_Speed(r);
			//double E	 = 1.0 / DM.mass * (v * v - vesc2);
			double E	 = 0.5 * DM.mass * (v * v - vesc2);

			//long int *DMinBin;//////////some problems
			int binsorder = (int)(r/(0.001 * rSun));    ///////////
			if(r < 2.0*rSun){
				DMinBin[binsorder] = DMinBin[binsorder]+1;  ///////////
				//std::cout<<binsorder<<'\t'<<DMinBin[binsorder]<<'\n';
			}
			
			  //////////
			//std::cout<< "E" << E << "\t";
			// 旧输出格式：10 x float64 = 80 bytes（含序号、能量、半径，已注释）
			// double record[10] = {
			// 	(double)num_data_points,
			// 	In_Units(event.time, sec),
			// 	In_Units(event.position[0], km),
			// 	In_Units(event.position[1], km),
			// 	In_Units(event.position[2], km),
			// 	In_Units(event.velocity[0], km / sec),
			// 	In_Units(event.velocity[1], km / sec),
			// 	In_Units(event.velocity[2], km / sec),
			// 	In_Units(E, eV),
			// 	In_Units(event.Radius(), km)
			// };

			// 输出格式：6 x float32 = 24 bytes
			// Columns: time[s], r[km], vx[km/s], vy[km/s], vz[km/s], E[eV]
			float record[6] = {
				(float)In_Units(event.time, sec),
				(float)In_Units(r, km),
				(float)In_Units(event.velocity[0], km / sec),
				(float)In_Units(event.velocity[1], km / sec),
				(float)In_Units(event.velocity[2], km / sec),
				(float)In_Units(E, eV)
			};
			f.write(reinterpret_cast<const char*>(record), sizeof(record));

			// 标记是否出现负能量（被捕获的标志）
			if(E <= 0)
				trajectory_has_negative_energy = true;
		}

		// Check for scatterings and reflection
		bool scattering = false;
		bool reflection = false;
		if(r_after < rSun)
		{ 
			// 边界检查：确保速度为正
			if(v_after < 0.0)
			{
				std::cerr << "Warning: Negative velocity detected (v = " << v_after << ") at r = " << r_after << ", skipping scattering calculation." << std::endl;
				break;
			}
			double total_rate	 = solar_model.Total_DM_Scattering_Rate(DM, r_after, v_after);
			double time_step_max = (total_rate > 0.0) ? (0.1 / total_rate) : (1e30);
			// 太阳内部：仅受散射率约束，无硬性步长上限
			// 当散射率低时，RK45 自适应步长可自由增长（数秒甚至更大）
			if(particle_propagator.time_step > time_step_max)
				particle_propagator.time_step = time_step_max;
			minus_log_xi -= particle_propagator.time_step * total_rate;
			if(minus_log_xi < 0.0)
				scattering = true;
		}
		else
		{
			// 太阳外部：不干预 RK45 自适应步长，让它自由增长
			// （椭圆轨道远日点附近步长可增长到 ~秒级，大幅减少步数）
			// 仅检查粒子是否逃逸
			if(r_before < maximum_distance && r_after > maximum_distance && v_after > solar_model.Local_Escape_Speed(r_after))
				reflection = true;
		}

		if(reflection || scattering)
			success = true;
	}

	//std::cout << "timesteps number" << time_steps << '\n'; ////////
	//std::cout << "if success" << success << "\n";
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

void Trajectory_Simulator::Toggle_Trajectory_Saving(unsigned int max_trajectories)
{
	saved_trajectories	   = 0;
	saved_trajectories_max = max_trajectories;
	save_trajectories	   = !save_trajectories;
	saved_trajectories_captured = 0;
	saved_trajectories_not_captured  = 0;
}

void Trajectory_Simulator::Fix_PRNG_Seed(int fixed_seed)
{
	PRNG.seed(fixed_seed);
}

Trajectory_Result Trajectory_Simulator::Simulate(const Event& initial_condition, obscura::DM_Particle& DM, unsigned int mpi_rank) //////
{
	////////////
	//std::string path = TOP_LEVEL_DIR "results_" + std::to_string(In_Units(DM.mass,GeV)) + "_" + std::to_string(In_Units(DM.Sigma_Proton(),cm*cm)) + "/"; 我改变了这里的路径，使得输出文件存在/scratch/s1/kennyng/DaMaSCUS_OUT/下
	std::string path = g_top_level_dir + "results_" + std::to_string(log10(In_Units(DM.mass,GeV))) + "_" + std::to_string(log10(In_Units(DM.Sigma_Proton(),cm*cm))) + "/";
	if (!directoryExists(path))
	{
		if (mkdir(path.c_str(), 0755) != 0)
		{
			std::cerr << "Error: Unable to create directory: " << path << std::endl;
		}
	}
	////////////
	std::ofstream f;
	std::string trajectory_filepath;
	if(save_trajectories && saved_trajectories < saved_trajectories_max)
	{
		saved_trajectories++;

		// std::cout << "trajectory is" << saved_trajectories << "\n";  //////////////
		//std::string path = TOP_LEVEL_DIR "results/";
		trajectory_filepath = path + "trajectory_" + std::to_string(saved_trajectories) + "_task" + std::to_string(mpi_rank) + ".dat";
		f.open(trajectory_filepath, std::ios::binary);
	}
	Event current_event						= initial_condition;
	long unsigned int number_of_scatterings = 0;
	
	//bool abcd = false; ///////////
	//static long int* DMinBin = new long int[50]; /////////////
	for(int i=0; i<2000; i++){
		DMinBin[i] = 0; //////////
	}

	// 初始化步数记录的自校准状态
	if(g_recording_step_override > 0)
	{
		recording_step_interval  = g_recording_step_override;
		step_interval_calibrated = true;   // 手动指定，跳过校准
	}
	else
	{
		recording_step_interval  = 10;     // 校准前默认: 每10步记录一次
		step_interval_calibrated = false;
	}
	calibration_peri_count     = 0;
	peri_step_1                = 0;
	calib_prev_r               = 0.0;
	calib_prev_prev_r          = 0.0;

	static bool recording_printed = false;
	if(!recording_printed)
	{
		double mass_gev = In_Units(DM.mass, GeV);
		double sigma_cm2 = In_Units(DM.Sigma_Proton(), cm * cm);
		if(g_recording_step_override > 0)
			std::cout << "[Recording] m=" << mass_gev << " GeV, sigma=" << sigma_cm2
			          << " cm^2 -> step-based (N_step=" << recording_step_interval
			          << ", override)" << std::endl;
		else
			std::cout << "[Recording] m=" << mass_gev << " GeV, sigma=" << sigma_cm2
			          << " cm^2 -> step-based (N_step_init=" << recording_step_interval
			          << ", auto-calibrate on)" << std::endl;
		recording_printed = true;
	}

	num_data_points = 0; ////////
	trajectory_has_negative_energy = false;
	while(Propagate_Freely(current_event, DM, f) && number_of_scatterings < maximum_scatterings)
	{
		
		if(current_event.Radius() < rSun)
		{
			Scatter(current_event, DM);
			number_of_scatterings++;
		}
		else
			break;
	}
	if(save_trajectories)
	{
		f.close();
		// 统计捕获与未捕获的轨迹文件数（不删除文件）
		if(!trajectory_filepath.empty())
		{
			if(trajectory_has_negative_energy)
				saved_trajectories_captured++;
			else
				saved_trajectories_not_captured++;
		}
	}	

	
	//std::cout << "trajectory is" << saved_trajectories << "," << "number of scatterings" << number_of_scatterings << "\n"; ///////
	/*if(abcd == false) 
	    std::cout << "false" << "\n";
	else if (abcd == true)
	{
		std::cout << "true" << "\n";
	}
	else
	    std::cout << "something wrong" << "\n";*/

	return Trajectory_Result(initial_condition, current_event, number_of_scatterings);
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

void Free_Particle_Propagator::Runge_Kutta_45_Step(double mass)
{
	// RK coefficients:
	double k_r[6];
	double k_v[6];
	double k_p[6];

	k_r[0] = time_step * dr_dt(v_radial);
	k_v[0] = time_step * dv_dt(radius, mass);
	k_p[0] = time_step * dphi_dt(radius);

	k_r[1] = time_step * dr_dt(v_radial + k_v[0] / 4.0);
	k_v[1] = time_step * dv_dt(radius + k_r[0] / 4.0, mass);
	// k_p[1]=	dt*dphi_dt(radius+k_r[0]/4.0,J);

	k_r[2] = time_step * dr_dt(v_radial + 3.0 / 32.0 * k_v[0] + 9.0 / 32.0 * k_v[1]);
	k_v[2] = time_step * dv_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1], mass);
	k_p[2] = time_step * dphi_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1]);

	k_r[3] = time_step * dr_dt(v_radial + 1932.0 / 2197.0 * k_v[0] - 7200.0 / 2197.0 * k_v[1] + 7296.0 / 2197.0 * k_v[2]);
	k_v[3] = time_step * dv_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2], mass);
	k_p[3] = time_step * dphi_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2]);

	k_r[4] = time_step * dr_dt(v_radial + 439.0 / 216.0 * k_v[0] - 8.0 * k_v[1] + 3680.0 / 513.0 * k_v[2] - 845.0 / 4104.0 * k_v[3]);
	k_v[4] = time_step * dv_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3], mass);
	k_p[4] = time_step * dphi_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3]);

	k_r[5] = time_step * dr_dt(v_radial - 8.0 / 27.0 * k_v[0] + 2.0 * k_v[1] - 3544.0 / 2565.0 * k_v[2] + 1859.0 / 4104.0 * k_v[3] - 11.0 / 40.0 * k_v[4]);
	k_v[5] = time_step * dv_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4], mass);
	k_p[5] = time_step * dphi_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4]);

	// New values with Runge Kutta 4 and Runge Kutta 5
	double radius_4	  = radius + 25.0 / 216.0 * k_r[0] + 1408.0 / 2565.0 * k_r[2] + 2197.0 / 4101.0 * k_r[3] - 1.0 / 5.0 * k_r[4];
	double v_radial_4 = v_radial + 25.0 / 216.0 * k_v[0] + 1408.0 / 2565.0 * k_v[2] + 2197.0 / 4101.0 * k_v[3] - 1.0 / 5.0 * k_v[4];
	double phi_4	  = phi + 25.0 / 216.0 * k_p[0] + 1408.0 / 2565.0 * k_p[2] + 2197.0 / 4101.0 * k_p[3] - 1.0 / 5.0 * k_p[4];

	double radius_5	  = radius + 16.0 / 135.0 * k_r[0] + 6656.0 / 12825.0 * k_r[2] + 28561.0 / 56430.0 * k_r[3] - 9.0 / 50.0 * k_r[4] + 2.0 / 55.0 * k_r[5];
	double v_radial_5 = v_radial + 16.0 / 135.0 * k_v[0] + 6656.0 / 12825.0 * k_v[2] + 28561.0 / 56430.0 * k_v[3] - 9.0 / 50.0 * k_v[4] + 2.0 / 55.0 * k_v[5];
	double phi_5	  = phi + 16.0 / 135.0 * k_p[0] + 6656.0 / 12825.0 * k_p[2] + 28561.0 / 56430.0 * k_p[3] - 9.0 / 50.0 * k_p[4] + 2.0 / 55.0 * k_p[5];

	// Error and adapting the time step (B: 栈数组替代堆分配vector)
	double errors[3] = {fabs(radius_5 - radius_4), fabs(v_radial_5 - v_radial_4), fabs(phi_5 - phi_4)};
	double delta = 1e30;
	for(int i = 0; i < 3; i++)
	{
		double d = 0.84 * pow(error_tolerances[i] / errors[i], 0.25);
		if(d < delta) delta = d;
	}
	double time_step_new = delta * time_step;

	// Check if errors fall below the tolerance
	// 自适应步长方法：根据误差调整 time_step
	if(errors[0] < error_tolerances[0] && errors[1] < error_tolerances[1] && errors[2] < error_tolerances[2])
	{
		time	  = time + time_step;
		radius	  = radius_4;
		// 边界检查：防止半径变为负数（数值误差）
		if(radius < 0.0)
			radius = 0.0;
		v_radial  = v_radial_4;
		phi		  = phi_4;
		time_step = time_step_new;
	}
	else
	{
		time_step = time_step_new;
		Runge_Kutta_45_Step(mass);
	}
	// 旧方法：固定步长，不根据误差调整（已注释）
	// time	  = time + time_step;
	// radius	  = radius_4;
	// if(radius < 0.0)
	// 	radius = 0.0;
	// v_radial  = v_radial_4;
	// phi		  = phi_4;
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

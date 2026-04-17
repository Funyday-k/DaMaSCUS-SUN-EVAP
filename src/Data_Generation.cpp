#include "Data_Generation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <sys/stat.h>
#include <vector>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/Astronomy.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

Simulation_Data::Simulation_Data(unsigned int sample_size, unsigned int max_trajectories, double u_min, unsigned int iso_rings)
: minimum_speed_threshold(u_min), isoreflection_rings(iso_rings),
  number_of_trajectories(0), number_of_free_particles(0), number_of_reflected_particles(0), number_of_captured_particles(0),
  average_number_of_scatterings(0.0), computing_time(0.0), early_stopped(false),
  number_of_data_points(std::vector<unsigned long int>(iso_rings, 0)),
  data(iso_rings, std::vector<libphysica::DataPoint>())
{
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    target_captured_per_rank = (sample_size + mpi_processes - 1) / mpi_processes;  // ceil division
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

void Simulation_Data::Configure(double initial_radius, unsigned int min_scattering, long int max_scattering, unsigned long int max_free_steps)
{
	initial_and_final_radius      = initial_radius;
	minimum_number_of_scatterings = min_scattering;
	maximum_number_of_scatterings = max_scattering;
	maximum_free_time_steps       = max_free_steps;
}

void Simulation_Data::Generate_Data(obscura::DM_Particle& DM, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, SnapshotConfig snapshot_cfg, unsigned int fixed_seed)
{
	auto time_start = std::chrono::system_clock::now();

	// Configure the simulator
	Trajectory_Simulator simulator(solar_model, maximum_free_time_steps, maximum_number_of_scatterings, initial_and_final_radius);
	if(fixed_seed != 0)
		simulator.Fix_PRNG_Seed(fixed_seed);

	// Snapshot configuration
	if(snapshot_cfg.enabled)
	{
		simulator.snapshot_config = snapshot_cfg;
		std::string base_path = g_top_level_dir + "results_" + std::to_string(log10(In_Units(DM.mass, GeV))) + "_" + std::to_string(log10(In_Units(DM.Sigma_Proton(), cm * cm))) + "/";
		struct stat info;
		if(stat(base_path.c_str(), &info) != 0)
			mkdir(base_path.c_str(), 0755);
		simulator.snapshot_output_dir = base_path;
	}

	unsigned long int local_captured = 0;
	unsigned long int local_total = 0;
	early_stopped = false;

	while(local_captured < target_captured_per_rank && local_total < max_trajectories_per_rank)
	{
		Event IC = Initial_Conditions(halo_model, solar_model, simulator.PRNG);
		Hyperbolic_Kepler_Shift(IC, initial_and_final_radius);
		Trajectory_Result trajectory = simulator.Simulate(IC, DM, mpi_rank);

		local_total++;
		number_of_trajectories++;
		average_number_of_scatterings = 1.0 / number_of_trajectories * ((number_of_trajectories - 1) * average_number_of_scatterings + trajectory.number_of_scatterings);

		if(trajectory.bincount.is_captured)
		{
			number_of_captured_particles++;
			local_captured++;

			// Accumulate captured bincount
			for(int b = 0; b < NUM_BINS; b++)
			{
				captured_dt_hist[b]   += trajectory.bincount.dt_hist[b];
				captured_v2dt_hist[b] += trajectory.bincount.v2dt_hist[b];
				captured_dt_sq_hist[b]   += trajectory.bincount.dt_hist[b] * trajectory.bincount.dt_hist[b];
				captured_v2dt_sq_hist[b] += trajectory.bincount.v2dt_hist[b] * trajectory.bincount.v2dt_hist[b];
			}

			// Record evaporation time
			EvaporationRecord rec;
			rec.trajectory_id = number_of_trajectories;
			rec.t_evap = trajectory.bincount.t_last_negative - trajectory.bincount.t_first_negative;
			rec.truncated = trajectory.bincount.truncated;
			evaporation_records.push_back(rec);
		}
		else
		{
			// Accumulate non-captured bincount
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
				// Keep reflection data for compatibility (if needed)
				Hyperbolic_Kepler_Shift(trajectory.final_event, 1.0 * AU);
				double v_final = trajectory.final_event.Speed();
				if(trajectory.number_of_scatterings >= minimum_number_of_scatterings && v_final > KDE_boundary_correction_factor * minimum_speed_threshold)
				{
					unsigned int isoreflection_ring = (isoreflection_rings == 1) ? 0 : trajectory.final_event.Isoreflection_Ring(obscura::Sun_Velocity(), isoreflection_rings);
					data[isoreflection_ring].push_back(libphysica::DataPoint(v_final));
				}
			}
		}

		// Progress bar (periodic)
		if(local_total % 100 == 0 && mpi_rank == 0)
		{
			double progress = std::min(1.0, 1.0 * local_captured / target_captured_per_rank);
			double time_elapsed = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - time_start).count();
			libphysica::Print_Progress_Bar(progress, 0, 44, time_elapsed);
		}
	}

	if(local_total >= max_trajectories_per_rank && local_captured < target_captured_per_rank)
		early_stopped = true;

	auto time_end  = std::chrono::system_clock::now();
	computing_time = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	if(mpi_rank == 0)
	{
		libphysica::Print_Progress_Bar(1.0, 0, 44, computing_time);
		std::cout << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	Perform_MPI_Reductions();
}

void Simulation_Data::Perform_MPI_Reductions()
{
	average_number_of_scatterings *= number_of_trajectories;
	MPI_Allreduce(MPI_IN_PLACE, &number_of_trajectories, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_free_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_reflected_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_captured_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &average_number_of_scatterings, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	average_number_of_scatterings /= number_of_trajectories;

	// Reduce bincount histograms across all ranks
	MPI_Allreduce(MPI_IN_PLACE, captured_dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, captured_v2dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, not_captured_dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, not_captured_v2dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, captured_dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, captured_v2dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, not_captured_dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, not_captured_v2dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// Reduce early_stopped flag (any rank early stopped => global flag)
	int local_es = early_stopped ? 1 : 0;
	int global_es = 0;
	MPI_Allreduce(&local_es, &global_es, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	early_stopped = (global_es > 0);

	// Gather evaporation records from all ranks
	// First gather counts
	int local_evap_count = evaporation_records.size();
	std::vector<int> evap_counts(mpi_processes);
	MPI_Allgather(&local_evap_count, 1, MPI_INT, evap_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

	int total_evap = std::accumulate(evap_counts.begin(), evap_counts.end(), 0);

	// Pack local evaporation data: (trajectory_id, t_evap, truncated)
	// Use doubles for MPI transfer
	std::vector<double> local_evap_data(local_evap_count * 3);
	for(int i = 0; i < local_evap_count; i++)
	{
		local_evap_data[3*i]     = static_cast<double>(evaporation_records[i].trajectory_id);
		local_evap_data[3*i + 1] = evaporation_records[i].t_evap;
		local_evap_data[3*i + 2] = evaporation_records[i].truncated ? 1.0 : 0.0;
	}

	std::vector<int> recv_counts(mpi_processes), displacements(mpi_processes);
	for(int j = 0; j < mpi_processes; j++)
	{
		recv_counts[j] = evap_counts[j] * 3;
		displacements[j] = (j == 0) ? 0 : displacements[j-1] + recv_counts[j-1];
	}

	std::vector<double> global_evap_data(total_evap * 3);
	MPI_Allgatherv(local_evap_data.data(), local_evap_count * 3, MPI_DOUBLE,
	               global_evap_data.data(), recv_counts.data(), displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);

	// Unpack into evaporation_records
	evaporation_records.clear();
	evaporation_records.resize(total_evap);
	for(int i = 0; i < total_evap; i++)
	{
		evaporation_records[i].trajectory_id = static_cast<unsigned long int>(global_evap_data[3*i]);
		evaporation_records[i].t_evap        = global_evap_data[3*i + 1];
		evaporation_records[i].truncated     = (global_evap_data[3*i + 2] > 0.5);
	}

	MPI_Allreduce(MPI_IN_PLACE, &computing_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

void Simulation_Data::Write_Output_Files(const std::string& output_dir, obscura::DM_Particle& DM)
{
	if(mpi_rank != 0)
		return;

	// Create output directory
	struct stat info;
	if(stat(output_dir.c_str(), &info) != 0)
		mkdir(output_dir.c_str(), 0755);

	double mass_gev = In_Units(DM.mass, GeV);
	double sigma_cm2 = In_Units(DM.Sigma_Proton(), cm * cm);

	auto write_header = [&](std::ofstream& f) {
		f << "# DM_mass_GeV = " << std::scientific << std::setprecision(6) << mass_gev << "\n";
		f << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(6) << sigma_cm2 << "\n";
		f << "# total_trajectories = " << number_of_trajectories << "\n";
		f << "# captured_particles = " << number_of_captured_particles << "\n";
		double p_cap = Capture_Ratio();
		f << "# capture_rate = " << std::fixed << std::setprecision(8) << p_cap << "\n";
		// Capture rate error (Wilson interval, 95% CL)
		{
			double N = static_cast<double>(number_of_trajectories);
			double p = p_cap;
			double z = 1.96;
			double sigma_p = (N > 0) ? sqrt(p * (1.0 - p) / N) : 0.0;
			f << "# capture_rate_err = " << std::fixed << std::setprecision(8) << sigma_p << "\n";
			if(N > 0)
			{
				double denom = 1.0 + z * z / N;
				double center = p + z * z / (2.0 * N);
				double spread = z * sqrt(p * (1.0 - p) / N + z * z / (4.0 * N * N));
				double ci_lower = (center - spread) / denom;
				double ci_upper = (center + spread) / denom;
				if(ci_lower < 0.0) ci_lower = 0.0;
				if(ci_upper > 1.0) ci_upper = 1.0;
				f << "# capture_rate_CI_95_lower = " << std::fixed << std::setprecision(8) << ci_lower << "\n";
				f << "# capture_rate_CI_95_upper = " << std::fixed << std::setprecision(8) << ci_upper << "\n";
			}
		}
		if(early_stopped)
			f << "# EARLY_STOP: max_trajectories reached\n";
	};

	// 1. Captured bincount
	{
		std::ofstream f(output_dir + "/captured_bincount.txt");
		write_header(f);
		f << "# bin_index  Sigma_dt[s]  Sigma_v2dt[km2/s]  err_dt[s]  err_v2dt[km2/s]\n";
		double N_cap = static_cast<double>(number_of_captured_particles);
		for(int b = 0; b < NUM_BINS; b++)
		{
			double err_dt = 0.0, err_v2dt = 0.0;
			if(N_cap > 1)
			{
				double mean_dt = captured_dt_hist[b] / N_cap;
				double var_dt  = captured_dt_sq_hist[b] / N_cap - mean_dt * mean_dt;
				if(var_dt < 0.0) var_dt = 0.0;
				err_dt = sqrt(N_cap * var_dt);

				double mean_v2dt = captured_v2dt_hist[b] / N_cap;
				double var_v2dt  = captured_v2dt_sq_hist[b] / N_cap - mean_v2dt * mean_v2dt;
				if(var_v2dt < 0.0) var_v2dt = 0.0;
				err_v2dt = sqrt(N_cap * var_v2dt);
			}
			f << b << "\t" << std::scientific << std::setprecision(10)
			  << captured_dt_hist[b] << "\t" << captured_v2dt_hist[b]
			  << "\t" << err_dt << "\t" << err_v2dt << "\n";
		}
		f.close();
	}

	// 2. Not-captured bincount
	{
		std::ofstream f(output_dir + "/not_captured_bincount.txt");
		write_header(f);
		f << "# bin_index  Sigma_dt[s]  Sigma_v2dt[km2/s]  err_dt[s]  err_v2dt[km2/s]\n";
		double N_nc = static_cast<double>(number_of_trajectories - number_of_captured_particles);
		for(int b = 0; b < NUM_BINS; b++)
		{
			double err_dt = 0.0, err_v2dt = 0.0;
			if(N_nc > 1)
			{
				double mean_dt = not_captured_dt_hist[b] / N_nc;
				double var_dt  = not_captured_dt_sq_hist[b] / N_nc - mean_dt * mean_dt;
				if(var_dt < 0.0) var_dt = 0.0;
				err_dt = sqrt(N_nc * var_dt);

				double mean_v2dt = not_captured_v2dt_hist[b] / N_nc;
				double var_v2dt  = not_captured_v2dt_sq_hist[b] / N_nc - mean_v2dt * mean_v2dt;
				if(var_v2dt < 0.0) var_v2dt = 0.0;
				err_v2dt = sqrt(N_nc * var_v2dt);
			}
			f << b << "\t" << std::scientific << std::setprecision(10)
			  << not_captured_dt_hist[b] << "\t" << not_captured_v2dt_hist[b]
			  << "\t" << err_dt << "\t" << err_v2dt << "\n";
		}
		f.close();
	}

	// 3. Evaporation summary
	{
		std::ofstream f(output_dir + "/evaporation_summary.txt");
		write_header(f);
		f << "# trajectory_id  t_evap[s]  truncated(0/1)\n";
		for(const auto& rec : evaporation_records)
			f << rec.trajectory_id << "\t" << std::scientific << std::setprecision(10) << rec.t_evap << "\t" << (rec.truncated ? 1 : 0) << "\n";
		f.close();
	}
}

double Simulation_Data::Free_Ratio() const
{
	return 1.0 * number_of_free_particles / number_of_trajectories;
}
double Simulation_Data::Capture_Ratio() const
{
	return 1.0 * number_of_captured_particles / number_of_trajectories;
}
double Simulation_Data::Reflection_Ratio(int isoreflection_ring) const
{
	if(isoreflection_ring < 0)
		return 1.0 * number_of_reflected_particles / number_of_trajectories;
	else
		return 1.0 * data[isoreflection_ring].size() / number_of_trajectories;
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

void Simulation_Data::Print_Summary(unsigned int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR
				  << "Simulation data summary" << std::endl
				  << std::endl
				  << "Results:" << std::endl
				  << "Simulated trajectories:\t\t" << number_of_trajectories << std::endl
				  << "Average # of scatterings:\t" << libphysica::Round(average_number_of_scatterings) << std::endl
				  << "Free particles [%]:\t\t" << libphysica::Round(100.0 * Free_Ratio()) << std::endl
				  << "Reflected particles [%]:\t" << libphysica::Round(100.0 * Reflection_Ratio()) << std::endl
				  << "Captured particles [%]:\t\t" << libphysica::Round(100.0 * Capture_Ratio()) << std::endl
				  << "Captured count:\t\t\t" << number_of_captured_particles << std::endl;

		// Capture rate error (Wilson interval)
		{
			double N = static_cast<double>(number_of_trajectories);
			double p = Capture_Ratio();
			double z = 1.96;
			double sigma_p = (N > 0) ? sqrt(p * (1.0 - p) / N) : 0.0;
			std::cout << "Capture rate error (1σ):\t" << std::fixed << std::setprecision(6) << sigma_p << std::endl;
			if(N > 0)
			{
				double denom = 1.0 + z * z / N;
				double center = p + z * z / (2.0 * N);
				double spread = z * sqrt(p * (1.0 - p) / N + z * z / (4.0 * N * N));
				double ci_lower = (center - spread) / denom;
				double ci_upper = (center + spread) / denom;
				if(ci_lower < 0.0) ci_lower = 0.0;
				if(ci_upper > 1.0) ci_upper = 1.0;
				std::cout << "Capture rate 95% CI:\t\t[" << std::fixed << std::setprecision(6) << ci_lower << ", " << ci_upper << "]" << std::endl;
			}
		}

		if(early_stopped)
			std::cout << "*** EARLY STOP: max_trajectories reached ***" << std::endl;

		// Evaporation time median (non-truncated only)
		std::vector<double> non_truncated_evap;
		for(const auto& rec : evaporation_records)
		{
			if(!rec.truncated)
				non_truncated_evap.push_back(rec.t_evap);
		}
		if(!non_truncated_evap.empty())
		{
			std::sort(non_truncated_evap.begin(), non_truncated_evap.end());
			double median;
			size_t n = non_truncated_evap.size();
			if(n % 2 == 0)
				median = 0.5 * (non_truncated_evap[n/2 - 1] + non_truncated_evap[n/2]);
			else
				median = non_truncated_evap[n/2];
			std::cout << "Evaporation time median [s]:\t" << std::scientific << std::setprecision(4) << median << " (" << non_truncated_evap.size() << " non-truncated)" << std::endl;
		}
		else
		{
			std::cout << "Evaporation time median:\tN/A (all truncated or no captures)" << std::endl;
		}

		std::cout << std::endl
				  << "Trajectory rate [1/s]:\t\t" << libphysica::Round(1.0 * number_of_trajectories / computing_time) << std::endl
				  << "Capture rate [1/s]:\t\t" << libphysica::Round(1.0 * number_of_captured_particles / computing_time) << std::endl
				  << "Simulation time:\t\t" << libphysica::Time_Display(computing_time) << std::endl;

		std::cout << SEPARATOR << std::endl;
	}
}

}	// namespace DaMaSCUS_SUN

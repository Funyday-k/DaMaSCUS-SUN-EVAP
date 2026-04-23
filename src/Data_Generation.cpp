#include "Data_Generation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/Astronomy.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

namespace
{
struct Rank_Snapshot_State
{
	uint64_t run_id = 0;
	int32_t snapshot_index = 0;
	int32_t rank = 0;
	int32_t done = 0;
	int32_t trajectory_in_progress = 0;
	uint64_t local_captured = 0;
	uint64_t local_total = 0;
	uint64_t current_trajectory_id = 0;
	double rank_elapsed_wall_sec = 0.0;
	double current_trajectory_wall_sec = 0.0;
	double current_trajectory_physical_sec = 0.0;
	std::array<double, NUM_BINS> captured_dt_hist{};
	std::array<double, NUM_BINS> captured_v2dt_hist{};
};

bool Path_Exists(const std::string& path)
{
	struct stat info;
	return stat(path.c_str(), &info) == 0;
}

void Ensure_Directory_Exists(const std::string& directory)
{
	struct stat info;
	if(stat(directory.c_str(), &info) != 0)
		mkdir(directory.c_str(), 0755);
}

long long Snapshot_Time_Label_Seconds(int snapshot_index, double interval_seconds)
{
	return static_cast<long long>(std::llround(snapshot_index * interval_seconds));
}

std::string Snapshot_Text_File_Path(const std::string& snapshot_root, int snapshot_index, double interval_seconds)
{
	return snapshot_root + "snapshot_" + std::to_string(Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds)) + "s.txt";
}

std::string Rank_Snapshot_Checkpoint_Path(const std::string& rank_snapshot_dir, int rank, int snapshot_index, double interval_seconds)
{
	return rank_snapshot_dir + "snapshot_" + std::to_string(Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds)) + "s_rank" + std::to_string(rank) + ".bin";
}

std::string Rank_Snapshot_Final_Path(const std::string& rank_snapshot_dir, int rank)
{
	return rank_snapshot_dir + "rank" + std::to_string(rank) + "_final.bin";
}

bool Write_Rank_Snapshot_State(const std::string& path, const Rank_Snapshot_State& state)
{
	std::string tmp_path = path + ".tmp." + std::to_string(getpid()) + "." + std::to_string(state.rank);
	std::ofstream file(tmp_path, std::ios::binary | std::ios::trunc);
	if(!file.is_open())
		return false;

	file.write(reinterpret_cast<const char*>(&state.run_id), sizeof(state.run_id));
	file.write(reinterpret_cast<const char*>(&state.snapshot_index), sizeof(state.snapshot_index));
	file.write(reinterpret_cast<const char*>(&state.rank), sizeof(state.rank));
	file.write(reinterpret_cast<const char*>(&state.done), sizeof(state.done));
	file.write(reinterpret_cast<const char*>(&state.trajectory_in_progress), sizeof(state.trajectory_in_progress));
	file.write(reinterpret_cast<const char*>(&state.local_captured), sizeof(state.local_captured));
	file.write(reinterpret_cast<const char*>(&state.local_total), sizeof(state.local_total));
	file.write(reinterpret_cast<const char*>(&state.current_trajectory_id), sizeof(state.current_trajectory_id));
	file.write(reinterpret_cast<const char*>(&state.rank_elapsed_wall_sec), sizeof(state.rank_elapsed_wall_sec));
	file.write(reinterpret_cast<const char*>(&state.current_trajectory_wall_sec), sizeof(state.current_trajectory_wall_sec));
	file.write(reinterpret_cast<const char*>(&state.current_trajectory_physical_sec), sizeof(state.current_trajectory_physical_sec));
	file.write(reinterpret_cast<const char*>(state.captured_dt_hist.data()), NUM_BINS * sizeof(double));
	file.write(reinterpret_cast<const char*>(state.captured_v2dt_hist.data()), NUM_BINS * sizeof(double));
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

bool Read_Rank_Snapshot_State(const std::string& path, uint64_t expected_run_id, Rank_Snapshot_State& state)
{
	std::ifstream file(path, std::ios::binary);
	if(!file.is_open())
		return false;

	file.read(reinterpret_cast<char*>(&state.run_id), sizeof(state.run_id));
	file.read(reinterpret_cast<char*>(&state.snapshot_index), sizeof(state.snapshot_index));
	file.read(reinterpret_cast<char*>(&state.rank), sizeof(state.rank));
	file.read(reinterpret_cast<char*>(&state.done), sizeof(state.done));
	file.read(reinterpret_cast<char*>(&state.trajectory_in_progress), sizeof(state.trajectory_in_progress));
	file.read(reinterpret_cast<char*>(&state.local_captured), sizeof(state.local_captured));
	file.read(reinterpret_cast<char*>(&state.local_total), sizeof(state.local_total));
	file.read(reinterpret_cast<char*>(&state.current_trajectory_id), sizeof(state.current_trajectory_id));
	file.read(reinterpret_cast<char*>(&state.rank_elapsed_wall_sec), sizeof(state.rank_elapsed_wall_sec));
	file.read(reinterpret_cast<char*>(&state.current_trajectory_wall_sec), sizeof(state.current_trajectory_wall_sec));
	file.read(reinterpret_cast<char*>(&state.current_trajectory_physical_sec), sizeof(state.current_trajectory_physical_sec));
	file.read(reinterpret_cast<char*>(state.captured_dt_hist.data()), NUM_BINS * sizeof(double));
	file.read(reinterpret_cast<char*>(state.captured_v2dt_hist.data()), NUM_BINS * sizeof(double));

	if(!file)
		return false;

	return state.run_id == expected_run_id;
}

std::string Format_Rank_Status(const Rank_Snapshot_State& state)
{
	std::ostringstream stream;
	stream << std::fixed << std::setprecision(3);
	if(state.done)
	{
		stream << "done";
	}
	else
	{
		stream << "running";
	}

	stream << " (captured=" << state.local_captured
	       << ", completed=" << state.local_total
	       << ", rank_wall=" << state.rank_elapsed_wall_sec << "s";

	if(state.trajectory_in_progress)
	{
		stream << ", current_traj=" << state.current_trajectory_id
		       << ", traj_wall=" << state.current_trajectory_wall_sec << "s"
		       << ", traj_physical=" << state.current_trajectory_physical_sec << "s";
	}

	stream << ")";
	return stream.str();
}

bool Try_Write_Merged_Snapshot(const std::string& snapshot_root, const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes, uint64_t run_id)
{
	std::array<double, NUM_BINS> sum_dt{}, sum_v2dt{};
	std::vector<Rank_Snapshot_State> rank_states;
	rank_states.reserve(mpi_processes);
	const long long snapshot_time_label = Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds);

	for(int rank = 0; rank < mpi_processes; rank++)
	{
		Rank_Snapshot_State state;
		const std::string checkpoint_path = Rank_Snapshot_Checkpoint_Path(rank_snapshot_dir, rank, snapshot_index, interval_seconds);
		if(!Read_Rank_Snapshot_State(checkpoint_path, run_id, state))
		{
			const std::string final_path = Rank_Snapshot_Final_Path(rank_snapshot_dir, rank);
			if(!Read_Rank_Snapshot_State(final_path, run_id, state) || !state.done || state.rank_elapsed_wall_sec > snapshot_time_label)
				return false;
		}

		rank_states.push_back(state);
		for(int bin = 0; bin < NUM_BINS; bin++)
		{
			sum_dt[bin] += state.captured_dt_hist[bin];
			sum_v2dt[bin] += state.captured_v2dt_hist[bin];
		}
	}

	const std::string output_path = Snapshot_Text_File_Path(snapshot_root, snapshot_index, interval_seconds);
	const std::string tmp_path = output_path + ".tmp." + std::to_string(getpid()) + "." + std::to_string(snapshot_index);
	std::ofstream file(tmp_path, std::ios::trunc);
	if(!file.is_open())
		return false;

	file << "# Cumulative captured bincount snapshot\n";
	file << "# snapshot_target_wall_time_s = " << snapshot_time_label << "\n";
	file << "# snapshot_interval_s = " << std::fixed << std::setprecision(3) << interval_seconds << "\n";
	file << "# ranks_ready = " << mpi_processes << " / " << mpi_processes << "\n";
	file << "#\n";
	file << "# Rank status:\n";
	for(const Rank_Snapshot_State& state : rank_states)
		file << "#   rank " << state.rank << ": " << Format_Rank_Status(state) << "\n";
	file << "#\n";
	file << "# bin_index  cap_dt[s]  cap_v2dt[km2/s]\n";
	for(int bin = 0; bin < NUM_BINS; bin++)
		file << bin << "\t" << std::scientific << std::setprecision(10) << sum_dt[bin] << "\t" << sum_v2dt[bin] << "\n";
	file.close();

	if(!file.good())
	{
		std::remove(tmp_path.c_str());
		return false;
	}

	if(std::rename(tmp_path.c_str(), output_path.c_str()) != 0)
	{
		std::remove(tmp_path.c_str());
		return false;
	}

	return true;
}
}

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

    // Initialize computation time statistics
    total_wall_time_captured = 0.0;
    total_wall_time_not_captured = 0.0;
    wall_time_hist_captured.fill(0);
    wall_time_hist_not_captured.fill(0);
    wall_time_overflow_captured = 0;
    wall_time_overflow_not_captured = 0;
    total_rk45_steps_captured = 0;
    total_rk45_steps_not_captured = 0;
    step_count_hist_captured.fill(0);
    step_count_hist_not_captured.fill(0);
    step_count_overflow_captured = 0;
    step_count_overflow_not_captured = 0;
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
	unsigned long int local_captured = 0;
	unsigned long int local_total = 0;

	// Configure the simulator
	Trajectory_Simulator simulator(solar_model, maximum_free_time_steps, maximum_number_of_scatterings, initial_and_final_radius);
	if(fixed_seed != 0)
		simulator.Fix_PRNG_Seed(fixed_seed);

	// Snapshot configuration
	const double snapshot_interval = (snapshot_cfg.interval_seconds > 0.0) ? snapshot_cfg.interval_seconds : 60.0;
	std::string snapshot_root;
	std::string rank_snapshot_dir;
	uint64_t snapshot_run_id = 0;
	int next_snapshot_index = 1;
	if(snapshot_cfg.enabled)
	{
		snapshot_root = g_top_level_dir + "results_" + std::to_string(log10(In_Units(DM.mass, GeV))) + "_" + std::to_string(log10(In_Units(DM.Sigma_Proton(), cm * cm))) + "/snapshot/";
		rank_snapshot_dir = snapshot_root + "rank_snapshot/";
		Ensure_Directory_Exists(g_top_level_dir + "results_" + std::to_string(log10(In_Units(DM.mass, GeV))) + "_" + std::to_string(log10(In_Units(DM.Sigma_Proton(), cm * cm))) + "/");
		Ensure_Directory_Exists(snapshot_root);
		Ensure_Directory_Exists(rank_snapshot_dir);
		if(mpi_rank == 0)
			snapshot_run_id = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
		MPI_Bcast(&snapshot_run_id, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	}

	auto elapsed_since_start = [&]()
	{
		return 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - time_start).count();
	};

	auto build_rank_snapshot_state = [&](bool done)
	{
		Rank_Snapshot_State state;
		state.run_id = snapshot_run_id;
		state.rank = mpi_rank;
		state.done = done ? 1 : 0;
		state.trajectory_in_progress = (!done && simulator.Trajectory_In_Progress()) ? 1 : 0;
		state.local_captured = static_cast<uint64_t>(local_captured);
		state.local_total = static_cast<uint64_t>(local_total);
		state.current_trajectory_id = state.trajectory_in_progress ? static_cast<uint64_t>(simulator.Current_Trajectory_ID()) : 0;
		state.rank_elapsed_wall_sec = done ? computing_time : elapsed_since_start();
		state.current_trajectory_wall_sec = state.trajectory_in_progress ? simulator.Current_Trajectory_Wall_Time_Seconds() : 0.0;
		state.current_trajectory_physical_sec = state.trajectory_in_progress ? simulator.Current_Trajectory_Physical_Time_Seconds() : 0.0;
		state.captured_dt_hist = captured_dt_hist;
		state.captured_v2dt_hist = captured_v2dt_hist;
		return state;
	};

	auto try_write_ready_snapshots = [&](int max_snapshot_index)
	{
		for(int snapshot_index = 1; snapshot_index <= max_snapshot_index; snapshot_index++)
			Try_Write_Merged_Snapshot(snapshot_root, rank_snapshot_dir, snapshot_index, snapshot_interval, mpi_processes, snapshot_run_id);
	};

	auto publish_checkpoint_snapshots = [&]()
	{
		if(!snapshot_cfg.enabled)
			return;

		double elapsed = elapsed_since_start();
		while(elapsed >= next_snapshot_index * snapshot_interval)
		{
			Rank_Snapshot_State state = build_rank_snapshot_state(false);
			state.snapshot_index = next_snapshot_index;
			Write_Rank_Snapshot_State(Rank_Snapshot_Checkpoint_Path(rank_snapshot_dir, mpi_rank, next_snapshot_index, snapshot_interval), state);
			try_write_ready_snapshots(next_snapshot_index);
			next_snapshot_index++;
		}
	};
	early_stopped = false;
	int last_progress_milestone = -1;
	if(snapshot_cfg.enabled)
		simulator.Set_Snapshot_Progress_Callback([&](const Trajectory_Simulator&)
		{
			publish_checkpoint_snapshots();
		});

	while(local_captured < target_captured_per_rank && local_total < max_trajectories_per_rank)
	{
		Event IC = Initial_Conditions(halo_model, solar_model, simulator.PRNG);
		Hyperbolic_Kepler_Shift(IC, initial_and_final_radius);

		auto traj_t0 = std::chrono::high_resolution_clock::now();
		Trajectory_Result trajectory = simulator.Simulate(IC, DM, mpi_rank);
		auto traj_t1 = std::chrono::high_resolution_clock::now();
		double traj_wall_sec = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(traj_t1 - traj_t0).count();

		local_total++;
		number_of_trajectories++;
		average_number_of_scatterings = 1.0 / number_of_trajectories * ((number_of_trajectories - 1) * average_number_of_scatterings + trajectory.number_of_scatterings);

		// Helper lambdas for log-histogram binning
		auto wall_time_bin = [](double t) -> int {
			if(t <= 0.0) return -1;
			return static_cast<int>((log10(t) - WALL_TIME_LOG_MIN) * 10.0);
		};
		auto step_count_bin = [](unsigned long int s) -> int {
			if(s == 0) return -1;
			return static_cast<int>((log10(static_cast<double>(s)) - STEP_COUNT_LOG_MIN) * 10.0);
		};

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

			// Computation time & step count statistics (captured)
			total_wall_time_captured += traj_wall_sec;
			total_rk45_steps_captured += trajectory.total_rk45_steps;
			int wb = wall_time_bin(traj_wall_sec);
			if(wb >= 0 && wb < WALL_TIME_BINS) wall_time_hist_captured[wb]++;
			else if(wb >= WALL_TIME_BINS) wall_time_overflow_captured++;
			int sb = step_count_bin(trajectory.total_rk45_steps);
			if(sb >= 0 && sb < STEP_COUNT_BINS) step_count_hist_captured[sb]++;
			else if(sb >= STEP_COUNT_BINS) step_count_overflow_captured++;
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

			// Computation time & step count statistics (not captured)
			total_wall_time_not_captured += traj_wall_sec;
			total_rk45_steps_not_captured += trajectory.total_rk45_steps;
			int wb = wall_time_bin(traj_wall_sec);
			if(wb >= 0 && wb < WALL_TIME_BINS) wall_time_hist_not_captured[wb]++;
			else if(wb >= WALL_TIME_BINS) wall_time_overflow_not_captured++;
			int sb = step_count_bin(trajectory.total_rk45_steps);
			if(sb >= 0 && sb < STEP_COUNT_BINS) step_count_hist_not_captured[sb]++;
			else if(sb >= STEP_COUNT_BINS) step_count_overflow_not_captured++;

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

		// Progress bar (every 20%)
		if(mpi_rank == 0)
		{
			double progress = std::min(1.0, 1.0 * local_captured / target_captured_per_rank);
			int milestone = static_cast<int>(progress * 5);  // 0=0%,1=20%,...,5=100%
			if(milestone > last_progress_milestone)
			{
				last_progress_milestone = milestone;
				double time_elapsed = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - time_start).count();
				libphysica::Print_Progress_Bar(progress, 0, 44, time_elapsed);
			}
		}

		publish_checkpoint_snapshots();
	}

	if(local_total >= max_trajectories_per_rank && local_captured < target_captured_per_rank)
		early_stopped = true;

	auto time_end  = std::chrono::system_clock::now();
	computing_time = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();

	publish_checkpoint_snapshots();

	// Final snapshot state: keep per-rank binary under snapshot/rank_snapshot and merge on interval boundaries.
	if(snapshot_cfg.enabled)
	{
		Rank_Snapshot_State final_state = build_rank_snapshot_state(true);
		final_state.snapshot_index = next_snapshot_index - 1;
		Write_Rank_Snapshot_State(Rank_Snapshot_Final_Path(rank_snapshot_dir, mpi_rank), final_state);

		int final_snapshot_index = std::max(1, static_cast<int>(std::ceil(computing_time / snapshot_interval)));
		try_write_ready_snapshots(final_snapshot_index);
	}

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

	// Reduce computation time & step count statistics
	MPI_Allreduce(MPI_IN_PLACE, &total_wall_time_captured, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &total_wall_time_not_captured, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, wall_time_hist_captured.data(), WALL_TIME_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, wall_time_hist_not_captured.data(), WALL_TIME_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &wall_time_overflow_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &wall_time_overflow_not_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &total_rk45_steps_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &total_rk45_steps_not_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, step_count_hist_captured.data(), STEP_COUNT_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, step_count_hist_not_captured.data(), STEP_COUNT_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &step_count_overflow_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &step_count_overflow_not_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
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

	// 4. Computation time & RK45 step count summary
	{
		std::ofstream f(output_dir + "/computation_time_summary.txt");
		write_header(f);

		double mean_wall_cap = (number_of_captured_particles > 0) ? total_wall_time_captured / number_of_captured_particles : 0.0;
		double mean_wall_nc  = (number_of_trajectories > number_of_captured_particles) ? total_wall_time_not_captured / (number_of_trajectories - number_of_captured_particles) : 0.0;
		double mean_steps_cap = (number_of_captured_particles > 0) ? static_cast<double>(total_rk45_steps_captured) / number_of_captured_particles : 0.0;
		double mean_steps_nc  = (number_of_trajectories > number_of_captured_particles) ? static_cast<double>(total_rk45_steps_not_captured) / (number_of_trajectories - number_of_captured_particles) : 0.0;

		f << "# total_wall_time_captured_sec = " << std::scientific << std::setprecision(6) << total_wall_time_captured << "\n";
		f << "# total_wall_time_not_captured_sec = " << std::scientific << std::setprecision(6) << total_wall_time_not_captured << "\n";
		f << "# mean_wall_time_captured_sec = " << std::scientific << std::setprecision(6) << mean_wall_cap << "\n";
		f << "# mean_wall_time_not_captured_sec = " << std::scientific << std::setprecision(6) << mean_wall_nc << "\n";
		f << "# total_rk45_steps_captured = " << total_rk45_steps_captured << "\n";
		f << "# total_rk45_steps_not_captured = " << total_rk45_steps_not_captured << "\n";
		f << "# mean_rk45_steps_captured = " << std::scientific << std::setprecision(6) << mean_steps_cap << "\n";
		f << "# mean_rk45_steps_not_captured = " << std::scientific << std::setprecision(6) << mean_steps_nc << "\n";
		f << "#\n";

		// Wall-clock time histogram
		f << "# [Wall-clock time histogram]\n";
		f << "# bin_index  log10_t_lower  log10_t_upper  count_captured  count_not_captured\n";
		for(int b = 0; b < WALL_TIME_BINS; b++)
		{
			double lo = WALL_TIME_LOG_MIN + b * 0.1;
			double hi = lo + 0.1;
			f << b << "\t" << std::fixed << std::setprecision(1) << lo << "\t" << hi
			  << "\t" << wall_time_hist_captured[b] << "\t" << wall_time_hist_not_captured[b] << "\n";
		}
		f << "# overflow_captured = " << wall_time_overflow_captured << "\n";
		f << "# overflow_not_captured = " << wall_time_overflow_not_captured << "\n";
		f << "#\n";

		// RK45 step count histogram
		f << "# [RK45 step count histogram]\n";
		f << "# bin_index  log10_steps_lower  log10_steps_upper  count_captured  count_not_captured\n";
		for(int b = 0; b < STEP_COUNT_BINS; b++)
		{
			double lo = STEP_COUNT_LOG_MIN + b * 0.1;
			double hi = lo + 0.1;
			f << b << "\t" << std::fixed << std::setprecision(1) << lo << "\t" << hi
			  << "\t" << step_count_hist_captured[b] << "\t" << step_count_hist_not_captured[b] << "\n";
		}
		f << "# overflow_captured = " << step_count_overflow_captured << "\n";
		f << "# overflow_not_captured = " << step_count_overflow_not_captured << "\n";
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

		// Per-trajectory computation time & RK45 step count summary
		std::cout << std::endl
				  << "Captured wall-clock time:\t" << libphysica::Time_Display(total_wall_time_captured) << std::endl
				  << "Not-captured wall-clock time:\t" << libphysica::Time_Display(total_wall_time_not_captured) << std::endl;
		if(number_of_captured_particles > 0)
			std::cout << "Mean captured traj time:\t\t" << std::scientific << std::setprecision(4) << total_wall_time_captured / number_of_captured_particles << " s" << std::endl;
		if(number_of_trajectories > number_of_captured_particles)
			std::cout << "Mean not-captured traj time:\t" << std::scientific << std::setprecision(4) << total_wall_time_not_captured / (number_of_trajectories - number_of_captured_particles) << " s" << std::endl;
		std::cout << "Total RK45 steps (captured):\t" << total_rk45_steps_captured << std::endl
				  << "Total RK45 steps (not captured):\t" << total_rk45_steps_not_captured << std::endl;
		if(number_of_captured_particles > 0)
			std::cout << "Mean RK45 steps/traj (captured):\t" << std::scientific << std::setprecision(4) << static_cast<double>(total_rk45_steps_captured) / number_of_captured_particles << std::endl;
		if(number_of_trajectories > number_of_captured_particles)
			std::cout << "Mean RK45 steps/traj (not cap):\t" << std::scientific << std::setprecision(4) << static_cast<double>(total_rk45_steps_not_captured) / (number_of_trajectories - number_of_captured_particles) << std::endl;

		std::cout << SEPARATOR << std::endl;
	}
}

}	// namespace DaMaSCUS_SUN

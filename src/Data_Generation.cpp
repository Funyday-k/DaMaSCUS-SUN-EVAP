#include "Data_Generation.hpp"

#include <dirent.h>
#include <errno.h>
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
constexpr int SNAPSHOT_WALL_TIME_BINS = 120;
constexpr double SNAPSHOT_WALL_TIME_LOG_MIN = -4.0;
constexpr int SNAPSHOT_STEP_COUNT_BINS = 140;
constexpr double SNAPSHOT_STEP_COUNT_LOG_MIN = 0.0;

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
	std::array<double, NUM_BINS> captured_dt_sq_hist{};
	std::array<double, NUM_BINS> captured_v2dt_sq_hist{};
	std::array<double, NUM_BINS> not_captured_dt_hist{};
	std::array<double, NUM_BINS> not_captured_v2dt_hist{};
	std::array<double, NUM_BINS> not_captured_dt_sq_hist{};
	std::array<double, NUM_BINS> not_captured_v2dt_sq_hist{};
	double total_wall_time_captured = 0.0;
	double total_wall_time_not_captured = 0.0;
	uint64_t total_rk45_steps_captured = 0;
	uint64_t total_rk45_steps_not_captured = 0;
	std::array<uint64_t, SNAPSHOT_WALL_TIME_BINS> wall_time_hist_captured{};
	std::array<uint64_t, SNAPSHOT_WALL_TIME_BINS> wall_time_hist_not_captured{};
	uint64_t wall_time_overflow_captured = 0;
	uint64_t wall_time_overflow_not_captured = 0;
	std::array<uint64_t, SNAPSHOT_STEP_COUNT_BINS> step_count_hist_captured{};
	std::array<uint64_t, SNAPSHOT_STEP_COUNT_BINS> step_count_hist_not_captured{};
	uint64_t step_count_overflow_captured = 0;
	uint64_t step_count_overflow_not_captured = 0;
	std::vector<EvaporationRecord> evaporation_records;
};

struct Snapshot_Report_State
{
	long long snapshot_time_label = 0;
	double snapshot_interval_seconds = 0.0;
	uint64_t total_trajectories = 0;
	uint64_t captured_particles = 0;
	std::array<double, NUM_BINS> captured_dt_hist{};
	std::array<double, NUM_BINS> captured_v2dt_hist{};
	std::array<double, NUM_BINS> captured_dt_sq_hist{};
	std::array<double, NUM_BINS> captured_v2dt_sq_hist{};
	std::array<double, NUM_BINS> not_captured_dt_hist{};
	std::array<double, NUM_BINS> not_captured_v2dt_hist{};
	std::array<double, NUM_BINS> not_captured_dt_sq_hist{};
	std::array<double, NUM_BINS> not_captured_v2dt_sq_hist{};
	double total_wall_time_captured = 0.0;
	double total_wall_time_not_captured = 0.0;
	uint64_t total_rk45_steps_captured = 0;
	uint64_t total_rk45_steps_not_captured = 0;
	std::array<uint64_t, SNAPSHOT_WALL_TIME_BINS> wall_time_hist_captured{};
	std::array<uint64_t, SNAPSHOT_WALL_TIME_BINS> wall_time_hist_not_captured{};
	uint64_t wall_time_overflow_captured = 0;
	uint64_t wall_time_overflow_not_captured = 0;
	std::array<uint64_t, SNAPSHOT_STEP_COUNT_BINS> step_count_hist_captured{};
	std::array<uint64_t, SNAPSHOT_STEP_COUNT_BINS> step_count_hist_not_captured{};
	uint64_t step_count_overflow_captured = 0;
	uint64_t step_count_overflow_not_captured = 0;
	std::vector<EvaporationRecord> evaporation_records;
	std::vector<Rank_Snapshot_State> rank_states;
};

struct Snapshot_Load_Diagnostics
{
	int ready_ranks = 0;
	std::vector<std::string> per_rank_status;
};

std::string Format_Rank_Status(const Rank_Snapshot_State& state);

template<typename T, size_t N>
bool Find_Nonzero_Bin_Range(const std::array<T, N>& captured_hist, const std::array<T, N>& not_captured_hist, int& first_bin, int& last_bin)
{
	first_bin = 0;
	while(first_bin < static_cast<int>(N) && captured_hist[first_bin] == 0 && not_captured_hist[first_bin] == 0)
		first_bin++;

	last_bin = static_cast<int>(N) - 1;
	while(last_bin >= first_bin && captured_hist[last_bin] == 0 && not_captured_hist[last_bin] == 0)
		last_bin--;

	return first_bin <= last_bin;
}

template<typename T>
void Write_Binary_Value(std::ofstream& file, const T& value)
{
	file.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<typename T, size_t N>
void Write_Binary_Array(std::ofstream& file, const std::array<T, N>& values)
{
	file.write(reinterpret_cast<const char*>(values.data()), N * sizeof(T));
}

template<typename T>
void Read_Binary_Value(std::ifstream& file, T& value)
{
	file.read(reinterpret_cast<char*>(&value), sizeof(T));
}

template<typename T, size_t N>
void Read_Binary_Array(std::ifstream& file, std::array<T, N>& values)
{
	file.read(reinterpret_cast<char*>(values.data()), N * sizeof(T));
}

double Capture_Ratio_From_Counts(uint64_t total_trajectories, uint64_t captured_particles)
{
	return (total_trajectories > 0) ? static_cast<double>(captured_particles) / static_cast<double>(total_trajectories) : 0.0;
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

void Write_Report_Header(std::ofstream& file, double mass_gev, double sigma_cm2, uint64_t total_trajectories, uint64_t captured_particles, bool early_stopped, long long snapshot_time_label = -1, double snapshot_interval_seconds = 0.0)
{
	if(snapshot_time_label >= 0)
	{
		file << "# snapshot_target_wall_time_s = " << snapshot_time_label << "\n";
		file << "# snapshot_interval_s = " << std::fixed << std::setprecision(3) << snapshot_interval_seconds << "\n";
	}

	file << "# DM_mass_GeV = " << std::scientific << std::setprecision(6) << mass_gev << "\n";
	file << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(6) << sigma_cm2 << "\n";
	file << "# total_trajectories = " << total_trajectories << "\n";
	file << "# captured_particles = " << captured_particles << "\n";

	double p_cap = Capture_Ratio_From_Counts(total_trajectories, captured_particles);
	file << "# capture_rate = " << std::fixed << std::setprecision(8) << p_cap << "\n";

	{
		double N = static_cast<double>(total_trajectories);
		double p = p_cap;
		double z = 1.96;
		double sigma_p = (N > 0.0) ? sqrt(p * (1.0 - p) / N) : 0.0;
		file << "# capture_rate_err = " << std::fixed << std::setprecision(8) << sigma_p << "\n";
		if(N > 0.0)
		{
			double denom = 1.0 + z * z / N;
			double center = p + z * z / (2.0 * N);
			double spread = z * sqrt(p * (1.0 - p) / N + z * z / (4.0 * N * N));
			double ci_lower = (center - spread) / denom;
			double ci_upper = (center + spread) / denom;
			if(ci_lower < 0.0) ci_lower = 0.0;
			if(ci_upper > 1.0) ci_upper = 1.0;
			file << "# capture_rate_CI_95_lower = " << std::fixed << std::setprecision(8) << ci_lower << "\n";
			file << "# capture_rate_CI_95_upper = " << std::fixed << std::setprecision(8) << ci_upper << "\n";
		}
	}

	if(early_stopped)
		file << "# EARLY_STOP: max_trajectories reached\n";
}

std::string Format_Physical_Time_Scientific(double physical_time_sec)
{
	std::ostringstream stream;

	if(!std::isfinite(physical_time_sec))
		stream << "nan";
	else
		stream << std::uppercase << std::scientific << std::setprecision(3) << physical_time_sec;

	return stream.str();
}

bool Path_Exists(const std::string& path)
{
	struct stat info;
	return stat(path.c_str(), &info) == 0;
}

std::string Join_Path(const std::string& directory, const std::string& name)
{
	if(directory.empty() || directory.back() == '/')
		return directory + name;
	return directory + "/" + name;
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

std::string Snapshot_Status_File_Path(const std::string& snapshot_root, int snapshot_index, double interval_seconds)
{
	return snapshot_root + "snapshot_" + std::to_string(Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds)) + "s_status.txt";
}

std::string Snapshot_Report_File_Path(const std::string& snapshot_root, int snapshot_index, double interval_seconds, const std::string& suffix)
{
	return snapshot_root + "snapshot_" + std::to_string(Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds)) + "s_" + suffix + ".txt";
}

std::string Rank_Snapshot_Checkpoint_Path(const std::string& rank_snapshot_dir, int rank, int snapshot_index, double interval_seconds)
{
	return rank_snapshot_dir + "snapshot_" + std::to_string(Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds)) + "s_rank" + std::to_string(rank) + ".bin";
}

std::string Rank_Snapshot_Final_Path(const std::string& rank_snapshot_dir, int rank)
{
	return rank_snapshot_dir + "rank" + std::to_string(rank) + "_final.bin";
}

void Cleanup_Snapshot_Checkpoints(const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes)
{
	for(int rank = 0; rank < mpi_processes; rank++)
		std::remove(Rank_Snapshot_Checkpoint_Path(rank_snapshot_dir, rank, snapshot_index, interval_seconds).c_str());
}

void Cleanup_Final_Snapshot_States(const std::string& rank_snapshot_dir, int mpi_processes)
{
	for(int rank = 0; rank < mpi_processes; rank++)
		std::remove(Rank_Snapshot_Final_Path(rank_snapshot_dir, rank).c_str());

	rmdir(rank_snapshot_dir.c_str());
}

bool Write_Text_File_Atomically(const std::string& path, int unique_tag, const std::function<void(std::ofstream&)>& writer)
{
	std::string tmp_path = path + ".tmp." + std::to_string(getpid()) + "." + std::to_string(unique_tag);
	std::ofstream file(tmp_path, std::ios::trunc);
	if(!file.is_open())
		return false;

	writer(file);
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

bool Write_Rank_Snapshot_State(const std::string& path, const Rank_Snapshot_State& state)
{
	std::string tmp_path = path + ".tmp." + std::to_string(getpid()) + "." + std::to_string(state.rank);
	std::ofstream file(tmp_path, std::ios::binary | std::ios::trunc);
	if(!file.is_open())
		return false;

	Write_Binary_Value(file, state.run_id);
	Write_Binary_Value(file, state.snapshot_index);
	Write_Binary_Value(file, state.rank);
	Write_Binary_Value(file, state.done);
	Write_Binary_Value(file, state.trajectory_in_progress);
	Write_Binary_Value(file, state.local_captured);
	Write_Binary_Value(file, state.local_total);
	Write_Binary_Value(file, state.current_trajectory_id);
	Write_Binary_Value(file, state.rank_elapsed_wall_sec);
	Write_Binary_Value(file, state.current_trajectory_wall_sec);
	Write_Binary_Value(file, state.current_trajectory_physical_sec);
	Write_Binary_Array(file, state.captured_dt_hist);
	Write_Binary_Array(file, state.captured_v2dt_hist);
	Write_Binary_Array(file, state.captured_dt_sq_hist);
	Write_Binary_Array(file, state.captured_v2dt_sq_hist);
	Write_Binary_Array(file, state.not_captured_dt_hist);
	Write_Binary_Array(file, state.not_captured_v2dt_hist);
	Write_Binary_Array(file, state.not_captured_dt_sq_hist);
	Write_Binary_Array(file, state.not_captured_v2dt_sq_hist);
	Write_Binary_Value(file, state.total_wall_time_captured);
	Write_Binary_Value(file, state.total_wall_time_not_captured);
	Write_Binary_Value(file, state.total_rk45_steps_captured);
	Write_Binary_Value(file, state.total_rk45_steps_not_captured);
	Write_Binary_Array(file, state.wall_time_hist_captured);
	Write_Binary_Array(file, state.wall_time_hist_not_captured);
	Write_Binary_Value(file, state.wall_time_overflow_captured);
	Write_Binary_Value(file, state.wall_time_overflow_not_captured);
	Write_Binary_Array(file, state.step_count_hist_captured);
	Write_Binary_Array(file, state.step_count_hist_not_captured);
	Write_Binary_Value(file, state.step_count_overflow_captured);
	Write_Binary_Value(file, state.step_count_overflow_not_captured);

	uint64_t evaporation_count = static_cast<uint64_t>(state.evaporation_records.size());
	Write_Binary_Value(file, evaporation_count);
	for(const auto& rec : state.evaporation_records)
	{
		uint64_t trajectory_id = static_cast<uint64_t>(rec.trajectory_id);
		uint8_t truncated = rec.truncated ? 1 : 0;
		Write_Binary_Value(file, trajectory_id);
		Write_Binary_Value(file, rec.t_evap);
		Write_Binary_Value(file, truncated);
	}
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

	Read_Binary_Value(file, state.run_id);
	Read_Binary_Value(file, state.snapshot_index);
	Read_Binary_Value(file, state.rank);
	Read_Binary_Value(file, state.done);
	Read_Binary_Value(file, state.trajectory_in_progress);
	Read_Binary_Value(file, state.local_captured);
	Read_Binary_Value(file, state.local_total);
	Read_Binary_Value(file, state.current_trajectory_id);
	Read_Binary_Value(file, state.rank_elapsed_wall_sec);
	Read_Binary_Value(file, state.current_trajectory_wall_sec);
	Read_Binary_Value(file, state.current_trajectory_physical_sec);
	Read_Binary_Array(file, state.captured_dt_hist);
	Read_Binary_Array(file, state.captured_v2dt_hist);
	Read_Binary_Array(file, state.captured_dt_sq_hist);
	Read_Binary_Array(file, state.captured_v2dt_sq_hist);
	Read_Binary_Array(file, state.not_captured_dt_hist);
	Read_Binary_Array(file, state.not_captured_v2dt_hist);
	Read_Binary_Array(file, state.not_captured_dt_sq_hist);
	Read_Binary_Array(file, state.not_captured_v2dt_sq_hist);
	Read_Binary_Value(file, state.total_wall_time_captured);
	Read_Binary_Value(file, state.total_wall_time_not_captured);
	Read_Binary_Value(file, state.total_rk45_steps_captured);
	Read_Binary_Value(file, state.total_rk45_steps_not_captured);
	Read_Binary_Array(file, state.wall_time_hist_captured);
	Read_Binary_Array(file, state.wall_time_hist_not_captured);
	Read_Binary_Value(file, state.wall_time_overflow_captured);
	Read_Binary_Value(file, state.wall_time_overflow_not_captured);
	Read_Binary_Array(file, state.step_count_hist_captured);
	Read_Binary_Array(file, state.step_count_hist_not_captured);
	Read_Binary_Value(file, state.step_count_overflow_captured);
	Read_Binary_Value(file, state.step_count_overflow_not_captured);

	uint64_t evaporation_count = 0;
	Read_Binary_Value(file, evaporation_count);
	state.evaporation_records.clear();
	state.evaporation_records.reserve(static_cast<size_t>(evaporation_count));
	for(uint64_t i = 0; i < evaporation_count; i++)
	{
		uint64_t trajectory_id = 0;
		double t_evap = 0.0;
		uint8_t truncated = 0;
		Read_Binary_Value(file, trajectory_id);
		Read_Binary_Value(file, t_evap);
		Read_Binary_Value(file, truncated);
		EvaporationRecord rec;
		rec.trajectory_id = static_cast<unsigned long int>(trajectory_id);
		rec.t_evap = t_evap;
		rec.truncated = (truncated != 0);
		state.evaporation_records.push_back(rec);
	}

	if(!file)
		return false;

	return state.run_id == expected_run_id;
}

bool Write_Snapshot_Status_File(const std::string& snapshot_root, int snapshot_index, double interval_seconds, int caller_rank, const Snapshot_Load_Diagnostics& diagnostics, bool merged)
{
	return Write_Text_File_Atomically(Snapshot_Status_File_Path(snapshot_root, snapshot_index, interval_seconds), snapshot_index * 1000 + caller_rank, [&](std::ofstream& file)
	{
		file << "# snapshot_status = " << (merged ? "merged" : "waiting") << "\n";
		file << "# snapshot_target_wall_time_s = " << Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds) << "\n";
		file << "# snapshot_interval_s = " << std::fixed << std::setprecision(3) << interval_seconds << "\n";
		file << "# attempted_by_rank = " << caller_rank << "\n";
		file << "# ready_ranks = " << diagnostics.ready_ranks << " / " << diagnostics.per_rank_status.size() << "\n";
		file << "#\n";
		file << "# rank_index  source_or_wait_reason\n";
		for(size_t rank = 0; rank < diagnostics.per_rank_status.size(); rank++)
			file << rank << "\t" << diagnostics.per_rank_status[rank] << "\n";
	});
}

void Accumulate_Snapshot_Report_State(Snapshot_Report_State& report, const Rank_Snapshot_State& state)
{
	report.total_trajectories += state.local_total;
	report.captured_particles += state.local_captured;
	report.total_wall_time_captured += state.total_wall_time_captured;
	report.total_wall_time_not_captured += state.total_wall_time_not_captured;
	report.total_rk45_steps_captured += state.total_rk45_steps_captured;
	report.total_rk45_steps_not_captured += state.total_rk45_steps_not_captured;
	report.wall_time_overflow_captured += state.wall_time_overflow_captured;
	report.wall_time_overflow_not_captured += state.wall_time_overflow_not_captured;
	report.step_count_overflow_captured += state.step_count_overflow_captured;
	report.step_count_overflow_not_captured += state.step_count_overflow_not_captured;

	for(int bin = 0; bin < NUM_BINS; bin++)
	{
		report.captured_dt_hist[bin] += state.captured_dt_hist[bin];
		report.captured_v2dt_hist[bin] += state.captured_v2dt_hist[bin];
		report.captured_dt_sq_hist[bin] += state.captured_dt_sq_hist[bin];
		report.captured_v2dt_sq_hist[bin] += state.captured_v2dt_sq_hist[bin];
		report.not_captured_dt_hist[bin] += state.not_captured_dt_hist[bin];
		report.not_captured_v2dt_hist[bin] += state.not_captured_v2dt_hist[bin];
		report.not_captured_dt_sq_hist[bin] += state.not_captured_dt_sq_hist[bin];
		report.not_captured_v2dt_sq_hist[bin] += state.not_captured_v2dt_sq_hist[bin];
	}

	for(int bin = 0; bin < SNAPSHOT_WALL_TIME_BINS; bin++)
	{
		report.wall_time_hist_captured[bin] += state.wall_time_hist_captured[bin];
		report.wall_time_hist_not_captured[bin] += state.wall_time_hist_not_captured[bin];
	}

	for(int bin = 0; bin < SNAPSHOT_STEP_COUNT_BINS; bin++)
	{
		report.step_count_hist_captured[bin] += state.step_count_hist_captured[bin];
		report.step_count_hist_not_captured[bin] += state.step_count_hist_not_captured[bin];
	}

	report.evaporation_records.insert(report.evaporation_records.end(), state.evaporation_records.begin(), state.evaporation_records.end());
}

bool Load_Snapshot_Report_State(const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes, uint64_t run_id, Snapshot_Report_State& report, Snapshot_Load_Diagnostics* diagnostics = NULL)
{
	report = Snapshot_Report_State();
	report.snapshot_time_label = Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds);
	report.snapshot_interval_seconds = interval_seconds;
	report.rank_states.reserve(mpi_processes);
	if(diagnostics != NULL)
	{
		diagnostics->ready_ranks = 0;
		diagnostics->per_rank_status.assign(mpi_processes, "pending");
	}

	bool all_ranks_ready = true;

	for(int rank = 0; rank < mpi_processes; rank++)
	{
		Rank_Snapshot_State state;
		const std::string checkpoint_path = Rank_Snapshot_Checkpoint_Path(rank_snapshot_dir, rank, snapshot_index, interval_seconds);
		const bool checkpoint_exists = Path_Exists(checkpoint_path);
		if(Read_Rank_Snapshot_State(checkpoint_path, run_id, state))
		{
			Accumulate_Snapshot_Report_State(report, state);
			report.rank_states.push_back(state);
			if(diagnostics != NULL)
			{
				diagnostics->ready_ranks++;
				diagnostics->per_rank_status[rank] = "checkpoint";
			}
			continue;
		}

		const std::string final_path = Rank_Snapshot_Final_Path(rank_snapshot_dir, rank);
		const bool final_exists = Path_Exists(final_path);
		if(Read_Rank_Snapshot_State(final_path, run_id, state) && state.done && state.rank_elapsed_wall_sec <= report.snapshot_time_label)
		{
			Accumulate_Snapshot_Report_State(report, state);
			report.rank_states.push_back(state);
			if(diagnostics != NULL)
			{
				diagnostics->ready_ranks++;
				std::ostringstream status;
				status << "final(done_at=" << std::fixed << std::setprecision(3) << state.rank_elapsed_wall_sec << "s)";
				diagnostics->per_rank_status[rank] = status.str();
			}
			continue;
		}

		all_ranks_ready = false;
		if(diagnostics != NULL)
		{
			std::ostringstream status;
			if(checkpoint_exists)
				status << "checkpoint_unreadable_or_stale";
			else
				status << "no_checkpoint";

			status << "; ";
			if(!final_exists)
				status << "no_final";
			else if(!Read_Rank_Snapshot_State(final_path, run_id, state))
				status << "final_unreadable_or_stale";
			else if(!state.done)
				status << "final_not_done";
			else if(state.rank_elapsed_wall_sec > report.snapshot_time_label)
				status << "final_after_target(done_at=" << std::fixed << std::setprecision(3) << state.rank_elapsed_wall_sec << "s)";
			else
				status << "final_unknown";

			diagnostics->per_rank_status[rank] = status.str();
		}
	}

	return all_ranks_ready;
}

bool Write_Snapshot_Report_Files(const std::string& snapshot_root, int snapshot_index, double interval_seconds, double mass_gev, double sigma_cm2, const Snapshot_Report_State& report)
{
	if(!Write_Text_File_Atomically(Snapshot_Text_File_Path(snapshot_root, snapshot_index, interval_seconds), snapshot_index, [&](std::ofstream& file)
	{
		file << "# Cumulative snapshot report\n";
		Write_Report_Header(file, mass_gev, sigma_cm2, report.total_trajectories, report.captured_particles, false, report.snapshot_time_label, report.snapshot_interval_seconds);

		double mean_wall_cap = (report.captured_particles > 0) ? report.total_wall_time_captured / static_cast<double>(report.captured_particles) : 0.0;
		double not_captured_particles = static_cast<double>(report.total_trajectories - report.captured_particles);
		double captured_particles = static_cast<double>(report.captured_particles);
		double mean_wall_nc = (not_captured_particles > 0.0) ? report.total_wall_time_not_captured / not_captured_particles : 0.0;
		double mean_steps_cap = (report.captured_particles > 0) ? static_cast<double>(report.total_rk45_steps_captured) / captured_particles : 0.0;
		double mean_steps_nc = (not_captured_particles > 0.0) ? static_cast<double>(report.total_rk45_steps_not_captured) / not_captured_particles : 0.0;

		file << "# total_wall_time_captured_sec = " << std::scientific << std::setprecision(6) << report.total_wall_time_captured << "\n";
		file << "# total_wall_time_not_captured_sec = " << std::scientific << std::setprecision(6) << report.total_wall_time_not_captured << "\n";
		file << "# mean_wall_time_captured_sec = " << std::scientific << std::setprecision(6) << mean_wall_cap << "\n";
		file << "# mean_wall_time_not_captured_sec = " << std::scientific << std::setprecision(6) << mean_wall_nc << "\n";
		file << "# total_rk45_steps_captured = " << report.total_rk45_steps_captured << "\n";
		file << "# total_rk45_steps_not_captured = " << report.total_rk45_steps_not_captured << "\n";
		file << "# mean_rk45_steps_captured = " << std::scientific << std::setprecision(6) << mean_steps_cap << "\n";
		file << "# mean_rk45_steps_not_captured = " << std::scientific << std::setprecision(6) << mean_steps_nc << "\n";
		file << "# ranks_ready = " << report.rank_states.size() << " / " << report.rank_states.size() << "\n";
		file << "#\n";
		file << "# Rank status:\n";
		for(const Rank_Snapshot_State& state : report.rank_states)
			file << "#   rank " << state.rank << ": " << Format_Rank_Status(state) << "\n";
		file << "#\n";
		file << "# bin_index  cap_dt[s]  cap_v2dt[km2/s]  cap_err_dt[s]  cap_err_v2dt[km2/s]  not_cap_dt[s]  not_cap_v2dt[km2/s]  not_cap_err_dt[s]  not_cap_err_v2dt[km2/s]\n";
		for(int bin = 0; bin < NUM_BINS; bin++)
		{
			double cap_err_dt = Snapshot_Bin_Error(report.captured_dt_hist[bin], report.captured_dt_sq_hist[bin], captured_particles);
			double cap_err_v2dt = Snapshot_Bin_Error(report.captured_v2dt_hist[bin], report.captured_v2dt_sq_hist[bin], captured_particles);
			double not_cap_err_dt = Snapshot_Bin_Error(report.not_captured_dt_hist[bin], report.not_captured_dt_sq_hist[bin], not_captured_particles);
			double not_cap_err_v2dt = Snapshot_Bin_Error(report.not_captured_v2dt_hist[bin], report.not_captured_v2dt_sq_hist[bin], not_captured_particles);

			file << bin << "\t" << std::scientific << std::setprecision(10)
			     << report.captured_dt_hist[bin] << "\t" << report.captured_v2dt_hist[bin]
			     << "\t" << cap_err_dt << "\t" << cap_err_v2dt
			     << "\t" << report.not_captured_dt_hist[bin] << "\t" << report.not_captured_v2dt_hist[bin]
			     << "\t" << not_cap_err_dt << "\t" << not_cap_err_v2dt << "\n";
		}
	}))
		return false;

	if(!Write_Text_File_Atomically(Snapshot_Report_File_Path(snapshot_root, snapshot_index, interval_seconds, "captured_bincount"), snapshot_index, [&](std::ofstream& file)
	{
		Write_Report_Header(file, mass_gev, sigma_cm2, report.total_trajectories, report.captured_particles, false, report.snapshot_time_label, report.snapshot_interval_seconds);
		file << "# bin_index  Sigma_dt[s]  Sigma_v2dt[km2/s]  err_dt[s]  err_v2dt[km2/s]\n";
		double N_cap = static_cast<double>(report.captured_particles);
		for(int bin = 0; bin < NUM_BINS; bin++)
		{
			double err_dt = Snapshot_Bin_Error(report.captured_dt_hist[bin], report.captured_dt_sq_hist[bin], N_cap);
			double err_v2dt = Snapshot_Bin_Error(report.captured_v2dt_hist[bin], report.captured_v2dt_sq_hist[bin], N_cap);

			file << bin << "\t" << std::scientific << std::setprecision(10)
			     << report.captured_dt_hist[bin] << "\t" << report.captured_v2dt_hist[bin]
			     << "\t" << err_dt << "\t" << err_v2dt << "\n";
		}
	}))
		return false;

	if(!Write_Text_File_Atomically(Snapshot_Report_File_Path(snapshot_root, snapshot_index, interval_seconds, "not_captured_bincount"), snapshot_index, [&](std::ofstream& file)
	{
		Write_Report_Header(file, mass_gev, sigma_cm2, report.total_trajectories, report.captured_particles, false, report.snapshot_time_label, report.snapshot_interval_seconds);
		file << "# bin_index  Sigma_dt[s]  Sigma_v2dt[km2/s]  err_dt[s]  err_v2dt[km2/s]\n";
		double N_nc = static_cast<double>(report.total_trajectories - report.captured_particles);
		for(int bin = 0; bin < NUM_BINS; bin++)
		{
			double err_dt = Snapshot_Bin_Error(report.not_captured_dt_hist[bin], report.not_captured_dt_sq_hist[bin], N_nc);
			double err_v2dt = Snapshot_Bin_Error(report.not_captured_v2dt_hist[bin], report.not_captured_v2dt_sq_hist[bin], N_nc);

			file << bin << "\t" << std::scientific << std::setprecision(10)
			     << report.not_captured_dt_hist[bin] << "\t" << report.not_captured_v2dt_hist[bin]
			     << "\t" << err_dt << "\t" << err_v2dt << "\n";
		}
	}))
		return false;

	return true;
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
		       << ", traj_physical_s=" << Format_Physical_Time_Scientific(state.current_trajectory_physical_sec);
	}
	else if(!state.done)
	{
		stream << ", next_traj=" << (state.local_total + 1);
	}

	stream << ")";
	return stream.str();
}

bool Try_Write_Merged_Snapshot(const std::string& snapshot_root, const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes, uint64_t run_id, double mass_gev, double sigma_cm2, int caller_rank)
{
	if(Path_Exists(Snapshot_Text_File_Path(snapshot_root, snapshot_index, interval_seconds)))
	{
		Cleanup_Snapshot_Checkpoints(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes);
		return true;
	}

	Snapshot_Report_State report;
	Snapshot_Load_Diagnostics diagnostics;
	if(!Load_Snapshot_Report_State(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes, run_id, report, &diagnostics))
	{
		Write_Snapshot_Status_File(snapshot_root, snapshot_index, interval_seconds, caller_rank, diagnostics, false);
		return false;
	}

	if(!Write_Snapshot_Report_Files(snapshot_root, snapshot_index, interval_seconds, mass_gev, sigma_cm2, report))
		return false;

	Write_Snapshot_Status_File(snapshot_root, snapshot_index, interval_seconds, caller_rank, diagnostics, true);
	Cleanup_Snapshot_Checkpoints(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes);
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
	const double snapshot_mass_gev = In_Units(DM.mass, GeV);
	const double snapshot_sigma_cm2 = In_Units(DM.Sigma_Proton(), cm * cm);
	std::string snapshot_root;
	std::string rank_snapshot_dir;
	uint64_t snapshot_run_id = 0;
	int next_snapshot_index = 1;
	if(snapshot_cfg.enabled)
	{
		snapshot_root = g_top_level_dir + "results_" + std::to_string(log10(In_Units(DM.mass, GeV))) + "_" + std::to_string(log10(In_Units(DM.Sigma_Proton(), cm * cm))) + "/snapshot/";
		rank_snapshot_dir = snapshot_root + "rank_snapshot/";
		Ensure_Directory_Exists(g_top_level_dir + "results_" + std::to_string(log10(In_Units(DM.mass, GeV))) + "_" + std::to_string(log10(In_Units(DM.Sigma_Proton(), cm * cm))) + "/");
		if(mpi_rank == 0)
		{
			Ensure_Directory_Exists(snapshot_root);
			Clear_Directory_Contents(snapshot_root);
			Ensure_Directory_Exists(snapshot_root);
			Ensure_Directory_Exists(rank_snapshot_dir);
		}
		MPI_Barrier(MPI_COMM_WORLD);
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
		state.captured_dt_sq_hist = captured_dt_sq_hist;
		state.captured_v2dt_sq_hist = captured_v2dt_sq_hist;
		state.not_captured_dt_hist = not_captured_dt_hist;
		state.not_captured_v2dt_hist = not_captured_v2dt_hist;
		state.not_captured_dt_sq_hist = not_captured_dt_sq_hist;
		state.not_captured_v2dt_sq_hist = not_captured_v2dt_sq_hist;
		state.total_wall_time_captured = total_wall_time_captured;
		state.total_wall_time_not_captured = total_wall_time_not_captured;
		state.total_rk45_steps_captured = static_cast<uint64_t>(total_rk45_steps_captured);
		state.total_rk45_steps_not_captured = static_cast<uint64_t>(total_rk45_steps_not_captured);
		std::copy(wall_time_hist_captured.begin(), wall_time_hist_captured.end(), state.wall_time_hist_captured.begin());
		std::copy(wall_time_hist_not_captured.begin(), wall_time_hist_not_captured.end(), state.wall_time_hist_not_captured.begin());
		state.wall_time_overflow_captured = static_cast<uint64_t>(wall_time_overflow_captured);
		state.wall_time_overflow_not_captured = static_cast<uint64_t>(wall_time_overflow_not_captured);
		std::copy(step_count_hist_captured.begin(), step_count_hist_captured.end(), state.step_count_hist_captured.begin());
		std::copy(step_count_hist_not_captured.begin(), step_count_hist_not_captured.end(), state.step_count_hist_not_captured.begin());
		state.step_count_overflow_captured = static_cast<uint64_t>(step_count_overflow_captured);
		state.step_count_overflow_not_captured = static_cast<uint64_t>(step_count_overflow_not_captured);
		state.evaporation_records = evaporation_records;
		return state;
	};

	auto try_write_ready_snapshots = [&](int max_snapshot_index)
	{
		bool all_snapshots_merged = true;
		for(int snapshot_index = 1; snapshot_index <= max_snapshot_index; snapshot_index++)
		{
			if(!Try_Write_Merged_Snapshot(snapshot_root, rank_snapshot_dir, snapshot_index, snapshot_interval, mpi_processes, snapshot_run_id, snapshot_mass_gev, snapshot_sigma_cm2, mpi_rank))
				all_snapshots_merged = false;
		}
		return all_snapshots_merged;
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
		if(try_write_ready_snapshots(final_snapshot_index))
			Cleanup_Final_Snapshot_States(rank_snapshot_dir, mpi_processes);
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
		Write_Report_Header(f, mass_gev, sigma_cm2, number_of_trajectories, number_of_captured_particles, early_stopped);
	};

	std::remove((output_dir + "/captured_bincount.txt").c_str());
	std::remove((output_dir + "/not_captured_bincount.txt").c_str());

	// 1. Merged bincount
	{
		std::ofstream f(output_dir + "/bincount.txt");
		write_header(f);
		f << "# bin_index  cap_dt[s]  cap_v2dt[km2/s]  cap_err_dt[s]  cap_err_v2dt[km2/s]  not_cap_dt[s]  not_cap_v2dt[km2/s]  not_cap_err_dt[s]  not_cap_err_v2dt[km2/s]\n";
		double N_cap = static_cast<double>(number_of_captured_particles);
		double N_nc = static_cast<double>(number_of_trajectories - number_of_captured_particles);
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

	// 2. Evaporation summary
	{
		std::ofstream f(output_dir + "/evaporation_summary.txt");
		write_header(f);
		f << "# trajectory_id  t_evap[s]  truncated(0/1)\n";
		for(const auto& rec : evaporation_records)
			f << rec.trajectory_id << "\t" << std::scientific << std::setprecision(10) << rec.t_evap << "\t" << (rec.truncated ? 1 : 0) << "\n";
		f.close();
	}

	// 3. Computation time & RK45 step count summary
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
		int wall_first_bin = 0;
		int wall_last_bin = -1;
		bool has_wall_bins = Find_Nonzero_Bin_Range(wall_time_hist_captured, wall_time_hist_not_captured, wall_first_bin, wall_last_bin);
		int wall_output_first_bin = has_wall_bins ? wall_first_bin - 1 : 0;
		f << "# [Wall-clock time histogram]\n";
		if(has_wall_bins)
		{
			double observed_log10_t_min = WALL_TIME_LOG_MIN + 0.1 * wall_output_first_bin;
			double observed_log10_t_max = WALL_TIME_LOG_MIN + 0.1 * (wall_last_bin + 1);
			f << "# occupied_log10_t_min = " << std::fixed << std::setprecision(1) << observed_log10_t_min << "\n";
			f << "# occupied_log10_t_max = " << std::fixed << std::setprecision(1) << observed_log10_t_max << "\n";
		}
		else
		{
			f << "# occupied_log10_t_min = nan\n";
			f << "# occupied_log10_t_max = nan\n";
		}
		f << "# bin_index  log10_t_lower  log10_t_upper  count_captured  count_not_captured\n";
		for(int b = wall_output_first_bin, compact_index = 0; has_wall_bins && b <= wall_last_bin; b++, compact_index++)
		{
			double lo = WALL_TIME_LOG_MIN + b * 0.1;
			double hi = lo + 0.1;
			f << compact_index << "\t" << std::fixed << std::setprecision(1) << lo << "\t" << hi
			  << "\t" << ((b >= 0) ? wall_time_hist_captured[b] : 0)
			  << "\t" << ((b >= 0) ? wall_time_hist_not_captured[b] : 0) << "\n";
		}
		f << "# overflow_captured = " << wall_time_overflow_captured << "\n";
		f << "# overflow_not_captured = " << wall_time_overflow_not_captured << "\n";
		f << "#\n";

		// RK45 step count histogram
		int step_first_bin = 0;
		int step_last_bin = -1;
		bool has_step_bins = Find_Nonzero_Bin_Range(step_count_hist_captured, step_count_hist_not_captured, step_first_bin, step_last_bin);
		f << "# [RK45 step count histogram]\n";
		if(has_step_bins)
		{
			double observed_log10_steps_min = STEP_COUNT_LOG_MIN + 0.1 * step_first_bin;
			double observed_log10_steps_max = STEP_COUNT_LOG_MIN + 0.1 * (step_last_bin + 1);
			f << "# occupied_log10_steps_min = " << std::fixed << std::setprecision(1) << observed_log10_steps_min << "\n";
			f << "# occupied_log10_steps_max = " << std::fixed << std::setprecision(1) << observed_log10_steps_max << "\n";
		}
		else
		{
			f << "# occupied_log10_steps_min = nan\n";
			f << "# occupied_log10_steps_max = nan\n";
		}
		f << "# bin_index  log10_steps_lower  log10_steps_upper  count_captured  count_not_captured\n";
		for(int b = step_first_bin, compact_index = 0; has_step_bins && b <= step_last_bin; b++, compact_index++)
		{
			double lo = STEP_COUNT_LOG_MIN + b * 0.1;
			double hi = lo + 0.1;
			f << compact_index << "\t" << std::fixed << std::setprecision(1) << lo << "\t" << hi
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

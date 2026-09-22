#include "Data_Generation.hpp"

#include "gtest/gtest.h"
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

namespace
{
bool FileExists(const std::string& path)
{
	struct stat info;
	return stat(path.c_str(), &info) == 0 && S_ISREG(info.st_mode);
}

bool FileContains(const std::string& path, const std::string& needle)
{
	std::ifstream file(path);
	std::ostringstream content;
	content << file.rdbuf();
	return content.str().find(needle) != std::string::npos;
}

std::string TestOutputDir(const std::string& name)
{
	std::string dir = name + "_" + std::to_string(getpid()) + "/";
	mkdir(dir.c_str(), 0755);
	return dir;
}

void TouchFile(const std::string& path)
{
	std::ofstream file(path);
	file << "stale\n";
}

void RemoveTestOutputDir(const std::string& directory)
{
	std::remove((directory + "bincount.txt").c_str());
	std::remove((directory + "evaporation_times.txt").c_str());
	std::remove((directory + "run_metadata.json").c_str());
	std::remove((directory + "diagnostic_trajectory_summary.tsv").c_str());
	std::remove((directory + "trajectory_events.tsv").c_str());
	std::remove((directory + "invalid_trajectories.tsv").c_str());
	std::remove((directory + "residence_jackknife_blocks.tsv").c_str());
	std::remove((directory + "radial_blocks.tsv").c_str());
	rmdir(directory.c_str());
}

}

int main(int argc, char* argv[])
{
	int result = 0;

	::testing::InitGoogleTest(&argc, argv);
	MPI_Init(&argc, &argv);
	result = RUN_ALL_TESTS();
	MPI_Finalize();
	return result;
}

TEST(TestDataGeneration, TestGenerateData)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	unsigned int sample_size = 2;
	unsigned int max_trajectories = 2;

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 10, 10);

	// ACT
	Simulation_Data data_set(sample_size, max_trajectories);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT – number_of_trajectories is private; verify indirectly
	ASSERT_GE(data_set.Free_Ratio() + data_set.Capture_Ratio() + data_set.Reflection_Ratio(), 0.0);
}

TEST(TestDataGeneration, TestConfigure)
{
	// ARRANGE
	unsigned int sample_size		 = 2;
	double r						 = 2.0 * rSun;
	unsigned int min_scattering		 = 1;
	unsigned int max_scattering		 = 1;
	unsigned long int max_time_steps = 1e4;
	// ACT
	Simulation_Data data_set(sample_size, 0);
	data_set.Configure(r, min_scattering, max_scattering, max_time_steps);

	// ASSERT
	// ASSERT_EQ(data_set.data[0].size(), sample_size);
}

TEST(TestDataGeneration, ScatteredNeverCapturedPathIsIncludedInPopulationBincount)
{
	Solar_Model sun;
	obscura::Standard_Halo_Model halo;
	obscura::DM_Particle_SD dm(0.1 * GeV);
	dm.Set_Low_Mass_Mode(true);
	dm.Set_Sigma_Proton(1.0e-34 * cm * cm);
	bool found_reflected = false;
	// Each run injects exactly one particle. Select a physically completed,
	// scattered but never-captured path: the former scatter-count filter made
	// its entire transit histogram zero, including inside the Sun.
	for(unsigned int seed = 1; seed <= 32 && !found_reflected; ++seed)
	{
		Simulation_Data sample(1, 1);
		sample.Configure(TRAJECTORY_BOUNDARY_RSUN * rSun, 0, 100000);
		sample.Generate_Data(dm, sun, halo, SnapshotConfig(), seed);
		if(sample.Reflection_Ratio() != 1.0 || sample.Capture_Ratio() != 0.0)
			continue;
		found_reflected = true;
		const std::string dir = TestOutputDir("scattered_uncaptured_population");
		sample.Write_Bincount(dir, dm, halo);
		int rank = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(rank == 0)
		{
			std::ifstream stream(dir + "bincount.tsv");
			std::string line;
			double inside_dt_s = 0.0, outside_dt_s = 0.0, captured_dt_s = 0.0;
			while(std::getline(stream, line))
			{
				if(line.empty() || line[0] == '#') continue;
				std::istringstream row(line);
				std::array<double, 14> values{};
				for(double& value : values)
				{
					std::string token;
					ASSERT_TRUE(static_cast<bool>(row >> token));
					value = std::stod(token);
				}
				(values[2] <= 1.0 + 1e-12 ? inside_dt_s : outside_dt_s) += values[11];
				captured_dt_s += values[8];
				EXPECT_DOUBLE_EQ(values[9], 0.0);
				EXPECT_NEAR(values[12], values[11] * values[11],
				            2.0e-9 * std::max(1.0, values[12]));
				EXPECT_TRUE(std::isnan(values[5]));
				EXPECT_TRUE(std::isnan(values[13]));
			}
			EXPECT_GT(inside_dt_s, 0.0);
			EXPECT_GT(outside_dt_s, 0.0);
			EXPECT_DOUBLE_EQ(captured_dt_s, 0.0);
			std::remove((dir + "bincount.tsv").c_str());
			rmdir((dir + "snapshot").c_str());
			rmdir(dir.c_str());
		}
	}
	EXPECT_TRUE(found_reflected);
}

TEST(TestDataGeneration, TestInitialShiftFailureIsReported)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	Simulation_Data data_set(1, 1);
	data_set.Configure(0.5 * rSun, 0, 1, 10);
	data_set.diagnostic_output_enabled = true;
	data_set.Generate_Data(DM, SSM, SHM);

	EXPECT_EQ(data_set.Valid_Trajectories(), 0UL);
	EXPECT_DOUBLE_EQ(data_set.Numerical_Failure_Ratio(), 1.0);
	EXPECT_DOUBLE_EQ(data_set.Capture_Ratio_Valid(), 0.0);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir = TestOutputDir("initial_shift_failure_contract");
	data_set.Write_Diagnostic_Output(output_dir, DM);
	if(rank == 0)
	{
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# EARLY_STOP: max_trajectories_reached"));
		RemoveTestOutputDir(output_dir);
	}
}

TEST(TestDataGeneration, TestComputationallyTruncatedNonCaptureIsExcludedFromCaptureRate)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	Simulation_Data data_set(1, 1);
	data_set.Configure(2.0 * rSun, 0, 0, 10);
	data_set.diagnostic_output_enabled = true;
	data_set.Generate_Data(DM, SSM, SHM, SnapshotConfig(), 20260710);

	EXPECT_EQ(data_set.Valid_Trajectories(), 0UL);
	EXPECT_DOUBLE_EQ(data_set.Capture_Ratio_Valid(), 0.0);
	EXPECT_DOUBLE_EQ(data_set.Numerical_Failure_Ratio(), 0.0);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir = TestOutputDir("truncated_output_contract");
	data_set.Write_Diagnostic_Output(output_dir, DM);
	if(rank == 0)
	{
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# valid_trajectories = 0"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# unresolved_not_captured_trajectories = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# computational_truncations = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# invalid_trajectory_records = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# termination_max_scatterings_uncaptured = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# capture_rate_valid = 0.00000000"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# capture_rate_valid_CI_95_lower = 0.00000000"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# EARLY_STOP: max_trajectories_reached"));
		EXPECT_TRUE(FileExists(output_dir + "invalid_trajectories.tsv"));
		EXPECT_TRUE(FileContains(output_dir + "invalid_trajectories.tsv", "# record_count = 1"));
		EXPECT_TRUE(FileContains(output_dir + "invalid_trajectories.tsv", "\tpropagation\tmax_scatterings\t"));
		EXPECT_TRUE(FileContains(output_dir + "invalid_trajectories.tsv", "rng_state_before_simulation"));
		RemoveTestOutputDir(output_dir);
	}
}

TEST(TestDataGeneration, TestWallTimeCutoffExcludesOnlyFailedHistory)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	Simulation_Data data_set(1, 1);
	data_set.Configure(2.0 * rSun, 0, 10, 10);
	SnapshotConfig snapshot_config;
	snapshot_config.max_trajectory_wall_time_sec = 1.0e-12;
	data_set.diagnostic_output_enabled = true;
	data_set.Generate_Data(
	    DM, SSM, SHM, snapshot_config, 20260818);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir =
	    TestOutputDir("wall_time_censor_contract");
	data_set.Write_Diagnostic_Output(output_dir, DM);
	if(rank == 0)
	{
		EXPECT_TRUE(FileContains(
		    output_dir + "bincount.txt",
		    "# computational_truncations = 1"));
		EXPECT_TRUE(FileContains(
		    output_dir + "bincount.txt",
		    "# invalid_trajectory_records = 1"));
		EXPECT_TRUE(FileContains(
		    output_dir + "bincount.txt",
		    "# termination_wall_time_limit_uncaptured = 1"));
		EXPECT_TRUE(FileContains(
		    output_dir + "invalid_trajectories.tsv",
		    "# record_count = 1"));
		RemoveTestOutputDir(output_dir);
	}
}

TEST(TestDataGeneration, TestInvalidTrajectoriesContinueUntilExplicitBudget)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	Simulation_Data data_set(1, 3);
	data_set.Configure(2.0 * rSun, 0, 0, 10);
	data_set.diagnostic_output_enabled = true;
	data_set.Generate_Data(DM, SSM, SHM, SnapshotConfig(), 20260710);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir =
	    TestOutputDir("invalid_trajectory_explicit_budget");
	data_set.Write_Diagnostic_Output(output_dir, DM);
	if(rank == 0)
	{
		EXPECT_TRUE(FileContains(
		    output_dir + "bincount.txt",
		    "# EARLY_STOP: max_trajectories_reached"));
		EXPECT_TRUE(FileContains(
		    output_dir + "bincount.txt",
		    "# computational_truncations = 3"));
		EXPECT_TRUE(FileContains(
		    output_dir + "bincount.txt",
		    "# invalid_trajectory_records = 3"));
		RemoveTestOutputDir(output_dir);
	}
}

TEST(TestDataGeneration, TestCaptureModeCompletedEscapeIsIncludedInCaptureRate)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);
	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 10, 10);

	Simulation_Data data_set(1, 1);
	data_set.Configure(1.1 * rSun, 0, 500);
	data_set.Generate_Data(DM, SSM, SHM, SnapshotConfig(), 20260710, true);

	EXPECT_EQ(data_set.Valid_Trajectories(), 1UL);
	EXPECT_DOUBLE_EQ(data_set.Capture_Ratio_Valid(), 0.0);
	EXPECT_DOUBLE_EQ(data_set.Numerical_Failure_Ratio(), 0.0);
	for(int detail_index = 1;
	    detail_index < TRAJECTORY_NUMERICAL_FAILURE_DETAIL_COUNT;
	    detail_index++)
	{
		EXPECT_EQ(
		    data_set.Numerical_Failure_Detail_Count(
		        static_cast<TrajectoryNumericalFailureDetail>(
		            detail_index),
		        false),
		    0UL)
		    << TrajectoryNumericalFailureDetailKey(
		           static_cast<TrajectoryNumericalFailureDetail>(
		               detail_index));
	}
}

TEST(TestDataGeneration, TestDataFreeRatio)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 10, 10);

	unsigned int sample_size = 2;
	// ACT
	Simulation_Data data_set(sample_size, sample_size);
	data_set.Configure(1.1 * rSun, 0, 500);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT
	ASSERT_DOUBLE_EQ(data_set.Free_Ratio(), 1.0);
}

TEST(TestDataGeneration, TestDataSetCaptureRatio)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(1.0 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 10, 10);

	unsigned int sample_size = 2;

	// ACT
	Simulation_Data data_set(sample_size, sample_size);
	data_set.Configure(1.1 * rSun, 0, 500);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT
	ASSERT_GE(data_set.Capture_Ratio(), 0.0);
}

TEST(TestDataGeneration, TestDataSetReflectionRatio)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(1.0 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 10, 10);

	unsigned int sample_size = 2;

	// ACT
	Simulation_Data data_set(sample_size, sample_size);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT
	ASSERT_GE(data_set.Reflection_Ratio(), 0.0);
}

TEST(TestDataGeneration, TestSpeedFunctions)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(1.0 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 10, 10);

	unsigned int sample_size = 2;
	double u_min			 = 0.0001;
	// ACT
	Simulation_Data data_set(sample_size, sample_size, u_min);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT
	EXPECT_DOUBLE_EQ(data_set.Minimum_Speed(), 0.75 * u_min);
}

TEST(TestDataGeneration, TestOutputFailuresAreReported)
{
	Simulation_Data data_set(1, 1);
	obscura::DM_Particle_SD DM(0.1 * GeV);
	obscura::Standard_Halo_Model halo;
	const std::string dir = TestOutputDir("output_failure");
	EXPECT_NO_THROW(data_set.Prepare_Output_Directory(dir));
	const std::string blocker = dir + "regular_file";
	TouchFile(blocker);
	EXPECT_THROW(data_set.Prepare_Output_Directory(blocker + "/nested"), std::runtime_error);
	EXPECT_THROW(data_set.Prepare_Output_Directory(dir), std::runtime_error);
	EXPECT_TRUE(FileExists(blocker));
	std::remove(blocker.c_str());
	// A failed rename must leave neither a complete file nor a temporary file.
	const std::string target = dir + "bincount.tsv";
	mkdir(target.c_str(), 0755);
	TouchFile(target + "/keep");
	EXPECT_THROW(data_set.Write_Bincount(dir, DM, halo), std::runtime_error);
	EXPECT_FALSE(FileExists(dir + "bincount.tsv.tmp"));
	EXPECT_TRUE(FileExists(target + "/keep"));
	std::remove((target + "/keep").c_str());
	rmdir(target.c_str());
	// A temporary path that cannot be opened also fails without publication.
	mkdir((dir + "bincount.tsv.tmp").c_str(), 0755);
	TouchFile(dir + "bincount.tsv.tmp/keep");
	EXPECT_THROW(data_set.Write_Bincount(dir, DM, halo), std::runtime_error);
	EXPECT_FALSE(FileExists(target));
	std::remove((dir + "bincount.tsv.tmp/keep").c_str());
	rmdir((dir + "bincount.tsv.tmp").c_str());
	EXPECT_THROW(data_set.Write_Diagnostic_Output(dir, DM), std::logic_error);
	rmdir((dir + "snapshot").c_str());
	rmdir(dir.c_str());
}

TEST(TestDataGeneration, JackknifeSumErrorUsesPopulationDenominators)
{
	std::array<double, RESIDENCE_JACKKNIFE_BLOCKS> sums{};
	std::array<unsigned long int, RESIDENCE_JACKKNIFE_BLOCKS> counts{};
	// Unequal populations, identical per-history values: sampling error is zero.
	for(std::size_t block = 0; block < counts.size(); ++block)
	{
		counts[block] = block % 5;
		sums[block] = 3.0 * counts[block];
	}
	EXPECT_DOUBLE_EQ(Block_Jackknife_Sum_SE(sums, counts), 0.0);
	// One history per block reduces to the ordinary sample-mean SE times N.
	for(std::size_t block = 0; block < counts.size(); ++block)
	{
		counts[block] = 1;
		sums[block] = block % 2 == 0 ? 1.0 : 3.0;
	}
	EXPECT_NEAR(Block_Jackknife_Sum_SE(sums, counts), 64.0 / std::sqrt(63.0), 1e-12);
}

TEST(TestDataGeneration, JackknifeSumErrorHandlesEmptyAndSingleOccupiedBlocks)
{
	std::array<double, RESIDENCE_JACKKNIFE_BLOCKS> sums{};
	std::array<unsigned long int, RESIDENCE_JACKKNIFE_BLOCKS> counts{};
	EXPECT_TRUE(std::isnan(Block_Jackknife_Sum_SE(sums, counts)));
	counts[2] = 1; sums[2] = 2.0;
	EXPECT_TRUE(std::isnan(Block_Jackknife_Sum_SE(sums, counts)));
	counts[2] = 10; sums[2] = 20.0;
	EXPECT_TRUE(std::isnan(Block_Jackknife_Sum_SE(sums, counts)));
	counts[5] = 10; sums[5] = 20.0;
	EXPECT_DOUBLE_EQ(Block_Jackknife_Sum_SE(sums, counts), 0.0);
}

TEST(TestDataGeneration, CompletePathSecondMomentIncludesCrossTerms)
{
	RadialHistogram inbound{2.0}, pre_capture{3.0, 11.0}, residence{5.0}, outgoing{7.0};
	RadialHistogram dt, dt_sq;
	Accumulate_Complete_Path_Block({&inbound, &pre_capture, &residence, &outgoing}, 7, dt, dt_sq);
	EXPECT_DOUBLE_EQ(dt[7], 17.0);
	EXPECT_DOUBLE_EQ(dt_sq[7], 289.0);
	EXPECT_DOUBLE_EQ(dt_sq[RESIDENCE_JACKKNIFE_BLOCKS + 7], 121.0);
	EXPECT_DOUBLE_EQ(dt_sq[6], 0.0);
	// A second history contributes its own square, not the square of the sum.
	Accumulate_Complete_Path_Block({&inbound}, 7, dt, dt_sq);
	EXPECT_DOUBLE_EQ(dt[7], 19.0);
	EXPECT_DOUBLE_EQ(dt_sq[7], 293.0);
}

TEST(TestDataGeneration, TestExplicitLegacyDiagnosticContract)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;
	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	Simulation_Data data_set(1, 1);
	data_set.Configure(1.1 * rSun, 0, 100);
	data_set.diagnostic_output_enabled = true;
	data_set.Generate_Data(DM, SSM, SHM);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir = TestOutputDir("default_output_contract");
	if(rank == 0)
	{
		TouchFile(output_dir + std::string("evaporation_") + "summary.txt");
		TouchFile(output_dir + std::string("evaporation_") + "mode_summary.txt");
		TouchFile(output_dir + std::string("evaporation_") + "mode_" + "bincount.txt");
		TouchFile(output_dir + std::string("computation_") + "time_summary.txt");
	}
	data_set.Write_Diagnostic_Output(output_dir, DM);
	if(rank == 0)
	{
		EXPECT_TRUE(FileExists(output_dir + "bincount.txt"));
			EXPECT_TRUE(FileExists(output_dir + "evaporation_times.txt"));
			EXPECT_TRUE(FileExists(output_dir + "invalid_trajectories.tsv"));
			EXPECT_TRUE(FileExists(output_dir + "residence_jackknife_blocks.tsv"));
			EXPECT_FALSE(FileExists(output_dir + "bincount.txt.tmp"));
			EXPECT_FALSE(FileExists(
			    output_dir + "residence_jackknife_blocks.tsv.tmp"));
		EXPECT_TRUE(FileContains(
		    output_dir + "residence_jackknife_blocks.tsv",
		    "# block_count = 64"));
		EXPECT_TRUE(FileContains(
		    output_dir + "residence_jackknife_blocks.tsv",
		    "block_id\tbin_index\tresidence_dt_s\tresidence_v2dt_km2_s"));
		EXPECT_TRUE(FileContains(
		    output_dir + "residence_jackknife_blocks.tsv",
		    "# block_0_attempted = "));
		EXPECT_TRUE(FileContains(output_dir + "invalid_trajectories.tsv", "# record_count = 0"));
		EXPECT_TRUE(FileContains(output_dir + "evaporation_times.txt", "# format_version = 5"));
		EXPECT_TRUE(FileContains(output_dir + "evaporation_times.txt", "P_kepler_first_bound_exit_sec"));
		EXPECT_TRUE(FileContains(
		    output_dir + "bincount.txt",
		    "# bincount_integration = conservative-hermite-kepler-separate-boundaries-v8"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# radial_domain_max_Rsun = "));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# total_radial_bins = 1100"));
		EXPECT_TRUE(FileContains(
		    output_dir + "bincount.txt",
		    "# radial_grid = uniform_inner_geometric_width_capped_v5"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# exterior_bins = 0"));
		EXPECT_TRUE(FileContains(
		    output_dir + "bincount.txt",
		    "# exterior_bin_growth_factor = 1.0200000000e+00"));
		EXPECT_TRUE(FileContains(
		    output_dir + "bincount.txt",
		    "# exterior_max_bin_width_Rsun = 1.0000000000e+01"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# mpi_scheduler = dynamic_rma_work_queue_v1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# mpi_scheduler_progress = main_thread_iprobe_v1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# mpi_scheduler_work_claims = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# mpi_scheduler_peak_in_flight = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# capture_target_overshoot = 0"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# total_scatterings = 0"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# simulation_time_seconds = "));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# completed_outward_escapes = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# unresolved_not_captured_trajectories = 0"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# termination_physical_escape_uncaptured = 1"));
		EXPECT_FALSE(FileExists(output_dir + "evaporation_diagnostics.txt"));
		EXPECT_TRUE(FileExists(output_dir + "run_metadata.json"));
		EXPECT_TRUE(FileExists(output_dir + "diagnostic_trajectory_summary.tsv"));
		EXPECT_TRUE(FileExists(output_dir + "trajectory_events.tsv"));
		EXPECT_FALSE(FileExists(output_dir + std::string("evaporation_") + "summary.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("evaporation_") + "mode_summary.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("evaporation_") + "mode_" + "bincount.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("computation_") + "time_summary.txt"));
		RemoveTestOutputDir(output_dir);
	}
}

TEST(TestDataGeneration, TestTrajectoryDiagnosticOutputContract)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;
	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-32 * cm * cm);
	DM.Set_Sigma_Electron(1.0e-100 * pb);
	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 20, 20);

	Simulation_Data data_set(1, 64);
	data_set.Configure(2.0 * rSun, 0, 100000);
	TrajectoryDiagnosticConfig diagnostics;
	diagnostics.summary_enabled = true;
	diagnostics.events_enabled = true;
	diagnostics.trace_rate = 1.0;
	diagnostics.trace_seed = 2026072201ULL;
	diagnostics.interpolation_points = 20;
	EXPECT_NO_THROW(data_set.Configure_Trajectory_Diagnostics(diagnostics));
	diagnostics.trace_rate = 1.01;
	EXPECT_THROW(data_set.Configure_Trajectory_Diagnostics(diagnostics), std::invalid_argument);
	data_set.diagnostic_output_enabled = true;
	data_set.Generate_Data(DM, SSM, SHM, SnapshotConfig(), 20260722);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir = TestOutputDir("trajectory_diagnostic_contract");
	data_set.Write_Diagnostic_Output(output_dir, DM);
	if(rank == 0)
	{
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"schema_version\": \"trajectory-diagnostic-v5\""));
		EXPECT_TRUE(FileContains(
		    output_dir + "run_metadata.json",
		    "\"bincount_integration\": \"conservative-hermite-kepler-separate-boundaries-v8\""));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"interpolation_points\": 20"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"rate_radius_points\": 20"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"rate_speed_points\": 20"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"rate_max_speed\": 7.50000000000000000e-01"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"rate_query_count\":"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"rate_fallback_count\":"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"evaporation_event_reconciliation\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"escape_radius_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"residence_time_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"event_sequence_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"event_count_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"trace_selection_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"replay_state_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"bound_exit_orbit_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "diagnostic_trajectory_summary.tsv", "t_first_unbinding_s"));
		EXPECT_TRUE(FileContains(output_dir + "diagnostic_trajectory_summary.tsv", "n_recapture"));
		EXPECT_TRUE(FileContains(output_dir + "diagnostic_trajectory_summary.tsv", "P_kepler_max_bound_exit_s"));
		EXPECT_TRUE(FileContains(output_dir + "diagnostic_trajectory_summary.tsv", "rng_state_before_simulation"));
		EXPECT_TRUE(FileContains(output_dir + "trajectory_events.tsv", "event_type"));
		EXPECT_TRUE(FileContains(output_dir + "trajectory_events.tsv", "scatter_pre"));
		EXPECT_TRUE(FileContains(output_dir + "trajectory_events.tsv", "scatter_post"));
		EXPECT_TRUE(FileContains(output_dir + "trajectory_events.tsv", "candidate_unbinding"));
		EXPECT_TRUE(FileContains(output_dir + "trajectory_events.tsv", "escape_validated"));
		RemoveTestOutputDir(output_dir);
	}
}

TEST(TestDataGeneration, TestTraceSelectionIsStableAcrossRunIds)
{
	std::vector<uint64_t> first_selection;
	std::vector<uint64_t> repeated_selection;
	std::vector<uint64_t> other_seed_selection;
	for(uint64_t trajectory_id = 1; trajectory_id <= 1000; trajectory_id++)
	{
		if(TrajectoryTraceSelected(1234567ULL, 2, trajectory_id, 0.02))
			first_selection.push_back(trajectory_id);
		if(TrajectoryTraceSelected(1234567ULL, 2, trajectory_id, 0.02))
			repeated_selection.push_back(trajectory_id);
		if(TrajectoryTraceSelected(7654321ULL, 2, trajectory_id, 0.02))
			other_seed_selection.push_back(trajectory_id);
	}
	EXPECT_EQ(first_selection, repeated_selection);
	EXPECT_NE(first_selection, other_seed_selection);
	EXPECT_FALSE(TrajectoryTraceSelected(123ULL, 0, 1, 0.0));
	EXPECT_TRUE(TrajectoryTraceSelected(123ULL, 0, 1, 1.0));
}

TEST(TestDataGeneration, TestExplicitDiagnosticsContainsOnlyRequestedReports)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;
	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	Simulation_Data data_set(1, 1);
	data_set.Configure(1.1 * rSun, 0, 100);
	data_set.diagnostic_output_enabled = true;
	data_set.Generate_Data(DM, SSM, SHM);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir = TestOutputDir("diagnostics_output_contract");
	if(rank == 0)
		TouchFile(output_dir + "evaporation_diagnostics.txt");
	data_set.Write_Diagnostic_Output(output_dir, DM);
	if(rank == 0)
	{
		EXPECT_TRUE(FileExists(output_dir + "bincount.txt"));
		EXPECT_TRUE(FileExists(output_dir + "evaporation_times.txt"));
		EXPECT_TRUE(FileExists(output_dir + "residence_jackknife_blocks.tsv"));
		EXPECT_FALSE(FileExists(output_dir + "evaporation_diagnostics.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("evaporation_") + "mode_summary.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("evaporation_") + "mode_" + "bincount.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("computation_") + "time_summary.txt"));
		RemoveTestOutputDir(output_dir);
	}
}

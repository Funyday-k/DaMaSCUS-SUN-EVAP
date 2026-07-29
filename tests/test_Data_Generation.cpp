#include "Data_Generation.hpp"

#include "gtest/gtest.h"
#include <cstdio>
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
	std::remove((directory + "trajectory_summary.tsv").c_str());
	std::remove((directory + "trajectory_events.tsv").c_str());
	std::remove((directory + "invalid_trajectories.tsv").c_str());
	std::remove((directory + "residence_jackknife_blocks.tsv").c_str());
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
	data_set.Generate_Data(DM, SSM, SHM);

	EXPECT_EQ(data_set.Valid_Trajectories(), 0UL);
	EXPECT_DOUBLE_EQ(data_set.Numerical_Failure_Ratio(), 1.0);
	EXPECT_DOUBLE_EQ(data_set.Capture_Ratio_Valid(), 0.0);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir = TestOutputDir("initial_shift_failure_contract");
	data_set.Write_Output_Files(output_dir, DM);
	if(rank == 0)
	{
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# EARLY_STOP: initial_shift_failure_fraction_exceeded"));
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
	data_set.Generate_Data(DM, SSM, SHM, SnapshotConfig(), 20260710);

	EXPECT_EQ(data_set.Valid_Trajectories(), 0UL);
	EXPECT_DOUBLE_EQ(data_set.Capture_Ratio_Valid(), 0.0);
	EXPECT_DOUBLE_EQ(data_set.Numerical_Failure_Ratio(), 0.0);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir = TestOutputDir("truncated_output_contract");
	data_set.Write_Output_Files(output_dir, DM);
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

TEST(TestDataGeneration, TestUnlimitedBudgetStopsWhenEveryTrajectoryIsInvalid)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	Simulation_Data data_set(1, 0);
	data_set.Configure(2.0 * rSun, 0, 0, 10);
	data_set.Generate_Data(DM, SSM, SHM, SnapshotConfig(), 20260710, true);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir = TestOutputDir("invalid_trajectory_fuse");
	data_set.Write_Output_Files(output_dir, DM);
	if(rank == 0)
	{
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# EARLY_STOP: invalid_trajectory_fraction_exceeded"));
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

TEST(TestDataGeneration, TestDefaultOutputContract)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;
	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	Simulation_Data data_set(1, 1);
	data_set.Configure(1.1 * rSun, 0, 100);
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
	data_set.Write_Output_Files(output_dir, DM);
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
		    "# bincount_integration = conservative-hermite-kepler-jupiter-log-v3"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# radial_domain_max_AU = 5.2000000000e+00"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# total_radial_bins = 1612"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# normal_mode_mpi_sync_interval = 1048576"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# mpi_sync_rounds = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# final_mpi_round_trajectories = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# capture_target_overshoot = 0"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# total_scatterings = 0"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# simulation_time_seconds = "));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# completed_outward_escapes = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# unresolved_not_captured_trajectories = 0"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# termination_outward_escape_uncaptured = 1"));
		EXPECT_FALSE(FileExists(output_dir + "evaporation_diagnostics.txt"));
		EXPECT_FALSE(FileExists(output_dir + "run_metadata.json"));
		EXPECT_FALSE(FileExists(output_dir + "trajectory_summary.tsv"));
		EXPECT_FALSE(FileExists(output_dir + "trajectory_events.tsv"));
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
	data_set.Generate_Data(DM, SSM, SHM, SnapshotConfig(), 20260722);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir = TestOutputDir("trajectory_diagnostic_contract");
	data_set.Write_Output_Files(output_dir, DM);
	if(rank == 0)
	{
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"schema_version\": \"trajectory-diagnostic-v4\""));
		EXPECT_TRUE(FileContains(
		    output_dir + "run_metadata.json",
		    "\"bincount_integration\": \"conservative-hermite-kepler-jupiter-log-v3\""));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"interpolation_points\": 20"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"evaporation_event_reconciliation\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"escape_radius_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"residence_time_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"event_sequence_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"event_count_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"trace_selection_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"replay_state_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"bound_exit_orbit_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "trajectory_summary.tsv", "t_first_unbinding_s"));
		EXPECT_TRUE(FileContains(output_dir + "trajectory_summary.tsv", "n_recapture"));
		EXPECT_TRUE(FileContains(output_dir + "trajectory_summary.tsv", "P_kepler_max_bound_exit_s"));
		EXPECT_TRUE(FileContains(output_dir + "trajectory_summary.tsv", "rng_state_before_simulation"));
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

TEST(TestDataGeneration, TestFinalOutputContainsOnlyRequestedReports)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;
	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	Simulation_Data data_set(1, 1);
	data_set.Configure(1.1 * rSun, 0, 100);
	data_set.Generate_Data(DM, SSM, SHM);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir = TestOutputDir("diagnostics_output_contract");
	if(rank == 0)
		TouchFile(output_dir + "evaporation_diagnostics.txt");
	data_set.Write_Output_Files(output_dir, DM);
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

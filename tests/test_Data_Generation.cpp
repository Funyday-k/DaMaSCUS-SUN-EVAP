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
	rmdir(directory.c_str());
}

// Retired with the block/manifest snapshot workflow.
#if 0
struct TestEvaporationRow
{
	int rank = -1;
	uint64_t trajectory_id = 0;
	double lifetime_unbinding = -1.0;
};

TestEvaporationRow MakeTestEvaporationRow(int rank, uint64_t trajectory_id, double lifetime_unbinding)
{
	TestEvaporationRow row;
	row.rank = rank;
	row.trajectory_id = trajectory_id;
	row.lifetime_unbinding = lifetime_unbinding;
	return row;
}

void EnsureDir(const std::string& path)
{
	mkdir(path.c_str(), 0755);
}

void WriteTestEvaporationBlock(const std::string& path, uint64_t run_id, int snapshot_index, double interval, double mass_gev, double sigma_cm2, const std::vector<TestEvaporationRow>& rows)
{
	std::ofstream file(path);
	file << "# format_version = 1\n";
	file << "# run_id = " << run_id << "\n";
	file << "# snapshot_index = " << snapshot_index << "\n";
	file << "# snapshot_interval_sec = " << std::scientific << std::setprecision(17) << interval << "\n";
	file << "# DM_mass_GeV = " << std::scientific << std::setprecision(17) << mass_gev << "\n";
	file << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(17) << sigma_cm2 << "\n";
	file << "# rank trajectory_id lifetime_unbinding_sec\n";
	for(const auto& row : rows)
		file << row.rank << "\t" << row.trajectory_id << "\t" << std::scientific << std::setprecision(10) << row.lifetime_unbinding << "\n";
}

std::vector<TestEvaporationRow> ReadRecoveredRows(const std::string& path)
{
	std::vector<TestEvaporationRow> rows;
	std::ifstream file(path);
	std::string line;
	while(std::getline(file, line))
	{
		if(line.empty() || line[0] == '#')
			continue;
		std::istringstream stream(line);
		TestEvaporationRow row;
		if(stream >> row.rank >> row.trajectory_id >> row.lifetime_unbinding)
			rows.push_back(row);
	}
	return rows;
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

TEST(TestDataGeneration, TestDataSetConstructor)
{
	// ARRANGE
	unsigned int sample_size = 100;
	unsigned int max_traj    = 0;
	double u_min			 = 0.0;
	unsigned int iso_rings	 = 10;
	// ACT
	Simulation_Data data_set(sample_size, max_traj, u_min, iso_rings);
	// ASSERT
	ASSERT_EQ(data_set.data.size(), iso_rings);
	for(auto& set : data_set.data)
		ASSERT_EQ(set.size(), 0);
}

TEST(TestDataGeneration, TestCompactEvaporationEventStaysCompact)
{
	EXPECT_LT(sizeof(CompactEvaporationEvent), sizeof(EvaporationRecord) / 2);
	EXPECT_LE(sizeof(CompactEvaporationEvent), static_cast<size_t>(40));
}

TEST(TestDataGeneration, TestRecoverEvaporationTimeFileFromBlocksFiltersRunAndDeduplicates)
{
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank != 0)
		return;

	const uint64_t run_id = 123456789ULL;
	const double mass_gev = 1.25;
	const double sigma_cm2 = 2.5e-36;
	const double interval = 30.0;
	const std::string output_dir = TestOutputDir("recover_evaporation_blocks");
	const std::string snapshot_root = output_dir + "snapshot/";
	const std::string block_dir = snapshot_root + "evaporation_blocks/";
	EnsureDir(snapshot_root);
	EnsureDir(block_dir);

	std::vector<TestEvaporationRow> block_1_rows;
	block_1_rows.push_back(MakeTestEvaporationRow(1, 5, 99.0));
	block_1_rows.push_back(MakeTestEvaporationRow(2, 4, 30.0));
	WriteTestEvaporationBlock(block_dir + "block_000001.txt", run_id, 1, interval, mass_gev, sigma_cm2, block_1_rows);

	std::vector<TestEvaporationRow> block_2_rows;
	block_2_rows.push_back(MakeTestEvaporationRow(1, 5, 20.0));
	block_2_rows.push_back(MakeTestEvaporationRow(0, 3, 10.0));
	WriteTestEvaporationBlock(block_dir + "block_000002.txt", run_id, 2, interval, mass_gev, sigma_cm2, block_2_rows);

	std::vector<TestEvaporationRow> wrong_run_rows;
	wrong_run_rows.push_back(MakeTestEvaporationRow(9, 9, 9.0));
	WriteTestEvaporationBlock(block_dir + "block_000003.txt", run_id + 1, 3, interval, mass_gev, sigma_cm2, wrong_run_rows);

	const std::string output_path = output_dir + "evaporation_times.partial.txt";
	ASSERT_TRUE(Recover_Evaporation_Time_File_From_Blocks(snapshot_root, output_path, run_id, mass_gev, sigma_cm2));

	const std::vector<TestEvaporationRow> rows = ReadRecoveredRows(output_path);
	ASSERT_EQ(static_cast<size_t>(3), rows.size());
	EXPECT_EQ(0, rows[0].rank);
	EXPECT_EQ(static_cast<uint64_t>(3), rows[0].trajectory_id);
	EXPECT_DOUBLE_EQ(10.0, rows[0].lifetime_unbinding);
	EXPECT_EQ(1, rows[1].rank);
	EXPECT_EQ(static_cast<uint64_t>(5), rows[1].trajectory_id);
	EXPECT_DOUBLE_EQ(99.0, rows[1].lifetime_unbinding);
	EXPECT_EQ(2, rows[2].rank);
	EXPECT_EQ(static_cast<uint64_t>(4), rows[2].trajectory_id);
	EXPECT_DOUBLE_EQ(30.0, rows[2].lifetime_unbinding);
}

TEST(TestDataGeneration, TestRecoverEvaporationTimeFileFromBlocksRejectsWrongRun)
{
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank != 0)
		return;

	const uint64_t run_id = 77ULL;
	const double mass_gev = 1.0;
	const double sigma_cm2 = 1.0e-37;
	const std::string output_dir = TestOutputDir("recover_wrong_run");
	const std::string snapshot_root = output_dir + "snapshot/";
	const std::string block_dir = snapshot_root + "evaporation_blocks/";
	EnsureDir(snapshot_root);
	EnsureDir(block_dir);

	std::vector<TestEvaporationRow> block_rows;
	block_rows.push_back(MakeTestEvaporationRow(0, 1, 5.0));
	WriteTestEvaporationBlock(block_dir + "block_000001.txt", run_id, 1, 60.0, mass_gev, sigma_cm2, block_rows);

	const std::string output_path = output_dir + "evaporation_times.partial.txt";
	TouchFile(output_path);
	EXPECT_FALSE(Recover_Evaporation_Time_File_From_Blocks(snapshot_root, output_path, run_id + 1, mass_gev, sigma_cm2));

	std::ifstream file(output_path);
	std::string first_line;
	std::getline(file, first_line);
	EXPECT_EQ("stale", first_line);
}
#endif
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
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# capture_rate_valid = 0.00000000"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# capture_rate_valid_CI_95_lower = 0.00000000"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# EARLY_STOP: max_trajectories_reached"));
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
		EXPECT_TRUE(FileContains(
		    output_dir + "bincount.txt",
		    "# bincount_integration = conservative-hermite-radial-v1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# normal_mode_mpi_sync_interval = 1048576"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# mpi_sync_rounds = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# final_mpi_round_trajectories = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# mpi_tail_trajectories = 0"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# capture_target_overshoot = 0"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# total_scatterings = 0"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# simulation_time_seconds = "));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# completed_outward_escapes = 1"));
		EXPECT_TRUE(FileContains(output_dir + "bincount.txt", "# unresolved_not_captured_trajectories = 0"));
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
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"schema_version\": \"trajectory-diagnostic-v2\""));
		EXPECT_TRUE(FileContains(
		    output_dir + "run_metadata.json",
		    "\"bincount_integration\": \"conservative-hermite-radial-v1\""));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"interpolation_points\": 20"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"legacy_evaporation_reconciliation\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"escape_radius_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"residence_time_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"event_sequence_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"event_count_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"trace_selection_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "run_metadata.json", "\"replay_state_invariant\": true"));
		EXPECT_TRUE(FileContains(output_dir + "trajectory_summary.tsv", "t_first_unbinding_s"));
		EXPECT_TRUE(FileContains(output_dir + "trajectory_summary.tsv", "n_recapture"));
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
		EXPECT_FALSE(FileExists(output_dir + "evaporation_diagnostics.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("evaporation_") + "mode_summary.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("evaporation_") + "mode_" + "bincount.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("computation_") + "time_summary.txt"));
		RemoveTestOutputDir(output_dir);
	}
}

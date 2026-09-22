#include "Parameter_Scan.hpp"

#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

#include "libphysica/Natural_Units.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

int main(int argc, char* argv[])
{
	int result = 0;

	::testing::InitGoogleTest(&argc, argv);
	MPI_Init(&argc, &argv);
	result = RUN_ALL_TESTS();
	MPI_Finalize();
	return result;
}

TEST(TestParameterScan, TestConfiguration)
{
	// ARRANGE
	// ACT & ASSERT
	Configuration cfg(PROJECT_DIR "tests/config_unittest.cfg", 1);
	// ASSERT
	EXPECT_TRUE(cfg.compute_halo_constraints);
	EXPECT_EQ(cfg.sample_size, 50);
	EXPECT_EQ(cfg.fixed_seed, 20260710u);
	EXPECT_DOUBLE_EQ(cfg.cross_section_min, 1.0e-35 * cm * cm);
	EXPECT_DOUBLE_EQ(cfg.cross_section_max, 1.0e-32 * cm * cm);
	EXPECT_EQ(cfg.cross_sections, 5);
	EXPECT_EQ(cfg.interpolation_points, 150);
	EXPECT_EQ(cfg.rate_radius_points, 150u);
	EXPECT_EQ(cfg.rate_speed_points, 150u);
	EXPECT_DOUBLE_EQ(cfg.rate_max_speed, 0.75);
	EXPECT_EQ(cfg.isoreflection_rings, 3);
	EXPECT_EQ(g_top_level_dir, "./unit_test_output/");
	EXPECT_TRUE(cfg.snapshot_config.enabled);
	EXPECT_DOUBLE_EQ(cfg.snapshot_config.interval_seconds, 10.0);
	const std::string physical_config = cfg.Physical_Configuration_JSON();
	EXPECT_NE(physical_config.find("\"DM_mass\":"), std::string::npos);
	EXPECT_NE(physical_config.find("\"DM_relative_couplings\":["), std::string::npos);
	EXPECT_NE(physical_config.find("\"DM_distribution\":\"SHM\""), std::string::npos);
	EXPECT_NE(physical_config.find("\"SHM_vObserver\":["), std::string::npos);
	EXPECT_EQ(physical_config.find("\"sample_size\":"), std::string::npos);
}

TEST(TestParameterScan, TestMinimalCaptureConfigurationDefaults)
{
	Configuration cfg(PROJECT_DIR "tests/config_capture_minimal.cfg", 1);

	EXPECT_EQ(cfg.run_mode, "Capture");
	EXPECT_TRUE(cfg.capture_mode);
	EXPECT_EQ(cfg.isoreflection_rings, 1);
	EXPECT_EQ(cfg.interpolation_points, 0);
	EXPECT_EQ(cfg.rate_radius_points, 0u);
	EXPECT_EQ(cfg.rate_speed_points, 0u);
	EXPECT_DOUBLE_EQ(cfg.rate_max_speed, 0.75);
	EXPECT_DOUBLE_EQ(cfg.cross_section_min, 0.0);
	EXPECT_DOUBLE_EQ(cfg.cross_section_max, 0.0);
	EXPECT_EQ(cfg.cross_sections, 0);
	EXPECT_FALSE(cfg.compute_halo_constraints);
	EXPECT_FALSE(cfg.perform_full_scan);
	EXPECT_EQ(cfg.fixed_seed, 0u);
}

TEST(TestParameterScan, TestRectangularRateGridConfiguration)
{
	const std::string path = "/tmp/damascus_rate_grid_" + std::to_string(getpid()) + ".cfg";
	{
		std::ifstream source(PROJECT_DIR "tests/config_unittest.cfg");
		ASSERT_TRUE(source.good());
		std::ofstream destination(path);
		destination << source.rdbuf()
		            << "\nrate_radius_points = 1000;\n"
		            << "rate_speed_points = 256;\n"
		            << "rate_max_speed = 0.02;\n";
		ASSERT_TRUE(destination.good());
	}
	Configuration cfg(path, 1);
	EXPECT_EQ(cfg.interpolation_points, 150u);
	EXPECT_EQ(cfg.rate_radius_points, 1000u);
	EXPECT_EQ(cfg.rate_speed_points, 256u);
	EXPECT_DOUBLE_EQ(cfg.rate_max_speed, 0.02);
	std::remove(path.c_str());
}

TEST(TestParameterScan, TestRejectsInvalidExplicitRateGridConfiguration)
{
	auto expect_invalid = [](const std::string& suffix, const std::string& settings) {
		const std::string path = "/tmp/damascus_invalid_rate_grid_"
		                       + std::to_string(getpid()) + "_" + suffix + ".cfg";
		{
			std::ifstream source(PROJECT_DIR "tests/config_capture_minimal.cfg");
			EXPECT_TRUE(source.good());
			std::ofstream destination(path);
			destination << source.rdbuf() << '\n' << settings;
			EXPECT_TRUE(destination.good());
		}
		EXPECT_THROW(Configuration cfg(path, 1), std::invalid_argument);
		std::remove(path.c_str());
	};

	expect_invalid("one_one",
	               "rate_radius_points = 1;\nrate_speed_points = 1;\n");
	expect_invalid("one_zero",
	               "rate_radius_points = 1;\nrate_speed_points = 0;\n");
	expect_invalid("zero_one",
	               "rate_radius_points = 0;\nrate_speed_points = 1;\n");
	expect_invalid("excessive_speed",
	               "rate_radius_points = 2;\nrate_speed_points = 2;\nrate_max_speed = 0.750001;\n");
}

TEST(TestParameterScan, TestLegacySinglePointRateGridStillDisablesInterpolation)
{
	const std::string path = "/tmp/damascus_legacy_rate_grid_"
	                       + std::to_string(getpid()) + ".cfg";
	{
		std::ifstream source(PROJECT_DIR "tests/config_capture_minimal.cfg");
		ASSERT_TRUE(source.good());
		std::ofstream destination(path);
		destination << source.rdbuf() << "\ninterpolation_points = 1;\n";
		ASSERT_TRUE(destination.good());
	}
	Configuration cfg(path, 1);
	EXPECT_EQ(cfg.interpolation_points, 1u);
	EXPECT_EQ(cfg.rate_radius_points, 1u);
	EXPECT_EQ(cfg.rate_speed_points, 1u);
	std::remove(path.c_str());
}

TEST(TestParameterScan, TestConfigurationSummary)
{
	// ARRANGE
	Configuration cfg(PROJECT_DIR "tests/config_unittest.cfg", 1);
	// ACT & ASSERT
	cfg.Print_Summary(0);
	// ASSERT
}

TEST(TestParameterScan, TestRejectsInvalidGridDefinitions)
{
	EXPECT_THROW(Parameter_Scan({}, {1.0 * cm * cm}, "empty_mass", 1), std::invalid_argument);
	EXPECT_THROW(Parameter_Scan({1.0 * MeV}, {}, "empty_coupling", 1), std::invalid_argument);
	EXPECT_THROW(Parameter_Scan({std::numeric_limits<double>::quiet_NaN()}, {1.0 * cm * cm}, "nan_mass", 1), std::invalid_argument);
	EXPECT_THROW(Parameter_Scan({1.0 * MeV}, {1.0 * cm * cm}, "bad_cl", 1, 0, 1.0), std::invalid_argument);
}

TEST(TestParameterScan, AtomicGridWriteRoundTripsAndPreservesOldFileOnFailure)
{
	const std::string root =
	    "/tmp/damascus_atomic_grid_" + std::to_string(getpid()) + "/";
	const std::string path = root + "P_Values_Grid.txt";
	ASSERT_EQ(0, mkdir(root.c_str(), 0755));

	const std::vector<std::vector<double>> expected{
	    {0.12345678901234566, -1.0}, {1.0, std::numeric_limits<double>::min()}};
	ASSERT_TRUE(Write_P_Value_Grid_Atomically(path, expected));
	{
		std::ifstream file(path);
		ASSERT_TRUE(file.good());
		double value = 0.0;
		for(const auto& row : expected)
			for(double expected_value : row)
			{
				ASSERT_TRUE(static_cast<bool>(file >> value));
				EXPECT_DOUBLE_EQ(expected_value, value);
			}
	}

	{
		std::ofstream file(path, std::ios::out | std::ios::trunc);
		ASSERT_TRUE(file.good());
		file << "stable resume grid\n";
	}
	const std::string tmp_path = path + ".tmp." + std::to_string(getpid());
	ASSERT_EQ(0, mkdir(tmp_path.c_str(), 0755));
	EXPECT_FALSE(Write_P_Value_Grid_Atomically(path, {{0.5}}));
	{
		std::ifstream file(path);
		ASSERT_TRUE(file.good());
		std::string contents;
		std::getline(file, contents);
		EXPECT_EQ("stable resume grid", contents);
	}

	ASSERT_EQ(0, rmdir(tmp_path.c_str()));
	ASSERT_EQ(0, std::remove(path.c_str()));
	ASSERT_EQ(0, rmdir(root.c_str()));
}

TEST(TestParameterScan, TestCriticalProbabilityDoesNotStallSquareTrace)
{
	const std::string previous_output_root = g_top_level_dir;
	const std::string root = "/tmp/damascus_sta_boundary_" + std::to_string(getpid()) + "/";
	const std::string results = root + "results/";
	const std::string run = results + "critical_probability/";
	mkdir(root.c_str(), 0755);
	mkdir(results.c_str(), 0755);
	mkdir(run.c_str(), 0755);

	const double certainty = 0.95;
	{
		std::ofstream file(run + "P_Values_Grid.txt");
		ASSERT_TRUE(file.good());
		file << std::setprecision(std::numeric_limits<double>::max_digits10)
		     << (1.0 - certainty) << "\n";
	}

	g_top_level_dir = root;
	Parameter_Scan scan({1.0 * MeV}, {1.0e-35 * cm * cm},
	                    "critical_probability", 1, 0, certainty, 1);
	std::vector<std::vector<double>> curve;
	EXPECT_NO_THROW(curve = scan.Limit_Curve());
	EXPECT_TRUE(curve.empty());
	g_top_level_dir = previous_output_root;

	std::remove((run + "P_Values_Grid.txt").c_str());
	rmdir(run.c_str());
	rmdir(results.c_str());
	rmdir(root.c_str());
}

TEST(TestParameterScan, TransportControlsAreTypedAndNeverSilentlyIgnored)
{
    std::ifstream source(PROJECT_DIR "tests/config_capture_minimal.cfg");
    const std::string text((std::istreambuf_iterator<char>(source)),std::istreambuf_iterator<char>());
    const std::string path="/tmp/damascus_transport_config_"+std::to_string(getpid())+".cfg";
    for(const auto& literal : {"3000", "3000.0", "3000L"}) {
        { std::ofstream out(path); out << text << "\nouter_removal_radius_rsun = " << literal << ";\n"; }
        Configuration cfg(path,1);
        EXPECT_DOUBLE_EQ(cfg.outer_removal_radius_rsun,3000.0);
        EXPECT_FALSE(cfg.diagnostic_mode);
    }
    for(const auto& invalid : {"outer_removal_radius_rsun = \"3000\";", "outer_removal_radius_rsun = 0.005;", "outer_boundary_radius_au = 1100;", "production_mode = true;", "production_mode = false;", "production_mode = 1;", "trajectory_events_enabled = false;", "thermal_validation_mode = \"false\";"}) {
        { std::ofstream out(path); out << text << "\n" << invalid << "\n"; }
        EXPECT_THROW(Configuration cfg(path,1),std::invalid_argument);
    }
    std::remove(path.c_str());
}

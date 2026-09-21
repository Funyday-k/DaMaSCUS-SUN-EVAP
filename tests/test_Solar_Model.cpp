#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mpi.h>
#include <random>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/DM_Particle_Standard.hpp"

#include "Solar_Model.hpp"

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

TEST(TestSolarModel, TestSolarNucleus)
{
	// ARRANGE
	double rho	= 1.505e+02 * gram / cm / cm / cm;
	double f_H1 = 0.36209;
	double n_H1 = f_H1 * rho / mProton;
	Solar_Model SSM;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(SSM.target_isotopes[0].Number_Density(0.0), n_H1);
	ASSERT_DOUBLE_EQ(SSM.target_isotopes[0].Number_Density(0.001 * rSun), n_H1);
	ASSERT_DOUBLE_EQ(SSM.target_isotopes[0].Number_Density(0.00150 * rSun), n_H1);
}

TEST(TestSolarModel, TestConstructor)
{
	// ARRANGE
	// ACT
	Solar_Model SSM;
	// ASSERT
	ASSERT_EQ(SSM.name, "Standard Solar Model AGSS09");
	ASSERT_NE(SSM.target_isotopes.size(), 0);
}

TEST(TestSolarModel, TestMass)
{
	// ARRANGE
	Solar_Model SSM;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(SSM.Mass(0), 0.0);
	ASSERT_DOUBLE_EQ(SSM.Mass(rSun), mSun);
	ASSERT_DOUBLE_EQ(SSM.Mass(0.00150 * rSun), 0.0000004 * mSun);
}

TEST(TestSolarModel, TestTemperature)
{
	// ARRANGE
	Solar_Model SSM;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(SSM.Temperature(0), 1.549e+07 * Kelvin);
	EXPECT_NEAR(SSM.Temperature(rSun), 5800.0 * Kelvin, 1.0e-20);
	EXPECT_DOUBLE_EQ(SSM.Temperature(0.96950 * rSun), 1.544e+05 * Kelvin);
}

TEST(TestSolarModel, TestLocalEscapeSpeed)
{
	// ARRANGE
	Solar_Model SSM;
	// ACT & ASSERT
	ASSERT_NEAR(SSM.Local_Escape_Speed(0), 1384.13 * km / sec, 0.1 * km / sec);
	ASSERT_DOUBLE_EQ(SSM.Local_Escape_Speed(rSun), sqrt(2.0 * G_Newton * mSun / rSun));
	ASSERT_DOUBLE_EQ(SSM.Local_Escape_Speed(2.0 * rSun), sqrt(2.0 * G_Newton * mSun / 2.0 / rSun));
}

TEST(TestSolarModel, TestNumberDensityNucleus)
{
	// ARRANGE
	Solar_Model SSM;
	double r1 = 0.0 * rSun;
	double r2 = 0.5 * rSun;
	double r3 = 1.5 * rSun;
	// ACT & ASSERT
	for(unsigned int i = 0; i < SSM.target_isotopes.size(); i++)
	{
		EXPECT_DOUBLE_EQ(SSM.Number_Density_Nucleus(r1, i), SSM.target_isotopes[i].Number_Density(r1));
		EXPECT_DOUBLE_EQ(SSM.Number_Density_Nucleus(r2, i), SSM.target_isotopes[i].Number_Density(r2));
		EXPECT_DOUBLE_EQ(SSM.Number_Density_Nucleus(r3, i), 0.0);
	}
}

TEST(TestSolarModel, TestNumberDensityElectron)
{
	// ARRANGE
	Solar_Model SSM;
	double r		 = 0.5 * rSun;
	double nElectron = 0.0;
	for(auto& isotope : SSM.target_isotopes)
		nElectron += isotope.Z * isotope.Number_Density(r);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(SSM.Number_Density_Electron(r), nElectron);
	EXPECT_DOUBLE_EQ(SSM.Number_Density_Electron(1.5 * rSun), 0.0);
}

TEST(TestSolarModel, TestDMScatteringRateElectron)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::DM_Particle_SI DM;
	DM.Set_Sigma_Electron(pb);
	double v_DM = 1e-3;
	double r1	= 0.5 * rSun;
	double r2	= 1.5 * rSun;
	// ACT & ASSERT
	EXPECT_GE(SSM.DM_Scattering_Rate_Electron(DM, r1, v_DM), 0.0);
	EXPECT_DOUBLE_EQ(SSM.DM_Scattering_Rate_Electron(DM, r2, v_DM), 0.0);
	DM.Set_Sigma_Electron(0.0);
	EXPECT_DOUBLE_EQ(SSM.DM_Scattering_Rate_Electron(DM, r1, v_DM), 0.0);
}

TEST(TestSolarModel, TestDMScatteringRateNucleus)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::DM_Particle_SI DM;
	DM.Set_Sigma_Proton(pb);
	double v_DM = 1e-3;
	double r1	= 0.5 * rSun;
	double r2	= 1.5 * rSun;
	// ACT & ASSERT
	EXPECT_GE(SSM.DM_Scattering_Rate_Nucleus(DM, r1, v_DM, 0), 0.0);
	EXPECT_DOUBLE_EQ(SSM.DM_Scattering_Rate_Nucleus(DM, r2, v_DM, 0), 0.0);
	DM.Unfix_Coupling_Ratios();
	DM.Set_Sigma_Proton(0.0);
	EXPECT_DOUBLE_EQ(SSM.DM_Scattering_Rate_Nucleus(DM, r1, v_DM, 0), 0.0);
}

TEST(TestSolarModel, TestTotalDMScatteringRate)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::DM_Particle_SI DM;
	DM.Set_Sigma_Proton(pb);
	double v_DM	 = 1e-3;
	double r1	 = 0.5 * rSun;
	double r2	 = 1.5 * rSun;
	double total = SSM.DM_Scattering_Rate_Electron(DM, r1, v_DM);
	for(unsigned int i = 0; i < SSM.target_isotopes.size(); i++)
		total += SSM.DM_Scattering_Rate_Nucleus(DM, r1, v_DM, i);
	// ACT & ASSERT
	EXPECT_GE(SSM.Total_DM_Scattering_Rate(DM, r1, v_DM), 0.0);
	EXPECT_DOUBLE_EQ(SSM.Total_DM_Scattering_Rate(DM, r2, v_DM), 0.0);
	EXPECT_DOUBLE_EQ(SSM.Total_DM_Scattering_Rate(DM, r1, v_DM), total);
	DM.Unfix_Coupling_Ratios();
	DM.Set_Sigma_Proton(0.0);
	DM.Set_Sigma_Neutron(0.0);
	DM.Set_Sigma_Electron(0.0);
	EXPECT_DOUBLE_EQ(SSM.Total_DM_Scattering_Rate(DM, r1, v_DM), 0.0);
}

TEST(TestSolarModel, TestTotalDMScatteringRateInterpolation)
{
	// ARRANGE
	int fixed_seed = 998;
	std::mt19937 PRNG(fixed_seed);
	Solar_Model SSM;
	obscura::DM_Particle_SI DM(0.01);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(pb);
	int trials		 = 500;
	double tolerance = 0.1;
	// ACT
	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 1000, 1000);
	// ASSERT
	for(int i = 0; i < trials; i++)
	{
		double r			 = libphysica::Sample_Uniform(PRNG, 0, rSun);
		double w			 = libphysica::Sample_Uniform(PRNG, 0, 0.3);
		double correct_value = SSM.Total_DM_Scattering_Rate_Computed(DM, r, w);
		double max_deviation = tolerance * correct_value;
		ASSERT_NEAR(SSM.Total_DM_Scattering_Rate(DM, r, w), correct_value, max_deviation);
	}
}

TEST(TestSolarModel, TestTotalDMScatteringRateRegularGridNodes)
{
	Solar_Model SSM;
	obscura::DM_Particle_SI DM(0.01);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(pb);
	const unsigned int radius_points = 31;
	const unsigned int speed_points = 29;
	const double max_speed = 0.02;
	const unsigned int radius_indices[] = {0, 7, 15, 30};
	const unsigned int speed_indices[] = {0, 5, 14, 28};

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, radius_points, speed_points, max_speed);
	EXPECT_EQ(SSM.Scattering_Rate_Interpolation_Radius_Points(), radius_points);
	EXPECT_EQ(SSM.Scattering_Rate_Interpolation_Speed_Points(), speed_points);
	EXPECT_DOUBLE_EQ(SSM.Scattering_Rate_Interpolation_Max_Speed(), max_speed);
	for(const auto radius_index : radius_indices)
	{
		const double radius = rSun * radius_index / (radius_points - 1);
		for(const auto speed_index : speed_indices)
		{
			const double speed = max_speed * speed_index / (speed_points - 1);
			const double expected = SSM.Total_DM_Scattering_Rate_Computed(DM, radius, speed);
			const double tolerance = std::max(1.0e-30, 1.0e-10 * std::fabs(expected));
			EXPECT_NEAR(SSM.Total_DM_Scattering_Rate(DM, radius, speed), expected, tolerance);
		}
	}

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 1, speed_points);
	EXPECT_EQ(SSM.Scattering_Rate_Interpolation_Radius_Points(), 0u);
	EXPECT_EQ(SSM.Scattering_Rate_Interpolation_Speed_Points(), 0u);
	const double radius = 0.37 * rSun;
	const double speed = 0.21;
	EXPECT_DOUBLE_EQ(SSM.Total_DM_Scattering_Rate(DM, radius, speed),
	                 SSM.Total_DM_Scattering_Rate_Computed(DM, radius, speed));
}

TEST(TestSolarModel, TestScatteringRateFallbackIsCountedWithoutWarning)
{
	Solar_Model SSM;
	obscura::DM_Particle_SD DM(0.1 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Fix_Coupling_Ratio(1.0, 0.0);
	DM.Set_Sigma_Proton(1.0e-32 * cm * cm);
	DM.Set_Sigma_Electron(1.0e-80 * cm * cm);
	const double max_speed = 0.02;
	const double query_speed = 1.1 * max_speed;
	const double radius = 0.4 * rSun;

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 31, 17, max_speed);
	const double expected = SSM.Total_DM_Scattering_Rate_Computed(DM, radius, query_speed);
	testing::internal::CaptureStderr();
	const double actual = SSM.Total_DM_Scattering_Rate(DM, radius, query_speed);
	const std::string stderr_output = testing::internal::GetCapturedStderr();

	EXPECT_DOUBLE_EQ(actual, expected);
	EXPECT_TRUE(stderr_output.empty());
	EXPECT_EQ(SSM.Scattering_Rate_Query_Count(), 1u);
	EXPECT_EQ(SSM.Scattering_Rate_Fallback_Count(), 1u);
	EXPECT_DOUBLE_EQ(SSM.Scattering_Rate_Fallback_Fraction(), 1.0);
	EXPECT_DOUBLE_EQ(SSM.Maximum_Scattering_Rate_Query_Speed(), query_speed);
}

TEST(TestSolarModel, TestRectangularGridAccuracyAcrossRepresentativeMasses)
{
	const double masses[] = {0.01, 0.1, 1.0};
	const double max_speeds[] = {0.08, 0.02, 0.01};
	const double radius_fractions[] = {0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99};
	const double speed_fractions[] = {0.0, 0.1, 0.25, 0.5, 0.75, 0.95};

	for(unsigned int mass_index = 0; mass_index < 3; mass_index++)
	{
		Solar_Model SSM;
		obscura::DM_Particle_SD DM(masses[mass_index] * GeV);
		DM.Set_Low_Mass_Mode(true);
		DM.Fix_Coupling_Ratio(1.0, 0.0);
		DM.Set_Sigma_Proton(1.0e-32 * cm * cm);
		DM.Set_Sigma_Electron(1.0e-80 * cm * cm);
		SSM.Interpolate_Total_DM_Scattering_Rate(
		    DM, 1000, 256, max_speeds[mass_index]);

		double weighted_absolute_error = 0.0;
		double total_exact_rate = 0.0;
		for(const double radius_fraction : radius_fractions)
			for(const double speed_fraction : speed_fractions)
			{
				const double radius = radius_fraction * rSun;
				const double speed = speed_fraction * max_speeds[mass_index];
				const double exact = SSM.Total_DM_Scattering_Rate_Computed(DM, radius, speed);
				const double interpolated = SSM.Total_DM_Scattering_Rate(DM, radius, speed);
				weighted_absolute_error += std::fabs(interpolated - exact);
				total_exact_rate += exact;
			}
		ASSERT_GT(total_exact_rate, 0.0);
		EXPECT_LT(weighted_absolute_error / total_exact_rate, 1.0e-3)
		    << "mass_GeV=" << masses[mass_index];
		EXPECT_EQ(SSM.Scattering_Rate_Fallback_Count(), 0u);
	}
}

TEST(TestSolarModel, TestTotalDMScatteringRateRegularGridSDProton)
{
	Solar_Model SSM;
	obscura::DM_Particle_SD DM(4.0 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Fix_Coupling_Ratio(1.0, 0.0);
	DM.Set_Sigma_Proton(1.0e-36 * cm * cm);
	const unsigned int grid_points = 17;

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, grid_points, grid_points);
	const double radius = rSun * 9.0 / (grid_points - 1);
	const double speed = 0.75 * 11.0 / (grid_points - 1);
	const double expected = SSM.Total_DM_Scattering_Rate_Computed(DM, radius, speed);
	const double tolerance = std::max(1.0e-30, 1.0e-10 * std::fabs(expected));
	EXPECT_NEAR(SSM.Total_DM_Scattering_Rate(DM, radius, speed), expected, tolerance);
}

TEST(TestSolarModel, TestNonFiniteInterpolationCoordinatesAreRejectedBeforeIndexing)
{
	Solar_Model SSM;
	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(pb);
	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 5, 5);

	const double nan = std::numeric_limits<double>::quiet_NaN();
	EXPECT_TRUE(std::isnan(SSM.Total_DM_Scattering_Rate(DM, nan, 1.0e-3)));
	EXPECT_TRUE(std::isnan(SSM.Total_DM_Scattering_Rate_Interpolated(DM, 0.5 * rSun, nan)));
}

TEST(TestSolarModel, TestPrintSummary)
{
	// ARRANGE
	Solar_Model SSM;

	// ACT & ASSERT
	SSM.Print_Summary();
}

// 3. Thermal average of relative speed between a particle of speed v_DM and a solar thermal target.
TEST(TestSolarModel, TestThermalAveragedRelativeSpeed)
{
	// ARRANGE
	double temperature = 1.0e7 * Kelvin;
	double mass_target = mProton;
	double v_DM		   = 1e-3;
	// ACT & ASSERT
	EXPECT_NEAR(Thermal_Averaged_Relative_Speed(0.01 * Kelvin, mass_target, v_DM), v_DM, 1e-10);
	EXPECT_NEAR(Thermal_Averaged_Relative_Speed(temperature, mass_target, v_DM), 0.0017928, 1.0e-7);
	const double zero_speed_average = Thermal_Averaged_Relative_Speed(temperature, mass_target, 0.0);
	EXPECT_TRUE(std::isfinite(zero_speed_average));
	EXPECT_GT(zero_speed_average, 0.0);
	EXPECT_TRUE(std::isnan(Thermal_Averaged_Relative_Speed(0.0, mass_target, v_DM)));
	EXPECT_TRUE(std::isnan(Thermal_Averaged_Relative_Speed(temperature, mass_target, -v_DM)));
}

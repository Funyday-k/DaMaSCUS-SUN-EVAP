#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"

#include "Dark_Photon.hpp"
#include "Simulation_Trajectory.hpp"
#include "Simulation_Utilities.hpp"
#include "Solar_Model.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

namespace
{
struct OrbitError
{
	double position;
	double velocity;
};

struct OrbitBincount
{
	std::array<double, NUM_BINS> dt{};
	std::array<double, NUM_BINS> v2dt{};
};

void AddBincountContributions(
	OrbitBincount& histogram,
	const std::vector<BincountContribution>& contributions)
{
	for(const BincountContribution& contribution : contributions)
	{
		ASSERT_GE(contribution.bin, 0);
		ASSERT_LT(contribution.bin, NUM_BINS);
		histogram.dt[contribution.bin] += contribution.dt_sec;
		histogram.v2dt[contribution.bin] += contribution.v2dt_km2_per_sec;
	}
}

OrbitBincount EllipticOutboundBincount(double step, double first_step)
{
	const double semi_major_axis = 1.5 * rSun;
	const double eccentricity = 1.0 / 3.0;
	const double mu = G_Newton * mSun;
	const double mean_motion = std::sqrt(mu / std::pow(semi_major_axis, 3.0));
	const double duration = M_PI / mean_motion;
	const double periapsis = semi_major_axis * (1.0 - eccentricity);
	const double periapsis_speed = std::sqrt(
	    mu * (1.0 + eccentricity) / (semi_major_axis * (1.0 - eccentricity)));
	Event initial(
	    0.0,
	    libphysica::Vector({periapsis, 0.0, 0.0}),
	    libphysica::Vector({0.0, periapsis_speed, 0.0}));
	Free_Particle_Propagator propagator(initial);
	OrbitBincount histogram;
	std::vector<BincountContribution> contributions;
	unsigned int steps = 0;
	while(propagator.Current_Time() < duration * (1.0 - 1.0e-14))
	{
		const Event before = propagator.Event_In_3D();
		const double requested_step = (steps == 0 && first_step > 0.0) ? first_step : step;
		propagator.time_step = std::min(
		    requested_step, duration - propagator.Current_Time());
		EXPECT_TRUE(propagator.Runge_Kutta_45_Step(mSun));
		const Event after = propagator.Event_In_3D();
		Compute_Bincount_Interval_Contributions(before, after, contributions);
		AddBincountContributions(histogram, contributions);
		steps++;
		if(steps > 10000)
		{
			ADD_FAILURE() << "Elliptic bincount validation exceeded its step guard";
			break;
		}
	}
	return histogram;
}

OrbitBincount LegacyLeftEndpointEllipticBincount(double step)
{
	const double semi_major_axis = 1.5 * rSun;
	const double eccentricity = 1.0 / 3.0;
	const double mu = G_Newton * mSun;
	const double mean_motion = std::sqrt(mu / std::pow(semi_major_axis, 3.0));
	const double duration = M_PI / mean_motion;
	const double periapsis = semi_major_axis * (1.0 - eccentricity);
	const double periapsis_speed = std::sqrt(
	    mu * (1.0 + eccentricity) / (semi_major_axis * (1.0 - eccentricity)));
	Free_Particle_Propagator propagator(Event(
	    0.0,
	    libphysica::Vector({periapsis, 0.0, 0.0}),
	    libphysica::Vector({0.0, periapsis_speed, 0.0})));
	OrbitBincount histogram;
	while(propagator.Current_Time() < duration * (1.0 - 1.0e-14))
	{
		const Event before = propagator.Event_In_3D();
		propagator.time_step = std::min(
		    step, duration - propagator.Current_Time());
		EXPECT_TRUE(propagator.Runge_Kutta_45_Step(mSun));
		const Event after = propagator.Event_In_3D();
		const double radius_km = In_Units(before.Radius(), km);
		const double speed_km_s = In_Units(before.Speed(), km / sec);
		const double dt_sec = In_Units(after.time - before.time, sec);
		if(radius_km >= 0.0 && radius_km < BIN_MAX_KM)
		{
			const int bin = static_cast<int>(radius_km / BIN_WIDTH_KM);
			histogram.dt[bin] += dt_sec;
			histogram.v2dt[bin] += speed_km_s * speed_km_s * dt_sec;
		}
	}
	return histogram;
}

OrbitBincount ExactEllipticOutboundBincount()
{
	const double semi_major_axis_km = 1.5 * R_SUN_KM;
	const double eccentricity = 1.0 / 3.0;
	const double mu_km3_s2 = In_Units(
	    G_Newton * mSun, km * km * km / sec / sec);
	const double mean_motion_s =
	    std::sqrt(mu_km3_s2 / std::pow(semi_major_axis_km, 3.0));
	OrbitBincount histogram;
	for(int bin = 1000; bin < NUM_BINS; bin++)
	{
		const double radius_low_km =
		    std::max(R_SUN_KM, static_cast<double>(bin) * BIN_WIDTH_KM);
		const double radius_high_km =
		    std::min(BIN_MAX_KM, static_cast<double>(bin + 1) * BIN_WIDTH_KM);
		if(!(radius_high_km > radius_low_km))
			continue;
		const double cosine_low = std::max(
		    -1.0, std::min(1.0,
		                  (1.0 - radius_low_km / semi_major_axis_km) / eccentricity));
		const double cosine_high = std::max(
		    -1.0, std::min(1.0,
		                  (1.0 - radius_high_km / semi_major_axis_km) / eccentricity));
		const double anomaly_low = std::acos(cosine_low);
		const double anomaly_high = std::acos(cosine_high);
		const double mean_anomaly_low =
		    anomaly_low - eccentricity * std::sin(anomaly_low);
		const double mean_anomaly_high =
		    anomaly_high - eccentricity * std::sin(anomaly_high);
		histogram.dt[bin] =
		    (mean_anomaly_high - mean_anomaly_low) / mean_motion_s;
		histogram.v2dt[bin] =
		    mu_km3_s2 / (semi_major_axis_km * mean_motion_s)
		    * (anomaly_high - anomaly_low
		       + eccentricity * (std::sin(anomaly_high) - std::sin(anomaly_low)));
	}
	return histogram;
}

double RelativeL1(
	const std::array<double, NUM_BINS>& actual,
	const std::array<double, NUM_BINS>& expected)
{
	double absolute_error = 0.0;
	double expected_total = 0.0;
	for(int bin = 0; bin < NUM_BINS; bin++)
	{
		absolute_error += std::fabs(actual[bin] - expected[bin]);
		expected_total += std::fabs(expected[bin]);
	}
	return absolute_error / expected_total;
}

OrbitError EllipticOrbitError(double step, double duration)
{
	const double semi_major_axis = 2.0 * rSun;
	const double eccentricity = 0.3;
	const double mu = G_Newton * mSun;
	const double mean_motion = std::sqrt(mu / std::pow(semi_major_axis, 3.0));
	const double periapsis = semi_major_axis * (1.0 - eccentricity);
	const double periapsis_speed = std::sqrt(
	    mu * (1.0 + eccentricity) / (semi_major_axis * (1.0 - eccentricity)));
	Event initial(
	    0.0,
	    libphysica::Vector({periapsis, 0.0, 0.0}),
	    libphysica::Vector({0.0, periapsis_speed, 0.0}));
	Free_Particle_Propagator propagator(initial);

	const unsigned int nominal_steps = static_cast<unsigned int>(std::ceil(duration / step));
	unsigned int steps = 0;
	while(propagator.Current_Time() < duration * (1.0 - 1.0e-14))
	{
		// Hold the step fixed for this convergence experiment. Production uses
		// the same RKF45 stepper with adaptive error control. If the requested
		// cap is rejected, keep the accepted smaller step and continue to the
		// common physical end time.
		propagator.time_step = std::min(step, duration - propagator.Current_Time());
		EXPECT_TRUE(propagator.Runge_Kutta_45_Step(mSun));
		steps++;
		if(steps > 4u * nominal_steps)
		{
			ADD_FAILURE() << "RKF45 required too many rejected/subdivided steps";
			break;
		}
	}

	const Event final = propagator.Event_In_3D();
	const double mean_anomaly = mean_motion * duration;
	double eccentric_anomaly = mean_anomaly;
	for(int iteration = 0; iteration < 20; iteration++)
	{
		const double residual = eccentric_anomaly
		                        - eccentricity * std::sin(eccentric_anomaly)
		                        - mean_anomaly;
		const double derivative = 1.0 - eccentricity * std::cos(eccentric_anomaly);
		eccentric_anomaly -= residual / derivative;
	}
	const double vertical_factor = std::sqrt(1.0 - eccentricity * eccentricity);
	const double anomaly_rate = mean_motion
	                            / (1.0 - eccentricity * std::cos(eccentric_anomaly));
	const libphysica::Vector expected_position(
	    {semi_major_axis * (std::cos(eccentric_anomaly) - eccentricity),
	     semi_major_axis * vertical_factor * std::sin(eccentric_anomaly),
	     0.0});
	const libphysica::Vector expected_velocity(
	    {-semi_major_axis * std::sin(eccentric_anomaly) * anomaly_rate,
	     semi_major_axis * vertical_factor * std::cos(eccentric_anomaly) * anomaly_rate,
	     0.0});

	EXPECT_NEAR(final.time, duration, 1.0e-12 * duration);
	return {(final.position - expected_position).Norm(),
	        (final.velocity - expected_velocity).Norm()};
}

double UniformKSDistance(std::vector<double> samples)
{
	std::sort(samples.begin(), samples.end());
	double distance = 0.0;
	const double count = static_cast<double>(samples.size());
	for(std::size_t i = 0; i < samples.size(); i++)
	{
		const double lower_empirical = static_cast<double>(i) / count;
		const double upper_empirical = static_cast<double>(i + 1) / count;
		distance = std::max(distance, std::fabs(samples[i] - lower_empirical));
		distance = std::max(distance, std::fabs(upper_empirical - samples[i]));
	}
	return distance;
}

} // namespace

TEST(PhysicsValidation, PointMassEllipticOrbitShowsFourthOrderConvergence)
{
	const double duration = 4000.0 * sec;
	const OrbitError coarse = EllipticOrbitError(400.0 * sec, duration);
	const OrbitError medium = EllipticOrbitError(200.0 * sec, duration);
	const OrbitError fine = EllipticOrbitError(100.0 * sec, duration);

	// A fourth-order global error should fall by about 2^4 when the fixed
	// step is halved. The lower bound leaves room for roundoff and the polar
	// reconstruction while still detecting a loss of RK order.
	EXPECT_GT(coarse.position / medium.position, 8.0);
	EXPECT_GT(medium.position / fine.position, 8.0);
	EXPECT_GT(coarse.velocity / medium.velocity, 8.0);
	EXPECT_GT(medium.velocity / fine.velocity, 8.0);
	EXPECT_LT(In_Units(fine.position, km), 0.01);
	EXPECT_LT(In_Units(fine.velocity, km / sec), 1.0e-5);
}

TEST(PhysicsValidation, BincountMatchesAnalyticKeplerShellResidenceTimes)
{
	const OrbitBincount exact = ExactEllipticOutboundBincount();
	const OrbitBincount coarse =
	    EllipticOutboundBincount(400.0 * sec, 0.0);
	const OrbitBincount shifted =
	    EllipticOutboundBincount(400.0 * sec, 137.0 * sec);
	const OrbitBincount fine =
	    EllipticOutboundBincount(100.0 * sec, 0.0);
	const OrbitBincount legacy =
	    LegacyLeftEndpointEllipticBincount(400.0 * sec);

	EXPECT_LT(RelativeL1(coarse.dt, exact.dt), 5.0e-3);
	EXPECT_LT(RelativeL1(shifted.dt, exact.dt), 5.0e-3);
	EXPECT_LT(RelativeL1(fine.dt, exact.dt), 3.0e-3);
	EXPECT_LT(RelativeL1(coarse.v2dt, exact.v2dt), 5.0e-3);
	EXPECT_LT(RelativeL1(shifted.v2dt, exact.v2dt), 5.0e-3);
	EXPECT_LT(RelativeL1(fine.v2dt, exact.v2dt), 3.0e-3);
	EXPECT_LT(RelativeL1(coarse.dt, shifted.dt), 5.0e-3);
	EXPECT_GT(RelativeL1(legacy.dt, exact.dt), 5.0e-1);
}

TEST(PhysicsValidation, SolarProfilesObeyBasicPhysicalStructure)
{
	Solar_Model solar_model;
	const unsigned int samples = 1000;
	double previous_mass = solar_model.Mass(0.0);
	double previous_escape_speed = solar_model.Local_Escape_Speed(0.0);

	for(unsigned int i = 0; i <= samples; i++)
	{
		const double radius = rSun * static_cast<double>(i) / samples;
		const double enclosed_mass = solar_model.Mass(radius);
		const double temperature = solar_model.Temperature(radius);
		const double density = solar_model.Mass_Density(radius);
		const double escape_speed = solar_model.Local_Escape_Speed(radius);

		EXPECT_TRUE(std::isfinite(enclosed_mass));
		EXPECT_TRUE(std::isfinite(temperature));
		EXPECT_TRUE(std::isfinite(density));
		EXPECT_TRUE(std::isfinite(escape_speed));
		EXPECT_GE(enclosed_mass, previous_mass);
		EXPECT_GT(temperature, 0.0);
		EXPECT_GE(density, 0.0);
		EXPECT_LE(escape_speed, previous_escape_speed);

		previous_mass = enclosed_mass;
		previous_escape_speed = escape_speed;
	}

	EXPECT_DOUBLE_EQ(solar_model.Mass(0.0), 0.0);
	EXPECT_DOUBLE_EQ(solar_model.Mass(rSun), mSun);
	EXPECT_DOUBLE_EQ(
	    solar_model.Local_Escape_Speed(3.0 * rSun),
	    std::sqrt(2.0 * G_Newton * mSun / (3.0 * rSun)));
}

TEST(PhysicsValidation, InitialImpactParameterIsUniformInArea)
{
	Solar_Model solar_model;
	obscura::Standard_Halo_Model halo_model;
	std::mt19937 prng(20260714u);
	const unsigned int sample_count = 8000;
	std::vector<double> normalized_area;
	normalized_area.reserve(sample_count);

	double first_moment = 0.0;
	double second_moment = 0.0;
	for(unsigned int i = 0; i < sample_count; i++)
	{
		const Event event = Initial_Conditions(halo_model, solar_model, prng);
		const double maximum_angular_momentum = rSun * std::sqrt(
		    event.Speed() * event.Speed()
		    + std::pow(solar_model.Local_Escape_Speed(rSun), 2.0));
		const double area_coordinate = std::pow(
		    event.Angular_Momentum() / maximum_angular_momentum, 2.0);

		ASSERT_TRUE(std::isfinite(area_coordinate));
		ASSERT_GE(area_coordinate, 0.0);
		ASSERT_LE(area_coordinate, 1.0 + 1.0e-12);
		normalized_area.push_back(std::min(1.0, area_coordinate));
		first_moment += area_coordinate;
		second_moment += area_coordinate * area_coordinate;
	}

	first_moment /= sample_count;
	second_moment /= sample_count;
	EXPECT_NEAR(first_moment, 0.5, 0.01);
	EXPECT_NEAR(second_moment, 1.0 / 3.0, 0.01);
	// 1.75/sqrt(N) is a conservative deterministic threshold near the 0.5%
	// Kolmogorov-Smirnov critical value for a truly uniform sample.
	EXPECT_LT(UniformKSDistance(normalized_area), 1.75 / std::sqrt(sample_count));
}

TEST(PhysicsValidation, DarkPhotonAngleSamplerMatchesItsCDF)
{
	DM_Particle_Dark_Photon dark_matter(0.1 * GeV, 1.0e-36 * cm * cm);
	dark_matter.Set_Low_Mass_Mode(true);
	dark_matter.Set_FormFactor_DM("General", 10.0 * keV);
	std::mt19937 prng(314159u);
	const double speed = 1.0e-3;
	const unsigned int sample_count = 10000;
	std::vector<double> probability_integral_transform;
	probability_integral_transform.reserve(sample_count);

	for(unsigned int i = 0; i < sample_count; i++)
	{
		const double cosine = dark_matter.Sample_Scattering_Angle_Electron(prng, speed);
		ASSERT_GE(cosine, -1.0);
		ASSERT_LE(cosine, 1.0);
		probability_integral_transform.push_back(
		    dark_matter.CDF_Scattering_Angle_Electron(cosine, speed));
	}

	EXPECT_LT(
	    UniformKSDistance(probability_integral_transform),
	    1.75 / std::sqrt(sample_count));

	// Independently integrate the PDF with a midpoint rule.
	const unsigned int integration_bins = 20000;
	double normalization = 0.0;
	const double width = 2.0 / integration_bins;
	for(unsigned int i = 0; i < integration_bins; i++)
	{
		const double cosine = -1.0 + (i + 0.5) * width;
		normalization += width * dark_matter.PDF_Scattering_Angle_Electron(cosine, speed);
	}
	EXPECT_NEAR(normalization, 1.0, 1.0e-6);
}

TEST(PhysicsValidation, ThermalRelativeSpeedMatchesMaxwellLimitAtRest)
{
	const double temperature = 1.2e7 * Kelvin;
	const double expected_mean_speed = std::sqrt(8.0 * temperature / (M_PI * mProton));
	EXPECT_NEAR(
	    Thermal_Averaged_Relative_Speed(temperature, mProton, 0.0),
	    expected_mean_speed,
	    1.0e-12 * expected_mean_speed);
}

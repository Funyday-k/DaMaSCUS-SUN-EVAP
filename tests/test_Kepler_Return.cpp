#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <random>

#include "Simulation_Trajectory.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

namespace
{
libphysica::Vector EccentricityVector(const Event& event)
{
	return event.velocity.Cross(event.position.Cross(event.velocity)) / (G_Newton * mSun)
	    - event.position / event.Radius();
}

// Construct the initial and final states from orbital elements, independently
// of the propagation routine's reconstruction from Cartesian invariants.
void CheckEllipticReturn(double periapsis, double apoapsis,
                        const libphysica::Vector& axis_x,
                        const libphysica::Vector& axis_y)
{
	const double radius = 1.2 * rSun;
	const double mu = G_Newton * mSun;
	const double a = 0.5 * (periapsis + apoapsis);
	const double e = (apoapsis - periapsis) / (apoapsis + periapsis);
	const double cos_e = (1.0 - radius / a) / e;
	const double sin_e = std::sqrt(1.0 - cos_e * cos_e);
	const double anomaly = std::acos(cos_e);
	const double beta = std::sqrt(1.0 - e * e);
	const double n = std::sqrt(mu / (a * a * a));
	const double velocity_scale = n * a / (1.0 - e * cos_e);
	Event outward(3.0 * sec,
	    a * ((cos_e - e) * axis_x + beta * sin_e * axis_y),
	    velocity_scale * (-sin_e * axis_x + beta * cos_e * axis_y));
	const libphysica::Vector expected_position =
	    a * ((cos_e - e) * axis_x - beta * sin_e * axis_y);
	const libphysica::Vector expected_velocity =
	    velocity_scale * (sin_e * axis_x + beta * cos_e * axis_y);

	BoundKeplerExteriorArc arc;
	ASSERT_TRUE(Compute_Bound_Kepler_Exterior_Arc(outward, arc));
	ASSERT_FALSE(arc.outer_domain_removed);
	const Event& inbound = arc.terminal_event;
	EXPECT_LT((inbound.position - expected_position).Norm() / radius, 2.0e-9);
	EXPECT_LT((inbound.velocity - expected_velocity).Norm() / outward.Speed(), 2.0e-9);
	EXPECT_LT((EccentricityVector(inbound) - EccentricityVector(outward)).Norm(), 2.0e-12);
	EXPECT_LT((inbound.position.Cross(inbound.velocity)
	           - outward.position.Cross(outward.velocity)).Norm()
	          / outward.Angular_Momentum(), 2.0e-10);
	EXPECT_NEAR(inbound.Speed() / outward.Speed(), 1.0, 2.0e-12);
	EXPECT_LT(inbound.position.Dot(inbound.velocity), 0.0);
	const double expected_time = In_Units(
	    (2.0 * M_PI - 2.0 * anomaly + 2.0 * e * sin_e) / n, sec);
	EXPECT_NEAR(arc.elapsed_time_sec / expected_time, 1.0, 2.0e-9);
}
}

TEST(KeplerReturn, CartesianStateMatchesIndependentOrbitalElements)
{
	std::mt19937 generator(20260910);
	std::uniform_real_distribution<double> unit(-1.0, 1.0);
	std::uniform_real_distribution<double> uniform(0.0, 1.0);
	for(unsigned int sample = 0; sample < 1200; sample++)
	{
		SCOPED_TRACE(sample);
		const libphysica::Vector axis_x =
		    libphysica::Vector({unit(generator), unit(generator), unit(generator)}).Normalized();
		const libphysica::Vector reference = (std::fabs(axis_x[2]) < 0.9)
		    ? libphysica::Vector({0.0, 0.0, 1.0})
		    : libphysica::Vector({1.0, 0.0, 0.0});
		const libphysica::Vector axis_y = axis_x.Cross(reference).Normalized();
		CheckEllipticReturn((0.02 + 1.1 * uniform(generator)) * rSun,
		                   (1.3 + 18.7 * uniform(generator)) * rSun, axis_x, axis_y);
	}
}

TEST(KeplerReturn, EccentricAndNearlyCircularOrbitsKeepTheirApsidalDirection)
{
	const libphysica::Vector axis_x({1.0, 0.0, 0.0});
	const libphysica::Vector axis_y({0.0, 1.0, 0.0});
	CheckEllipticReturn(1.199 * rSun, 1.201 * rSun, axis_x, axis_y);
	CheckEllipticReturn(1.0e-7 * rSun, 4.0 * rSun, axis_x, axis_y);
}

TEST(KeplerReturn, DegenerateRadialAndUnboundStatesAreRejected)
{
	const libphysica::Vector position({1.2 * rSun, 0.0, 0.0});
	BoundKeplerExteriorArc arc;
	EXPECT_FALSE(Compute_Bound_Kepler_Exterior_Arc(
	    Event(0.0, position, libphysica::Vector({100.0 * km / sec, 0.0, 0.0})), arc));
	EXPECT_FALSE(Compute_Bound_Kepler_Exterior_Arc(
	    Event(0.0, position, libphysica::Vector({1000.0 * km / sec, 0.0, 0.0})), arc));
	EXPECT_FALSE(Compute_Bound_Kepler_Exterior_Arc(
	    Event(0.0, position, libphysica::Vector({-100.0 * km / sec, 200.0 * km / sec, 0.0})), arc));
}

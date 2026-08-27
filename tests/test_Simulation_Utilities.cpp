#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "obscura/Astronomy.hpp"
#include "obscura/DM_Halo_Models.hpp"

#include "Simulation_Trajectory.hpp"
#include "Simulation_Utilities.hpp"
#include "Solar_Model.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

namespace
{
double Kepler_Energy(const Event& event)
{
	const double radius = event.Radius();
	const double speed = event.Speed();
	return 0.5 * speed * speed - G_Newton * mSun / radius;
}

double Radial_Velocity_For_Test(const Event& event)
{
	return event.position.Dot(event.velocity) / event.Radius();
}

void Expect_Relative_Near(double value, double reference, double rel_tol)
{
	const double scale = std::max(std::fabs(reference), 1.0e-300);
	EXPECT_NEAR(value, reference, rel_tol * scale);
}

}	// namespace

// 1. Event class
TEST(TestSimulationUtilities, TestEventConstructor)
{
	// ARRANGE
	Event default_event;
	double t = 7 * sec;
	libphysica::Vector x({1 * meter, 2 * meter, 3 * meter});
	libphysica::Vector v({4 * meter / sec, 5 * meter / sec, 6 * meter / sec});
	Event event(t, x, v);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(default_event.time, 0.0);
	for(int i = 0; i < 3; i++)
		EXPECT_DOUBLE_EQ(default_event.position[i], libphysica::Vector({0, 0, 0})[i]);
	for(int i = 0; i < 3; i++)
		EXPECT_DOUBLE_EQ(default_event.velocity[i], libphysica::Vector({0, 0, 0})[i]);
	EXPECT_DOUBLE_EQ(In_Units(event.time, sec), 7.0);
	for(int i = 0; i < 3; i++)
		EXPECT_DOUBLE_EQ(event.position[i], x[i]);
	for(int i = 0; i < 3; i++)
		EXPECT_DOUBLE_EQ(event.velocity[i], v[i]);
}

TEST(TestSimulationUtilities, TestRadius)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector x({1, 2, 3});
	libphysica::Vector v({4, 5, 6});
	Event event(t, x, v);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(event.Radius(), sqrt(14));
}

TEST(TestSimulationUtilities, TestSpeed)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector x({1, 2, 3});
	libphysica::Vector v({4, 5, 6});
	Event event(t, x, v);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(event.Speed(), sqrt(77));
}

TEST(TestSimulationUtilities, TestAngularMomentum)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector x({1, 2, 3});
	libphysica::Vector v({4, 5, 6});
	Event event(t, x, v);
	double J = x.Cross(v).Norm();
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(event.Angular_Momentum(), J);
}

TEST(TestSimulationUtilities, TestAsymptoticSpeedSqr)
{
	// ARRANGE
	Solar_Model SSM;
	double r	= 5.0 * rSun;
	double u	= 100.0 * km / sec;
	double vesc = SSM.Local_Escape_Speed(r);
	double v	= sqrt(u * u + vesc * vesc);

	libphysica::Vector x({0, r, 0});
	libphysica::Vector vel({v, 0, 0});
	Event event(0, x, vel);

	double tol = 1.0e-10;
	// ACT & ASSERT
	ASSERT_NEAR(event.Asymptotic_Speed_Sqr(SSM), u * u, tol);
}

TEST(TestSimulationUtilities, TestIsoreflectionAngle)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector x({1, 2, 3});
	libphysica::Vector v({4, 5, 6});
	Event event(t, x, v);

	libphysica::Vector vSun = obscura::Sun_Velocity();
	double phi				= acos(x.Normalized() * vSun.Normalized());
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(event.Isoreflection_Angle(vSun), phi);
}

TEST(TestSimulationUtilities, TestIsoreflectionRing)
{
	// ARRANGE
	libphysica::Vector vSun = obscura::Sun_Velocity();
	unsigned int rings		= 3;
	double t				= 0;
	libphysica::Vector x({1, 2, 3});
	libphysica::Vector v({4, 5, 6});
	Event event(t, x, v);

	// ACT & ASSERT
	event.position = vSun;
	EXPECT_EQ(event.Isoreflection_Ring(vSun, rings), 0);
	event.position = -1.0 * vSun;
	EXPECT_EQ(event.Isoreflection_Ring(vSun, rings), 2);
	event.position = libphysica::Vector({vSun[1], -vSun[0], 0});
	EXPECT_EQ(event.Isoreflection_Ring(vSun, rings), 1);
}

TEST(TestSimulationUtilities, TestInUnits)
{
	// ARRANGE
	double t = 7 * sec;
	libphysica::Vector x({1 * km, 2 * km, 3 * km});
	libphysica::Vector v({4 * km / sec, 5 * km / sec, 6 * km / sec});
	Event event(t, x, v);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(event.In_Units(km, minute).time, 7.0 / 60);
	EXPECT_DOUBLE_EQ(event.In_Units(rSun, sec).position[0], km / rSun);
	EXPECT_DOUBLE_EQ(event.In_Units(km, sec).position[1], 2);
	EXPECT_DOUBLE_EQ(event.In_Units(km, sec).position[2], 3);
	EXPECT_DOUBLE_EQ(event.In_Units(km, sec).velocity[0], 4);
	EXPECT_DOUBLE_EQ(event.In_Units(km, sec).velocity[1], 5);
	EXPECT_DOUBLE_EQ(event.In_Units(km, sec).velocity[2], 6);
}

// 2. Generator of initial conditions

TEST(TestSimulationUtilities, TestInitialConditions)
{
	// ARRANGE
	std::random_device rd;
	std::mt19937 PRNG;
	PRNG.seed(3);

	Solar_Model SSM;

	obscura::Standard_Halo_Model SHM;
	// SHM.Set_Observer_Velocity(libphysica::Vector({0, 0, 0}));

	double R_distance	= 1000 * AU;
	unsigned int trials = 1000;
	// ACT & ASSERT
	for(unsigned int i = 0; i < trials; i++)
	{
		Event IC	= Initial_Conditions(SHM, SSM, PRNG);
		double r	= IC.Radius();
		double v	= IC.Speed();
		double vesc = SSM.Local_Escape_Speed(r);
		double Jmax = rSun * sqrt(v * v + SSM.Local_Escape_Speed(rSun) * SSM.Local_Escape_Speed(rSun));
		ASSERT_GE(IC.Radius(), R_distance);
		ASSERT_GE(IC.Speed(), vesc);
		ASSERT_LE(IC.Angular_Momentum(), Jmax);
		ASSERT_GE(IC.Angular_Momentum(), 0.0);
	}
}

TEST(TestSimulationUtilities, TestInitialSpeedPDFIsFiniteAtZero)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;
	EXPECT_DOUBLE_EQ(PDF_Initial_Speed(0.0, SHM, SSM), 0.0);
	EXPECT_TRUE(std::isfinite(PDF_Initial_Speed(1.0e-12 * km / sec, SHM, SSM)));
	EXPECT_THROW(PDF_Initial_Speed(-1.0, SHM, SSM), std::invalid_argument);
}

// 3. Analytically propagate a particle at event on a hyperbolic Kepler orbit to a radius R (without passing the periapsis)
TEST(TestSimulationUtilities, TestHyperbolicKeplerShift)
{
	int fixed_seed = 123;
	std::mt19937 PRNG(fixed_seed);

	Solar_Model SSM;

	obscura::Standard_Halo_Model SHM;
	SHM.Set_Observer_Velocity(libphysica::Vector({0, 0, 0}));

	int trials = 100;
	for(int k = 0; k < trials; k++)
	{
		Event IC = Initial_Conditions(SHM, SSM, PRNG);
		const double energy_initial = Kepler_Energy(IC);
		const double angular_momentum_initial = IC.Angular_Momentum();

		ASSERT_TRUE(Hyperbolic_Kepler_Shift(IC, 5.0 * rSun));
		ASSERT_NEAR(IC.Radius(), 5.0 * rSun, 1.0e-10 * rSun);
		EXPECT_LT(Radial_Velocity_For_Test(IC), 0.0);
		Expect_Relative_Near(Kepler_Energy(IC), energy_initial, 1.0e-10);
		Expect_Relative_Near(IC.Angular_Momentum(), angular_momentum_initial, 1.0e-10);

		ASSERT_TRUE(Hyperbolic_Kepler_Shift(IC, rSun));
		ASSERT_NEAR(IC.Radius(), rSun, 1.0e-10 * rSun);
		EXPECT_LT(Radial_Velocity_For_Test(IC), 0.0);
		Expect_Relative_Near(Kepler_Energy(IC), energy_initial, 1.0e-10);
		Expect_Relative_Near(IC.Angular_Momentum(), angular_momentum_initial, 1.0e-10);
		for(int i = 0; i < 3; i++)
		{
			EXPECT_TRUE(std::isfinite(IC.position[i]));
			EXPECT_TRUE(std::isfinite(IC.velocity[i]));
		}
	}
}

TEST(TestSimulationUtilities, TestHyperbolicKeplerShiftReturnsFalseForInvalidBranches)
{
	Event inside_sun(0.0, libphysica::Vector({0.5 * rSun, 0.0, 0.0}), libphysica::Vector({-1000.0 * km / sec, 0.0, 0.0}));
	EXPECT_FALSE(Hyperbolic_Kepler_Shift(inside_sun, 2.0 * rSun));

	Event inbound_outward_request(0.0, libphysica::Vector({5.0 * rSun, 0.0, 0.0}), libphysica::Vector({-1000.0 * km / sec, 0.0, 0.0}));
	Event inbound_outward_original = inbound_outward_request;
	EXPECT_FALSE(Hyperbolic_Kepler_Shift(inbound_outward_request, 6.0 * rSun));
	for(int i = 0; i < 3; i++)
	{
		EXPECT_DOUBLE_EQ(inbound_outward_request.position[i], inbound_outward_original.position[i]);
		EXPECT_DOUBLE_EQ(inbound_outward_request.velocity[i], inbound_outward_original.velocity[i]);
	}

	Event bound_orbit(0.0, libphysica::Vector({2.0 * rSun, 0.0, 0.0}), libphysica::Vector({-1.0 * km / sec, 0.0, 0.0}));
	EXPECT_FALSE(Hyperbolic_Kepler_Shift(bound_orbit, rSun));
}

// 4. Equiareal isodetection rings
TEST(TestSimulationUtilities, TestIsoreflectionRingAngles)
{
	// ARRANGE
	std::vector<double> two_rings	  = {90 * deg, 180 * deg};
	std::vector<double> hundred_rings = {36.8699 * deg, 53.1301 * deg, 66.4218 * deg, 78.463 * deg, 90. * deg, 101.537 * deg, 113.578 * deg, 126.87 * deg, 143.13 * deg, 180. * deg};

	double tol = 1.0e-3 * deg;
	// ACT & ASSERT
	for(unsigned int i = 0; i < two_rings.size(); i++)
		EXPECT_NEAR(Isoreflection_Ring_Angles(2)[i], two_rings[i], tol);
	for(unsigned int i = 0; i < hundred_rings.size(); i++)
		EXPECT_NEAR(Isoreflection_Ring_Angles(10)[i], hundred_rings[i], tol);
	EXPECT_THROW(Isoreflection_Ring_Angles(0), std::invalid_argument);
}

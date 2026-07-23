#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"

#include "Simulation_Trajectory.hpp"
#include "Simulation_Utilities.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

namespace
{
struct AggregatedBincount
{
	std::array<double, NUM_BINS> dt{};
	std::array<double, NUM_BINS> v2dt{};
};

void AddContributions(
	AggregatedBincount& aggregate,
	const std::vector<BincountContribution>& contributions)
{
	for(const BincountContribution& contribution : contributions)
	{
		ASSERT_GE(contribution.bin, 0);
		ASSERT_LT(contribution.bin, NUM_BINS);
		aggregate.dt[contribution.bin] += contribution.dt_sec;
		aggregate.v2dt[contribution.bin] += contribution.v2dt_km2_per_sec;
	}
}

AggregatedBincount ComputeContributions(const Event& before, const Event& after)
{
	std::vector<BincountContribution> contributions;
	Compute_Bincount_Interval_Contributions(before, after, contributions);
	AggregatedBincount aggregate;
	AddContributions(aggregate, contributions);
	return aggregate;
}

double SumBins(const std::array<double, NUM_BINS>& values)
{
	double sum = 0.0;
	for(double value : values)
		sum += value;
	return sum;
}
}

// 1. Result of one trajectory
TEST(TestSimulationTrajectory, TestTrajectoryResultConstructor)
{
	// ARRANGE
	double t_i = 0;
	libphysica::Vector r_i({1, 2, 3});
	libphysica::Vector v_i({4, 5, 6});
	double t_f = 7;
	libphysica::Vector r_f({8, 9, 10});
	libphysica::Vector v_f({11, 12, 13});
	Event initial_event(t_i, r_i, v_i);
	Event final_event(t_f, r_f, v_f);
	unsigned int nscat = 14;
	// ACT
	Trajectory_Result result(initial_event, final_event, nscat, TrajectoryBincount{});
	// ASSERT
	EXPECT_EQ(result.initial_event.time, t_i);
	EXPECT_EQ(result.initial_event.position, r_i);
	EXPECT_EQ(result.initial_event.velocity, v_i);
	EXPECT_EQ(result.final_event.time, t_f);
	EXPECT_EQ(result.final_event.position, r_f);
	EXPECT_EQ(result.final_event.velocity, v_f);
	EXPECT_EQ(result.number_of_scatterings, nscat);
}

TEST(TestSimulationTrajectory, TestParticleReflected)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector r_in({0, 0.5 * rSun, 0});
	libphysica::Vector r_out({0, 1.5 * rSun, 0});
	libphysica::Vector v_1({0, 0, 0});
	libphysica::Vector v_2({0, 1000 * km / sec, 0});
	Event initial_event;
	Event not_reflected_1(t, r_in, v_1);
	Event not_reflected_2(t, r_in, v_2);
	Event not_reflected_3(t, r_out, v_1);
	Event reflected_1(t, r_out, v_2);
	TrajectoryBincount escaped_bincount;
	escaped_bincount.termination_reason = TrajectoryTerminationReason::OutwardEscape;

	// ACT & ASSERT
	EXPECT_FALSE(Trajectory_Result(initial_event, not_reflected_1, 1, TrajectoryBincount{}).Particle_Reflected());
	EXPECT_FALSE(Trajectory_Result(initial_event, not_reflected_2, 1, TrajectoryBincount{}).Particle_Reflected());
	EXPECT_FALSE(Trajectory_Result(initial_event, not_reflected_3, 1, TrajectoryBincount{}).Particle_Reflected());
	EXPECT_FALSE(Trajectory_Result(initial_event, reflected_1, 1, TrajectoryBincount{}).Particle_Reflected());
	EXPECT_TRUE(Trajectory_Result(initial_event, reflected_1, 1, escaped_bincount).Particle_Reflected());
	EXPECT_FALSE(Trajectory_Result(initial_event, reflected_1, 0, escaped_bincount).Particle_Reflected());
}

TEST(TestSimulationTrajectory, TestParticleFree)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector r({0, 0.5 * rSun, 0});
	libphysica::Vector v({0, 1000 * km / sec, 0});
	Event initial_event;
	Event final_event(t, r, v);
	TrajectoryBincount escaped_bincount;
	escaped_bincount.termination_reason = TrajectoryTerminationReason::OutwardEscape;

	// ACT & ASSERT
	EXPECT_FALSE(Trajectory_Result(initial_event, final_event, 1, TrajectoryBincount{}).Particle_Free());
	EXPECT_FALSE(Trajectory_Result(initial_event, final_event, 0, TrajectoryBincount{}).Particle_Free());
	EXPECT_TRUE(Trajectory_Result(initial_event, final_event, 0, escaped_bincount).Particle_Free());
}

TEST(TestSimulationTrajectory, TestParticleCaptured)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector r_in({0, 0.5 * rSun, 0});
	libphysica::Vector r_out({0, 1.5 * rSun, 0});
	libphysica::Vector v_1({0, 0, 0});
	libphysica::Vector v_2({0, 1000 * km / sec, 0});
	Event initial_event;
	Event captured_1(t, r_in, v_1);
	Event not_captured_2(t, r_in, v_2);
	Event captured_3(t, r_out, v_1);
	Event not_captured(t, r_out, v_2);
	Solar_Model solar_model;

	// ACT & ASSERT
	EXPECT_TRUE(Trajectory_Result(initial_event, captured_1, 0, TrajectoryBincount{}).Particle_Captured(solar_model));
	EXPECT_FALSE(Trajectory_Result(initial_event, not_captured_2, 0, TrajectoryBincount{}).Particle_Captured(solar_model));
	EXPECT_TRUE(Trajectory_Result(initial_event, captured_3, 0, TrajectoryBincount{}).Particle_Captured(solar_model));
	EXPECT_FALSE(Trajectory_Result(initial_event, not_captured, 0, TrajectoryBincount{}).Particle_Captured(solar_model));
}

// TEST(TestSimulationTrajectory, TestResultPrintSummary)
// {
// 	//ARRANGE
// 	//ACT
// 	//ASSERT
// }

// 2. Simulator

TEST(TestSimulationTrajectory, TestDefaultTrajectoryWallTimeIsUnlimited)
{
	SnapshotConfig cfg;
	EXPECT_FALSE(cfg.enabled);
	EXPECT_DOUBLE_EQ(cfg.max_trajectory_wall_time_sec, 0.0);

	Solar_Model SSM;
	Trajectory_Simulator simulator(SSM);
	EXPECT_DOUBLE_EQ(simulator.max_trajectory_wall_time_sec, 0.0);
}

TEST(TestSimulationTrajectory, TestSimulatorRejectsInvalidNumericalBounds)
{
	Solar_Model SSM;
	EXPECT_THROW(Trajectory_Simulator(SSM, 100, 100, 0.0), std::invalid_argument);
	EXPECT_THROW(Trajectory_Simulator(SSM, 100, 100, std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
}

TEST(TestSimulationTrajectory, TestTrajectoryBincountSurvivalDefaults)
{
	TrajectoryBincount bincount;
	EXPECT_FALSE(bincount.is_captured);
	EXPECT_FALSE(bincount.event_observed);
	EXPECT_FALSE(bincount.boundary_escape_observed);
	EXPECT_TRUE(bincount.survival_valid);
	EXPECT_FALSE(bincount.numerically_invalid_escape);
	EXPECT_FALSE(bincount.truncated);
	EXPECT_DOUBLE_EQ(bincount.t_first_negative, -1.0);
	EXPECT_DOUBLE_EQ(bincount.t_last_negative, -1.0);
	EXPECT_DOUBLE_EQ(bincount.t_capture, -1.0);
	EXPECT_DOUBLE_EQ(bincount.t_last_bound, -1.0);
	EXPECT_DOUBLE_EQ(bincount.t_termination, -1.0);
	EXPECT_TRUE(std::isnan(bincount.t_final_unbinding_scatter));
	EXPECT_TRUE(std::isnan(bincount.t_first_unbinding_scatter));
	EXPECT_TRUE(std::isnan(bincount.t_boundary_escape));
	EXPECT_DOUBLE_EQ(bincount.max_free_energy_drift_eV, 0.0);
	EXPECT_DOUBLE_EQ(bincount.max_free_energy_drift_rel, 0.0);
	EXPECT_EQ(bincount.number_of_scatterings, 0UL);
	EXPECT_EQ(bincount.number_of_bound_to_unbound, 0UL);
	EXPECT_EQ(bincount.number_of_recaptures, 0UL);
	EXPECT_EQ(bincount.number_of_integrator_steps_after_capture, 0UL);
	EXPECT_TRUE(std::isnan(bincount.min_energy_after_capture_eV));
	EXPECT_EQ(TRAJECTORY_TERMINATION_REASON_COUNT, 11);
	EXPECT_EQ(static_cast<int>(TrajectoryTerminationReason::EnergyDriftEscape), 10);
}

TEST(TestSimulationTrajectory, TestBincountDepositsLinearRadialIntervalAcrossEveryBin)
{
	const int first_bin = 1500;
	const double start_radius_km = (static_cast<double>(first_bin) + 0.25) * BIN_WIDTH_KM;
	const double speed_km_s = BIN_WIDTH_KM;
	const double duration_sec = 10.5;
	Event before(
	    0.0,
	    libphysica::Vector({start_radius_km * km, 0.0, 0.0}),
	    libphysica::Vector({speed_km_s * km / sec, 0.0, 0.0}));
	Event after(
	    duration_sec * sec,
	    libphysica::Vector({
	        (start_radius_km + speed_km_s * duration_sec) * km, 0.0, 0.0}),
	    before.velocity);

	const AggregatedBincount aggregate = ComputeContributions(before, after);
	EXPECT_NEAR(aggregate.dt[first_bin], 0.75, 1.0e-11);
	for(int bin = first_bin + 1; bin < first_bin + 10; bin++)
		EXPECT_NEAR(aggregate.dt[bin], 1.0, 1.0e-11);
	EXPECT_NEAR(aggregate.dt[first_bin + 10], 0.75, 1.0e-11);
	EXPECT_NEAR(SumBins(aggregate.dt), duration_sec, 1.0e-10);
	EXPECT_NEAR(
	    SumBins(aggregate.v2dt),
	    speed_km_s * speed_km_s * duration_sec,
	    1.0e-9 * speed_km_s * speed_km_s * duration_sec);
}

TEST(TestSimulationTrajectory, TestBincountLinearDepositionIsStepPhaseInvariant)
{
	const double start_radius_km = 1300.37 * BIN_WIDTH_KM;
	const double speed_km_s = 0.8 * BIN_WIDTH_KM;
	const double duration_sec = 17.25;
	const double split_sec = 4.137;
	auto event_at = [&](double time_sec)
	{
		return Event(
		    time_sec * sec,
		    libphysica::Vector({
		        (start_radius_km + speed_km_s * time_sec) * km, 0.0, 0.0}),
		    libphysica::Vector({speed_km_s * km / sec, 0.0, 0.0}));
	};

	const AggregatedBincount one_step =
	    ComputeContributions(event_at(0.0), event_at(duration_sec));
	AggregatedBincount split_steps;
	std::vector<BincountContribution> contributions;
	Compute_Bincount_Interval_Contributions(
	    event_at(0.0), event_at(split_sec), contributions);
	AddContributions(split_steps, contributions);
	Compute_Bincount_Interval_Contributions(
	    event_at(split_sec), event_at(duration_sec), contributions);
	AddContributions(split_steps, contributions);

	for(int bin = 0; bin < NUM_BINS; bin++)
	{
		EXPECT_NEAR(one_step.dt[bin], split_steps.dt[bin], 1.0e-10);
		EXPECT_NEAR(
		    one_step.v2dt[bin],
		    split_steps.v2dt[bin],
		    1.0e-8 * std::max(1.0, one_step.v2dt[bin]));
	}
}

TEST(TestSimulationTrajectory, TestBincountClipsOutwardIntervalAtTwoSolarRadii)
{
	const double speed_km_s = BIN_WIDTH_KM;
	const double start_radius_km = 1998.5 * BIN_WIDTH_KM;
	const double duration_sec = 3.0;
	Event before(
	    0.0,
	    libphysica::Vector({start_radius_km * km, 0.0, 0.0}),
	    libphysica::Vector({speed_km_s * km / sec, 0.0, 0.0}));
	Event after(
	    duration_sec * sec,
	    libphysica::Vector({
	        (start_radius_km + speed_km_s * duration_sec) * km, 0.0, 0.0}),
	    before.velocity);

	const AggregatedBincount aggregate = ComputeContributions(before, after);
	EXPECT_NEAR(aggregate.dt[1998], 0.5, 1.0e-11);
	EXPECT_NEAR(aggregate.dt[1999], 1.0, 1.0e-11);
	EXPECT_NEAR(SumBins(aggregate.dt), 1.5, 1.0e-10);
}

TEST(TestSimulationTrajectory, TestSurvivalInvalidTerminationReasons)
{
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::EnergyDriftEscape));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::NumericalFailure));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::NonFiniteState));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::SpeedLimit));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::Unknown));

	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::WallTimeLimit));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::MaxFreeSteps));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::MaxScatterings));
	EXPECT_FALSE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::OutwardEscape));
}

TEST(TestSimulationTrajectory, TestScatter)
{
	// ARRANGE
	obscura::DM_Particle_SI DM(0.5 * GeV);
	DM.Set_Sigma_Proton(pb);
	double time = 0.0;
	libphysica::Vector r({0.25 * rSun, 0.25 * rSun, 0.25 * rSun});
	libphysica::Vector vel({100 * km / sec, 100 * km / sec, 100 * km / sec});
	Event event(time, r, vel);
	Solar_Model SSM;
	Trajectory_Simulator simulator(SSM);
	// ACT & ASSERT
	int trials = 10;
	for(int i = 0; i < trials; i++)
	{
		libphysica::Vector v_ini = event.velocity;
		simulator.Scatter(event, DM);
		for(int j = 0; j < 3; j++)
			EXPECT_NE(event.velocity[j], v_ini[j]);
	}
}

TEST(TestSimulationTrajectory, TestScatterWithZAxisVelocityIsFinite)
{
	obscura::DM_Particle_SI DM(0.5 * GeV);
	DM.Set_Sigma_Proton(pb);
	DM.Set_Sigma_Electron(pb);

	Solar_Model SSM;
	Trajectory_Simulator simulator(SSM);
	std::vector<libphysica::Vector> velocities;
	velocities.push_back(libphysica::Vector({0.0, 0.0, 100.0 * km / sec}));
	velocities.push_back(libphysica::Vector({0.0, 0.0, -100.0 * km / sec}));
	velocities.push_back(libphysica::Vector({1.0e-15 * km / sec, 0.0, 100.0 * km / sec}));

	for(const auto& velocity : velocities)
	{
		Event event(0.0, libphysica::Vector({0.25 * rSun, 0.0, 0.0}), velocity);
		ASSERT_NO_THROW(simulator.Scatter(event, DM));
		for(int component = 0; component < 3; component++)
			EXPECT_TRUE(std::isfinite(event.velocity[component]));
	}
}

TEST(TestSimulationTrajectory, TestMaxFreeStepsTerminationReason)
{
	obscura::DM_Particle_SI DM(0.5 * GeV);
	DM.Set_Sigma_Proton(0.1 * pb);
	Solar_Model SSM;
	Trajectory_Simulator simulator(SSM, 0, 10, 2.0 * rSun);

	const double radius = 0.5 * rSun;
	const double unbound_speed = 1.1 * SSM.Local_Escape_Speed(radius);
	Event IC(0.0, libphysica::Vector({radius, 0.0, 0.0}), libphysica::Vector({0.0, unbound_speed, 0.0}));
	Trajectory_Result result = simulator.Simulate(IC, DM, 0);

	EXPECT_EQ(result.bincount.termination_reason, TrajectoryTerminationReason::MaxFreeSteps);
	EXPECT_FALSE(result.bincount.truncated);
	EXPECT_FALSE(result.bincount.survival_valid);
	EXPECT_FALSE(result.Particle_Free());
	EXPECT_FALSE(result.Particle_Reflected());
}

TEST(TestSimulationTrajectory, TestUncapturedBoundFreeFlightTerminatesAsNumericalFailure)
{
	obscura::DM_Particle_SI DM(0.5 * GeV);
	DM.Set_Sigma_Proton(1.0e-40 * cm * cm);
	Solar_Model SSM;
	Trajectory_Simulator simulator(SSM, 100, 10, 2.0 * rSun);

	const double radius = 1.5 * rSun;
	const double bound_speed = 0.5 * SSM.Local_Escape_Speed(radius);
	Event IC(
		0.0,
		libphysica::Vector({radius, 0.0, 0.0}),
		libphysica::Vector({bound_speed, 0.0, 0.0}));

	testing::internal::CaptureStderr();
	Trajectory_Result result = simulator.Simulate(IC, DM, 0);
	const std::string stderr_output = testing::internal::GetCapturedStderr();

	EXPECT_EQ(result.bincount.termination_reason, TrajectoryTerminationReason::NumericalFailure);
	EXPECT_TRUE(stderr_output.empty());
	EXPECT_FALSE(result.bincount.is_captured);
	EXPECT_FALSE(result.bincount.survival_valid);
	EXPECT_EQ(result.number_of_scatterings, 0UL);
	EXPECT_DOUBLE_EQ(result.bincount.t_termination, 0.0);
}

TEST(TestSimulationTrajectory, TestSimulate)
{
	// ARRANGE
	obscura::DM_Particle_SI DM(0.5 * GeV);
	DM.Set_Sigma_Proton(0.1 * pb);
	Solar_Model SSM;
	Trajectory_Simulator simulator(SSM);
	obscura::Standard_Halo_Model SHM;
	// ACT & ASSERT
	int trials = 5;
	for(int i = 0; i < trials; i++)
	{
		Event IC = Initial_Conditions(SHM, SSM, simulator.PRNG);
		ASSERT_TRUE(Hyperbolic_Kepler_Shift(IC, 1.5 * rSun));
		Trajectory_Result result = simulator.Simulate(IC, DM, 0);
		if(result.Particle_Reflected() || result.Particle_Free())
			ASSERT_NEAR(result.final_event.Radius(), simulator.maximum_distance, 1.0e-10 * rSun);
		else
			ASSERT_GT(result.number_of_scatterings, 0);
	}
}

TEST(TestSimulationTrajectory, TestSerializedPRNGStateReplaysTrajectoryExactly)
{
	obscura::DM_Particle_SI DM(0.5 * GeV);
	DM.Set_Sigma_Proton(0.1 * pb);
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	Trajectory_Simulator original(SSM, 1000000, 1000, 2.0 * rSun);
	original.Fix_PRNG_Seed(20260722);
	Event IC = Initial_Conditions(SHM, SSM, original.PRNG);
	ASSERT_TRUE(Hyperbolic_Kepler_Shift(IC, 1.5 * rSun));
	const std::string state_before_simulation = original.Serialize_PRNG_State();
	original.Enable_Diagnostic_Trace(true);
	Trajectory_Result first = original.Simulate(IC, DM, 0);

	Trajectory_Simulator replay(SSM, 1000000, 1000, 2.0 * rSun);
	ASSERT_NO_THROW(replay.Restore_PRNG_State(state_before_simulation));
	replay.Enable_Diagnostic_Trace(true);
	Trajectory_Result second = replay.Simulate(IC, DM, 0);

	EXPECT_EQ(first.number_of_scatterings, second.number_of_scatterings);
	EXPECT_EQ(first.bincount.termination_reason, second.bincount.termination_reason);
	EXPECT_DOUBLE_EQ(first.bincount.t_capture, second.bincount.t_capture);
	if(std::isnan(first.bincount.t_final_unbinding_scatter))
		EXPECT_TRUE(std::isnan(second.bincount.t_final_unbinding_scatter));
	else
		EXPECT_DOUBLE_EQ(first.bincount.t_final_unbinding_scatter, second.bincount.t_final_unbinding_scatter);
	EXPECT_DOUBLE_EQ(first.final_event.time, second.final_event.time);
	for(size_t component = 0; component < 3; component++)
	{
		EXPECT_DOUBLE_EQ(first.final_event.position[component], second.final_event.position[component]);
		EXPECT_DOUBLE_EQ(first.final_event.velocity[component], second.final_event.velocity[component]);
	}
	ASSERT_EQ(first.diagnostic_events.size(), second.diagnostic_events.size());
	for(size_t index = 0; index < first.diagnostic_events.size(); index++)
	{
		EXPECT_EQ(first.diagnostic_events[index].event_type, second.diagnostic_events[index].event_type);
		EXPECT_DOUBLE_EQ(first.diagnostic_events[index].t_s, second.diagnostic_events[index].t_s);
		EXPECT_EQ(first.diagnostic_events[index].scatter_index, second.diagnostic_events[index].scatter_index);
		EXPECT_EQ(first.diagnostic_events[index].step_index, second.diagnostic_events[index].step_index);
	}
}

TEST(TestSimulationTrajectory, TestSimulatorPrintSummary)
{
	// ARRANGE
	obscura::DM_Particle_SI DM(0.5 * GeV);
	DM.Set_Sigma_Proton(0.1 * pb);
	Solar_Model SSM;
	Trajectory_Simulator simulator(SSM);
	obscura::Standard_Halo_Model SHM;
	// ACT
	Event IC = Initial_Conditions(SHM, SSM, simulator.PRNG);
	ASSERT_TRUE(Hyperbolic_Kepler_Shift(IC, 1.5 * rSun));
	Trajectory_Result result = simulator.Simulate(IC, DM, 0);
	// ASSERT
	result.Print_Summary(SSM);
}

// 3. Equation of motion solution with Runge-Kutta-Fehlberg

TEST(TestSimulationTrajectory, TestRungeKuttaStep)
{
	// ASSERT
	double t = 0;
	libphysica::Vector r({0.25 * rSun, 0.5 * rSun, 0.25 * rSun});
	libphysica::Vector v({km / sec, 1000 * km / sec, km / sec});
	Event event(t, r, v);
	Free_Particle_Propagator propagator(event);
	double initial_timestep = propagator.time_step;
	Solar_Model SSM;
	// ACT
	propagator.Runge_Kutta_45_Step(SSM);
	// // ASSERT
	EXPECT_NE(propagator.time_step, initial_timestep);
	EXPECT_GT(propagator.Current_Time(), event.time);
	EXPECT_GT(propagator.Current_Radius(), event.Radius());
	EXPECT_LT(propagator.Current_Speed(), event.Speed());
}

TEST(TestSimulationTrajectory, TestRungeKuttaZeroErrorStepGrowthIsBounded)
{
	// ASSERT
	double t = 0.0;
	libphysica::Vector r({17.0 * km, 0.0, 0.0});
	libphysica::Vector v({0.0, km / sec, 0.0});
	Event event(t, r, v);
	Free_Particle_Propagator propagator(event);
	double initial_timestep = propagator.time_step;

	// ACT
	for(int i = 0; i < 3; i++)
		propagator.Runge_Kutta_45_Step(0.0);

	// ASSERT
	EXPECT_TRUE(std::isfinite(propagator.Current_Time()));
	EXPECT_TRUE(std::isfinite(propagator.time_step));
	EXPECT_LT(propagator.Current_Time(), 100.0 * initial_timestep);
	EXPECT_LE(propagator.time_step, 64.0 * initial_timestep);
	EXPECT_GT(propagator.Current_Radius(), 0.0);
	EXPECT_GT(propagator.Current_Speed(), 0.0);
}

// Regression guard: before the non-finite-error fix, a RK stage that produced
// NaN errors caused RK45_Next_Step_Size to grow the step geometrically and the
// inner while(!accepted) loop spun forever, hanging one MPI rank and deadlocking
// the snapshot / MPI_Barrier path (observed on 0.5 GeV runs at the 300 s snapshot).
// This test seeds the propagator with a NaN mass so every derivative becomes NaN,
// then checks that the step returns in bounded time and never pollutes public state
// with non-finite values.
TEST(TestSimulationTrajectory, TestRungeKuttaHandlesNonFiniteDerivativesWithoutHanging)
{
	double t = 0.0;
	libphysica::Vector r({17.0 * km, 0.0, 0.0});
	libphysica::Vector v({0.0, km / sec, 0.0});
	Event event(t, r, v);
	Free_Particle_Propagator propagator(event);

	// ACT: a NaN mass makes every dv/dt NaN, hence every stage and every error NaN.
	const double nan_mass = std::nan("");
	bool step_ok = true;
	for(int i = 0; i < 5; i++)
		step_ok = propagator.Runge_Kutta_45_Step(nan_mass);

	// ASSERT: propagator must still expose finite public state and not have
	// exploded the time step to infinity.
	EXPECT_FALSE(step_ok);
	EXPECT_TRUE(std::isfinite(propagator.Current_Time()));
	EXPECT_TRUE(std::isfinite(propagator.time_step));
	EXPECT_TRUE(std::isfinite(propagator.Current_Radius()));
	EXPECT_TRUE(std::isfinite(propagator.Current_Speed()));
	EXPECT_GE(propagator.Current_Radius(), 0.0);
}

TEST(TestSimulationTrajectory, TestPropagatorCurrentRadius)
{
	// ASSERT
	double t	  = 13 * sec;
	double radius = 17 * km;
	libphysica::Vector r({radius, 0, 0});
	libphysica::Vector vel;
	Event event(t, r, vel);
	Free_Particle_Propagator propagator(event);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(propagator.Current_Radius(), radius);
}

TEST(TestSimulationTrajectory, TestPropagatorCurrentSpeed)
{
	// ASSERT
	double t = 13 * sec;
	double v = 123.0 * km / sec;
	libphysica::Vector r;
	libphysica::Vector vel({0, v, 0});
	Event event(t, r, vel);
	Free_Particle_Propagator propagator(event);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(propagator.Current_Speed(), v);
}

TEST(TestSimulationTrajectory, TestPropagatorEventIn3D)
{
	// ASSERT
	double t = 13 * sec;
	libphysica::Vector r({km, 0.5 * rSun, km});
	libphysica::Vector v({km / sec, 1000 * km / sec, km / sec});
	Event event(t, r, v);
	Free_Particle_Propagator propagator(event);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(propagator.Event_In_3D().time, t);
	for(int i = 0; i < 3; i++)
		EXPECT_DOUBLE_EQ(propagator.Event_In_3D().position[i], r[i]);
	for(int i = 0; i < 3; i++)
		EXPECT_DOUBLE_EQ(propagator.Event_In_3D().velocity[i], v[i]);
}

TEST(TestSimulationTrajectory, TestPropagatorEventIn3DRadialTrajectory)
{
	double t = 13 * sec;
	libphysica::Vector r({2.0 * rSun, 0.0, 0.0});
	libphysica::Vector v({-100.0 * km / sec, 0.0, 0.0});
	Event event(t, r, v);
	Free_Particle_Propagator propagator(event);

	Event recovered = propagator.Event_In_3D();

	EXPECT_DOUBLE_EQ(recovered.time, t);
	for(int i = 0; i < 3; i++)
	{
		EXPECT_TRUE(std::isfinite(recovered.position[i]));
		EXPECT_TRUE(std::isfinite(recovered.velocity[i]));
		EXPECT_DOUBLE_EQ(recovered.position[i], r[i]);
		EXPECT_DOUBLE_EQ(recovered.velocity[i], v[i]);
	}
}

TEST(TestSimulationTrajectory, TestPropagatorEventIn3DZeroRadius)
{
	double t = 13 * sec;
	double v0 = 123.0 * km / sec;
	libphysica::Vector r({0.0, 0.0, 0.0});
	libphysica::Vector v({0.0, v0, 0.0});
	Event event(t, r, v);
	Free_Particle_Propagator propagator(event);

	Event recovered = propagator.Event_In_3D();

	EXPECT_DOUBLE_EQ(recovered.time, t);
	for(int i = 0; i < 3; i++)
	{
		EXPECT_TRUE(std::isfinite(recovered.position[i]));
		EXPECT_TRUE(std::isfinite(recovered.velocity[i]));
		EXPECT_DOUBLE_EQ(recovered.position[i], r[i]);
		EXPECT_DOUBLE_EQ(recovered.velocity[i], v[i]);
	}
}

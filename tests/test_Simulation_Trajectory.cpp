#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"

#include "Simulation_Trajectory.hpp"
#include "Solar_Model.hpp"
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


TEST(TestSimulationTrajectory, RuntimeRadialGridUsesEarthScaleForBincount)
{
        constexpr double R_EARTH_KM = 6371.0;

        Bincount_Radial_Grid earth_grid(R_EARTH_KM);

        // 500.5 inner-bin widths corresponds to 0.5005 R_earth. It lies
        // strictly inside Earth bin 500, avoiding any bin-edge ambiguity.
        const double radius =
            500.5 * earth_grid.Inner_Bin_Width_Km() * km;

        const Event before(
            0.0,
            libphysica::Vector({radius, 0.0, 0.0}),
            libphysica::Vector({0.0, 0.0, 0.0}));

        const Event after(
            1.0 * sec,
            libphysica::Vector({radius, 0.0, 0.0}),
            libphysica::Vector({0.0, 0.0, 0.0}));

        std::vector<BincountContribution> earth_contributions;
        Compute_Bincount_Interval_Contributions(
            before,
            after,
            earth_grid,
            earth_contributions);

        ASSERT_EQ(1u, earth_contributions.size());
        EXPECT_EQ(500, earth_contributions.front().bin);
        EXPECT_DOUBLE_EQ(1.0, earth_contributions.front().dt_sec);
        EXPECT_DOUBLE_EQ(
            0.0,
            earth_contributions.front().v2dt_km2_per_sec);

        // The legacy no-grid overload remains Sun-scaled for compatibility.
        std::vector<BincountContribution> legacy_contributions;
        Compute_Bincount_Interval_Contributions(
            before,
            after,
            legacy_contributions);

        ASSERT_EQ(1u, legacy_contributions.size());
        EXPECT_EQ(4, legacy_contributions.front().bin);
        EXPECT_DOUBLE_EQ(1.0, legacy_contributions.front().dt_sec);
}

TEST(TestSimulationTrajectory, SnapshotProgressPublicationHasStepAndWallClockBounds)
{
	const double wall_interval = SnapshotProgressPublishWallIntervalSeconds();
	ASSERT_GT(wall_interval, 0.0);
	EXPECT_FALSE(SnapshotProgressPublishDue(0, 0.0, false));
	EXPECT_FALSE(SnapshotProgressPublishDue(511, wall_interval * 0.5, false));
	EXPECT_TRUE(SnapshotProgressPublishDue(512, 0.0, false));
	EXPECT_TRUE(SnapshotProgressPublishDue(1, wall_interval, false));
	EXPECT_TRUE(SnapshotProgressPublishDue(0, 0.0, true));
	EXPECT_FALSE(SnapshotProgressPublishDue(
	    0, std::numeric_limits<double>::quiet_NaN(), false));
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
	EXPECT_FALSE(bincount.outer_domain_removed);
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
	EXPECT_EQ(bincount.number_of_bound_exterior_arcs, 0UL);
	EXPECT_TRUE(std::isnan(bincount.first_bound_exit_kepler_period_sec));
	EXPECT_TRUE(std::isnan(bincount.last_bound_exit_kepler_period_sec));
	EXPECT_TRUE(std::isnan(bincount.max_bound_exit_kepler_period_sec));
	EXPECT_TRUE(std::isnan(bincount.first_bound_exit_exterior_time_sec));
	EXPECT_TRUE(std::isnan(bincount.last_bound_exit_exterior_time_sec));
	EXPECT_TRUE(std::isnan(bincount.max_bound_exit_exterior_time_sec));
	EXPECT_TRUE(std::isnan(bincount.min_energy_after_capture_eV));
	EXPECT_EQ(TRAJECTORY_TERMINATION_REASON_COUNT, 12);
	EXPECT_EQ(static_cast<int>(TrajectoryTerminationReason::EnergyDriftEscape), 10);
	EXPECT_EQ(static_cast<int>(TrajectoryTerminationReason::OuterDomainRemoval), 11);
	EXPECT_EQ(
	    bincount.numerical_failure_detail,
	    TrajectoryNumericalFailureDetail::None);
	EXPECT_TRUE(std::isnan(bincount.failure_energy_before_step_eV));
	EXPECT_TRUE(std::isnan(bincount.failure_energy_after_step_eV));
	EXPECT_TRUE(std::isnan(bincount.failure_energy_at_boundary_eV));
	EXPECT_TRUE(std::isnan(bincount.failure_reference_energy_eV));
	EXPECT_TRUE(std::isnan(bincount.failure_boundary_vr_km_s));
	EXPECT_EQ(
	    std::string(TrajectoryNumericalFailureDetailKey(
	        TrajectoryNumericalFailureDetail::
	            BoundaryEnergyMismatch)),
	    "boundary_energy_mismatch");
}

TEST(TestSimulationTrajectory, HermiteBoundaryLocatorSelectsOutwardCrossing)
{
	const double target_radius = 100.0 * km;
	Event before(
	    0.0 * sec,
	    libphysica::Vector({90.0 * km, 0.0, 0.0}),
	    libphysica::Vector({2.0 * km / sec, 0.0, 0.0}));
	Event after(
	    100.0 * sec,
	    libphysica::Vector({110.0 * km, 0.0, 0.0}),
	    libphysica::Vector({-1.0 * km / sec, 0.0, 0.0}));
	Event crossing;

	ASSERT_TRUE(Find_First_Outward_Hermite_Radius_Crossing(
	    before, after, target_radius, crossing));
	EXPECT_NEAR(
	    In_Units(crossing.Radius(), km),
	    In_Units(target_radius, km),
	    1.0e-10);
	EXPECT_GT(
	    In_Units(
	        crossing.position.Dot(crossing.velocity)
	        / crossing.Radius(),
	        km / sec),
	    0.0);
	EXPECT_GT(crossing.time, before.time);
	EXPECT_LT(crossing.time, after.time);
}

TEST(TestSimulationTrajectory, TestBincountDepositsLinearRadialIntervalAcrossEveryBin)
{
	const int first_bin = 750;
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

TEST(TestSimulationTrajectory, TestBincountDepositsSingleBinIntervalConservatively)
{
	const int bin = 100;
	const double duration_sec = 1.0;
	const double start_radius_km = (static_cast<double>(bin) + 0.2) * BIN_WIDTH_KM;
	const double speed_km_s = 0.6 * BIN_WIDTH_KM;
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
	EXPECT_NEAR(aggregate.dt[bin], duration_sec, 1.0e-12);
	EXPECT_NEAR(SumBins(aggregate.dt), duration_sec, 1.0e-12);
	EXPECT_NEAR(
	    SumBins(aggregate.v2dt),
	    speed_km_s * speed_km_s * duration_sec,
	    1.0e-12 * speed_km_s * speed_km_s * duration_sec);
}

TEST(TestSimulationTrajectory, TestBincountSingleBinFastPathRejectsCurvedBoundaryCrossing)
{
	const int endpoint_bin = 100;
	const double duration_sec = 1.0;
	const double endpoint_radius_km =
	    (static_cast<double>(endpoint_bin) + 0.5) * BIN_WIDTH_KM;
	const double endpoint_speed_km_s = 4.0 * BIN_WIDTH_KM;
	Event before(
	    0.0,
	    libphysica::Vector({endpoint_radius_km * km, 0.0, 0.0}),
	    libphysica::Vector({endpoint_speed_km_s * km / sec, 0.0, 0.0}));
	Event after(
	    duration_sec * sec,
	    before.position,
	    libphysica::Vector({-endpoint_speed_km_s * km / sec, 0.0, 0.0}));

	const AggregatedBincount aggregate = ComputeContributions(before, after);
	const double expected_outer_fraction = std::sqrt(0.5);
	EXPECT_NEAR(aggregate.dt[endpoint_bin], 1.0 - expected_outer_fraction, 1.0e-3);
	EXPECT_NEAR(aggregate.dt[endpoint_bin + 1], expected_outer_fraction, 1.0e-3);
	EXPECT_NEAR(SumBins(aggregate.dt), duration_sec, 1.0e-10);
	EXPECT_NEAR(
	    SumBins(aggregate.v2dt),
	    endpoint_speed_km_s * endpoint_speed_km_s / 3.0,
	    1.0e-9 * endpoint_speed_km_s * endpoint_speed_km_s);
}

TEST(TestSimulationTrajectory, TestBincountLinearDepositionIsStepPhaseInvariant)
{
	const double start_radius_km = 800.37 * BIN_WIDTH_KM;
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

TEST(TestSimulationTrajectory, TestNumericalBincountClipsAtKeplerBoundary)
{
	const double speed_km_s = BIN_WIDTH_KM;
	const double start_radius_km = 1098.5 * BIN_WIDTH_KM;
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
	EXPECT_NEAR(aggregate.dt[1098], 0.5, 1.0e-11);
	EXPECT_NEAR(aggregate.dt[1099], 1.0, 1.0e-11);
	EXPECT_NEAR(SumBins(aggregate.dt), 1.5, 1.0e-10);
}

TEST(TestSimulationTrajectory, CompactRadialGridReachesConfiguredCutoff)
{
	EXPECT_EQ(TOTAL_BINS, static_cast<std::size_t>(NUM_BINS) + EXTERIOR_BINS);
	EXPECT_DOUBLE_EQ(BincountBinLowerKm(0), 0.0);
	EXPECT_DOUBLE_EQ(BincountBinLowerKm(NUM_BINS), BIN_MAX_KM);
	EXPECT_NEAR(
	    BincountBinUpperKm(NUM_BINS) - BincountBinLowerKm(NUM_BINS),
	    BIN_WIDTH_KM,
	    1.0e-12 * BIN_WIDTH_KM);
	for(std::size_t exterior_bin = 0;
	    exterior_bin + 1 < EXTERIOR_FIRST_CAPPED_BIN
	    && exterior_bin + 2 < EXTERIOR_BINS;
	    exterior_bin++)
	{
		const std::size_t bin = static_cast<std::size_t>(NUM_BINS) + exterior_bin;
		const double width = BincountBinUpperKm(bin) - BincountBinLowerKm(bin);
		const double next_width =
		    BincountBinUpperKm(bin + 1) - BincountBinLowerKm(bin + 1);
		EXPECT_NEAR(
		    next_width / width,
		    EXTERIOR_BIN_GROWTH_FACTOR,
		    1.0e-11);
	}
	if(EXTERIOR_BINS > EXTERIOR_FIRST_CAPPED_BIN + 2)
	{
		const std::size_t first_capped_bin =
		    static_cast<std::size_t>(NUM_BINS) + EXTERIOR_FIRST_CAPPED_BIN;
		EXPECT_NEAR(
		    BincountBinUpperKm(first_capped_bin) - BincountBinLowerKm(first_capped_bin),
		    EXTERIOR_MAX_BIN_WIDTH_KM,
		    1.0e-12 * EXTERIOR_MAX_BIN_WIDTH_KM);
		EXPECT_NEAR(
		    BincountBinUpperKm(first_capped_bin + 1) - BincountBinLowerKm(first_capped_bin + 1),
		    EXTERIOR_MAX_BIN_WIDTH_KM,
		    1.0e-12 * EXTERIOR_MAX_BIN_WIDTH_KM);
	}
	EXPECT_NEAR(
	    BincountBinUpperKm(TOTAL_BINS - 1),
	    RADIAL_DOMAIN_MAX_KM,
	    1.0e-12 * RADIAL_DOMAIN_MAX_KM);
	EXPECT_EQ(BincountBinIndexKm(BIN_MAX_KM), NUM_BINS);
	EXPECT_EQ(BincountBinIndexKm(std::nextafter(RADIAL_DOMAIN_MAX_KM, 0.0)),
	          static_cast<int>(TOTAL_BINS - 1));
	EXPECT_EQ(BincountBinIndexKm(RADIAL_DOMAIN_MAX_KM), -1);
	for(std::size_t bin = 0; bin < TOTAL_BINS; bin++)
	{
		EXPECT_LT(BincountBinLowerKm(bin), BincountBinUpperKm(bin));
		EXPECT_EQ(BincountBinIndexKm(BincountBinLowerKm(bin)), static_cast<int>(bin));
	}
	EXPECT_LE(
	    BincountBinUpperKm(TOTAL_BINS - 1) - BincountBinLowerKm(TOTAL_BINS - 1),
	    EXTERIOR_MAX_BIN_WIDTH_KM);
}

TEST(TestSimulationTrajectory, BoundKeplerExteriorArcUsesGeometricWidthGrid)
{
	const double semi_major_axis = 1.5 * rSun;
	const double eccentricity = 1.0 / 3.0;
	const double boundary_radius = TRAJECTORY_BOUNDARY_RSUN * rSun;
	const double mu = G_Newton * mSun;
	const double angular_momentum =
	    std::sqrt(mu * semi_major_axis * (1.0 - eccentricity * eccentricity));
	const double tangential_speed = angular_momentum / boundary_radius;
	const double total_speed =
	    std::sqrt(mu * (2.0 / boundary_radius - 1.0 / semi_major_axis));
	const double radial_speed =
	    std::sqrt(total_speed * total_speed - tangential_speed * tangential_speed);
	Event outward(
	    3.0 * sec,
	    libphysica::Vector({boundary_radius, 0.0, 0.0}),
	    libphysica::Vector({radial_speed, tangential_speed, 0.0}));

	BoundKeplerExteriorArc arc;
	ASSERT_TRUE(Compute_Bound_Kepler_Exterior_Arc(outward, arc));
	EXPECT_FALSE(arc.outer_domain_removed);
	EXPECT_NEAR(arc.apoapsis_km, 2.0 * R_SUN_KM, 1.0e-10 * R_SUN_KM);
	EXPECT_EQ(arc.dt_hist.size(), TOTAL_BINS);
	EXPECT_EQ(arc.v2dt_hist.size(), arc.dt_hist.size());
	EXPECT_NEAR(
	    In_Units(arc.terminal_event.Radius(), rSun),
	    TRAJECTORY_BOUNDARY_RSUN,
	    1.0e-12);
	EXPECT_LT(
	    arc.terminal_event.position.Dot(arc.terminal_event.velocity),
	    0.0);
	EXPECT_NEAR(
	    In_Units(arc.terminal_event.Speed(), km / sec),
	    In_Units(outward.Speed(), km / sec),
	    1.0e-10);
	EXPECT_NEAR(
	    std::accumulate(arc.dt_hist.begin(), arc.dt_hist.end(), 0.0),
	    arc.elapsed_time_sec,
	    1.0e-10 * arc.elapsed_time_sec);
	const double expected_period_sec = In_Units(
	    2.0 * M_PI * std::sqrt(
	        semi_major_axis * semi_major_axis * semi_major_axis / mu),
	    sec);
	EXPECT_NEAR(
	    arc.kepler_period_sec,
	    expected_period_sec,
	    1.0e-12 * expected_period_sec);
	EXPECT_GT(arc.kepler_period_sec, arc.elapsed_time_sec);
	EXPECT_EQ(
	    std::accumulate(arc.dt_hist.begin(), arc.dt_hist.begin() + NUM_BINS, 0.0),
	    0.0);
	EXPECT_GT(
	    std::accumulate(arc.dt_hist.begin() + NUM_BINS, arc.dt_hist.end(), 0.0),
	    0.0);
	EXPECT_GT(
	    std::accumulate(arc.v2dt_hist.begin(), arc.v2dt_hist.end(), 0.0),
	    0.0);
}
TEST(TestSimulationTrajectory, BoundKeplerExteriorArcUsesRuntimeEarthRadialGrid)
{
        // Reuse a known-valid bound solar Kepler orbit. Only the histogram
        // grid changes here: the runtime Earth grid must determine the
        // shell layout and vector size.
        const double semi_major_axis = 1.5 * rSun;
        const double eccentricity = 1.0 / 3.0;
        const double boundary_radius = TRAJECTORY_BOUNDARY_RSUN * rSun;
        const double mu = G_Newton * mSun;
        const double angular_momentum =
            std::sqrt(mu * semi_major_axis * (1.0 - eccentricity * eccentricity));
        const double tangential_speed = angular_momentum / boundary_radius;
        const double total_speed =
            std::sqrt(mu * (2.0 / boundary_radius - 1.0 / semi_major_axis));
        const double radial_speed =
            std::sqrt(total_speed * total_speed - tangential_speed * tangential_speed);

        Event outward(
            3.0 * sec,
            libphysica::Vector({boundary_radius, 0.0, 0.0}),
            libphysica::Vector({radial_speed, tangential_speed, 0.0}));

		constexpr double earth_radius_km = 6371.0;

        const Bincount_Radial_Grid earth_grid(
            earth_radius_km,
            RADIAL_DOMAIN_MAX_AU);

        BoundKeplerExteriorArc earth_arc;
        ASSERT_TRUE(Compute_Bound_Kepler_Exterior_Arc(
            outward,
            earth_grid,
            earth_arc));

        EXPECT_FALSE(earth_arc.outer_domain_removed);
        EXPECT_EQ(earth_arc.dt_hist.size(), earth_grid.Bin_Count());
        EXPECT_EQ(earth_arc.v2dt_hist.size(), earth_grid.Bin_Count());
        EXPECT_NE(earth_arc.dt_hist.size(), TOTAL_BINS);

        const double boundary_radius_km =
            In_Units(boundary_radius, km);
        const int first_bin_index =
            earth_grid.Bin_Index_Km(boundary_radius_km);

        ASSERT_GE(
            first_bin_index,
            static_cast<int>(earth_grid.Inner_Bin_Count()));
        ASSERT_LT(
            static_cast<std::size_t>(first_bin_index),
            earth_grid.Bin_Count());

        // This distinguishes the Earth-scaled grid from the legacy Sun grid.
        EXPECT_NE(
            first_bin_index,
            BincountBinIndexKm(boundary_radius_km));

        const std::size_t first_bin =
            static_cast<std::size_t>(first_bin_index);

        EXPECT_EQ(
            std::accumulate(
                earth_arc.dt_hist.begin(),
                earth_arc.dt_hist.begin() + first_bin,
                0.0),
            0.0);

        EXPECT_GT(earth_arc.dt_hist[first_bin], 0.0);
        EXPECT_GT(
            std::accumulate(
                earth_arc.dt_hist.begin(),
                earth_arc.dt_hist.end(),
                0.0),
            0.0);

        EXPECT_NEAR(
            std::accumulate(
                earth_arc.dt_hist.begin(),
                earth_arc.dt_hist.end(),
                0.0),
            earth_arc.elapsed_time_sec,
            1.0e-10 * earth_arc.elapsed_time_sec);

        EXPECT_GT(
            std::accumulate(
                earth_arc.v2dt_hist.begin(),
                earth_arc.v2dt_hist.end(),
                0.0),
            0.0);
}

TEST(TestSimulationTrajectory, BoundKeplerExteriorArcFlagsApoapsisBeyondConfiguredCutoff)
{
	const double periapsis = 0.5 * rSun;
	const double apoapsis = 6.0 * AU;
	const double semi_major_axis = 0.5 * (periapsis + apoapsis);
	const double eccentricity = (apoapsis - periapsis) / (apoapsis + periapsis);
	const double boundary_radius = TRAJECTORY_BOUNDARY_RSUN * rSun;
	const double mu = G_Newton * mSun;
	const double angular_momentum =
	    std::sqrt(mu * semi_major_axis * (1.0 - eccentricity * eccentricity));
	const double tangential_speed = angular_momentum / boundary_radius;
	const double total_speed =
	    std::sqrt(mu * (2.0 / boundary_radius - 1.0 / semi_major_axis));
	const double radial_speed =
	    std::sqrt(total_speed * total_speed - tangential_speed * tangential_speed);
	Event outward(
	    3.0 * sec,
	    libphysica::Vector({boundary_radius, 0.0, 0.0}),
	    libphysica::Vector({radial_speed, tangential_speed, 0.0}));

	BoundKeplerExteriorArc arc;
	ASSERT_TRUE(Compute_Bound_Kepler_Exterior_Arc(outward, arc));
	EXPECT_TRUE(arc.outer_domain_removed);
	EXPECT_GT(arc.apoapsis_km, RADIAL_DOMAIN_MAX_KM);
	EXPECT_GT(arc.kepler_period_sec, arc.elapsed_time_sec);
	EXPECT_EQ(arc.dt_hist.size(), TOTAL_BINS);
	EXPECT_NEAR(
	    In_Units(arc.terminal_event.Radius(), km),
	    RADIAL_DOMAIN_MAX_KM,
	    1.0e-10 * RADIAL_DOMAIN_MAX_KM);
	EXPECT_GT(
	    arc.terminal_event.position.Dot(arc.terminal_event.velocity),
	    0.0);
	EXPECT_NEAR(
	    std::accumulate(arc.dt_hist.begin(), arc.dt_hist.end(), 0.0),
	    arc.elapsed_time_sec,
	    1.0e-9 * arc.elapsed_time_sec);
	EXPECT_GT(
	    std::accumulate(arc.v2dt_hist.begin(), arc.v2dt_hist.end(), 0.0),
	    0.0);
}

TEST(TestSimulationTrajectory, TestSurvivalInvalidTerminationReasons)
{
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::EnergyDriftEscape));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::NumericalFailure));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::NonFiniteState));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::SpeedLimit));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::Unknown));

	EXPECT_FALSE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::WallTimeLimit));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::MaxFreeSteps));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::MaxScatterings));
	EXPECT_FALSE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::OutwardEscape));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason::OuterDomainRemoval));
}

TEST(TestSimulationTrajectory, TestResidenceBincountRejectsOnlyNumericalTerminations)
{
	EXPECT_TRUE(TrajectoryTerminationInvalidatesResidenceBincount(
	    TrajectoryTerminationReason::NumericalFailure));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesResidenceBincount(
	    TrajectoryTerminationReason::NonFiniteState));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesResidenceBincount(
	    TrajectoryTerminationReason::SpeedLimit));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesResidenceBincount(
	    TrajectoryTerminationReason::EnergyDriftEscape));
	EXPECT_TRUE(TrajectoryTerminationInvalidatesResidenceBincount(
	    TrajectoryTerminationReason::Unknown));

	EXPECT_FALSE(TrajectoryTerminationInvalidatesResidenceBincount(
	    TrajectoryTerminationReason::OutwardEscape));
	EXPECT_FALSE(TrajectoryTerminationInvalidatesResidenceBincount(
	    TrajectoryTerminationReason::WallTimeLimit));
	EXPECT_FALSE(TrajectoryTerminationInvalidatesResidenceBincount(
	    TrajectoryTerminationReason::MaxFreeSteps));
	EXPECT_FALSE(TrajectoryTerminationInvalidatesResidenceBincount(
	    TrajectoryTerminationReason::MaxScatterings));
	EXPECT_FALSE(TrajectoryTerminationInvalidatesResidenceBincount(
	    TrajectoryTerminationReason::OuterDomainRemoval));
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
	EXPECT_FALSE(result.bincount.outer_domain_removed);
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
	EXPECT_EQ(
	    result.bincount.numerical_failure_detail,
	    TrajectoryNumericalFailureDetail::
	        UncapturedBoundMismatch);
	EXPECT_TRUE(stderr_output.empty());
	EXPECT_FALSE(result.bincount.is_captured);
	EXPECT_FALSE(result.bincount.survival_valid);
	EXPECT_EQ(result.number_of_scatterings, 0UL);
	EXPECT_DOUBLE_EQ(result.bincount.t_termination, 0.0);
}

TEST(TestSimulationTrajectory, TestUncapturedOutwardBoundaryStateTerminatesAsEscape)
{
    obscura::DM_Particle_SI DM(0.5 * GeV);

    // Use negligible but nonzero cross sections. This particle begins outside
    // the trajectory boundary and must terminate before any scattering query.
    // DM_Particle_SI does not permit sigma_proton = 0 with fixed fn/fp.
    DM.Set_Sigma_Proton(1.0e-100 * pb);
    DM.Set_Sigma_Electron(1.0e-100 * pb);

    Solar_Model SSM;
    const double boundary = 2.0 * rSun;
    Trajectory_Simulator simulator(SSM, 100, 10, boundary);

    const double radius = 1.01 * boundary;
    const double unbound_outward_speed =
        1.1 * SSM.Local_Escape_Speed(radius);

    Event IC(
        0.0,
        libphysica::Vector({radius, 0.0, 0.0}),
        libphysica::Vector({unbound_outward_speed, 0.0, 0.0}));

    Trajectory_Result result = simulator.Simulate(IC, DM, 0);

    EXPECT_EQ(
        result.bincount.termination_reason,
        TrajectoryTerminationReason::OutwardEscape);
    EXPECT_EQ(result.number_of_scatterings, 0UL);
    EXPECT_FALSE(result.bincount.is_captured);
    EXPECT_TRUE(result.bincount.survival_valid);
    EXPECT_FALSE(result.bincount.outer_domain_removed);
    EXPECT_TRUE(result.Particle_Free());
    EXPECT_FALSE(result.Particle_Reflected());
    EXPECT_DOUBLE_EQ(result.bincount.t_termination, 0.0);
    EXPECT_GT(result.final_event.Radius(), simulator.maximum_distance);
    EXPECT_GT(
        result.final_event.position.Dot(result.final_event.velocity),
        0.0);
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

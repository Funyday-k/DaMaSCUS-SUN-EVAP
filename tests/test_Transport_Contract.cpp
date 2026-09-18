#include "gtest/gtest.h"
#include "Simulation_Trajectory.hpp"
#include "obscura/DM_Particle_Standard.hpp"
#include "obscura/DM_Halo_Models.hpp"
#include <cmath>
#include <numeric>
#include <random>
using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

namespace {
Event ellipse(double e, double r, double peri) {
    const double a=peri/(1-e), mu=G_Newton*mSun;
    const double c=(1-r/a)/e, s=std::sqrt(1-c*c), beta=std::sqrt(1-e*e);
    const double n=std::sqrt(mu/(a*a*a));
    return Event(0.0,libphysica::Vector({a*(c-e),a*beta*s,0.0}),
        libphysica::Vector({-n*a*s/(1-e*c),n*a*beta*c/(1-e*c),0.0}));
}
double sum(const RadialHistogram& x){return std::accumulate(x.begin(),x.end(),0.0);}
}
TEST(TransportContract, GridRetainsNativeGeometryAndClipsFinalShell) {
    for(double outer:{1100.0,3000.0,10000.0,DEFAULT_OUTER_BOUNDARY_AU*AU_KM/R_SUN_KM}) {
        auto edges=BuildRadialGrid(outer*R_SUN_KM);
        EXPECT_DOUBLE_EQ(edges.front(),0.0); EXPECT_DOUBLE_EQ(edges.back(),outer*R_SUN_KM);
        for(std::size_t b=0;b+1<edges.size();++b) {
            EXPECT_LT(edges[b],edges[b+1]); EXPECT_DOUBLE_EQ(edges[b],BincountBinLowerKm(b));
            if(b<NUM_BINS) EXPECT_NEAR(edges[b+1]-edges[b],BIN_WIDTH_KM,1e-6);
        }
    }
    EXPECT_THROW(BuildRadialGrid(-1),std::invalid_argument);
}
TEST(TransportContract, EllipsePeriodAndVelocityMomentClosure) {
    for(double e:{0.5,0.9,0.99,0.999,0.9999}) {
        SCOPED_TRACE(e); const double peri=0.5*rSun, radius=1.2*rSun, a=peri/(1-e);
        const auto out=ellipse(e,radius,peri); BoundKeplerExteriorArc arc;
        ASSERT_TRUE(Compute_Bound_Kepler_Exterior_Arc(out,arc));
        const double n=std::sqrt(G_Newton*mSun/(a*a*a)), E=std::acos((1-radius/a)/e);
        const double inner_t=In_Units(2*(E-e*std::sin(E))/n,sec);
        const double inner_v2=In_Units(2*G_Newton*mSun/(a*n)*(E+e*std::sin(E)),km*km/sec);
        const double period=In_Units(2*M_PI/n,sec);
        const double full_v2=In_Units(2*M_PI*G_Newton*mSun/(a*n),km*km/sec);
        EXPECT_NEAR((sum(arc.dt_hist)+inner_t)/period,1.0,2e-8);
        EXPECT_NEAR((sum(arc.v2dt_hist)+inner_v2)/full_v2,1.0,2e-8);
        EXPECT_NEAR(sum(arc.dt_hist)/arc.elapsed_time_sec,1.0,2e-8);
    }
}
TEST(TransportContract, RemovalIsOneWayAndPreservesEnergyAngularMomentum) {
    const double peri=0.5*rSun, apo=4000*rSun, e=(apo-peri)/(apo+peri);
    const auto out=ellipse(e,1.2*rSun,peri); BoundKeplerExteriorArc arc;
    ASSERT_TRUE(Compute_Bound_Kepler_Exterior_Arc(out,arc,1100*R_SUN_KM));
    ASSERT_TRUE(arc.outer_removed); EXPECT_NEAR(In_Units(arc.terminal_event.Radius(),rSun),1100,1e-7);
    EXPECT_GT(arc.terminal_event.position.Dot(arc.terminal_event.velocity),0);
    EXPECT_NEAR(sum(arc.dt_hist)/arc.elapsed_time_sec,1.0,2e-8);
    EXPECT_NEAR(arc.terminal_event.Angular_Momentum()/out.Angular_Momentum(),1.0,2e-8);
    EXPECT_NEAR((arc.terminal_event.Speed()*arc.terminal_event.Speed()-2*G_Newton*mSun/arc.terminal_event.Radius()) /
                (out.Speed()*out.Speed()-2*G_Newton*mSun/out.Radius()),1.0,2e-8);
    BoundKeplerExteriorArc complete; ASSERT_TRUE(Compute_Bound_Kepler_Exterior_Arc(out,complete));
    EXPECT_LT(arc.elapsed_time_sec,complete.elapsed_time_sec/2);
    EXPECT_FALSE(TrajectoryTerminationInvalidatesResidenceBincount(TrajectoryTerminationReason::OuterDomainRemoval));
}
TEST(TransportContract, WideBoundOrbitStopsAtTheShared1100AUSurface) {
    const double peri=0.5*rSun, apo=2000*AU;
    const double e=(apo-peri)/(apo+peri);
    const Event out=ellipse(e,1.2*rSun,peri);
    BoundKeplerExteriorArc arc;
    ASSERT_TRUE(Compute_Bound_Kepler_Exterior_Arc(out,arc,DEFAULT_OUTER_BOUNDARY_AU*AU_KM));
    EXPECT_TRUE(arc.outer_removed);
    EXPECT_NEAR(In_Units(arc.terminal_event.Radius(),AU),DEFAULT_OUTER_BOUNDARY_AU,1e-7);
    EXPECT_GT(arc.terminal_event.position.Dot(arc.terminal_event.velocity),0.0);
    EXPECT_NEAR(sum(arc.dt_hist)/arc.elapsed_time_sec,1.0,1e-10);
    EXPECT_GT(arc.dt_hist.back(),0.0);
    EXPECT_GT(sum(arc.v2dt_hist),0.0);
}
TEST(TransportContract, HyperbolicMomentsAgreeWithRadialQuadrature) {
    const auto out=Event(0.0,libphysica::Vector({1.2*rSun,0.0,0.0}),libphysica::Vector({600*km/sec,150*km/sec,0.0}));
    BoundKeplerExteriorArc arc; ASSERT_TRUE(Compute_Unbound_Kepler_Exterior_Arc(out,1100*R_SUN_KM,arc));
    const double mu=In_Units(G_Newton*mSun,km*km*km/(sec*sec)), r0=1.2*R_SUN_KM;
    const double v=In_Units(out.Speed(),km/sec), h=In_Units(out.Angular_Momentum(),km*km/sec), energy=v*v/2-mu/r0;
    double dt=0,dv2=0; const int N=200000; const double dy=std::log(1100/1.2)/N;
    for(int j=0;j<N;++j) {
        const double r=r0*std::exp((j+0.5)*dy), v2=2*(energy+mu/r), vr=std::sqrt(v2-h*h/(r*r));
        dt+=r*dy/vr; dv2+=v2*r*dy/vr;
    }
    EXPECT_NEAR(sum(arc.dt_hist)/dt,1.0,1e-8); EXPECT_NEAR(sum(arc.v2dt_hist)/dv2,1.0,1e-8);
    EXPECT_NEAR(In_Units(arc.terminal_event.time,sec)/dt,1.0,1e-8);
}
TEST(TransportContract, Incoming1100AUPathClosesAtTheSolarSurface) {
    Solar_Model sun;
    obscura::Standard_Halo_Model halo;
    std::mt19937 prng(20260918u);
    for(int i=0;i<32;++i) {
        const Event sampled=Initial_Conditions(halo,sun,prng);
        EXPECT_NEAR(In_Units(sampled.Radius(),AU),DEFAULT_OUTER_BOUNDARY_AU,1e-9);
        Event surface=sampled;
        ASSERT_TRUE(Hyperbolic_Kepler_Shift(surface,rSun));
        EXPECT_NEAR(In_Units(surface.Radius(),rSun),1.0,1e-10);
        EXPECT_LT(surface.position.Dot(surface.velocity),0.0);
        BoundKeplerExteriorArc inbound;
        ASSERT_TRUE(Compute_Unbound_Kepler_Exterior_Arc(
            Event(0.0,surface.position,(-1.0)*surface.velocity),
            DEFAULT_OUTER_BOUNDARY_AU*AU_KM,inbound));
        EXPECT_NEAR(sum(inbound.dt_hist)/inbound.elapsed_time_sec,1.0,1e-10);
        EXPECT_GT(inbound.dt_hist[BincountBinIndexKm(R_SUN_KM)],0.0);
        EXPECT_GT(inbound.dt_hist.back(),0.0);
        EXPECT_NEAR(In_Units(inbound.terminal_event.Radius(),AU),DEFAULT_OUTER_BOUNDARY_AU,1e-9);
        EXPECT_NEAR(inbound.terminal_event.Angular_Momentum()/sampled.Angular_Momentum(),1.0,1e-9);
    }
}
TEST(TransportContract, SolarPotentialMatchesExteriorAtTheSurface) {
    Solar_Model sun;
    const double analytic=2*G_Newton*mSun/rSun;
    EXPECT_NEAR(std::pow(sun.Local_Escape_Speed(rSun),2)/analytic,1.0,1e-10);
    EXPECT_NEAR(std::pow(sun.Local_Escape_Speed((1-1e-11)*rSun),2)/analytic,1.0,1e-9);
    EXPECT_NEAR(sun.Mass(rSun)/mSun,1.0,1e-10);
}

TEST(TransportContract, CaptureAndTransportShareFirstCollisionAccuracy) {
    // Low positive incident energy exposed the former coarse capture-only
    // collision-location approximation. Compare the same first-capture state
    // in the two workflows while allowing transport to continue afterwards.
    Solar_Model sun;
    obscura::DM_Particle_SD dm(0.01*GeV);
    dm.Set_Sigma_Proton(1e-32*cm*cm);
    const Event initial(0.0,
        libphysica::Vector({-742895.7231921622*km,183547.0308320915*km,-7389.508476071721*km}),
        libphysica::Vector({115.6738867389554*km/sec,-529.697120864836*km/sec,230.2163702333374*km/sec}));
    unsigned captures=0;
    for(unsigned seed=1;seed<=8;++seed) {
        SCOPED_TRACE(seed);
        Trajectory_Simulator capture(sun,100000,128,1.1*rSun), transport(sun,100000,128,1.1*rSun);
        capture.Enable_Capture_Mode(true); capture.Fix_PRNG_Seed(seed); transport.Fix_PRNG_Seed(seed);
        const auto a=capture.Simulate(initial,dm,0), b=transport.Simulate(initial,dm,0);
        if(b.bincount.is_captured) {
            double inside=0,outside=0;
            for(std::size_t bin=0;bin<b.bincount.dt_hist.size();++bin)
                (BincountBinUpperKm(bin)<=R_SUN_KM*(1+1e-12) ? inside : outside) += b.bincount.dt_hist[bin];
            EXPECT_NEAR(b.bincount.time_inside_sun_after_capture_sec,inside,1e-10*std::max(1.0,inside));
            EXPECT_NEAR(b.bincount.time_outside_sun_after_capture_sec,outside,1e-10*std::max(1.0,outside));
        }
        EXPECT_NE(a.bincount.termination_reason,TrajectoryTerminationReason::NumericalFailure);
        EXPECT_EQ(a.bincount.is_captured,b.bincount.is_captured);
        if(a.bincount.is_captured) {
            ++captures;
            EXPECT_DOUBLE_EQ(a.bincount.t_capture,b.bincount.t_capture);
            EXPECT_DOUBLE_EQ(a.bincount.r_first_negative_km,b.bincount.r_first_negative_km);
            EXPECT_DOUBLE_EQ(a.bincount.E_first_negative_eV,b.bincount.E_first_negative_eV);
        }
    }
    EXPECT_GT(captures,0u);
}

TEST(TransportContract, NearlyParabolicMomentsRemainAccurateAtFiniteCutoff) {
    // Integrate dr/vr directly, independently of eccentric/hyperbolic anomalies.
    const double r0=1.2*R_SUN_KM, r1=1100*R_SUN_KM;
    const double mu=In_Units(G_Newton*mSun,km*km*km/(sec*sec)), vt=150.0;
    for(double offset:{-1e-4,-1e-8,-1e-10,-1e-12,1e-4,1e-8,1e-10,1e-12}) {
        SCOPED_TRACE(offset);
        const double vr0=std::sqrt(2*mu/r0*(1+offset)-vt*vt);
        const Event out(0.0,libphysica::Vector({r0*km,0.0,0.0}),libphysica::Vector({vr0*km/sec,vt*km/sec,0.0}));
        BoundKeplerExteriorArc arc;
        ASSERT_TRUE(offset<0 ? Compute_Bound_Kepler_Exterior_Arc(out,arc,r1)
                             : Compute_Unbound_Kepler_Exterior_Arc(out,r1,arc));
        const double v=In_Units(out.Speed(),km/sec), h=In_Units(out.Angular_Momentum(),km*km/sec);
        const double energy=v*v/2-mu/r0;
        const int N=100000; const double dy=std::log(r1/r0)/N;
        double dt=0,v2dt=0;
        for(int j=0;j<N;++j) {
            const double r=r0*std::exp((j+0.5)*dy), v2=2*(energy+mu/r);
            const double step=r*dy/std::sqrt(v2-h*h/(r*r));
            dt+=step; v2dt+=v2*step;
        }
        EXPECT_NEAR(sum(arc.dt_hist)/dt,1.0,1e-8);
        EXPECT_NEAR(sum(arc.v2dt_hist)/v2dt,1.0,1e-8);
        EXPECT_NEAR(arc.elapsed_time_sec/dt,1.0,1e-8);
        EXPECT_NEAR(In_Units(arc.terminal_event.Radius(),km)/r1,1.0,1e-8);
    }
}

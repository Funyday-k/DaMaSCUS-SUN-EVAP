#include "Simulation_Trajectory.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/Astronomy.hpp"
#include "Snapshot_Shared_State.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

namespace
{
constexpr double RK45_MIN_STEP_FACTOR = 0.1;
constexpr double RK45_MAX_STEP_FACTOR = 4.0;
// 单次 Runge_Kutta_45_Step 内层 while(!accepted) 的最大重试次数。
// 正常物理场景下通常在 100 次内收敛；保留额外裕量，同时避免一次病态
// 步长重试长时间占住外层，使轨迹墙钟保护无法及时获得控制权。
constexpr int RK45_MAX_INNER_RETRIES = 256;
constexpr double FREE_ENERGY_DRIFT_REL_E_SCALE_EV = 1.0e-30;
constexpr double BOUNDARY_ENERGY_ABSOLUTE_TOLERANCE_EV_AT_10_MEV = 1.0e-4;
constexpr double BOUNDARY_ENERGY_RELATIVE_TOLERANCE = 1.0e-8;
constexpr unsigned int MAX_BOUNDARY_REFINEMENT_RETRIES = 32;
constexpr int BOUNDARY_HERMITE_SCAN_INTERVALS = 64;
constexpr int BOUNDARY_HERMITE_BISECTION_ITERATIONS = 80;
constexpr double MAX_OPTICAL_DEPTH_STEP = 0.05;
constexpr double CAPTURE_MODE_MAX_OPTICAL_DEPTH_STEP = 0.10;
constexpr double OPTICAL_DEPTH_RELATIVE_TOLERANCE = 1.0e-2;
constexpr double OPTICAL_DEPTH_ABSOLUTE_TOLERANCE = 1.0e-12 * MAX_OPTICAL_DEPTH_STEP;
constexpr unsigned int MAX_OPTICAL_DEPTH_RETRIES = 100;
constexpr std::size_t MAX_OPTICAL_DEPTH_PIECES = 4;
constexpr std::size_t CAPTURE_MODE_OPTICAL_DEPTH_PIECES = 2;
constexpr unsigned long int TARGET_VELOCITY_MAX_REJECTION_ATTEMPTS = 10000UL;
// Publish at either bound: the step limit keeps fast trajectories inexpensive,
// while the wall-clock limit keeps slow trajectories visible to snapshots.
constexpr unsigned long int SNAPSHOT_PUBLISH_STEP_INTERVAL = 512UL;
constexpr double SNAPSHOT_PUBLISH_WALL_INTERVAL_SEC = 0.25;
constexpr double BINCOUNT_DENSE_POSITION_TOLERANCE_KM = 2.0e-3 * BIN_WIDTH_KM;
constexpr int BINCOUNT_DENSE_MAX_RECURSION = 20;
constexpr double BINCOUNT_GAUSS_LEGENDRE_OFFSET = 0.77459666924148337704;

struct OpticalDepthPiece
{
	double start;
	double end;
	double rate_start;
	double rate_end;
};

double RK45_Absolute_Max_Time_Step()
{
	return 1.0e6 * sec;  // 约 11.6 天；防止太阳外弱加速度区域步长增长到 inf。
}

double RK45_Sanitized_Time_Step(double step)
{
	if(!std::isfinite(step) || step <= 0.0)
		return 0.1 * sec;
	return std::min(step, RK45_Absolute_Max_Time_Step());
}

double Free_Propagation_Time_Step_Cap(double radius, double speed, double maximum_distance)
{
	double cap = RK45_Absolute_Max_Time_Step();
	const double safe_speed = std::max(std::fabs(speed), 1.0e-12 * km / sec);

	// 避免单步跨越过大的径向距离，否则可能直接跳过 1.1R_sun 边界并把 r 推到非物理值。
	const double crossing_scale = std::max(0.25 * maximum_distance, 10.0 * km);
	cap = std::min(cap, crossing_scale / safe_speed);

	// 太阳外近似 Kepler 运动，步长不应大于局域动力学时间的一个固定比例。
	const double safe_radius = std::max(radius, 1.0 * km);
	const double dynamical_time = sqrt(safe_radius * safe_radius * safe_radius / (G_Newton * mSun));
	if(std::isfinite(dynamical_time) && dynamical_time > 0.0)
		cap = std::min(cap, 0.1 * dynamical_time);

	return std::max(cap, 1.0e-8 * sec);
}

double Clamp_Unit_Interval(double value)
{
	return std::max(0.0, std::min(1.0, value));
}

double Clamp_Cosine(double value)
{
	return std::max(-1.0, std::min(1.0, value));
}

double Positive_Unit_Uniform(std::mt19937& PRNG)
{
	return std::max(libphysica::Sample_Uniform(PRNG, 0.0, 1.0),
	                std::numeric_limits<double>::min());
}

libphysica::Vector Unit_Vector_At_Angle_From_Axis(double cos_theta, double phi, const libphysica::Vector& axis)
{
	if(!std::isfinite(cos_theta) || cos_theta < -1.0 - 1.0e-12 || cos_theta > 1.0 + 1.0e-12
	   || !std::isfinite(phi))
		throw std::runtime_error("Unit_Vector_At_Angle_From_Axis(): angle is non-finite or outside its domain.");
	const double axis_norm = axis.Norm();
	if(!std::isfinite(axis_norm) || axis_norm <= 0.0)
		throw std::runtime_error("Unit_Vector_At_Angle_From_Axis(): axis is zero or non-finite.");

	libphysica::Vector e3 = axis / axis_norm;
	libphysica::Vector reference = (std::fabs(e3[2]) < 0.9)
	                             ? libphysica::Vector({0.0, 0.0, 1.0})
	                             : libphysica::Vector({1.0, 0.0, 0.0});
	libphysica::Vector e1 = reference.Cross(e3).Normalized();
	libphysica::Vector e2 = e3.Cross(e1).Normalized();

	cos_theta = Clamp_Cosine(cos_theta);
	const double sin_theta = sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
	return cos_theta * e3 + sin_theta * cos(phi) * e1 + sin_theta * sin(phi) * e2;
}

void Build_Stable_Orthonormal_Basis(const libphysica::Vector& axis, libphysica::Vector& e1, libphysica::Vector& e2, libphysica::Vector& e3)
{
	const double axis_norm = axis.Norm();
	if(!std::isfinite(axis_norm) || axis_norm <= 0.0)
		throw std::runtime_error("Build_Stable_Orthonormal_Basis(): axis is zero or non-finite.");

	e1 = axis / axis_norm;
	libphysica::Vector reference = (std::fabs(e1[2]) < 0.9)
	                             ? libphysica::Vector({0.0, 0.0, 1.0})
	                             : libphysica::Vector({1.0, 0.0, 0.0});
	e2 = reference.Cross(e1).Normalized();
	e3 = e1.Cross(e2).Normalized();
}

double Radial_Velocity(const Event& event)
{
	const double radius = event.Radius();
	if(radius <= 0.0)
		return 0.0;
	return event.position.Dot(event.velocity) / radius;
}

bool Bound_Kepler_Return_At_Same_Radius(const Event& outward_event, Event& inbound_event,
	                                    double& return_time, double& apoapsis_radius)
{
	// Bound particles outside the Sun can spend very long times on Kepler arcs.
	// Fast-forward only this scatter-free exterior segment; this is not a
	// termination or censoring condition, the trajectory continues inbound.
	const double radius = outward_event.Radius();
	const double speed = outward_event.Speed();
	const double radial_velocity = Radial_Velocity(outward_event);
	if(!std::isfinite(radius) || !std::isfinite(speed) || !std::isfinite(radial_velocity)
	   || radius <= rSun || speed <= 0.0 || radial_velocity <= 0.0)
		return false;

	const double mu = G_Newton * mSun;
	const double specific_energy = 0.5 * speed * speed - mu / radius;
	if(!std::isfinite(specific_energy) || specific_energy >= 0.0)
		return false;

	const double angular_momentum = outward_event.Angular_Momentum();
	if(!std::isfinite(angular_momentum))
		return false;

	const double semi_major_axis = -mu / (2.0 * specific_energy);
	if(!std::isfinite(semi_major_axis) || semi_major_axis <= 0.0)
		return false;

	double eccentricity_sqr = 1.0 + 2.0 * specific_energy * angular_momentum * angular_momentum / (mu * mu);
	if(!std::isfinite(eccentricity_sqr) || eccentricity_sqr < 0.0)
		return false;
	eccentricity_sqr = std::max(0.0, eccentricity_sqr);
	const double eccentricity = sqrt(eccentricity_sqr);
	if(eccentricity <= 1.0e-14 || eccentricity >= 1.0)
		return false;

	const double cos_eccentric_anomaly = Clamp_Cosine((1.0 - radius / semi_major_axis) / eccentricity);
	const double eccentric_anomaly = acos(cos_eccentric_anomaly);
	const double sin_eccentric_anomaly = sqrt(std::max(0.0, 1.0 - cos_eccentric_anomaly * cos_eccentric_anomaly));
	const double mean_motion = sqrt(mu / (semi_major_axis * semi_major_axis * semi_major_axis));
	if(!std::isfinite(mean_motion) || mean_motion <= 0.0)
		return false;

	return_time = (2.0 * M_PI - 2.0 * eccentric_anomaly + 2.0 * eccentricity * sin_eccentric_anomaly) / mean_motion;
	if(!std::isfinite(return_time) || return_time <= 0.0)
		return false;
	apoapsis_radius = semi_major_axis * (1.0 + eccentricity);
	if(!std::isfinite(apoapsis_radius) || apoapsis_radius < radius)
		return false;

	const libphysica::Vector radial_unit = outward_event.position / radius;
	inbound_event = outward_event;
	inbound_event.time += return_time;
	inbound_event.velocity = outward_event.velocity - 2.0 * radial_velocity * radial_unit;
	return inbound_event.time > outward_event.time && std::isfinite(inbound_event.Speed());
}

// Overwrite an Event in place. Assigning a libphysica::Vector allocates twice
// (its assignment operator takes and returns by value), so the propagation loop
// copies components into buffers it already owns.
void Copy_Event_Into(Event& destination, const Event& source)
{
	if(destination.position.Size() != 3)
		destination.position = libphysica::Vector(3);
	if(destination.velocity.Size() != 3)
		destination.velocity = libphysica::Vector(3);
	destination.time = source.time;
	for(unsigned int component = 0; component < 3; component++)
	{
		destination.position[component] = source.position[component];
		destination.velocity[component] = source.velocity[component];
	}
}

Event Interpolate_Event(const Event& before, const Event& after, double fraction)
{
	fraction = Clamp_Unit_Interval(fraction);
	libphysica::Vector position(3);
	libphysica::Vector velocity(3);
	for(unsigned int component = 0; component < 3; component++)
	{
		position[component] = before.position[component]
		                    + fraction * (after.position[component] - before.position[component]);
		velocity[component] = before.velocity[component]
		                    + fraction * (after.velocity[component] - before.velocity[component]);
	}
	return Event(before.time + fraction * (after.time - before.time), position, velocity);
}

Event Hermite_Interpolate_Event(const Event& before, const Event& after, double fraction)
{
	fraction = Clamp_Unit_Interval(fraction);
	const double dt = after.time - before.time;
	const double s = fraction;
	const double s2 = s * s;
	const double s3 = s2 * s;
	const double h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
	const double h10 = s3 - 2.0 * s2 + s;
	const double h01 = -2.0 * s3 + 3.0 * s2;
	const double h11 = s3 - s2;
	const double dh00 = 6.0 * s2 - 6.0 * s;
	const double dh10 = 3.0 * s2 - 4.0 * s + 1.0;
	const double dh01 = -6.0 * s2 + 6.0 * s;
	const double dh11 = 3.0 * s2 - 2.0 * s;

	libphysica::Vector position(3);
	libphysica::Vector velocity(3);
	for(unsigned int component = 0; component < 3; component++)
	{
		position[component] =
		    h00 * before.position[component]
		    + h10 * dt * before.velocity[component]
		    + h01 * after.position[component]
		    + h11 * dt * after.velocity[component];
		if(dt > 0.0)
		{
			velocity[component] =
			    dh00 / dt * before.position[component]
			    + dh10 * before.velocity[component]
			    + dh01 / dt * after.position[component]
			    + dh11 * after.velocity[component];
		}
		else
			velocity[component] = before.velocity[component];
	}
	return Event(before.time + fraction * dt, position, velocity);
}

double Hermite_Radius_Residual(
	const Event& before,
	const Event& after,
	double fraction,
	double target_radius)
{
	return Hermite_Interpolate_Event(before, after, fraction).Radius()
	     - target_radius;
}

double Boundary_Energy_Tolerance_eV(
	double reference_energy_eV,
	const obscura::DM_Particle& DM)
{
	const double mass_scale =
	    std::max(In_Units(DM.mass, MeV) / 10.0, 1.0e-12);
	const double absolute_tolerance =
	    BOUNDARY_ENERGY_ABSOLUTE_TOLERANCE_EV_AT_10_MEV * mass_scale;
	return absolute_tolerance
	     + BOUNDARY_ENERGY_RELATIVE_TOLERANCE
	       * std::max(std::fabs(reference_energy_eV), absolute_tolerance);
}

double Radius_At_Fraction(const Event& before, const Event& after, double fraction)
{
	fraction = Clamp_Unit_Interval(fraction);
	double radius_squared = 0.0;
	for(unsigned int component = 0; component < 3; component++)
	{
		const double value = before.position[component]
		                   + fraction * (after.position[component] - before.position[component]);
		radius_squared += value * value;
	}
	return sqrt(std::max(0.0, radius_squared));
}

bool Find_First_Outward_Hermite_Radius_Crossing_Impl(
	const Event& before,
	const Event& after,
	double target_radius,
	Event& crossing)
{
	if(!(after.time > before.time)
	   || !std::isfinite(target_radius)
	   || target_radius <= 0.0)
		return false;

	double left = 0.0;
	double residual_left =
	    Hermite_Radius_Residual(before, after, left, target_radius);
	for(int interval = 1; interval <= BOUNDARY_HERMITE_SCAN_INTERVALS; interval++)
	{
		const double right =
		    static_cast<double>(interval)
		    / static_cast<double>(BOUNDARY_HERMITE_SCAN_INTERVALS);
		const double residual_right =
		    Hermite_Radius_Residual(before, after, right, target_radius);
		const bool bracketed =
		    std::isfinite(residual_left)
		    && std::isfinite(residual_right)
		    && ((residual_left <= 0.0 && residual_right >= 0.0)
		        || (residual_left >= 0.0 && residual_right <= 0.0));
		if(bracketed)
		{
			double lower = left;
			double upper = right;
			double lower_residual = residual_left;
			for(int iteration = 0;
			    iteration < BOUNDARY_HERMITE_BISECTION_ITERATIONS;
			    iteration++)
			{
				const double middle = 0.5 * (lower + upper);
				const double middle_residual =
				    Hermite_Radius_Residual(
				        before, after, middle, target_radius);
				if(!std::isfinite(middle_residual))
					return false;
				if((lower_residual <= 0.0 && middle_residual <= 0.0)
				   || (lower_residual >= 0.0 && middle_residual >= 0.0))
				{
					lower = middle;
					lower_residual = middle_residual;
				}
				else
					upper = middle;
			}
			Event candidate =
			    Hermite_Interpolate_Event(
			        before, after, 0.5 * (lower + upper));
			const double candidate_radius = candidate.Radius();
			if(std::isfinite(candidate_radius)
			   && candidate_radius > 0.0)
				candidate.position =
				    target_radius / candidate_radius * candidate.position;
			if(std::isfinite(candidate.Radius())
			   && std::isfinite(candidate.Speed())
			   && Radial_Velocity(candidate) > 0.0)
			{
				crossing = candidate;
				return true;
			}
		}
		left = right;
		residual_left = residual_right;
	}
	return false;
}

bool Surface_Crossing_Fractions(const Event& before, const Event& after, double& first, double& second)
{
	double a = 0.0;
	double b = 0.0;
	double c = -rSun * rSun;
	for(unsigned int component = 0; component < 3; component++)
	{
		const double displacement = after.position[component] - before.position[component];
		a += displacement * displacement;
		b += 2.0 * before.position[component] * displacement;
		c += before.position[component] * before.position[component];
	}
	if(a <= 0.0)
		return false;

	const double discriminant = b * b - 4.0 * a * c;
	if(!std::isfinite(discriminant) || discriminant < 0.0)
		return false;

	const double sqrt_discriminant = std::sqrt(std::max(0.0, discriminant));
	first = (-b - sqrt_discriminant) / (2.0 * a);
	second = (-b + sqrt_discriminant) / (2.0 * a);
	if(first > second)
		std::swap(first, second);
	return true;
}

bool Fraction_In_Unit_Interval(double fraction)
{
	constexpr double tolerance = 1.0e-12;
	return fraction >= -tolerance && fraction <= 1.0 + tolerance;
}

double Find_Surface_Crossing_Fraction(const Event& before, const Event& after, double lower, double upper)
{
	double r_lower = Radius_At_Fraction(before, after, lower) - rSun;
	for(int iteration = 0; iteration < 60; iteration++)
	{
		const double mid = 0.5 * (lower + upper);
		const double r_mid = Radius_At_Fraction(before, after, mid) - rSun;
		if((r_lower <= 0.0 && r_mid <= 0.0) || (r_lower >= 0.0 && r_mid >= 0.0))
		{
			lower = mid;
			r_lower = r_mid;
		}
		else
			upper = mid;
	}
	return Clamp_Unit_Interval(0.5 * (lower + upper));
}

bool Solar_Interior_Fraction_Interval(const Event& before, const Event& after, double r_before, double r_after, double& start, double& end)
{
	const bool before_inside = r_before < rSun;
	const bool after_inside = r_after < rSun;
	double first_crossing = 0.0;
	double second_crossing = 0.0;
	const bool has_crossings = Surface_Crossing_Fractions(before, after, first_crossing, second_crossing);
	if(before_inside && after_inside)
	{
		start = 0.0;
		end = 1.0;
		return true;
	}
	if(!before_inside && after_inside)
	{
		if(has_crossings && Fraction_In_Unit_Interval(first_crossing))
			start = Clamp_Unit_Interval(first_crossing);
		else if(has_crossings && Fraction_In_Unit_Interval(second_crossing))
			start = Clamp_Unit_Interval(second_crossing);
		else
			start = Find_Surface_Crossing_Fraction(before, after, 0.0, 1.0);
		end = 1.0;
		return end > start;
	}
	if(before_inside && !after_inside)
	{
		start = 0.0;
		if(has_crossings && Fraction_In_Unit_Interval(second_crossing))
			end = Clamp_Unit_Interval(second_crossing);
		else if(has_crossings && Fraction_In_Unit_Interval(first_crossing))
			end = Clamp_Unit_Interval(first_crossing);
		else
			end = Find_Surface_Crossing_Fraction(before, after, 0.0, 1.0);
		return end > start;
	}

	if(has_crossings && Fraction_In_Unit_Interval(first_crossing) && Fraction_In_Unit_Interval(second_crossing))
	{
		start = Clamp_Unit_Interval(first_crossing);
		end = Clamp_Unit_Interval(second_crossing);
		return end > start;
	}
	return false;
}

double Scattering_Rate_At_Fraction(const Event& before, const Event& after, double fraction, Solar_Model& solar_model, obscura::DM_Particle& DM)
{
	fraction = Clamp_Unit_Interval(fraction);
	double radius_squared = 0.0;
	double speed_squared = 0.0;
	for(unsigned int component = 0; component < 3; component++)
	{
		const double position = before.position[component]
		                      + fraction * (after.position[component] - before.position[component]);
		const double velocity = before.velocity[component]
		                      + fraction * (after.velocity[component] - before.velocity[component]);
		radius_squared += position * position;
		speed_squared += velocity * velocity;
	}
	const double radius = sqrt(std::max(0.0, radius_squared));
	if(radius >= rSun)
		return 0.0;
	return solar_model.Total_DM_Scattering_Rate(DM, radius, sqrt(std::max(0.0, speed_squared)));
}

bool Build_Optical_Depth_Pieces(const Event& before,
                                const Event& after,
                                double interior_start,
                                double interior_end,
                                double actual_dt,
                                Solar_Model& solar_model,
                                obscura::DM_Particle& DM,
                                bool use_cached_start_rate,
                                double cached_start_rate,
                                std::size_t requested_piece_count,
                                std::array<OpticalDepthPiece, MAX_OPTICAL_DEPTH_PIECES>& pieces,
                                std::size_t& piece_count,
                                double& tau_two_piece,
                                double& tau_four_piece,
                                bool& invalid_rate)
{
	piece_count = 0;
	tau_two_piece = 0.0;
	tau_four_piece = 0.0;
	invalid_rate = false;
	if(!(interior_end > interior_start) || actual_dt <= 0.0)
		return false;

	requested_piece_count = std::max<std::size_t>(1, std::min(requested_piece_count, MAX_OPTICAL_DEPTH_PIECES));
	const double span = interior_end - interior_start;
	std::array<double, MAX_OPTICAL_DEPTH_PIECES + 1> fractions;
	std::array<double, MAX_OPTICAL_DEPTH_PIECES + 1> rates;
	for(std::size_t i = 0; i <= requested_piece_count; i++)
		fractions[i] = interior_start + static_cast<double>(i) / static_cast<double>(requested_piece_count) * span;

	for(std::size_t i = 0; i <= requested_piece_count; i++)
	{
		if(i == 0 && use_cached_start_rate)
			rates[i] = cached_start_rate;
		else
			rates[i] = Scattering_Rate_At_Fraction(before, after, fractions[i], solar_model, DM);

		if(!std::isfinite(rates[i]) || rates[i] < 0.0)
		{
			invalid_rate = true;
			return true;
		}
	}

	const auto trapezoid_tau = [&](std::size_t start_index, std::size_t end_index)
	{
		const double piece_dt = actual_dt * (fractions[end_index] - fractions[start_index]);
		return 0.5 * (rates[start_index] + rates[end_index]) * piece_dt;
	};

	if(requested_piece_count >= 4)
		tau_two_piece = trapezoid_tau(0, 2) + trapezoid_tau(2, 4);
	else
	{
		for(std::size_t i = 0; i < requested_piece_count; i++)
			tau_two_piece += trapezoid_tau(i, i + 1);
	}
	for(std::size_t i = 0; i < requested_piece_count; i++)
	{
		if(fractions[i + 1] > fractions[i])
		{
			pieces[piece_count++] = OpticalDepthPiece{fractions[i], fractions[i + 1], rates[i], rates[i + 1]};
			tau_four_piece += trapezoid_tau(i, i + 1);
		}
	}

	return piece_count > 0;
}

double Solve_Linear_Rate_Optical_Depth_Fraction(double rate_start, double rate_end, double interval_dt, double target_tau)
{
	if(target_tau <= 0.0)
		return 0.0;
	if(!std::isfinite(target_tau) || interval_dt <= 0.0)
		return 1.0;
	if(!std::isfinite(rate_start) || !std::isfinite(rate_end) || rate_start < 0.0 || rate_end < 0.0)
		return 1.0;

	const double total_tau = 0.5 * (rate_start + rate_end) * interval_dt;
	if(!std::isfinite(total_tau) || total_tau <= 0.0)
		return 1.0;
	if(target_tau >= total_tau)
		return 1.0;

	const double target = target_tau / interval_dt;
	const double a = 0.5 * (rate_end - rate_start);
	const double b = rate_start;
	if(std::fabs(a) < 1.0e-300)
		return (b > 0.0) ? Clamp_Unit_Interval(target / b) : 1.0;

	const double discriminant = b * b + 4.0 * a * target;
	if(std::isfinite(discriminant) && discriminant >= 0.0)
	{
		const double sqrt_discriminant = sqrt(discriminant);
		const double roots[2] = {
		    (-b + sqrt_discriminant) / (2.0 * a),
		    (-b - sqrt_discriminant) / (2.0 * a)
		};
		for(double root : roots)
		{
			if(std::isfinite(root) && root >= 0.0 && root <= 1.0)
				return Clamp_Unit_Interval(root);
		}
	}

	double low = 0.0;
	double high = 1.0;
	for(int iteration = 0; iteration < 30; iteration++)
	{
		const double mid = 0.5 * (low + high);
		const double tau_mid = interval_dt * (rate_start * mid + 0.5 * (rate_end - rate_start) * mid * mid);
		if(tau_mid < target_tau)
			low = mid;
		else
			high = mid;
	}
	return Clamp_Unit_Interval(0.5 * (low + high));
}

bool RK45_Errors_Within_Tolerance(const double errors[3], const double tolerances[3])
{
	for(int i = 0; i < 3; i++)
	{
		if(!std::isfinite(errors[i]) || errors[i] > tolerances[i])
			return false;
	}

	return true;
}

double RK45_Next_Step_Size(double current_step, const double errors[3], const double tolerances[3])
{
	// 若任何一个 error 为 NaN/Inf，说明 RK 中间阶段出现数值崩溃。
	// 此时必须让步长"缩小"而不是"放大"，否则接下来会继续崩，外层
	// while(!accepted) 会陷入无限循环（见 DaMaSCUS-SUN-EVAP#rank1-stuck）。
	for(int i = 0; i < 3; i++)
	{
		if(!std::isfinite(errors[i]))
			return RK45_MIN_STEP_FACTOR * current_step;
	}

	double factor = RK45_MAX_STEP_FACTOR;

	for(int i = 0; i < 3; i++)
	{
		if(errors[i] <= 0.0)
			continue;

		const double candidate = 0.84 * pow(tolerances[i] / errors[i], 0.25);
		if(std::isfinite(candidate))
			factor = std::min(factor, candidate);
	}

	factor = std::max(RK45_MIN_STEP_FACTOR, std::min(factor, RK45_MAX_STEP_FACTOR));
	return RK45_Sanitized_Time_Step(factor * current_step);
}

using Cartesian3 = std::array<double, 3>;

Cartesian3 Cartesian_Add(const Cartesian3& lhs, const Cartesian3& rhs)
{
	Cartesian3 result;
	for(std::size_t component = 0; component < result.size(); component++)
		result[component] = lhs[component] + rhs[component];
	return result;
}

Cartesian3 Cartesian_Subtract(const Cartesian3& lhs, const Cartesian3& rhs)
{
	Cartesian3 result;
	for(std::size_t component = 0; component < result.size(); component++)
		result[component] = lhs[component] - rhs[component];
	return result;
}

Cartesian3 Cartesian_Scale(double scale, const Cartesian3& value)
{
	Cartesian3 result;
	for(std::size_t component = 0; component < result.size(); component++)
		result[component] = scale * value[component];
	return result;
}

double Cartesian_Dot(const Cartesian3& lhs, const Cartesian3& rhs)
{
	double result = 0.0;
	for(std::size_t component = 0; component < lhs.size(); component++)
		result += lhs[component] * rhs[component];
	return result;
}

double Cartesian_Norm(const Cartesian3& value)
{
	return sqrt(std::max(0.0, Cartesian_Dot(value, value)));
}

bool Cartesian_Finite(const Cartesian3& value)
{
	for(double component : value)
		if(!std::isfinite(component))
			return false;
	return true;
}

class HermiteBincountDepositor
{
  public:
	HermiteBincountDepositor(
		const Event& before,
		const Event& after,
		std::vector<BincountContribution>& contributions)
	: dt_sec_(In_Units(after.time - before.time, sec)),
	  contributions_(contributions)
	{
		for(std::size_t component = 0; component < 3; component++)
		{
			position_before_km_[component] = In_Units(before.position[component], km);
			position_after_km_[component] = In_Units(after.position[component], km);
			velocity_before_km_s_[component] = In_Units(before.velocity[component], km / sec);
			velocity_after_km_s_[component] = In_Units(after.velocity[component], km / sec);
		}
		if(std::isfinite(dt_sec_) && dt_sec_ > 0.0)
		{
			// Cache the cubic Hermite curve and its derivative in power-basis
			// form. This is algebraically identical to evaluating h00..h11,
			// but avoids rebuilding the basis at every recursion probe.
			for(std::size_t component = 0; component < 3; component++)
			{
				position_constant_km_[component] = position_before_km_[component];
				position_linear_km_[component] = dt_sec_ * velocity_before_km_s_[component];
				position_quadratic_km_[component] =
				    -3.0 * position_before_km_[component]
				    -2.0 * dt_sec_ * velocity_before_km_s_[component]
				    +3.0 * position_after_km_[component]
				    -dt_sec_ * velocity_after_km_s_[component];
				position_cubic_km_[component] =
				    2.0 * position_before_km_[component]
				    +dt_sec_ * velocity_before_km_s_[component]
				    -2.0 * position_after_km_[component]
				    +dt_sec_ * velocity_after_km_s_[component];
				velocity_constant_km_s_[component] = velocity_before_km_s_[component];
				velocity_linear_km_s_[component] =
				    2.0 * position_quadratic_km_[component] / dt_sec_;
				velocity_quadratic_km_s_[component] =
				    3.0 * position_cubic_km_[component] / dt_sec_;
			}
		}
	}

	void Deposit()
	{
		if(!std::isfinite(dt_sec_) || dt_sec_ <= 0.0
		   || !Cartesian_Finite(position_before_km_) || !Cartesian_Finite(position_after_km_)
		   || !Cartesian_Finite(velocity_before_km_s_) || !Cartesian_Finite(velocity_after_km_s_)
		   || !Cartesian_Finite(position_constant_km_) || !Cartesian_Finite(position_linear_km_)
		   || !Cartesian_Finite(position_quadratic_km_) || !Cartesian_Finite(position_cubic_km_)
		   || !Cartesian_Finite(velocity_constant_km_s_) || !Cartesian_Finite(velocity_linear_km_s_)
		   || !Cartesian_Finite(velocity_quadratic_km_s_))
			return;
		Deposit_Dense_Interval(0.0, 1.0, 0);
	}

  private:
	double dt_sec_;
	Cartesian3 position_before_km_;
	Cartesian3 position_after_km_;
	Cartesian3 velocity_before_km_s_;
	Cartesian3 velocity_after_km_s_;
	Cartesian3 position_constant_km_{};
	Cartesian3 position_linear_km_{};
	Cartesian3 position_quadratic_km_{};
	Cartesian3 position_cubic_km_{};
	Cartesian3 velocity_constant_km_s_{};
	Cartesian3 velocity_linear_km_s_{};
	Cartesian3 velocity_quadratic_km_s_{};
	std::vector<BincountContribution>& contributions_;

	Cartesian3 Position(double fraction) const
	{
		fraction = Clamp_Unit_Interval(fraction);
		Cartesian3 result;
		for(std::size_t component = 0; component < result.size(); component++)
		{
			result[component] = position_constant_km_[component]
			    + fraction * (position_linear_km_[component]
			    + fraction * (position_quadratic_km_[component]
			    + fraction * position_cubic_km_[component]));
		}
		return result;
	}

	Cartesian3 Velocity(double fraction) const
	{
		fraction = Clamp_Unit_Interval(fraction);
		Cartesian3 result;
		for(std::size_t component = 0; component < result.size(); component++)
		{
			result[component] = velocity_constant_km_s_[component]
			    + fraction * (velocity_linear_km_s_[component]
			    + fraction * velocity_quadratic_km_s_[component]);
		}
		return result;
	}

	double Position_Linearity_Error(double start, double end) const
	{
		const Cartesian3 position_start = Position(start);
		const Cartesian3 position_end = Position(end);
		double maximum_error = 0.0;
		const double probes[3] = {0.25, 0.5, 0.75};
		for(double probe : probes)
		{
			const double fraction = start + probe * (end - start);
			const Cartesian3 dense_position = Position(fraction);
			const Cartesian3 chord_position = Cartesian_Add(
			    Cartesian_Scale(1.0 - probe, position_start),
			    Cartesian_Scale(probe, position_end));
			maximum_error = std::max(
			    maximum_error,
			    Cartesian_Norm(Cartesian_Subtract(dense_position, chord_position)));
		}
		return maximum_error;
	}

	double Integrate_Speed_Squared(double start, double end) const
	{
		if(!(end > start))
			return 0.0;
		const double midpoint = 0.5 * (start + end);
		const double half_width = 0.5 * (end - start);
		const double offset = half_width * BINCOUNT_GAUSS_LEGENDRE_OFFSET;
		const Cartesian3 velocity_left = Velocity(midpoint - offset);
		const Cartesian3 velocity_midpoint = Velocity(midpoint);
		const Cartesian3 velocity_right = Velocity(midpoint + offset);
		const double integral_over_fraction = half_width * (
		    (5.0 / 9.0) * Cartesian_Dot(velocity_left, velocity_left)
		    + (8.0 / 9.0) * Cartesian_Dot(velocity_midpoint, velocity_midpoint)
		    + (5.0 / 9.0) * Cartesian_Dot(velocity_right, velocity_right));
		return dt_sec_ * integral_over_fraction;
	}

	void Append_Contribution(int bin, double start, double end)
	{
		if(bin < 0 || bin >= NUM_BINS || !(end > start))
			return;
		BincountContribution contribution;
		contribution.bin = bin;
		contribution.dt_sec = dt_sec_ * (end - start);
		contribution.v2dt_km2_per_sec = Integrate_Speed_Squared(start, end);
		if(!std::isfinite(contribution.dt_sec) || contribution.dt_sec <= 0.0
		   || !std::isfinite(contribution.v2dt_km2_per_sec)
		   || contribution.v2dt_km2_per_sec < 0.0)
			return;
		if(!contributions_.empty() && contributions_.back().bin == contribution.bin)
		{
			contributions_.back().dt_sec += contribution.dt_sec;
			contributions_.back().v2dt_km2_per_sec += contribution.v2dt_km2_per_sec;
		}
		else
			contributions_.push_back(contribution);
	}

	void Append_Chord_Piece(
		double dense_start,
		double dense_end,
		const Cartesian3& chord_start,
		const Cartesian3& chord_delta,
		double chord_piece_start,
		double chord_piece_end)
	{
		if(!(chord_piece_end > chord_piece_start))
			return;
		const double chord_midpoint = 0.5 * (chord_piece_start + chord_piece_end);
		const Cartesian3 position_midpoint = Cartesian_Add(
		    chord_start, Cartesian_Scale(chord_midpoint, chord_delta));
		const double radius_midpoint = Cartesian_Norm(position_midpoint);
		int bin = -1;
		if(std::isfinite(radius_midpoint)
		   && radius_midpoint >= 0.0 && radius_midpoint < BIN_MAX_KM)
			bin = static_cast<int>(radius_midpoint / BIN_WIDTH_KM);
		const double dense_piece_start =
		    dense_start + (dense_end - dense_start) * chord_piece_start;
		const double dense_piece_end =
		    dense_start + (dense_end - dense_start) * chord_piece_end;
		Append_Contribution(bin, dense_piece_start, dense_piece_end);
	}

	bool Chord_Boundary_Root(
		const Cartesian3& chord_start,
		const Cartesian3& chord_delta,
		double quadratic_a,
		double quadratic_b,
		double target_radius,
		double chord_fraction_start,
		double chord_fraction_end,
		double& root) const
	{
		const double quadratic_c =
		    Cartesian_Dot(chord_start, chord_start) - target_radius * target_radius;
		const double discriminant =
		    quadratic_b * quadratic_b - 4.0 * quadratic_a * quadratic_c;
		if(!std::isfinite(discriminant) || discriminant < 0.0 || quadratic_a <= 0.0)
			return false;
		const double sqrt_discriminant = sqrt(std::max(0.0, discriminant));
		const double roots[2] = {
		    (-quadratic_b - sqrt_discriminant) / (2.0 * quadratic_a),
		    (-quadratic_b + sqrt_discriminant) / (2.0 * quadratic_a)
		};
		for(double candidate : roots)
		{
			if(candidate > chord_fraction_start + 1.0e-12
			   && candidate < chord_fraction_end - 1.0e-12)
			{
				root = candidate;
				return true;
			}
		}
		return false;
	}

	void Deposit_Monotonic_Chord(
		double dense_start,
		double dense_end,
		const Cartesian3& chord_start,
		const Cartesian3& chord_delta,
		double chord_fraction_start,
		double chord_fraction_end)
	{
		const Cartesian3 position_start = Cartesian_Add(
		    chord_start, Cartesian_Scale(chord_fraction_start, chord_delta));
		const Cartesian3 position_end = Cartesian_Add(
		    chord_start, Cartesian_Scale(chord_fraction_end, chord_delta));
		const double radius_start = Cartesian_Norm(position_start);
		const double radius_end = Cartesian_Norm(position_end);
		const double radius_min = std::min(radius_start, radius_end);
		const double radius_max = std::max(radius_start, radius_end);
		const int first_boundary = std::max(
		    1, static_cast<int>(floor(radius_min / BIN_WIDTH_KM)) + 1);
		const int last_boundary = std::min(
		    NUM_BINS, static_cast<int>(floor(radius_max / BIN_WIDTH_KM)));
		const double quadratic_a = Cartesian_Dot(chord_delta, chord_delta);
		const double quadratic_b = 2.0 * Cartesian_Dot(chord_start, chord_delta);
		const int boundary_step = (radius_end >= radius_start) ? 1 : -1;
		int boundary = (boundary_step > 0) ? first_boundary : last_boundary;
		const int boundary_end = (boundary_step > 0) ? last_boundary : first_boundary;
		double previous_fraction = chord_fraction_start;
		// Radius is monotonic on this chord segment, so its shell crossings
		// already occur in radial-boundary order. Emit them directly instead
		// of allocating and sorting a temporary fraction vector.
		while((boundary_step > 0 && boundary <= boundary_end)
		      || (boundary_step < 0 && boundary >= boundary_end))
		{
			const double target_radius = static_cast<double>(boundary) * BIN_WIDTH_KM;
			if(target_radius > radius_min + 1.0e-10 * BIN_WIDTH_KM
			   && target_radius < radius_max - 1.0e-10 * BIN_WIDTH_KM)
			{
				double root = 0.0;
				if(Chord_Boundary_Root(
				       chord_start,
				       chord_delta,
				       quadratic_a,
				       quadratic_b,
				       target_radius,
				       chord_fraction_start,
				       chord_fraction_end,
				       root))
				{
					Append_Chord_Piece(
					    dense_start,
					    dense_end,
					    chord_start,
					    chord_delta,
					    previous_fraction,
					    root);
					previous_fraction = root;
				}
			}
			boundary += boundary_step;
		}
		Append_Chord_Piece(
		    dense_start,
		    dense_end,
		    chord_start,
		    chord_delta,
		    previous_fraction,
		    chord_fraction_end);
	}

	void Deposit_Chord(double start, double end)
	{
		const Cartesian3 chord_start = Position(start);
		const Cartesian3 chord_end = Position(end);
		const Cartesian3 chord_delta = Cartesian_Subtract(chord_end, chord_start);
		const double chord_length_sqr = Cartesian_Dot(chord_delta, chord_delta);
		if(!std::isfinite(chord_length_sqr))
			return;
		if(chord_length_sqr <= 0.0)
		{
			const double radius = Cartesian_Norm(chord_start);
			const int bin = (std::isfinite(radius) && radius >= 0.0 && radius < BIN_MAX_KM)
			              ? static_cast<int>(radius / BIN_WIDTH_KM)
			              : -1;
			Append_Contribution(bin, start, end);
			return;
		}

		const double closest_approach =
		    -Cartesian_Dot(chord_start, chord_delta) / chord_length_sqr;
		if(closest_approach > 1.0e-12 && closest_approach < 1.0 - 1.0e-12)
		{
			Deposit_Monotonic_Chord(
			    start,
			    end,
			    chord_start,
			    chord_delta,
			    0.0,
			    closest_approach);
			Deposit_Monotonic_Chord(
			    start,
			    end,
			    chord_start,
			    chord_delta,
			    closest_approach,
			    1.0);
		}
		else
			Deposit_Monotonic_Chord(
			    start, end, chord_start, chord_delta, 0.0, 1.0);
	}

	bool Try_Deposit_Single_Bin(double start, double end)
	{
		if(!(end > start))
			return true;
		const Cartesian3 control_start = Position(start);
		const Cartesian3 control_end = Position(end);
		const Cartesian3 velocity_start = Velocity(start);
		const Cartesian3 velocity_end = Velocity(end);
		const double interval_dt_sec = dt_sec_ * (end - start);
		const Cartesian3 control_inner_start = Cartesian_Add(
		    control_start, Cartesian_Scale(interval_dt_sec / 3.0, velocity_start));
		const Cartesian3 control_inner_end = Cartesian_Subtract(
		    control_end, Cartesian_Scale(interval_dt_sec / 3.0, velocity_end));
		const Cartesian3 controls[4] = {
		    control_start, control_inner_start, control_inner_end, control_end
		};
		const Cartesian3 midpoint_position = Position(0.5 * (start + end));
		const double midpoint_radius = Cartesian_Norm(midpoint_position);
		if(!std::isfinite(midpoint_radius)
		   || !Cartesian_Finite(control_start) || !Cartesian_Finite(control_inner_start)
		   || !Cartesian_Finite(control_inner_end) || !Cartesian_Finite(control_end))
			return false;

		// A cubic Hermite segment is the cubic Bezier curve defined by these
		// four controls. Its convex hull bounds the maximum norm, while the
		// projection onto any unit vector is bounded by the control-point
		// projections. Together these give a conservative proof that the
		// entire dense curve stays inside one radial bin.
		Cartesian3 radial_direction{};
		if(midpoint_radius > 0.0)
			radial_direction = Cartesian_Scale(1.0 / midpoint_radius, midpoint_position);

		double maximum_control_radius = 0.0;
		double minimum_radial_projection = std::numeric_limits<double>::infinity();
		for(const Cartesian3& control : controls)
		{
			maximum_control_radius =
			    std::max(maximum_control_radius, Cartesian_Norm(control));
			if(midpoint_radius > 0.0)
				minimum_radial_projection =
				    std::min(minimum_radial_projection, Cartesian_Dot(radial_direction, control));
		}
		const double margin =
		    BINCOUNT_DENSE_POSITION_TOLERANCE_KM + 1.0e-10 * BIN_WIDTH_KM;
		if(midpoint_radius >= BIN_MAX_KM)
			return midpoint_radius > 0.0
			    && minimum_radial_projection >= BIN_MAX_KM + margin;

		const int bin = static_cast<int>(midpoint_radius / BIN_WIDTH_KM);
		const double lower_radius = static_cast<double>(bin) * BIN_WIDTH_KM;
		const double upper_radius = static_cast<double>(bin + 1) * BIN_WIDTH_KM;
		if(maximum_control_radius >= upper_radius - margin)
			return false;
		if(lower_radius > 0.0
		   && (midpoint_radius <= 0.0
		       || minimum_radial_projection <= lower_radius + margin))
			return false;
		Append_Contribution(bin, start, end);
		return true;
	}

	void Deposit_Dense_Interval(double start, double end, int depth)
	{
		if(Try_Deposit_Single_Bin(start, end))
			return;
		const double linearity_error = Position_Linearity_Error(start, end);
		if(std::isfinite(linearity_error)
		   && linearity_error > BINCOUNT_DENSE_POSITION_TOLERANCE_KM
		   && depth < BINCOUNT_DENSE_MAX_RECURSION)
		{
			const double midpoint = 0.5 * (start + end);
			Deposit_Dense_Interval(start, midpoint, depth + 1);
			Deposit_Dense_Interval(midpoint, end, depth + 1);
			return;
		}
		Deposit_Chord(start, end);
	}
};
}

bool Find_First_Outward_Hermite_Radius_Crossing(
	const Event& before,
	const Event& after,
	double target_radius,
	Event& crossing)
{
	return Find_First_Outward_Hermite_Radius_Crossing_Impl(
	    before, after, target_radius, crossing);
}

void Compute_Bincount_Interval_Contributions(
	const Event& before,
	const Event& after,
	std::vector<BincountContribution>& contributions)
{
	contributions.clear();
	HermiteBincountDepositor(before, after, contributions).Deposit();
}

double BincountBinLowerKm(std::size_t bin)
{
	if(bin <= static_cast<std::size_t>(NUM_BINS))
		return static_cast<double>(bin) * BIN_WIDTH_KM;
	if(bin > TOTAL_BINS)
		return std::numeric_limits<double>::quiet_NaN();
	const std::size_t exterior_edge = bin - static_cast<std::size_t>(NUM_BINS);
	const double log_span = log(RADIAL_DOMAIN_MAX_KM / BIN_MAX_KM);
	return BIN_MAX_KM
	     * exp(log_span * static_cast<double>(exterior_edge)
	           / static_cast<double>(EXTERIOR_LOG_BINS));
}

double BincountBinUpperKm(std::size_t bin)
{
	return BincountBinLowerKm(bin + 1);
}

int BincountBinIndexKm(double radius_km)
{
	if(!std::isfinite(radius_km) || radius_km < 0.0
	   || radius_km >= RADIAL_DOMAIN_MAX_KM)
		return -1;
	if(radius_km < BIN_MAX_KM)
		return std::min(
		    NUM_BINS - 1,
		    std::max(0, static_cast<int>(floor(radius_km / BIN_WIDTH_KM))));

	const double log_span = log(RADIAL_DOMAIN_MAX_KM / BIN_MAX_KM);
	const double coordinate =
	    log(radius_km / BIN_MAX_KM) / log_span
	    * static_cast<double>(EXTERIOR_LOG_BINS);
	const std::size_t exterior_bin = std::min<std::size_t>(
	    EXTERIOR_LOG_BINS - 1,
	    static_cast<std::size_t>(std::max(0.0, floor(coordinate))));
	return NUM_BINS + static_cast<int>(exterior_bin);
}

bool Compute_Bound_Kepler_Exterior_Arc(
	const Event& outward_event,
	BoundKeplerExteriorArc& arc)
{
	Event inbound_event;
	double return_time = 0.0;
	double apoapsis_radius = 0.0;
	if(!Bound_Kepler_Return_At_Same_Radius(
	       outward_event, inbound_event, return_time, apoapsis_radius))
		return false;

	const double radius_km = In_Units(outward_event.Radius(), km);
	const double speed_km_s = In_Units(outward_event.Speed(), km / sec);
	const double angular_momentum_km2_s =
	    In_Units(outward_event.Angular_Momentum(), km * km / sec);
	const double mu_km3_s2 =
	    In_Units(G_Newton * mSun, km * km * km / (sec * sec));
	const double specific_energy_km2_s2 =
	    0.5 * speed_km_s * speed_km_s - mu_km3_s2 / radius_km;
	const double semi_major_axis_km =
	    -mu_km3_s2 / (2.0 * specific_energy_km2_s2);
	double eccentricity_squared =
	    1.0 + 2.0 * specific_energy_km2_s2
	        * angular_momentum_km2_s * angular_momentum_km2_s
	        / (mu_km3_s2 * mu_km3_s2);
	if(!std::isfinite(radius_km) || !std::isfinite(semi_major_axis_km)
	   || !std::isfinite(eccentricity_squared) || eccentricity_squared < 0.0)
		return false;
	const double eccentricity = sqrt(std::max(0.0, eccentricity_squared));
	const double mean_motion_s_inv = sqrt(
	    mu_km3_s2
	    / (semi_major_axis_km * semi_major_axis_km * semi_major_axis_km));
	if(!(eccentricity > 1.0e-14 && eccentricity < 1.0)
	   || !std::isfinite(mean_motion_s_inv) || mean_motion_s_inv <= 0.0)
		return false;

	const double apoapsis_km = In_Units(apoapsis_radius, km);
	if(!std::isfinite(apoapsis_km) || apoapsis_km < radius_km)
		return false;
	arc = BoundKeplerExteriorArc();
	arc.terminal_event = inbound_event;
	arc.elapsed_time_sec = In_Units(return_time, sec);
	arc.kepler_period_sec = 2.0 * M_PI / mean_motion_s_inv;
	arc.apoapsis_km = apoapsis_km;
	arc.outer_domain_removed =
	    apoapsis_km > RADIAL_DOMAIN_MAX_KM * (1.0 + 1.0e-12);
	const double integration_limit_km =
	    arc.outer_domain_removed ? RADIAL_DOMAIN_MAX_KM : apoapsis_km;
	const double pass_factor = arc.outer_domain_removed ? 1.0 : 2.0;
	if(arc.outer_domain_removed)
	{
		const double cos_e_start = Clamp_Cosine(
		    (1.0 - radius_km / semi_major_axis_km) / eccentricity);
		const double cos_e_terminal = Clamp_Cosine(
		    (1.0 - RADIAL_DOMAIN_MAX_KM / semi_major_axis_km) / eccentricity);
		const double e_start = acos(cos_e_start);
		const double e_terminal = acos(cos_e_terminal);
		const double sin_e_start =
		    sqrt(std::max(0.0, 1.0 - cos_e_start * cos_e_start));
		const double sin_e_terminal =
		    sqrt(std::max(0.0, 1.0 - cos_e_terminal * cos_e_terminal));
		arc.elapsed_time_sec =
		    ((e_terminal - eccentricity * sin_e_terminal)
		     - (e_start - eccentricity * sin_e_start))
		    / mean_motion_s_inv;

		const double radius = outward_event.Radius();
		const double mu = G_Newton * mSun;
		const libphysica::Vector radial_unit = outward_event.position / radius;
		const libphysica::Vector angular_momentum_vector =
		    outward_event.position.Cross(outward_event.velocity);
		const double angular_momentum = angular_momentum_vector.Norm();
		const double semi_latus_rectum =
		    angular_momentum * angular_momentum / mu;
		const libphysica::Vector eccentricity_vector =
		    outward_event.velocity.Cross(angular_momentum_vector) / mu
		    - radial_unit;
		const double eccentricity_natural = eccentricity_vector.Norm();
		if(!std::isfinite(angular_momentum) || angular_momentum <= 0.0
		   || !std::isfinite(semi_latus_rectum) || semi_latus_rectum <= 0.0
		   || !std::isfinite(eccentricity_natural)
		   || eccentricity_natural <= 0.0)
			return false;

		libphysica::Vector axis_x =
		    eccentricity_vector / eccentricity_natural;
		const libphysica::Vector axis_z =
		    angular_momentum_vector / angular_momentum;
		libphysica::Vector axis_y = axis_z.Cross(axis_x);
		const double axis_y_norm = axis_y.Norm();
		if(!std::isfinite(axis_y_norm) || axis_y_norm <= 0.0)
			return false;
		axis_y = axis_y / axis_y_norm;
		axis_x = axis_y.Cross(axis_z).Normalized();

		const double terminal_radius = RADIAL_DOMAIN_MAX_KM * km;
		const double cos_true_anomaly = Clamp_Cosine(
		    (semi_latus_rectum / terminal_radius - 1.0)
		    / eccentricity_natural);
		const double true_anomaly = acos(cos_true_anomaly);
		arc.terminal_event = outward_event;
		arc.terminal_event.time += arc.elapsed_time_sec * sec;
		arc.terminal_event.position =
		    terminal_radius
		    * (cos(true_anomaly) * axis_x + sin(true_anomaly) * axis_y);
		arc.terminal_event.velocity =
		    sqrt(mu / semi_latus_rectum)
		    * (-sin(true_anomaly) * axis_x
		       + (eccentricity_natural + cos(true_anomaly)) * axis_y);
		if(!std::isfinite(arc.elapsed_time_sec)
		   || arc.elapsed_time_sec <= 0.0
		   || !std::isfinite(arc.terminal_event.Radius())
		   || !std::isfinite(arc.terminal_event.Speed())
		   || Radial_Velocity(arc.terminal_event) <= 0.0)
			return false;
	}

	const int first_bin_index = BincountBinIndexKm(radius_km);
	if(first_bin_index < NUM_BINS)
		return false;
	for(std::size_t bin = static_cast<std::size_t>(first_bin_index);
	    bin < TOTAL_BINS; bin++)
	{
		const double lower_radius_km =
		    std::max(radius_km, BincountBinLowerKm(bin));
		const double upper_radius_km =
		    std::min(integration_limit_km, BincountBinUpperKm(bin));
		if(!(upper_radius_km > lower_radius_km))
			break;

		const double cos_e_lower = Clamp_Cosine(
		    (1.0 - lower_radius_km / semi_major_axis_km) / eccentricity);
		const double cos_e_upper = Clamp_Cosine(
		    (1.0 - upper_radius_km / semi_major_axis_km) / eccentricity);
		const double e_lower = acos(cos_e_lower);
		const double e_upper = acos(cos_e_upper);
		const double sin_e_lower = sqrt(std::max(0.0, 1.0 - cos_e_lower * cos_e_lower));
		const double sin_e_upper = sqrt(std::max(0.0, 1.0 - cos_e_upper * cos_e_upper));
		const double delta_mean_anomaly =
		    (e_upper - eccentricity * sin_e_upper)
		    - (e_lower - eccentricity * sin_e_lower);
		const double delta_v2_primitive =
		    (e_upper - e_lower)
		    + eccentricity * (sin_e_upper - sin_e_lower);
		arc.dt_hist[bin] =
		    pass_factor * delta_mean_anomaly / mean_motion_s_inv;
		arc.v2dt_hist[bin] =
		    pass_factor * mu_km3_s2
		    / (semi_major_axis_km * mean_motion_s_inv)
		    * delta_v2_primitive;
		if(!std::isfinite(arc.dt_hist[bin]) || arc.dt_hist[bin] < 0.0
		   || !std::isfinite(arc.v2dt_hist[bin]) || arc.v2dt_hist[bin] < 0.0)
			return false;
	}
	return std::isfinite(arc.elapsed_time_sec) && arc.elapsed_time_sec > 0.0
	    && std::isfinite(arc.kepler_period_sec) && arc.kepler_period_sec > 0.0;
}

const char* TrajectoryDiagnosticEventTypeKey(TrajectoryDiagnosticEventType type)
{
	switch(type)
	{
		case TrajectoryDiagnosticEventType::Capture: return "capture";
		case TrajectoryDiagnosticEventType::ScatterPre: return "scatter_pre";
		case TrajectoryDiagnosticEventType::ScatterPost: return "scatter_post";
		case TrajectoryDiagnosticEventType::CandidateUnbinding: return "candidate_unbinding";
		case TrajectoryDiagnosticEventType::Recapture: return "recapture";
		case TrajectoryDiagnosticEventType::SunExit: return "sun_exit";
		case TrajectoryDiagnosticEventType::SunEntry: return "sun_entry";
		case TrajectoryDiagnosticEventType::EscapeValidated: return "escape_validated";
		case TrajectoryDiagnosticEventType::Censored: return "censored";
		case TrajectoryDiagnosticEventType::NumericalFailure: return "numerical_failure";
		default: return "unknown";
	}
}

const char* TrajectoryNumericalFailureDetailKey(
	TrajectoryNumericalFailureDetail detail)
{
	switch(detail)
	{
		case TrajectoryNumericalFailureDetail::None:
			return "none";
		case TrajectoryNumericalFailureDetail::RK45ToleranceFailure:
			return "rk45_tolerance_failure";
		case TrajectoryNumericalFailureDetail::NonFinitePropagationState:
			return "non_finite_propagation_state";
		case TrajectoryNumericalFailureDetail::InvalidScatteringRate:
			return "invalid_scattering_rate";
		case TrajectoryNumericalFailureDetail::OpticalDepthRetryExhausted:
			return "optical_depth_retry_exhausted";
		case TrajectoryNumericalFailureDetail::UncapturedBoundMismatch:
			return "uncaptured_bound_mismatch";
		case TrajectoryNumericalFailureDetail::BoundaryCrossingLocalizationFailure:
			return "boundary_crossing_localization_failure";
		case TrajectoryNumericalFailureDetail::BoundaryDirectionMismatch:
			return "boundary_direction_mismatch";
		case TrajectoryNumericalFailureDetail::BoundaryEnergyMismatch:
			return "boundary_energy_mismatch";
		case TrajectoryNumericalFailureDetail::KeplerArcConstructionFailure:
			return "kepler_arc_construction_failure";
		default:
			return "unknown";
	}
}

double RK45PositionToleranceKm() { return In_Units(1.0 * km, km); }
double RK45VelocityToleranceKmPerSec() { return In_Units(1.0e-3 * km / sec, km / sec); }
double RK45PhaseTolerance() { return 1.0e-7; }
double RK45AbsoluteMaxStepSec() { return In_Units(RK45_Absolute_Max_Time_Step(), sec); }
double NormalModeMaxOpticalDepthStep() { return MAX_OPTICAL_DEPTH_STEP; }
double OpticalDepthRelativeTolerance() { return OPTICAL_DEPTH_RELATIVE_TOLERANCE; }
const char* BincountIntegrationScheme() { return "conservative-hermite-kepler-jupiter-log-v3"; }
double BincountDensePositionToleranceKm() { return BINCOUNT_DENSE_POSITION_TOLERANCE_KM; }
double SnapshotProgressPublishWallIntervalSeconds() { return SNAPSHOT_PUBLISH_WALL_INTERVAL_SEC; }
bool SnapshotProgressPublishDue(
    unsigned long int accepted_steps_since_publish, double wall_seconds_since_publish, bool force)
{
	return force || accepted_steps_since_publish >= SNAPSHOT_PUBLISH_STEP_INTERVAL
	       || (std::isfinite(wall_seconds_since_publish)
	           && wall_seconds_since_publish >= SNAPSHOT_PUBLISH_WALL_INTERVAL_SEC);
}

// 1. Result of one trajectory
bool TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason reason)
{
	switch(reason)
	{
		case TrajectoryTerminationReason::EnergyDriftEscape:
		case TrajectoryTerminationReason::NumericalFailure:
		case TrajectoryTerminationReason::NonFiniteState:
		case TrajectoryTerminationReason::SpeedLimit:
		case TrajectoryTerminationReason::WallTimeLimit:
		case TrajectoryTerminationReason::MaxFreeSteps:
		case TrajectoryTerminationReason::MaxScatterings:
		case TrajectoryTerminationReason::OuterDomainRemoval:
		case TrajectoryTerminationReason::Unknown:
			return true;
		default:
			return false;
	}
}

Trajectory_Result::Trajectory_Result(const Event& event_ini, const Event& event_final, unsigned long int nScat,
	                                 TrajectoryBincount bc, std::vector<TrajectoryDiagnosticEvent> events)
: initial_event(event_ini), final_event(event_final), number_of_scatterings(nScat), bincount(std::move(bc)),
  diagnostic_events(std::move(events))
{
}

bool Trajectory_Result::Particle_Reflected() const
{
	if(bincount.termination_reason != TrajectoryTerminationReason::OutwardEscape)
		return false;
	double r	= final_event.Radius();
	double vesc = sqrt(2 * G_Newton * mSun / r);
	return r > rSun && final_event.Speed() > vesc && number_of_scatterings > 0;
}

bool Trajectory_Result::Particle_Free() const
{
	return bincount.termination_reason == TrajectoryTerminationReason::OutwardEscape
	    && number_of_scatterings == 0;
}

bool Trajectory_Result::Particle_Captured(Solar_Model& solar_model) const
{
	double r	= final_event.Radius();
	double vesc = solar_model.Local_Escape_Speed(r);
	return final_event.Speed() < vesc;
}

void Trajectory_Result::Print_Summary(Solar_Model& solar_model, unsigned int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR
				  << "Trajectory result summary" << std::endl
				  << std::endl
				  << "Number of scatterings:\t" << number_of_scatterings << std::endl
				  << "Simulation time [days]:\t" << libphysica::Round(In_Units(final_event.time, day)) << std::endl
				  << "Final radius [rSun]:\t" << libphysica::Round(In_Units(final_event.Radius(), rSun)) << std::endl
				  << "Final speed [km/sec]:\t" << libphysica::Round(In_Units(final_event.Speed(), km / sec)) << std::endl
				  << "Free particle:\t\t[" << (Particle_Free() ? "x" : " ") << "]" << std::endl
				  << "Captured:\t\t[" << (Particle_Captured(solar_model) ? "x" : " ") << "]" << std::endl
				  << "Reflection:\t\t[" << (Particle_Reflected() ? "x" : " ") << "]";

		if(Particle_Reflected())
		{
			double u_i_sqr = initial_event.Asymptotic_Speed_Sqr(solar_model);
			double u_f_sqr = final_event.Asymptotic_Speed_Sqr(solar_model);
			// 检查渐近速度平方是否为正
			if(u_i_sqr > 0.0 && u_f_sqr > 0.0)
			{
				double u_i = sqrt(u_i_sqr);
				double u_f = sqrt(u_f_sqr);
				std::cout << "\t(ratio u_f/u_i = " << libphysica::Round(u_f / u_i) << ")" << std::endl;
			}
			else
				std::cout << "\t(asymptotic speeds unavailable)" << std::endl;
		}
		else
			std::cout << std::endl;
		std::cout << SEPARATOR << std::endl;
	}
}

// 2. Simulator
Trajectory_Simulator::Trajectory_Simulator(const Solar_Model& model, unsigned long int max_time_steps, unsigned long int max_scatterings, double max_distance)
: solar_model(model), has_previous_bincount_event(false),
  free_flight_reference_energy_eV(std::numeric_limits<double>::quiet_NaN()), current_ballistic_energy_drift_eV(0.0),
  current_physical_bound_state(false), terminate_on_capture(false), diagnostic_trace_enabled(false),
  diagnostic_event_index(0), diagnostic_scatter_index(0), diagnostic_step_index(0), last_scatter_target_index(-2),
  snapshot_recorder(nullptr), trajectory_in_progress(false), track_trajectory_wall_time(false), accumulated_snapshot_overhead_sec(0.0),
  steps_since_snapshot_publish(0),
  maximum_time_steps(max_time_steps), maximum_scatterings(max_scatterings), maximum_distance(max_distance), current_mpi_rank(0), current_trajectory_id(0)
{
	if(!std::isfinite(maximum_distance) || maximum_distance <= 0.0)
		throw std::invalid_argument("Trajectory_Simulator(): maximum distance must be finite and positive.");
	std::random_device rd;
	PRNG.seed(rd());
	rate_nuclei_cache.resize(solar_model.target_isotopes.size());
	bincount_contribution_cache.reserve(256);
}

void Trajectory_Simulator::Set_Snapshot_Recorder(SnapshotRecorder* recorder)
{
	snapshot_recorder = recorder;
}

void Trajectory_Simulator::Enable_Capture_Mode(bool enabled)
{
	terminate_on_capture = enabled;
}

void Trajectory_Simulator::Enable_Diagnostic_Trace(bool enabled)
{
	diagnostic_trace_enabled = enabled;
}

void Trajectory_Simulator::Restore_PRNG_State(const std::string& serialized_state)
{
	std::istringstream stream(serialized_state);
	stream >> PRNG;
	stream >> std::ws;
	if(!stream.eof())
		throw std::invalid_argument("Restore_PRNG_State(): invalid std::mt19937 state.");
}

std::string Trajectory_Simulator::Serialize_PRNG_State() const
{
	std::ostringstream stream;
	stream << PRNG;
	return stream.str();
}

bool Trajectory_Simulator::Trajectory_In_Progress() const
{
	return trajectory_in_progress;
}

unsigned long int Trajectory_Simulator::Current_Trajectory_ID() const
{
	return current_trajectory_id;
}

double Trajectory_Simulator::Current_Trajectory_Wall_Time_Seconds() const
{
	if(!trajectory_in_progress || !track_trajectory_wall_time)
		return 0.0;

	const double total_wall_time_sec = 1.0e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - current_trajectory_wall_start).count();
	return std::max(0.0, total_wall_time_sec - accumulated_snapshot_overhead_sec);
}

void Trajectory_Simulator::Accumulate_Snapshot_Overhead(const std::chrono::steady_clock::time_point& operation_start)
{
	if(!trajectory_in_progress || !track_trajectory_wall_time)
		return;
	accumulated_snapshot_overhead_sec += 1.0e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(
		std::chrono::steady_clock::now() - operation_start).count();
}

const TrajectoryBincount& Trajectory_Simulator::Current_Trajectory_Bincount() const
{
	return current_bincount;
}

void Trajectory_Simulator::Accumulate_Bincount_Interval(
	const Event& before,
	const Event& after,
	double simulated_time_sec)
{
	const double dt_sec = In_Units(after.time - before.time, sec);
	if(!std::isfinite(dt_sec) || dt_sec <= 0.0)
		return;
	if(!current_bincount.is_captured)
		return;
	Compute_Bincount_Interval_Contributions(before, after, bincount_contribution_cache);
	for(const BincountContribution& contribution : bincount_contribution_cache)
	{
		if(contribution.bin < 0
		   || static_cast<std::size_t>(contribution.bin) >= current_bincount.dt_hist.size())
			continue;
		current_bincount.dt_hist[contribution.bin] += contribution.dt_sec;
		current_bincount.v2dt_hist[contribution.bin] += contribution.v2dt_km2_per_sec;
	}

	const double r_before_km = In_Units(before.Radius(), km);
	double interior_start = 0.0;
	double interior_end = 0.0;
	double inside_fraction = 0.0;
	if(Solar_Interior_Fraction_Interval(before, after, before.Radius(), after.Radius(), interior_start, interior_end))
		inside_fraction = std::max(0.0, std::min(1.0, interior_end - interior_start));
	current_bincount.time_inside_sun_after_capture_sec += inside_fraction * dt_sec;
	current_bincount.time_outside_sun_after_capture_sec += (1.0 - inside_fraction) * dt_sec;
	const double max_radius_km = std::max(r_before_km, In_Units(after.Radius(), km));
	if(!std::isfinite(current_bincount.max_radius_after_capture_km)
	   || max_radius_km > current_bincount.max_radius_after_capture_km)
		current_bincount.max_radius_after_capture_km = max_radius_km;
}

void Trajectory_Simulator::Maybe_Publish_Snapshot_Progress(double simulated_time_sec, bool force)
{
	if(snapshot_recorder == nullptr)
		return;
	steps_since_snapshot_publish++;
	const auto now = std::chrono::steady_clock::now();
	const double wall_seconds_since_publish =
	    std::chrono::duration<double>(now - last_snapshot_publish_wall_time).count();
	if(!SnapshotProgressPublishDue(steps_since_snapshot_publish, wall_seconds_since_publish, force))
		return;
	steps_since_snapshot_publish = 0;
	last_snapshot_publish_wall_time = now;
	// current_bincount already holds exactly the accumulation the shared state
	// used to build step by step, so publishing it wholesale is equivalent.
	const auto snapshot_operation_start = now;
	snapshot_recorder->PublishCurrentTrajectoryProgress(
	    current_bincount.dt_hist, current_bincount.v2dt_hist, simulated_time_sec);
	Accumulate_Snapshot_Overhead(snapshot_operation_start);
}

void Trajectory_Simulator::Reset_Bincount_Anchor(const Event& event)
{
	if(terminate_on_capture)
		return;
	Copy_Event_Into(previous_bincount_event, event);
	has_previous_bincount_event = true;
}

double Trajectory_Simulator::Capture_Energy_eV(double radius, double speed, obscura::DM_Particle& DM)
{
	double vesc = solar_model.Local_Escape_Speed(radius);
	double E = 0.5 * DM.mass * (speed * speed - vesc * vesc);
	return In_Units(E, eV);
}

void Trajectory_Simulator::Record_Diagnostic_Event(
	TrajectoryDiagnosticEventType type,
	const Event& event,
	obscura::DM_Particle& DM,
	const char* target_species)
{
	if(!diagnostic_trace_enabled)
		return;
	const auto diagnostic_operation_start = std::chrono::steady_clock::now();
	TrajectoryDiagnosticEvent record;
	record.rank = static_cast<int>(current_mpi_rank);
	record.trajectory_id = current_trajectory_id;
	record.event_index = diagnostic_event_index++;
	record.scatter_index = diagnostic_scatter_index;
	record.step_index = diagnostic_step_index;
	record.event_type = type;
	record.t_s = In_Units(event.time, sec);
	record.r_km = In_Units(event.Radius(), km);
	record.vr_km_s = In_Units(Radial_Velocity(event), km / sec);
	record.speed_km_s = In_Units(event.Speed(), km / sec);
	for(size_t component = 0; component < 3; component++)
	{
		record.position_km[component] = In_Units(event.position[component], km);
		record.velocity_km_s[component] = In_Units(event.velocity[component], km / sec);
	}
	record.energy_eV = Capture_Energy_eV(event.Radius(), event.Speed(), DM);
	record.angular_momentum_km2_s = In_Units(event.Angular_Momentum(), km * km / sec);
	record.is_bound = record.energy_eV < 0.0 ? 1 : 0;
	record.inside_sun = event.Radius() < rSun ? 1 : 0;
	record.candidate_active = std::isfinite(current_bincount.t_final_unbinding_scatter) ? 1 : 0;
	std::strncpy(record.target_species, target_species == nullptr ? "" : target_species,
	             sizeof(record.target_species) - 1);
	record.target_species[sizeof(record.target_species) - 1] = '\0';
	record.ballistic_energy_drift_eV = current_ballistic_energy_drift_eV;
	current_diagnostic_events.push_back(record);
	Accumulate_Snapshot_Overhead(diagnostic_operation_start);
}

void Trajectory_Simulator::Record_Surface_Crossing_Events(
	const Event& before,
	const Event& after,
	obscura::DM_Particle& DM)
{
	if(!diagnostic_trace_enabled || !current_bincount.is_captured)
		return;
	double first = 0.0;
	double second = 0.0;
	if(!Surface_Crossing_Fractions(before, after, first, second))
		return;
	const double fractions[2] = {first, second};
	for(double fraction : fractions)
	{
		if(fraction <= 1.0e-12 || fraction > 1.0 + 1.0e-12)
			continue;
		Event crossing = Interpolate_Event(before, after, fraction);
		const bool entering = Radial_Velocity(crossing) < 0.0;
		Record_Diagnostic_Event(entering ? TrajectoryDiagnosticEventType::SunEntry
		                                 : TrajectoryDiagnosticEventType::SunExit,
		                        crossing, DM);
		current_diagnostic_events.back().inside_sun = entering ? 1 : 0;
	}
}

bool Trajectory_Simulator::Update_Capture_State(double radius, double speed, double time, obscura::DM_Particle& DM, bool allow_new_capture)
{
	if(terminate_on_capture && !allow_new_capture && !current_bincount.is_captured)
		return false;

	double E_eV = Capture_Energy_eV(radius, speed, DM);
	double dE_from_prev_eV = std::numeric_limits<double>::quiet_NaN();
	if(std::isfinite(previous_capture_energy_eV))
		dE_from_prev_eV = E_eV - previous_capture_energy_eV;
	previous_capture_energy_eV = E_eV;
	bool negative_energy = E_eV < 0.0;
	double t_now_sec = In_Units(time, sec);
	if(current_bincount.is_captured)
	{
		if(!std::isfinite(current_bincount.min_energy_after_capture_eV)
		   || E_eV < current_bincount.min_energy_after_capture_eV)
			current_bincount.min_energy_after_capture_eV = E_eV;
		const double radius_km = In_Units(radius, km);
		if(!std::isfinite(current_bincount.max_radius_after_capture_km)
		   || radius_km > current_bincount.max_radius_after_capture_km)
			current_bincount.max_radius_after_capture_km = radius_km;
	}

	if(!allow_new_capture)
	{
		if(current_bincount.is_captured && std::isfinite(free_flight_reference_energy_eV))
		{
			double drift_eV = std::fabs(E_eV - free_flight_reference_energy_eV);
			if(std::isfinite(drift_eV) && drift_eV > current_bincount.max_free_energy_drift_eV)
				current_bincount.max_free_energy_drift_eV = drift_eV;
			current_ballistic_energy_drift_eV = std::isfinite(drift_eV) ? drift_eV : 0.0;
			double scale_eV = std::max(std::fabs(free_flight_reference_energy_eV), FREE_ENERGY_DRIFT_REL_E_SCALE_EV);
			double drift_rel = drift_eV / scale_eV;
			if(std::isfinite(drift_rel) && drift_rel > current_bincount.max_free_energy_drift_rel)
				current_bincount.max_free_energy_drift_rel = drift_rel;
		}

		if(current_physical_bound_state && current_bincount.is_captured)
		{
			current_bincount.t_last_bound = t_now_sec;
		}
		return false;
	}

	bool was_bound = current_physical_bound_state;
	current_physical_bound_state = negative_energy;
	free_flight_reference_energy_eV = E_eV;
	current_ballistic_energy_drift_eV = 0.0;

	if(negative_energy)
	{
		if(!current_bincount.is_captured)
		{
			current_bincount.is_captured = true;
			current_bincount.t_capture = t_now_sec;
			current_bincount.t_last_bound = t_now_sec;
			current_bincount.r_first_negative_km = In_Units(radius, km);
			current_bincount.E_first_negative_eV = E_eV;
			current_bincount.dE_first_negative_from_prev_eV = dE_from_prev_eV;
			current_bincount.min_energy_after_capture_eV = E_eV;
			current_bincount.max_radius_after_capture_km = In_Units(radius, km);
			if(snapshot_recorder != nullptr)
			{
				const auto snapshot_operation_start = std::chrono::steady_clock::now();
				snapshot_recorder->MarkCurrentCaptured(true);
				Accumulate_Snapshot_Overhead(snapshot_operation_start);
				Maybe_Publish_Snapshot_Progress(t_now_sec, true);
			}
		}
		else if(!was_bound)
		{
			current_bincount.number_of_recaptures++;
			current_bincount.t_final_unbinding_scatter = std::numeric_limits<double>::quiet_NaN();
			current_bincount.r_final_unbinding_km = std::numeric_limits<double>::quiet_NaN();
			current_bincount.E_final_unbinding_eV = std::numeric_limits<double>::quiet_NaN();
		}
		current_bincount.t_last_bound = t_now_sec;
		return true;
	}

	if(was_bound && current_bincount.is_captured)
	{
		current_bincount.number_of_bound_to_unbound++;
		if(!std::isfinite(current_bincount.t_first_unbinding_scatter))
		{
			current_bincount.t_first_unbinding_scatter = t_now_sec;
			current_bincount.r_first_unbinding_km = In_Units(radius, km);
			current_bincount.E_first_unbinding_eV = E_eV;
		}
		current_bincount.t_final_unbinding_scatter = t_now_sec;
		current_bincount.r_final_unbinding_km = In_Units(radius, km);
		current_bincount.E_final_unbinding_eV = E_eV;
	}

	return false;
}

TrajectoryTerminationReason Trajectory_Simulator::Propagate_Freely(Event& current_event, obscura::DM_Particle& DM)
{
	if(!current_physical_bound_state
	   && current_event.Radius() >= maximum_distance
	   && Radial_Velocity(current_event) > 0.0)
	{
		const double boundary_energy_eV =
		    Capture_Energy_eV(
		        current_event.Radius(),
		        current_event.Speed(),
		        DM);
		const double energy_tolerance_eV =
		    Boundary_Energy_Tolerance_eV(
		        free_flight_reference_energy_eV,
		        DM);
		if(std::isfinite(boundary_energy_eV)
		   && boundary_energy_eV >= -energy_tolerance_eV)
			return TrajectoryTerminationReason::OutwardEscape;
	}

	auto abort_if_uncaptured_bound = [&](const Event& event)
	{
		if(current_bincount.is_captured)
			return false;

		const double energy_eV = Capture_Energy_eV(event.Radius(), event.Speed(), DM);
		if(std::isfinite(energy_eV) && energy_eV >= 0.0)
			return false;

		current_bincount.numerical_failure_detail =
		    TrajectoryNumericalFailureDetail::UncapturedBoundMismatch;
		current_bincount.failure_energy_after_step_eV = energy_eV;
		current_bincount.failure_reference_energy_eV =
		    free_flight_reference_energy_eV;
		current_event = event;
		return true;
	};

	// A trajectory can become physically bound only at a scattering, where
	// Update_Capture_State(..., allow_new_capture=true) records the capture.
	// Entering free flight already bound while still marked uncaptured is a
	// numerical-state mismatch; without this guard it can repeat bound Kepler
	// returns up to maximum_time_steps and stall an entire MPI batch.
	if(abort_if_uncaptured_bound(current_event))
		return TrajectoryTerminationReason::NumericalFailure;

	// Keep each free-propagation segment near t=0. Otherwise a captured particle
	// with a large absolute lifetime can make small RK45 increments disappear in
	// double precision, producing dt=0 and an apparent hang.
	double time_origin = current_event.time;
	Event local_initial_event = current_event;
	local_initial_event.time = 0.0;
	Free_Particle_Propagator particle_propagator(local_initial_event);
	double minus_log_xi = -log(Positive_Unit_Uniform(PRNG));
	TrajectoryTerminationReason outcome = TrajectoryTerminationReason::MaxFreeSteps;
	unsigned long int time_steps = 0;
	unsigned long int step_attempts = 0;
	unsigned int optical_depth_retries = 0;
	unsigned int boundary_refinement_retries = 0;
	const bool fast_capture_mode = terminate_on_capture;
	const double optical_depth_step_limit = fast_capture_mode ? CAPTURE_MODE_MAX_OPTICAL_DEPTH_STEP : MAX_OPTICAL_DEPTH_STEP;
	const std::size_t optical_depth_piece_target = fast_capture_mode ? CAPTURE_MODE_OPTICAL_DEPTH_PIECES : MAX_OPTICAL_DEPTH_PIECES;

	// Reused per-step Event buffers. Constructing or assigning an Event allocates,
	// because libphysica::Vector owns a std::vector and its assignment operator
	// takes and returns by value, so the loop below fills these in place.
	Event event_before;
	Event event_after;
	Event accepted_event;
	Event boundary_event;
	Event commit_absolute_event;

	auto to_absolute_event = [&](const Event& local_event)
	{
		Event absolute_event = local_event;
		absolute_event.time += time_origin;
		return absolute_event;
	};

	auto reset_propagator_at_absolute_event = [&](const Event& absolute_event)
	{
		time_origin = absolute_event.time;
		Event local_event = absolute_event;
		local_event.time = 0.0;
		particle_propagator = Free_Particle_Propagator(local_event);
	};

	auto abort_if_wall_time_exceeded = [&](const char* phase)
	{
		if(max_trajectory_wall_time_sec <= 0.0)
			return false;

		double traj_wall = Current_Trajectory_Wall_Time_Seconds();
		if(traj_wall <= max_trajectory_wall_time_sec)
			return false;

		std::cerr << "\nWarning in Propagate_Freely(): trajectory wall-time budget exceeded (rank "
		          << current_mpi_rank << ", traj " << current_trajectory_id
		          << ", phase=" << phase << ", step_attempts=" << step_attempts
		          << ", traj_wall=" << traj_wall << "s > " << max_trajectory_wall_time_sec
		          << "s). Aborting trajectory." << std::endl;
		current_event = to_absolute_event(particle_propagator.Event_In_3D());
		return true;
	};

	auto commit_accepted_event = [&](const Event& local_accepted_event)
	{
		Event& absolute_accepted_event = commit_absolute_event;
		Copy_Event_Into(absolute_accepted_event, local_accepted_event);
		absolute_accepted_event.time += time_origin;
		const double simulated_time_sec = In_Units(absolute_accepted_event.time, sec);
		diagnostic_step_index++;
		if(!terminate_on_capture)
		{
			if(has_previous_bincount_event)
			{
				Accumulate_Bincount_Interval(previous_bincount_event, absolute_accepted_event, simulated_time_sec);
				Record_Surface_Crossing_Events(previous_bincount_event, absolute_accepted_event, DM);
			}
			Copy_Event_Into(previous_bincount_event, absolute_accepted_event);
			has_previous_bincount_event = true;
		}

		Maybe_Publish_Snapshot_Progress(simulated_time_sec, false);
		if(current_bincount.is_captured)
			current_bincount.number_of_integrator_steps_after_capture++;
		return Update_Capture_State(absolute_accepted_event.Radius(), absolute_accepted_event.Speed(), absolute_accepted_event.time, DM, false);
	};

	while(time_steps < maximum_time_steps && outcome == TrajectoryTerminationReason::MaxFreeSteps)
	{
		step_attempts++;
		if(abort_if_wall_time_exceeded("before_step"))
			return TrajectoryTerminationReason::WallTimeLimit;

		particle_propagator.Fill_Event_In_3D(event_before);
		double r_before = particle_propagator.Current_Radius();
		double v_before = particle_propagator.Current_Speed();
		double rate_before = 0.0;
		if(r_before >= rSun)
		{
			const double step_cap = Free_Propagation_Time_Step_Cap(r_before, v_before, maximum_distance);
			particle_propagator.time_step = std::min(RK45_Sanitized_Time_Step(particle_propagator.time_step), step_cap);
			const double radial_velocity_before =
			    (r_before > 0.0) ? event_before.position.Dot(event_before.velocity) / r_before : 0.0;
			if(r_before < maximum_distance && radial_velocity_before > 0.0)
			{
				const double boundary_time_estimate =
				    (maximum_distance - r_before) / radial_velocity_before;
				if(std::isfinite(boundary_time_estimate)
				   && boundary_time_estimate > 0.0)
				{
					// Slightly overshoot the estimated first outward crossing so
					// the accepted step brackets it, while preventing one step
					// from also crossing the exterior turning point.
					const double boundary_step_cap =
					    std::max(1.05 * boundary_time_estimate, 1.0e-8 * sec);
					particle_propagator.time_step =
					    std::min(particle_propagator.time_step, boundary_step_cap);
				}
			}
			if(r_before > rSun && radial_velocity_before < 0.0)
			{
				const double surface_time_estimate = (r_before - rSun) / (-radial_velocity_before);
				if(std::isfinite(surface_time_estimate) && surface_time_estimate > 0.0)
					particle_propagator.time_step =
					    std::min(particle_propagator.time_step, std::max(surface_time_estimate, 1.0e-8 * sec));
			}
		}
		else
		{
			particle_propagator.time_step = RK45_Sanitized_Time_Step(particle_propagator.time_step);
			rate_before = solar_model.Total_DM_Scattering_Rate(DM, r_before, v_before);
			if(!std::isfinite(rate_before) || rate_before < 0.0)
				{
					std::cerr << "\nWarning in Propagate_Freely(): invalid pre-step scattering rate (rank "
					          << current_mpi_rank << ", traj " << current_trajectory_id
					          << ", rate=" << rate_before << "). Marking trajectory as numerically failed." << std::endl;
					current_event = to_absolute_event(particle_propagator.Event_In_3D());
					current_bincount.numerical_failure_detail =
					    TrajectoryNumericalFailureDetail::InvalidScatteringRate;
					return TrajectoryTerminationReason::NumericalFailure;
				}
			if(rate_before > 0.0)
			{
				particle_propagator.time_step = std::min(particle_propagator.time_step, optical_depth_step_limit / rate_before);
			}
		}

		const Free_Particle_Propagator::Scalar_State propagator_state_before = particle_propagator.Save_Scalar_State();
		double t_before = particle_propagator.Current_Time();
		bool rk_step_ok = particle_propagator.Runge_Kutta_45_Step(solar_model);
		double actual_dt = particle_propagator.Current_Time() - t_before;
		double r_after = particle_propagator.Current_Radius();
		double v_after = particle_propagator.Current_Speed();
		particle_propagator.Fill_Event_In_3D(event_after);

		if(!rk_step_ok)
		{
				std::cerr << "\nWarning in Propagate_Freely(): RK45 step failed tolerance after retry limits (rank "
				          << current_mpi_rank << ", traj " << current_trajectory_id
				          << "). Marking trajectory as numerically failed." << std::endl;
			current_event = to_absolute_event(particle_propagator.Event_In_3D());
			current_bincount.numerical_failure_detail =
			    TrajectoryNumericalFailureDetail::RK45ToleranceFailure;
			current_bincount.failure_attempted_step_s =
			    In_Units(propagator_state_before.time_step, sec);
			current_bincount.failure_accepted_step_s =
			    In_Units(actual_dt, sec);
			return TrajectoryTerminationReason::NumericalFailure;
		}

			if(!std::isfinite(r_after) || !std::isfinite(v_after) || !std::isfinite(actual_dt) || actual_dt <= 0.0)
			{
				// 出现 NaN/Inf，轨迹已无法继续。中止这条轨迹，避免 rank 被单条坏轨迹卡住。
				std::cerr << "\nWarning in Propagate_Freely(): non-finite state (rank " << current_mpi_rank
				          << ", traj " << current_trajectory_id << ", r=" << r_after
				          << ", v=" << v_after << ", dt=" << actual_dt << "). Aborting trajectory." << std::endl;
			current_event = to_absolute_event(particle_propagator.Event_In_3D());
			current_bincount.numerical_failure_detail =
			    TrajectoryNumericalFailureDetail::NonFinitePropagationState;
			current_bincount.failure_attempted_step_s =
			    In_Units(propagator_state_before.time_step, sec);
			current_bincount.failure_accepted_step_s =
			    In_Units(actual_dt, sec);
			return TrajectoryTerminationReason::NonFiniteState;
		}

		if(v_after > v_max)
			{
				std::cerr << "\nWarning in Propagate_Freely(): DM speed exceeds the maximum of v_max = " << v_max << std::endl
						  << "\tAbort simulation." << std::endl;
				current_event = to_absolute_event(particle_propagator.Event_In_3D());
				return TrajectoryTerminationReason::SpeedLimit;
			}

		if(abort_if_wall_time_exceeded("after_rk45_step"))
		{
			commit_accepted_event(event_after);
			current_event = to_absolute_event(event_after);
			return TrajectoryTerminationReason::WallTimeLimit;
		}

		// Check for scatterings and reflection
		bool scattering = false;
		bool reflection = false;
		const double minus_log_xi_before_step = minus_log_xi;
		Copy_Event_Into(accepted_event, event_after);
		double interior_start = 0.0;
		double interior_end = 0.0;
		if(Solar_Interior_Fraction_Interval(event_before, event_after, r_before, r_after, interior_start, interior_end))
		{
			if(r_after < rSun && v_after < 0.0)
			{
				std::cerr << "Warning: Negative velocity detected (v = " << v_after << ") at r = " << r_after << ", skipping scattering calculation." << std::endl;
				break;
			}
			std::array<OpticalDepthPiece, MAX_OPTICAL_DEPTH_PIECES> optical_pieces;
			std::size_t optical_piece_count = 0;
			double tau_two_piece = 0.0;
			double delta_tau = 0.0;
			bool invalid_rate = false;
			const bool use_cached_start_rate = r_before < rSun && interior_start <= 1.0e-12;
				Build_Optical_Depth_Pieces(event_before,
				                           event_after,
				                           interior_start,
				                           interior_end,
				                           actual_dt,
				                           solar_model,
				                           DM,
				                           use_cached_start_rate,
				                           rate_before,
				                           optical_depth_piece_target,
				                           optical_pieces,
				                           optical_piece_count,
				                           tau_two_piece,
				                           delta_tau,
				                           invalid_rate);
			if(invalid_rate || !std::isfinite(delta_tau) || delta_tau < 0.0)
				{
					std::cerr << "\nWarning in Propagate_Freely(): invalid scattering rate (rank "
					          << current_mpi_rank << ", traj " << current_trajectory_id
					          << "). Marking trajectory as numerically failed." << std::endl;
					current_event = to_absolute_event(event_after);
					current_bincount.numerical_failure_detail =
					    TrajectoryNumericalFailureDetail::InvalidScatteringRate;
					return TrajectoryTerminationReason::NumericalFailure;
				}

				const double tau_error_scale = std::max(delta_tau, OPTICAL_DEPTH_ABSOLUTE_TOLERANCE);
				const double tau_relative_error = std::fabs(delta_tau - tau_two_piece) / tau_error_scale;
				const bool reject_for_optical_depth = delta_tau > optical_depth_step_limit;
				const bool reject_for_tau_accuracy = !fast_capture_mode && tau_relative_error > OPTICAL_DEPTH_RELATIVE_TOLERANCE;
				if(reject_for_optical_depth || reject_for_tau_accuracy)
				{
					optical_depth_retries++;
					if(optical_depth_retries > MAX_OPTICAL_DEPTH_RETRIES)
					{
							std::cerr << "\nWarning in Propagate_Freely(): optical-depth step rejected more than "
							          << MAX_OPTICAL_DEPTH_RETRIES << " consecutive times (rank "
							          << current_mpi_rank << ", traj " << current_trajectory_id
							          << ", delta_tau=" << delta_tau
							          << ", relative_error=" << tau_relative_error
							          << "). Marking trajectory as numerically failed." << std::endl;
							current_event = to_absolute_event(event_before);
							current_bincount.numerical_failure_detail =
							    TrajectoryNumericalFailureDetail::OpticalDepthRetryExhausted;
							return TrajectoryTerminationReason::NumericalFailure;
						}

					particle_propagator.Restore_Scalar_State(propagator_state_before);
					double factor = RK45_MAX_STEP_FACTOR;
					if(delta_tau > optical_depth_step_limit && delta_tau > 0.0)
						factor = std::min(factor, 0.8 * optical_depth_step_limit / delta_tau);
					if(reject_for_tau_accuracy && tau_relative_error > 0.0)
						factor = std::min(factor, 0.8 * std::sqrt(OPTICAL_DEPTH_RELATIVE_TOLERANCE / tau_relative_error));
					factor = std::max(RK45_MIN_STEP_FACTOR, std::min(factor, 1.0));
					particle_propagator.time_step = RK45_Sanitized_Time_Step(factor * actual_dt);
					continue;
				}

			if(delta_tau >= minus_log_xi)
			{
				double remaining_tau = minus_log_xi;
				for(std::size_t piece_index = 0; piece_index < optical_piece_count; piece_index++)
				{
					const OpticalDepthPiece& piece = optical_pieces[piece_index];
					const double piece_dt = actual_dt * (piece.end - piece.start);
					const double piece_tau = 0.5 * (piece.rate_start + piece.rate_end) * piece_dt;
					if(piece_tau >= remaining_tau)
					{
						const double local_fraction =
						    Solve_Linear_Rate_Optical_Depth_Fraction(piece.rate_start, piece.rate_end, piece_dt, remaining_tau);
						const double fraction = piece.start + local_fraction * (piece.end - piece.start);
						accepted_event = Interpolate_Event(event_before, event_after, fraction);
						scattering = true;
						outcome = TrajectoryTerminationReason::Scatter;
						break;
					}
					remaining_tau -= piece_tau;
				}
				if(!scattering)
				{
					accepted_event = Interpolate_Event(event_before, event_after, interior_end);
					scattering = true;
					outcome = TrajectoryTerminationReason::Scatter;
				}
			}
			else
				minus_log_xi -= delta_tau;
		}

		if(!scattering)
		{
			Copy_Event_Into(boundary_event, event_after);
			const bool endpoint_brackets_matching_boundary =
			    r_before <= maximum_distance * (1.0 + 1.0e-12)
			    && r_after >= maximum_distance;
			if(endpoint_brackets_matching_boundary)
			{
				if(!Find_First_Outward_Hermite_Radius_Crossing(
				       event_before,
				       event_after,
				       maximum_distance,
				       boundary_event))
				{
					boundary_refinement_retries++;
					if(boundary_refinement_retries
					   <= MAX_BOUNDARY_REFINEMENT_RETRIES)
					{
						particle_propagator.Restore_Scalar_State(
						    propagator_state_before);
						particle_propagator.time_step =
						    RK45_Sanitized_Time_Step(
						        std::max(
						            0.5 * actual_dt,
						            1.0e-8 * sec));
						minus_log_xi = minus_log_xi_before_step;
						continue;
					}
					current_bincount.numerical_failure_detail =
					    TrajectoryNumericalFailureDetail::
					        BoundaryCrossingLocalizationFailure;
					current_bincount.failure_energy_before_step_eV =
					    Capture_Energy_eV(
					        event_before.Radius(),
					        event_before.Speed(),
					        DM);
					current_bincount.failure_energy_after_step_eV =
					    Capture_Energy_eV(
					        event_after.Radius(),
					        event_after.Speed(),
					        DM);
					current_bincount.failure_reference_energy_eV =
					    free_flight_reference_energy_eV;
					current_bincount.failure_attempted_step_s =
					    In_Units(
					        propagator_state_before.time_step,
					        sec);
					current_bincount.failure_accepted_step_s =
					    In_Units(actual_dt, sec);
					current_event = to_absolute_event(event_after);
					return TrajectoryTerminationReason::NumericalFailure;
				}
				const double boundary_energy_eV =
				    Capture_Energy_eV(
				        boundary_event.Radius(),
				        boundary_event.Speed(),
				        DM);
				const double reference_energy_eV =
				    free_flight_reference_energy_eV;
				const double energy_tolerance_eV =
				    Boundary_Energy_Tolerance_eV(
				        reference_energy_eV,
				        DM);
				const bool boundary_energy_consistent =
				    std::isfinite(boundary_energy_eV)
				    && std::isfinite(reference_energy_eV)
				    && (current_physical_bound_state
				        ? boundary_energy_eV
				              <= energy_tolerance_eV
				        : boundary_energy_eV
				              >= -energy_tolerance_eV);
				if(!boundary_energy_consistent)
				{
					boundary_refinement_retries++;
					if(boundary_refinement_retries
					   <= MAX_BOUNDARY_REFINEMENT_RETRIES)
					{
						particle_propagator.Restore_Scalar_State(
						    propagator_state_before);
						particle_propagator.time_step =
						    RK45_Sanitized_Time_Step(
						        std::max(
						            0.5 * actual_dt,
						            1.0e-8 * sec));
						minus_log_xi = minus_log_xi_before_step;
						continue;
					}
					current_bincount.numerical_failure_detail =
					    TrajectoryNumericalFailureDetail::
					        BoundaryEnergyMismatch;
					current_bincount.failure_energy_before_step_eV =
					    Capture_Energy_eV(
					        event_before.Radius(),
					        event_before.Speed(),
					        DM);
					current_bincount.failure_energy_after_step_eV =
					    Capture_Energy_eV(
					        event_after.Radius(),
					        event_after.Speed(),
					        DM);
					current_bincount.failure_energy_at_boundary_eV =
					    boundary_energy_eV;
					current_bincount.failure_reference_energy_eV =
					    reference_energy_eV;
					current_bincount.failure_boundary_vr_km_s =
					    In_Units(
					        Radial_Velocity(boundary_event),
					        km / sec);
					current_bincount.failure_attempted_step_s =
					    In_Units(
					        propagator_state_before.time_step,
					        sec);
					current_bincount.failure_accepted_step_s =
					    In_Units(actual_dt, sec);
					std::cerr
					    << "\nWarning in Propagate_Freely(): "
					    << "boundary energy refinement exhausted (rank "
					    << current_mpi_rank
					    << ", traj "
					    << current_trajectory_id
					    << ", E_ref="
					    << reference_energy_eV
					    << " eV, E_boundary="
					    << boundary_energy_eV
					    << " eV, tolerance="
					    << energy_tolerance_eV
					    << " eV, vr_boundary="
					    << current_bincount.failure_boundary_vr_km_s
					    << " km/s)."
					    << std::endl;
					current_event =
					    to_absolute_event(boundary_event);
					return TrajectoryTerminationReason::NumericalFailure;
				}
				boundary_refinement_retries = 0;

				if(!current_physical_bound_state)
				{
					accepted_event = boundary_event;
					reflection = true;
					outcome =
					    TrajectoryTerminationReason::OutwardEscape;
				}
				else
				{
				const Event absolute_boundary_event = to_absolute_event(boundary_event);
				if(abort_if_uncaptured_bound(absolute_boundary_event))
					return TrajectoryTerminationReason::NumericalFailure;

				BoundKeplerExteriorArc kepler_arc;
				if(Compute_Bound_Kepler_Exterior_Arc(boundary_event, kepler_arc))
				{
					time_steps++;
					optical_depth_retries = 0;
					commit_accepted_event(boundary_event);
					if(current_bincount.is_captured)
					{
						current_bincount.number_of_bound_exterior_arcs++;
						if(current_bincount.number_of_bound_exterior_arcs == 1)
						{
							current_bincount.first_bound_exit_kepler_period_sec =
							    kepler_arc.kepler_period_sec;
							current_bincount.first_bound_exit_exterior_time_sec =
							    kepler_arc.elapsed_time_sec;
						}
						current_bincount.last_bound_exit_kepler_period_sec =
						    kepler_arc.kepler_period_sec;
						current_bincount.last_bound_exit_exterior_time_sec =
						    kepler_arc.elapsed_time_sec;
						if(!std::isfinite(current_bincount.max_bound_exit_kepler_period_sec)
						   || kepler_arc.kepler_period_sec
						      > current_bincount.max_bound_exit_kepler_period_sec)
							current_bincount.max_bound_exit_kepler_period_sec =
							    kepler_arc.kepler_period_sec;
						if(!std::isfinite(current_bincount.max_bound_exit_exterior_time_sec)
						   || kepler_arc.elapsed_time_sec
						      > current_bincount.max_bound_exit_exterior_time_sec)
							current_bincount.max_bound_exit_exterior_time_sec =
							    kepler_arc.elapsed_time_sec;
						for(std::size_t bin = NUM_BINS; bin < TOTAL_BINS; bin++)
						{
							current_bincount.dt_hist[bin] += kepler_arc.dt_hist[bin];
							current_bincount.v2dt_hist[bin] += kepler_arc.v2dt_hist[bin];
						}
						current_bincount.time_outside_sun_after_capture_sec +=
						    kepler_arc.elapsed_time_sec;
					}
					current_event = to_absolute_event(kepler_arc.terminal_event);
					if(kepler_arc.outer_domain_removed)
					{
						if(current_bincount.is_captured)
						{
							if(!std::isfinite(current_bincount.max_radius_after_capture_km)
							   || RADIAL_DOMAIN_MAX_KM
							      > current_bincount.max_radius_after_capture_km)
								current_bincount.max_radius_after_capture_km =
								    RADIAL_DOMAIN_MAX_KM;
							current_bincount.t_last_bound = In_Units(current_event.time, sec);
						}
						Maybe_Publish_Snapshot_Progress(
						    In_Units(current_event.time, sec), true);
						return TrajectoryTerminationReason::OuterDomainRemoval;
					}
					if(current_bincount.is_captured)
					{
						const double apoapsis_km = kepler_arc.apoapsis_km;
						if(!std::isfinite(current_bincount.max_radius_after_capture_km)
						   || apoapsis_km > current_bincount.max_radius_after_capture_km)
							current_bincount.max_radius_after_capture_km = apoapsis_km;
						current_bincount.t_last_bound = In_Units(current_event.time, sec);
					}
					Reset_Bincount_Anchor(current_event);
					// Synchronize the in-progress state after the large analytic
					// jump. Snapshot files themselves remain on fixed wall cadence.
					Maybe_Publish_Snapshot_Progress(In_Units(current_event.time, sec), true);
					reset_propagator_at_absolute_event(current_event);
					continue;
				}
				current_bincount.numerical_failure_detail =
				    TrajectoryNumericalFailureDetail::
				        KeplerArcConstructionFailure;
				current_bincount.failure_energy_at_boundary_eV =
				    boundary_energy_eV;
				current_bincount.failure_reference_energy_eV =
				    reference_energy_eV;
				current_bincount.failure_boundary_vr_km_s =
				    In_Units(
				        Radial_Velocity(boundary_event),
				        km / sec);
				current_event = absolute_boundary_event;
				return TrajectoryTerminationReason::NumericalFailure;
				}
			}
		}

		time_steps++;
		optical_depth_retries = 0;
		const bool captured_now = commit_accepted_event(accepted_event);
		Copy_Event_Into(current_event, accepted_event);
		current_event.time += time_origin;
		if(abort_if_uncaptured_bound(current_event))
			return TrajectoryTerminationReason::NumericalFailure;
		if(terminate_on_capture && captured_now)
			return TrajectoryTerminationReason::CaptureMode;

		if(reflection || scattering)
			break;

	}

	return outcome;
}


int Trajectory_Simulator::Sample_Target(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(r > rSun)
	{
		throw std::runtime_error("Sample_Target(): r > rSun.");
	}
	else
	{
		// C: 复用预分配的 rate_nuclei_cache，避免每次散射事件堆分配
		for(unsigned int i = 0; i < solar_model.target_isotopes.size(); i++)
		{
			rate_nuclei_cache[i] = solar_model.DM_Scattering_Rate_Nucleus(DM, r, DM_speed, i);
			if(!std::isfinite(rate_nuclei_cache[i]) || rate_nuclei_cache[i] < 0.0)
				throw std::runtime_error("Sample_Target(): nucleus scattering rate is negative or non-finite.");
		}
		double rate_electron = solar_model.DM_Scattering_Rate_Electron(DM, r, DM_speed);
		if(!std::isfinite(rate_electron) || rate_electron < 0.0)
			throw std::runtime_error("Sample_Target(): electron scattering rate is negative or non-finite.");
		double total_rate	 = std::accumulate(rate_nuclei_cache.begin(), rate_nuclei_cache.end(), rate_electron);
		if(!std::isfinite(total_rate) || total_rate <= 0.0)
			throw std::runtime_error("Sample_Target(): total scattering rate is non-positive or non-finite.");

		double xi = libphysica::Sample_Uniform(PRNG);
		// Electron
		double sum = rate_electron / total_rate;
		if(sum > xi)
			return -1;
		// Nuclei
		for(unsigned int i = 0; i < solar_model.target_isotopes.size(); i++)
		{
			sum += rate_nuclei_cache[i] / total_rate;
			if(sum > xi || i == solar_model.target_isotopes.size() - 1)
				return i;
		}
		throw std::runtime_error("Sample_Target(): no target could be sampled.");
	}
}

libphysica::Vector Trajectory_Simulator::Sample_Target_Velocity(double temperature, double target_mass, const libphysica::Vector& vel_DM)
{
	// Sampling algorithm taken from Romano & Walsh, "An improved target velocity sampling algorithm for free gas elastic scattering"
	if(!std::isfinite(temperature) || temperature <= 0.0
	   || !std::isfinite(target_mass) || target_mass <= 0.0)
		throw std::runtime_error("Sample_Target_Velocity(): temperature and target mass must be finite and positive.");
	double kappa = sqrt(target_mass / 2.0 / temperature);
	double vDM	 = vel_DM.Norm();
	if(!std::isfinite(kappa) || kappa <= 0.0 || !std::isfinite(vDM) || vDM <= 0.0)
		throw std::runtime_error("Sample_Target_Velocity(): DM speed is non-positive or non-finite.");
	// 1. Sample target speed vT and mu = cos alpha
	double y  = kappa * vDM;
	double x  = y;
	double mu = 1.0;
	unsigned long int rejection_attempts = 0;
	auto relative_speed_ratio = [&]() {
		const double x_sqr = x * x;
		const double y_sqr = y * y;
		const double relative_speed_sqr_raw = x_sqr + y_sqr - 2.0 * x * y * mu;
		const double denominator = x + y;
		const double roundoff_scale = std::max(1.0, x_sqr + y_sqr + 2.0 * std::fabs(x * y));
		if(!std::isfinite(x) || x < 0.0 || !std::isfinite(y) || y <= 0.0
		   || !std::isfinite(mu) || mu < -1.0 || mu > 1.0
		   || !std::isfinite(relative_speed_sqr_raw)
		   || relative_speed_sqr_raw < -1.0e-12 * roundoff_scale
		   || !std::isfinite(denominator) || denominator <= 0.0)
			throw std::runtime_error("Sample_Target_Velocity(): rejection state is invalid.");
		const double relative_speed_sqr = std::max(0.0, relative_speed_sqr_raw);
		const double ratio = sqrt(relative_speed_sqr) / denominator;
		if(!std::isfinite(ratio) || ratio < 0.0 || ratio > 1.0 + 1.0e-12)
			throw std::runtime_error("Sample_Target_Velocity(): rejection ratio is invalid.");
		return std::min(1.0, ratio);
	};
	while(relative_speed_ratio() < libphysica::Sample_Uniform(PRNG, 0.0, 1.0))
	{
		rejection_attempts++;
		if(rejection_attempts > TARGET_VELOCITY_MAX_REJECTION_ATTEMPTS)
			throw std::runtime_error("Sample_Target_Velocity(): rejection sampling exceeded maximum attempts.");
		mu			= libphysica::Sample_Uniform(PRNG, -1.0, 1.0);
		double xi_1 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
		if(xi_1 < 2.0 / (sqrt(M_PI) * y + 2.0))
		{
			double xi_2 = Positive_Unit_Uniform(PRNG);
			double xi_3 = Positive_Unit_Uniform(PRNG);
			double z	= -log(xi_2) - log(xi_3);
			x			= sqrt(z);
		}
		else
		{
			double xi_2 = Positive_Unit_Uniform(PRNG);
			double xi_3 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
			double xi_4 = Positive_Unit_Uniform(PRNG);
			double z	= -log(xi_2) - pow(cos(M_PI / 2.0 * xi_3), 2.0) * log(xi_4);
			x			= sqrt(z);
		}
	}
	// 2. Construct the target velocity vel_T
	double vT		 = x / kappa;
	double cos_theta = mu;
	double phi		 = libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	libphysica::Vector unit_vector_T = Unit_Vector_At_Angle_From_Axis(cos_theta, phi, vel_DM);

	libphysica::Vector velocity = vT * unit_vector_T;
	if(!std::isfinite(velocity.Norm()))
		throw std::runtime_error("Sample_Target_Velocity(): sampled target velocity is non-finite.");
	return velocity;
}

libphysica::Vector Trajectory_Simulator::New_DM_Velocity(double cos_scattering_angle, double DM_mass, double target_mass, libphysica::Vector& vel_DM, libphysica::Vector& vel_target)
{
	if(!std::isfinite(cos_scattering_angle)
	   || cos_scattering_angle < -1.0 - 1.0e-12 || cos_scattering_angle > 1.0 + 1.0e-12
	   || !std::isfinite(DM_mass) || DM_mass <= 0.0
	   || !std::isfinite(target_mass) || target_mass <= 0.0
	   || !std::isfinite(vel_DM.Norm()) || !std::isfinite(vel_target.Norm()))
		throw std::runtime_error("New_DM_Velocity(): masses, velocities, or scattering angle are invalid.");
	cos_scattering_angle = Clamp_Cosine(cos_scattering_angle);
	// Construction of n, the unit vector pointing into the direction of vfinal.
	double phi			 = libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	libphysica::Vector n = Unit_Vector_At_Angle_From_Axis(cos_scattering_angle, phi, vel_DM);

	double relative_speed = (vel_target - vel_DM).Norm();

	libphysica::Vector velocity = target_mass * relative_speed / (target_mass + DM_mass) * n + (DM_mass * vel_DM + target_mass * vel_target) / (target_mass + DM_mass);
	if(!std::isfinite(velocity.Norm()))
		throw std::runtime_error("New_DM_Velocity(): final DM velocity is non-finite.");
	return velocity;
}

void Trajectory_Simulator::Scatter(Event& current_event, obscura::DM_Particle& DM)
{
	double r = current_event.Radius();
	double v = current_event.Speed();
	// 1. Find target properties.
	int target_index = Sample_Target(DM, r, v);
	last_scatter_target_index = target_index;

	double target_mass;
	if(target_index == -1)
		target_mass = mElectron;
	else
		target_mass = solar_model.target_isotopes[target_index].mass;

	libphysica::Vector vel_target = Sample_Target_Velocity(solar_model.Temperature(r), target_mass, current_event.velocity);

	// 2. Sample the scattering angle
	double cos_alpha = (target_index == -1) ? DM.Sample_Scattering_Angle_Electron(PRNG, v, r) : DM.Sample_Scattering_Angle_Nucleus(PRNG, solar_model.target_isotopes[target_index], v, r);

	// 3. Construct the final DM velocity
	current_event.velocity = New_DM_Velocity(cos_alpha, DM.mass, target_mass, current_event.velocity, vel_target);
}

void Trajectory_Simulator::Fix_PRNG_Seed(unsigned int fixed_seed)
{
	PRNG.seed(fixed_seed);
}

Trajectory_Result Trajectory_Simulator::Simulate(const Event& initial_condition, obscura::DM_Particle& DM, unsigned int mpi_rank)
{
	if(!std::isfinite(initial_condition.time) || initial_condition.time < 0.0
	   || !std::isfinite(initial_condition.Radius()) || initial_condition.Radius() <= 0.0
	   || !std::isfinite(initial_condition.Speed()) || initial_condition.Speed() <= 0.0
	   || !std::isfinite(DM.mass) || DM.mass <= 0.0
	   || !std::isfinite(max_trajectory_wall_time_sec) || max_trajectory_wall_time_sec < 0.0)
		throw std::invalid_argument("Trajectory_Simulator::Simulate(): initial state, DM mass, or wall-time limit is invalid.");
	Event current_event = initial_condition;
	long unsigned int number_of_scatterings = 0;

	// Initialize per-trajectory bincount
	current_bincount = TrajectoryBincount();
	Reset_Bincount_Anchor(current_event);
	previous_capture_energy_eV = Capture_Energy_eV(current_event.Radius(), current_event.Speed(), DM);
	free_flight_reference_energy_eV = previous_capture_energy_eV;
	current_ballistic_energy_drift_eV = 0.0;
	// Injection is physically unbound by construction. A substantially negative
	// initial boundary energy is handled by the uncaptured-bound mismatch guard;
	// a sub-tolerance sign flip must not reclassify the injected particle.
	current_physical_bound_state = false;
	current_mpi_rank = mpi_rank;
	current_trajectory_id++;
	diagnostic_event_index = 0;
	diagnostic_scatter_index = 0;
	diagnostic_step_index = 0;
	last_scatter_target_index = -2;
	current_diagnostic_events.clear();
	trajectory_in_progress = true;
	steps_since_snapshot_publish = 0;
	track_trajectory_wall_time = snapshot_recorder != nullptr || max_trajectory_wall_time_sec > 0.0;
	if(track_trajectory_wall_time)
		current_trajectory_wall_start = std::chrono::steady_clock::now();
	last_snapshot_publish_wall_time = std::chrono::steady_clock::now();
	accumulated_snapshot_overhead_sec = 0.0;
	if(snapshot_recorder != nullptr)
	{
		const auto snapshot_operation_start = std::chrono::steady_clock::now();
		snapshot_recorder->BeginTrajectory(current_trajectory_id, In_Units(current_event.time, sec));
		Accumulate_Snapshot_Overhead(snapshot_operation_start);
	}

	TrajectoryTerminationReason termination_reason = TrajectoryTerminationReason::Unknown;
	while(number_of_scatterings < maximum_scatterings)
	{
		TrajectoryTerminationReason propagation_reason = Propagate_Freely(current_event, DM);
		if(propagation_reason == TrajectoryTerminationReason::Scatter && current_event.Radius() < rSun)
		{
			diagnostic_scatter_index = number_of_scatterings + 1;
			const size_t scatter_pre_event_index = current_diagnostic_events.size();
			Record_Diagnostic_Event(TrajectoryDiagnosticEventType::ScatterPre, current_event, DM);
			const bool was_captured = current_bincount.is_captured;
			const unsigned long int bound_to_unbound_before = current_bincount.number_of_bound_to_unbound;
			const unsigned long int recaptures_before = current_bincount.number_of_recaptures;
			try
			{
				Scatter(current_event, DM);
			}
			catch(const std::exception& error)
			{
				std::cerr << "\nWarning in Trajectory_Simulator::Simulate(): " << error.what()
				          << " Marking trajectory as numerically failed." << std::endl;
				termination_reason = TrajectoryTerminationReason::NumericalFailure;
				break;
			}
			number_of_scatterings++;
			std::string target_species = "unknown";
			if(last_scatter_target_index == -1)
				target_species = "electron";
			else if(last_scatter_target_index >= 0
			        && static_cast<size_t>(last_scatter_target_index) < solar_model.target_isotopes.size())
				target_species = solar_model.target_isotopes[static_cast<size_t>(last_scatter_target_index)].name;
			if(diagnostic_trace_enabled && scatter_pre_event_index < current_diagnostic_events.size())
			{
				std::strncpy(current_diagnostic_events[scatter_pre_event_index].target_species,
				             target_species.c_str(), sizeof(current_diagnostic_events[scatter_pre_event_index].target_species) - 1);
			}
			if(snapshot_recorder != nullptr)
			{
				const auto snapshot_operation_start = std::chrono::steady_clock::now();
				snapshot_recorder->UpdateCurrentScatterings(number_of_scatterings);
				Accumulate_Snapshot_Overhead(snapshot_operation_start);
			}
			bool captured_after_scatter = Update_Capture_State(current_event.Radius(), current_event.Speed(), current_event.time, DM, true);
			Record_Diagnostic_Event(TrajectoryDiagnosticEventType::ScatterPost, current_event, DM, target_species.c_str());
			if(!was_captured && current_bincount.is_captured)
				Record_Diagnostic_Event(TrajectoryDiagnosticEventType::Capture, current_event, DM, target_species.c_str());
			if(current_bincount.number_of_bound_to_unbound > bound_to_unbound_before)
				Record_Diagnostic_Event(TrajectoryDiagnosticEventType::CandidateUnbinding, current_event, DM, target_species.c_str());
			if(current_bincount.number_of_recaptures > recaptures_before)
				Record_Diagnostic_Event(TrajectoryDiagnosticEventType::Recapture, current_event, DM, target_species.c_str());
			Reset_Bincount_Anchor(current_event);
			if(terminate_on_capture && captured_after_scatter)
			{
				termination_reason = TrajectoryTerminationReason::CaptureMode;
				break;
			}
		}
		else
		{
			termination_reason = propagation_reason;
			break;
		}
	}

	if(termination_reason == TrajectoryTerminationReason::Unknown && number_of_scatterings >= maximum_scatterings)
		termination_reason = TrajectoryTerminationReason::MaxScatterings;

	trajectory_in_progress = false;
	track_trajectory_wall_time = false;
	current_bincount.termination_reason = termination_reason;
	current_bincount.number_of_scatterings = number_of_scatterings;
	current_bincount.t_termination = In_Units(current_event.time, sec);

	if(current_bincount.is_captured)
	{
		double r_final = current_event.Radius();
		double v_final = current_event.Speed();
		double vesc_final = solar_model.Local_Escape_Speed(r_final);
		double E_final = 0.5 * DM.mass * (v_final * v_final - vesc_final * vesc_final);
		double E_final_eV = In_Units(E_final, eV);
		if(termination_reason == TrajectoryTerminationReason::OutwardEscape && E_final_eV >= 0.0)
		{
			if(std::isfinite(current_bincount.t_final_unbinding_scatter))
			{
				current_bincount.boundary_escape_observed = true;
				current_bincount.t_boundary_escape = current_bincount.t_termination;
				current_bincount.r_boundary_escape_km = In_Units(r_final, km);
				current_bincount.vr_boundary_escape_km_s = In_Units(current_event.position.Dot(current_event.velocity) / r_final, km / sec);
				current_bincount.E_boundary_escape_eV = E_final_eV;
				current_bincount.event_observed = true;
				Record_Diagnostic_Event(TrajectoryDiagnosticEventType::EscapeValidated, current_event, DM);
			}
			else
			{
				termination_reason = TrajectoryTerminationReason::EnergyDriftEscape;
				current_bincount.termination_reason = termination_reason;
				current_bincount.numerical_failure_detail =
				    TrajectoryNumericalFailureDetail::BoundaryEnergyMismatch;
				current_bincount.failure_energy_at_boundary_eV =
				    E_final_eV;
				current_bincount.failure_reference_energy_eV =
				    free_flight_reference_energy_eV;
				current_bincount.failure_boundary_vr_km_s =
				    In_Units(
				        current_event.position.Dot(current_event.velocity)
				        / r_final,
				        km / sec);
				current_bincount.survival_valid = false;
				current_bincount.numerically_invalid_escape = true;
				current_bincount.boundary_escape_observed = false;
				current_bincount.event_observed = false;
			}
		}
	}

	if(TrajectoryTerminationInvalidatesSurvival(termination_reason))
		current_bincount.survival_valid = false;

	if(current_bincount.is_captured && !current_bincount.event_observed)
	{
		const bool numerical_failure = termination_reason == TrajectoryTerminationReason::NumericalFailure
		                            || termination_reason == TrajectoryTerminationReason::NonFiniteState
		                            || termination_reason == TrajectoryTerminationReason::SpeedLimit
		                            || termination_reason == TrajectoryTerminationReason::EnergyDriftEscape
		                            || termination_reason == TrajectoryTerminationReason::Unknown;
		Record_Diagnostic_Event(numerical_failure ? TrajectoryDiagnosticEventType::NumericalFailure
		                                          : TrajectoryDiagnosticEventType::Censored,
		                        current_event, DM);
	}

	if(current_bincount.is_captured)
	{
		if(!current_bincount.event_observed)
		{
			current_bincount.t_final_unbinding_scatter = std::numeric_limits<double>::quiet_NaN();
			if(termination_reason != TrajectoryTerminationReason::CaptureMode)
				current_bincount.survival_valid = false;
		}
		current_bincount.outer_domain_removed =
		    termination_reason == TrajectoryTerminationReason::OuterDomainRemoval;
	}

	return Trajectory_Result(initial_condition, current_event, number_of_scatterings, current_bincount,
	                         std::move(current_diagnostic_events));
}  

// 3. Equation of motion solution with Runge-Kutta-Fehlberg
Free_Particle_Propagator::Free_Particle_Propagator(const Event& event)
{
	time = event.time;
	radius = event.Radius();
	phi = 0.0;
	radial_mode = false;

	const double speed = event.Speed();
	const double radius_eps = 1.0e-12 * km;
	const double angular_scale = std::max(radius * speed, 1.0 * km * km / sec);
	const double angular_eps = 1.0e-12 * angular_scale;
	const libphysica::Vector angular_momentum_vector = event.position.Cross(event.velocity);
	const double angular_momentum_norm = angular_momentum_vector.Norm();

	if(std::isfinite(radius) && radius > radius_eps)
	{
		axis_x = event.position / radius;
		v_radial = event.position.Dot(event.velocity) / radius;
	}
	else if(std::isfinite(speed) && speed > 0.0)
	{
		axis_x = event.velocity / speed;
		v_radial = speed;
		radius = 0.0;
		radial_mode = true;
	}
	else
	{
		axis_x = libphysica::Vector({1.0, 0.0, 0.0});
		v_radial = 0.0;
		radius = 0.0;
		radial_mode = true;
	}

	if(!radial_mode && std::isfinite(angular_momentum_norm) && angular_momentum_norm > angular_eps)
	{
		axis_z = angular_momentum_vector / angular_momentum_norm;
		axis_y = axis_z.Cross(axis_x);
		const double axis_y_norm = axis_y.Norm();
		if(std::isfinite(axis_y_norm) && axis_y_norm > 0.0)
		{
			axis_y = axis_y / axis_y_norm;
			axis_x = axis_y.Cross(axis_z).Normalized();
			angular_momentum = angular_momentum_norm;
		}
		else
			radial_mode = true;
	}
	else
		radial_mode = true;

	if(radial_mode)
	{
		Build_Stable_Orthonormal_Basis(axis_x, axis_x, axis_y, axis_z);
		angular_momentum = 0.0;
	}

	// 3. Error tolerances (fixed-size array, no heap allocation)
	error_tolerances[0] = 1.0 * km;
	error_tolerances[1] = 1.0e-3 * km / sec;
	error_tolerances[2] = 1.0e-7;
}

double Free_Particle_Propagator::dr_dt(double v)
{
	return v;
}

double Free_Particle_Propagator::dv_dt(double r, double mass)
{
	// 防止 RK45 中间阶段 r_i ≤ 0（向内飞越球心时的数值越界）
	// 导致离心项 L²/r³ 符号翻转从而使误差发散、步长无限缩小。
	// 1 km 远小于任何真实轨迹的转折半径（~10⁴ km），不影响物理结果。
	if(r < 1.0 * km) r = 1.0 * km;
	return angular_momentum * angular_momentum / r / r / r - G_Newton * mass / r / r;
}

double Free_Particle_Propagator::dphi_dt(double r)
{
	if(radial_mode || angular_momentum == 0.0 || r <= 0.0)
		return 0.0;
	return angular_momentum / r / r;
}

bool Free_Particle_Propagator::Runge_Kutta_45_Step(Solar_Model& solar_model)
{
	time_step = RK45_Sanitized_Time_Step(time_step);
	bool accepted = false;
	int inner_iter = 0;
	while(!accepted)
	{
	inner_iter++;
	// RK coefficients:
	double k_r[6];
	double k_v[6];
	double k_p[6];
	double r_i;  // intermediate radius for each RK stage

	// Stage 0
	k_r[0] = time_step * dr_dt(v_radial);
	k_v[0] = time_step * dv_dt(radius, solar_model.Mass(radius));
	k_p[0] = time_step * dphi_dt(radius);

	// Stage 1
	r_i = radius + k_r[0] / 4.0;
	k_r[1] = time_step * dr_dt(v_radial + k_v[0] / 4.0);
	k_v[1] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	// k_p[1]=	dt*dphi_dt(radius+k_r[0]/4.0,J);

	// Stage 2
	r_i = radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1];
	k_r[2] = time_step * dr_dt(v_radial + 3.0 / 32.0 * k_v[0] + 9.0 / 32.0 * k_v[1]);
	k_v[2] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	k_p[2] = time_step * dphi_dt(r_i);

	// Stage 3
	r_i = radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2];
	k_r[3] = time_step * dr_dt(v_radial + 1932.0 / 2197.0 * k_v[0] - 7200.0 / 2197.0 * k_v[1] + 7296.0 / 2197.0 * k_v[2]);
	k_v[3] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	k_p[3] = time_step * dphi_dt(r_i);

	// Stage 4
	r_i = radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3];
	k_r[4] = time_step * dr_dt(v_radial + 439.0 / 216.0 * k_v[0] - 8.0 * k_v[1] + 3680.0 / 513.0 * k_v[2] - 845.0 / 4104.0 * k_v[3]);
	k_v[4] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	k_p[4] = time_step * dphi_dt(r_i);

	// Stage 5
	r_i = radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4];
	k_r[5] = time_step * dr_dt(v_radial - 8.0 / 27.0 * k_v[0] + 2.0 * k_v[1] - 3544.0 / 2565.0 * k_v[2] + 1859.0 / 4104.0 * k_v[3] - 11.0 / 40.0 * k_v[4]);
	k_v[5] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	k_p[5] = time_step * dphi_dt(r_i);

	// New values with Runge Kutta 4 and Runge Kutta 5
	double radius_4	  = radius + 25.0 / 216.0 * k_r[0] + 1408.0 / 2565.0 * k_r[2] + 2197.0 / 4104.0 * k_r[3] - 1.0 / 5.0 * k_r[4];
	double v_radial_4 = v_radial + 25.0 / 216.0 * k_v[0] + 1408.0 / 2565.0 * k_v[2] + 2197.0 / 4104.0 * k_v[3] - 1.0 / 5.0 * k_v[4];
	double phi_4	  = phi + 25.0 / 216.0 * k_p[0] + 1408.0 / 2565.0 * k_p[2] + 2197.0 / 4104.0 * k_p[3] - 1.0 / 5.0 * k_p[4];

	double radius_5	  = radius + 16.0 / 135.0 * k_r[0] + 6656.0 / 12825.0 * k_r[2] + 28561.0 / 56430.0 * k_r[3] - 9.0 / 50.0 * k_r[4] + 2.0 / 55.0 * k_r[5];
	double v_radial_5 = v_radial + 16.0 / 135.0 * k_v[0] + 6656.0 / 12825.0 * k_v[2] + 28561.0 / 56430.0 * k_v[3] - 9.0 / 50.0 * k_v[4] + 2.0 / 55.0 * k_v[5];
	double phi_5	  = phi + 16.0 / 135.0 * k_p[0] + 6656.0 / 12825.0 * k_p[2] + 28561.0 / 56430.0 * k_p[3] - 9.0 / 50.0 * k_p[4] + 2.0 / 55.0 * k_p[5];

	// Error and adapting the time step (B: 栈数组替代堆分配vector)
	double errors[3] = {fabs(radius_5 - radius_4), fabs(v_radial_5 - v_radial_4), fabs(phi_5 - phi_4)};
	double time_step_new = RK45_Next_Step_Size(time_step, errors, error_tolerances);

	// Check if errors fall below the tolerance
	// 自适应步长方法：根据误差调整 time_step
	if(RK45_Errors_Within_Tolerance(errors, error_tolerances))
	{
		time	  = time + time_step;
		radius	  = radius_4;
		// 边界检查：防止半径变为负数（数值误差）
		if(radius < 0.0)
			radius = 0.0;
		v_radial  = v_radial_4;
		phi		  = phi_4;
		time_step = time_step_new;
		accepted  = true;
		return true;
	}
	else
	{
		// 绝对最小步长保护：防止步长无限缩小导致死循环
		// 当粒子轨迹中间阶段 r_i 为负值时（角动量 L≠0 且接近 r≈0），
		// dv_dt = L²/r³ 会发散，误差无法满足容差，步长每次乘 0.1 直到下溢。
		// 强制接受：确保内层循环始终能返回，snapshot callback 才能被触发。
		static const double abs_min_step = 1.0e-8 * sec;  // 10 纳秒物理时间下限
		const bool reached_min_step  = time_step_new < abs_min_step;
		const bool reached_iter_cap  = inner_iter >= RK45_MAX_INNER_RETRIES;
		if(reached_min_step || reached_iter_cap)
		{
			// NaN 守卫：若本轮 radius_4/v_radial_4/phi_4 已经是非有限数，
			// 绝不能把 NaN 写回传播器状态，否则后续的每一步都会继续崩。
			// 这时保留原状态，只把时间推进 abs_min_step，由外层 Propagate_Freely
			// 的 is_finite 检查来中止本条轨迹。
		const bool state_ok =
		    std::isfinite(radius_4) && std::isfinite(v_radial_4) && std::isfinite(phi_4);
		time = time + abs_min_step;
		if(state_ok)
		{
			radius = radius_4;
			if(radius < 0.0)
				radius = 0.0;
			v_radial = v_radial_4;
			phi      = phi_4;
		}
		time_step = abs_min_step;
		accepted  = true;
		return false;
	}
	else
	{
		time_step = time_step_new;
		// Loop continues with smaller time_step (was recursive call before Q1 fix)
	}
	}
	}   // end while(!accepted)
	return true;
}

// Overload with constant mass (for Kepler orbits outside the sun / testing)
bool Free_Particle_Propagator::Runge_Kutta_45_Step(double constant_mass)
{
	time_step = RK45_Sanitized_Time_Step(time_step);
	bool accepted = false;
	int inner_iter = 0;
	while(!accepted)
	{
	inner_iter++;
	double k_r[6];
	double k_v[6];
	double k_p[6];

	k_r[0] = time_step * dr_dt(v_radial);
	k_v[0] = time_step * dv_dt(radius, constant_mass);
	k_p[0] = time_step * dphi_dt(radius);

	k_r[1] = time_step * dr_dt(v_radial + k_v[0] / 4.0);
	k_v[1] = time_step * dv_dt(radius + k_r[0] / 4.0, constant_mass);

	k_r[2] = time_step * dr_dt(v_radial + 3.0 / 32.0 * k_v[0] + 9.0 / 32.0 * k_v[1]);
	k_v[2] = time_step * dv_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1], constant_mass);
	k_p[2] = time_step * dphi_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1]);

	k_r[3] = time_step * dr_dt(v_radial + 1932.0 / 2197.0 * k_v[0] - 7200.0 / 2197.0 * k_v[1] + 7296.0 / 2197.0 * k_v[2]);
	k_v[3] = time_step * dv_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2], constant_mass);
	k_p[3] = time_step * dphi_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2]);

	k_r[4] = time_step * dr_dt(v_radial + 439.0 / 216.0 * k_v[0] - 8.0 * k_v[1] + 3680.0 / 513.0 * k_v[2] - 845.0 / 4104.0 * k_v[3]);
	k_v[4] = time_step * dv_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3], constant_mass);
	k_p[4] = time_step * dphi_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3]);

	k_r[5] = time_step * dr_dt(v_radial - 8.0 / 27.0 * k_v[0] + 2.0 * k_v[1] - 3544.0 / 2565.0 * k_v[2] + 1859.0 / 4104.0 * k_v[3] - 11.0 / 40.0 * k_v[4]);
	k_v[5] = time_step * dv_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4], constant_mass);
	k_p[5] = time_step * dphi_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4]);

	double radius_4	  = radius + 25.0 / 216.0 * k_r[0] + 1408.0 / 2565.0 * k_r[2] + 2197.0 / 4104.0 * k_r[3] - 1.0 / 5.0 * k_r[4];
	double v_radial_4 = v_radial + 25.0 / 216.0 * k_v[0] + 1408.0 / 2565.0 * k_v[2] + 2197.0 / 4104.0 * k_v[3] - 1.0 / 5.0 * k_v[4];
	double phi_4	  = phi + 25.0 / 216.0 * k_p[0] + 1408.0 / 2565.0 * k_p[2] + 2197.0 / 4104.0 * k_p[3] - 1.0 / 5.0 * k_p[4];

	double radius_5	  = radius + 16.0 / 135.0 * k_r[0] + 6656.0 / 12825.0 * k_r[2] + 28561.0 / 56430.0 * k_r[3] - 9.0 / 50.0 * k_r[4] + 2.0 / 55.0 * k_r[5];
	double v_radial_5 = v_radial + 16.0 / 135.0 * k_v[0] + 6656.0 / 12825.0 * k_v[2] + 28561.0 / 56430.0 * k_v[3] - 9.0 / 50.0 * k_v[4] + 2.0 / 55.0 * k_v[5];
	double phi_5	  = phi + 16.0 / 135.0 * k_p[0] + 6656.0 / 12825.0 * k_p[2] + 28561.0 / 56430.0 * k_p[3] - 9.0 / 50.0 * k_p[4] + 2.0 / 55.0 * k_p[5];

	double errors[3] = {fabs(radius_5 - radius_4), fabs(v_radial_5 - v_radial_4), fabs(phi_5 - phi_4)};
	double time_step_new = RK45_Next_Step_Size(time_step, errors, error_tolerances);

	if(RK45_Errors_Within_Tolerance(errors, error_tolerances))
	{
		time	  = time + time_step;
		radius	  = radius_4;
		if(radius < 0.0)
			radius = 0.0;
		v_radial  = v_radial_4;
		phi		  = phi_4;
		time_step = time_step_new;
		accepted  = true;
		return true;
	}
	else
	{
		static const double abs_min_step = 1.0e-8 * sec;  // 10 纳秒物理时间下限（同 Solar Model 版本）
		const bool reached_min_step  = time_step_new < abs_min_step;
		const bool reached_iter_cap  = inner_iter >= RK45_MAX_INNER_RETRIES;
		if(reached_min_step || reached_iter_cap)
		{
		const bool state_ok =
		    std::isfinite(radius_4) && std::isfinite(v_radial_4) && std::isfinite(phi_4);
		time = time + abs_min_step;
		if(state_ok)
		{
			radius = radius_4;
			if(radius < 0.0)
				radius = 0.0;
			v_radial = v_radial_4;
			phi      = phi_4;
		}
		time_step = abs_min_step;
		accepted  = true;
		return false;
	}
	else
	{
		time_step = time_step_new;
	}
	}
	}
	return true;
}

double Free_Particle_Propagator::Current_Time()
{
	return time;
}

double Free_Particle_Propagator::Current_Radius()
{
	// 边界检查：确保返回非负半径
	return (radius < 0.0) ? 0.0 : radius;
}

double Free_Particle_Propagator::Current_Speed()
{
	// 边界检查：确保 radius 非负
	double r = (radius < 0.0) ? 0.0 : radius;
	if(r == 0 || angular_momentum == 0)
		return fabs(v_radial);  // 返回速度绝对值
	else
		return sqrt(v_radial * v_radial + angular_momentum * angular_momentum / r / r);
}

Free_Particle_Propagator::Scalar_State Free_Particle_Propagator::Save_Scalar_State() const
{
	Scalar_State state;
	state.time		= time;
	state.radius	= radius;
	state.phi		= phi;
	state.v_radial	= v_radial;
	state.time_step = time_step;
	return state;
}

void Free_Particle_Propagator::Restore_Scalar_State(const Scalar_State& state)
{
	time	  = state.time;
	radius	  = state.radius;
	phi		  = state.phi;
	v_radial  = state.v_radial;
	time_step = state.time_step;
}

void Free_Particle_Propagator::Fill_Event_In_3D(Event& event) const
{
	if(event.position.Size() != 3)
		event.position = libphysica::Vector(3);
	if(event.velocity.Size() != 3)
		event.velocity = libphysica::Vector(3);

	const double safe_radius = (std::isfinite(radius) && radius > 0.0) ? radius : 0.0;
	event.time				 = time;
	if(radial_mode || angular_momentum == 0.0 || safe_radius == 0.0)
	{
		for(unsigned int component = 0; component < 3; component++)
		{
			event.position[component] = safe_radius * axis_x[component];
			event.velocity[component] = v_radial * axis_x[component];
		}
		return;
	}

	const double cos_phi = cos(phi);
	const double sin_phi = sin(phi);
	const double tangential_speed = angular_momentum / safe_radius;
	const double velocity_x = v_radial * cos_phi - tangential_speed * sin_phi;
	const double velocity_y = v_radial * sin_phi + tangential_speed * cos_phi;
	for(unsigned int component = 0; component < 3; component++)
	{
		event.position[component] = safe_radius * (cos_phi * axis_x[component] + sin_phi * axis_y[component]);
		event.velocity[component] = velocity_x * axis_x[component] + velocity_y * axis_y[component];
	}
}

Event Free_Particle_Propagator::Event_In_3D()
{
	Event event;
	Fill_Event_In_3D(event);
	return event;
}

}	// namespace DaMaSCUS_SUN

#include "Simulation_Utilities.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/DM_Halo_Models.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

namespace
{
constexpr double KEPLER_DOMAIN_TOLERANCE = 1.0e-10;

double Clamp_Cosine(double value)
{
	return std::max(-1.0, std::min(1.0, value));
}

bool Report_Hyperbolic_Kepler_Shift_Failure(const std::string& message)
{
	static unsigned int warning_count = 0;
	if(warning_count < 10)
		std::cerr << "Warning in Hyperbolic_Kepler_Shift(): " << message << std::endl;
	else if(warning_count == 10)
		std::cerr << "Warning in Hyperbolic_Kepler_Shift(): additional failures suppressed." << std::endl;
	warning_count++;
	return false;
}

bool Finite_Positive(double value)
{
	return std::isfinite(value) && value > 0.0;
}

}	// namespace

// 1. Event Class
Event::Event()
: time(0.0), position({0, 0, 0}), velocity({0, 0, 0})
{
}

Event::Event(double t, const libphysica::Vector& pos, const libphysica::Vector& vel)
: time(t), position(pos), velocity(vel)
{
}

double Event::Radius() const
{
	return position.Norm();
}

double Event::Speed() const
{
	return velocity.Norm();
}

double Event::Angular_Momentum() const
{
	return position.Cross(velocity).Norm();
}

double Event::Asymptotic_Speed_Sqr(Celestial_Model& solar_model) const
{
	double r	 = Radius();
	double v	 = Speed();
	double vesc	 = solar_model.Local_Escape_Speed(r);
	double u_sqr = v * v - vesc * vesc;
	return u_sqr;
}

double Event::Isoreflection_Angle(const libphysica::Vector& vel_sun) const
{
	const double position_norm = position.Norm();
	const double sun_speed = vel_sun.Norm();
	if(!Finite_Positive(position_norm) || !Finite_Positive(sun_speed))
		throw std::runtime_error("Event::Isoreflection_Angle(): position and solar velocity must be finite and non-zero.");
	const double cosine = position.Dot(vel_sun) / position_norm / sun_speed;
	if(!std::isfinite(cosine))
		throw std::runtime_error("Event::Isoreflection_Angle(): angle cosine is non-finite.");
	return acos(Clamp_Cosine(cosine));
}

int Event::Isoreflection_Ring(const libphysica::Vector& vel_sun, unsigned int number_of_rings) const
{
	double theta					= Isoreflection_Angle(vel_sun);
	std::vector<double> ring_angles = Isoreflection_Ring_Angles(number_of_rings);
	for(unsigned int ring = 0; ring < ring_angles.size(); ring++)
		if(theta <= ring_angles[ring])
			return ring;
	std::cerr << "Error in Event::Isoreflection_Ring(): Angle = " << theta << " out of bound." << std::endl;
	std::exit(EXIT_FAILURE);
}

Event Event::In_Units(double unit_distance, double unit_time) const
{
	return Event(libphysica::natural_units::In_Units(time, unit_time), libphysica::natural_units::In_Units(position, unit_distance), libphysica::natural_units::In_Units(velocity, unit_distance / unit_time));
}

// Overload <<
std::ostream& operator<<(std::ostream& output, const Event& event)
{
	return output << "{"
				  << event.time
				  << ","
				  << event.position
				  << ","
				  << event.velocity
				  << "}";
}

// 2. Generator of initial conditions
double PDF_Initial_Speed(double v, obscura::DM_Distribution& halo_model, Celestial_Model& solar_model)
{
	if(!std::isfinite(v) || v < 0.0)
		throw std::invalid_argument("PDF_Initial_Speed(): speed must be finite and non-negative.");
	// The physical speed PDF is O(v^2), so the flux-weighted distribution has
	// a finite zero limit despite the apparent v_esc^2/v factor.
	if(v == 0.0)
		return 0.0;

	double v_esc			 = solar_model.Local_Escape_Speed(g_body_radius);
	double v_average		 = halo_model.Average_Speed();
	double v_inverse_average = halo_model.Eta_Function(0.0);
	const double normalization = v_average + v_esc * v_esc * v_inverse_average;
	const double speed_pdf = halo_model.PDF_Speed(v);
	if(!Finite_Positive(normalization) || !std::isfinite(speed_pdf) || speed_pdf < 0.0)
		throw std::runtime_error("PDF_Initial_Speed(): halo normalization or speed PDF is invalid.");
	const double pdf = (speed_pdf * v + (speed_pdf / v) * v_esc * v_esc) / normalization;
	if(!std::isfinite(pdf) || pdf < 0.0)
		throw std::runtime_error("PDF_Initial_Speed(): flux-weighted speed PDF is invalid.");
	return pdf;
}

// Conditional pdf for cos_theta given a speed value v
double PDF_Cos_Theta(double cos_theta, double v, obscura::DM_Distribution& halo_model)
{
	if(!std::isfinite(cos_theta) || cos_theta < -1.0 || cos_theta > 1.0)
		throw std::invalid_argument("PDF_Cos_Theta(): cosine must be finite and lie in [-1, 1].");
	if(!Finite_Positive(v))
		throw std::invalid_argument("PDF_Cos_Theta(): speed must be finite and positive.");
	const double speed_pdf = halo_model.PDF_Speed(v);
	if(!Finite_Positive(speed_pdf))
		throw std::runtime_error("PDF_Cos_Theta(): halo speed PDF must be finite and positive.");
	double normalization = 2.0 * M_PI * v * v / speed_pdf;

	// Construct velocity vector with angle theta to the Sun's velocity
	auto* standard_halo = dynamic_cast<obscura::Standard_Halo_Model*>(&halo_model);
	if(standard_halo == nullptr)
		throw std::invalid_argument("PDF_Cos_Theta(): initial-condition sampling requires a Standard_Halo_Model.");
	libphysica::Vector vel_sun = standard_halo->Get_Observer_Velocity();
	const double observer_speed = vel_sun.Norm();
	if(!std::isfinite(observer_speed))
		throw std::runtime_error("PDF_Cos_Theta(): observer velocity is non-finite.");
	const libphysica::Vector angular_axis = (observer_speed > 0.0)
	                                           ? vel_sun
	                                           : libphysica::Vector({0.0, 0.0, 1.0});
	libphysica::Vector vel	   = libphysica::Spherical_Coordinates(v, acos(cos_theta), 0.0, angular_axis);

	const double pdf = normalization * halo_model.PDF_Velocity(vel);
	if(!std::isfinite(pdf) || pdf < 0.0)
		throw std::runtime_error("PDF_Cos_Theta(): conditional angular PDF is invalid.");
	return pdf;
}

Event Initial_Conditions(obscura::DM_Distribution& halo_model, Celestial_Model& solar_model, std::mt19937& PRNG)
{
	auto* standard_halo = dynamic_cast<obscura::Standard_Halo_Model*>(&halo_model);
	if(standard_halo == nullptr)
		throw std::invalid_argument("Initial_Conditions(): only Standard_Halo_Model distributions are supported.");
	// 1. Initial velocity
	// 1.1. Sample initial speed u asymptotically far from the Sun.
	std::function<double(double)> pdf_v = [&halo_model, &solar_model](double v) {
		return PDF_Initial_Speed(v, halo_model, solar_model);
	};
	double u = libphysica::Rejection_Sampling(pdf_v, halo_model.Minimum_DM_Speed(), halo_model.Maximum_DM_Speed(), 1200.0, PRNG);
	if(!Finite_Positive(u))
		throw std::runtime_error("Initial_Conditions(): sampled asymptotic speed is invalid.");

	// 1.2. Sample cos(theta) where theta is the angle between v and v_sun.
	std::function<double(double)> pdf_cos_theta = [u, &halo_model](double cos_theta) {
		return PDF_Cos_Theta(cos_theta, u, halo_model);
	};
	// Determine domain of cos theta
	libphysica::Vector vel_sun = standard_halo->Get_Observer_Velocity();
	double v_sun			   = vel_sun.Norm();
	if(!std::isfinite(v_sun))
		throw std::runtime_error("Initial_Conditions(): observer velocity is non-finite.");
	double v_gal			   = halo_model.Maximum_DM_Speed() - v_sun;
	double cos_theta_max	   = (v_sun == 0.0)
	                              ? 1.0
	                              : std::min(1.0, (v_gal * v_gal - v_sun * v_sun - u * u) / (2.0 * u * v_sun));
	if(!std::isfinite(cos_theta_max) || cos_theta_max < -1.0)
		throw std::runtime_error("Initial_Conditions(): angular sampling domain is empty or non-finite.");

	double y_max	 = PDF_Cos_Theta(-1.0, u, halo_model);
	double cos_theta = libphysica::Rejection_Sampling(pdf_cos_theta, -1.0, cos_theta_max, y_max, PRNG);
	// double cos_theta = libphysica::Sample_Uniform(PRNG,-1.0,1.0);// to test isotropic initial conditions

	// 1.3. Construct velocity vector
	double phi							= libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	const libphysica::Vector angular_axis = (v_sun > 0.0)
	                                           ? vel_sun
	                                           : libphysica::Vector({0.0, 0.0, 1.0});
	libphysica::Vector initial_velocity = libphysica::Spherical_Coordinates(u, acos(Clamp_Cosine(cos_theta)), phi, angular_axis);
	if(!Finite_Positive(initial_velocity.Norm()))
		throw std::runtime_error("Initial_Conditions(): sampled velocity is zero or non-finite.");

	// 1.4. Blue-shift the speed
	double asymptotic_distance = 1000.0 * AU;
	double vesc_asymptotic	   = solar_model.Local_Escape_Speed(asymptotic_distance);
	double v				   = sqrt(u * u + vesc_asymptotic * vesc_asymptotic);
	if(!std::isfinite(vesc_asymptotic) || vesc_asymptotic < 0.0 || !Finite_Positive(v))
		throw std::runtime_error("Initial_Conditions(): asymptotic escape speed is invalid.");
	initial_velocity		   = v * initial_velocity.Normalized();

	// 2. Initial position
	// 2.1 Find the maximum impact parameter such that the particle still hits the Sun.
	double v_esc				= solar_model.Local_Escape_Speed(g_body_radius);
	double impact_parameter_max = sqrt(u * u + v_esc * v_esc) / v * g_body_radius;
	if(!std::isfinite(v_esc) || v_esc < 0.0 || !Finite_Positive(impact_parameter_max))
		throw std::runtime_error("Initial_Conditions(): solar escape speed or impact-parameter bound is invalid.");
	libphysica::Vector e_z		= (-1.0) * initial_velocity.Normalized();
	const libphysica::Vector reference = (std::fabs(e_z[2]) < 0.9)
	                                         ? libphysica::Vector({0.0, 0.0, 1.0})
	                                         : libphysica::Vector({1.0, 0.0, 0.0});
	libphysica::Vector e_x = reference.Cross(e_z).Normalized();
	libphysica::Vector e_y = e_z.Cross(e_x).Normalized();

	// 2.2 Find a random point in the plane.
	double phi_disk						= libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	double xi							= libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
	double impact_parameter				= sqrt(xi) * impact_parameter_max;
	libphysica::Vector initial_position = asymptotic_distance * e_z + impact_parameter * (cos(phi_disk) * e_x + sin(phi_disk) * e_y);
	if(!Finite_Positive(initial_position.Norm()))
		throw std::runtime_error("Initial_Conditions(): sampled position is zero or non-finite.");

	return Event(0.0, initial_position, initial_velocity);
}

// 3. Analytically propagate a particle at event on a hyperbolic Kepler orbit to a radius R (without passing the periapsis)
bool Hyperbolic_Kepler_Shift(Event& event, double R_final)
{
	const double mu = G_Newton * g_body_mass;
	const double R_initial = event.Radius();
	const double speed = event.Speed();

	if(R_final < g_body_radius || R_initial < g_body_radius)
		return Report_Hyperbolic_Kepler_Shift_Failure("orbits inside the Sun cannot be described analytically.");
	if(!Finite_Positive(R_final) || !Finite_Positive(R_initial) || !Finite_Positive(speed))
		return Report_Hyperbolic_Kepler_Shift_Failure("initial or final radius/speed is non-finite or non-positive.");
	if(std::fabs(R_final - R_initial) <= KEPLER_DOMAIN_TOLERANCE * std::max(R_initial, R_final))
		return true;

	const double radial_velocity = event.position.Dot(event.velocity) / R_initial;
	if(!std::isfinite(radial_velocity) || radial_velocity == 0.0)
		return Report_Hyperbolic_Kepler_Shift_Failure("cannot choose a forward Kepler branch at zero or non-finite radial velocity.");

	const double radial_direction = (radial_velocity < 0.0) ? -1.0 : 1.0;
	if(radial_direction < 0.0 && R_final > R_initial)
		return Report_Hyperbolic_Kepler_Shift_Failure("requested outward shift while particle is on the inbound branch.");
	if(radial_direction > 0.0 && R_final < R_initial)
		return Report_Hyperbolic_Kepler_Shift_Failure("requested inward shift while particle is on the outbound branch.");

	const double energy = 0.5 * speed * speed - mu / R_initial;
	if(!Finite_Positive(energy))
		return Report_Hyperbolic_Kepler_Shift_Failure("orbit is not hyperbolic.");

	const libphysica::Vector r_hat = event.position / R_initial;
	const libphysica::Vector h_vec = event.position.Cross(event.velocity);
	const double h = h_vec.Norm();
	const double h_scale = R_initial * speed;

	if(!std::isfinite(h) || h <= KEPLER_DOMAIN_TOLERANCE * h_scale)
	{
		const double final_speed_sqr = 2.0 * (energy + mu / R_final);
		if(!Finite_Positive(final_speed_sqr))
			return Report_Hyperbolic_Kepler_Shift_Failure("radial branch has non-finite or non-positive final speed.");
		event.position = R_final * r_hat;
		event.velocity = radial_direction * sqrt(final_speed_sqr) * r_hat;
		return true;
	}

	const double p = h * h / mu;
	const libphysica::Vector eccentricity_vector = event.velocity.Cross(h_vec) / mu - r_hat;
	const double eccentricity = eccentricity_vector.Norm();
	if(!Finite_Positive(p) || !std::isfinite(eccentricity) || eccentricity <= 1.0)
		return Report_Hyperbolic_Kepler_Shift_Failure("invalid hyperbolic orbital elements.");

	const double cos_theta_final_raw = (p / R_final - 1.0) / eccentricity;
	if(!std::isfinite(cos_theta_final_raw)
	   || cos_theta_final_raw > 1.0 + KEPLER_DOMAIN_TOLERANCE
	   || cos_theta_final_raw < -1.0 - KEPLER_DOMAIN_TOLERANCE)
		return Report_Hyperbolic_Kepler_Shift_Failure("target radius is not reachable on this Kepler orbit.");

	libphysica::Vector axis_x = eccentricity_vector / eccentricity;
	libphysica::Vector axis_z = h_vec / h;
	libphysica::Vector axis_y = axis_z.Cross(axis_x);
	const double axis_y_norm = axis_y.Norm();
	if(!Finite_Positive(axis_y_norm))
		return Report_Hyperbolic_Kepler_Shift_Failure("failed to construct orbital basis.");
	axis_y = axis_y / axis_y_norm;
	axis_x = axis_y.Cross(axis_z).Normalized();

	const double theta_final_abs = acos(Clamp_Cosine(cos_theta_final_raw));
	const double theta_final = radial_direction * theta_final_abs;

	event.position = R_final * (cos(theta_final) * axis_x + sin(theta_final) * axis_y);
	event.velocity = sqrt(mu / p) * (-sin(theta_final) * axis_x + (eccentricity + cos(theta_final)) * axis_y);
	return true;
}

// 4. Equiareal isodetection rings
std::vector<double> Isoreflection_Ring_Angles(unsigned int number_of_rings)
{
	if(number_of_rings == 0)
		throw std::invalid_argument("Isoreflection_Ring_Angles(): number_of_rings must be positive.");
	std::vector<double> thetas;
	double theta = 0;
	for(unsigned int i = 0; i < number_of_rings - 1; i++)
	{
		theta = acos(cos(theta) - 2.0 / number_of_rings);
		thetas.push_back(theta);
	}
	thetas.push_back(180. * deg);
	return thetas;
}

}	// namespace DaMaSCUS_SUN

#include "Celestial_Model.hpp"

#include "libphysica/Natural_Units.hpp"

namespace DaMaSCUS_SUN
{
// Default to the Sun so pre-existing Sun runs are unchanged.
double g_body_radius = libphysica::natural_units::rSun;
double g_body_mass   = libphysica::natural_units::mSun;
}	// namespace DaMaSCUS_SUN

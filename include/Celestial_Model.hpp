#ifndef __Celestial_Model_hpp_
#define __Celestial_Model_hpp_

#include <string>

#include "obscura/DM_Particle.hpp"
#include "obscura/Target_Nucleus.hpp"

namespace DaMaSCUS_SUN
{

// Abstract base class for a celestial body (Sun, Earth, ...).
// The simulation engine talks to this interface so that the body's
// geometry (Radius) and physics are supplied polymorphically instead
// of being hard-coded to the Sun.
class Celestial_Model
{
  public:
    virtual ~Celestial_Model() = default;

    virtual const std::string& Name() const = 0;

    // Geometry
    virtual double Radius() const     = 0;
    virtual double Total_Mass() const = 0;

    // Bulk profiles
    virtual double Mass(double r)                          = 0;
    virtual double Mass_Density(double r)                  = 0;
    virtual double Temperature(double r)                   = 0;
    virtual double Local_Escape_Speed(double r)            = 0;
    virtual double Debye_Screening_Scale_Squared(double r) = 0;

    // Targets
    virtual unsigned int Target_Count() const                                  = 0;
    virtual const obscura::Isotope& Target_Isotope(unsigned int index) const   = 0;
    virtual double Number_Density_Nucleus(double r, unsigned int nucleus_index) = 0;
    virtual double Number_Density_Electron(double r)                           = 0;

    // DM scattering
    virtual double DM_Scattering_Rate_Electron(obscura::DM_Particle& DM, double r, double DM_speed)                            = 0;
    virtual double DM_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double r, double DM_speed, unsigned int nucleus_index) = 0;
    virtual double Total_DM_Scattering_Rate(obscura::DM_Particle& DM, double r, double DM_speed)                               = 0;
    virtual void Interpolate_Total_DM_Scattering_Rate(obscura::DM_Particle& DM, unsigned int N_radius, unsigned int N_speed)   = 0;
};


// Global current-body parameters (set at startup per target_body).
// Default to Sun so existing Sun runs are byte-for-byte unchanged.
extern double g_body_radius;
extern double g_body_mass;

}	// namespace DaMaSCUS_SUN

#endif

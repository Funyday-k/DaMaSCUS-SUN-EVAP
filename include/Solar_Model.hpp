#ifndef __Solar_Model_hpp_
#define __Solar_Model_hpp_

#include <string>

#include "libphysica/Linear_Algebra.hpp"
#include "libphysica/Numerics.hpp"

#include "obscura/DM_Particle.hpp"
#include "obscura/Target_Nucleus.hpp"

namespace DaMaSCUS_SUN
{

// Locate the AGSS09 input table for build-tree and installed executions.
// DAMASCUS_SUN_SOLAR_MODEL (file) and DAMASCUS_SUN_DATA_DIR (directory)
// provide explicit runtime overrides.
std::string Locate_Solar_Model_Data_File(const std::string& executable_hint = "");

// 1. Nuclear targets in the Sun
class Solar_Isotope : public obscura::Isotope
{
  private:
	libphysica::Interpolation number_density;

  public:
	Solar_Isotope(const obscura::Isotope& isotope, const std::vector<std::vector<double>>& density_table, double abundance = 1.0);

	double Number_Density(double r);
};

// 2. Solar model
class Solar_Model
{
  private:
	libphysica::Interpolation mass, temperature, local_escape_speed_squared, mass_density;

	// Auxiliary functions for the data import
	std::vector<std::vector<double>> raw_data;
	void Import_Raw_Data(const std::string& data_file);
	std::vector<std::vector<double>> Create_Interpolation_Table(unsigned int row) const;
	std::vector<std::vector<double>> Create_Escape_Speed_Table();
	std::vector<std::vector<double>> Create_Number_Density_Table(unsigned int target, double mass) const;
	std::vector<std::vector<double>> Create_Number_Density_Table_Electron();

	// Solar electrons
	libphysica::Interpolation number_density_electron;

	// Interpolation of total scattering rate
	bool using_interpolated_rate;
	std::vector<double> rate_grid;
	unsigned int rate_grid_radius_points;
	unsigned int rate_grid_speed_points;
	double rate_grid_inverse_radius_step;
	double rate_grid_inverse_speed_step;
	double rate_grid_max_speed;

  public:
	std::string name;
	std::vector<Solar_Isotope> target_isotopes;

	explicit Solar_Model(const std::string& data_file = "");

	double Mass(double r);
	double Mass_Density(double r);
	double Temperature(double r);
	double Local_Escape_Speed(double r);
	double Debye_Screening_Scale_Squared(double r);

	double Number_Density_Nucleus(double r, unsigned int nucleus_index);
	double Number_Density_Electron(double r);

	double DM_Scattering_Rate_Electron(obscura::DM_Particle& DM, double r, double DM_speed);
	double DM_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double r, double DM_speed, unsigned int nucleus_index);

	double Total_DM_Scattering_Rate(obscura::DM_Particle& DM, double r, double DM_speed);
	double Total_DM_Scattering_Rate_Computed(obscura::DM_Particle& DM, double r, double DM_speed);

	double Total_DM_Scattering_Rate_Interpolated(obscura::DM_Particle& DM, double r, double DM_speed);
	void Interpolate_Total_DM_Scattering_Rate(obscura::DM_Particle& DM, unsigned int N_radius, unsigned int N_speed);
	unsigned int Scattering_Rate_Interpolation_Radius_Points() const;
	unsigned int Scattering_Rate_Interpolation_Speed_Points() const;

	void Print_Summary(int mpi_rank = 0) const;
};

// 3. Thermal average of relative speed between a particle of speed v_DM and a solar thermal target.
extern double Thermal_Averaged_Relative_Speed(double temperature, double mass_target, double v_DM);
}	// namespace DaMaSCUS_SUN
#endif

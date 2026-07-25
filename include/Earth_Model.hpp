#ifndef __Earth_Model_hpp_
#define __Earth_Model_hpp_

#include <string>
#include <vector>

#include "libphysica/Numerics.hpp"

#include "Celestial_Model.hpp"

namespace DaMaSCUS_SUN
{

class Earth_Isotope : public obscura::Isotope
{
  private:
        std::vector<double> layer_mass_fractions;

  public:
        Earth_Isotope(const obscura::Isotope& isotope, const std::vector<double>& mass_fractions_by_layer);
        double Mass_Fraction(unsigned int layer_index) const;
};


class Earth_Model : public Celestial_Model
{
  private:
	std::string model_file;
	std::vector<std::vector<double>> raw_data;
	libphysica::Interpolation mass, temperature, local_escape_speed_squared, mass_density;
	bool using_interpolated_rate;
	libphysica::Interpolation_2D rate_interpolation;
	std::vector<double> layer_outer_radii;

	std::vector<std::vector<double>> raw_element_ppm;
    std::vector<double> layer_boundaries;
	std::vector<double> layer_radii;       // 新增：儲存每一層的半徑


	void Import_Raw_Data(const std::string& path);
	std::vector<std::vector<double>> Create_Interpolation_Table(unsigned int column) const;
	std::vector<std::vector<double>> Create_Escape_Speed_Table();
	void Initialize_Composition();
	unsigned int Layer_Index(double r) const;

  public:
	std::string name;
	std::vector<Earth_Isotope> target_isotopes;

	Earth_Model();
	explicit Earth_Model(const std::string& path);

	const std::string& Name() const override;
	const std::string& Model_File() const;
	double Radius() const override;
	double Total_Mass() const override;

	double Mass(double r) override;
	double Mass_Density(double r) override;
	double Temperature(double r) override;
	double Local_Escape_Speed(double r) override;
	double Debye_Screening_Scale_Squared(double r) override;

	unsigned int Target_Count() const override;
	const obscura::Isotope& Target_Isotope(unsigned int index) const override;
	double Number_Density_Nucleus(double r, unsigned int nucleus_index) override;
	double Number_Density_Electron(double r) override;

	double DM_Scattering_Rate_Electron(obscura::DM_Particle& DM, double r, double DM_speed) override;
	double DM_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double r, double DM_speed, unsigned int nucleus_index) override;
	double Total_DM_Scattering_Rate(obscura::DM_Particle& DM, double r, double DM_speed) override;
	double Total_DM_Scattering_Rate_Computed(obscura::DM_Particle& DM, double r, double DM_speed);
	double Total_DM_Scattering_Rate_Interpolated(obscura::DM_Particle& DM, double r, double DM_speed);
	void Interpolate_Total_DM_Scattering_Rate(obscura::DM_Particle& DM, unsigned int N_radius, unsigned int N_speed) override;

	void Print_Summary(int mpi_rank = 0);
};

} // namespace DaMaSCUS_SUN

#endif
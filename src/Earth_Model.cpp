#include "Earth_Model.hpp"
#include "Solar_Model.hpp"  // for Thermal_Averaged_Relative_Speed

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <sstream>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

#include "version.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

namespace
{
std::vector<double> Normalize_Layer_Fractions(std::vector<double> fractions)
{
	double sum = 0.0;
	for(double fraction : fractions)
		sum += fraction;
	if(sum <= 0.0)
		return fractions;
	for(double& fraction : fractions)
		fraction /= sum;
	return fractions;
}
}

Earth_Isotope::Earth_Isotope(const obscura::Isotope& isotope, const std::vector<double>& mass_fractions_by_layer)
: Isotope(isotope), layer_mass_fractions(mass_fractions_by_layer)
{
}

double Earth_Isotope::Mass_Fraction(unsigned int layer_index) const
{
	if(layer_index >= layer_mass_fractions.size())
		return 0.0;
	return layer_mass_fractions[layer_index];
}

Earth_Model::Earth_Model()
: Earth_Model(PROJECT_DIR "/data/earth_prem.dat")
{
}

Earth_Model::Earth_Model(const std::string& path)
: model_file(path), using_interpolated_rate(false), name("Simplified Earth PREM smoke model")
{
	Import_Raw_Data(model_file);
	mass = libphysica::Interpolation(Create_Interpolation_Table(1));
	mass_density = libphysica::Interpolation(Create_Interpolation_Table(2));
	temperature = libphysica::Interpolation(Create_Interpolation_Table(3));
	local_escape_speed_squared = libphysica::Interpolation(Create_Escape_Speed_Table());
	Initialize_Composition();
}

void Earth_Model::Import_Raw_Data(const std::string& path)
{
    std::ifstream input(path.c_str());
    if(!input)
    {
        std::cerr << "Error in Earth_Model::Import_Raw_Data(): Could not open " << path << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    double previous_radius = -1.0;
    raw_data.clear();
    raw_element_ppm.clear();  // 清空元素數據

    while(std::getline(input, line))
    {
        size_t comment = line.find('#');
        if(comment != std::string::npos)
            line = line.substr(0, comment);
        std::istringstream stream(line);
        double r_km = 0.0;
        double enclosed_mass_kg = 0.0;
        double density_g_cm3 = 0.0;
        double temperature_K = 0.0;
        if(!(stream >> r_km >> enclosed_mass_kg >> density_g_cm3 >> temperature_K))
            continue;

        // 讀取所有剩餘的數值（元素數據）
        std::vector<double> element_ppm;
        double value;
        while(stream >> value)
            element_ppm.push_back(value);

        if(r_km < 0.0 || enclosed_mass_kg < 0.0 || density_g_cm3 < 0.0 || temperature_K <= 0.0 || r_km <= previous_radius)
        {
            std::cerr << "Error in Earth_Model::Import_Raw_Data(): Invalid row in " << path << ": " << line << std::endl;
            std::exit(EXIT_FAILURE);
        }

        previous_radius = r_km;
        raw_data.push_back({r_km * km, enclosed_mass_kg * kg, density_g_cm3 * gram / cm / cm / cm, temperature_K * Kelvin});
        raw_element_ppm.push_back(element_ppm);  // 儲存元素數據
    }

    if(raw_data.size() < 2 || std::fabs(raw_data.front()[0]) > 1.0e-12 * km || raw_data.back()[0] < 6371.0 * km)
    {
        std::cerr << "Error in Earth_Model::Import_Raw_Data(): " << path << " must cover 0 <= r_km <= 6371." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // 除錯輸出
    int mpi_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if(mpi_rank == 0)
    {
        std::cout << ">>> Import_Raw_Data: Loaded " << raw_data.size() << " rows." << std::endl;
        std::cout << ">>> raw_element_ppm size: " << raw_element_ppm.size() << std::endl;
        if(!raw_element_ppm.empty())
            std::cout << ">>> First row has " << raw_element_ppm[0].size() << " elements." << std::endl;
    }
}

std::vector<std::vector<double>> Earth_Model::Create_Interpolation_Table(unsigned int column) const
{
	std::vector<std::vector<double>> table(raw_data.size(), std::vector<double>(2, 0.0));
	for(unsigned int i = 0; i < raw_data.size(); i++)
	{
		table[i][0] = raw_data[i][0];
		table[i][1] = raw_data[i][column];
	}
	return table;
}

std::vector<std::vector<double>> Earth_Model::Create_Escape_Speed_Table()
{
	std::vector<std::vector<double>> table(raw_data.size(), std::vector<double>(2, 0.0));
	const double body_radius = Radius();
	const double body_mass = Total_Mass();
	for(unsigned int i = 0; i < raw_data.size(); i++)
	{
		double r = raw_data[i][0];
		auto integrand = [this](double x) {
			if(x <= 0.0)
				return 0.0;
			return Mass(x) / x / x;
		};
		double integral = (r < body_radius) ? libphysica::Integrate(integrand, r, body_radius) : body_mass * (1.0 / r - 1.0 / body_radius);
		table[i][0] = r;
		table[i][1] = 2.0 * G_Newton * body_mass / body_radius * (1.0 + body_radius / body_mass * integral);
	}
	return table;
}

void Earth_Model::Initialize_Composition()
{
    // 如果沒有元素數據，使用硬編碼（向後相容）
    if(raw_element_ppm.empty() || raw_element_ppm[0].empty())
    {
        // ... 保留原有的 8 種硬編碼 ...
        layer_outer_radii = {3480.0 * km, 6346.0 * km, Radius()};
        std::vector<std::vector<double>> layers;
        layers.push_back(Normalize_Layer_Fractions({0.0, 0.0, 0.0, 0.85, 0.05, 0.10, 0.0, 0.0}));
        layers.push_back(Normalize_Layer_Fractions({0.440, 0.228, 0.210, 0.063, 0.0, 0.0, 0.025, 0.034}));
        layers.push_back(Normalize_Layer_Fractions({0.466, 0.021, 0.277, 0.050, 0.001, 0.003, 0.036, 0.081}));
        // ... 後續硬編碼 ...
        return;
    }

    // ============================================
    // 動態讀取所有元素
    // ============================================
    
    int num_elements = raw_element_ppm[0].size();
    int num_layers = raw_data.size();
    
    // 設定層邊界
    layer_outer_radii.clear();
    for(size_t i = 0; i < raw_data.size(); i++)
    {
        layer_outer_radii.push_back(raw_data[i][0]);
    }
    
    // 元素 (Z, A) 映射 - 涵蓋所有可能的元素
    // 由於 raw_element_ppm 只儲存數值，沒有名稱，
    // 我們用索引來對應，但需要知道 Z 和 A
    // 這裡使用通用的方法：Z = 元素索引 + 1（近似）
    // 更好的方法：從 CSV 標題解析
    
    // 手動定義 78 種元素的 (Z, A)（根據 CSV 順序）
    // 這些數值需要與 CSV 中的元素順序對應
    std::vector<std::pair<int, int>> element_ZA = {
        {1, 1},   // H
        {2, 4},   // He
        {3, 7},   // Li
        {4, 9},   // Be
        {5, 11},  // B
        {6, 12},  // C
        {7, 14},  // N
        {8, 16},  // O
        {9, 19},  // F
        {10, 20}, // Ne
        {11, 23}, // Na
        {12, 24}, // Mg
        {13, 27}, // Al
        {14, 28}, // Si
        {15, 31}, // P
        {16, 32}, // S
        {17, 35}, // Cl
        {18, 40}, // Ar
        {19, 39}, // K
        {20, 40}, // Ca
        {21, 45}, // Sc
        {22, 48}, // Ti
        {23, 51}, // V
        {24, 52}, // Cr
        {25, 55}, // Mn
        {26, 56}, // Fe
        {27, 59}, // Co
        {28, 58}, // Ni
        {29, 63}, // Cu
        {30, 64}, // Zn
        {31, 69}, // Ga
        {32, 74}, // Ge
        {33, 75}, // As
        {34, 80}, // Se
        {35, 79}, // Br
        {36, 84}, // Kr
        {37, 85}, // Rb
        {38, 88}, // Sr
        {39, 89}, // Y
        {40, 90}, // Zr
        {41, 93}, // Nb
        {42, 98}, // Mo
        {44, 101}, // Ru
        {45, 103}, // Rh
        {46, 106}, // Pd
        {47, 107}, // Ag
        {48, 112}, // Cd
        {49, 115}, // In
        {50, 118}, // Sn
        {51, 121}, // Sb
        {52, 130}, // Te
        {53, 127}, // I
        {54, 131}, // Xe
        {55, 133}, // Cs
        {56, 137}, // Ba
        {57, 138}, // La
        {58, 140}, // Ce
        {59, 141}, // Pr
        {60, 142}, // Nd
        {62, 147}, // Sm
        {63, 153}, // Eu
        {64, 157}, // Gd
        {65, 159}, // Tb
        {66, 163}, // Dy
        {67, 165}, // Ho
        {68, 166}, // Er
        {69, 169}, // Tm
        {70, 174}, // Yb
        {71, 175}, // Lu
        {72, 178}, // Hf
        {73, 181}, // Ta
        {74, 184}, // W
        {75, 187}, // Re
        {76, 190}, // Os
        {77, 193}, // Ir
        {78, 195}, // Pt
        {79, 197}, // Au
        {80, 200}, // Hg
        {81, 203}, // Tl
        {82, 207}, // Pb
        {83, 209}, // Bi
        {90, 232}, // Th
        {92, 238}, // U
    };
    
    obscura::Import_Nuclear_Data();
    target_isotopes.clear();
    
    // 對每個元素建立同位素
    int num_to_load = std::min(num_elements, (int)element_ZA.size());
    for(int e = 0; e < num_to_load; e++)
    {
        int Z = element_ZA[e].first;
        int A = element_ZA[e].second;
        
        try {
            obscura::Isotope isotope = obscura::Get_Isotope(Z, A);
            isotope.abundance = 1.0;
            
            std::vector<double> mass_fractions;
            for(int l = 0; l < num_layers; l++)
            {
                double ppm = (e < raw_element_ppm[l].size()) ? raw_element_ppm[l][e] : 0.0;
                mass_fractions.push_back(ppm / 1e6);
            }
            target_isotopes.push_back(Earth_Isotope(isotope, mass_fractions));
        } catch(...) {
            // 如果同位素不存在，跳過
            continue;
        }
    }
    
    // 除錯輸出
    int mpi_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if(mpi_rank == 0)
    {
        std::cout << "========================================" << std::endl;
        std::cout << "Earth Composition: " << num_layers << " layers, " 
                  << target_isotopes.size() << " elements loaded." << std::endl;
        std::cout << "========================================" << std::endl;
    }
}


unsigned int Earth_Model::Layer_Index(double r) const
{
    if(r < 0.0)
        r = 0.0;
    
    for(unsigned int layer = 0; layer < layer_outer_radii.size(); layer++)
    {
        if(r <= layer_outer_radii[layer])
            return layer;
    }
    return static_cast<unsigned int>(layer_outer_radii.size()) - 1;
}

const std::string& Earth_Model::Name() const
{
	return name;
}

const std::string& Earth_Model::Model_File() const
{
	return model_file;
}

double Earth_Model::Radius() const
{
	return rEarth;
}

double Earth_Model::Total_Mass() const
{
	return mEarth;
}

double Earth_Model::Mass(double r)
{
	if(r < 0.0)
		r = 0.0;
	if(r > Radius())
		return Total_Mass();
	return mass(r);
}

double Earth_Model::Mass_Density(double r)
{
	if(r < 0.0)
		r = 0.0;
	if(r > Radius())
		return 0.0;
	return mass_density(r);
}

double Earth_Model::Temperature(double r)
{
	if(r < 0.0)
		r = 0.0;
	if(r > Radius())
		r = Radius();
	return temperature(r);
}

double Earth_Model::Local_Escape_Speed(double r)
{
	if(r < 0.0)
		r = 0.0;
	if(r > Radius())
		return sqrt(2.0 * G_Newton * Total_Mass() / r);
	return sqrt(local_escape_speed_squared(r));
}

double Earth_Model::Debye_Screening_Scale_Squared(double r)
{
	if(r > Radius())
		return 0.0;
	double T = Temperature(r);
	if(T <= 0.0)
		return 0.0;
	double debye_scale_2 = 4.0 * M_PI * aEM / T * Number_Density_Electron(r);
	for(unsigned int i = 0; i < Target_Count(); i++)
		debye_scale_2 += 4.0 * M_PI * aEM / T * Number_Density_Nucleus(r, i) * Target_Isotope(i).Z;
	return debye_scale_2;
}

unsigned int Earth_Model::Target_Count() const
{
	return static_cast<unsigned int>(target_isotopes.size());
}

const obscura::Isotope& Earth_Model::Target_Isotope(unsigned int index) const
{
	if(index >= target_isotopes.size())
	{
		std::cerr << "Error in Earth_Model::Target_Isotope(): Index = " << index << " is out of bound (number of targets: " << target_isotopes.size() << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	return target_isotopes[index];
}

double Earth_Model::Number_Density_Nucleus(double r, unsigned int nucleus_index)
{
	if(nucleus_index >= target_isotopes.size())
	{
		std::cerr << "Error in Earth_Model::Number_Density_Nucleus(): Index = " << nucleus_index << " is out of bound (number of targets: " << target_isotopes.size() << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(r < 0.0)
		r = 0.0;
	if(r > Radius())
		return 0.0;
	unsigned int layer = Layer_Index(r);
	return target_isotopes[nucleus_index].Mass_Fraction(layer) * Mass_Density(r) / target_isotopes[nucleus_index].mass;
}

double Earth_Model::Number_Density_Electron(double r)
{
	if(r < 0.0)
		r = 0.0;
	if(r > Radius())
		return 0.0;
	double density = 0.0;
	for(unsigned int i = 0; i < Target_Count(); i++)
		density += Target_Isotope(i).Z * Number_Density_Nucleus(r, i);
	return density;
}

double Earth_Model::DM_Scattering_Rate_Electron(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(r < 0.0)
		r = 0.0;
	if(DM_speed < 0.0)
		DM_speed = 0.0;
	if(r > Radius())
		return 0.0;
	double v_rel = Thermal_Averaged_Relative_Speed(Temperature(r), mElectron, DM_speed);
	return Number_Density_Electron(r) * DM.Sigma_Total_Electron(DM_speed) * v_rel;
}

double Earth_Model::DM_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double r, double DM_speed, unsigned int nucleus_index)
{
	if(nucleus_index >= target_isotopes.size())
	{
		std::cerr << "Error in Earth_Model::DM_Scattering_Rate_Nucleus(): Index = " << nucleus_index << " is out of bound (number of targets: " << target_isotopes.size() << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(r < 0.0)
		r = 0.0;
	if(DM_speed < 0.0)
		DM_speed = 0.0;
	if(r > Radius())
		return 0.0;
	double m_target = target_isotopes[nucleus_index].mass;
	double v_rel = Thermal_Averaged_Relative_Speed(Temperature(r), m_target, DM_speed);
	return Number_Density_Nucleus(r, nucleus_index) * DM.Sigma_Total_Nucleus(target_isotopes[nucleus_index], DM_speed, r) * v_rel;
}

double Earth_Model::Total_DM_Scattering_Rate(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(r < 0.0)
		r = 0.0;
	if(DM_speed < 0.0)
		DM_speed = 0.0;
	if(r > Radius())
		return 0.0;
	if(using_interpolated_rate && DM_speed < rate_interpolation.domain[1][1])
		return Total_DM_Scattering_Rate_Interpolated(DM, r, DM_speed);
	return Total_DM_Scattering_Rate_Computed(DM, r, DM_speed);
}

double Earth_Model::Total_DM_Scattering_Rate_Computed(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(r < 0.0)
		r = 0.0;
	if(DM_speed < 0.0)
		DM_speed = 0.0;
	if(r > Radius())
		return 0.0;
	double total_rate = DM_Scattering_Rate_Electron(DM, r, DM_speed);
	for(unsigned int i = 0; i < Target_Count(); i++)
		total_rate += DM_Scattering_Rate_Nucleus(DM, r, DM_speed, i);
	return total_rate;
}

double Earth_Model::Total_DM_Scattering_Rate_Interpolated(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(r > Radius() || r < 0.0)
		return 0.0;
	if(DM_speed < 0.0)
		return 0.0;
	if(DM_speed > rate_interpolation.domain[1][1])
		return Total_DM_Scattering_Rate_Computed(DM, r, DM_speed);
	return rate_interpolation(r, DM_speed);
}

void Earth_Model::Interpolate_Total_DM_Scattering_Rate(obscura::DM_Particle& DM, unsigned int N_radius, unsigned int N_speed)
{
	// Interpolation_2D must not be constructed with a 1x1 or 2x2 grid.
	// Use direct scattering-rate evaluation for undersized grids.
	if(N_radius < 3 || N_speed < 3)
	{
	        using_interpolated_rate = false;
	        return;
	}

	int mpi_processes = 1;
	int mpi_rank = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	using_interpolated_rate = true;
	double vMax = 0.75;
	unsigned int local_N_radius = std::ceil(1.0 * N_radius / mpi_processes);
	unsigned int global_N_radius = mpi_processes * local_N_radius;

	std::vector<double> global_radii = libphysica::Linear_Space(0, Radius(), global_N_radius);
	std::vector<double> local_radii(local_N_radius, 0.0);
	MPI_Scatter(global_radii.data(), local_N_radius, MPI_DOUBLE, local_radii.data(), local_N_radius, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	std::vector<double> speeds = libphysica::Linear_Space(0, vMax, N_speed);
	std::vector<double> local_rates;
	std::vector<double> global_rates(N_speed * global_N_radius, 0.0);
	for(auto& radius : local_radii)
		for(auto& speed : speeds)
			local_rates.push_back(Total_DM_Scattering_Rate_Computed(DM, radius, speed));
	MPI_Allgather(local_rates.data(), local_N_radius * N_speed, MPI_DOUBLE, global_rates.data(), local_N_radius * N_speed, MPI_DOUBLE, MPI_COMM_WORLD);

	std::vector<std::vector<double>> rates;
	int i = 0;
	for(auto& radius : global_radii)
		for(auto& speed : speeds)
			rates.push_back({radius, speed, global_rates[i++]});
	rate_interpolation = libphysica::Interpolation_2D(rates);
}

void Earth_Model::Print_Summary(int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR
		          << "Earth model:\t\t" << name << std::endl
		          << "Model file:\t\t" << model_file << std::endl
		          << "Radius [km]:\t\t" << In_Units(Radius(), km) << std::endl
		          << "Mass [kg]:\t\t" << In_Units(Total_Mass(), kg) << std::endl
		          << "Surface escape speed [km/s]:\t" << In_Units(Local_Escape_Speed(Radius()), km / sec) << std::endl
		          << "Nuclear targets:\t" << target_isotopes.size() << std::endl
		          << SEPARATOR;
	}
}

} // namespace DaMaSCUS_SUN
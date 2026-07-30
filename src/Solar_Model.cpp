#include "Solar_Model.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <mpi.h>
#include <sstream>
#include <stdexcept>

#if defined(__APPLE__)
#include <mach-o/dyld.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

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

bool File_Is_Readable(const std::string& path)
{
	std::ifstream input(path.c_str());
	return input.good();
}

std::string Join_Path(const std::string& directory, const std::string& path)
{
	if(directory.empty())
		return path;
	if(directory[directory.size() - 1] == '/')
		return directory + path;
	return directory + "/" + path;
}

std::string Directory_Name(const std::string& path)
{
	const std::string::size_type separator = path.find_last_of('/');
	if(separator == std::string::npos)
		return "";
	if(separator == 0)
		return "/";
	return path.substr(0, separator);
}

std::string Running_Executable_Path(const std::string& executable_hint)
{
#if defined(__APPLE__)
	uint32_t size = 0;
	_NSGetExecutablePath(NULL, &size);
	if(size > 0)
	{
		std::vector<char> path(size, '\0');
		if(_NSGetExecutablePath(path.data(), &size) == 0)
			return std::string(path.data());
	}
#elif defined(__linux__)
	std::vector<char> path(4096, '\0');
	const ssize_t length = readlink("/proc/self/exe", path.data(), path.size() - 1);
	if(length > 0)
	{
		path[static_cast<std::size_t>(length)] = '\0';
		return std::string(path.data());
	}
#endif
	return executable_hint;
}

std::runtime_error Missing_Data_Error(const std::vector<std::string>& attempted_paths)
{
	std::ostringstream message;
	message << "Unable to locate model_agss09.dat. Checked:";
	for(std::vector<std::string>::const_iterator path = attempted_paths.begin(); path != attempted_paths.end(); ++path)
		message << "\n  - " << *path;
	message << "\nSet DAMASCUS_SUN_SOLAR_MODEL to the data file, or "
	        << "DAMASCUS_SUN_DATA_DIR to its containing directory.";
	return std::runtime_error(message.str());
}

} // namespace

std::string Locate_Solar_Model_Data_File(const std::string& executable_hint)
{
	std::vector<std::string> attempted_paths;

	const char* explicit_file = std::getenv("DAMASCUS_SUN_SOLAR_MODEL");
	if(explicit_file != NULL && explicit_file[0] != '\0')
	{
		const std::string path(explicit_file);
		if(File_Is_Readable(path))
			return path;
		attempted_paths.push_back(path + " (DAMASCUS_SUN_SOLAR_MODEL)");
		throw Missing_Data_Error(attempted_paths);
	}

	const char* explicit_directory = std::getenv("DAMASCUS_SUN_DATA_DIR");
	if(explicit_directory != NULL && explicit_directory[0] != '\0')
	{
		const std::string path = Join_Path(explicit_directory, "model_agss09.dat");
		if(File_Is_Readable(path))
			return path;
		attempted_paths.push_back(path + " (DAMASCUS_SUN_DATA_DIR)");
		throw Missing_Data_Error(attempted_paths);
	}

	const std::string executable_path = Running_Executable_Path(executable_hint);
	std::string executable_directory = Directory_Name(executable_path);
	// Multi-config generators may add one configuration directory below the
	// normal build output. Check a small number of ancestors while retaining the
	// standard <prefix>/bin -> <prefix>/share lookup as the first candidate.
	for(unsigned int level = 0; level < 3 && !executable_directory.empty(); ++level)
	{
		const std::string installed_path = Join_Path(executable_directory, INSTALL_DATA_PATH_FROM_BINDIR);
		attempted_paths.push_back(installed_path);
		if(File_Is_Readable(installed_path))
			return installed_path;
		const std::string parent_directory = Directory_Name(executable_directory);
		if(parent_directory == executable_directory)
			break;
		executable_directory = parent_directory;
	}

	const std::string working_directory_path = "data/model_agss09.dat";
	attempted_paths.push_back(working_directory_path);
	if(File_Is_Readable(working_directory_path))
		return working_directory_path;

	throw Missing_Data_Error(attempted_paths);
}

// 1. Nuclear targets in the Sun
Solar_Isotope::Solar_Isotope(const obscura::Isotope& isotope, const std::vector<std::vector<double>>& density_table, double abundance)
: Isotope(isotope), number_density(libphysica::Interpolation(density_table))
{
	number_density.Multiply(abundance);
}

double Solar_Isotope::Number_Density(double r)
{
	if(!std::isfinite(r))
		return std::numeric_limits<double>::quiet_NaN();
	// 边界检查：防止负半径
	if(r < 0.0)
		r = 0.0;
	if(r > rSun)
		return 0.0;
	else
		return number_density(r);
}

// 2. Solar model
// Auxiliary functions for the data import
void Solar_Model::Import_Raw_Data(const std::string& data_file)
{
	std::vector<double> units = {mSun, rSun, Kelvin, gram / cm / cm / cm, dyne / cm / cm, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	raw_data				  = libphysica::Import_Table(data_file, units, 20);
	// Add a first and last line to ensure a full domain of the interpolations.
	std::vector<double> first_line(35, 0.0);
	std::vector<double> last_line(35, 0.0);
	for(unsigned int i = 0; i < 35; i++)
	{
		first_line[i] = (i == 0 || i == 1) ? 0.0 : raw_data.front()[i];
		if(i == 0)
			last_line[i] = mSun;
		else if(i == 1)
			last_line[i] = rSun;
		else if(i == 2)
			last_line[i] = 5800 * Kelvin;	// temperature of the photosphere (http://solar-center.stanford.edu/vitalstats.html)
		else if(i == 3)
			last_line[i] = 1.0e-9 * gram / cm / cm / cm;   // mass density of the photosphere (http://solar-center.stanford.edu/vitalstats.html)
		else
			last_line[i] = raw_data.back()[i];
	}
	raw_data.insert(raw_data.begin(), first_line);
	raw_data.push_back(last_line);
}

std::vector<std::vector<double>> Solar_Model::Create_Interpolation_Table(unsigned int row) const
{
	std::vector<std::vector<double>> table(raw_data.size(), std::vector<double>(2, 0.0));
	for(unsigned int i = 0; i < table.size(); i++)
	{
		table[i][0] = raw_data[i][1];
		table[i][1] = raw_data[i][row - 1];
	}
	return table;
}

std::vector<std::vector<double>> Solar_Model::Create_Escape_Speed_Table()
{
	std::vector<std::vector<double>> table_vesc(raw_data.size(), std::vector<double>(2, 0.0));
	for(unsigned int i = 0; i < raw_data.size(); i++)
	{
		double r	   = raw_data[i][1];
		auto integrand = [this](double x) {
			if(x == 0)
				return 0.0;
			else
				return Mass(x) / x / x;
		};
		double integral	 = (r < rSun) ? libphysica::Integrate(integrand, r, rSun) : mSun * (1.0 / r - 1.0 / rSun);
		table_vesc[i][0] = r;
		table_vesc[i][1] = 2.0 * G_Newton * mSun / rSun * (1.0 + rSun / mSun * integral);
	}
	return table_vesc;
}

std::vector<std::vector<double>> Solar_Model::Create_Number_Density_Table(unsigned int target, double mass) const
{
	std::vector<std::vector<double>> table;
	for(unsigned int i = 0; i < raw_data.size(); i++)
	{
		double r   = raw_data[i][1];
		double rho = raw_data[i][3];
		double f   = raw_data[i][6 + target];
		double n   = f * rho / mass;
		table.push_back({r, n});
	}
	return table;
}

std::vector<std::vector<double>> Solar_Model::Create_Number_Density_Table_Electron()
{
	std::vector<std::vector<double>> table(raw_data.size(), std::vector<double>(2, 0.0));
	for(auto& isotope : target_isotopes)
	{
		for(unsigned int i = 0; i < raw_data.size(); i++)
		{
			double r	= raw_data[i][1];
			double n	= isotope.Number_Density(r);
			table[i][0] = r;
			table[i][1] += isotope.Z * n;
		}
	}
	return table;
}

Solar_Model::Solar_Model(const std::string& data_file)
: using_interpolated_rate(false),
  rate_grid_radius_points(0), rate_grid_speed_points(0),
  rate_grid_inverse_radius_step(0.0), rate_grid_inverse_speed_step(0.0),
  rate_grid_max_speed(0.0), name("Standard Solar Model AGSS09")
{
	Import_Raw_Data(data_file.empty() ? Locate_Solar_Model_Data_File() : data_file);

	// Interpolate tables.
	mass					   = libphysica::Interpolation(Create_Interpolation_Table(1));
	temperature				   = libphysica::Interpolation(Create_Interpolation_Table(3));
	local_escape_speed_squared = libphysica::Interpolation(Create_Escape_Speed_Table());
	mass_density			   = libphysica::Interpolation(Create_Interpolation_Table(4));

	// Nuclear abundances
	obscura::Import_Nuclear_Data();
	std::vector<int> Zs				  = {1, 2, 2, 6, 6, 7, 7, 8, 8, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
	std::vector<double> As			  = {1.0, 4.0, 3.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0};
	std::vector<int> included_targets = libphysica::Range(Zs.size());
	for(auto& target_index : included_targets)
	{
		int Z = Zs[target_index];
		if(target_index < 10)
		{
			double A				 = As[target_index];
			obscura::Isotope isotope = obscura::Get_Isotope(Z, A);
			isotope.abundance		 = 1.0;
			target_isotopes.push_back(Solar_Isotope(isotope, Create_Number_Density_Table(target_index, isotope.mass)));
		}
		else
		{
			obscura::Nucleus nucleus = obscura::Get_Nucleus(Z);
			for(const auto& isotope : nucleus.isotopes)
				target_isotopes.push_back(Solar_Isotope(isotope, Create_Number_Density_Table(target_index, isotope.mass), isotope.abundance));
		}
	}
	// Electron number density
	number_density_electron = libphysica::Interpolation(Create_Number_Density_Table_Electron());
}

const std::string& Solar_Model::Name() const
{
    return name;
}
double Solar_Model::Radius() const
{
    return rSun;
}
double Solar_Model::Total_Mass() const
{
    return mSun;
}
unsigned int Solar_Model::Target_Count() const
{
    return target_isotopes.size();
}
const obscura::Isotope& Solar_Model::Target_Isotope(unsigned int index) const
{
    return target_isotopes[index];
}

double Solar_Model::Mass(double r)
{
	if(!std::isfinite(r))
		return std::numeric_limits<double>::quiet_NaN();
	// 边界检查：防止负半径
	if(r < 0.0)
		r = 0.0;
	if(r > rSun)
		return mSun;
	else
		return mass(r);
}

double Solar_Model::Mass_Density(double r)
{
	if(!std::isfinite(r))
		return std::numeric_limits<double>::quiet_NaN();
	// 边界检查：防止负半径
	if(r < 0.0)
		r = 0.0;
	if(r > rSun)
		return 0.0;
	else
		return mass_density(r);
}

double Solar_Model::Temperature(double r)
{
	if(!std::isfinite(r))
		return std::numeric_limits<double>::quiet_NaN();
	// 边界检查：防止负半径
	if(r < 0.0)
		r = 0.0;
	return temperature(r);
}

double Solar_Model::Local_Escape_Speed(double r)
{
	if(!std::isfinite(r))
		return std::numeric_limits<double>::quiet_NaN();
	// 边界检查：防止负半径导致插值错误
	if(r < 0.0)
	{
		std::cerr << "Warning Solar_Model::Local_Escape_Speed(): Negative radius r = " << r << ", using r = 0." << std::endl;
		r = 0.0;
	}
	if(r > rSun)
		return sqrt(2 * G_Newton * mSun / r);
	else
		return sqrt(local_escape_speed_squared(r));
}

double Solar_Model::Debye_Screening_Scale_Squared(double r)
{
	if(!std::isfinite(r))
		return std::numeric_limits<double>::quiet_NaN();
	if(r <= rSun)
	{
		double T			 = Temperature(r);
		double debye_scale_2 = 4.0 * M_PI * aEM / T * Number_Density_Electron(r);
		for(unsigned int i = 0; i < target_isotopes.size(); i++)
			debye_scale_2 += 4.0 * M_PI * aEM / T * Number_Density_Nucleus(r, i) * target_isotopes[i].Z;
		return debye_scale_2;
	}
	else
	{
		std::cerr << "Error in Solar_Model::Debye_Screening_Scale(): r/rSun = " << r / rSun << " is outside the Sun." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

double Solar_Model::Number_Density_Nucleus(double r, unsigned int nucleus_index)
{
	if(nucleus_index >= target_isotopes.size())
	{
		std::cerr << "Error in Solar_Model::Number_Density_Nucleus(): Index = " << nucleus_index << " is out of bound (number of targets: " << target_isotopes.size() << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
		return target_isotopes[nucleus_index].Number_Density(r);
}

double Solar_Model::Number_Density_Electron(double r)
{
	if(!std::isfinite(r))
		return std::numeric_limits<double>::quiet_NaN();
	// 边界检查：防止负半径
	if(r < 0.0)
		r = 0.0;
	if(r > rSun)
		return 0.0;
	else
		return number_density_electron(r);
}

double Solar_Model::DM_Scattering_Rate_Electron(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(!std::isfinite(r) || !std::isfinite(DM_speed))
		return std::numeric_limits<double>::quiet_NaN();
	// 边界检查：防止负半径或负速度
	if(r < 0.0)
		r = 0.0;
	if(DM_speed < 0.0)
		DM_speed = 0.0;
	if(r > rSun)
		return 0.0;
	else
	{
		double v_rel = Thermal_Averaged_Relative_Speed(Temperature(r), mElectron, DM_speed);
		return Number_Density_Electron(r) * DM.Sigma_Total_Electron(DM_speed) * v_rel;
	}
}

double Solar_Model::DM_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double r, double DM_speed, unsigned int nucleus_index)
{
	if(nucleus_index >= target_isotopes.size())
	{
		std::cerr << "Error in Solar_Model::Number_Density_Nucleus(): Index = " << nucleus_index << " is out of bound (number of targets: " << target_isotopes.size() << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(!std::isfinite(r) || !std::isfinite(DM_speed))
		return std::numeric_limits<double>::quiet_NaN();
	// 边界检查：防止负半径或负速度
	if(r < 0.0)
		r = 0.0;
	if(DM_speed < 0.0)
		DM_speed = 0.0;
	if(r > rSun)
		return 0.0;
	else
	{
		double m_target = target_isotopes[nucleus_index].mass;
		double v_rel	= Thermal_Averaged_Relative_Speed(Temperature(r), m_target, DM_speed);
		return Number_Density_Nucleus(r, nucleus_index) * DM.Sigma_Total_Nucleus(target_isotopes[nucleus_index], DM_speed, r) * v_rel;
	}
}

double Solar_Model::Total_DM_Scattering_Rate(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(!std::isfinite(r) || !std::isfinite(DM_speed))
		return std::numeric_limits<double>::quiet_NaN();
	// 边界检查：防止负半径或负速度
	if(r < 0.0)
		r = 0.0;
	if(DM_speed < 0.0)
		DM_speed = 0.0;
	if(using_interpolated_rate && DM_speed <= rate_grid_max_speed)
		return Total_DM_Scattering_Rate_Interpolated(DM, r, DM_speed);
	else
	{
		if(using_interpolated_rate)
			std::cerr << "Warning Solar_Model::Total_DM_Scattering_Rate(): DM speed is out of bound (vDM = " << DM_speed << ")\n\tScattering rate must be computed on the fly." << std::endl;
		return Total_DM_Scattering_Rate_Computed(DM, r, DM_speed);
	}
}

double Solar_Model::Total_DM_Scattering_Rate_Computed(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(!std::isfinite(r) || !std::isfinite(DM_speed))
		return std::numeric_limits<double>::quiet_NaN();
	// 边界检查：防止负半径或负速度
	if(r < 0.0)
		r = 0.0;
	if(DM_speed < 0.0)
		DM_speed = 0.0;
	if(r > rSun)
		return 0.0;
	else
	{
		double total_rate = DM_Scattering_Rate_Electron(DM, r, DM_speed);
		for(unsigned int i = 0; i < target_isotopes.size(); i++)
			total_rate += DM_Scattering_Rate_Nucleus(DM, r, DM_speed, i);
		return total_rate;
	}
}

double Solar_Model::Total_DM_Scattering_Rate_Interpolated(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(!std::isfinite(r) || !std::isfinite(DM_speed))
		return std::numeric_limits<double>::quiet_NaN();
	const std::size_t expected_grid_size = static_cast<std::size_t>(rate_grid_radius_points)
	                                       * static_cast<std::size_t>(rate_grid_speed_points);
	if(!using_interpolated_rate || rate_grid_radius_points < 2 || rate_grid_speed_points < 2
	   || rate_grid.size() != expected_grid_size
	   || !std::isfinite(rate_grid_inverse_radius_step) || rate_grid_inverse_radius_step <= 0.0
	   || !std::isfinite(rate_grid_inverse_speed_step) || rate_grid_inverse_speed_step <= 0.0)
		return Total_DM_Scattering_Rate_Computed(DM, r, DM_speed);
	// 边界检查：防止超出插值域
	if(r > rSun || r < 0.0)
		return 0.0;
	if(DM_speed < 0.0)
	{
		std::cerr << "Warning Solar_Model::Total_DM_Scattering_Rate_Interpolated(): Negative DM speed (" << DM_speed << "), returning 0." << std::endl;
		return 0.0;
	}
	// 确保速度不超过插值上限
	if(DM_speed > rate_grid_max_speed)
	{
		std::cerr << "Warning Solar_Model::Total_DM_Scattering_Rate_Interpolated(): DM speed (" << DM_speed << ") exceeds interpolation domain, using computed rate." << std::endl;
		return Total_DM_Scattering_Rate_Computed(DM, r, DM_speed);
	}

	const double radius_coordinate = r * rate_grid_inverse_radius_step;
	const double speed_coordinate = DM_speed * rate_grid_inverse_speed_step;
	if(!std::isfinite(radius_coordinate) || !std::isfinite(speed_coordinate)
	   || radius_coordinate < 0.0 || speed_coordinate < 0.0)
		return Total_DM_Scattering_Rate_Computed(DM, r, DM_speed);
	const unsigned int radius_index = std::min(
	    static_cast<unsigned int>(radius_coordinate), rate_grid_radius_points - 2);
	const unsigned int speed_index = std::min(
	    static_cast<unsigned int>(speed_coordinate), rate_grid_speed_points - 2);
	const double radius_fraction = radius_coordinate - radius_index;
	const double speed_fraction = speed_coordinate - speed_index;
	const size_t row_stride = rate_grid_speed_points;
	const size_t lower_offset = static_cast<size_t>(radius_index) * row_stride + speed_index;
	const size_t upper_offset = lower_offset + row_stride;

	const double lower_rate = rate_grid[lower_offset]
	                        + speed_fraction * (rate_grid[lower_offset + 1] - rate_grid[lower_offset]);
	const double upper_rate = rate_grid[upper_offset]
	                        + speed_fraction * (rate_grid[upper_offset + 1] - rate_grid[upper_offset]);
	const double interpolated_rate = lower_rate + radius_fraction * (upper_rate - lower_rate);
	return std::isfinite(interpolated_rate)
	           ? std::max(0.0, interpolated_rate)
	           : std::numeric_limits<double>::quiet_NaN();
}

void Solar_Model::Interpolate_Total_DM_Scattering_Rate(obscura::DM_Particle& DM, unsigned int N_radius, unsigned int N_speed)
{
	if(N_radius < 2 || N_speed < 2)
	{
		using_interpolated_rate = false;
		rate_grid.clear();
		rate_grid_radius_points = 0;
		rate_grid_speed_points = 0;
		rate_grid_inverse_radius_step = 0.0;
		rate_grid_inverse_speed_step = 0.0;
		rate_grid_max_speed = 0.0;
	}
	else
	{
		int mpi_processes, mpi_rank;
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

		const double vMax = 0.75;
		const unsigned int process_count = static_cast<unsigned int>(mpi_processes);
		const unsigned int base_radius_count = N_radius / process_count;
		const unsigned int remainder = N_radius % process_count;
		const unsigned int local_N_radius = base_radius_count + (static_cast<unsigned int>(mpi_rank) < remainder ? 1u : 0u);
		const unsigned int local_radius_offset = static_cast<unsigned int>(mpi_rank) * base_radius_count
		                                         + std::min(static_cast<unsigned int>(mpi_rank), remainder);

		const unsigned long long total_rate_count = static_cast<unsigned long long>(N_radius) * N_speed;
		if(total_rate_count > static_cast<unsigned long long>(std::numeric_limits<int>::max()))
			throw std::invalid_argument("scattering-rate interpolation grid is too large for MPI_Allgatherv");

		std::vector<int> receive_counts(mpi_processes, 0);
		std::vector<int> displacements(mpi_processes, 0);
		for(int rank = 0; rank < mpi_processes; rank++)
		{
			const unsigned int rank_index = static_cast<unsigned int>(rank);
			const unsigned int rank_radius_count = base_radius_count + (rank_index < remainder ? 1u : 0u);
			const unsigned int rank_radius_offset = rank_index * base_radius_count + std::min(rank_index, remainder);
			receive_counts[rank] = static_cast<int>(static_cast<unsigned long long>(rank_radius_count) * N_speed);
			displacements[rank] = static_cast<int>(static_cast<unsigned long long>(rank_radius_offset) * N_speed);
		}

		// Each rank computes its exact contiguous slice of the requested grid. Uneven
		// slices avoid padding N_radius to a multiple of the MPI process count.
		const std::vector<double> speeds = libphysica::Linear_Space(0, vMax, N_speed);
		std::vector<double> local_rates;
		local_rates.reserve(static_cast<size_t>(local_N_radius) * N_speed);
		for(unsigned int local_radius_index = 0; local_radius_index < local_N_radius; local_radius_index++)
		{
			const unsigned int global_radius_index = local_radius_offset + local_radius_index;
			const double radius = rSun * static_cast<double>(global_radius_index) / static_cast<double>(N_radius - 1);
			for(const auto speed : speeds)
				local_rates.push_back(Total_DM_Scattering_Rate_Computed(DM, radius, speed));
		}

		std::vector<double> global_rates(static_cast<size_t>(N_speed) * N_radius, 0.0);
		MPI_Allgatherv(local_rates.empty() ? nullptr : local_rates.data(),
		               receive_counts[mpi_rank], MPI_DOUBLE,
		               global_rates.data(), receive_counts.data(), displacements.data(),
		               MPI_DOUBLE, MPI_COMM_WORLD);

		// The table is regular and radius-major, so queries can use direct indices
		// instead of two generic binary searches and nested vector lookups.
		rate_grid.swap(global_rates);
		rate_grid_radius_points = N_radius;
		rate_grid_speed_points = N_speed;
		rate_grid_inverse_radius_step = static_cast<double>(N_radius - 1) / rSun;
		rate_grid_inverse_speed_step = static_cast<double>(N_speed - 1) / vMax;
		rate_grid_max_speed = vMax;
		using_interpolated_rate = true;
	}
}

unsigned int Solar_Model::Scattering_Rate_Interpolation_Radius_Points() const
{
	return rate_grid_radius_points;
}

unsigned int Solar_Model::Scattering_Rate_Interpolation_Speed_Points() const
{
	return rate_grid_speed_points;
}

void Solar_Model::Print_Summary(int mpi_rank) const
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR
				  << "Solar model:\t\t" << name << std::endl
				  << "Nuclear targets:\t" << target_isotopes.size() << std::endl
				  << std::endl
				  << "Isotope\tZ\tA\tAbund.[%]\tSpin\t<sp>\t<sn>"
				  << SEPARATOR_LINE;
		for(auto& isotope : target_isotopes)
			isotope.Print_Summary(mpi_rank);
		std::cout << SEPARATOR;
	}
}

double Thermal_Averaged_Relative_Speed(double temperature, double mass_target, double v_DM)
{
	if(!std::isfinite(temperature) || temperature <= 0.0
	   || !std::isfinite(mass_target) || mass_target <= 0.0
	   || !std::isfinite(v_DM) || v_DM < 0.0)
		return std::numeric_limits<double>::quiet_NaN();

	const double kappa = sqrt(mass_target / (2.0 * temperature));
	if(!std::isfinite(kappa) || kappa <= 0.0)
		return std::numeric_limits<double>::quiet_NaN();
	if(v_DM == 0.0)
		return 2.0 / sqrt(M_PI) / kappa;

	const double x = kappa * v_DM;
	if(!std::isfinite(x))
		return v_DM;
	if(x < 1.0e-4)
	{
		const double x2 = x * x;
		return 2.0 / sqrt(M_PI) / kappa
		       * (1.0 + x2 / 3.0 - x2 * x2 / 30.0 + x2 * x2 * x2 / 210.0);
	}
	if(x > 26.0)
		return v_DM * (1.0 + 0.5 / (x * x));

	const double x2 = x * x;
	return v_DM * ((1.0 + 2.0 * x2) * erf(x) / (2.0 * x2)
	               + exp(-x2) / (sqrt(M_PI) * x));
}

}	// namespace DaMaSCUS_SUN

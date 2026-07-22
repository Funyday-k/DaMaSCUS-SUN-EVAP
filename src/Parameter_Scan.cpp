#include "Parameter_Scan.hpp"

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <libconfig.h++>
#include <limits>
#include <mpi.h>
#include <set>
#include <stdexcept>
#include <string>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "Dark_Photon.hpp"
#include "Data_Generation.hpp"
#include "Reflection_Spectrum.hpp"

std::string g_top_level_dir;  // 全局输出目录，从config文件读取
unsigned int g_max_trajectories = 0;  // 最大轨迹数安全阀；0 表示不限制

namespace DaMaSCUS_SUN
{

using namespace libconfig;
using namespace libphysica::natural_units;

namespace
{

double Checked_Probability(double p, const std::string& context)
{
	if(!std::isfinite(p) || p < 0.0 || p > 1.0)
		throw std::runtime_error(context + ": p-value must be finite and lie in [0, 1].");
	return p;
}

double Log10_Probability(double p, const std::string& context)
{
	p = Checked_Probability(p, context);
	return (p < 1.0e-100) ? -100.0 : log10(p);
}

std::size_t STA_Step_Limit(std::size_t rows, std::size_t columns)
{
	if(rows == 0 || columns == 0)
		throw std::invalid_argument("square tracing requires a non-empty parameter grid");
	if(rows > std::numeric_limits<std::size_t>::max() / columns)
		throw std::length_error("parameter grid is too large for square tracing");

	const std::size_t cells = rows * columns;
	if(cells > (std::numeric_limits<std::size_t>::max() - 16) / 16)
		throw std::length_error("parameter grid is too large for a bounded square trace");
	return std::max<std::size_t>(32, 16 * cells + 16);
}

} // namespace

// 1. Configuration class for input file, which extends the obscura::Configuration class.

Configuration::Configuration(std::string cfg_filename, int MPI_rank)
{
	cfg_file	 = cfg_filename;
	results_path = "./";

	// 1. Read the cfg file.
	Read_Config_File();

	// 2. Find the run ID, create a folder and copy the cfg file.
	Initialize_Result_Folder(MPI_rank);

	// 3. DM particle
	Construct_DM_Particle();

	// 4. DM Distribution
	Construct_DM_Distribution();

	// 5. DM-detection experiment
	Construct_DM_Detector();

	// 6. Computation of exclusion limits
	Initialize_Parameters();

	// 7. DaMaSCUS specific parameters
	Import_Parameter_Scan_Parameter();
}

void Configuration::Import_Parameter_Scan_Parameter()
{
	try
	{
		const int configured_sample_size = config.lookup("sample_size");
		if(configured_sample_size <= 0)
		{
			std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): 'sample_size' must be positive." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		sample_size = static_cast<unsigned int>(configured_sample_size);
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'sample_size' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	fixed_seed = 0;
	try
	{
		const int configured_seed = config.lookup("fixed_seed");
		if(configured_seed < 0)
		{
			std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): 'fixed_seed' must be non-negative." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		fixed_seed = static_cast<unsigned int>(configured_seed);
	}
	catch(const SettingNotFoundException& nfex)
	{
	}
	try
	{
		run_mode = config.lookup("run_mode").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'run_mode' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	capture_mode = (run_mode == "Capture");
	try
	{
		bool cfg_capture_mode = config.lookup("capture_mode");
		capture_mode = cfg_capture_mode;
	}
	catch(const SettingNotFoundException& nfex)
	{
	}
	if(run_mode == "Capture")
		capture_mode = true;
	const bool parameter_scan_mode = (run_mode == "Parameter scan");
	try
	{
		const int configured_rings = config.lookup("isoreflection_rings");
		if(configured_rings <= 0)
		{
			std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): 'isoreflection_rings' must be positive." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		isoreflection_rings = static_cast<unsigned int>(configured_rings);
	}
	catch(const SettingNotFoundException& nfex)
	{
		isoreflection_rings = 1;
	}
	try
	{
		const int configured_interpolation_points = config.lookup("interpolation_points");
		if(configured_interpolation_points < 0)
		{
			std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): 'interpolation_points' must be non-negative." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		interpolation_points = static_cast<unsigned int>(configured_interpolation_points);
	}
	catch(const SettingNotFoundException& nfex)
	{
		interpolation_points = 0;
	}
	try
	{
		cross_section_min = config.lookup("cross_section_min");
		cross_section_min *= cm * cm;
	}
	catch(const SettingNotFoundException& nfex)
	{
		if(parameter_scan_mode)
		{
			std::cerr << "No 'cross_section_min' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		cross_section_min = 0.0;
	}

	try
	{
		cross_section_max = config.lookup("cross_section_max");
		cross_section_max *= cm * cm;
	}
	catch(const SettingNotFoundException& nfex)
	{
		if(parameter_scan_mode)
		{
			std::cerr << "No 'cross_section_max' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		cross_section_max = 0.0;
	}

	try
	{
		const int configured_cross_sections = config.lookup("cross_sections");
		if(parameter_scan_mode && configured_cross_sections <= 0)
		{
			std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): 'cross_sections' must be positive in parameter-scan mode." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		if(configured_cross_sections < 0)
		{
			std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): 'cross_sections' must be non-negative." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		cross_sections = static_cast<unsigned int>(configured_cross_sections);
	}
	catch(const SettingNotFoundException& nfex)
	{
		if(parameter_scan_mode)
		{
			std::cerr << "No 'cross_sections' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		cross_sections = 0;
	}

	try
	{
		compute_halo_constraints = config.lookup("compute_halo_constraints");
	}
	catch(const SettingNotFoundException& nfex)
	{
		if(parameter_scan_mode)
		{
			std::cerr << "No 'compute_halo_constraints' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		compute_halo_constraints = false;
	}
	try
	{
		perform_full_scan = config.lookup("perform_full_scan");
	}
	catch(const SettingNotFoundException& nfex)
	{
		if(parameter_scan_mode)
		{
			std::cerr << "No 'perform_full_scan' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		perform_full_scan = false;
	}
	try
	{
		g_top_level_dir = config.lookup("output_dir").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		g_top_level_dir = "./";
	}
	if(g_top_level_dir.empty())
		g_top_level_dir = "./";
	else if(g_top_level_dir[g_top_level_dir.size() - 1] != '/')
		g_top_level_dir += "/";
	try
	{
		int mt = config.lookup("max_trajectories");
		if(mt < 0)
		{
			std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): 'max_trajectories' must be non-negative." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		max_trajectories = static_cast<unsigned int>(mt);
		g_max_trajectories = max_trajectories;
	}
	catch(const SettingNotFoundException& nfex)
	{
		// Optional: no implicit trajectory cap. sample_size is the target
		// number of captured particles, so a sample_size-dependent cap can
		// silently terminate a low-capture-rate run before that target is met.
		max_trajectories = 0;
		g_max_trajectories = max_trajectories;
	}

	// Maximum number of scatterings/collisions per trajectory. Optional in cfg.
	maximum_number_of_scatterings = DEFAULT_MAXIMUM_SCATTERINGS;
	try
	{
		config.lookup("maximum_number_of_scatterings");
		unsigned long long max_scatterings = 0;
		if(config.lookupValue("maximum_number_of_scatterings", max_scatterings))
		{
			if(max_scatterings > std::numeric_limits<unsigned long int>::max())
			{
				std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): 'maximum_number_of_scatterings' is too large." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			maximum_number_of_scatterings = static_cast<unsigned long int>(max_scatterings);
		}
		else
		{
			double max_scatterings_float = 0.0;
			if(!config.lookupValue("maximum_number_of_scatterings", max_scatterings_float) || !std::isfinite(max_scatterings_float) || max_scatterings_float < 0.0 || std::floor(max_scatterings_float) != max_scatterings_float || max_scatterings_float > static_cast<double>(std::numeric_limits<unsigned long int>::max()))
			{
				std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): 'maximum_number_of_scatterings' must be a non-negative integer." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			maximum_number_of_scatterings = static_cast<unsigned long int>(max_scatterings_float);
		}
	}
	catch(const SettingNotFoundException& nfex)
	{
	}

	// Snapshot configuration (optional)
	try
	{
		bool snap_enabled = config.lookup("snapshot_enabled");
		snapshot_config.enabled = snap_enabled;
	}
	catch(const SettingNotFoundException& nfex)
	{
		snapshot_config.enabled = false;
	}
	if(snapshot_config.enabled)
	{
		try
		{
			double interval = 0.0;
			if(!config.lookupValue("snapshot_interval", interval))
			{
				long long integer_interval = 0;
				if(!config.lookupValue("snapshot_interval", integer_interval))
				{
					std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): "
					          << "'snapshot_interval' must be a numeric value." << std::endl;
					std::exit(EXIT_FAILURE);
				}
				interval = static_cast<double>(integer_interval);
			}
			if(!IsValidSnapshotIntervalSeconds(interval))
			{
				std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): "
				          << "'snapshot_interval' must be a positive integer number of seconds." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			snapshot_config.interval_seconds = interval;
		}
		catch(const SettingNotFoundException& nfex)
		{
			snapshot_config.interval_seconds = 60.0;  // default: 60 seconds
		}
	}
	// Per-trajectory wall-time budget (seconds). 0 = unlimited. Optional in cfg.
	try
	{
		double wall_budget = config.lookup("max_trajectory_wall_time_sec");
		if(!std::isfinite(wall_budget) || wall_budget < 0.0)
		{
			std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): 'max_trajectory_wall_time_sec' must be finite and non-negative." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		snapshot_config.max_trajectory_wall_time_sec = wall_budget;
	}
	catch(const SettingNotFoundException& nfex)
	{
		snapshot_config.max_trajectory_wall_time_sec = 0.0;  // default: unlimited
	}

	trajectory_diagnostic_config = TrajectoryDiagnosticConfig();
	trajectory_diagnostic_config.interpolation_points = interpolation_points;
	try
	{
		trajectory_diagnostic_config.summary_enabled = config.lookup("trajectory_summary_enabled");
	}
	catch(const SettingNotFoundException& nfex)
	{
	}
	try
	{
		trajectory_diagnostic_config.events_enabled = config.lookup("trajectory_events_enabled");
	}
	catch(const SettingNotFoundException& nfex)
	{
	}
	try
	{
		trajectory_diagnostic_config.trace_rate = config.lookup("trajectory_trace_rate");
	}
	catch(const SettingNotFoundException& nfex)
	{
		trajectory_diagnostic_config.trace_rate = trajectory_diagnostic_config.events_enabled ? 0.02 : 0.0;
	}
	try
	{
		const long long trace_seed = config.lookup("trajectory_trace_seed");
		if(trace_seed < 0)
		{
			std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): 'trajectory_trace_seed' must be a non-negative integer." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		trajectory_diagnostic_config.trace_seed = static_cast<uint64_t>(trace_seed);
	}
	catch(const SettingNotFoundException& nfex)
	{
		trajectory_diagnostic_config.trace_seed = static_cast<uint64_t>(fixed_seed);
	}
	if(!std::isfinite(trajectory_diagnostic_config.trace_rate)
	   || trajectory_diagnostic_config.trace_rate < 0.0
	   || trajectory_diagnostic_config.trace_rate > 1.0)
	{
		std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): 'trajectory_trace_rate' must lie in [0, 1]." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(trajectory_diagnostic_config.events_enabled)
		trajectory_diagnostic_config.summary_enabled = true;

	if(run_mode != "Parameter point" && run_mode != "Parameter scan" && run_mode != "Custom" && run_mode != "Capture")
	{
		std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): Run mode " << run_mode << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(parameter_scan_mode)
	{
		const int configured_mass_points = config.lookup("constraints_masses");
		if(!std::isfinite(constraints_certainty) || constraints_certainty <= 0.0 || constraints_certainty >= 1.0
		   || !std::isfinite(constraints_mass_min) || constraints_mass_min <= 0.0
		   || !std::isfinite(constraints_mass_max) || constraints_mass_max <= constraints_mass_min
		   || configured_mass_points <= 0
		   || !std::isfinite(cross_section_min) || cross_section_min <= 0.0
		   || !std::isfinite(cross_section_max) || cross_section_max <= cross_section_min)
		{
			std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): invalid or non-finite parameter-scan bounds." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		constraints_masses = static_cast<unsigned int>(configured_mass_points);
	}
}

void Configuration::Construct_DM_Particle()
{
	double DM_mass, DM_spin, DM_fraction;
	bool DM_light;
	// 3.1 General properties
	try
	{
		DM_mass = config.lookup("DM_mass");
		DM_mass *= GeV;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_mass' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	try
	{
		DM_spin = config.lookup("DM_spin");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_spin' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	try
	{
		DM_fraction = config.lookup("DM_fraction");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_fraction' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	try
	{
		DM_light = config.lookup("DM_light");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_light' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(!std::isfinite(DM_mass) || DM_mass <= 0.0
	   || !std::isfinite(DM_spin) || DM_spin < 0.0
	   || !std::isfinite(DM_fraction) || DM_fraction <= 0.0 || DM_fraction > 1.0)
	{
		std::cerr << "Error in Configuration::Construct_DM_Particle(): DM mass, spin, and fraction must be finite and physical." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// 3.2 DM interactions
	std::string DM_interaction;
	try
	{
		DM_interaction = config.lookup("DM_interaction").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_interaction' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// 3.2.1 SI and SD
	if(DM_interaction == "SI" || DM_interaction == "SD")
		Configuration::Construct_DM_Particle_Standard(DM_interaction);
	else if(DM_interaction == "Dark photon")
		Configuration::Construct_DM_Particle_Dark_Photon();
	else
	{
		std::cerr << "Error in DaMaSCUS_SUN::Configuration::Construct_DM_Particle(): 'DM_interaction' setting " << DM_interaction << " in configuration file not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	DM->Set_Mass(DM_mass);
	DM->Set_Spin(DM_spin);
	DM->Set_Fractional_Density(DM_fraction);
	DM->Set_Low_Mass_Mode(DM_light);
	const double primary_cross_section = DM->Get_Interaction_Parameter("Nuclei");
	const double electron_cross_section = DM->Get_Interaction_Parameter("Electrons");
	if(!std::isfinite(primary_cross_section) || primary_cross_section < 0.0
	   || !std::isfinite(electron_cross_section) || electron_cross_section < 0.0)
	{
		std::cerr << "Error in Configuration::Construct_DM_Particle(): DM cross sections must be finite and non-negative." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void Configuration::Construct_DM_Particle_Dark_Photon()
{
	DM = new DM_Particle_Dark_Photon();

	// DM form factor
	std::string DM_form_factor;
	double DM_mediator_mass = -1.0;
	try
	{
		DM_form_factor = config.lookup("DM_form_factor").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "Error in DaMaSCUS_SUN::Configuration::Construct_DM_Particle_DP(): No 'DM_form_factor' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(DM_form_factor == "General")
	{
		try
		{
			DM_mediator_mass = config.lookup("DM_mediator_mass");
			DM_mediator_mass *= MeV;
		}
		catch(const SettingNotFoundException& nfex)
		{
			std::cerr << "Error in Configuration::Construct_DM_Particle_DP(): No 'DM_mediator_mass' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		if(!std::isfinite(DM_mediator_mass) || DM_mediator_mass <= 0.0)
		{
			std::cerr << "Error in Configuration::Construct_DM_Particle_DP(): 'DM_mediator_mass' must be finite and positive for the General form factor." << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
	dynamic_cast<DM_Particle_Dark_Photon*>(DM)->Set_FormFactor_DM(DM_form_factor, DM_mediator_mass);
	double DM_cross_section_electron;
	try
	{
		DM_cross_section_electron = config.lookup("DM_cross_section_electron");
		DM_cross_section_electron *= cm * cm;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "Error in Configuration::Construct_DM_Particle_DP(): No 'DM_cross_section_electron' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	DM->Set_Interaction_Parameter(DM_cross_section_electron, "Electrons");
}

void Configuration::Print_Summary(int mpi_rank)
{
	if(mpi_rank == 0)
	{
		Print_Summary_Base(mpi_rank);
		std::cout << "DaMaSCUS-SUN options" << std::endl
				  << std::endl
				  << "\tRun mode:\t\t\t" << run_mode << std::endl
				  << "\tCapture mode:\t\t\t" << (capture_mode ? "[x]" : "[ ]") << std::endl
				  << "\tSample size:\t\t\t" << sample_size << std::endl
				  << "\tFixed PRNG seed:\t\t" << (fixed_seed == 0 ? "random" : std::to_string(fixed_seed)) << std::endl
				  << "\tMax scatterings/traj:\t\t" << maximum_number_of_scatterings << std::endl
				  << "\tSc. rate interpolation:\t\t" << ((interpolation_points > 0) ? "[x] (Grid: " + std::to_string(interpolation_points) + "×" + std::to_string(interpolation_points) + ")" : "[ ]") << std::endl;
		if(run_mode == "Parameter point" && isoreflection_rings > 1)
			std::cout << "\tIsoreflection rings:\t\t" << isoreflection_rings << std::endl;
		else if(run_mode == "Parameter scan")
			std::cout
				<< "\tCross section (min) [cm^2]:\t" << libphysica::Round(In_Units(cross_section_min, cm * cm)) << std::endl
				<< "\tCross section (max) [cm^2]:\t" << libphysica::Round(In_Units(cross_section_max, cm * cm)) << std::endl
				<< "\tCross section steps:\t\t" << cross_sections << std::endl;
		std::cout << SEPARATOR << std::endl;
	}
}

double Compute_p_Value(unsigned int sample_size, obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, unsigned int rate_interpolation_points, int mpi_rank, unsigned long int max_scatterings, SnapshotConfig snapshot_config, unsigned int fixed_seed)
{
	double u_min = detector.Minimum_DM_Speed(DM);

	solar_model.Interpolate_Total_DM_Scattering_Rate(DM, rate_interpolation_points, rate_interpolation_points);
	Simulation_Data data_set(sample_size, g_max_trajectories, u_min);
	data_set.Configure(2.0 * rSun, 1, max_scatterings);
	data_set.Generate_Data(DM, solar_model, halo_model, snapshot_config, fixed_seed, false);
	data_set.Print_Summary(mpi_rank);
	Reflection_Spectrum spectrum(data_set, solar_model, halo_model, DM.mass);
	double p = detector.P_Value(DM, spectrum);
	p = Checked_Probability(p, "Compute_p_Value()");
	return (p < 1.0e-100) ? 0.0 : p;
}

// 2. 	Class to perform parameter scans in the (m_DM, sigma)-plane to search for equal-p-value contours.
Parameter_Scan::Parameter_Scan(const std::vector<double>& masses, const std::vector<double>& coupl, std::string ID, unsigned int samplesize, unsigned int interpolation_points, double CL, unsigned long int max_scatterings)
: DM_masses(masses), couplings(coupl), sample_size(samplesize), scattering_rate_interpolation_points(interpolation_points), maximum_number_of_scatterings(max_scatterings), snapshot_config(), fixed_seed(0), certainty_level(CL)
{
	if(DM_masses.empty() || couplings.empty())
		throw std::invalid_argument("Parameter_Scan requires non-empty mass and coupling grids");
	if(DM_masses.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())
	   || couplings.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
		throw std::length_error("Parameter_Scan grid dimensions exceed the supported index range");
	if(!std::isfinite(certainty_level) || certainty_level <= 0.0 || certainty_level >= 1.0)
		throw std::invalid_argument("Parameter_Scan certainty level must be finite and lie strictly between 0 and 1");
	for(const double mass : DM_masses)
		if(!std::isfinite(mass) || mass <= 0.0)
			throw std::invalid_argument("Parameter_Scan masses must be finite and positive");
	for(const double coupling : couplings)
		if(!std::isfinite(coupling) || coupling <= 0.0)
			throw std::invalid_argument("Parameter_Scan couplings must be finite and positive");

	results_path = g_top_level_dir + "results/" + ID + "/";
	std::sort(DM_masses.begin(), DM_masses.end());
	std::sort(couplings.begin(), couplings.end());
	p_value_grid = std::vector<std::vector<double>>(couplings.size(), std::vector<double>(DM_masses.size(), -1.0));

	// Try to import previous results from an incomplete run
	Import_P_Values();
}

Parameter_Scan::Parameter_Scan(Configuration& config)
: Parameter_Scan(libphysica::Log_Space(config.constraints_mass_min, config.constraints_mass_max, config.constraints_masses), libphysica::Log_Space(config.cross_section_min, config.cross_section_max, config.cross_sections), config.ID, config.sample_size, config.interpolation_points, config.constraints_certainty, config.maximum_number_of_scatterings)
{
	snapshot_config = config.snapshot_config;
	fixed_seed = config.fixed_seed;
}

void Parameter_Scan::Import_P_Values()
{
	// Import p-values if a corresponding file exists and the grid dimensions fit.
	// CAREFUL: Changes in the grid's mininum/maximum mass/cross section will not be detected at this point.
	std::string filepath = results_path + "P_Values_Grid.txt";
	if(libphysica::File_Exists(filepath))
	{
		std::vector<std::vector<double>> imported_table = libphysica::Import_Table(filepath);
		bool dimensions_match = imported_table.size() == p_value_grid.size();
		if(dimensions_match)
			for(const auto& imported_row : imported_table)
				if(imported_row.size() != DM_masses.size())
				{
					dimensions_match = false;
					break;
				}
		if(dimensions_match)
		{
			for(const auto& imported_row : imported_table)
				for(const double p : imported_row)
					if(p != -1.0)
						Checked_Probability(p, "Parameter_Scan::Import_P_Values()");
			p_value_grid = imported_table;
		}
	}
}

bool Parameter_Scan::STA_Point_On_Grid(int row, int column)
{
	return row >= 0 && column >= 0 && row < couplings.size() && column < DM_masses.size();
}

void Parameter_Scan::STA_Go_Forward(int& row, int& column, std::string& STA_direction)
{
	if(STA_direction == "N")
		row++;
	else if(STA_direction == "E")
		column++;
	else if(STA_direction == "S")
		row--;
	else if(STA_direction == "W")
		column--;
}

void Parameter_Scan::STA_Go_Left(int& row, int& column, std::string& STA_direction)
{
	if(STA_direction == "N")
	{
		column--;
		STA_direction = "W";
	}
	else if(STA_direction == "E")
	{
		row++;
		STA_direction = "N";
	}
	else if(STA_direction == "S")
	{
		column++;
		STA_direction = "E";
	}
	else if(STA_direction == "W")
	{
		row--;
		STA_direction = "S";
	}
}

void Parameter_Scan::STA_Go_Right(int& row, int& column, std::string& STA_direction)
{
	if(STA_direction == "N")
	{
		column++;
		STA_direction = "E";
	}
	else if(STA_direction == "E")
	{
		row--;
		STA_direction = "S";
	}
	else if(STA_direction == "S")
	{
		column--;
		STA_direction = "W";
	}
	else if(STA_direction == "W")
	{
		row++;
		STA_direction = "N";
	}
}

void Parameter_Scan::STA_Fill_Gaps()
{
	for(unsigned int row = 0; row < couplings.size(); row++)
	{
		bool excluded = false;
		for(unsigned int column = 0; column < DM_masses.size(); column++)
		{
			double p = p_value_grid[row][column];
			if(p < 0)
				p_value_grid[row][column] = excluded ? 0.0 : 1.0;
			else if(Checked_Probability(p, "Parameter_Scan::STA_Fill_Gaps()") < 1.0 - certainty_level)
				excluded = true;
			else
				excluded = false;
		}
	}
}

std::vector<double> Parameter_Scan::Find_Contour_Point(int row, int column, int row_previous, int column_previous, double p_critical)
{
	if(row_previous == row)
	{
		double sigma = couplings[row];
		int column_1 = (DM_masses[column] < DM_masses[column_previous]) ? column : column_previous;
		int column_2 = (column_1 == column) ? column_previous : column;
		double x_1	 = DM_masses[column_1];
		double x_2	 = DM_masses[column_2];
		double p_1	 = p_value_grid[row][column_1];
		double p_2	 = p_value_grid[row][column_2];
		double y_1	 = Log10_Probability(p_1, "Parameter_Scan::Find_Contour_Point()");
		double y_2	 = Log10_Probability(p_2, "Parameter_Scan::Find_Contour_Point()");
		double y	 = log10(p_critical);
		const double denominator = y_2 - y_1;
		if(!std::isfinite(denominator) || denominator == 0.0)
			throw std::runtime_error("Parameter_Scan::Find_Contour_Point(): degenerate p-values cannot be interpolated");
		double x	 = (x_2 - x_1) * (y - y_1) / denominator + x_1;
		if(!std::isfinite(x))
			throw std::runtime_error("Parameter_Scan::Find_Contour_Point(): interpolation produced a non-finite mass");
		return {x, sigma};
	}
	else
	{
		double mDM = DM_masses[column];
		int row_1  = (couplings[row] < couplings[row_previous]) ? row : row_previous;
		int row_2  = (row_1 == row) ? row_previous : row;
		double x_1 = couplings[row_1];
		double x_2 = couplings[row_2];
		double p_1 = p_value_grid[row_1][column];
		double p_2 = p_value_grid[row_2][column];
		double y_1 = Log10_Probability(p_1, "Parameter_Scan::Find_Contour_Point()");
		double y_2 = Log10_Probability(p_2, "Parameter_Scan::Find_Contour_Point()");
		double y   = log10(p_critical);
		const double denominator = y_2 - y_1;
		if(!std::isfinite(denominator) || denominator == 0.0)
			throw std::runtime_error("Parameter_Scan::Find_Contour_Point(): degenerate p-values cannot be interpolated");
		double x   = (x_2 - x_1) * (y - y_1) / denominator + x_1;
		if(!std::isfinite(x))
			throw std::runtime_error("Parameter_Scan::Find_Contour_Point(): interpolation produced a non-finite coupling");
		return {mDM, x};
	}
}

std::vector<std::vector<double>> Parameter_Scan::Limit_Curve()
{
	std::vector<std::vector<double>> limit_curve;
	std::string STA_direction = "W";
	int row					  = couplings.size() - 1;
	int column				  = DM_masses.size() - 1;
	int row_previous = -10, column_previous = -10;
	double p_previous = 10.0;
	double p_critical = 1.0 - certainty_level;
	std::vector<int> first_excluded_point;
	unsigned int first_excluded_point_visits = 0;
	const std::size_t maximum_steps = STA_Step_Limit(couplings.size(), DM_masses.size());
	std::size_t steps = 0;
	while(first_excluded_point_visits < 2)
	{
		if(++steps > maximum_steps)
			throw std::runtime_error("Parameter_Scan::Limit_Curve(): square tracing exceeded its bounded step count");
		double p = STA_Point_On_Grid(row, column)
		             ? Checked_Probability(p_value_grid[row][column], "Parameter_Scan::Limit_Curve()")
		             : 1.0;
		// Abort if no point in the upper row can be excluded:
		if(first_excluded_point.empty() && p >= p_critical && row == static_cast<int>(couplings.size()) - 1 && column == 0)
			break;
		// Save the first excluded point and count how often we re-visit that point
		if(p < p_critical && first_excluded_point.empty())
			first_excluded_point = {row, column};
		if(!first_excluded_point.empty() && row == first_excluded_point[0] && column == first_excluded_point[1])
			first_excluded_point_visits++;
		// Interpolate at the boundary to find the point where p == p_critical
		if(STA_Point_On_Grid(row, column) && STA_Point_On_Grid(row_previous, column_previous)
		   && ((p < p_critical) != (p_previous < p_critical)))
		{
			std::vector<double> contour_point = Find_Contour_Point(row, column, row_previous, column_previous, p_critical);
			if(limit_curve.empty() || limit_curve.back() != contour_point)
				limit_curve.push_back(contour_point);
		}
		p_previous		= p;
		row_previous	= row;
		column_previous = column;
		if(first_excluded_point.empty())
			STA_Go_Forward(row, column, STA_direction);
		else if(p < p_critical)
			STA_Go_Left(row, column, STA_direction);
		else
			STA_Go_Right(row, column, STA_direction);
	}
	std::reverse(limit_curve.begin(), limit_curve.end());
	return limit_curve;
}

void Parameter_Scan::Perform_STA_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank)
{
	Import_P_Values();
	double mDM_original		 = DM.mass;
	double coupling_original = DM.Get_Interaction_Parameter(detector.Target_Particles());

	std::vector<int> first_excluded_point;

	int counter				  = 0;
	std::string STA_direction = "W";
	int row					  = couplings.size() - 1;
	int column				  = DM_masses.size() - 1;
	double p_critical		  = 1.0 - certainty_level;

	int first_excluded_point_counter = 0;
	const std::size_t maximum_steps = STA_Step_Limit(couplings.size(), DM_masses.size());
	std::size_t steps = 0;
	while(first_excluded_point_counter < 2)
	{
		if(++steps > maximum_steps)
			throw std::runtime_error("Parameter_Scan::Perform_STA_Scan(): square tracing exceeded its bounded step count");
		MPI_Barrier(MPI_COMM_WORLD);
		double p;
		if(!STA_Point_On_Grid(row, column))
			p = 1.0;
		else if(p_value_grid[row][column] >= 0)
			p = p_value_grid[row][column];
		else
		{
			DM.Set_Interaction_Parameter(couplings[row], detector.Target_Particles());
			DM.Set_Mass(DM_masses[column]);
			double u_min = detector.Minimum_DM_Speed(DM);
			if(mpi_rank == 0)
				std::cout << std::endl
						  << ++counter << ")\t"
						  << "m_DM [MeV]:\t" << libphysica::Round(In_Units(DM.mass, MeV)) << "\t\t"
						  << "sigma_p [cm2]:\t" << libphysica::Round(In_Units(DM.Get_Interaction_Parameter("Nuclei"), cm * cm)) << std::endl
						  << "\tu_min [km/sec]:\t" << libphysica::Round(In_Units(u_min, km / sec)) << "\t\t"
						  << "sigma_e [cm2]:\t" << libphysica::Round(In_Units(DM.Get_Interaction_Parameter("Electrons"), cm * cm)) << std::endl
						  << std::endl;
			Print_Grid(mpi_rank, row, column);
			MPI_Barrier(MPI_COMM_WORLD);

			p = Compute_p_Value(sample_size, DM, detector, solar_model, halo_model, scattering_rate_interpolation_points, mpi_rank, maximum_number_of_scatterings, snapshot_config, fixed_seed);

			p_value_grid[row][column] = p;
			libphysica::Export_Table(results_path + "P_Values_Grid.txt", p_value_grid);
			if(mpi_rank == 0)
			{
				std::cout << std::endl
						  << std::endl;
				libphysica::Print_Box("p = " + std::to_string(libphysica::Round(p)), 1);
			}
		}
		p = Checked_Probability(p, "Parameter_Scan::Perform_STA_Scan()");
		// If the upper row does not contain excluded points, we abort.
		if(first_excluded_point.empty() && p >= p_critical && row == static_cast<int>(couplings.size()) - 1 && column == 0)
			break;
		// Check if we arrived back at the first excluded point
		if(first_excluded_point.empty() && p < p_critical)
			first_excluded_point = {row, column};
		if(!first_excluded_point.empty() && row == first_excluded_point[0] && column == first_excluded_point[1])
			first_excluded_point_counter++;

		// Go to the next parameter point
		if(first_excluded_point.empty())
			STA_Go_Forward(row, column, STA_direction);
		else if(p < p_critical)
			STA_Go_Left(row, column, STA_direction);
		else
			STA_Go_Right(row, column, STA_direction);
	}
	STA_Fill_Gaps();
	Print_Grid(mpi_rank);
	libphysica::Export_Table(results_path + "P_Values_Grid.txt", p_value_grid);
	DM.Set_Mass(mDM_original);
	DM.Set_Interaction_Parameter(coupling_original, detector.Target_Particles());
}

void Parameter_Scan::Perform_Full_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank)
{
	Import_P_Values();
	double mDM_original		 = DM.mass;
	double coupling_original = DM.Get_Interaction_Parameter(detector.Target_Particles());

	double p_critical					  = 1.0 - certainty_level;
	unsigned int counter				  = 0;
	for(unsigned int i = 0; i < couplings.size(); i++)
	{
		int row			   = couplings.size() - 1 - i;
		bool row_exclusion = false;
		for(unsigned int j = 0; j < DM_masses.size(); j++)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			int column = DM_masses.size() - 1 - j;
			DM.Set_Mass(DM_masses[column]);
			DM.Set_Interaction_Parameter(couplings[row], detector.Target_Particles());
			double u_min = detector.Minimum_DM_Speed(DM);
			double p;
			if(p_value_grid[row][column] >= 0)
				p = p_value_grid[row][column];
			else
			{
				if(mpi_rank == 0)
					std::cout << std::endl
							  << ++counter << ")\t"
							  << "m_DM [MeV]:\t" << libphysica::Round(In_Units(DM.mass, MeV)) << "\t\t"
							  << "sigma_p [cm2]:\t" << libphysica::Round(In_Units(DM.Get_Interaction_Parameter("Nuclei"), cm * cm)) << std::endl
							  << "\tu_min [km/sec]:\t" << libphysica::Round(In_Units(u_min, km / sec)) << "\t\t"
							  << "sigma_e [cm2]:\t" << libphysica::Round(In_Units(DM.Get_Interaction_Parameter("Electrons"), cm * cm)) << std::endl
							  << std::endl;
				Print_Grid(mpi_rank, row, column);
				MPI_Barrier(MPI_COMM_WORLD);

				p = Compute_p_Value(sample_size, DM, detector, solar_model, halo_model, scattering_rate_interpolation_points, mpi_rank, maximum_number_of_scatterings, snapshot_config, fixed_seed);

				p_value_grid[row][column] = p;
				libphysica::Export_Table(results_path + "P_Values_Grid.txt", p_value_grid);
				if(mpi_rank == 0)
				{
					std::cout << std::endl
							  << std::endl;
					libphysica::Print_Box("p = " + std::to_string(libphysica::Round(p)), 1);
				}
			}
			p = Checked_Probability(p, "Parameter_Scan::Perform_Full_Scan()");

			if(p < p_critical)
				row_exclusion			 = true;
		}
		if(!row_exclusion)
		{
			for(int k = 0; k < row; k++)
				for(unsigned int j = 0; j < DM_masses.size(); j++)
					p_value_grid[k][j] = 1.0;
			break;
		}
	}
	libphysica::Export_Table(results_path + "P_Values_Grid.txt", p_value_grid);

	DM.Set_Mass(mDM_original);
	DM.Set_Interaction_Parameter(coupling_original, detector.Target_Particles());
}

void Parameter_Scan::Print_Grid(int mpi_rank, int marker_row, int marker_column)
{
	if(mpi_rank == 0)
	{
		double p_critical = 1.0 - certainty_level;
		std::cout << "\t┌";
		for(unsigned int column = 0; column < DM_masses.size(); column++)
			std::cout << "─";
		std::cout << "┐ σ [cm^2]" << std::endl;
		for(unsigned int row = 0; row < couplings.size(); row++)
		{
			std::cout << "\t┤";
			for(unsigned int column = 0; column < DM_masses.size(); column++)
			{
				double p = p_value_grid[couplings.size() - 1 - row][column];
				if(row == couplings.size() - 1 - marker_row && column == marker_column)
					std::cout << "¤";
				else if(p < 0.0)
					std::cout << "·";
				else if(p < p_critical)
					std::cout << "█";
				else if(p > p_critical)
					std::cout << "░";
			}
			std::cout << "├";	// ((row == 0 || row == couplings.size() - 1) ? "├" : "│");
			if(row == 0)
				std::cout << " " << libphysica::Round(In_Units(couplings.back(), cm * cm));
			else if(row == couplings.size() - 1)
				std::cout << " " << libphysica::Round(In_Units(couplings.front(), cm * cm));
			std::cout << std::endl;
		}
		std::cout << "\t└";
		for(unsigned int column = 0; column < DM_masses.size(); column++)
			std::cout << ((column == 0 || column == DM_masses.size() - 1) ? "┬" : "─");
		if(DM_masses.size() > 5)
		{
			std::cout << "┘ mDM[MeV]\n\t ";
			for(unsigned int i = 0; i < DM_masses.size() - 1; i++)
				std::cout
					<< " ";
			std::cout << libphysica::Round(In_Units(DM_masses.back(), MeV))
					  << "\r\t " << libphysica::Round(In_Units(DM_masses.front(), MeV)) << std::endl;
		}
		std::cout << std::endl;
	}
}

void Parameter_Scan::Export_Results(int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::vector<std::vector<double>> table;
		for(unsigned int i = 0; i < DM_masses.size(); i++)
			for(unsigned int j = 0; j < couplings.size(); j++)
				table.push_back({DM_masses[i], couplings[j], p_value_grid[j][i]});
		libphysica::Export_Table(results_path + "P_Values_List.txt", table, {GeV, cm * cm, 1.0});
		int CL										   = std::round(100.0 * certainty_level);
		std::vector<std::vector<double>> limit_contour = Limit_Curve();
		libphysica::Export_Table(results_path + "Reflection_Limit_" + std::to_string(CL) + ".txt", limit_contour, {GeV, cm * cm});
	}
}

}	// namespace DaMaSCUS_SUN

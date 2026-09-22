#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>	 // for strlen
#include <exception>
#include <iostream>
#include <memory>
#include <mpi.h>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "Data_Generation.hpp"
#include "Parameter_Scan.hpp"
#include "Solar_Model.hpp"
#include "version.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

int main(int argc, char* argv[])
{
	int mpi_thread_provided = MPI_THREAD_SINGLE;
	if(MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_thread_provided) != MPI_SUCCESS)
	{
		std::cerr << "Error: MPI initialization failed." << std::endl;
		return 1;
	}
	int mpi_processes, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	if(argc < 2 || argc > 4)
	{
		if(mpi_rank == 0)
			std::cerr << "Usage: " << argv[0] << " <config.cfg> [--diagnostic | --thermal-validation]" << std::endl;
		MPI_Finalize();
		return 1;
	}

	bool diagnostic = false, thermal_validation = false;
	for(int argument = 2; argument < argc; ++argument)
	{
		const std::string option = argv[argument];
		if(option == "--diagnostic") diagnostic = true;
		else if(option == "--thermal-validation") thermal_validation = true;
		else
		{
			if(mpi_rank == 0) std::cerr << "Unknown option: " << option << std::endl;
			MPI_Finalize(); return 1;
		}
	}
	// Keep library/configuration progress away from the machine-readable Capture record.
	struct LogStreamScope
	{
		std::streambuf* original = std::cout.rdbuf(std::cerr.rdbuf());
		~LogStreamScope() { std::cout.rdbuf(original); }
	} log_stream;

	// Initial terminal output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	if(mpi_rank == 0)
		std::cout << "[Started on " << ctime_start << "]" << std::endl
				  << PROJECT_NAME << "-" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
				  << DAMASCUS_SUN_LOGO
				  << std::endl
				  << "MPI processes:\t" << mpi_processes << std::endl;

	// Configuration parameters
	std::unique_ptr<Configuration> configuration;
	try
	{
		configuration.reset(new Configuration(argv[1], mpi_rank, diagnostic, thermal_validation));
	}
	catch(const std::exception& error)
	{
		if(mpi_rank == 0) std::cerr << "Error: " << error.what() << std::endl;
		MPI_Finalize(); return 1;
	}
	Configuration& cfg = *configuration;
	if(!cfg.capture_mode) std::cout.rdbuf(log_stream.original);
	if(cfg.diagnostic_mode) g_top_level_dir += "diagnostics/";
	if(cfg.snapshot_config.enabled && mpi_thread_provided < MPI_THREAD_FUNNELED)
	{
		if(mpi_rank == 0)
			std::cerr << "Warning: MPI implementation does not provide MPI_THREAD_FUNNELED; "
			          << "heartbeat snapshot is disabled. Final MPI-reduced outputs are unaffected." << std::endl;
		cfg.snapshot_config.enabled = false;
	}
	std::string solar_model_data_file;
	try
	{
		solar_model_data_file = Locate_Solar_Model_Data_File(argv[0]);
	}
	catch(const std::exception& error)
	{
		if(mpi_rank == 0)
			std::cerr << "Error: " << error.what() << std::endl;
		MPI_Finalize();
		return 1;
	}
	Solar_Model SSM(solar_model_data_file);
	cfg.Print_Summary(mpi_rank);
	MPI_Barrier(MPI_COMM_WORLD);
	////////////////////////////////////////////////////////////////////////

	// Generate data for one parameter point specified in the configuration file.
	if(cfg.run_mode == "Parameter point" || cfg.run_mode == "Capture")
	{
		double u_min = 0.0;
		Simulation_Data data_set(cfg.sample_size, cfg.max_trajectories, u_min, cfg.isoreflection_rings);
		data_set.physical_config_json = cfg.Physical_Configuration_JSON();
		data_set.physical_config_header = cfg.Physical_Configuration_Header();
		data_set.diagnostic_output_enabled = cfg.diagnostic_mode;
		data_set.Configure(TRAJECTORY_BOUNDARY_RSUN * rSun, 1, cfg.maximum_number_of_scatterings);
		data_set.Configure_Trajectory_Diagnostics(cfg.trajectory_diagnostic_config);
		const std::string output_prefix = "results_";
		const std::string output_path = g_top_level_dir + output_prefix + std::to_string(log10(In_Units(cfg.DM->mass, GeV))) + "_" + std::to_string(log10(In_Units(cfg.DM->Sigma_Proton(), cm * cm))) + "/";
		// All ranks take the same failure path, including errors that occur only
		// on rank zero while opening, flushing, or publishing the final files.
		auto root_output_succeeded = [&](auto action) {
			int success = 1;
			if(mpi_rank == 0)
			{
				try
				{
					action();
				}
				catch(const std::exception& error)
				{
					std::cerr << "Error: " << error.what() << std::endl;
					success = 0;
				}
			}
			MPI_Bcast(&success, 1, MPI_INT, 0, MPI_COMM_WORLD);
			return success != 0;
		};
		if(!cfg.capture_mode && !root_output_succeeded([&]() {
			data_set.Prepare_Output_Directory(output_path);
		}))
		{
			MPI_Finalize();
			return 1;
		}
		if(mpi_rank == 0)
			std::cout << (cfg.capture_mode ? "Generate data in CAPTURE MODE..." : "Generate data...") << std::endl
					  << "\tm_DM [MeV]:\t" << libphysica::Round(In_Units(cfg.DM->mass, MeV)) << "\t\t"
					  << "sigma_p [cm2]:\t" << libphysica::Round(In_Units(cfg.DM->Get_Interaction_Parameter("Nuclei"), cm * cm)) << std::endl
					  << "\tu_min [km/sec]:\t" << libphysica::Round(In_Units(u_min, km / sec)) << "\t\t"
					  << "sigma_e [cm2]:\t" << libphysica::Round(In_Units(cfg.DM->Get_Interaction_Parameter("Electrons"), cm * cm)) << std::endl
					  << std::endl;
		SSM.Interpolate_Total_DM_Scattering_Rate(
		    *cfg.DM, cfg.rate_radius_points, cfg.rate_speed_points, cfg.rate_max_speed);

		data_set.outer_removal_radius_rsun = cfg.outer_removal_radius_rsun;
		data_set.interpolation_points = cfg.interpolation_points;
		data_set.thermal_shape_run = cfg.thermal_validation_mode;
		data_set.Generate_Data(*cfg.DM, SSM, *cfg.DM_distr, cfg.snapshot_config, cfg.fixed_seed, cfg.capture_mode);
		if(cfg.capture_mode)
		{
			data_set.Print_Capture_Mode_Summary(mpi_rank);
			if(!root_output_succeeded([&]() {
				std::cout.rdbuf(log_stream.original);
				data_set.Print_Capture_Result_JSON(*cfg.DM, *cfg.DM_distr);
			})) { MPI_Finalize(); return 1; }
			const bool target_reached = data_set.Target_Reached();
			MPI_Finalize();
			return target_reached ? 0 : 2;
		}
		data_set.Print_Summary(mpi_rank);

		if(cfg.diagnostic_mode)
		{
			if(!root_output_succeeded([&]() { data_set.Write_Diagnostic_Output(output_path, *cfg.DM); }))
			{ MPI_Finalize(); return 1; }
		}
		else
		{
			if(!root_output_succeeded([&]() { data_set.Write_Bincount(output_path, *cfg.DM, *cfg.DM_distr); }))
			{ MPI_Finalize(); return 1; }
			if(!data_set.Target_Reached())
			{
				if(mpi_rank == 0) std::cerr << "Transport target not reached within the configured trajectory budget. See bincount.tsv." << std::endl;
				MPI_Finalize(); return 2;
			}
		}

	}
	// Perform a parameter scan to compute exclusion limits
	else if(cfg.run_mode == "Parameter scan")
	{
		if(mpi_rank == 0 && cfg.compute_halo_constraints)
		{
			std::cout << "Compute halo constraints for " << cfg.DM_detector->name << ":" << std::endl;
			double mDM_min								= cfg.DM_detector->Minimum_DM_Mass(*cfg.DM, *cfg.DM_distr);
			std::vector<double> DM_masses				= libphysica::Log_Space(mDM_min, GeV, 100);
			std::vector<std::vector<double>> halo_limit = cfg.DM_detector->Upper_Limit_Curve(*cfg.DM, *cfg.DM_distr, DM_masses, cfg.constraints_certainty);
			int CL										= std::round(100.0 * cfg.constraints_certainty);
			libphysica::Export_Table(g_top_level_dir + "results/" + cfg.ID + "/Halo_Limit_" + std::to_string(CL) + ".txt", halo_limit, {GeV, cm * cm});
		}
		Parameter_Scan scan(cfg);
		if(cfg.perform_full_scan)
			scan.Perform_Full_Scan(*cfg.DM, *cfg.DM_detector, SSM, *cfg.DM_distr, mpi_rank);
		else
			scan.Perform_STA_Scan(*cfg.DM, *cfg.DM_detector, SSM, *cfg.DM_distr, mpi_rank);
		scan.Export_Results(mpi_rank);
		if(mpi_rank == 0)
		{
			int CL = std::round(100.0 * cfg.constraints_certainty);
			std::cout << "\nFinal reflection constraints (" << CL << "% CL)" << std::endl;
			scan.Print_Grid(mpi_rank);
		}
	}
	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	MPI_Barrier(MPI_COMM_WORLD);
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	if(mpi_rank == 0)
		std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	MPI_Finalize();
	return 0;
}

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

#include "libphysica/Natural_Units.hpp"

#include "Parameter_Scan.hpp"
#include "Simulation_Trajectory.hpp"
#include "Solar_Model.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

namespace
{
std::vector<std::string> Split_Tabs(const std::string& line)
{
	std::vector<std::string> fields;
	std::string field;
	std::istringstream stream(line);
	while(std::getline(stream, field, '\t'))
		fields.push_back(field);
	return fields;
}

double Parse_Double(const std::string& value, const std::string& field)
{
	char* end = nullptr;
	const double parsed = std::strtod(value.c_str(), &end);
	if(end == value.c_str() || *end != '\0')
		throw std::runtime_error(
		    "invalid floating-point value for " + field + ": " + value);
	return parsed;
}

uint64_t Parse_UInt64(const std::string& value, const std::string& field)
{
	char* end = nullptr;
	const unsigned long long parsed = std::strtoull(value.c_str(), &end, 10);
	if(end == value.c_str() || *end != '\0')
		throw std::runtime_error(
		    "invalid integer value for " + field + ": " + value);
	return static_cast<uint64_t>(parsed);
}

std::string Decode_RNG_State(std::string state)
{
	std::replace(state.begin(), state.end(), ',', ' ');
	return state;
}

std::map<std::string, std::string> Read_Ledger_Row(
	const std::string& path,
	int requested_rank,
	uint64_t requested_trajectory_id)
{
	std::ifstream input(path);
	if(!input)
		throw std::runtime_error("cannot open invalid trajectory ledger: " + path);

	std::vector<std::string> columns;
	std::string line;
	while(std::getline(input, line))
	{
		if(line.empty() || line[0] == '#')
			continue;
		const std::vector<std::string> fields = Split_Tabs(line);
		if(columns.empty())
		{
			columns = fields;
			continue;
		}
		if(fields.size() != columns.size())
			throw std::runtime_error("invalid ledger row column count");
		std::map<std::string, std::string> row;
		for(size_t index = 0; index < columns.size(); index++)
			row[columns[index]] = fields[index];
		if(static_cast<int>(Parse_UInt64(row.at("rank"), "rank"))
		       == requested_rank
		   && Parse_UInt64(row.at("trajectory_id"), "trajectory_id")
		          == requested_trajectory_id)
			return row;
	}
	throw std::runtime_error(
	    "requested rank/trajectory_id was not found in the ledger");
}

Event Event_From_Ledger(const std::map<std::string, std::string>& row)
{
	libphysica::Vector position(3);
	libphysica::Vector velocity(3);
	const char* position_fields[3] = {
	    "initial_x_km", "initial_y_km", "initial_z_km"};
	const char* velocity_fields[3] = {
	    "initial_vx_km_s", "initial_vy_km_s", "initial_vz_km_s"};
	for(size_t component = 0; component < 3; component++)
	{
		position[component] =
		    Parse_Double(
		        row.at(position_fields[component]),
		        position_fields[component])
		    * km;
		velocity[component] =
		    Parse_Double(
		        row.at(velocity_fields[component]),
		        velocity_fields[component])
		    * km / sec;
	}
	return Event(
	    Parse_Double(row.at("initial_time_s"), "initial_time_s") * sec,
	    position,
	    velocity);
}
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int mpi_rank = 0;
	int mpi_processes = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	int result_code = 0;
	try
	{
		if(argc != 5)
			throw std::runtime_error(
			    std::string("usage: ") + argv[0]
			    + " <config.cfg> <invalid_trajectories.tsv> <rank> <trajectory_id>");
		if(mpi_processes != 1)
			throw std::runtime_error(
			    "invalid-trajectory replay requires exactly one MPI process");

		const int requested_rank =
		    static_cast<int>(Parse_UInt64(argv[3], "rank"));
		const uint64_t requested_trajectory_id =
		    Parse_UInt64(argv[4], "trajectory_id");
		const auto row =
		    Read_Ledger_Row(
		        argv[2], requested_rank, requested_trajectory_id);

		Configuration cfg(argv[1], 0);
		const std::string solar_model_data_file =
		    Locate_Solar_Model_Data_File(argv[0]);
		Solar_Model solar_model(solar_model_data_file);
		solar_model.Interpolate_Total_DM_Scattering_Rate(
		    *cfg.DM,
		    cfg.interpolation_points,
		    cfg.interpolation_points);

		Trajectory_Simulator simulator(
		    solar_model,
		    DEFAULT_MAXIMUM_FREE_TIME_STEPS,
		    cfg.maximum_number_of_scatterings,
		    TRAJECTORY_BOUNDARY_RSUN * g_body_radius);
		simulator.current_trajectory_id =
		    requested_trajectory_id > 0
		    ? requested_trajectory_id - 1
		    : 0;
		simulator.Restore_PRNG_State(
		    Decode_RNG_State(
		        row.at("rng_state_before_simulation")));
		simulator.Enable_Diagnostic_Trace(true);
		const Event initial_event = Event_From_Ledger(row);
		const Trajectory_Result replay =
		    simulator.Simulate(initial_event, *cfg.DM, requested_rank);

		const double final_radius = replay.final_event.Radius();
		const double final_radial_velocity =
		    final_radius > 0.0
		    ? replay.final_event.position.Dot(replay.final_event.velocity)
		          / final_radius
		    : 0.0;
		const double final_escape_speed =
		    solar_model.Local_Escape_Speed(final_radius);
		const double final_energy =
		    0.5 * cfg.DM->mass
		    * (replay.final_event.Speed() * replay.final_event.Speed()
		       - final_escape_speed * final_escape_speed);

		std::cout
		    << std::scientific
		    << std::setprecision(17)
		    << "replay_rank=" << requested_rank << "\n"
		    << "replay_trajectory_id=" << requested_trajectory_id << "\n"
		    << "original_termination_reason="
		    << row.at("termination_reason") << "\n"
		    << "original_failure_detail="
		    << row.at("numerical_failure_detail") << "\n"
		    << "replay_termination_reason="
		    << static_cast<int>(replay.bincount.termination_reason)
		    << "\n"
		    << "replay_failure_detail="
		    << TrajectoryNumericalFailureDetailKey(
		           replay.bincount.numerical_failure_detail)
		    << "\n"
		    << "replay_is_captured="
		    << (replay.bincount.is_captured ? 1 : 0) << "\n"
		    << "replay_survival_valid="
		    << (replay.bincount.survival_valid ? 1 : 0) << "\n"
		    << "replay_event_observed="
		    << (replay.bincount.event_observed ? 1 : 0) << "\n"
		    << "replay_scatterings="
		    << replay.number_of_scatterings << "\n"
		    << "replay_final_r_Rbody="
		    << In_Units(final_radius, g_body_radius) << "\n"
		    << "replay_final_vr_km_s="
		    << In_Units(final_radial_velocity, km / sec) << "\n"
		    << "replay_final_energy_eV="
		    << In_Units(final_energy, eV) << "\n"
		    << "replay_diagnostic_events="
		    << replay.diagnostic_events.size() << "\n";
	}
	catch(const std::exception& error)
	{
		if(mpi_rank == 0)
			std::cerr << "Replay error: " << error.what() << std::endl;
		result_code = 1;
	}
	MPI_Finalize();
	return result_code;
}

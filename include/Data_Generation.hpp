#ifndef __Data_Generation_hpp_
#define __Data_Generation_hpp_

#include <vector>
#include <string>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/DM_Particle.hpp"

#include "Simulation_Trajectory.hpp"

namespace DaMaSCUS_SUN
{

// Evaporation time record for a single captured trajectory
struct EvaporationRecord
{
	unsigned long int trajectory_id;
	double t_evap;    // seconds
	bool truncated;   // true if last step had E <= 0
};

class Simulation_Data
{
  private:
	// Configuration
	unsigned int target_captured_per_rank;   // ceil(sample_size / N_ranks)
	unsigned long int max_trajectories_per_rank;
	double initial_and_final_radius = 2.0 * libphysica::natural_units::rSun;
	unsigned int minimum_number_of_scatterings = 1;
	long int maximum_number_of_scatterings = 1e11;
	unsigned long int maximum_free_time_steps = 1e12;

	// Results
	unsigned long int number_of_trajectories;
	unsigned long int number_of_free_particles;
	unsigned long int number_of_reflected_particles;
	unsigned long int number_of_captured_particles;
	double average_number_of_scatterings;
	double computing_time;
	bool early_stopped;

	// Aggregated bincount histograms
	std::array<double, NUM_BINS> captured_dt_hist;
	std::array<double, NUM_BINS> captured_v2dt_hist;
	std::array<double, NUM_BINS> not_captured_dt_hist;
	std::array<double, NUM_BINS> not_captured_v2dt_hist;

	// Per-bin sum of squares for error estimation
	std::array<double, NUM_BINS> captured_dt_sq_hist;      // Σ (per-traj dt)²
	std::array<double, NUM_BINS> captured_v2dt_sq_hist;    // Σ (per-traj v²dt)²
	std::array<double, NUM_BINS> not_captured_dt_sq_hist;
	std::array<double, NUM_BINS> not_captured_v2dt_sq_hist;

	// Evaporation records
	std::vector<EvaporationRecord> evaporation_records;

	// MPI
	int mpi_rank, mpi_processes;
	void Perform_MPI_Reductions();

	// For reflection spectrum compatibility (kept but not actively used in new logic)
	unsigned int isoreflection_rings;
	double minimum_speed_threshold;
	std::vector<unsigned long int> number_of_data_points;
	double KDE_boundary_correction_factor = 0.75;

  public:
	std::vector<std::vector<libphysica::DataPoint>> data;

	Simulation_Data(unsigned int sample_size, unsigned int max_trajectories, double u_min = 0.0, unsigned int iso_rings = 1);

	void Configure(double initial_radius, unsigned int min_scattering, long int max_scattering, unsigned long int max_free_steps = 1e12);

	void Generate_Data(obscura::DM_Particle& DM, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, SnapshotConfig snapshot_cfg = SnapshotConfig(), unsigned int fixed_seed = 0);

	// Output files
	void Write_Output_Files(const std::string& output_dir, obscura::DM_Particle& DM);

	double Free_Ratio() const;
	double Capture_Ratio() const;
	double Reflection_Ratio(int isoreflection_ring = -1) const;

	double Minimum_Speed() const;
	double Lowest_Speed(unsigned int iso_ring = 0) const;
	double Highest_Speed(unsigned int iso_ring = 0) const;

	void Print_Summary(unsigned int mpi_rank = 0);
};
}	// namespace DaMaSCUS_SUN
#endif

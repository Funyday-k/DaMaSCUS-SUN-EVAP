#ifndef __Radial_Binning_hpp_
#define __Radial_Binning_hpp_

#include <cstddef>
#include <vector>

namespace DaMaSCUS_SUN
{

constexpr double BINCOUNT_RADIAL_GRID_AU_KM = 1.495978707e8;

// Runtime radial grid for the residence-time bincount.
//
// Default construction rules:
//   - uniform inner grid: dr = 0.001 * R_body, through 1.1 * R_body;
//   - exterior grid: starts with the inner width, grows geometrically,
//     and is capped at 10 * R_body;
//   - final shell is shortened so the domain ends exactly at domain_max_au.
class Bincount_Radial_Grid
{
  public:
    explicit Bincount_Radial_Grid(
        double body_radius_km,
        double domain_max_au = 1.0,
        double inner_bin_width_fraction = 0.001,
        double inner_extent_fraction = 1.1,
        double exterior_growth_factor = 1.02,
        double exterior_max_bin_width_fraction = 10.0);

    double Body_Radius_Km() const;
    double Inner_Bin_Width_Km() const;
    double Inner_Extent_Km() const;
    double Domain_Max_Km() const;
    double Exterior_Growth_Factor() const;
    double Exterior_Max_Bin_Width_Km() const;

    std::size_t Inner_Bin_Count() const;
    std::size_t Exterior_Bin_Count() const;
    std::size_t Bin_Count() const;

    const std::vector<double>& Bin_Edges_Km() const;

    double Bin_Lower_Km(std::size_t bin) const;
    double Bin_Upper_Km(std::size_t bin) const;

    // Returns -1 outside [0, Domain_Max_Km()).
    int Bin_Index_Km(double radius_km) const;

  private:
    double body_radius_km_;
    double inner_bin_width_km_;
    double inner_extent_km_;
    double domain_max_km_;
    double exterior_growth_factor_;
    double exterior_max_bin_width_km_;
    std::size_t inner_bin_count_;
    std::vector<double> bin_edges_km_;
};

} // namespace DaMaSCUS_SUN

#endif

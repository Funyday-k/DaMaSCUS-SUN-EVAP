#include "Radial_Binning.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace DaMaSCUS_SUN
{

Bincount_Radial_Grid::Bincount_Radial_Grid(
    double body_radius_km,
    double domain_max_au,
    double inner_bin_width_fraction,
    double inner_extent_fraction,
    double exterior_growth_factor,
    double exterior_max_bin_width_fraction)
: body_radius_km_(body_radius_km),
  inner_bin_width_km_(0.0),
  inner_extent_km_(0.0),
  domain_max_km_(0.0),
  exterior_growth_factor_(exterior_growth_factor),
  exterior_max_bin_width_km_(0.0),
  inner_bin_count_(0)
{
    if(!std::isfinite(body_radius_km_) || body_radius_km_ <= 0.0)
        throw std::invalid_argument(
            "Bincount_Radial_Grid: body radius must be finite and positive.");

    if(!std::isfinite(domain_max_au) || domain_max_au <= 0.0)
        throw std::invalid_argument(
            "Bincount_Radial_Grid: domain maximum must be finite and positive.");

    if(!std::isfinite(inner_bin_width_fraction)
       || inner_bin_width_fraction <= 0.0
       || !std::isfinite(inner_extent_fraction)
       || inner_extent_fraction <= 0.0)
        throw std::invalid_argument(
            "Bincount_Radial_Grid: inner-grid parameters must be finite and positive.");

    if(!std::isfinite(exterior_growth_factor_)
       || exterior_growth_factor_ < 1.0
       || !std::isfinite(exterior_max_bin_width_fraction)
       || exterior_max_bin_width_fraction <= 0.0)
        throw std::invalid_argument(
            "Bincount_Radial_Grid: exterior-grid parameters are invalid.");

    inner_bin_width_km_ =
        inner_bin_width_fraction * body_radius_km_;

    inner_extent_km_ =
        inner_extent_fraction * body_radius_km_;

    domain_max_km_ =
        domain_max_au * BINCOUNT_RADIAL_GRID_AU_KM;

    exterior_max_bin_width_km_ =
        exterior_max_bin_width_fraction * body_radius_km_;

    if(!(domain_max_km_ > inner_extent_km_))
        throw std::invalid_argument(
            "Bincount_Radial_Grid: outer domain must exceed the inner-grid extent.");

    const double raw_inner_bin_count =
        inner_extent_fraction / inner_bin_width_fraction;

    const double rounded_inner_bin_count =
        std::round(raw_inner_bin_count);

    if(std::fabs(raw_inner_bin_count - rounded_inner_bin_count)
       > 1.0e-12 * std::max(1.0, std::fabs(raw_inner_bin_count)))
        throw std::invalid_argument(
            "Bincount_Radial_Grid: inner extent must be an integer number of bins.");

    inner_bin_count_ =
        static_cast<std::size_t>(rounded_inner_bin_count);

    if(inner_bin_count_ == 0)
        throw std::invalid_argument(
            "Bincount_Radial_Grid: inner grid must contain at least one bin.");

    bin_edges_km_.reserve(inner_bin_count_ + 512);

    for(std::size_t edge = 0; edge <= inner_bin_count_; edge++)
    {
        if(edge == inner_bin_count_)
            bin_edges_km_.push_back(inner_extent_km_);
        else
            bin_edges_km_.push_back(
                static_cast<double>(edge) * inner_bin_width_km_);
    }

    double edge_km  = inner_extent_km_;
    double width_km = inner_bin_width_km_;

    while(edge_km < domain_max_km_)
    {
        width_km = std::min(
            width_km,
            exterior_max_bin_width_km_);

        const double next_edge_km = std::min(
            domain_max_km_,
            edge_km + width_km);

        if(!(next_edge_km > edge_km))
            throw std::runtime_error(
                "Bincount_Radial_Grid: exterior grid failed to advance.");

        bin_edges_km_.push_back(next_edge_km);
        edge_km = next_edge_km;

        width_km = std::min(
            exterior_max_bin_width_km_,
            width_km * exterior_growth_factor_);
    }

    // Keep the outer edge exact despite accumulated floating-point rounding.
    bin_edges_km_.back() = domain_max_km_;
}

double Bincount_Radial_Grid::Body_Radius_Km() const
{
    return body_radius_km_;
}

double Bincount_Radial_Grid::Inner_Bin_Width_Km() const
{
    return inner_bin_width_km_;
}

double Bincount_Radial_Grid::Inner_Extent_Km() const
{
    return inner_extent_km_;
}

double Bincount_Radial_Grid::Domain_Max_Km() const
{
    return domain_max_km_;
}

double Bincount_Radial_Grid::Exterior_Growth_Factor() const
{
    return exterior_growth_factor_;
}

double Bincount_Radial_Grid::Exterior_Max_Bin_Width_Km() const
{
    return exterior_max_bin_width_km_;
}

std::size_t Bincount_Radial_Grid::Inner_Bin_Count() const
{
    return inner_bin_count_;
}

std::size_t Bincount_Radial_Grid::Exterior_Bin_Count() const
{
    return Bin_Count() - Inner_Bin_Count();
}

std::size_t Bincount_Radial_Grid::Bin_Count() const
{
    return bin_edges_km_.size() - 1;
}

const std::vector<double>& Bincount_Radial_Grid::Bin_Edges_Km() const
{
    return bin_edges_km_;
}

double Bincount_Radial_Grid::Bin_Lower_Km(std::size_t bin) const
{
    if(bin > Bin_Count())
        return std::numeric_limits<double>::quiet_NaN();

    return bin_edges_km_[bin];
}

double Bincount_Radial_Grid::Bin_Upper_Km(std::size_t bin) const
{
    if(bin >= Bin_Count())
        return std::numeric_limits<double>::quiet_NaN();

    return bin_edges_km_[bin + 1];
}

int Bincount_Radial_Grid::Bin_Index_Km(double radius_km) const
{
    if(!std::isfinite(radius_km)
       || radius_km < 0.0
       || radius_km >= domain_max_km_)
        return -1;

    const std::vector<double>::const_iterator upper =
        std::upper_bound(
            bin_edges_km_.begin(),
            bin_edges_km_.end(),
            radius_km);

    return static_cast<int>(
        std::distance(bin_edges_km_.begin(), upper) - 1);
}

} // namespace DaMaSCUS_SUN

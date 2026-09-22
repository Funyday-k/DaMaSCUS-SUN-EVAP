#ifndef DAMASCUS_OUTPUT_STATISTICS_HPP
#define DAMASCUS_OUTPUT_STATISTICS_HPP

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace DaMaSCUS_SUN
{
constexpr std::size_t RESIDENCE_JACKKNIFE_BLOCKS = 64;
constexpr int OUTPUT_STATISTIC_DIGITS = 4;
constexpr int OUTPUT_PARAMETER_DIGITS = 12;
using OutputBlockCounts = std::array<unsigned long, RESIDENCE_JACKKNIFE_BLOCKS>;
using OutputMoments = std::array<long double, 7>;
using OutputReplicates = std::array<long double, RESIDENCE_JACKKNIFE_BLOCKS>;

// Decimal input parameters / geometry omit insignificant trailing zeros.
// Statistical values and errors are rounded only when serialized.
std::string Format_Output_Number(double value, int significant_digits = OUTPUT_STATISTIC_DIGITS);

// sigma [cm^2] = coefficient * 10^(-exponent); exponent is always an integer.
void Write_Cross_Section_Header(std::ostream& output, const std::string& species, double sigma_cm2);

// Fixed-N error of a sum from per-history first and second moments.
// Independent of how histories were assigned to jackknife blocks.
double Moment_Sum_SE(long double sum, long double square_sum, unsigned long count);

// k_B T [eV], then residence/ever/never pair times and ever-never cross time [s^2].
// temperature_factor = m_chi[eV] / (3 c[km/s]^2).
std::array<long double, 5> Derived_Radial_Observables(
    const OutputMoments& moments, unsigned long captured, unsigned long never,
    long double temperature_factor);

// All 64 hash blocks participate. A missing replicate makes this covariance NaN.
long double Block_Jackknife_Covariance(const OutputReplicates& a, const OutputReplicates& b);

// Integral_0^r r'^2 <|D-r'|^-2>_Omega dr', in cm. Requires 0 <= r < D.
long double Detector_Radial_Integral(long double radius_cm, long double distance_cm);

// Append compact derived statistics, selected integral covariances and radial rows.
// Block histograms are bin-major; moments are res S, res Q, res v2S, ever S/Q, never S/Q.
void Write_Compact_Radial_Statistics(
    std::ostream& output, const std::vector<double>& edges_km,
    const std::array<const std::vector<double>*, 7>& block_histograms,
    const std::vector<double>& residence_v2_square_sum,
    const OutputBlockCounts& captured_counts, const OutputBlockCounts& never_counts,
    long double temperature_factor);
}
#endif

#include "Output_Statistics.hpp"
#include "Simulation_Trajectory.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>

namespace DaMaSCUS_SUN
{
namespace
{
constexpr long double missing = std::numeric_limits<long double>::quiet_NaN();

long double Nonnegative_Difference(long double difference, long double scale)
{
    if(difference < -128 * std::numeric_limits<double>::epsilon() * scale)
        throw std::runtime_error("inconsistent first and second output moments");
    return std::max(0.0L, difference);
}

long double Pair_Time(long double sum, long double squares, unsigned long count)
{
    if(count < 2) return missing;
    return Nonnegative_Difference(sum * sum - squares, std::max(sum * sum, squares))
           / (static_cast<long double>(count) * (count - 1));
}

struct RadialStatistics
{
    OutputMoments totals{};
    std::array<long double, 5> central{};
    std::array<OutputReplicates, 5> deleted{};
};

struct Window
{
    const char* name;
    long double low_rsun, high_rsun;
    bool flux;
};

// Density is constant within each output shell. Clip boundary shells by volume,
// or by the exact spherical detector kernel for the gamma event flux.
long double Window_Weight(const Window& window, long double low_km, long double high_km)
{
    const long double low = std::max(low_km, window.low_rsun * R_SUN_KM) * 1e5L;
    const long double high = std::min(high_km, window.high_rsun * R_SUN_KM) * 1e5L;
    if(high <= low) return 0.0L;
    const long double a = low_km * 1e5L, b = high_km * 1e5L;
    const long double four_pi_over_three = 4 * std::acos(-1.0L) / 3;
    const long double volume = four_pi_over_three * (b - a) * (b*b + a*b + a*a);
    const long double integral = window.flux
        ? Detector_Radial_Integral(high, AU_KM * 1e5L) - Detector_Radial_Integral(low, AU_KM * 1e5L)
        : four_pi_over_three * (high - low) * (high*high + low*high + low*low);
    return integral / (volume * volume);
}

void Write_Number(std::ostream& output, long double number)
{
    const double stored = static_cast<double>(number);
    if(std::isinf(stored)) throw std::runtime_error("nonfinite derived output statistic");
    output << stored;
}
}

double Moment_Sum_SE(long double sum, long double square_sum, unsigned long count)
{
    if(count < 2) return static_cast<double>(missing);
    const long double n_q = count * square_sum;
    return static_cast<double>(std::sqrt(Nonnegative_Difference(
        n_q - sum * sum, std::max(n_q, sum * sum)) / (count - 1)));
}

std::array<long double, 5> Derived_Radial_Observables(
    const OutputMoments& m, unsigned long captured, unsigned long never, long double temperature_factor)
{
    return {{m[0] > 0 ? temperature_factor * m[2] / m[0] : missing,
             Pair_Time(m[0], m[1], captured), Pair_Time(m[3], m[4], captured),
             Pair_Time(m[5], m[6], never),
             captured > 0 && never > 0 ? (m[3] / captured) * (m[5] / never) : missing}};
}

long double Block_Jackknife_Covariance(const OutputReplicates& a, const OutputReplicates& b)
{
    for(std::size_t i = 0; i < a.size(); ++i)
        if(!std::isfinite(a[i]) || !std::isfinite(b[i])) return missing;
    const long double mean_a = std::accumulate(a.begin(), a.end(), 0.0L) / a.size();
    const long double mean_b = std::accumulate(b.begin(), b.end(), 0.0L) / b.size();
    long double covariance = 0;
    for(std::size_t i = 0; i < a.size(); ++i)
        covariance += (a[i] - mean_a) * (b[i] - mean_b);
    return covariance * (a.size() - 1) / a.size();
}

long double Detector_Radial_Integral(long double radius_cm, long double distance_cm)
{
    if(!(radius_cm >= 0 && distance_cm > radius_cm))
        throw std::invalid_argument("detector kernel requires 0 <= radius < distance");
    const long double x = radius_cm / distance_cm;
    if(x < 1e-3L)
    {
        const long double x2 = x * x;
        return distance_cm * x * x2 * (1.0L/3 + x2 * (1.0L/15 + x2 * (1.0L/35 + x2/63)));
    }
    return distance_cm * (x/2 + (x*x - 1)/4 * (std::log1p(x) - std::log1p(-x)));
}

void Write_Compact_Radial_Statistics(
    std::ostream& output, const std::vector<double>& edges_km,
    const std::array<const std::vector<double>*, 7>& histograms,
    const std::vector<double>& v2_square_sum,
    const OutputBlockCounts& captured_counts, const OutputBlockCounts& never_counts,
    long double temperature_factor)
{
    const unsigned long captured = std::accumulate(captured_counts.begin(), captured_counts.end(), 0UL);
    const unsigned long never = std::accumulate(never_counts.begin(), never_counts.end(), 0UL);
    const long double outer_rsun = edges_km.back() / R_SUN_KM;
    const std::array<Window, 6> windows{{
        {"solar_interior", 0, 1, false}, {"within_10_rsun", 0, 10, false},
        {"within_30_rsun", 0, 30, false}, {"solar_exterior", 1, outer_rsun, false},
        {"full_domain", 0, outer_rsun, false}, {"gamma_1_30_rsun_1au", 1, 30, true}}};
    const std::array<const char*, 4> components{{"residence_pair", "ever_pair", "never_pair", "ever_never_cross"}};
    constexpr std::size_t integral_count = 24;
    std::array<long double, integral_count> integrals{};
    std::array<OutputReplicates, integral_count> deleted_integrals{};
    std::vector<RadialStatistics> radial(edges_km.size() - 1);
    auto value = [](const std::vector<double>& histogram, std::size_t index) {
        return index < histogram.size() ? histogram[index] : 0.0;
    };
    for(std::size_t bin = 0; bin < radial.size(); ++bin)
    {
        auto& statistics = radial[bin];
        std::array<OutputMoments, RESIDENCE_JACKKNIFE_BLOCKS> blocks{};
        for(std::size_t moment = 0; moment < histograms.size(); ++moment)
            for(std::size_t block = 0; block < blocks.size(); ++block)
            {
                const double contribution = value(*histograms[moment], bin * blocks.size() + block);
                if(!std::isfinite(contribution) || contribution < 0)
                    throw std::runtime_error("nonfinite or negative bincount moment");
                blocks[block][moment] = contribution;
                statistics.totals[moment] += contribution;
            }
        statistics.central = Derived_Radial_Observables(statistics.totals, captured, never, temperature_factor);
        for(std::size_t block = 0; block < blocks.size(); ++block)
        {
            OutputMoments remaining{};
            // Sum retained blocks directly. Subtracting a dominant block from
            // the total can destroy the much smaller remaining second moment.
            for(std::size_t retained = 0; retained < blocks.size(); ++retained)
                if(retained != block)
                    for(std::size_t moment = 0; moment < remaining.size(); ++moment)
                        remaining[moment] += blocks[retained][moment];
            const auto deleted = Derived_Radial_Observables(remaining,
                captured - captured_counts[block], never - never_counts[block], temperature_factor);
            for(std::size_t observable = 0; observable < deleted.size(); ++observable)
                statistics.deleted[observable][block] = deleted[observable];
        }
        for(std::size_t window = 0; window < windows.size(); ++window)
        {
            const long double weight = Window_Weight(windows[window], edges_km[bin], edges_km[bin + 1]);
            if(weight == 0) continue; // Undefined populations outside a window do not affect it.
            for(std::size_t component = 0; component < components.size(); ++component)
            {
                const auto index = window * components.size() + component;
                integrals[index] += weight * statistics.central[component + 1];
                for(std::size_t block = 0; block < blocks.size(); ++block)
                    deleted_integrals[index][block] += weight * statistics.deleted[component + 1][block];
            }
        }
    }
    output << "# derived_observables_version = 1\n"
           << "# temperature_definition = kBT_eV=m_chi_eV/(3*c_km_s^2)*residence_v2dt/residence_dt\n"
           << "# pair_time_definition = (S*S-Q)/(N*(N-1)); residence/ever use N_residence_samples; never uses N_never_captured\n"
           << "# cross_time_definition = (ever_S/N_residence_samples)*(never_S/N_never_captured)\n"
           << "# nonlinear_covariance_definition = (63/64)*sum_b((theta_b-mean(theta))*(phi_b-mean(phi))); delete counts and all moments together\n"
           << "# nonlinear_uncertainty_undefined = nan if any of 64 deletions leaves an undefined observable; no blocks dropped\n"
           << "# integration_assumption = constant density within each bin; windows clipped to simulated radial domain\n"
           << "# rate_kernel = overlap_volume_cm3/bin_volume_cm3^2\n"
           << "# gamma_kernel = integral_overlap(r^2*angular_mean(1/distance^2)*dr)/bin_volume_cm3^2; isotropic emission; D=1_AU\n"
           << "# gamma_observer_distance_cm = " << AU_KM * 1e5 << '\n'
           << "# gamma_scope = annihilation event flux; no photon yield, spectrum, attenuation, or detector response\n"
           << "# physical_integral = 0.5*sigma_v_cm3_s*C_geom_per_s^2*(p2*A_ever+q2*A_never+2*pq*A_cross); residence-only uses p2*A_residence\n"
           << "# capture_factors = independent Capture: p2=Nc*(Nc-1)/(N*(N-1)); q2=Nu*(Nu-1)/(N*(N-1)); pq=Nc*Nu/(N*(N-1))\n"
           << "# integrated_covariance_scope = all listed components and windows, including cross-bin correlations; diagonal=SE^2\n"
           << "# integration_window_columns = name effective_r_low_rsun effective_r_high_rsun kernel\n"
           << "# integrated_stat_columns = window component units estimate SE\n"
           << "# integrated_stat_count = " << integral_count << '\n';
    for(std::size_t window = 0; window < windows.size(); ++window)
    {
        const auto& w = windows[window];
        output << "# integration_window_" << window << " = " << w.name << ' '
               << static_cast<double>(std::min(w.low_rsun, outer_rsun)) << ' '
               << static_cast<double>(std::min(w.high_rsun, outer_rsun)) << ' '
               << (w.flux ? "gamma_flux" : "rate") << '\n';
        for(std::size_t component = 0; component < components.size(); ++component)
        {
            const auto index = window * components.size() + component;
            output << "# integrated_stat_" << index << " = " << w.name << ' ' << components[component]
                   << (w.flux ? " s2_cm-5 " : " s2_cm-3 ");
            Write_Number(output, integrals[index]);
            output << ' ';
            Write_Number(output, std::sqrt(Block_Jackknife_Covariance(deleted_integrals[index], deleted_integrals[index])));
            output << '\n';
        }
    }
    for(std::size_t i = 0; i < integral_count; ++i)
        for(std::size_t j = i + 1; j < integral_count; ++j)
        {
            output << "# integrated_cov_" << i << '_' << j << " = ";
            Write_Number(output, Block_Jackknife_Covariance(deleted_integrals[i], deleted_integrals[j]));
            output << '\n';
        }
    output << "# columns = bin r_low_rsun r_high_rsun captured_residence_dt_sum_s captured_residence_dt_sq_sum_s2 captured_residence_dt_se_s captured_residence_v2dt_sum_km2_s captured_residence_v2dt_se_km2_s ever_captured_path_dt_sum_s ever_captured_path_dt_sq_sum_s2 ever_captured_path_dt_se_s never_captured_path_dt_sum_s never_captured_path_dt_sq_sum_s2 never_captured_path_dt_se_s captured_residence_v2dt_sq_sum_km4_s2 temperature_kBT_eV temperature_kBT_se_eV";
    for(const auto* component : components)
        output << ' ' << component << "_s2 " << component << "_se_s2";
    for(std::size_t i = 0; i < components.size(); ++i)
        for(std::size_t j = i + 1; j < components.size(); ++j)
            output << ' ' << components[i] << '_' << components[j] << "_cov_s4";
    output << '\n';
    for(std::size_t bin = 0; bin < radial.size(); ++bin)
    {
        const auto& statistics = radial[bin];
        const auto& t = statistics.totals;
        const double v2q = value(v2_square_sum, bin);
        if(!std::isfinite(v2q) || v2q < 0) throw std::runtime_error("invalid velocity second moment");
        output << bin << '\t' << edges_km[bin]/R_SUN_KM << '\t' << edges_km[bin + 1]/R_SUN_KM;
        for(const auto number : {t[0], t[1], static_cast<long double>(Moment_Sum_SE(t[0], t[1], captured)),
             t[2], static_cast<long double>(Moment_Sum_SE(t[2], v2q, captured)),
             t[3], t[4], static_cast<long double>(Moment_Sum_SE(t[3], t[4], captured)),
             t[5], t[6], static_cast<long double>(Moment_Sum_SE(t[5], t[6], never)), static_cast<long double>(v2q)})
        {
            output << '\t'; Write_Number(output, number);
        }
        for(std::size_t observable = 0; observable < statistics.central.size(); ++observable)
        {
            output << '\t'; Write_Number(output, statistics.central[observable]);
            output << '\t'; Write_Number(output, std::sqrt(Block_Jackknife_Covariance(statistics.deleted[observable], statistics.deleted[observable])));
        }
        for(std::size_t i = 1; i < statistics.central.size(); ++i)
            for(std::size_t j = i + 1; j < statistics.central.size(); ++j)
            {
                output << '\t'; Write_Number(output, Block_Jackknife_Covariance(statistics.deleted[i], statistics.deleted[j]));
            }
        output << '\n';
    }
}
}

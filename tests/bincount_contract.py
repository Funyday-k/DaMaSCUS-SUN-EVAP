"""Independent checks of compact moments, derived observables and integral covariance."""
from __future__ import annotations

import itertools
import math
from pathlib import Path

COMPONENTS = ('residence_pair', 'ever_pair', 'never_pair', 'ever_never_cross')
COLUMNS = (
    'bin r_low_rsun r_high_rsun captured_residence_dt_sum_s '
    'captured_residence_dt_sq_sum_s2 captured_residence_dt_se_s '
    'captured_residence_v2dt_sum_km2_s captured_residence_v2dt_se_km2_s '
    'ever_captured_path_dt_sum_s ever_captured_path_dt_sq_sum_s2 ever_captured_path_dt_se_s '
    'never_captured_path_dt_sum_s never_captured_path_dt_sq_sum_s2 never_captured_path_dt_se_s '
    'captured_residence_v2dt_sq_sum_km4_s2 temperature_kBT_eV temperature_kBT_se_eV'
).split()
for component in COMPONENTS:
    COLUMNS.extend((f'{component}_s2', f'{component}_se_s2'))
COLUMNS.extend(f'{a}_{b}_cov_s4' for a, b in itertools.combinations(COMPONENTS, 2))
MOMENT_COLUMNS = (3, 4, 6, 8, 9, 11, 12)
SE_COLUMNS = (5, 7, 10, 13)


def read_table(path: Path) -> tuple[dict[str, str], dict[int, tuple[int, ...]], list[list[float]]]:
    header: dict[str, str] = {}
    blocks: dict[int, tuple[int, ...]] = {}
    rows: list[list[float]] = []
    with path.open() as stream:
        for line in stream:
            if line.startswith('# '):
                key, separator, value = line[2:].strip().partition(' = ')
                if not separator:
                    continue
                if key == 'block_count':
                    block, *counts = map(int, value.split())
                    assert block not in blocks and len(counts) == 4
                    blocks[block] = tuple(counts)
                else:
                    assert key not in header, key
                    header[key] = value
            elif line.strip():
                rows.append(list(map(float, line.split())))
    assert set(blocks) == set(range(64))
    return header, blocks, rows


def read_bincount(path: Path) -> tuple[dict[str, str], dict[int, tuple[int, ...]], list[list[float]]]:
    header, blocks, rows = read_table(path)
    assert header['bincount_format_version'] == '3'
    assert header['derived_observables_version'] == '1'
    assert header['columns'].split() == COLUMNS
    assert header['output_significant_digits'] == '17'
    assert header['failed_history_policy'] == 'discard_and_replace'
    assert header['radial_covariance_available'] == 'selected_observables_and_integrals'
    assert header['capture_normalization_uncertainty_included'] == 'false'
    assert len(rows) == int(header['radial_bins'])
    assert all(len(row) == len(COLUMNS) == 31 for row in rows)
    assert path.stat().st_size < 1_500_000
    for index, key in enumerate(('N_injected', 'N_ever_captured', 'N_never_captured', 'N_residence_samples')):
        assert sum(counts[index] for counts in blocks.values()) == int(header[key])
    assert int(header['integrated_stat_count']) == 24
    outer = float(header['R_remove_rsun'])
    expected_windows = [('solar_interior', 0, 1, 'rate'), ('within_10_rsun', 0, 10, 'rate'),
                        ('within_30_rsun', 0, 30, 'rate'), ('solar_exterior', 1, outer, 'rate'),
                        ('full_domain', 0, outer, 'rate'), ('gamma_1_30_rsun_1au', 1, 30, 'gamma_flux')]
    for i, (name, low, high, kernel) in enumerate(expected_windows):
        actual_name, actual_low, actual_high, actual_kernel = header[f'integration_window_{i}'].split()
        assert (actual_name, actual_kernel) == (name, kernel)
        assert math.isclose(float(actual_low), min(low, outer), rel_tol=1e-14)
        assert math.isclose(float(actual_high), min(high, outer), rel_tol=1e-14)
    assert math.isclose(float(header['gamma_observer_distance_cm']), 14959787070000, rel_tol=2e-15)
    for i in range(24):
        name, component, units, estimate, error = header[f'integrated_stat_{i}'].split()
        assert component == COMPONENTS[i % 4]
        assert units == ('s2_cm-5' if i >= 20 else 's2_cm-3')
        assert name == header[f'integration_window_{i // 4}'].split()[0]
        assert not math.isinf(float(estimate)) and not math.isinf(float(error))
        for j in range(i + 1, 24):
            assert not math.isinf(float(header[f'integrated_cov_{i}_{j}']))
    return header, blocks, rows


def close(actual: float, expected: float, *, scale: float = 0.0) -> None:
    if math.isnan(expected):
        assert math.isnan(actual), (actual, expected)
    else:
        assert math.isclose(actual, expected, rel_tol=2e-8, abs_tol=2e-13 * scale), (actual, expected)


def sum_se(total: float, square_total: float, count: int) -> float:
    return math.sqrt(max(0, count * square_total - total * total) / (count - 1)) if count >= 2 else math.nan


def observables(m: list[float], captured: int, never: int, factor: float) -> list[float]:
    def pair(s: float, q: float, n: int) -> float:
        return max(0, s*s-q) / (n*(n-1)) if n >= 2 else math.nan
    return [factor * m[2] / m[0] if m[0] else math.nan,
            pair(m[0], m[1], captured), pair(m[3], m[4], captured), pair(m[5], m[6], never),
            (m[3] / captured) * (m[5] / never) if captured and never else math.nan]


def observable_scales(m: list[float], captured: int, never: int, factor: float) -> list[float]:
    """Use the pre-cancellation scale when checking S^2-Q near zero."""
    result = observables(m, captured, never, factor)
    for i, total, count in ((1, m[0], captured), (2, m[3], captured), (3, m[5], never)):
        result[i] = total * total / (count * (count - 1)) if count >= 2 else 0
    return [abs(x) if math.isfinite(x) else 0 for x in result]


def covariance(a: list[float], b: list[float]) -> float:
    if any(not math.isfinite(x) for x in a + b):
        return math.nan
    mean_a, mean_b = math.fsum(a)/64, math.fsum(b)/64
    return 63/64 * math.fsum((x - mean_a)*(y - mean_b) for x, y in zip(a, b))


def validate_moments(path: Path) -> dict[str, str]:
    header, blocks, rows = read_bincount(path)
    excluded = 0
    for injected, captured, never, residence in blocks.values():
        assert 0 <= residence <= captured and captured + never <= injected
        excluded += injected - residence - never
    assert excluded == int(header['N_excluded_trajectories'])
    assert excluded == int(header['numerical_failures']) + int(header['computational_failures'])
    captured, never = int(header['N_residence_samples']), int(header['N_never_captured'])
    factor = float(header['m_chi_GeV']) * 1e9 / (3 * 299792.458**2)
    previous_high = 0.0
    for index, row in enumerate(rows):
        assert row[0] == index
        assert all(math.isfinite(row[j]) and row[j] >= 0 for j in (0, 1, 2, *MOMENT_COLUMNS, 14))
        assert row[2] > row[1] and row[1] == previous_high
        previous_high = row[2]
        for first, second, error, count in ((3, 4, 5, captured), (6, 14, 7, captured), (8, 9, 10, captured), (11, 12, 13, never)):
            total, square_total = row[first], row[second]
            tolerance = 1e-10 * max(1.0, total * total, count * square_total)
            assert square_total <= total * total + tolerance
            assert total * total <= count * square_total + tolerance
            if count == 0:
                assert total == square_total == 0.0
            close(row[error], sum_se(total, square_total, count), scale=abs(total))
        moments = [row[j] for j in MOMENT_COLUMNS]
        expected = observables(moments, captured, never, factor)
        scales = observable_scales(moments, captured, never, factor)
        for i, number in enumerate(expected):
            close(row[15 + 2*i], number, scale=scales[i])
        assert not any(math.isinf(value) for value in row)
    assert math.isclose(previous_high, float(header['R_remove_rsun']), rel_tol=1e-12)
    return header


def validate_completed(path: Path) -> dict[str, str]:
    header = validate_moments(path)
    assert header['target_reached'] == 'true'
    assert int(header['N_residence_samples']) == int(header['requested_captured'])
    return header


def integration_weights(row: list[float], header: dict[str, str]) -> list[float]:
    rsun = float(header['R_sun_km']) * 1e5
    low, high = row[1] * rsun, row[2] * rsun
    volume = 4 * math.pi/3 * (high**3 - low**3)
    distance = float(header['gamma_observer_distance_cm'])
    # Independent convergent series for the detector integral; r/D <= 30 R_sun/AU.
    def detector(r: float) -> float:
        x = r / distance
        return distance * math.fsum(x**(2*k+3)/((2*k+1)*(2*k+3)) for k in range(20))
    weights = []
    for i in range(6):
        _, begin, end, kernel = header[f'integration_window_{i}'].split()
        a, b = max(low, float(begin)*rsun), min(high, float(end)*rsun)
        integral = 0 if b <= a else (detector(b)-detector(a) if kernel == 'gamma_flux' else 4*math.pi/3*(b**3-a**3))
        weights.append(integral / volume**2)
    return weights


def validate_against_blocks(compact_path: Path, diagnostic_path: Path) -> None:
    """Recompute nonlinear replicates and all 24 integral covariances from raw blocks."""
    header, counts, rows = read_bincount(compact_path)
    diagnostic, diagnostic_counts, full = read_table(diagnostic_path)
    assert diagnostic['radial_blocks_format_version'] == '1'
    assert counts == diagnostic_counts
    bins = int(header['radial_bins'])
    assert len(full) == 64 * bins
    for index, block_row in enumerate(full):
        assert len(block_row) == 11
        assert block_row[:2] == [index // bins, index % bins]
    captured = sum(counts[b][3] for b in range(64))
    never = sum(counts[b][2] for b in range(64))
    factor = float(header['m_chi_GeV']) * 1e9 / (3 * 299792.458**2)
    integrated = [0.0] * 24
    integrated_scales = [0.0] * 24
    integrated_deleted = [[0.0] * 64 for _ in range(24)]
    for bin_index, row in enumerate(rows):
        bin_blocks = [full[block * bins + bin_index][4:] for block in range(64)]
        totals = [math.fsum(b[moment] for b in bin_blocks) for moment in range(7)]
        for column, total in zip(MOMENT_COLUMNS, totals):
            close(row[column], total, scale=abs(total))
        central = observables(totals, captured, never, factor)
        scales = observable_scales(totals, captured, never, factor)
        deleted = [observables([math.fsum(bin_blocks[other][m] for other in range(64) if other != b) for m in range(7)],
                    captured-counts[b][3], never-counts[b][2], factor) for b in range(64)]
        transposed = [list(values) for values in zip(*deleted)]
        for i in range(5):
            close(row[15+2*i], central[i], scale=scales[i])
            close(row[16+2*i], math.sqrt(covariance(transposed[i], transposed[i])), scale=scales[i])
        for column, (i, j) in enumerate(itertools.combinations(range(1, 5), 2), 25):
            close(row[column], covariance(transposed[i], transposed[j]), scale=scales[i]*scales[j])
        for window, weight in enumerate(integration_weights(row, header)):
            if weight == 0:
                continue
            for component in range(4):
                index = window * 4 + component
                integrated[index] += weight * central[component+1]
                integrated_scales[index] += weight * scales[component+1]
                for b in range(64):
                    integrated_deleted[index][b] += weight * deleted[b][component+1]
    for i in range(24):
        *_, estimate, error = header[f'integrated_stat_{i}'].split()
        close(float(estimate), integrated[i], scale=integrated_scales[i])
        close(float(error), math.sqrt(covariance(integrated_deleted[i], integrated_deleted[i])), scale=integrated_scales[i])
        for j in range(i+1, 24):
            close(float(header[f'integrated_cov_{i}_{j}']), covariance(integrated_deleted[i], integrated_deleted[j]), scale=integrated_scales[i]*integrated_scales[j])

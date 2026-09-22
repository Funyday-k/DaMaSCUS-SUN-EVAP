"""Read format-3 totals and validate radial moments and jackknife SE (s, km, R_sun)."""
from __future__ import annotations

import math
from pathlib import Path

COLUMNS = (
    'bin r_low_rsun r_high_rsun captured_residence_dt_sum_s '
    'captured_residence_dt_sq_sum_s2 captured_residence_dt_se_s '
    'captured_residence_v2dt_sum_km2_s captured_residence_v2dt_se_km2_s '
    'ever_captured_path_dt_sum_s ever_captured_path_dt_sq_sum_s2 ever_captured_path_dt_se_s '
    'never_captured_path_dt_sum_s never_captured_path_dt_sq_sum_s2 never_captured_path_dt_se_s'
).split()
MOMENT_COLUMNS = (3, 4, 6, 8, 9, 11, 12)
SE_COLUMNS = (5, 7, 10, 13)


def read_table(path: Path) -> tuple[dict[str, str], dict[int, tuple[int, ...]], list[list[float]]]:
    """Read comment metadata, block counts, and a numeric radial table."""
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
    """Require one row per bin, fourteen columns, and small production output."""
    header, blocks, rows = read_table(path)
    assert header['bincount_format_version'] == '3'
    assert header['columns'].split() == COLUMNS
    assert header['output_significant_digits'] == '10'
    assert header['jackknife_se_scale'] == 'total_sum_at_fixed_population_count'
    assert header['radial_covariance_available'] == 'false'
    assert header['capture_normalization_uncertainty_included'] == 'false'
    assert len(rows) == int(header['radial_bins'])
    assert all(len(row) == 14 for row in rows)
    assert path.stat().st_size < 600_000
    for index, key in enumerate(('N_injected', 'N_ever_captured', 'N_never_captured', 'N_residence_samples')):
        assert sum(counts[index] for counts in blocks.values()) == int(header[key])
    return header, blocks, rows


def validate_moments(path: Path) -> dict[str, str]:
    """Check conditional-population moments; permit NaN only for undefined SE."""
    header, blocks, rows = read_bincount(path)
    excluded = 0
    for injected, captured, never, residence in blocks.values():
        assert 0 <= residence <= captured and captured + never <= injected
        excluded += injected - residence - never
    assert excluded == int(header['N_excluded_trajectories'])
    assert excluded == int(header['numerical_failures']) + int(header['computational_failures'])
    captured, never = int(header['N_residence_samples']), int(header['N_never_captured'])
    previous_high = 0.0
    for index, row in enumerate(rows):
        assert row[0] == index
        assert all(math.isfinite(row[j]) and row[j] >= 0 for j in (0, 1, 2, *MOMENT_COLUMNS))
        assert row[2] > row[1] and row[1] == previous_high
        previous_high = row[2]
        for first, second, count in ((3, 4, captured), (8, 9, captured), (11, 12, never)):
            total, square_total = row[first], row[second]
            tolerance = 2e-9 * max(1.0, total * total, count * square_total)
            assert square_total <= total * total + tolerance
            assert total * total <= count * square_total + tolerance
            if count == 0:
                assert total == square_total == 0.0
            elif count == 1:
                assert math.isclose(square_total, total * total, rel_tol=2e-9, abs_tol=1e-12)
        for column, population, count_index in ((5, captured, 3), (7, captured, 3), (10, captured, 3), (13, never, 2)):
            available = population >= 2 and max(counts[count_index] for counts in blocks.values()) < population
            assert (math.isfinite(row[column]) and row[column] >= 0) if available else math.isnan(row[column])
    assert math.isclose(previous_high, float(header['R_remove_rsun']), rel_tol=1e-9)
    return header


def validate_completed(path: Path) -> dict[str, str]:
    """Require the complete captured target, including runs with failed histories."""
    header = validate_moments(path)
    assert header['target_reached'] == 'true'
    assert int(header['N_residence_samples']) == int(header['requested_captured'])
    return header


def validate_against_blocks(compact_path: Path, diagnostic_path: Path) -> None:
    """Independently recompute all seven totals and four SE from full-precision blocks."""
    header, counts, rows = read_bincount(compact_path)
    diagnostic, diagnostic_counts, full = read_table(diagnostic_path)
    assert diagnostic['radial_blocks_format_version'] == '1'
    assert counts == diagnostic_counts
    bins = int(header['radial_bins'])
    assert len(full) == 64 * bins
    for index, block_row in enumerate(full):
        assert len(block_row) == 11
        assert block_row[:2] == [index // bins, index % bins]
    for bin_index, row in enumerate(rows):
        bin_blocks = [full[block * bins + bin_index] for block in range(64)]
        for moment, column in enumerate(MOMENT_COLUMNS):
            total = math.fsum(block_row[4 + moment] for block_row in bin_blocks)
            assert math.isclose(row[column], total, rel_tol=6e-10, abs_tol=1e-12)
        for moment, column, count_index in ((0, 5, 3), (2, 7, 3), (3, 10, 3), (5, 13, 2)):
            population = sum(counts[block][count_index] for block in range(64))
            if population < 2 or max(counts[block][count_index] for block in range(64)) == population:
                assert math.isnan(row[column])
                continue
            total = math.fsum(block_row[4 + moment] for block_row in bin_blocks)
            deleted_means = [(total - bin_blocks[block][4 + moment]) / (population - counts[block][count_index]) for block in range(64)]
            center = math.fsum(deleted_means) / 64
            expected = population * math.sqrt(63 / 64 * math.fsum((mean - center) ** 2 for mean in deleted_means))
            assert math.isclose(row[column], expected, rel_tol=2e-8, abs_tol=2e-13 * max(1, abs(total)))

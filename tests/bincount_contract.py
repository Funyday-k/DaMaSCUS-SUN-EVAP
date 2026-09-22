"""Read and validate format-2 radial sufficient statistics (seconds, km, R_sun)."""
from __future__ import annotations

import math
from pathlib import Path

COLUMNS = (
    'block bin r_low_rsun r_high_rsun captured_residence_dt_sum_s '
    'captured_residence_dt_sq_sum_s2 captured_residence_v2dt_sum_km2_s '
    'ever_captured_path_dt_sum_s ever_captured_path_dt_sq_sum_s2 '
    'never_captured_path_dt_sum_s never_captured_path_dt_sq_sum_s2'
).split()


def read_bincount(path: Path) -> tuple[dict[str, str], dict[int, tuple[int, ...]], list[list[float]]]:
    """Read header, integer block counts, and the eleven-column radial table."""
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
                values = list(map(float, line.split()))
                assert len(values) == 11
                rows.append(values)
    assert header['bincount_format_version'] == '2'
    assert header['columns'].split() == COLUMNS
    assert set(blocks) == set(range(64))
    assert len(rows) == 64 * int(header['radial_bins'])
    for index, key in enumerate(('N_injected', 'N_ever_captured', 'N_never_captured', 'N_residence_samples')):
        assert sum(counts[index] for counts in blocks.values()) == int(header[key])
    return header, blocks, rows


def validate_moments(path: Path) -> dict[str, str]:
    """Check complete-history normalization and moment bounds, including runs with failures."""
    header, blocks, rows = read_bincount(path)
    assert int(header['numerical_failures']) >= 0 and int(header['computational_failures']) >= 0
    excluded = 0
    for injected, captured, never, residence in blocks.values():
        assert 0 <= residence <= captured and captured + never <= injected
        excluded += injected - residence - never
    assert excluded == int(header['N_excluded_trajectories'])
    assert excluded == int(header['numerical_failures']) + int(header['computational_failures'])
    bins = int(header['radial_bins'])
    previous_high: dict[int, float] = {}
    for index, row in enumerate(rows):
        block, bin_index = int(row[0]), int(row[1])
        assert row[0] == block == index // bins and row[1] == bin_index == index % bins
        assert all(math.isfinite(value) and value >= 0 for value in row)
        assert row[3] > row[2]
        assert row[2] == previous_high.get(block, 0.0)
        previous_high[block] = row[3]
        for first, second, count in ((4, 5, blocks[block][3]), (7, 8, blocks[block][3]), (9, 10, blocks[block][2])):
            total, square_total = row[first], row[second]
            tolerance = 1e-10 * max(1.0, total * total, count * square_total)
            assert square_total <= total * total + tolerance
            assert total * total <= count * square_total + tolerance
            if count == 0:
                assert total == square_total == 0.0
            elif count == 1:
                assert math.isclose(square_total, total * total, rel_tol=1e-12, abs_tol=1e-12)
    assert all(math.isclose(radius, float(header['R_remove_rsun']), rel_tol=1e-12) for radius in previous_high.values())
    return header


def validate_completed(path: Path) -> dict[str, str]:
    """Require the target count, without requiring zero failed histories."""
    header = validate_moments(path)
    assert header['target_reached'] == 'true'
    assert int(header['N_residence_samples']) == int(header['requested_captured'])
    return header

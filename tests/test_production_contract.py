#!/usr/bin/env python3
"""Exercise scientific output, rejection, snapshots, and explicit local diagnostics."""
from __future__ import annotations
import argparse
import csv
import io
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile

from bincount_contract import read_bincount, validate_completed, validate_moments, validate_against_blocks


def validate_snapshots(directory: Path, final_header: dict[str, str], ranks: int) -> None:
    """Validate cumulative reports and honest missed checkpoints on any ranks."""
    assert ranks > 0
    assert not list(directory.rglob('*_evaporation_times.txt'))
    snapshots = sorted(directory.glob('snapshot_*s.txt'),
                       key=lambda path: int(re.fullmatch(r'snapshot_(\d+)s\.txt', path.name).group(1)))
    # A valid run can finish before the first one-second heartbeat.
    for path in snapshots:
        header: dict[str, str] = {}
        rows: list[list[float]] = []
        for line in path.read_text().splitlines():
            if line.startswith('#'):
                key, separator, value = line[1:].partition('=')
                if separator:
                    header[key.strip()] = value.strip()
            elif line.strip():
                row = [float(value) for value in line.split()]
                assert len(row) == 11, (path, row)
                assert all(math.isfinite(value) and value >= 0.0 for value in row)
                rows.append(row)
        status = header['snapshot_status']
        assert status in ('merged', 'partial', 'missed_writer_delay'), (path, header)
        target = int(header['snapshot_target_wall_time_s'])
        interval = float(header['snapshot_interval_s'])
        assert target == int(re.fullmatch(r'snapshot_(\d+)s\.txt', path.name).group(1))
        assert math.isfinite(interval) and interval > 0.0 and interval.is_integer()
        assert target > 0 and target % int(interval) == 0
        assert int(header['snapshot_run_id']) >= 0
        if status == 'missed_writer_delay':
            # The writer emits only these fields: no reconstructed populations
            # or ready/missing rank lists exist for an unwritten checkpoint.
            assert set(header) == {'snapshot_status', 'snapshot_run_id', 'snapshot_target_wall_time_s',
                                   'snapshot_interval_s', 'actual_write_wall_time_s', 'note'}, (path, header)
            actual = float(header['actual_write_wall_time_s'])
            assert math.isfinite(actual) and actual > target
            assert header['note'] == 'heartbeat woke too late; no state rewind is available'
            assert not rows, path
            continue
        assert header['format_version'] == '12'
        ready = [int(value) for value in header['ready_ranks'].split(',')] if header['ready_ranks'] else []
        missing = [int(value) for value in header['missing_ranks'].split(',')] if header['missing_ranks'] else []
        assert ready == sorted(set(ready)) and missing == sorted(set(missing))
        assert set(ready).isdisjoint(missing)
        assert sorted(ready + missing) == list(range(ranks))
        if status == 'merged':
            assert ready == list(range(ranks)) and not missing
        else:
            assert ready and missing
            assert header['WARNING'] == 'partial snapshot; not all ranks reached this checkpoint'
        assert rows and [row[0] for row in rows] == list(range(len(rows)))
        assert header['radial_range'] == 'dynamic_complete_shells'
        assert header['exterior_grid'] == 'geometric_width_unbounded'
        assert math.isclose(rows[-1][2], float(header['R_output_rsun']), rel_tol=1e-10)
        assert float(header['R_output_rsun']) >= float(header['R_output_min_rsun'])
        assert header['R_path_reference_rsun'] == header['R_output_rsun']
        for first, second in zip(rows[1100:], rows[1101:]):
            assert math.isclose((second[2]-second[1])/(first[2]-first[1]), 1.02, rel_tol=1e-7)
        total = int(header['total_trajectories'])
        captured = int(header['captured_particles'])
        classified = int(header['valid_trajectories'])
        never = int(header['never_captured_bincount_samples'])
        assert 0 <= captured <= classified <= total <= int(final_header['N_injected'])
        assert captured <= int(final_header['N_ever_captured'])
        assert never == classified - captured
        assert never <= int(final_header['N_never_captured'])
        assert int(header['unresolved_not_captured_trajectories']) == total - classified
        provisional = int(header['in_progress_bincount_captured_samples'])
        residence = int(header['residence_bincount_samples'])
        assert 0 <= provisional <= len(ready)
        assert provisional <= residence <= captured + provisional
        assert int(header['completed_captured_residence_bincount_samples']) == residence - provisional
        assert int(header['captured_residence_bincount_samples']) == residence
        assert header['never_captured_path_population'] == 'completed_physically_classified_never_captured_outward_escape'
        assert header['in_progress_uncaptured_path_included'] == '0'
        if not never:
            assert all(row[7:] == [0.0] * 4 for row in rows), path


def main() -> None:
    """Run small configurations with the same scientific contract on one or more ranks."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--program', type=Path, required=True)
    parser.add_argument('--replay', type=Path, required=True)
    parser.add_argument('--mpiexec')
    parser.add_argument('--numproc-flag', default='-np')
    parser.add_argument('--ranks', type=int, default=1)
    args = parser.parse_args()
    template = (Path(__file__).resolve().parents[1] / 'examples/quickstart.cfg').read_text()
    launcher = [args.mpiexec, args.numproc_flag, str(args.ranks)] if args.mpiexec else []
    with tempfile.TemporaryDirectory(prefix='damascus-production-contract-') as tmp:
        root = Path(tmp)

        def run(name: str, options: tuple[str, ...] = (), **changes: str) -> tuple[subprocess.CompletedProcess, Path, Path]:
            """Launch an external cfg; physical input units are those of quickstart.cfg."""
            output = root / name
            settings = {'run_mode': '"Capture"', 'sample_size': '128', 'max_trajectories': '0',
                        'DM_cross_section_nucleon': '1.0e-32', 'fixed_seed': '20260922',
                        'max_trajectory_wall_time_sec': '0.0', 'output_dir': json.dumps(str(output) + '/'),
                        'snapshot_enabled': 'true', 'snapshot_interval': '1', **changes}
            text = template
            for key, value in settings.items():
                pattern = rf'(?m)^{key}\s*=.*?;'
                text = re.sub(pattern, f'{key} = {value};', text) if re.search(pattern, text) else text + f'\n{key} = {value};\n'
            config = root / f'{name}.cfg'
            config.write_text(text)
            result = subprocess.run(launcher + [str(args.program), str(config), *options], cwd=root,
                                    text=True, capture_output=True, timeout=60)
            return result, output, config

        def capture_record(result: subprocess.CompletedProcess) -> dict:
            """Require stdout to contain exactly one JSON record, with all logs on stderr."""
            lines = result.stdout.splitlines()
            assert len(lines) == 1 and lines[0].startswith('CAPTURE_RESULT_JSON='), result.stdout
            record = json.loads(lines[0].split('=', 1)[1])
            assert record['capture_result_schema'] == 3 and record['mpi_ranks'] == args.ranks
            assert len(record['blocks']) == 64
            assert record['N_inj'] == sum(row[0] for row in record['blocks'])
            assert record['N_capt'] == sum(row[1] for row in record['blocks'])
            return record

        result, output, _ = run('capture')
        assert result.returncode == 0, result.stderr
        assert capture_record(result)['N_inj'] == 128
        assert not output.exists()
        result, output, _ = run('failed_capture', max_trajectories='128',
                                max_trajectory_wall_time_sec='1e-12')
        assert result.returncode == 2, result.stderr
        failed = capture_record(result)
        assert failed['target_reached'] is False and failed['N_computational_failures'] == 128
        assert failed['N_attempted'] == 128 and failed['N_inj'] == 0
        assert failed['N_unclassified'] == 128 and not output.exists()

        for name, snapshot_enabled in (('transport', 'true'), ('no_snapshot', 'false')):
            result, output, _ = run(name, run_mode='"Parameter point"', sample_size='64', snapshot_enabled=snapshot_enabled)
            assert result.returncode == 0, result.stdout + result.stderr
            product = output / 'results_-2.000000_-32.000000'
            assert {path.name for path in product.iterdir()} == {'bincount.tsv', 'snapshot'}
            header = validate_completed(product / 'bincount.tsv')
            assert int(header['N_ever_captured']) == 64
            assert (product / 'snapshot').is_dir()
            assert not list((product / 'snapshot').rglob('*_evaporation_times.txt'))
            assert 'source_sha256' not in header and 'compiler' not in header
            assert json.loads(header['SHM_vObserver_km_s']) == [11.1, 232.2, 7.3]
            assert json.loads(header['DM_relative_couplings']) == [1, 1]
            assert header['DM_form_factor'] == '"Contact"' and 'DM_mediator_mass_MeV' not in header

        result, output, _ = run('failed_transport', run_mode='"Parameter point"', sample_size='4', max_trajectories='12', max_trajectory_wall_time_sec='1e-12')
        assert result.returncode == 2, result.stdout + result.stderr
        product = output / 'results_-2.000000_-32.000000'
        assert {path.name for path in product.iterdir()} == {'bincount.tsv', 'snapshot'}
        header, _, _ = read_bincount(product / 'bincount.tsv')
        assert header['target_reached'] == 'false' and int(header['computational_failures']) == 12
        assert int(header['N_never_captured']) == 0 and int(header['N_unclassified']) == 12
        assert int(header['N_injected']) == 12
        validate_moments(product / 'bincount.tsv')
        sentinel = product / 'snapshot' / 'keep'
        sentinel.write_text('previous snapshot')
        saved = (product / 'bincount.tsv').read_bytes()
        result, _, _ = run('failed_transport', run_mode='"Parameter point"')
        assert result.returncode != 0 and 'not empty' in result.stderr
        assert sentinel.read_text() == 'previous snapshot' and (product / 'bincount.tsv').read_bytes() == saved

        # A smaller initial extent must never truncate statistics. The same
        # streams retain all recorded shells, moments and completion decisions.
        result, output, _ = run('small_window', run_mode='"Parameter point"', sample_size='64',
                                output_radius_rsun='2.05', snapshot_enabled='false')
        assert result.returncode == 0, result.stdout + result.stderr
        small_path = output / 'results_-2.000000_-32.000000/bincount.tsv'
        small_header, small_blocks, small_rows = read_bincount(small_path)
        validate_completed(small_path)
        full_header, full_blocks, full_rows = read_bincount(root / 'no_snapshot/results_-2.000000_-32.000000/bincount.tsv')
        assert small_blocks == full_blocks
        for key in ('N_injected', 'N_ever_captured', 'N_never_captured', 'N_residence_samples',
                    'N_excluded_trajectories', 'numerical_failures', 'computational_failures', 'stop_reason'):
            assert small_header[key] == full_header[key], key
        assert float(small_header['R_output_rsun']) > 2.05
        assert len(small_rows) <= len(full_rows)
        assert float(small_header['R_path_reference_rsun']) == small_rows[-1][2]
        assert float(full_header['R_path_reference_rsun']) == full_rows[-1][2]
        assert all(a == b or (math.isnan(a) and math.isnan(b))
                   for small, full in zip(small_rows, full_rows) for a, b in zip(small, full))

        # This fixed stream includes a numerical failure before completing 1000
        # captures. It must continue and retain correct conditional denominators.
        result, output, _ = run('continue_numerical', run_mode='"Parameter point"', sample_size='1000',
                                max_trajectories='100000', DM_mass='0.1', DM_cross_section_nucleon='1e-36',
                                fixed_seed='23260931', rate_radius_points='1000', rate_speed_points='256',
                                rate_max_speed='0.02', snapshot_enabled='true')
        assert result.returncode == 0, result.stdout + result.stderr
        header = validate_completed(output / 'results_-1.000000_-36.000000/bincount.tsv')
        assert int(header['numerical_failures']) > 0
        assert int(header['N_never_captured']) > 0
        validate_snapshots(output / 'results_-1.000000_-36.000000/snapshot', header, args.ranks)

        # Captured-but-truncated histories cannot enter the residence denominator.
        result, output, _ = run('captured_truncation', run_mode='"Parameter point"', sample_size='64',
                                max_trajectories='128', maximum_number_of_scatterings='1', snapshot_enabled='false')
        assert result.returncode == 2, result.stdout + result.stderr
        header = validate_moments(output / 'results_-2.000000_-32.000000/bincount.tsv')
        assert int(header['N_ever_captured']) > int(header['N_residence_samples'])
        assert int(header['N_injected']) == 128

        # Explicit debug runs stay outside the normal scientific result directory.
        result, output, config = run('local', ('--diagnostic',), run_mode='"Parameter point"',
                                     sample_size='64', snapshot_enabled='false')
        assert result.returncode == 0, result.stdout + result.stderr
        local = output / 'diagnostics/results_-2.000000_-32.000000'
        assert (local / 'radial_blocks.tsv').is_file() and not (local / 'bincount.tsv').exists()
        if args.ranks == 1:
            validate_against_blocks(root / 'no_snapshot/results_-2.000000_-32.000000/bincount.tsv',
                                    local / 'radial_blocks.tsv')
        rows = list(csv.DictReader(io.StringIO((local / 'diagnostic_trajectory_summary.tsv').read_text()), delimiter='\t'))
        assert rows and any(row['rng_state_before_simulation'] for row in rows)
        complete = [row for row in rows if row['event_observed'] == '1']
        assert complete and len(complete) == int(full_header['N_residence_samples'])
        assert all(row['termination_reason'] == 'physical_escape' for row in complete)
        assert not any(row['termination_reason'] == 'outer_domain_removal' for row in rows)
        for row in complete:
            lifetime = float(row['lifetime_validated_escape_s'])
            physical_time = float(row['time_inside_sun_s']) + float(row['time_outside_sun_s'])
            assert math.isclose(physical_time, lifetime, rel_tol=1e-8, abs_tol=1e-5)
        # Full radial residence retains every bound exterior waiting interval.
        # Compare rounded production totals to full-precision lifetimes.
        recorded = math.fsum(row[3] for row in full_rows)
        lifetime_total = math.fsum(float(row['lifetime_validated_escape_s']) for row in complete)
        assert math.isclose(recorded, lifetime_total, rel_tol=2e-3)
        assert len((local / 'trajectory_events.tsv').read_text().splitlines()) > 1
        result, output, config = run('local_invalid', ('--diagnostic',), run_mode='"Parameter point"', sample_size='1',
                                     max_trajectories='1', maximum_number_of_scatterings='1', snapshot_enabled='false')
        assert result.returncode == 0, result.stdout + result.stderr
        ledger = output / 'diagnostics/results_-2.000000_-32.000000/invalid_trajectories.tsv'
        rows = list(csv.DictReader((line for line in ledger.read_text().splitlines() if not line.startswith('#')), delimiter='\t'))
        assert rows and rows[0]['rng_state_before_simulation']
        replay = subprocess.run([str(args.replay), str(config), str(ledger), rows[0]['rank'], rows[0]['trajectory_id']],
                                cwd=root, capture_output=True, text=True, timeout=20)
        assert replay.returncode == 0, replay.stdout + replay.stderr
        assert 'original_termination_reason=max_scatterings' in replay.stdout and 'replay_scatterings=1' in replay.stdout
        result, _, _ = run('obsolete', production_mode='false')
        assert result.returncode != 0 and 'no longer a cfg setting' in result.stderr
        print(f'{args.ranks} rank(s): scientific output, continuation after errors, snapshot preservation and explicit replay passed')


if __name__ == '__main__':
    main()

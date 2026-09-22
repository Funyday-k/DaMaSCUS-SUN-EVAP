#!/usr/bin/env python3
"""Exercise scientific output, rejection, snapshots, and explicit local diagnostics."""
from __future__ import annotations
import argparse
import csv
import io
import json
from pathlib import Path
import re
import subprocess
import tempfile

from bincount_contract import read_bincount, validate_completed, validate_moments


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
            assert record['capture_result_schema'] == 2 and record['mpi_ranks'] == args.ranks
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

        # This fixed stream includes a numerical failure before completing 1000
        # captures. It must continue and retain correct conditional denominators.
        result, output, _ = run('continue_numerical', run_mode='"Parameter point"', sample_size='1000',
                                max_trajectories='100000', DM_mass='0.1', DM_cross_section_nucleon='1e-36',
                                fixed_seed='23260931', rate_radius_points='1000', rate_speed_points='256',
                                rate_max_speed='0.02', snapshot_enabled='false')
        assert result.returncode == 0, result.stdout + result.stderr
        header = validate_completed(output / 'results_-1.000000_-36.000000/bincount.tsv')
        assert int(header['numerical_failures']) > 0

        # Captured-but-truncated histories cannot enter the residence denominator.
        result, output, _ = run('captured_truncation', run_mode='"Parameter point"', sample_size='64',
                                max_trajectories='128', maximum_number_of_scatterings='1', snapshot_enabled='false')
        assert result.returncode == 2, result.stdout + result.stderr
        header = validate_moments(output / 'results_-2.000000_-32.000000/bincount.tsv')
        assert int(header['N_ever_captured']) > int(header['N_residence_samples'])
        assert int(header['N_injected']) == 128

        # Explicit debug runs stay outside the normal scientific result directory.
        result, output, config = run('local', ('--diagnostic',), run_mode='"Parameter point"',
                                     sample_size='4', snapshot_enabled='false')
        assert result.returncode == 0, result.stdout + result.stderr
        local = output / 'diagnostics/results_-2.000000_-32.000000'
        rows = list(csv.DictReader(io.StringIO((local / 'diagnostic_trajectory_summary.tsv').read_text()), delimiter='\t'))
        assert rows and any(row['rng_state_before_simulation'] for row in rows)
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

"""Verify output preflight and late publish failures exit collectively without hanging."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile

from bincount_contract import read_bincount


def main() -> None:
    """Inject filesystem failures without changing trajectory physics or the MPI queue."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--program', required=True)
    parser.add_argument('--mpiexec')
    parser.add_argument('--numproc-flag', default='-np')
    parser.add_argument('--ranks', default='1')
    args = parser.parse_args()
    template = (Path(__file__).resolve().parents[1] / 'examples/quickstart.cfg').read_text()
    launcher = [args.mpiexec, args.numproc_flag, args.ranks] if args.mpiexec else []
    with tempfile.TemporaryDirectory(prefix='damascus-output-failure-') as tmp:
        root = Path(tmp)

        def config_for(name: str, **changes: str) -> tuple[Path, Path]:
            """Create an isolated run; wall-clock budget settings are in seconds."""
            output = root / name
            settings = {'output_dir': json.dumps(str(output) + '/'), **changes}
            text = template
            for key, value in settings.items():
                pattern = rf'(?m)^{key}\s*=.*?;'
                text = re.sub(pattern, f'{key} = {value};', text) if re.search(pattern, text) else text + f'\n{key} = {value};\n'
            config = root / f'{name}.cfg'
            config.write_text(text)
            return config, output

        config, output = config_for('preflight')
        output.write_text('regular file blocking the output directory')
        result = subprocess.run(launcher + [args.program, str(config)], cwd=root, text=True, capture_output=True, timeout=20)
        assert result.returncode != 0 and 'failed to create output directory' in result.stderr
        assert 'Generate data...' not in result.stdout

        # Create a rename blocker AFTER preflight: wait for its flushed progress line.
        # The modest interpolation grid leaves ample time to inject the failure.
        config, output = config_for('late_publish', rate_radius_points='1000', rate_speed_points='256', rate_max_speed='0.02')
        process = subprocess.Popen(launcher + [args.program, str(config)], cwd=root, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        prefix = []
        try:
            for line in process.stdout:
                prefix.append(line)
                if 'Generate data...' in line:
                    break
            else:
                raise AssertionError('Run never reached post-preflight progress: ' + ''.join(prefix))
            target = output / 'results_-2.000000_-32.000000/bincount.tsv'
            target.mkdir()
            (target / 'keep').write_text('prevent atomic rename')
            remainder, _ = process.communicate(timeout=20)
            log = ''.join(prefix) + remainder
            assert process.returncode != 0 and 'cannot publish bincount.tsv' in log, log
            assert '[Finished in' not in log and not target.with_suffix('.tsv.tmp').exists()
            assert (target / 'keep').is_file() and (target.parent / 'snapshot').is_dir()
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()

        config, output = config_for('failed_transport', sample_size='4', max_trajectories='12', max_trajectory_wall_time_sec='1e-12')
        result = subprocess.run(launcher + [args.program, str(config)], cwd=root, text=True, capture_output=True, timeout=20)
        assert result.returncode == 2, result.stdout + result.stderr
        header, _, _ = read_bincount(output / 'results_-2.000000_-32.000000/bincount.tsv')
        assert header['target_reached'] == 'false' and int(header['computational_failures']) == 12
        assert int(header['N_injected']) == 12

        config, output = config_for('capture_untouched', run_mode='"Capture"', DM_cross_section_nucleon='1.0e-80')
        output.write_text('do not touch')

        def capture() -> dict:
            """Return the sole Capture JSON result; existing output paths stay untouched."""
            result = subprocess.run(launcher + [args.program, str(config)], cwd=root, text=True, capture_output=True, timeout=20)
            assert result.returncode == 0, result.stderr
            assert len(result.stdout.splitlines()) == 1
            return json.loads(result.stdout.split('=', 1)[1])

        rate = capture()['C_geom_s_inv']
        config.write_text(config.read_text().replace('DM_fraction = 1.0;', 'DM_fraction = 0.25;'))
        scaled = capture()
        assert math.isclose(scaled['C_geom_s_inv'] / rate, .25, rel_tol=1e-13)
        assert output.read_text() == 'do not touch'
        print('Preflight, late publication, failed transport and stdout-only Capture passed')


if __name__ == '__main__':
    main()

"""Verify that preflight and late write failures exit collectively, without hanging."""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--program", required=True)
    parser.add_argument("--mpiexec")
    parser.add_argument("--numproc-flag", default="-np")
    parser.add_argument("--ranks", default="1")
    args = parser.parse_args()
    template = (Path(__file__).resolve().parents[1] / "examples" / "quickstart.cfg").read_text()
    with tempfile.TemporaryDirectory(prefix="damascus-output-failure-") as tmp:
        root = Path(tmp)
        for case in ("preflight", "late_publish"):
            output = root / case
            if case == "preflight":
                output.write_text("This regular file cannot be an output directory.\n")
            else:
                # A nonempty directory prevents the final evaporation file rename.
                block = output / "results_-2.000000_-80.000000" / "evaporation_times.txt"
                block.mkdir(parents=True)
                (block / "keep").write_text("block publication\n")
            config = root / f"{case}.cfg"
            text = re.sub(r'output_dir = ".*?";', f'output_dir = "{output}/";', template)
            text = text.replace("DM_cross_section_nucleon = 1.0e-28;",
                                "DM_cross_section_nucleon = 1.0e-80;")
            config.write_text(text + "\nfixed_seed = 20260910;\n")
            command = [args.program, str(config)]
            if args.mpiexec:
                command = [args.mpiexec, args.numproc_flag, args.ranks] + command
            result = subprocess.run(command, cwd=root, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    timeout=20)
            assert result.returncode != 0, result.stdout
            assert "[Finished in" not in result.stdout, result.stdout
            expected = "failed to create output directory" if case == "preflight" else "failed to write evaporation_times.txt"
            assert expected in result.stdout, result.stdout
            if case == "preflight":
                assert "Generate data..." not in result.stdout, result.stdout
            else:
                assert (block.parent / "bincount.txt").is_file(), result.stdout
            print(f"{case}: exit={result.returncode}; no successful-finish message")

        # Unlimited incident budgets must still stop after the first failed
        # production history, including when several MPI claims are in flight.
        output=root/'failed_production'
        text=re.sub(r'output_dir = ".*?";',f'output_dir = "{output}/";',template)
        text=text.replace('sample_size = 1;','sample_size = 4;')
        text=text.replace('max_trajectories = 1;','max_trajectories = 0;')
        text=text.replace('max_trajectory_wall_time_sec = 5.0;','max_trajectory_wall_time_sec = 1e-12;')
        config=root/'failed_production.cfg'
        config.write_text(text+'\nproduction_mode = true;\nfixed_seed = 20260910;\n')
        command=[args.program,str(config)]
        if args.mpiexec:
            command=[args.mpiexec,args.numproc_flag,args.ranks]+command
        result=subprocess.run(command,cwd=root,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
        assert result.returncode==2,result.stdout
        metadata=json.loads((output/'results_-2.000000_-28.000000/metadata.json').read_text())
        assert metadata['production_accepted'] is False
        assert metadata['N_computational_failures']>=1
        assert metadata['N_inj']<=int(args.ranks)
        print('failed production: stopped without replacement loop; exit=2')

        # Capture never touches output_dir, including on a rerun.
        output=root/'capture_output_is_a_file'
        output.write_text('do not touch')
        text=re.sub(r'output_dir = ".*?";', f'output_dir = "{output}/";',template)
        text=text.replace('run_mode = "Parameter point";','run_mode = "Capture";')
        text=text.replace('DM_cross_section_nucleon = 1.0e-28;','DM_cross_section_nucleon = 1.0e-80;')
        text+='\nfixed_seed = 20260910;\n'
        config=root/'rerun.cfg'; config.write_text(text)
        def run(path: Path) -> dict:
            command=[args.program,str(path)]
            if args.mpiexec:
                command=[args.mpiexec,args.numproc_flag,args.ranks]+command
            result=subprocess.run(command,cwd=root,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
            assert result.returncode==0,result.stdout
            records=[line.split('=',1)[1] for line in result.stdout.splitlines() if line.startswith('CAPTURE_RESULT_JSON=')]
            assert len(records)==1,result.stdout
            return json.loads(records[0])
        rate=run(config)['C_geom_s_inv']
        config.write_text(text.replace('DM_fraction = 1.0;','DM_fraction = 0.25;'))
        scaled=run(config)
        assert math.isclose(scaled['C_geom_s_inv']/rate,.25,rel_tol=1e-13)
        assert scaled['physical_config']['DM_fraction']==.25
        assert output.read_text()=='do not touch'
        print('capture rerun: stdout only; geometric rate scales with DM_fraction')


if __name__ == "__main__":
    main()

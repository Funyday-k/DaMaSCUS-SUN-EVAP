"""Verify that preflight and late write failures exit collectively, without hanging."""
import argparse
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


if __name__ == "__main__":
    main()

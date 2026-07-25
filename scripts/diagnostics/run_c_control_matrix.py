#!/usr/bin/env python3
"""Run the C-point interpolation/seed/MPI control matrix reproducibly."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


REPLACEMENTS = {
    "ID": r'ID\s*=\s*"[^"]*";',
    "sample_size": r"sample_size\s*=\s*\d+;",
    "fixed_seed": r"fixed_seed\s*=\s*\d+;",
    "interpolation_points": r"interpolation_points\s*=\s*\d+;",
    "output_dir": r'output_dir\s*=\s*"[^"]*";',
    "trajectory_events_enabled": r"trajectory_events_enabled\s*=\s*(?:true|false);",
    "trajectory_trace_rate": r"trajectory_trace_rate\s*=\s*[0-9.eE+-]+;",
    "trajectory_trace_seed": r"trajectory_trace_seed\s*=\s*\d+;",
}


def comma_ints(value: str) -> list[int]:
    values = [int(item.strip()) for item in value.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected at least one integer")
    return values


def replace_setting(text: str, key: str, replacement: str) -> str:
    updated, count = re.subn(REPLACEMENTS[key], replacement, text, count=1)
    if count != 1:
        raise RuntimeError(f"could not replace exactly one {key!r} setting")
    return updated


def build_config(
    template: str,
    run_id: str,
    sample_size: int,
    seed: int,
    grid: int,
    output_dir: Path,
    trace_rate: float,
) -> str:
    text = template
    values = {
        "ID": f'ID = "{run_id}";',
        "sample_size": f"sample_size = {sample_size};",
        "fixed_seed": f"fixed_seed = {seed};",
        "interpolation_points": f"interpolation_points = {grid};",
        "output_dir": f'output_dir = "{output_dir}/";',
        "trajectory_events_enabled": "trajectory_events_enabled = " + ("true;" if trace_rate > 0 else "false;"),
        "trajectory_trace_rate": f"trajectory_trace_rate = {trace_rate:.8f};",
        "trajectory_trace_seed": f"trajectory_trace_seed = {seed};",
    }
    for key, value in values.items():
        text = replace_setting(text, key, value)
    return text


def output_complete(output_dir: Path) -> bool:
    result_dirs = list(output_dir.glob("results_*/run_metadata.json"))
    return len(result_dirs) == 1 and (result_dirs[0].parent / "evaporation_times.txt").exists()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path, default=Path("tests/config_evaporation_secondary_peak_C_10000.cfg"))
    parser.add_argument("--binary", type=Path, default=Path("build/src/DaMaSCUS-SUN"))
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--output-root", type=Path, default=Path("evaporation_c_controls"))
    parser.add_argument("--grids", type=comma_ints, default=comma_ints("0,20,1000,2000"))
    parser.add_argument("--seeds", type=comma_ints, default=comma_ints("20260722,20260723,20260724"))
    parser.add_argument("--ranks", type=comma_ints, default=comma_ints("4"))
    parser.add_argument("--sample-size", type=int, default=10_000)
    parser.add_argument("--trace-rate", type=float, default=0.0)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    repo = Path.cwd().resolve()
    template_path = (repo / args.template).resolve()
    binary = (repo / args.binary).resolve()
    output_root = (repo / args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    config_dir = output_root / "configs"
    log_dir = output_root / "logs"
    config_dir.mkdir(exist_ok=True)
    log_dir.mkdir(exist_ok=True)
    template = template_path.read_text(encoding="utf-8")

    manifest = {
        "schema_version": "c-control-matrix-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "template": str(template_path),
        "binary": str(binary),
        "grids": args.grids,
        "seeds": args.seeds,
        "ranks": args.ranks,
        "sample_size": args.sample_size,
        "trace_rate": args.trace_rate,
        "runs": [],
    }

    for ranks in args.ranks:
        for grid in args.grids:
            for seed in args.seeds:
                run_id = f"c_grid_{grid:04d}_seed_{seed}_r{ranks}"
                output_dir = output_root / run_id
                config_path = config_dir / f"{run_id}.cfg"
                log_path = log_dir / f"{run_id}.log"
                config_text = build_config(
                    template, run_id, args.sample_size, seed, grid, output_dir, args.trace_rate
                )
                config_path.write_text(config_text, encoding="utf-8")
                command = [args.mpiexec, "--oversubscribe", "-n", str(ranks), str(binary), str(config_path)]
                entry = {
                    "run_id": run_id,
                    "grid": grid,
                    "seed": seed,
                    "mpi_size": ranks,
                    "config": str(config_path),
                    "output_dir": str(output_dir),
                    "log": str(log_path),
                    "command": command,
                }
                manifest["runs"].append(entry)

                if output_complete(output_dir) and not args.force:
                    entry["status"] = "skipped_complete"
                    print(f"SKIP {run_id}: complete output exists", flush=True)
                    continue
                print("RUN  " + " ".join(command), flush=True)
                if args.dry_run:
                    entry["status"] = "dry_run"
                    continue
                with log_path.open("w", encoding="utf-8") as log:
                    completed = subprocess.run(command, cwd=repo, stdout=log, stderr=subprocess.STDOUT)
                entry["returncode"] = completed.returncode
                entry["status"] = "complete" if completed.returncode == 0 and output_complete(output_dir) else "failed"
                print(f"DONE {run_id}: {entry['status']} (returncode={completed.returncode})", flush=True)
                (output_root / "matrix_manifest.json").write_text(
                    json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
                )
                if entry["status"] == "failed":
                    return completed.returncode or 1

    (output_root / "matrix_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

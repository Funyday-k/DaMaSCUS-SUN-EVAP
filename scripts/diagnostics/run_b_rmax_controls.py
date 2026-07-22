#!/usr/bin/env python3
"""Run B-point escape-radius controls from one frozen executable."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from run_c_control_matrix import build_config, output_complete


def comma_floats(value: str) -> list[float]:
    values = [float(item.strip()) for item in value.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected at least one radius")
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path, default=Path("tests/config_evaporation_secondary_peak_B_10000.cfg"))
    parser.add_argument("--binary", type=Path, default=Path("build/src/DaMaSCUS-SUN"))
    parser.add_argument("--output-root", type=Path, default=Path("evaporation_b_rmax_controls"))
    parser.add_argument("--radii", type=comma_floats, default=comma_floats("2,3,5"))
    parser.add_argument("--interpolation-points", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260722)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--sample-size", type=int, default=10000)
    parser.add_argument("--trace-rate", type=float, default=0.10)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    repo = Path.cwd().resolve()
    output_root = (repo / args.output_root).resolve()
    config_dir = output_root / "configs"
    log_dir = output_root / "logs"
    config_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    template = (repo / args.template).resolve().read_text(encoding="utf-8")
    binary = (repo / args.binary).resolve()
    runs: list[dict[str, object]] = []

    for radius in args.radii:
        radius_tag = str(radius).replace(".", "p")
        run_id = f"b_rmax_{radius_tag}_grid_{args.interpolation_points}_seed_{args.seed}_r{args.ranks}"
        output_dir = output_root / run_id
        config_path = config_dir / f"{run_id}.cfg"
        log_path = log_dir / f"{run_id}.log"
        config = build_config(
            template,
            run_id,
            args.sample_size,
            args.seed,
            args.interpolation_points,
            output_dir,
            args.trace_rate,
        )
        config, count = re.subn(
            r"R_escape_Rsun\s*=\s*[0-9.eE+-]+;",
            f"R_escape_Rsun = {radius:.17g}" + (".0;" if radius.is_integer() else ";"),
            config,
            count=1,
        )
        if count != 1:
            raise RuntimeError("template must contain exactly one R_escape_Rsun setting")
        config_path.write_text(config, encoding="utf-8")
        command = ["mpiexec", "--oversubscribe", "-n", str(args.ranks), str(binary), str(config_path)]
        entry: dict[str, object] = {
            "run_id": run_id,
            "R_escape_Rsun": radius,
            "interpolation_points": args.interpolation_points,
            "seed": args.seed,
            "mpi_size": args.ranks,
            "sample_size": args.sample_size,
            "trace_rate": args.trace_rate,
            "config": str(config_path),
            "output_dir": str(output_dir),
            "log": str(log_path),
            "command": command,
        }
        runs.append(entry)
        if output_complete(output_dir) and not args.force:
            entry["status"] = "skipped_complete"
            print(f"SKIP {run_id}: complete output exists", flush=True)
            continue
        print("RUN  " + " ".join(command), flush=True)
        with log_path.open("w", encoding="utf-8") as log:
            completed = subprocess.run(command, cwd=repo, stdout=log, stderr=subprocess.STDOUT)
        entry["returncode"] = completed.returncode
        entry["status"] = "complete" if completed.returncode == 0 and output_complete(output_dir) else "failed"
        print(f"DONE {run_id}: {entry['status']}", flush=True)
        if entry["status"] == "failed":
            break

    manifest = {
        "schema_version": "b-rmax-control-matrix-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "binary": str(binary),
        "runs": runs,
    }
    (output_root / "matrix_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return 0 if all(entry.get("status") in {"complete", "skipped_complete"} for entry in runs) else 1


if __name__ == "__main__":
    raise SystemExit(main())

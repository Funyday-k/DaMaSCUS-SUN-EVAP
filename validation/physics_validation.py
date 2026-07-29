#!/usr/bin/env python3
"""Run and evaluate a production-scale numerical-physics convergence matrix.

The fast analytic/invariant checks live in tests/test_Physics_Validation.cpp.
This script targets the slower question: do final observables remain compatible
when the scattering-rate interpolation grid and random seed are changed?
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import statistics
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


HEADER_PATTERN = re.compile(r"^#\s*([^=]+?)\s*=\s*(.*?)\s*$")


@dataclass
class RunMetrics:
    grid: int
    seed: int
    result_dir: str
    valid_trajectories: int
    captured_particles: int
    capture_rate: float
    capture_error: float
    capture_ci_lower: float
    capture_ci_upper: float
    numerical_failure_rate: float
    average_scatterings: float
    early_stop: Optional[str]
    evaporation_events: int
    complete_evaporation_fraction: Optional[float]
    evaporation_median_sec: Optional[float]
    evaporation_p90_sec: Optional[float]
    captured_radial_distribution: List[float]


def _number(headers: Dict[str, str], key: str, converter=float):
    if key not in headers:
        raise ValueError(f"Missing required output header: {key}")
    try:
        return converter(headers[key])
    except ValueError as error:
        raise ValueError(f"Invalid value for output header {key}: {headers[key]}") from error


def _quantile(values: Sequence[float], probability: float) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    coordinate = probability * (len(ordered) - 1)
    lower = int(math.floor(coordinate))
    upper = int(math.ceil(coordinate))
    fraction = coordinate - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def _parse_bincount(path: Path):
    headers: Dict[str, str] = {}
    early_stop = None
    captured_dt: List[float] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("# EARLY_STOP:"):
                    early_stop = line.split(":", 1)[1].strip()
                match = HEADER_PATTERN.match(line)
                if match:
                    headers[match.group(1).strip()] = match.group(2).strip()
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"Malformed bincount row in {path}: {line}")
            # v2 adds explicit lower/upper radial edges before cap_dt.
            captured_dt_index = 3 if len(fields) >= 11 else 1
            captured_dt.append(float(fields[captured_dt_index]))

    total = sum(captured_dt)
    radial_distribution = (
        [value / total for value in captured_dt] if total > 0.0 else []
    )
    return headers, early_stop, radial_distribution


def _parse_evaporation_times(path: Path) -> List[float]:
    lifetimes = []
    if not path.exists():
        return lifetimes
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"Malformed evaporation row in {path}: {line}")
            lifetime = float(fields[2])
            if math.isfinite(lifetime) and lifetime > 0.0:
                lifetimes.append(lifetime)
    return lifetimes


def load_run_metrics(result_dir: Path, grid: int, seed: int) -> RunMetrics:
    bincount_path = result_dir / "bincount.txt"
    headers, early_stop, radial_distribution = _parse_bincount(bincount_path)
    lifetimes = _parse_evaporation_times(result_dir / "evaporation_times.txt")
    captured_particles = _number(headers, "captured_particles", int)
    return RunMetrics(
        grid=grid,
        seed=seed,
        result_dir=str(result_dir),
        valid_trajectories=_number(headers, "valid_trajectories", int),
        captured_particles=captured_particles,
        capture_rate=_number(headers, "capture_rate_valid"),
        capture_error=_number(headers, "capture_rate_valid_err"),
        capture_ci_lower=_number(headers, "capture_rate_valid_CI_95_lower"),
        capture_ci_upper=_number(headers, "capture_rate_valid_CI_95_upper"),
        numerical_failure_rate=_number(headers, "numerical_failure_rate"),
        average_scatterings=_number(headers, "average_scatterings"),
        early_stop=early_stop,
        evaporation_events=len(lifetimes),
        complete_evaporation_fraction=(
            len(lifetimes) / captured_particles if captured_particles > 0 else None
        ),
        evaporation_median_sec=statistics.median(lifetimes) if lifetimes else None,
        evaporation_p90_sec=_quantile(lifetimes, 0.9),
        captured_radial_distribution=radial_distribution,
    )


def _check(name: str, passed: bool, value, limit, required: bool = True):
    return {
        "name": name,
        "status": "pass" if passed else ("fail" if required else "not_evaluated"),
        "value": value,
        "limit": limit,
    }


def quality_checks(metrics: RunMetrics, maximum_failure_rate: float):
    return [
        _check("no_early_stop", metrics.early_stop is None, metrics.early_stop, None),
        _check(
            "numerical_failure_rate",
            metrics.numerical_failure_rate <= maximum_failure_rate,
            metrics.numerical_failure_rate,
            maximum_failure_rate,
        ),
        _check(
            "has_classified_trajectories",
            metrics.valid_trajectories > 0,
            metrics.valid_trajectories,
            "> 0",
        ),
    ]


def _relative_difference(candidate: float, baseline: float, floor: float = 0.0):
    return abs(candidate - baseline) / max(abs(baseline), floor, 1.0e-300)


def _total_variation(candidate: Sequence[float], baseline: Sequence[float]):
    if not candidate or not baseline or len(candidate) != len(baseline):
        return None
    return 0.5 * sum(abs(left - right) for left, right in zip(candidate, baseline))


def _two_proportion_z(successes_a: int, trials_a: int, successes_b: int, trials_b: int):
    if trials_a <= 0 or trials_b <= 0:
        return None
    pooled = (successes_a + successes_b) / (trials_a + trials_b)
    variance = pooled * (1.0 - pooled) * (1.0 / trials_a + 1.0 / trials_b)
    difference = abs(successes_a / trials_a - successes_b / trials_b)
    if variance <= 0.0:
        return 0.0 if difference == 0.0 else math.inf
    return difference / math.sqrt(variance)


def compare_runs(
    baseline: RunMetrics,
    candidate: RunMetrics,
    maximum_scattering_relative_difference: float = 0.10,
    maximum_radial_total_variation: float = 0.10,
    maximum_evaporation_median_relative_difference: float = 0.20,
    maximum_evaporation_p90_relative_difference: float = 0.30,
    maximum_complete_evaporation_fraction_z: float = 3.0,
    minimum_evaporation_events: int = 30,
    require_evaporation_events: bool = False,
):
    checks = []
    intervals_overlap = max(baseline.capture_ci_lower, candidate.capture_ci_lower) <= min(
        baseline.capture_ci_upper, candidate.capture_ci_upper
    )
    checks.append(
        _check(
            "capture_95pct_intervals_overlap",
            intervals_overlap,
            [
                baseline.capture_ci_lower,
                baseline.capture_ci_upper,
                candidate.capture_ci_lower,
                candidate.capture_ci_upper,
            ],
            "overlap",
        )
    )

    evaporation_fraction_z = _two_proportion_z(
        baseline.evaporation_events,
        baseline.captured_particles,
        candidate.evaporation_events,
        candidate.captured_particles,
    )
    fraction_available = evaporation_fraction_z is not None
    checks.append(
        _check(
            "complete_evaporation_fraction_z",
            fraction_available
            and evaporation_fraction_z <= maximum_complete_evaporation_fraction_z,
            evaporation_fraction_z,
            maximum_complete_evaporation_fraction_z,
            required=fraction_available,
        )
    )

    scattering_difference = _relative_difference(
        candidate.average_scatterings, baseline.average_scatterings, floor=1.0
    )
    checks.append(
        _check(
            "average_scatterings_relative_difference",
            scattering_difference <= maximum_scattering_relative_difference,
            scattering_difference,
            maximum_scattering_relative_difference,
        )
    )
    radial_distance = _total_variation(
        candidate.captured_radial_distribution,
        baseline.captured_radial_distribution,
    )
    radial_available = radial_distance is not None
    checks.append(
        _check(
            "captured_radial_total_variation",
            radial_available and radial_distance <= maximum_radial_total_variation,
            radial_distance,
            maximum_radial_total_variation,
            required=radial_available,
        )
    )

    enough_evaporation_events = (
        baseline.evaporation_events >= minimum_evaporation_events
        and candidate.evaporation_events >= minimum_evaporation_events
        and baseline.evaporation_median_sec is not None
        and candidate.evaporation_median_sec is not None
    )
    evaporation_difference = None
    if enough_evaporation_events:
        evaporation_difference = _relative_difference(
            candidate.evaporation_median_sec, baseline.evaporation_median_sec
        )
    checks.append(
        _check(
            "evaporation_median_relative_difference",
            enough_evaporation_events
            and evaporation_difference <= maximum_evaporation_median_relative_difference,
            evaporation_difference,
            maximum_evaporation_median_relative_difference,
            required=require_evaporation_events or enough_evaporation_events,
        )
    )

    evaporation_p90_difference = None
    if enough_evaporation_events:
        evaporation_p90_difference = _relative_difference(
            candidate.evaporation_p90_sec, baseline.evaporation_p90_sec
        )
    checks.append(
        _check(
            "evaporation_p90_relative_difference",
            enough_evaporation_events
            and evaporation_p90_difference <= maximum_evaporation_p90_relative_difference,
            evaporation_p90_difference,
            maximum_evaporation_p90_relative_difference,
            required=require_evaporation_events or enough_evaporation_events,
        )
    )

    return {
        "seed": candidate.seed,
        "baseline_grid": baseline.grid,
        "candidate_grid": candidate.grid,
        "status": "fail" if any(item["status"] == "fail" for item in checks) else "pass",
        "checks": checks,
    }


def _format_config_value(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return json.dumps(value)
    return str(value)


def rewrite_config(source: str, overrides: Dict[str, object]) -> str:
    result = source
    for key, value in overrides.items():
        rendered = _format_config_value(value)
        pattern = re.compile(
            r"(?m)^(\s*" + re.escape(key) + r"\s*=\s*)([^;]*)(;)"
        )
        result, replacements = pattern.subn(
            lambda match: match.group(1) + rendered + match.group(3), result, count=1
        )
        if replacements == 0:
            if result and not result.endswith("\n"):
                result += "\n"
            result += f"{key} = {rendered};\n"
    return result


def _find_result_directory(output_root: Path) -> Path:
    matches = sorted(path.parent for path in output_root.rglob("bincount.txt"))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one bincount.txt under {output_root}, found {len(matches)}"
        )
    return matches[0]


def _run_one(args, source_config: str, grid: int, seed: int):
    run_dir = args.output / f"grid-{grid}_seed-{seed}"
    if run_dir.exists() and args.overwrite:
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    output_root = run_dir / "outputs"
    config_path = run_dir / "config.cfg"
    log_path = run_dir / "run.log"

    overrides: Dict[str, object] = {
        "run_mode": "Parameter point",
        "interpolation_points": grid,
        "fixed_seed": seed,
        "snapshot_enabled": False,
        "output_dir": str(output_root) + "/",
    }
    if args.sample_size is not None:
        overrides["sample_size"] = args.sample_size
    if args.max_trajectories is not None:
        overrides["max_trajectories"] = args.max_trajectories
    config_path.write_text(rewrite_config(source_config, overrides), encoding="utf-8")

    if args.reuse_existing and list(output_root.rglob("bincount.txt")):
        return load_run_metrics(_find_result_directory(output_root), grid, seed)
    if args.dry_run:
        return None
    if output_root.exists() and not args.overwrite:
        raise RuntimeError(
            f"Output already exists for grid={grid}, seed={seed}; use --reuse-existing or --overwrite"
        )

    command = [str(args.executable), str(config_path)]
    if args.ranks > 1:
        command = [args.mpiexec, "-n", str(args.ranks)] + command
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            command,
            cwd=str(run_dir),
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=args.timeout,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Simulation failed for grid={grid}, seed={seed}; see {log_path}"
        )
    return load_run_metrics(_find_result_directory(output_root), grid, seed)


def _all_pass(check_groups: Iterable[Iterable[dict]]) -> bool:
    return all(
        item["status"] != "fail" for group in check_groups for item in group
    )


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--grids", type=int, nargs="+", default=[0, 1000, 2000])
    parser.add_argument("--seeds", type=int, nargs="+", default=[271828, 314159, 161803])
    parser.add_argument("--sample-size", type=int)
    parser.add_argument("--max-trajectories", type=int)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--timeout", type=float, default=86400.0)
    parser.add_argument("--maximum-failure-rate", type=float, default=1.0e-3)
    parser.add_argument("--maximum-scattering-relative-difference", type=float, default=0.10)
    parser.add_argument("--maximum-radial-total-variation", type=float, default=0.10)
    parser.add_argument("--maximum-evaporation-median-relative-difference", type=float, default=0.20)
    parser.add_argument("--maximum-evaporation-p90-relative-difference", type=float, default=0.30)
    parser.add_argument("--maximum-complete-evaporation-fraction-z", type=float, default=3.0)
    parser.add_argument("--minimum-evaporation-events", type=int, default=30)
    parser.add_argument("--require-evaporation-events", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--reuse-existing", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.executable = args.executable.resolve()
    args.config = args.config.resolve()
    args.output = args.output.resolve()
    if 0 not in args.grids:
        raise SystemExit("--grids must include 0 as the direct-rate baseline")
    if args.ranks < 1:
        raise SystemExit("--ranks must be positive")
    if args.overwrite and args.reuse_existing:
        raise SystemExit("--overwrite and --reuse-existing are mutually exclusive")
    if not args.executable.is_file():
        raise SystemExit(f"Executable not found: {args.executable}")
    if not args.config.is_file():
        raise SystemExit(f"Configuration not found: {args.config}")

    args.output.mkdir(parents=True, exist_ok=True)
    source_config = args.config.read_text(encoding="utf-8")
    metrics = []
    for seed in args.seeds:
        for grid in args.grids:
            print(f"[physics-validation] grid={grid} seed={seed}", flush=True)
            run_metrics = _run_one(args, source_config, grid, seed)
            if run_metrics is not None:
                metrics.append(run_metrics)

    if args.dry_run:
        print(f"Prepared {len(args.grids) * len(args.seeds)} configurations under {args.output}")
        return 0

    quality = []
    for item in metrics:
        quality.append(
            {
                "grid": item.grid,
                "seed": item.seed,
                "checks": quality_checks(item, args.maximum_failure_rate),
            }
        )

    indexed = {(item.grid, item.seed): item for item in metrics}
    comparisons = []
    for seed in args.seeds:
        baseline = indexed[(0, seed)]
        for grid in args.grids:
            if grid == 0:
                continue
            comparisons.append(
                compare_runs(
                    baseline,
                    indexed[(grid, seed)],
                    args.maximum_scattering_relative_difference,
                    args.maximum_radial_total_variation,
                    args.maximum_evaporation_median_relative_difference,
                    args.maximum_evaporation_p90_relative_difference,
                    args.maximum_complete_evaporation_fraction_z,
                    args.minimum_evaporation_events,
                    args.require_evaporation_events,
                )
            )

    quality_pass = _all_pass(item["checks"] for item in quality)
    comparisons_pass = all(item["status"] == "pass" for item in comparisons)
    report = {
        "status": "pass" if quality_pass and comparisons_pass else "fail",
        "baseline_grid": 0,
        "metrics": [asdict(item) for item in metrics],
        "quality": quality,
        "comparisons": comparisons,
    }
    report_path = args.output / "physics_validation_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"Physics validation: {report['status'].upper()}")
    print(f"Report: {report_path}")
    for comparison in comparisons:
        print(
            f"  seed={comparison['seed']} grid={comparison['candidate_grid']}: "
            f"{comparison['status']}"
        )
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    sys.exit(main())

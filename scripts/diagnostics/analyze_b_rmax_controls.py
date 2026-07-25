#!/usr/bin/env python3
"""Analyze B-point Rmax controls and extract stratified replay samples."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp

from analyze_evaporation_reproduction import (
    detect_logtime_peaks,
    ecdf,
    histogram_kde,
    mode_fractions,
    parse_header_metrics,
    read_json,
    write_csv,
)


RUN_PATTERN = re.compile(r"^b_rmax_([0-9p]+)_grid_(\d+)_seed_(\d+)_r(\d+)$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def vector(values: np.ndarray | list[float]) -> str:
    return ";".join(f"{float(value):.9g}" for value in values)


def key(row: dict[str, str]) -> tuple[int, int]:
    return int(row["rank"]), int(row["trajectory_id"])


def discover(root: Path) -> list[dict[str, Any]]:
    runs: list[dict[str, Any]] = []
    for directory in sorted(root.iterdir()):
        match = RUN_PATTERN.match(directory.name)
        if not match:
            continue
        radius = float(match.group(1).replace("p", "."))
        result_dirs = list(directory.glob("results_*/trajectory_summary.tsv"))
        if len(result_dirs) != 1:
            continue
        result_dir = result_dirs[0].parent
        fields, rows = read_tsv(result_dirs[0])
        completed = [row for row in rows if row["event_observed"] == "1"]
        final_log = np.asarray([math.log10(float(row["lifetime_final_unbinding_s"])) for row in completed])
        validated_log = np.asarray([math.log10(float(row["lifetime_validated_escape_s"])) for row in completed])
        runs.append(
            {
                "radius": radius,
                "grid": int(match.group(2)),
                "seed": int(match.group(3)),
                "mpi_size": int(match.group(4)),
                "directory": directory,
                "result_dir": result_dir,
                "summary_fields": fields,
                "rows": rows,
                "completed": completed,
                "final_log": final_log,
                "validated_log": validated_log,
                "metadata": read_json(result_dir / "run_metadata.json"),
                "bincount": parse_header_metrics(result_dir / "bincount.txt"),
            }
        )
    return sorted(runs, key=lambda run: run["radius"])


def plot(output: Path, runs: list[dict[str, Any]]) -> None:
    colors = {2.0: "#2563EB", 3.0: "#D97706", 5.0: "#059669"}
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.0), dpi=170)
    final_ecdf, final_kde, validated_ecdf, tail_axis = axes.flat
    all_final = [run["final_log"] for run in runs]
    all_validated = [run["validated_log"] for run in runs]
    lo = min(float(np.quantile(values, 0.001)) for values in all_final) - 0.15
    hi = max(float(np.quantile(values, 0.999)) for values in all_validated) + 0.15
    density_grid = np.linspace(lo, hi, 520)
    for run in runs:
        radius = run["radius"]
        x, y = ecdf(run["final_log"])
        final_ecdf.plot(x, y, color=colors[radius], linewidth=1.8, label=f"Rmax={radius:g}")
        final_kde.plot(density_grid, histogram_kde(run["final_log"], density_grid), color=colors[radius], linewidth=1.8, label=f"Rmax={radius:g}")
        x, y = ecdf(run["validated_log"])
        validated_ecdf.plot(x, y, color=colors[radius], linewidth=1.8, label=f"Rmax={radius:g}")
        tails = 10 ** run["validated_log"] - 10 ** run["final_log"]
        positive = tails[tails > 0]
        tail_axis.hist(np.log10(positive), bins=70, density=True, histtype="step", color=colors[radius], linewidth=1.5, label=f"Rmax={radius:g}")
    final_ecdf.set(xlim=(lo, hi), ylim=(0, 1), xlabel="log10(final-unbinding lifetime / s)", ylabel="Cumulative fraction")
    final_ecdf.set_title("Final-unbinding ECDF", loc="left", fontweight="semibold")
    final_kde.set(xlim=(lo, hi), xlabel="log10(final-unbinding lifetime / s)", ylabel="Probability density")
    final_kde.set_title("Final-unbinding KDE", loc="left", fontweight="semibold")
    validated_ecdf.set(xlim=(lo, hi), ylim=(0, 1), xlabel="log10(validated-escape lifetime / s)", ylabel="Cumulative fraction")
    validated_ecdf.set_title("Validated-escape ECDF", loc="left", fontweight="semibold")
    tail_axis.set(xlabel="log10(Rmax travel time / s)", ylabel="Probability density")
    tail_axis.set_title("Time from final unbinding to Rmax", loc="left", fontweight="semibold")
    for axis in axes.flat:
        axis.legend(frameon=False, fontsize=8)
        axis.grid(axis="y", color="#CBD5E1", linewidth=0.65, alpha=0.7)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    fig.suptitle("B-point escape-radius controls", x=0.07, ha="left", fontsize=16, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output)
    plt.close(fig)


def extract_strata(run: dict[str, Any], output_dir: Path, per_stratum: int) -> list[dict[str, Any]]:
    peaks, boundaries = detect_logtime_peaks(run["final_log"])
    if len(peaks) < 2 or not boundaries:
        raise RuntimeError("B run is not bimodal under the configured peak rule")
    strata = [
        ("short_peak", float(peaks[0]), 0.25),
        ("valley", float(boundaries[0]), 0.18),
        ("long_peak", float(peaks[1]), 0.25),
    ]
    completed = run["completed"]
    traced = [
        row
        for row in completed
        if row["trace_written"] == "1"
        and row["rng_state_before_initial_conditions"]
        and row["rng_state_before_simulation"]
    ]
    selected: list[dict[str, str]] = []
    summary: list[dict[str, Any]] = []
    for name, center, half_width in strata:
        all_candidates = [row for row in completed if abs(math.log10(float(row["lifetime_final_unbinding_s"])) - center) <= half_width]
        trace_candidates = [row for row in traced if abs(math.log10(float(row["lifetime_final_unbinding_s"])) - center) <= half_width]
        trace_candidates.sort(key=lambda row: (abs(math.log10(float(row["lifetime_final_unbinding_s"])) - center), key(row)))
        current = trace_candidates[:per_stratum]
        for row in current:
            row["replay_stratum"] = name
        selected.extend(current)
        outside_fraction = np.asarray([
            float(row["time_outside_sun_s"]) / (float(row["time_inside_sun_s"]) + float(row["time_outside_sun_s"]))
            for row in current
        ])
        summary.append(
            {
                "stratum": name,
                "center_log10_s": center,
                "half_width_dex": half_width,
                "all_completed_candidates": len(all_candidates),
                "traced_candidates": len(trace_candidates),
                "selected_replays": len(current),
                "median_log10_lifetime_s": float(np.median([math.log10(float(row["lifetime_final_unbinding_s"])) for row in current])),
                "median_capture_radius_Rsun": float(np.median([float(row["r_capture_Rsun"]) for row in current])),
                "median_capture_energy_eV": float(np.median([float(row["E_capture_eV"]) for row in current])),
                "median_outside_fraction": float(np.median(outside_fraction)),
                "recapture_fraction": float(np.mean([int(row["n_recapture"]) > 0 for row in current])),
                "median_scatter_count": float(np.median([int(row["n_scatter_total"]) for row in current])),
            }
        )
    fields = ["replay_stratum"] + run["summary_fields"]
    write_tsv(output_dir / "stratified_replays.tsv", fields, selected)
    selected_keys = {key(row) for row in selected}
    event_fields, event_rows = read_tsv(run["result_dir"] / "trajectory_events.tsv")
    selected_events = [row for row in event_rows if key(row) in selected_keys]
    write_tsv(output_dir / "stratified_events.tsv", ["replay_stratum"] + event_fields, [
        {"replay_stratum": next(row["replay_stratum"] for row in selected if key(row) == key(event)), **event}
        for event in selected_events
    ])
    write_csv(output_dir / "strata_summary.csv", summary)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replays-per-stratum", type=int, default=100)
    args = parser.parse_args()
    runs = discover(args.control_root.resolve())
    if [run["radius"] for run in runs] != [2.0, 3.0, 5.0]:
        raise RuntimeError("expected exactly Rmax=2,3,5 completed controls")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    reference = runs[0]["final_log"]
    metric_rows: list[dict[str, Any]] = []
    quantile_rows: list[dict[str, Any]] = []
    for run in runs:
        peaks, boundaries = detect_logtime_peaks(run["final_log"])
        fractions = mode_fractions(run["final_log"], boundaries)
        ks = ks_2samp(run["final_log"], reference)
        tail = 10 ** run["validated_log"] - 10 ** run["final_log"]
        completed = run["completed"]
        bincount = run["bincount"]
        total = float(bincount["total_trajectories"])
        failures = float(bincount["numerical_failures"])
        boolean_invariants = [value for key_name, value in run["metadata"].items() if key_name.endswith("_invariant") or key_name.endswith("_invariants") or key_name == "legacy_evaporation_reconciliation"]
        metric_rows.append(
            {
                "R_escape_Rsun": run["radius"],
                "n_complete": len(completed),
                "captured_invalid": len(run["rows"]) - len(completed),
                "peak_log10_s": vector(peaks),
                "mode_fractions": vector(fractions),
                "mode_boundary_log10_s": vector(boundaries),
                "ks_d_vs_R2_final_unbinding": float(ks.statistic),
                "ks_p_vs_R2_final_unbinding": float(ks.pvalue),
                "median_Rmax_travel_time_s": float(np.median(tail)),
                "p90_Rmax_travel_time_s": float(np.quantile(tail, 0.9)),
                "mean_scatterings_completed": float(np.mean([int(row["n_scatter_total"]) for row in completed])),
                "numerical_failure_rate": failures / total,
                "traced_completed": sum(row["trace_written"] == "1" for row in completed),
                "all_metadata_invariants": all(boolean_invariants),
                "git_commit": run["metadata"].get("git_commit"),
                "git_dirty": run["metadata"].get("git_dirty"),
                "interpolation_points": run["metadata"].get("interpolation_points"),
            }
        )
        for quantile in (0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99):
            quantile_rows.append(
                {
                    "R_escape_Rsun": run["radius"],
                    "quantile": quantile,
                    "final_unbinding_log10_s": float(np.quantile(run["final_log"], quantile)),
                    "validated_escape_log10_s": float(np.quantile(run["validated_log"], quantile)),
                    "validated_minus_final_dex": float(np.quantile(run["validated_log"], quantile) - np.quantile(run["final_log"], quantile)),
                }
            )
    write_csv(output_dir / "rmax_metrics.csv", metric_rows)
    write_csv(output_dir / "rmax_quantiles.csv", quantile_rows)
    strata_summary = extract_strata(runs[0], output_dir, args.replays_per_stratum)
    plot(output_dir / "b_rmax_controls.png", runs)
    manifest = {
        "schema_version": "b-rmax-control-analysis-v1",
        "control_root": str(args.control_root.resolve()),
        "inputs": [
            {
                "R_escape_Rsun": run["radius"],
                "summary": str(run["result_dir"] / "trajectory_summary.tsv"),
                "summary_sha256": sha256(run["result_dir"] / "trajectory_summary.tsv"),
                "events": str(run["result_dir"] / "trajectory_events.tsv"),
                "events_sha256": sha256(run["result_dir"] / "trajectory_events.tsv"),
                "metadata": str(run["result_dir"] / "run_metadata.json"),
                "metadata_sha256": sha256(run["result_dir"] / "run_metadata.json"),
            }
            for run in runs
        ],
        "outputs": ["rmax_metrics.csv", "rmax_quantiles.csv", "strata_summary.csv", "stratified_replays.tsv", "stratified_events.tsv", "b_rmax_controls.png"],
        "strata": strata_summary,
    }
    with (output_dir / "analysis_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Saved B Rmax analysis to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

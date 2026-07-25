#!/usr/bin/env python3
"""Summarize C-point interpolation, seed, and MPI controls."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp

from analyze_evaporation_reproduction import (
    Comparison,
    QUANTILES,
    analyze,
    detect_logtime_peaks,
    ecdf,
    histogram_kde,
    load_v3,
    mode_fractions,
    resolve_lifetime_file,
    sha256,
    write_csv,
)


RUN_PATTERN = re.compile(r"^c_grid_(\d+)_seed_(\d+)_r(\d+)$")


def discover(root: Path) -> list[dict[str, Any]]:
    runs: list[dict[str, Any]] = []
    if not root.exists():
        return runs
    for directory in sorted(root.iterdir()):
        match = RUN_PATTERN.match(directory.name)
        if not match or not directory.is_dir():
            continue
        grid, seed, ranks = (int(value) for value in match.groups())
        lifetime_file = resolve_lifetime_file(directory)
        runs.append(
            {
                "grid": grid,
                "seed": seed,
                "mpi_size": ranks,
                "directory": directory.resolve(),
                "lifetime_file": lifetime_file,
            }
        )
    return runs


def mean(values: list[float]) -> float:
    return float(np.mean(values)) if values else float("nan")


def range_text(values: list[float]) -> str:
    return f"{min(values):.9g};{max(values):.9g}" if values else ""


def parse_vector(value: str) -> list[float]:
    return [float(item) for item in value.split(";") if item]


def plot_controls(
    output: Path,
    runs: list[dict[str, Any]],
    reference_file: Path,
    metric_rows: list[dict[str, Any]],
) -> None:
    reference_log = np.log10(load_v3(reference_file)[2])
    rank4 = [run for run in runs if run["mpi_size"] == 4]
    selected_seed = min(run["seed"] for run in rank4)
    selected = {run["grid"]: run for run in rank4 if run["seed"] == selected_seed}
    colors = {0: "#64748B", 20: "#C2410C", 1000: "#2563EB", 2000: "#059669"}
    grid_positions = {0: 0, 20: 1, 1000: 2, 2000: 3}

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.2), dpi=170)
    ecdf_axis, kde_axis, ks_axis, weight_axis = axes.flat
    reference_x, reference_y = ecdf(reference_log)
    ecdf_axis.plot(reference_x, reference_y, color="#111827", linewidth=2.1, label="Historical (30 ranks)")
    all_selected_logs = [reference_log]
    for grid in sorted(selected):
        log_values = np.log10(load_v3(selected[grid]["lifetime_file"])[2])
        all_selected_logs.append(log_values)
        x, y = ecdf(log_values)
        ecdf_axis.plot(x, y, color=colors[grid], linewidth=1.5, label=f"grid={grid}")
    grid20_log = np.log10(load_v3(selected[20]["lifetime_file"])[2])
    grid20_ks = ks_2samp(grid20_log, reference_log)
    ecdf_axis.axvline(float(grid20_ks.statistic_location), color=colors[20], linestyle="--", linewidth=1.1)
    ecdf_axis.text(
        0.03,
        0.94,
        f"grid=20 max ΔCDF: D={grid20_ks.statistic:.4f} at {grid20_ks.statistic_location:.3f}",
        transform=ecdf_axis.transAxes,
        va="top",
        color="#9A3412",
        fontsize=9,
    )
    lo = min(float(np.quantile(values, 0.001)) for values in all_selected_logs) - 0.15
    hi = max(float(np.quantile(values, 0.999)) for values in all_selected_logs) + 0.15
    ecdf_axis.set(xlim=(lo, hi), ylim=(0, 1), xlabel="log10(lifetime / s)", ylabel="Cumulative fraction")
    ecdf_axis.set_title(f"ECDF control (seed={selected_seed}, 4 ranks)", loc="left", fontweight="semibold")
    ecdf_axis.legend(frameon=False, fontsize=8, loc="lower right")

    density_grid = np.linspace(lo, hi, 520)
    kde_axis.plot(density_grid, histogram_kde(reference_log, density_grid), color="#111827", linewidth=2.1, label="Historical")
    for grid in sorted(selected):
        log_values = np.log10(load_v3(selected[grid]["lifetime_file"])[2])
        kde_axis.plot(density_grid, histogram_kde(log_values, density_grid), color=colors[grid], linewidth=1.55, label=f"grid={grid}")
    kde_axis.set(xlim=(lo, hi), xlabel="log10(lifetime / s)", ylabel="Probability density")
    kde_axis.set_title("KDE: 20×20 changes the mode balance", loc="left", fontweight="semibold")
    kde_axis.legend(frameon=False, fontsize=8)

    rank4_rows = [row for row in metric_rows if int(row["mpi_size"]) == 4]
    for grid in sorted({int(row["interpolation_points"]) for row in rank4_rows}):
        rows = [row for row in rank4_rows if int(row["interpolation_points"]) == grid]
        x = np.full(len(rows), grid_positions[grid], dtype=float)
        ks_axis.scatter(x, [float(row["ks_d"]) for row in rows], color=colors[grid], s=42, zorder=3)
    ks_axis.set_xticks(list(grid_positions.values()), ["0", "20", "1000", "2000"])
    ks_axis.set(xlabel="Interpolation points per axis", ylabel="KS D vs historical")
    ks_axis.set_title("Three independent seeds per grid", loc="left", fontweight="semibold")

    reference_peaks, reference_boundaries = detect_logtime_peaks(reference_log)
    reference_weight = float(mode_fractions(reference_log, reference_boundaries)[0])
    weight_axis.axhline(reference_weight, color="#111827", linestyle="--", linewidth=1.4, label=f"Historical {reference_weight:.3f}")
    for grid in sorted({int(row["interpolation_points"]) for row in rank4_rows}):
        rows = [row for row in rank4_rows if int(row["interpolation_points"]) == grid]
        weights = [parse_vector(str(row["mode_fractions"]))[0] for row in rows]
        x_value = grid_positions[grid]
        weight_axis.scatter(np.full(len(weights), x_value), weights, color=colors[grid], s=42, zorder=3)
    weight_axis.set_xticks(list(grid_positions.values()), ["0", "20", "1000", "2000"])
    weight_axis.set(xlabel="Interpolation points per axis", ylabel="Short-mode fraction")
    weight_axis.set_title("Short-lifetime mode weight", loc="left", fontweight="semibold")
    weight_axis.legend(frameon=False, fontsize=8)

    for axis in axes.flat:
        axis.grid(axis="y", color="#CBD5E1", linewidth=0.65, alpha=0.7)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    fig.suptitle("C-point interpolation, seed, and MPI controls", x=0.07, ha="left", fontsize=16, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-root", type=Path, required=True)
    parser.add_argument("--mpi-root", type=Path)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    reference = resolve_lifetime_file(args.reference)
    runs = discover(args.matrix_root.resolve())
    if args.mpi_root:
        runs.extend(discover(args.mpi_root.resolve()))
    if not runs:
        raise RuntimeError("no C control runs discovered")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    metric_rows: list[dict[str, Any]] = []
    quantile_rows: list[dict[str, Any]] = []
    for run in runs:
        label = f"C-g{run['grid']}-s{run['seed']}-r{run['mpi_size']}"
        row, current_quantiles, _ = analyze(Comparison(label, run["directory"], reference))
        metric_rows.append(row)
        for quantile_row in current_quantiles:
            quantile_row.update({"grid": run["grid"], "seed": run["seed"], "mpi_size": run["mpi_size"]})
        quantile_rows.extend(current_quantiles)
    metric_rows.sort(key=lambda row: (int(row["mpi_size"]), int(row["interpolation_points"]), int(row["base_seed"])))
    write_csv(output_dir / "run_metrics.csv", metric_rows)
    write_csv(output_dir / "quantile_shifts.csv", quantile_rows)

    grouped: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in metric_rows:
        if int(row["mpi_size"]) == 4:
            grouped[int(row["interpolation_points"])].append(row)
    grid_rows: list[dict[str, Any]] = []
    for grid, rows in sorted(grouped.items()):
        peaks = [parse_vector(str(row["peak_log10_s"])) for row in rows]
        weights = [parse_vector(str(row["mode_fractions"])) for row in rows]
        grid_rows.append(
            {
                "interpolation_points": grid,
                "seeds": len(rows),
                "mean_n_complete": mean([float(row["n_complete"]) for row in rows]),
                "ks_d_mean": mean([float(row["ks_d"]) for row in rows]),
                "ks_d_range": range_text([float(row["ks_d"]) for row in rows]),
                "peak1_log10_s_mean": mean([value[0] for value in peaks]),
                "peak2_log10_s_mean": mean([value[1] for value in peaks]),
                "short_mode_weight_mean": mean([value[0] for value in weights]),
                "short_mode_weight_range": range_text([value[0] for value in weights]),
                "mean_scatterings_completed": mean([float(row["mean_scatterings_completed"]) for row in rows]),
                "numerical_failure_rate_mean": mean([float(row["numerical_failure_rate"]) for row in rows]),
                "captured_invalid_fraction_mean": mean([float(row["captured_invalid_fraction"]) for row in rows]),
            }
        )
    write_csv(output_dir / "grid_summary.csv", grid_rows)

    quantile_summary: list[dict[str, Any]] = []
    for grid in sorted(grouped):
        for quantile in QUANTILES:
            values = [
                float(row["shift_dex"])
                for row in quantile_rows
                if int(row["mpi_size"]) == 4 and int(row["grid"]) == grid and float(row["quantile"]) == quantile
            ]
            quantile_summary.append(
                {
                    "interpolation_points": grid,
                    "quantile": quantile,
                    "mean_shift_dex": mean(values),
                    "shift_range_dex": range_text(values),
                }
            )
    write_csv(output_dir / "grid_quantile_summary.csv", quantile_summary)

    indexed = {(run["grid"], run["seed"], run["mpi_size"]): run for run in runs}
    pairwise_rows: list[dict[str, Any]] = []
    for seed in sorted({run["seed"] for run in runs}):
        left = indexed.get((20, seed, 4))
        right = indexed.get((1000, seed, 4))
        if left and right:
            left_log = np.log10(load_v3(left["lifetime_file"])[2])
            right_log = np.log10(load_v3(right["lifetime_file"])[2])
            result = ks_2samp(left_log, right_log)
            pairwise_rows.append({"comparison": "grid20_vs_grid1000", "seed": seed, "grid": "20_vs_1000", "left_mpi": 4, "right_mpi": 4, "ks_d": result.statistic, "ks_p_value": result.pvalue, "ks_location_log10_s": result.statistic_location})
    common = sorted(set((run["grid"], run["seed"]) for run in runs if run["mpi_size"] == 4) & set((run["grid"], run["seed"]) for run in runs if run["mpi_size"] == 30))
    for grid, seed in common:
        left_log = np.log10(load_v3(indexed[(grid, seed, 4)]["lifetime_file"])[2])
        right_log = np.log10(load_v3(indexed[(grid, seed, 30)]["lifetime_file"])[2])
        result = ks_2samp(left_log, right_log)
        pairwise_rows.append({"comparison": "mpi4_vs_mpi30", "seed": seed, "grid": grid, "left_mpi": 4, "right_mpi": 30, "ks_d": result.statistic, "ks_p_value": result.pvalue, "ks_location_log10_s": result.statistic_location})
    write_csv(output_dir / "pairwise_controls.csv", pairwise_rows)
    plot_controls(output_dir / "c_control_summary.png", runs, reference, metric_rows)

    manifest = {
        "schema_version": "c-control-analysis-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "reference": str(reference),
        "reference_sha256": sha256(reference),
        "inputs": [{"file": str(run["lifetime_file"]), "sha256": sha256(run["lifetime_file"])} for run in runs],
        "outputs": ["run_metrics.csv", "quantile_shifts.csv", "grid_summary.csv", "grid_quantile_summary.csv", "pairwise_controls.csv", "c_control_summary.png"],
    }
    with (output_dir / "analysis_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Saved C control summary to {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

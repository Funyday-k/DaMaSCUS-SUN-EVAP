#!/usr/bin/env python3
"""Reproduce evaporation-lifetime comparison tables and ECDF/KDE figures.

Example:
    python3 scripts/diagnostics/analyze_evaporation_reproduction.py \
      --comparison A,evaporation_secondary_peak_A_10000/results_-2.000000_-32.000000,/path/to/A/evaporation_times.txt \
      --comparison B,evaporation_secondary_peak_B_10000/results_-2.000000_-34.000000,/path/to/B/evaporation_times.txt \
      --output-dir analysis/abc_reproduction

Each comparison is ``LABEL,RUN_PATH,REFERENCE_FILE``. ``RUN_PATH`` may be a
result directory, an evaporation_times.txt file, or a parent containing one
``results_*/evaporation_times.txt`` file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.stats import ks_2samp


QUANTILES = (0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)


@dataclass(frozen=True)
class Comparison:
    label: str
    run_path: Path
    reference_path: Path


def parse_comparison(spec: str) -> Comparison:
    parts = spec.split(",", 2)
    if len(parts) != 3 or not parts[0].strip():
        raise argparse.ArgumentTypeError(
            "comparison must be LABEL,RUN_PATH,REFERENCE_FILE"
        )
    return Comparison(parts[0].strip(), Path(parts[1]), Path(parts[2]))


def resolve_lifetime_file(path: Path) -> Path:
    path = path.expanduser().resolve()
    if path.is_file():
        return path
    direct = path / "evaporation_times.txt"
    if direct.is_file():
        return direct
    candidates = sorted(path.glob("results_*/evaporation_times.txt"))
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"{path}: expected one results_*/evaporation_times.txt, found {len(candidates)}"
        )
    return candidates[0]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_v3(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.loadtxt(path, comments="#")
    if values.ndim == 1:
        values = values.reshape(1, -1)
    if values.shape[1] != 3:
        raise ValueError(f"{path}: expected rank, trajectory_id, lifetime columns")
    rank = values[:, 0].astype(np.int64)
    trajectory_id = values[:, 1].astype(np.int64)
    lifetime = values[:, 2].astype(float)
    valid = np.isfinite(lifetime) & (lifetime > 0.0)
    return rank[valid], trajectory_id[valid], lifetime[valid]


def detect_logtime_peaks(log_times: np.ndarray) -> tuple[np.ndarray, list[float]]:
    """Match the historical multimodality diagnostic peak rule."""
    if log_times.size < 30:
        return np.array([], dtype=float), []
    lo = float(np.floor(log_times.min() * 4.0) / 4.0)
    hi = float(np.ceil(log_times.max() * 4.0) / 4.0)
    if hi <= lo:
        hi = lo + 1.0
    bin_count = int(np.clip(round((hi - lo) / 0.12), 35, 95))
    counts, edges = np.histogram(log_times, bins=np.linspace(lo, hi, bin_count + 1))
    centers = 0.5 * (edges[:-1] + edges[1:])
    smooth = gaussian_filter1d(counts.astype(float), sigma=1.15, mode="nearest")
    if smooth.size < 3 or smooth.max() <= 0:
        return np.array([], dtype=float), []
    width = float(edges[1] - edges[0])
    peak_indices, _ = find_peaks(
        smooth,
        prominence=max(3.0, 0.055 * float(smooth.max())),
        distance=max(2, int(round(0.45 / width))),
    )
    if peak_indices.size:
        minimum_height = max(5.0, 0.075 * float(smooth.max()))
        peak_indices = peak_indices[smooth[peak_indices] >= minimum_height]
    peak_indices = peak_indices[np.argsort(centers[peak_indices])]
    boundaries: list[float] = []
    for left, right in zip(peak_indices[:-1], peak_indices[1:]):
        boundaries.append(float(centers[left + int(np.argmin(smooth[left : right + 1]))]))
    return centers[peak_indices], boundaries


def mode_fractions(log_times: np.ndarray, boundaries: list[float]) -> np.ndarray:
    labels = np.searchsorted(np.asarray(boundaries), log_times, side="right")
    count = max(1, len(boundaries) + 1)
    return np.asarray([np.mean(labels == index) for index in range(count)])


def parse_header_metrics(path: Path) -> dict[str, str]:
    metrics: dict[str, str] = {}
    if not path.is_file():
        return metrics
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            payload = line[1:].strip()
            if "=" in payload:
                key, value = payload.split("=", 1)
                metrics[key.strip()] = value.strip()
    return metrics


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def summary_metrics(path: Path) -> dict[str, float | int]:
    if not path.is_file():
        return {}
    completed_scatterings: list[int] = []
    captured_rows = 0
    invalid_rows = 0
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            captured_rows += 1
            if row.get("event_observed") == "1":
                completed_scatterings.append(int(row["n_scatter_total"]))
            else:
                invalid_rows += 1
    return {
        "captured_summary_rows": captured_rows,
        "captured_invalid_rows": invalid_rows,
        "captured_invalid_fraction": invalid_rows / captured_rows if captured_rows else math.nan,
        "mean_scatterings_completed": float(np.mean(completed_scatterings))
        if completed_scatterings
        else math.nan,
    }


def numeric(value: str | None, default: float = math.nan) -> float:
    try:
        return float(value) if value is not None else default
    except ValueError:
        return default


def format_vector(values: np.ndarray | list[float]) -> str:
    return ";".join(f"{float(value):.9g}" for value in values)


def ecdf(values: np.ndarray, maximum_points: int = 6000) -> tuple[np.ndarray, np.ndarray]:
    sorted_values = np.sort(values)
    if sorted_values.size <= maximum_points:
        indices = np.arange(sorted_values.size)
    else:
        indices = np.unique(np.linspace(0, sorted_values.size - 1, maximum_points).astype(int))
    return sorted_values[indices], (indices + 1) / sorted_values.size


def histogram_kde(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    edges = np.linspace(grid[0], grid[-1], grid.size + 1)
    density, _ = np.histogram(values, bins=edges, density=True)
    width = float(edges[1] - edges[0])
    bandwidth = float(np.std(values, ddof=1) * max(values.size, 2) ** (-0.2))
    smooth = gaussian_filter1d(density, sigma=max(0.5, bandwidth / width), mode="nearest")
    return smooth


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def analyze(comparison: Comparison) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    run_file = resolve_lifetime_file(comparison.run_path)
    reference_file = resolve_lifetime_file(comparison.reference_path)
    run_rank, run_id, run_lifetime = load_v3(run_file)
    _, _, reference_lifetime = load_v3(reference_file)
    run_log = np.log10(run_lifetime)
    reference_log = np.log10(reference_lifetime)
    run_peaks, run_boundaries = detect_logtime_peaks(run_log)
    reference_peaks, reference_boundaries = detect_logtime_peaks(reference_log)
    run_fractions = mode_fractions(run_log, run_boundaries)
    reference_fractions = mode_fractions(reference_log, reference_boundaries)
    ks = ks_2samp(run_log, reference_log)
    result_dir = run_file.parent
    metadata = read_json(result_dir / "run_metadata.json")
    bincount = parse_header_metrics(result_dir / "bincount.txt")
    extra = summary_metrics(result_dir / "trajectory_summary.tsv")
    keys = np.column_stack((run_rank, run_id))
    duplicate_keys = int(keys.shape[0] - np.unique(keys, axis=0).shape[0])
    simulated = numeric(bincount.get("total_trajectories"))
    failures = numeric(bincount.get("numerical_failures"))
    row: dict[str, Any] = {
        "label": comparison.label,
        "n_complete": int(run_lifetime.size),
        "reference_n_complete": int(reference_lifetime.size),
        "duplicate_rank_trajectory_keys": duplicate_keys,
        "peak_log10_s": format_vector(run_peaks),
        "reference_peak_log10_s": format_vector(reference_peaks),
        "mode_fractions": format_vector(run_fractions),
        "reference_mode_fractions": format_vector(reference_fractions),
        "mode_boundaries_log10_s": format_vector(run_boundaries),
        "ks_d": float(ks.statistic),
        "ks_p_value": float(ks.pvalue),
        "ks_location_log10_s": float(ks.statistic_location),
        "mean_scatterings_all_simulated": numeric(bincount.get("average_scatterings")),
        "mean_scatterings_completed": extra.get("mean_scatterings_completed", math.nan),
        "numerical_failure_count": int(failures) if math.isfinite(failures) else "",
        "numerical_failure_rate": failures / simulated
        if math.isfinite(failures) and math.isfinite(simulated) and simulated > 0
        else math.nan,
        "captured_invalid_rows": extra.get("captured_invalid_rows", ""),
        "captured_invalid_fraction": extra.get("captured_invalid_fraction", math.nan),
        "interpolation_points": metadata.get("interpolation_points", ""),
        "base_seed": metadata.get("base_seed", ""),
        "mpi_size": metadata.get("mpi_size", ""),
        "git_commit": metadata.get("git_commit", ""),
        "git_dirty": metadata.get("git_dirty", ""),
        "run_file": str(run_file),
        "reference_file": str(reference_file),
    }
    quantile_rows: list[dict[str, Any]] = []
    for quantile in QUANTILES:
        current = float(np.quantile(run_log, quantile))
        reference = float(np.quantile(reference_log, quantile))
        quantile_rows.append(
            {
                "label": comparison.label,
                "quantile": quantile,
                "run_log10_s": current,
                "reference_log10_s": reference,
                "shift_dex": current - reference,
            }
        )
    plot_data = {
        "label": comparison.label,
        "run_log": run_log,
        "reference_log": reference_log,
        "ks_d": float(ks.statistic),
        "ks_location": float(ks.statistic_location),
    }
    return row, quantile_rows, plot_data


def plot_comparisons(path: Path, plot_data: list[dict[str, Any]], title: str) -> None:
    rows = len(plot_data)
    fig, axes = plt.subplots(rows, 2, figsize=(12.5, max(4.2, 3.55 * rows)), squeeze=False, dpi=160)
    for row_index, data in enumerate(plot_data):
        run_log = data["run_log"]
        reference_log = data["reference_log"]
        lo = float(min(np.quantile(run_log, 0.001), np.quantile(reference_log, 0.001)))
        hi = float(max(np.quantile(run_log, 0.999), np.quantile(reference_log, 0.999)))
        padding = max(0.15, 0.04 * (hi - lo))
        lo -= padding
        hi += padding
        current_x, current_y = ecdf(run_log)
        reference_x, reference_y = ecdf(reference_log)
        ecdf_axis = axes[row_index, 0]
        ecdf_axis.plot(reference_x, reference_y, color="#52606D", linewidth=1.7, label="Historical reference")
        ecdf_axis.plot(current_x, current_y, color="#2B6CB0", linewidth=1.7, label="Current run")
        ecdf_axis.axvline(data["ks_location"], color="#B7791F", linestyle="--", linewidth=1.25)
        ecdf_axis.text(
            0.03,
            0.93,
            f"KS D={data['ks_d']:.4f} at {data['ks_location']:.3f}",
            transform=ecdf_axis.transAxes,
            ha="left",
            va="top",
            color="#805A16",
            fontsize=9,
        )
        ecdf_axis.set_xlim(lo, hi)
        ecdf_axis.set_ylim(0.0, 1.0)
        ecdf_axis.set_title(f"{data['label']} — empirical CDF", loc="left", fontsize=11, fontweight="semibold")
        ecdf_axis.set_xlabel("log10(final-unbinding lifetime / s)")
        ecdf_axis.set_ylabel("Cumulative fraction")
        ecdf_axis.legend(frameon=False, loc="lower right", fontsize=8.5)

        kde_axis = axes[row_index, 1]
        grid = np.linspace(lo, hi, 500)
        kde_axis.plot(grid, histogram_kde(reference_log, grid), color="#52606D", linewidth=1.9, label="Historical reference")
        kde_axis.plot(grid, histogram_kde(run_log, grid), color="#2B6CB0", linewidth=1.9, label="Current run")
        kde_axis.set_xlim(lo, hi)
        kde_axis.set_title(f"{data['label']} — Scott-bandwidth KDE", loc="left", fontsize=11, fontweight="semibold")
        kde_axis.set_xlabel("log10(final-unbinding lifetime / s)")
        kde_axis.set_ylabel("Probability density")
        kde_axis.legend(frameon=False, loc="upper right", fontsize=8.5)
        for axis in (ecdf_axis, kde_axis):
            axis.grid(axis="y", color="#CBD5E1", linewidth=0.65, alpha=0.65)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
            axis.spines["left"].set_color("#7B8794")
            axis.spines["bottom"].set_color("#7B8794")
            axis.tick_params(colors="#52606D", labelsize=8.5)
    fig.suptitle(title, x=0.08, ha="left", fontsize=15, fontweight="semibold", color="#1F2933")
    fig.text(
        0.99,
        0.01,
        "ECDF uses all rows; KDE display range is the combined 0.1%–99.9% interval.",
        ha="right",
        va="bottom",
        fontsize=8,
        color="#7B8794",
    )
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.07, top=0.92, hspace=0.45, wspace=0.25)
    fig.savefig(path)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", action="append", type=parse_comparison, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--title", default="Evaporation lifetime reproduction")
    args = parser.parse_args()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict[str, Any]] = []
    quantile_rows: list[dict[str, Any]] = []
    plot_data: list[dict[str, Any]] = []
    manifest_inputs: list[dict[str, Any]] = []
    for comparison in args.comparison:
        row, current_quantiles, current_plot = analyze(comparison)
        metric_rows.append(row)
        quantile_rows.extend(current_quantiles)
        plot_data.append(current_plot)
        run_file = resolve_lifetime_file(comparison.run_path)
        reference_file = resolve_lifetime_file(comparison.reference_path)
        manifest_inputs.append(
            {
                "label": comparison.label,
                "run_file": str(run_file),
                "run_sha256": sha256(run_file),
                "reference_file": str(reference_file),
                "reference_sha256": sha256(reference_file),
            }
        )
    write_csv(output_dir / "comparison_metrics.csv", metric_rows)
    write_csv(output_dir / "quantile_shifts.csv", quantile_rows)
    plot_comparisons(output_dir / "ecdf_kde_comparison.png", plot_data, args.title)
    manifest = {
        "schema_version": "evaporation-reproduction-analysis-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "numpy": np.__version__,
        "cwd": os.getcwd(),
        "inputs": manifest_inputs,
        "outputs": [
            "comparison_metrics.csv",
            "quantile_shifts.csv",
            "ecdf_kde_comparison.png",
        ],
    }
    with (output_dir / "analysis_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Saved analysis to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

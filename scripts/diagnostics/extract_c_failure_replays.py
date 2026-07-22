#!/usr/bin/env python3
"""Extract exact C-point failure replays and duration-matched successes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


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


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def key(row: dict[str, str]) -> tuple[int, int]:
    return int(row["rank"]), int(row["trajectory_id"])


def observed_log_duration(row: dict[str, str]) -> float:
    if row["event_observed"] == "1":
        duration = float(row["lifetime_validated_escape_s"])
    else:
        duration = float(row["t_censor_s"]) - float(row["t_capture_s"])
    if not math.isfinite(duration) or duration <= 0:
        raise ValueError(f"invalid observed duration for {key(row)}: {duration}")
    return math.log10(duration)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-summary", type=Path, required=True)
    parser.add_argument("--replay-result-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    original_summary = args.original_summary.resolve()
    replay_dir = args.replay_result_dir.resolve()
    replay_summary = replay_dir / "trajectory_summary.tsv"
    replay_events = replay_dir / "trajectory_events.tsv"
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    original_fields, original_rows = read_tsv(original_summary)
    replay_fields, replay_rows = read_tsv(replay_summary)
    original_failures = {key(row): row for row in original_rows if row["event_observed"] == "0"}
    replay_by_key = {key(row): row for row in replay_rows}
    missing = sorted(set(original_failures) - set(replay_by_key))
    if missing:
        raise RuntimeError(f"missing original failure keys in replay: {missing}")
    failures = [replay_by_key[item] for item in sorted(original_failures)]
    for row in failures:
        original = original_failures[key(row)]
        if row["termination_reason"] != original["termination_reason"]:
            raise RuntimeError(f"termination mismatch for {key(row)}")
        if row["max_scaled_ballistic_energy_drift"] != original["max_scaled_ballistic_energy_drift"]:
            raise RuntimeError(f"energy-drift mismatch for {key(row)}")
        if row["trace_written"] != "1" or not row["rng_state_before_initial_conditions"] or not row["rng_state_before_simulation"]:
            raise RuntimeError(f"missing replay state for {key(row)}")

    successes = [row for row in replay_rows if row["event_observed"] == "1"]
    available = {key(row): row for row in successes}
    pairs: list[dict[str, Any]] = []
    matched_successes: list[dict[str, str]] = []
    for failure in sorted(failures, key=observed_log_duration):
        failure_log = observed_log_duration(failure)
        matched_key, matched = min(
            available.items(), key=lambda item: abs(observed_log_duration(item[1]) - failure_log)
        )
        del available[matched_key]
        success_log = observed_log_duration(matched)
        matched_successes.append(matched)
        pairs.append(
            {
                "failure_rank": failure["rank"],
                "failure_trajectory_id": failure["trajectory_id"],
                "failure_censor_duration_s": float(failure["t_censor_s"]) - float(failure["t_capture_s"]),
                "failure_log10_censor_duration_s": failure_log,
                "failure_scatter_count": failure["n_scatter_total"],
                "failure_max_scaled_energy_drift": failure["max_scaled_ballistic_energy_drift"],
                "success_rank": matched["rank"],
                "success_trajectory_id": matched["trajectory_id"],
                "success_validated_lifetime_s": matched["lifetime_validated_escape_s"],
                "success_log10_lifetime_s": success_log,
                "success_scatter_count": matched["n_scatter_total"],
                "absolute_log10_duration_difference": abs(success_log - failure_log),
                "within_0p25_dex": int(abs(success_log - failure_log) <= 0.25),
            }
        )

    write_tsv(output_dir / "failure_replays.tsv", replay_fields, failures)
    write_tsv(output_dir / "matched_success_replays.tsv", replay_fields, matched_successes)
    write_tsv(output_dir / "matched_pairs.tsv", list(pairs[0]), pairs)

    selected_keys = {key(row) for row in failures + matched_successes}
    event_fields, event_rows = read_tsv(replay_events)
    selected_events = [row for row in event_rows if key(row) in selected_keys]
    write_tsv(output_dir / "selected_trajectory_events.tsv", event_fields, selected_events)

    differences = [float(row["absolute_log10_duration_difference"]) for row in pairs]
    manifest = {
        "schema_version": "c-failure-replay-extract-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "original_summary": str(original_summary),
        "original_summary_sha256": sha256(original_summary),
        "replay_summary": str(replay_summary),
        "replay_summary_sha256": sha256(replay_summary),
        "replay_events": str(replay_events),
        "replay_events_sha256": sha256(replay_events),
        "original_failure_count": len(original_failures),
        "exactly_reproduced_failure_count": len(failures),
        "termination_reason_exact_count": sum(
            replay_by_key[item]["termination_reason"] == original_failures[item]["termination_reason"]
            for item in original_failures
        ),
        "scaled_energy_drift_exact_count": sum(
            replay_by_key[item]["max_scaled_ballistic_energy_drift"]
            == original_failures[item]["max_scaled_ballistic_energy_drift"]
            for item in original_failures
        ),
        "failure_replay_state_count": sum(
            replay_by_key[item]["trace_written"] == "1"
            and bool(replay_by_key[item]["rng_state_before_initial_conditions"])
            and bool(replay_by_key[item]["rng_state_before_simulation"])
            for item in original_failures
        ),
        "matched_success_count": len(matched_successes),
        "matches_within_0p25_dex": sum(float(value) <= 0.25 for value in differences),
        "maximum_match_difference_dex": max(differences),
        "matching_definition": "greedy one-to-one nearest log10(success validated lifetime) to log10(failure censor duration)",
        "matching_caveat": "a censor duration is a lower bound, not the unobserved physical lifetime",
        "selected_event_rows": len(selected_events),
    }
    with (output_dir / "replay_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Extracted {len(failures)} failures and {len(matched_successes)} controls to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

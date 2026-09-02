#!/usr/bin/env python3
"""Cheap replay analysis for the recorded P3.II.delta frontier outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np


T_975_DF15 = 2.131449545559323
P3II_DELTA_BOUND = Fraction(1265833333, 5000000000)
P3II_RANGE_BOUND = Fraction(913600001, 3600000000)


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "primegaps.p3ii-delta-frontier-raw.v1":
        raise ValueError(f"unexpected schema in {path}")
    if "translated-simplex" not in payload.get("method", ""):
        raise ValueError(f"superseded estimator in {path}")
    return payload


def arrays(payload: dict[str, object]):
    results = payload["results"]
    x = np.asarray([item["A_max"] for item in results], dtype=float)
    replicates = np.asarray(
        [item["lambda_48_replicates"] for item in results], dtype=float
    )
    return results, x, replicates, replicates.mean(axis=1)


def paired_linear_root(x: np.ndarray, replicates: np.ndarray):
    center = 0.25361
    mean_values = replicates.mean(axis=1)
    slope, intercept = np.polyfit(x - center, mean_values, 1)
    root = center + (1.0 - intercept) / slope
    replicate_roots = []
    replicate_slopes = []
    for column in range(replicates.shape[1]):
        seed_slope, seed_intercept = np.polyfit(
            x - center, replicates[:, column], 1
        )
        replicate_slopes.append(float(seed_slope))
        replicate_roots.append(float(center + (1.0 - seed_intercept) / seed_slope))
    replicate_roots = np.asarray(replicate_roots)
    root_mean = float(replicate_roots.mean())
    root_se = float(replicate_roots.std(ddof=1) / math.sqrt(len(replicate_roots)))
    critical = T_975_DF15 if len(replicate_roots) == 16 else 1.96
    return {
        "fit_root": float(root),
        "paired_replicate_root_mean": root_mean,
        "paired_replicate_root_standard_error": root_se,
        "paired_replicate_root_95_interval": [
            root_mean - critical * root_se,
            root_mean + critical * root_se,
        ],
        "slope": float(slope),
        "replicate_slopes": replicate_slopes,
        "replicate_roots": replicate_roots.tolist(),
    }


def nearest(results, target: float):
    item = min(results, key=lambda row: abs(row["A_max"] - target))
    return {
        "A_max": item["A_max"],
        "lambda_48": item["lambda_48"],
        "standard_error": item["lambda_48_standard_error"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--curve", type=Path, required=True)
    parser.add_argument("--crossing", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    curve = load(args.curve)
    crossing = load(args.crossing)
    curve_results, curve_x, curve_reps, curve_mean = arrays(curve)
    crossing_results, crossing_x, crossing_reps, _crossing_mean = arrays(crossing)
    if np.any(np.diff(curve_x) <= 0) or np.any(np.diff(curve_mean) <= 0):
        raise ValueError("recorded frontier is not strictly increasing")
    if not (curve_mean[0] < 1 < curve_mean[-1]):
        raise ValueError("recorded frontier does not bracket one")
    if any(item["best_tolerance"] != 1e-13 for item in curve_results):
        raise ValueError("the maximizing bank tolerance changes on the frontier")

    crossing_fit = paired_linear_root(crossing_x, crossing_reps)
    root = crossing_fit["paired_replicate_root_mean"]
    p3_delta = float(P3II_DELTA_BOUND)
    p3_range = float(P3II_RANGE_BOUND)
    summary = {
        "schema": "primegaps.p3ii-delta-frontier-summary.v1",
        "status": "numerical-screening-not-certificate",
        "target": "P3.II.delta",
        "response": "lambda_48",
        "raw_inputs": [
            {"path": str(args.curve), "sha256": file_hash(args.curve)},
            {"path": str(args.crossing), "sha256": file_hash(args.crossing)},
        ],
        "frontier": [
            {
                "A_max": item["A_max"],
                "lambda_48": item["lambda_48"],
                "standard_error": item["lambda_48_standard_error"],
            }
            for item in curve_results
        ],
        "crossing_fit": crossing_fit,
        "constraint_transition": {
            "relaxed_constraint": "P3.II.delta",
            "P3.II.delta_binding_A": {
                "fraction": str(P3II_DELTA_BOUND),
                "decimal": p3_delta,
            },
            "next_constraint": "P3.II.range",
            "P3.II.range_binding_A": {
                "fraction": str(P3II_RANGE_BOUND),
                "decimal": p3_range,
            },
        },
        "relaxation_to_cross_one": {
            "from_published_A_0.253": root - 0.253,
            "beyond_P3.II.delta_binding_A": root - p3_delta,
            "fraction_of_available_one_constraint_interval": (
                (root - p3_delta) / (p3_range - p3_delta)
            ),
            "remaining_A_headroom_before_P3.II.range": p3_range - root,
        },
        "selected_points": {
            "published_A": nearest(curve_results, 0.253),
            "P3.II.delta_boundary": nearest(curve_results, p3_delta),
            "lower_crossing_bracket": nearest(crossing_results, 0.253606),
            "upper_crossing_bracket": nearest(crossing_results, 0.253610),
            "prior_search_endpoint": nearest(curve_results, 0.2537),
            "P3.II.range_boundary": nearest(curve_results, p3_range),
        },
        "verification": {
            "curve_strictly_increasing": True,
            "crossing_bracketed": True,
            "best_tolerance_constant": "1e-13",
            "common_seed_pairing_used_for_root_uncertainty": True,
            "point_count": len(curve_results),
            "crossing_point_count": len(crossing_results),
            "crossing_replicate_count": crossing_reps.shape[1],
        },
        "caveats": [
            "The curve is a degree-21 numerical vector-bank/QMC screen, not an exact rational certificate or a proof of the global optimum over all supports.",
            "The reported root interval is a randomized-QMC sampling interval for a local linear fit; it excludes test-function truncation and outer-optimization error.",
            "The exact analytic constraint endpoints use Proposition 3 with theorem epsilon=1e-10; the floating endpoint score is evaluated at the nearest binary64 geometry.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.output)
    print(json.dumps({
        "root": root,
        "root_95_interval": crossing_fit["paired_replicate_root_95_interval"],
        "slope": crossing_fit["slope"],
        "next_constraint": "P3.II.range",
        "next_constraint_A": p3_range,
    }, sort_keys=True))


if __name__ == "__main__":
    main()

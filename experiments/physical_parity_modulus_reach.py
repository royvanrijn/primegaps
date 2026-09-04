#!/usr/bin/env python3
"""Exploratory modulus-reach truncation of the k=39 physical full face.

This is a research-work script.  It imports the tracked physical parity model,
but replaces each erased-coordinate feature E_i F(Y) by

    E_{i,theta} F(Y) = integral F(Y +_i X)
        1_{rho_star (|Y| + |Y +_i X|) <= theta} dnu(X).

On the midpoint grid the conservative cell rule uses the two upper total
endpoints.  The returned quadratic form is sum_i ||E_{i,theta} F||^2 over the
same unrestricted full-face masks as Jfull.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys
from time import monotonic

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "experiments" / "physical_parity_viability.py"


def load_model_module():
    spec = importlib.util.spec_from_file_location("physical_parity_viability", SOURCE)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


parity = load_model_module()


class ReachModel(parity.FloatPhysicalCapModel):
    """Physical cap model with the erased edge cut at a modulus exponent."""

    def _reach_correlations(self, shell_index: int, theta: Fraction):
        shell = self.shells["outer"][shell_index]
        mask = self.outer_masks[shell_index]
        theta = Fraction(theta)
        # Upper endpoints are (q+k)h and (r+k-1)h, where q=r+a.
        maximum_index_sum = int(theta // (self.upstream.RHO_STAR * self.hq))
        answer = {}
        midpoint_powers = tuple(self.midpoints**exponent for exponent in range(7))
        survival = self.survival(shell.ceiling)
        base_fiber = self.root * survival
        for degree in range(7):
            radial = self.radial_powers[degree] * mask
            for exponent in range(7):
                values = np.zeros(self.n, dtype=float)
                fiber = base_fiber * midpoint_powers[exponent]
                for retained_index in range(self.n):
                    last_outer_index = min(
                        self.n - 1,
                        maximum_index_sum - retained_index - (2 * self.k - 1),
                    )
                    if last_outer_index < retained_index:
                        continue
                    outer = radial[retained_index : last_outer_index + 1]
                    erased = fiber[: last_outer_index - retained_index + 1]
                    values[retained_index] = float(np.dot(outer, erased))
                answer[degree, exponent] = values
        return answer

    def _reach_layer_features(self, theta: Fraction):
        shell_correlations = tuple(
            self._reach_correlations(index, theta)
            for index in range(len(self.shells["outer"]))
        )
        by_cap = {}
        for cap in self.caps:
            rows = {}
            for shell_index, shell in enumerate(self.shells["outer"]):
                if cap > shell.ceiling:
                    continue
                correlations = shell_correlations[shell_index]
                for basis_index, (signature, degree) in enumerate(self.descriptors):
                    for remaining, exponent, multiplicity in self.upstream.fiber_splits(
                        signature
                    ):
                        if remaining not in rows:
                            rows[remaining] = np.zeros(
                                (self.basis_size, self.n), dtype=float
                            )
                        rows[remaining][basis_index] += (
                            multiplicity * correlations[degree, exponent]
                        )
            by_cap[cap] = rows
        return by_cap

    def reach_face_matrix(self, theta: Fraction) -> np.ndarray:
        features_by_cap = self._reach_layer_features(theta)
        result = np.zeros((self.basis_size, self.basis_size), dtype=float)
        previous = None
        for layer, cap in enumerate(self.caps):
            features = features_by_cap[cap]
            full_mask = self.inner_allowed[2, layer]
            if not features or not np.any(full_mask):
                previous = cap
                continue
            differences = {}
            signatures = tuple(sorted(features, key=lambda value: (len(value), value)))
            for right_index, right_signature in enumerate(signatures):
                right = features[right_signature]
                for left_index in range(right_index + 1):
                    left_signature = signatures[left_index]
                    left = features[left_signature]
                    joined = tuple(sorted(left_signature + right_signature))
                    if joined not in differences:
                        values = self.moment(cap, self.k - 1, joined).copy()
                        if previous is not None:
                            values -= self.moment(previous, self.k - 1, joined)
                        scale = max(
                            float(np.max(np.abs(values))), np.finfo(float).tiny
                        )
                        if float(np.min(values)) < -2e-6 * scale:
                            raise ArithmeticError(
                                "a cap-layer moment difference became materially negative"
                            )
                        differences[joined] = np.maximum(values, 0.0)
                    moment = differences[joined]
                    weighted_right = right[:, full_mask] * moment[full_mask]
                    block = left[:, full_mask] @ weighted_right.T
                    result += block
                    if left_index != right_index:
                        result += block.T
            previous = cap
        normalization = self.k * self.h / self.Z
        return normalization * (result + result.T) / 2


def clean_eigen_result(value):
    return {key: item for key, item in value.items() if key != "vector"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--intervals", type=int, required=True)
    parser.add_argument("--theta", type=Fraction, nargs="+", required=True)
    parser.add_argument("--cutoff", type=float, default=1e-11)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--matrix-output", type=Path)
    args = parser.parse_args()

    upstream = parity.load_inputs()
    model = ReachModel(upstream, dimension=39, intervals=args.intervals)
    started = monotonic()
    metric = model.mass_matrix()
    rows = []
    matrices = {}
    for theta in sorted(set(args.theta)):
        item_started = monotonic()
        matrix = model.reach_face_matrix(theta)
        eigen = parity.spectral_generalized_maximum(matrix, metric, args.cutoff)
        score = float(upstream.RHO_STAR) * float(eigen["quotient"])
        rows.append(
            {
                "theta": str(theta),
                "theta_float": float(theta),
                "score": score,
                "eigen": clean_eigen_result(eigen),
                "elapsed_seconds": monotonic() - item_started,
            }
        )
        matrices[f"theta_{str(theta).replace('/', '_')}"] = matrix
        matrices[f"vector_{str(theta).replace('/', '_')}"] = np.asarray(eigen["vector"])
        print(f"theta={float(theta):.9f} score={score:.12f}", flush=True)

    full = model.face_matrices()["Jfull"]
    full_eigen = parity.spectral_generalized_maximum(full, metric, args.cutoff)
    theta_endpoint = 2 * upstream.RHO_STAR * upstream.OUTER_RADIUS
    endpoint_matrix = None
    for theta in sorted(set(args.theta)):
        if theta >= theta_endpoint:
            endpoint_matrix = matrices[f"theta_{str(theta).replace('/', '_')}"]
            break
    check = {
        "theta_endpoint": str(theta_endpoint),
        "theta_endpoint_float": float(theta_endpoint),
        "full_face_score": float(upstream.RHO_STAR) * float(full_eigen["quotient"]),
        "full_face_eigen": clean_eigen_result(full_eigen),
    }
    if endpoint_matrix is not None:
        scale = max(float(np.max(np.abs(full))), np.finfo(float).tiny)
        check["endpoint_matrix_max_abs_difference"] = float(
            np.max(np.abs(endpoint_matrix - full))
        )
        check["endpoint_matrix_relative_max_difference"] = (
            check["endpoint_matrix_max_abs_difference"] / scale
        )

    payload = {
        "schema": "primegaps.physical-parity-modulus-reach.v1",
        "status": "exploratory-float64-artificial-edge-truncation",
        "definition": {
            "operator": "sum_i || E_{i,theta} F ||^2 on the full face",
            "edge_condition": "rho_star * (outer_total_upper + retained_face_total_upper) <= theta",
            "grid_condition": "rho_star*h*(2*r+a+2*k-1) <= theta, q=r+a",
            "endpoint_convention": "conservative upper cell endpoints",
        },
        "dimension": 39,
        "intervals": args.intervals,
        "rho_star": str(upstream.RHO_STAR),
        "input_sha256": parity.file_hash(parity.DEFAULT_INPUT),
        "spectral_cutoff": args.cutoff,
        "rows": rows,
        "endpoint_check": check,
        "elapsed_seconds": monotonic() - started,
    }
    if args.matrix_output is not None:
        args.matrix_output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.matrix_output, I=metric, **matrices)
        payload["matrix_receipt"] = {
            "path": str(args.matrix_output),
            "sha256": parity.file_hash(args.matrix_output),
        }
    parity.atomic_json(args.output, payload)
    print(json.dumps(payload, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()

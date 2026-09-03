#!/usr/bin/env python3
"""Checkpointed boundary-only Arb certifier for one fixed J candidate."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
import multiprocessing
from pathlib import Path
import sys
import time

import gmpy2

from primegaps.fast_exact import arb_j, compiled_poly, fast_j


ROOT = Path(__file__).resolve().parents[1]
ORBIT_PATH = ROOT / "reproduction/240/symmetry-assembler-design/orbit_status_densities.py"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def file_hash(path):
    digest = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_payload(value):
    if isinstance(value, (tuple, list)):
        return [exact_payload(item) for item in value]
    numerator = value.numerator
    denominator = value.denominator
    if callable(numerator):
        numerator = numerator()
    if callable(denominator):
        denominator = denominator()
    return [str(int(numerator)), str(int(denominator))]


def cell_payload(cell):
    return {
        "large": cell.large,
        "shifted": cell.shifted,
        "left_large": cell.left_large,
        "right_large": cell.right_large,
        "left_legal": cell.left_legal,
        "right_legal": cell.right_legal,
        "left_limit": None if cell.left_limit is None else exact_payload(cell.left_limit),
        "right_limit": None if cell.right_limit is None else exact_payload(cell.right_limit),
        "kind": cell.kind,
        "cell": exact_payload(cell.cell),
    }


def value_hash(value):
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def rational_payload(value):
    return [str(value.numerator()), str(value.denominator())]


def load_candidate(verifier, path, dimension):
    if hasattr(verifier, "load_candidate"):
        return verifier.load_candidate(path, dimension)
    payload = json.loads(Path(path).read_text())
    degree = payload.get("degree")
    if (
        payload.get("schema") != "primegaps-stadlmann-rational-candidate-v1"
        or payload.get("k") != dimension
        or not isinstance(degree, int)
    ):
        raise ValueError("candidate schema, k, or degree mismatch")
    terms = []
    for item in payload.get("terms", ()):
        signature = tuple(int(value) for value in item["signature"])
        slack = int(item["slack_power"])
        coefficient = gmpy2.mpq(int(item["numerator"]), int(item["denominator"]))
        terms.append(verifier.Term(signature, slack, coefficient))
    if len(terms) != payload.get("basis_dimension"):
        raise ValueError("candidate term count does not match basis dimension")
    return degree, tuple(terms)


VERIFIER = None
ORBIT = None
FEATURE_GROUPS = None
PAIR_ROUTES = None
TARGETS = None
BACKEND = None
PRECISION = 0


def _cell_row(index, cell, target_densities, target_polynomials):
    started = time.perf_counter()
    contribution, statistics = arb_j.contract_cell_real_ball(
        target_polynomials,
        target_densities,
        cell,
        verifier=VERIFIER,
        polynomial_backend=BACKEND,
        precision=PRECISION,
    )
    lower, upper = contribution.endpoints()
    payload = cell_payload(cell)
    return {
        "cell_index": index,
        "cell_sha256": value_hash(payload),
        "cell": payload,
        "lower": rational_payload(lower.exact_rational()),
        "upper": rational_payload(upper.exact_rational()),
        "statistics": {
            **statistics,
            "target_polynomials": len(target_polynomials),
            "contraction_seconds": time.perf_counter() - started,
        },
        "seconds": time.perf_counter() - started,
    }


def status_worker(job):
    """Compile one density status and reuse it across all matching cells."""
    status, indexed_cells = job
    density_started = time.perf_counter()
    target_densities = fast_j.target_status_densities(
        TARGETS,
        status,
        common_dimension=DIMENSION - 1,
        delta=VERIFIER.DELTA,
        rational=gmpy2.mpq,
        orbit_size=ORBIT.monomial_symmetric_orbit_size,
    )
    density_seconds = time.perf_counter() - density_started
    active_routes = {
        pair: tuple(
            (target, structure)
            for target, structure in outputs
            if target in target_densities
        )
        for pair, outputs in PAIR_ROUTES.items()
    }
    active_routes = {pair: outputs for pair, outputs in active_routes.items() if outputs}
    cells_by_regime = {}
    for index, cell in indexed_cells:
        regime = arb_j.slice_regime_key(cell, verifier=VERIFIER)
        cells_by_regime.setdefault(regime, []).append((index, cell))
    rows = []
    polynomial_seconds = 0.0
    for regime_cells in cells_by_regime.values():
        representative = regime_cells[0][1]
        polynomial_started = time.perf_counter()
        polynomial_statistics = {}
        target_polynomials = arb_j.target_correction_polynomials(
            active_routes,
            FEATURE_GROUPS,
            representative,
            verifier=VERIFIER,
            polynomial_backend=BACKEND,
            statistics=polynomial_statistics,
        )
        current_polynomial_seconds = time.perf_counter() - polynomial_started
        polynomial_seconds += current_polynomial_seconds
        for regime_index, (index, cell) in enumerate(regime_cells):
            row = _cell_row(index, cell, target_densities, target_polynomials)
            row["statistics"].update({
                "active_pairs": len(active_routes),
                **polynomial_statistics,
                "polynomial_seconds": (
                    current_polynomial_seconds if regime_index == 0 else 0.0
                ),
                "reused_target_polynomials": regime_index != 0,
            })
            row["seconds"] += (
                current_polynomial_seconds if regime_index == 0 else 0.0
            )
            rows.append(row)
    return {
        "status": status,
        "density_seconds": density_seconds,
        "density_targets": len(target_densities),
        "active_pairs": len(active_routes),
        "slice_regimes": len(cells_by_regime),
        "polynomial_seconds": polynomial_seconds,
        "rows": rows,
    }


def read_rows(path, cells):
    rows = {}
    if not path.exists():
        return rows
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        row = json.loads(line)
        index = int(row["cell_index"])
        if index in rows or not 0 <= index < len(cells):
            raise ValueError(f"invalid/duplicate cell index on line {line_number}")
        if row["cell_sha256"] != value_hash(cell_payload(cells[index])):
            raise ValueError(f"cell hash mismatch on line {line_number}")
        rows[index] = row
    return rows


def bind_manifest(args, degree, targets, pairs, cells):
    payload = {
        "schema": "primegaps-boundary-J-arb-checkpoint-v4",
        "strategy": "boundary-cell-first-single-density-status-slice-regime-reuse",
        "k": args.k,
        "degree": degree,
        "precision_bits": args.precision,
        "stride": args.stride,
        "candidate_sha256": file_hash(args.candidate),
        "verifier_sha256": file_hash(args.verifier),
        "orbit_sha256": file_hash(ORBIT_PATH),
        "arb_j_sha256": file_hash(Path(arb_j.__file__)),
        "compiled_poly_sha256": file_hash(Path(compiled_poly.__file__)),
        "fast_j_sha256": file_hash(Path(fast_j.__file__)),
        "fast_i_sha256": file_hash(Path(fast_j.fast_i.__file__)),
        "runner_sha256": file_hash(Path(__file__)),
        "target_count": len(targets),
        "limited_targets": args.target_limit is not None,
        "pair_count": len(pairs),
        "boundary_cell_count": len(cells),
        "limited_cells": args.limit_cells is not None,
        "boundary_status_count": len({cell.density_status for cell in cells}),
        "cell_list_sha256": value_hash([cell_payload(cell) for cell in cells]),
    }
    path = args.output.with_suffix(args.output.suffix + ".manifest.json")
    if path.exists():
        if json.loads(path.read_text()) != payload:
            raise ValueError(f"checkpoint manifest mismatch: {path}")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(json.dumps(payload, sort_keys=True) + "\n")
        temporary.replace(path)


def run(args):
    global VERIFIER, ORBIT, FEATURE_GROUPS, PAIR_ROUTES, TARGETS
    global BACKEND, PRECISION, DIMENSION, MAXIMUM_OFFSET
    VERIFIER = load_module("boundary_arb_verifier", args.verifier)
    ORBIT = load_module("boundary_arb_orbit", ORBIT_PATH)
    if hasattr(VERIFIER, "configure_rational"):
        VERIFIER.configure_rational(gmpy2.mpq)
    else:
        VERIFIER.Fraction = gmpy2.mpq
        VERIFIER.kernel.Fraction = gmpy2.mpq
    degree, terms = load_candidate(VERIFIER, args.candidate, args.k)
    prep_started = time.perf_counter()
    FEATURE_GROUPS = VERIFIER.grouped_marginal_coefficients(terms)
    pair_groups = VERIFIER.grouped_signature_pairs(FEATURE_GROUPS, args.k)
    targets = tuple(pair_groups)
    if args.target_limit is not None:
        targets = targets[:args.target_limit]
    TARGETS, PAIR_ROUTES = fast_j.pair_routes(pair_groups, targets)
    DIMENSION = args.k
    MAXIMUM_OFFSET = int(VERIFIER.R // VERIFIER.DELTA)
    cells = tuple(arb_j.iter_boundary_cells(
        common_dimension=args.k - 1, verifier=VERIFIER
    ))
    if args.limit_cells is not None:
        cells = cells[:args.limit_cells]
    BACKEND = compiled_poly.ArbEncodedPolynomialBackend(
        precision=args.precision, stride=args.stride, rational=gmpy2.mpq
    )
    PRECISION = args.precision
    bind_manifest(args, degree, targets, PAIR_ROUTES, cells)
    completed = read_rows(args.output, cells)
    pending_by_status = {}
    for index, cell in enumerate(cells):
        if index not in completed:
            pending_by_status.setdefault(cell.density_status, []).append((index, cell))
    jobs = tuple(pending_by_status.items())
    pending_cell_count = sum(len(items) for _status, items in jobs)
    print(json.dumps({
        "degree": degree,
        "targets": len(targets),
        "pairs": len(PAIR_ROUTES),
        "boundary_cells": len(cells),
        "completed_cells": len(completed),
        "boundary_statuses": len({cell.density_status for cell in cells}),
        "pending_statuses": len(jobs),
        "pending_cells": pending_cell_count,
        "precision_bits": args.precision,
        "preparation_seconds": time.perf_counter() - prep_started,
    }), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with args.output.open("a") as stream, ProcessPoolExecutor(
        max_workers=args.workers, mp_context=multiprocessing.get_context("fork")
    ) as pool:
        futures = {pool.submit(status_worker, job): job[0] for job in jobs}
        new_cell_count = 0
        for status_count, future in enumerate(as_completed(futures), 1):
            result = future.result()
            for row in result["rows"]:
                stream.write(json.dumps(row, sort_keys=True) + "\n")
            stream.flush()
            new_cell_count += len(result["rows"])
            print(json.dumps({
                "completed_statuses": status_count,
                "newly_completed_cells": new_cell_count,
                "total_completed_cells": len(completed) + new_cell_count,
                "last_status": result["status"],
                "last_status_cells": len(result["rows"]),
                "last_density_seconds": result["density_seconds"],
                "last_density_targets": result["density_targets"],
                "last_active_pairs": result["active_pairs"],
                "last_cells_seconds": sum(row["seconds"] for row in result["rows"]),
                "elapsed_seconds": time.perf_counter() - started,
            }), flush=True)


def finalize(args):
    manifest_path = args.input.with_suffix(args.input.suffix + ".manifest.json")
    manifest = json.loads(manifest_path.read_text())
    lower = Fraction()
    upper = Fraction()
    indices = set()
    for line_number, line in enumerate(args.input.read_text().splitlines(), 1):
        row = json.loads(line)
        index = int(row["cell_index"])
        if index in indices:
            raise ValueError(f"duplicate cell index on line {line_number}")
        indices.add(index)
        lower += Fraction(int(row["lower"][0]), int(row["lower"][1]))
        upper += Fraction(int(row["upper"][0]), int(row["upper"][1]))
    expected_count = int(manifest["boundary_cell_count"])
    checkpoint_complete = indices == set(range(expected_count))
    complete = checkpoint_complete
    certification_eligible = (
        checkpoint_complete
        and manifest.get("limited_cells") is False
        and manifest.get("limited_targets") is False
    )
    result = {
        "schema": "primegaps-boundary-J-arb-result-v1",
        "cell_count": len(indices),
        "expected_cell_count": expected_count,
        "checkpoint_complete": checkpoint_complete,
        "complete": complete,
        "certification_eligible": certification_eligible,
        "correction_lower": [str(lower.numerator), str(lower.denominator)],
        "correction_upper": [str(upper.numerator), str(upper.denominator)],
        "correction_lower_decimal": f"{float(lower):.17g}",
        "correction_upper_decimal": f"{float(upper):.17g}",
        "interval_width_decimal": f"{float(upper - lower):.17g}",
        "input_sha256": file_hash(args.input),
        "manifest_sha256": file_hash(manifest_path),
    }
    unrestricted = None
    if args.unrestricted_j is not None:
        unrestricted = json.loads(args.unrestricted_j.read_text())
    if args.i_upper is not None:
        if unrestricted is None:
            raise ValueError("--i-upper requires --unrestricted-j")
        i_upper = json.loads(args.i_upper.read_text())
        if (
            unrestricted["candidate_sha256"] != manifest["candidate_sha256"]
            or i_upper["candidate_sha256"] != manifest["candidate_sha256"]
            or int(unrestricted["k"]) != int(manifest["k"])
            or int(i_upper["k"]) != int(manifest["k"])
        ):
            raise ValueError("candidate or k mismatch in exact gate inputs")
        normalized_unrestricted = Fraction(
            int(unrestricted["normalized_kJ_numerator"]),
            int(unrestricted["normalized_kJ_denominator"]),
        )
        normalized_i_upper = Fraction(
            int(i_upper["normalized_I_numerator"]),
            int(i_upper["normalized_I_denominator"]),
        )
        radius_u = Fraction(unrestricted["radius_U"])
        normalization = (
            int(manifest["k"]) * __import__("math").factorial(int(manifest["k"]))
            / radius_u ** int(manifest["k"])
        )
        normalized_lower = lower * normalization
        normalized_upper = upper * normalization
        required = normalized_i_upper - normalized_unrestricted
        legal_lower = normalized_unrestricted + normalized_lower
        legal_upper = normalized_unrestricted + normalized_upper
        result["exact_gate"] = {
            "unrestricted_j_sha256": file_hash(args.unrestricted_j),
            "i_upper_sha256": file_hash(args.i_upper),
            "normalized_correction_lower": exact_payload(normalized_lower),
            "normalized_correction_upper": exact_payload(normalized_upper),
            "normalized_correction_lower_decimal": f"{float(normalized_lower):.17g}",
            "normalized_correction_upper_decimal": f"{float(normalized_upper):.17g}",
            "required_normalized_correction": exact_payload(required),
            "required_normalized_correction_decimal": f"{float(required):.17g}",
            "normalized_legal_kJ_lower": exact_payload(legal_lower),
            "normalized_legal_kJ_upper": exact_payload(legal_upper),
            "normalized_i_upper": exact_payload(normalized_i_upper),
            "certified_strictly_above_one": (
                certification_eligible and legal_lower > normalized_i_upper
            ),
            "next_precision_bits_if_inconclusive": (
                None
                if certification_eligible and legal_lower > normalized_i_upper
                else 2 * int(manifest["precision_bits"])
            ),
        }
    if args.oracle_legal_result is not None:
        if args.unrestricted_j is None:
            raise ValueError("--oracle-legal-result requires --unrestricted-j")
        oracle = json.loads(args.oracle_legal_result.read_text())
        legal = Fraction(int(oracle["J"]["numerator"]), int(oracle["J"]["denominator"]))
        unrestricted_raw = Fraction(
            int(unrestricted["J_numerator"]), int(unrestricted["J_denominator"])
        )
        exact_correction = legal - unrestricted_raw
        result["oracle_check"] = {
            "oracle_sha256": file_hash(args.oracle_legal_result),
            "exact_correction": exact_payload(exact_correction),
            "interval_contains_exact": lower <= exact_correction <= upper,
        }
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.output)
    print(json.dumps(result, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    build = subparsers.add_parser("run")
    build.add_argument("--k", type=int, required=True)
    build.add_argument("--candidate", type=Path, required=True)
    build.add_argument("--verifier", type=Path, required=True)
    build.add_argument("--output", type=Path, required=True)
    build.add_argument("--workers", type=int, default=16)
    build.add_argument("--precision", type=int, default=128)
    build.add_argument("--stride", type=int, default=384)
    build.add_argument("--limit-cells", type=int)
    build.add_argument("--target-limit", type=int)
    finish = subparsers.add_parser("finalize")
    finish.add_argument("--input", type=Path, required=True)
    finish.add_argument("--output", type=Path, required=True)
    finish.add_argument("--unrestricted-j", type=Path)
    finish.add_argument("--i-upper", type=Path)
    finish.add_argument("--oracle-legal-result", type=Path)
    args = parser.parse_args()
    if args.command == "run":
        run(args)
    else:
        finalize(args)


if __name__ == "__main__":
    main()

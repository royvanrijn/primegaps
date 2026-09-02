#!/usr/bin/env python3
"""Parallel, checkpointed exact contraction for one fixed rational candidate."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import multiprocessing
from pathlib import Path
import sys
import time
import importlib.util

import gmpy2

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exact_symmetric_verifier as verifier
import exact_provenance as provenance

DP_PATH = HERE.parent / "symmetry-assembler-design" / "orbit_status_densities.py"
dp_spec = importlib.util.spec_from_file_location("orbit_status_dp", DP_PATH)
orbit_status_dp = importlib.util.module_from_spec(dp_spec)
assert dp_spec.loader is not None
dp_spec.loader.exec_module(orbit_status_dp)


verifier.Fraction = gmpy2.mpq
verifier.kernel.Fraction = gmpy2.mpq

FEATURE_GROUPS = None
PAIR_GROUPS = None
I_GROUPS = None
DIMENSION = 0


def _j_worker(signature):
    started = time.time()
    value = verifier.exact_j_signature_group(
        signature, PAIR_GROUPS[signature], FEATURE_GROUPS, DIMENSION
    )
    verifier.orbit_status_densities.cache_clear()
    return {
        "signature": list(signature),
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "seconds": time.time() - started,
    }


def _i_worker(signature):
    started = time.time()
    value = gmpy2.mpq(0)
    density = orbit_status_dp.orbit_status_densities(
        signature,
        k=DIMENSION,
        delta=verifier.DELTA,
        max_large=len(verifier.B),
        max_offset_count=int(verifier.U // verifier.DELTA),
    )
    orbit = orbit_status_dp.monomial_symmetric_orbit_size(signature, DIMENSION)
    for (large_count, shifted_count, large_power, small_power), density_coefficient in (
        orbit_status_dp.normalized_density_terms(density, DIMENSION)
    ):
        height = verifier.U - (large_count + shifted_count) * verifier.DELTA
        large_cap = (
            None if large_count == 0
            else verifier.B[large_count - 1] - large_count * verifier.DELTA
        )
        for slack, candidate_coefficient in I_GROUPS[signature].items():
            value += (
                orbit * density_coefficient * candidate_coefficient
                * verifier._radial_group_integral(
                    large_power, small_power, slack, height, large_cap
                )
            )
    return {
        "signature": list(signature),
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("I", "J"), required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--binding", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    global FEATURE_GROUPS, PAIR_GROUPS, I_GROUPS, DIMENSION
    DIMENSION = args.k
    terms = verifier.rational_terms_from_candidate(args.candidate, args.k)
    if args.kind == "J":
        FEATURE_GROUPS = verifier.grouped_marginal_coefficients(terms)
        PAIR_GROUPS = verifier.grouped_signature_pairs(FEATURE_GROUPS, args.k)
        signatures = list(PAIR_GROUPS)
        worker = _j_worker
    else:
        atoms = verifier.aggregate_i_atoms(terms, args.k)
        I_GROUPS = {}
        for (signature, slack), coefficient in atoms.items():
            I_GROUPS.setdefault(signature, {})[slack] = coefficient
        signatures = list(I_GROUPS)
        worker = _i_worker

    evaluator_files = (
        Path(verifier.__file__),
        Path(verifier.kernel.__file__),
        DP_PATH,
        Path(__file__),
        Path(provenance.__file__),
    )
    context, context_hash = provenance.build_context(
        kind=args.kind,
        k=args.k,
        terms=terms,
        signatures=signatures,
        candidate=args.candidate,
        binding=args.binding,
        evaluator_files=evaluator_files,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = provenance.manifest_path(args.output)
    if args.output.exists() and not manifest_path.exists():
        raise ValueError("refusing to resume an unbound legacy checkpoint")
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("context") != context or manifest.get("context_sha256") != context_hash:
            raise ValueError("checkpoint provenance does not match this exact run")
        if manifest.get("completed"):
            if manifest.get("records_sha256") != provenance.file_hash(args.output):
                raise ValueError("completed checkpoint record hash does not match")
            summary_path = args.output.with_suffix(".summary.json")
            if (
                not summary_path.exists()
                or manifest.get("summary_sha256") != provenance.file_hash(summary_path)
            ):
                raise ValueError("completed checkpoint summary hash does not match")
    else:
        manifest = {
            "schema": "primegaps-stadlmann-exact-group-manifest-v1",
            "context": context,
            "context_sha256": context_hash,
            "completed": False,
            "records_sha256": None,
        }
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    completed = {}
    expected_signatures = set(signatures)
    if args.output.exists():
        for line_number, line in enumerate(args.output.read_text().splitlines(), 1):
            record = json.loads(line)
            signature = tuple(record["signature"])
            if record.get("schema") != provenance.RECORD_SCHEMA:
                raise ValueError(f"unbound checkpoint record at line {line_number}")
            if record.get("context_sha256") != context_hash:
                raise ValueError(f"mixed checkpoint provenance at line {line_number}")
            if signature not in expected_signatures:
                raise ValueError(f"unexpected checkpoint signature at line {line_number}")
            if signature in completed:
                raise ValueError(f"duplicate checkpoint signature at line {line_number}")
            completed[signature] = record
    pending = [signature for signature in signatures if signature not in completed]
    print(
        json.dumps({
            "kind": args.kind, "k": args.k, "terms": len(terms),
            "groups": len(signatures), "completed": len(completed),
            "pending": len(pending), "workers": args.workers,
        }),
        flush=True,
    )
    started = time.time()
    with ProcessPoolExecutor(
        max_workers=args.workers, mp_context=multiprocessing.get_context("fork")
    ) as executor:
        futures = {executor.submit(worker, signature): signature for signature in pending}
        with args.output.open("a") as stream:
            for index, future in enumerate(as_completed(futures), 1):
                record = future.result()
                record.update({
                    "schema": provenance.RECORD_SCHEMA,
                    "kind": args.kind,
                    "k": args.k,
                    "context_sha256": context_hash,
                })
                stream.write(json.dumps(record, sort_keys=True) + "\n")
                stream.flush()
                completed[tuple(record["signature"])] = record
                if index % 25 == 0 or index == len(pending):
                    print(
                        json.dumps({
                            "newly_completed": index,
                            "total_completed": len(completed),
                            "groups": len(signatures),
                            "elapsed_seconds": time.time() - started,
                        }),
                        flush=True,
                    )
    total = gmpy2.mpq(0)
    for record in completed.values():
        total += gmpy2.mpq(int(record["numerator"]), int(record["denominator"]))
    summary = {
        "kind": args.kind,
        "k": args.k,
        "degree": 21,
        "term_count": len(terms),
        "group_count": len(signatures),
        "candidate_sha256": context["candidate_sha256"],
        "context_sha256": context_hash,
        "numerator": str(total.numerator),
        "denominator": str(total.denominator),
        "seconds_this_run": time.time() - started,
    }
    summary_path = args.output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    manifest.update({
        "completed": True,
        "records_sha256": provenance.file_hash(args.output),
        "summary_sha256": provenance.file_hash(summary_path),
    })
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({**summary, "numerator": str(total.numerator)[:40] + "...",
                      "denominator": str(total.denominator)[:40] + "..."}, indent=2), flush=True)


if __name__ == "__main__":
    main()

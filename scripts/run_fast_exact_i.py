#!/usr/bin/env python3
"""Checkpointed benchmark runner for the closed-zero exact I backend."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
import importlib.util
import json
import multiprocessing
from pathlib import Path
import sys
import time

import gmpy2


HERE = Path(__file__).resolve().parent
ROOT = next(
    parent for parent in (HERE, *HERE.parents)
    if (parent / "pyproject.toml").is_file()
)
FROZEN = ROOT / "reproduction" / "240" / "independent-reproducer"
DP_PATH = (
    ROOT / "reproduction" / "240" / "symmetry-assembler-design"
    / "orbit_status_densities.py"
)
sys.path.insert(0, str(FROZEN))

from primegaps.fast_exact import fast_i, moment_cache
import exact_symmetric_verifier as verifier

dp_spec = importlib.util.spec_from_file_location("fast_i_orbit_helpers", DP_PATH)
orbit_helpers = importlib.util.module_from_spec(dp_spec)
assert dp_spec.loader is not None
dp_spec.loader.exec_module(orbit_helpers)

verifier.Fraction = gmpy2.mpq
verifier.kernel.Fraction = gmpy2.mpq


DIMENSION = 0
GROUPS = None
POSITIVE_CACHE = {}
MOMENT_VALUES = None


def file_hash(path):
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def bind_manifest(args):
    payload = {
        "schema": "primegaps-fast-exact-I-checkpoint-v1",
        "k": args.k,
        "candidate_sha256": file_hash(args.candidate),
        "verifier_sha256": file_hash(FROZEN / "exact_symmetric_verifier.py"),
        "orbit_helpers_sha256": file_hash(DP_PATH),
        "fast_i_sha256": file_hash(Path(fast_i.__file__)),
        "runner_sha256": file_hash(Path(__file__)),
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
    return path


def worker(signature):
    started = time.perf_counter()
    fresh_moments = {}
    if MOMENT_VALUES is None:
        value = fast_i.signature_value(
            signature,
            GROUPS[signature],
            k=DIMENSION,
            delta=verifier.DELTA,
            total_cap=verifier.U,
            large_caps=verifier.B,
            rational=gmpy2.mpq,
            orbit_size=orbit_helpers.monomial_symmetric_orbit_size,
            radial_moment=verifier._radial_group_integral,
            positive_cache=POSITIVE_CACHE,
        )
    else:
        existing = MOMENT_VALUES.get(signature, {})
        missing = tuple(sorted(set(GROUPS[signature]) - set(existing)))
        if missing:
            fresh_moments = fast_i.signature_moments(
                signature,
                missing,
                k=DIMENSION,
                delta=verifier.DELTA,
                total_cap=verifier.U,
                large_caps=verifier.B,
                rational=gmpy2.mpq,
                orbit_size=orbit_helpers.monomial_symmetric_orbit_size,
                radial_moment=verifier._radial_group_integral,
                positive_cache=POSITIVE_CACHE,
            )
        all_moments = dict(existing)
        all_moments.update(fresh_moments)
        value = sum(
            (
                coefficient * all_moments[slack]
                for slack, coefficient in GROUPS[signature].items()
            ),
            gmpy2.mpq(0),
        )
    row = {
        "signature": list(signature),
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "seconds": time.perf_counter() - started,
    }
    if fresh_moments:
        row["_fresh_moments"] = fresh_moments
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--limit", type=int)
    parser.add_argument(
        "--moment-cache", type=Path,
        help="append/reuse candidate-independent exact I moments",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    bind_manifest(args)

    global DIMENSION, GROUPS, MOMENT_VALUES
    DIMENSION = args.k
    terms = verifier.rational_terms_from_candidate(args.candidate, args.k)
    atoms = verifier.aggregate_i_atoms(terms, args.k)
    GROUPS = {}
    for (signature, slack), coefficient in atoms.items():
        GROUPS.setdefault(signature, {})[slack] = coefficient
    cache = None
    if args.moment_cache is not None:
        context = {
            "schema": "primegaps-fast-exact-I-context-v1",
            "k": args.k,
            "delta": str(verifier.DELTA),
            "total_cap": str(verifier.U),
            "large_caps": [str(value) for value in verifier.B],
            "verifier_sha256": file_hash(FROZEN / "exact_symmetric_verifier.py"),
        }
        cache = moment_cache.IMomentCache(
            args.moment_cache, context=context, rational=gmpy2.mpq
        )
        MOMENT_VALUES = cache.values
    signatures = tuple(GROUPS)
    if args.limit is not None:
        signatures = signatures[: args.limit]

    completed = {}
    if args.output.exists():
        for line_number, line in enumerate(args.output.read_text().splitlines(), 1):
            row = json.loads(line)
            signature = tuple(row["signature"])
            if signature in completed:
                raise ValueError(f"duplicate signature on line {line_number}")
            completed[signature] = row
    unexpected = set(completed) - set(signatures)
    if unexpected:
        raise ValueError(f"unexpected checkpoint signatures: {sorted(unexpected)[:3]}")
    pending = [signature for signature in signatures if signature not in completed]
    print(json.dumps({
        "k": args.k,
        "groups": len(signatures),
        "completed": len(completed),
        "pending": len(pending),
        "workers": args.workers,
    }), flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with args.output.open("a") as stream, ProcessPoolExecutor(
        max_workers=args.workers,
        mp_context=multiprocessing.get_context("fork"),
    ) as pool:
        futures = {pool.submit(worker, signature): signature for signature in pending}
        for count, future in enumerate(as_completed(futures), 1):
            row = future.result()
            fresh_moments = row.pop("_fresh_moments", {})
            if cache is not None and fresh_moments:
                cache.append(tuple(row["signature"]), fresh_moments)
            stream.write(json.dumps(row, sort_keys=True) + "\n")
            stream.flush()
            if count % 100 == 0 or count == len(pending):
                print(json.dumps({
                    "newly_completed": count,
                    "total_completed": len(completed) + count,
                    "elapsed_seconds": time.perf_counter() - started,
                }), flush=True)


if __name__ == "__main__":
    main()

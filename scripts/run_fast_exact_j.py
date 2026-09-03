#!/usr/bin/env python3
"""Checkpointed pair-first exact J benchmark runner."""

from __future__ import annotations

import argparse
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
import importlib.util
import json
import multiprocessing
import os
from pathlib import Path
import sys
import tempfile
import time

import gmpy2


HERE = Path(__file__).resolve().parent
ROOT = next(
    parent for parent in (HERE, *HERE.parents)
    if (parent / "pyproject.toml").is_file()
)
FROZEN = ROOT / "reproduction" / "240" / "independent-reproducer"
DEFAULT_VERIFIER_PATH = FROZEN / "exact_symmetric_verifier.py"
DP_PATH = (
    ROOT / "reproduction" / "240" / "symmetry-assembler-design"
    / "orbit_status_densities.py"
)
sys.path.insert(0, str(FROZEN))

from primegaps.fast_exact import (
    compiled_poly,
    fast_i,
    fast_j,
    modular_exact,
    moment_cache,
)
import exact_symmetric_verifier as verifier

dp_spec = importlib.util.spec_from_file_location("fast_j_orbit_helpers", DP_PATH)
orbit_helpers = importlib.util.module_from_spec(dp_spec)
assert dp_spec.loader is not None
dp_spec.loader.exec_module(orbit_helpers)

verifier.Fraction = gmpy2.mpq
verifier.kernel.Fraction = gmpy2.mpq


DIMENSION = 0
FEATURE_GROUPS = None
PAIR_GROUPS = None
POSITIVE_CACHE = {}
POLYNOMIAL_BACKEND = None
MODULAR_PRIMES = ()
FUNCTIONAL_VALUES = None
DENSITY_STATUSES = None
CACHE_CONTEXT_HASH = None
CACHE_SPOOL_DIRECTORY = None


class BoundedLRUCache:
    def __init__(self, maximum_size):
        self.maximum_size = int(maximum_size)
        self.values = OrderedDict()

    def get(self, key):
        value = self.values.get(key)
        if value is not None:
            self.values.move_to_end(key)
        return value

    def __setitem__(self, key, value):
        self.values[key] = value
        self.values.move_to_end(key)
        if len(self.values) > self.maximum_size:
            self.values.popitem(last=False)


SLICE_CACHE = BoundedLRUCache(10_000)


def file_hash(path):
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_candidate(path, dimension):
    if hasattr(verifier, "load_candidate"):
        return verifier.load_candidate(path, dimension)
    payload = json.loads(path.read_text())
    degree = payload.get("degree")
    if (
        payload.get("schema") != "primegaps-stadlmann-rational-candidate-v1"
        or payload.get("k") != dimension
        or not isinstance(degree, int)
        or degree < 0
        or payload.get("radial_coordinate") != "(U-sum(t)) with U=521/2000"
    ):
        raise ValueError("candidate schema, k, degree, or convention mismatch")
    combined = {}
    for item in payload.get("terms", []):
        raw_signature = item.get("signature")
        slack = item.get("slack_power")
        if not isinstance(raw_signature, list):
            raise ValueError("candidate signature must be a JSON list")
        signature = tuple(raw_signature)
        if (
            signature != tuple(sorted(signature, reverse=True))
            or any(
                isinstance(value, bool)
                or not isinstance(value, int)
                or value <= 0
                or value % 2
                for value in signature
            )
            or isinstance(slack, bool)
            or not isinstance(slack, int)
            or slack < 0
            or sum(signature) + slack > degree
        ):
            raise ValueError("invalid candidate basis term")
        key = (signature, slack)
        if key in combined:
            raise ValueError("duplicate candidate basis term")
        numerator = item.get("numerator")
        denominator = item.get("denominator")
        if (
            isinstance(numerator, bool)
            or not isinstance(numerator, int)
            or isinstance(denominator, bool)
            or not isinstance(denominator, int)
            or denominator <= 0
        ):
            raise ValueError("invalid candidate rational coefficient")
        value = gmpy2.mpq(numerator, denominator)
        if not value:
            raise ValueError("zero candidate term")
        combined[key] = value
    if len(combined) != payload.get("basis_dimension"):
        raise ValueError("candidate term count/basis dimension mismatch")
    return degree, tuple(
        verifier.Term(signature, slack, coefficient)
        for (signature, slack), coefficient in sorted(combined.items())
    )


def bind_manifest(args, degree):
    backend = (
        "modular:" + ",".join(map(str, args.modular_prime))
        if args.modular_prime
        else ("flint-rational" if args.compiled else "python-dict")
    )
    payload = {
        "schema": "primegaps-fast-exact-J-checkpoint-v2",
        "k": args.k,
        "degree": degree,
        "backend": backend,
        "candidate_sha256": file_hash(args.candidate),
        "verifier_sha256": file_hash(args.verifier or DEFAULT_VERIFIER_PATH),
        "orbit_helpers_sha256": file_hash(DP_PATH),
        "fast_i_sha256": file_hash(Path(fast_i.__file__)),
        "fast_j_sha256": file_hash(Path(fast_j.__file__)),
        "compiled_poly_sha256": file_hash(Path(compiled_poly.__file__)),
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


def worker(targets):
    started = time.perf_counter()
    fresh_functionals = {}
    fresh_statuses = {}
    spool_path = None
    if FUNCTIONAL_VALUES is None:
        values = fast_j.evaluate_target_chunk(
            targets,
            dimension=DIMENSION,
            feature_groups=FEATURE_GROUPS,
            pair_groups=PAIR_GROUPS,
            verifier=verifier,
            orbit_size=orbit_helpers.monomial_symmetric_orbit_size,
            rational=gmpy2.mpq,
            positive_cache=POSITIVE_CACHE,
            polynomial_backend=POLYNOMIAL_BACKEND,
            slice_cache=SLICE_CACHE,
        )
    else:
        descriptor, raw_spool_path = tempfile.mkstemp(
            prefix=".j-functional-worker-",
            suffix=".jsonl",
            dir=CACHE_SPOOL_DIRECTORY,
            text=True,
        )
        spool_path = Path(raw_spool_path)
        try:
            with os.fdopen(descriptor, "w") as spool:
                def write_status(target, statuses):
                    record = moment_cache.j_density_status_record(
                        CACHE_CONTEXT_HASH, target, statuses
                    )
                    spool.write(json.dumps(record, sort_keys=True) + "\n")

                def write_functional(functional_id, moments):
                    record = moment_cache.j_functional_record(
                        CACHE_CONTEXT_HASH, functional_id, moments
                    )
                    spool.write(json.dumps(record, sort_keys=True) + "\n")

                values, fresh_functionals, fresh_statuses = (
                    fast_j.evaluate_target_chunk_cached(
                        targets,
                        dimension=DIMENSION,
                        feature_groups=FEATURE_GROUPS,
                        pair_groups=PAIR_GROUPS,
                        verifier=verifier,
                        orbit_size=orbit_helpers.monomial_symmetric_orbit_size,
                        rational=gmpy2.mpq,
                        functional_values=FUNCTIONAL_VALUES,
                        density_statuses=DENSITY_STATUSES,
                        positive_cache=POSITIVE_CACHE,
                        polynomial_backend=POLYNOMIAL_BACKEND,
                        slice_cache=SLICE_CACHE,
                        fresh_functional_sink=write_functional,
                        fresh_status_sink=write_status,
                    )
                )
        except BaseException:
            spool_path.unlink(missing_ok=True)
            raise
        if not spool_path.stat().st_size:
            spool_path.unlink()
            spool_path = None
    elapsed = time.perf_counter() - started
    per_target = elapsed / len(targets)
    rows = []
    for signature, value in values.items():
        row = {
            "signature": list(signature),
            "chunk_seconds_per_target": per_target,
        }
        if not MODULAR_PRIMES:
            row.update({
                "numerator": str(value.numerator),
                "denominator": str(value.denominator),
            })
        elif len(MODULAR_PRIMES) == 1:
            row.update({"prime": MODULAR_PRIMES[0], "residue": int(value)})
        else:
            row.update({
                "primes": list(MODULAR_PRIMES),
                "residues": [int(residue) for residue in value],
            })
        rows.append(row)
    return rows, fresh_functionals, fresh_statuses, spool_path


def main():
    global verifier
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--chunk-size", type=int, default=32)
    parser.add_argument("--limit", type=int)
    parser.add_argument(
        "--verifier", type=Path,
        help="exact support adapter; defaults to the frozen 240 verifier",
    )
    parser.add_argument("--compiled", action="store_true")
    parser.add_argument(
        "--modular-prime", type=int, action="append",
        help="repeat to batch several CRT primes in one geometry pass",
    )
    parser.add_argument(
        "--moment-cache", type=Path,
        help="append/reuse candidate-independent exact J functionals",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.verifier is not None:
        verifier = load_module("fast_j_support_adapter", args.verifier)
        if hasattr(verifier, "configure_rational"):
            verifier.configure_rational(gmpy2.mpq)
        else:
            verifier.Fraction = gmpy2.mpq
            verifier.kernel.Fraction = gmpy2.mpq
    if args.chunk_size < 1 or args.workers < 1:
        raise ValueError("workers and chunk size must be positive")
    if (
        args.modular_prime
        and not all(modular_exact.is_prime_64(p) for p in args.modular_prime)
    ):
        raise ValueError("--modular-prime must be prime")
    if args.compiled and args.modular_prime:
        raise ValueError("choose rational compiled or modular, not both")
    if args.moment_cache is not None and args.modular_prime:
        raise ValueError("the exact moment cache cannot be used modulo a prime")
    if args.modular_prime and len(set(args.modular_prime)) != len(
        args.modular_prime
    ):
        raise ValueError("modular primes must be distinct")

    global DIMENSION, FEATURE_GROUPS, PAIR_GROUPS, POLYNOMIAL_BACKEND
    global MODULAR_PRIMES, FUNCTIONAL_VALUES, DENSITY_STATUSES
    global CACHE_CONTEXT_HASH, CACHE_SPOOL_DIRECTORY
    DIMENSION = args.k
    degree, terms = load_candidate(args.candidate, args.k)
    bind_manifest(args, degree)
    FEATURE_GROUPS = verifier.grouped_marginal_coefficients(terms)
    PAIR_GROUPS = verifier.grouped_signature_pairs(FEATURE_GROUPS, args.k)
    if args.compiled:
        POLYNOMIAL_BACKEND = compiled_poly.FlintEncodedPolynomialBackend(
            stride=256, rational=gmpy2.mpq
        )
    elif args.modular_prime:
        MODULAR_PRIMES = tuple(args.modular_prime)
        backends = tuple(
            compiled_poly.FlintModularEncodedPolynomialBackend(prime, stride=256)
            for prime in MODULAR_PRIMES
        )
        POLYNOMIAL_BACKEND = (
            backends[0]
            if len(backends) == 1
            else compiled_poly.ProductPolynomialBackend(backends)
        )
    cache = None
    if args.moment_cache is not None:
        context = {
            "schema": "primegaps-fast-exact-J-context-v1",
            "functional_key_schema": "target-status-cell-v1",
            "k": args.k,
            "delta": str(verifier.DELTA),
            "residual_cap": str(verifier.R),
            "total_cap": str(verifier.U),
            "large_caps": [str(value) for value in verifier.B],
            "verifier_sha256": file_hash(
                args.verifier or DEFAULT_VERIFIER_PATH
            ),
            "orbit_helpers_sha256": file_hash(DP_PATH),
            "fast_i_sha256": file_hash(Path(fast_i.__file__)),
            "fast_j_sha256": file_hash(Path(fast_j.__file__)),
        }
        cache = moment_cache.JFunctionalCache(
            args.moment_cache, context=context, rational=gmpy2.mpq
        )
        FUNCTIONAL_VALUES = cache.values
        DENSITY_STATUSES = cache.density_statuses
        CACHE_CONTEXT_HASH = cache.context_hash
        CACHE_SPOOL_DIRECTORY = str(cache.path.parent)
    signatures = list(PAIR_GROUPS)
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
    expected = set(signatures)
    unexpected = set(completed) - expected
    if unexpected:
        raise ValueError(f"unexpected checkpoint signatures: {sorted(unexpected)[:3]}")
    pending = [signature for signature in signatures if signature not in completed]
    chunks = [
        tuple(pending[index : index + args.chunk_size])
        for index in range(0, len(pending), args.chunk_size)
    ]
    print(json.dumps({
        "k": args.k,
        "degree": degree,
        "groups": len(signatures),
        "completed": len(completed),
        "pending": len(pending),
        "chunks": len(chunks),
        "chunk_size": args.chunk_size,
        "workers": args.workers,
        "compiled": args.compiled,
        "modular_prime": args.modular_prime,
        "moment_cache": str(args.moment_cache) if args.moment_cache else None,
    }), flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with args.output.open("a") as stream, ProcessPoolExecutor(
        max_workers=args.workers,
        mp_context=multiprocessing.get_context("fork"),
    ) as pool:
        futures = {pool.submit(worker, chunk): chunk for chunk in chunks}
        newly_completed = 0
        for future in as_completed(futures):
            rows, fresh_functionals, fresh_statuses, spool_path = future.result()
            if cache is not None:
                if spool_path is not None:
                    cache.ingest(spool_path, retain=False)
                    Path(spool_path).unlink()
                else:
                    for target, statuses in fresh_statuses.items():
                        cache.append_density_statuses(target, statuses)
                    for functional_id, moments in fresh_functionals.items():
                        cache.append(functional_id, moments)
            for row in rows:
                stream.write(json.dumps(row, sort_keys=True) + "\n")
                completed[tuple(row["signature"])] = row
            stream.flush()
            newly_completed += len(rows)
            print(json.dumps({
                "newly_completed": newly_completed,
                "total_completed": len(completed),
                "elapsed_seconds": time.perf_counter() - started,
            }), flush=True)


if __name__ == "__main__":
    main()

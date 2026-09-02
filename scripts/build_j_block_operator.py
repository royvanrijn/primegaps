#!/usr/bin/env python3
"""Experimentally compile exact cached moments into floating J blocks.

Production-degree monomial moments are extremely ill-conditioned when converted
to float64 or long double.  This tool is retained for low-degree checks and
conditioning experiments; discovery operators should instead be assembled from
directly evaluated orthogonal features.
"""

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
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = next(parent for parent in (HERE, *HERE.parents) if (parent / "pyproject.toml").is_file())
RUNNER_PATH = HERE / "run_fast_exact_j.py"
spec = importlib.util.spec_from_file_location("j_block_source_runner", RUNNER_PATH)
runner = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = runner
spec.loader.exec_module(runner)

from primegaps.fast_exact import j_block, moment_cache
from primegaps.basis import symmetric_basis


COMPILE_STATE = None


def compile_worker(pair):
    left, right = pair
    state = COMPILE_STATE
    return pair, j_block.compile_signature_pair_block(
        left,
        right,
        left_keys=state["marginal_map"].feature_keys[left],
        right_keys=state["marginal_map"].feature_keys[right],
        pair_groups=state["pair_groups"],
        functional_values=state["functional_values"],
        density_statuses=state["density_statuses"],
        common_dimension=state["common_dimension"],
        verifier=runner.verifier,
        rational=gmpy2.mpq,
        dtype=state["dtype"],
        control_variate=state["control_variate"],
        radial_basis=state["radial_basis"],
        route_index=state["route_index"],
    )


def value_hash(value):
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def block_name(left, right):
    digest = value_hash([list(left), list(right)])
    return f"blocks/{digest}.npy"


def atomic_save_array(path, array):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        np.save(stream, array, allow_pickle=False)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--target-cache", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--dtype", choices=("float64", "longdouble"), default="float64")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--control-variate", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--radial-basis", choices=("monomial", "jacobi"), default="monomial"
    )
    parser.add_argument(
        "--accept-ill-conditioned-cache-conversion",
        action="store_true",
        help="acknowledge that production-degree floating blocks require validation",
    )
    args = parser.parse_args()
    if args.workers < 1:
        raise ValueError("workers must be positive")
    if not args.accept_ill_conditioned_cache_conversion:
        parser.error(
            "cached exact moments are not numerically safe at production degree; "
            "pass --accept-ill-conditioned-cache-conversion only for experiments"
        )

    header = json.loads(args.target_cache.read_text().splitlines()[0])
    context = header.get("context")
    if not isinstance(context, dict) or context.get("k") != args.k:
        raise ValueError("target cache context/k mismatch")
    cache = moment_cache.JFunctionalCache(
        args.target_cache, context=context, rational=gmpy2.mpq
    )
    degree, terms = runner.load_candidate(args.candidate, args.k)
    feature_groups = runner.verifier.grouped_marginal_coefficients(terms)
    pair_groups = runner.verifier.grouped_signature_pairs(feature_groups, args.k)
    route_index = j_block.signature_pair_route_index(pair_groups)
    operator_basis = (
        tuple((term.signature, term.slack) for term in terms)
        if args.radial_basis == "monomial"
        else tuple(
            (tuple(2 * value for value in signature), radial_degree)
            for signature, radial_degree in symmetric_basis(degree, args.k)
        )
    )
    marginal_map = j_block.MarginalMap.from_basis(operator_basis)
    signatures = tuple(marginal_map.feature_keys)
    pairs = tuple(
        (left, right)
        for left_index, left in enumerate(signatures)
        for right in signatures[left_index:]
        if (left, right) in route_index
    )
    dtype = np.float64 if args.dtype == "float64" else np.longdouble
    manifest = {
        "schema": "primegaps-J-block-operator-v1",
        "k": args.k,
        "degree": degree,
        "dtype": args.dtype,
        "control_variate": args.control_variate,
        "radial_basis": args.radial_basis,
        "numerically_validated": False,
        "warning": "floating conversion of exact monomial moments is ill-conditioned",
        "basis": [
            {"signature": list(signature), "slack": slack}
            for signature, slack in marginal_map.basis
        ],
        "expected_block_count": len(pairs),
        "target_cache_sha256": j_block.file_sha256(args.target_cache),
        "target_cache_context_sha256": cache.context_hash,
        "candidate_basis_sha256": value_hash([
            [list(signature), slack] for signature, slack in marginal_map.basis
        ]),
        "compiler_sha256": j_block.file_sha256(Path(__file__)),
        "j_block_sha256": j_block.file_sha256(Path(j_block.__file__)),
    }
    args.output_directory.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_directory / "manifest.json"
    if manifest_path.exists():
        if json.loads(manifest_path.read_text()) != manifest:
            raise ValueError("J block manifest mismatch")
    else:
        temporary = manifest_path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(manifest, sort_keys=True) + "\n")
        temporary.replace(manifest_path)

    index_path = args.output_directory / "index.jsonl"
    completed = {}
    if index_path.exists():
        for line_number, line in enumerate(index_path.read_text().splitlines(), 1):
            record = json.loads(line)
            key = (tuple(record["left_signature"]), tuple(record["right_signature"]))
            if key in completed:
                raise ValueError(f"duplicate block index row {line_number}")
            block_path = args.output_directory / record["path"]
            if j_block.file_sha256(block_path) != record["sha256"]:
                raise ValueError(f"block hash mismatch: {block_path}")
            completed[key] = record
    pending = [pair for pair in pairs if pair not in completed]
    if args.limit is not None:
        pending = pending[:args.limit]
    required_targets = {
        target
        for left, right in pending
        for target, _structure in route_index[(left, right)]
    }
    missing_statuses = required_targets - set(cache.density_statuses)
    if missing_statuses:
        raise ValueError(
            f"target cache lacks {len(missing_statuses)} required density-status indexes"
        )
    print(json.dumps({
        "blocks": len(pairs),
        "completed": len(completed),
        "pending_this_run": len(pending),
        "dtype": args.dtype,
    }), flush=True)
    started = time.perf_counter()
    global COMPILE_STATE
    COMPILE_STATE = {
        "marginal_map": marginal_map,
        "pair_groups": pair_groups,
        "route_index": route_index,
        "functional_values": cache.values,
        "density_statuses": cache.density_statuses,
        "common_dimension": args.k - 1,
        "dtype": dtype,
        "control_variate": args.control_variate,
        "radial_basis": args.radial_basis,
    }
    with index_path.open("a") as index_stream, ProcessPoolExecutor(
        max_workers=args.workers,
        mp_context=multiprocessing.get_context("fork"),
    ) as pool:
        futures = {pool.submit(compile_worker, pair): pair for pair in pending}
        for block_index, future in enumerate(as_completed(futures), 1):
            (left, right), block = future.result()
            relative_path = block_name(left, right)
            output_path = args.output_directory / relative_path
            atomic_save_array(output_path, block)
            record = {
                "left_signature": list(left),
                "right_signature": list(right),
                "path": relative_path,
                "sha256": j_block.file_sha256(output_path),
                "shape": list(block.shape),
            }
            index_stream.write(json.dumps(record, sort_keys=True) + "\n")
            index_stream.flush()
            print(json.dumps({
                "newly_completed": block_index,
                "total_completed": len(completed) + block_index,
                "elapsed_seconds": time.perf_counter() - started,
            }), flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Cheap replay/finalization of recorded exact I/J signature blocks."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import gmpy2

import exact_provenance as provenance
import exact_symmetric_verifier as verifier


HERE = Path(__file__).resolve().parent
DP_PATH = HERE.parent / "symmetry-assembler-design" / "orbit_status_densities.py"
RUNNER_PATH = HERE / "run_parallel_exact.py"


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_groups(
    path: Path, *, kind: str, k: int, candidate_hash: str, binding_hash: str
) -> tuple[gmpy2.mpq, set[tuple[int, ...]], dict]:
    manifest_path = provenance.manifest_path(path)
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != "primegaps-stadlmann-exact-group-manifest-v1":
        raise ValueError(f"unsupported group manifest: {manifest_path}")
    if not manifest.get("completed") or manifest.get("records_sha256") != file_hash(path):
        raise ValueError(f"group manifest is incomplete or does not bind {path}")
    context = manifest.get("context", {})
    if manifest.get("context_sha256") != provenance.value_hash(context):
        raise ValueError(f"invalid context hash in {manifest_path}")
    required = {
        "kind": kind,
        "k": k,
        "candidate_sha256": candidate_hash,
        "candidate_binding_sha256": binding_hash,
        "support_sha256": provenance.value_hash(provenance.SUPPORT),
    }
    for key, expected in required.items():
        if context.get(key) != expected:
            raise ValueError(f"group context {key} does not match: {manifest_path}")
    expected_evaluators = {
        path.name: file_hash(path)
        for path in (
            Path(verifier.__file__),
            Path(verifier.kernel.__file__),
            DP_PATH,
            RUNNER_PATH,
            Path(provenance.__file__),
        )
    }
    if context.get("evaluator_files") != expected_evaluators:
        raise ValueError(f"group context does not match frozen evaluator files: {manifest_path}")
    total = gmpy2.mpq(0)
    signatures = set()
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        record = json.loads(line)
        if record.get("schema") != provenance.RECORD_SCHEMA:
            raise ValueError(f"unbound group record at {path}:{line_number}")
        if record.get("kind") != kind or record.get("k") != k:
            raise ValueError(f"wrong k/kind at {path}:{line_number}")
        if record.get("context_sha256") != manifest["context_sha256"]:
            raise ValueError(f"mixed group context at {path}:{line_number}")
        signature = tuple(record["signature"])
        if signature in signatures:
            raise ValueError(f"duplicate signature at {path}:{line_number}: {signature}")
        signatures.add(signature)
        total += gmpy2.mpq(int(record["numerator"]), int(record["denominator"]))
    if len(signatures) != context.get("group_count"):
        raise ValueError(f"group count does not match manifest: {path}")
    if provenance.value_hash(provenance.signature_payload(signatures)) != context.get(
        "expected_signatures_sha256"
    ):
        raise ValueError(f"signature set does not match manifest: {path}")
    return total, signatures, manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--i-groups", type=Path, required=True)
    parser.add_argument("--j-groups", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--binding", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-groups", type=int, default=2714)
    args = parser.parse_args()
    candidate_hash = file_hash(args.candidate)
    binding_hash = file_hash(args.binding)
    i_value, i_signatures, i_manifest = load_groups(
        args.i_groups, kind="I", k=args.k,
        candidate_hash=candidate_hash, binding_hash=binding_hash
    )
    j_value, j_signatures, j_manifest = load_groups(
        args.j_groups, kind="J", k=args.k,
        candidate_hash=candidate_hash, binding_hash=binding_hash
    )
    i_count, j_count = len(i_signatures), len(j_signatures)
    if (i_count, j_count) != (args.expected_groups, args.expected_groups):
        raise ValueError(
            f"incomplete exact replay: I={i_count}, J={j_count}, "
            f"expected={args.expected_groups}"
        )
    if i_signatures != j_signatures:
        raise ValueError("I and J block files contain different signature sets")
    terms = verifier.rational_terms_from_candidate(args.candidate, args.k)
    expected_i = {
        signature for signature, _slack in verifier.aggregate_i_atoms(terms, args.k)
    }
    feature_groups = verifier.grouped_marginal_coefficients(terms)
    expected_j = set(verifier.grouped_signature_pairs(feature_groups, args.k))
    if i_signatures != expected_i or j_signatures != expected_j:
        raise ValueError("group files do not match signature sets derived from candidate")
    binding = json.loads(args.binding.read_text())
    if not binding.get("exact_termwise_match"):
        raise ValueError("candidate binding did not record an exact termwise match")
    if binding.get("k") != args.k or binding.get("candidate_sha256") != candidate_hash:
        raise ValueError("candidate binding does not match k/candidate")
    numerator = args.k * j_value
    difference = numerator - i_value
    quotient = numerator / i_value
    deficit = 1 - quotient
    gmpy2.get_context().precision = 256
    payload = {
        "schema": "primegaps-stadlmann-240-exact-replay-v1",
        "k": args.k,
        "degree": 21,
        "basis_dimension": 846,
        "basis": "P_(2 lambda)(t) times a radial polynomial of degree b; 2|lambda|+b<=21",
        "rational_candidate": {
            "path": str(args.candidate),
            "sha256": candidate_hash,
            "binding_path": str(args.binding),
            "binding_sha256": file_hash(args.binding),
        },
        "exact_group_inputs": {
            "I": {"path": str(args.i_groups), "sha256": file_hash(args.i_groups), "count": i_count},
            "J": {"path": str(args.j_groups), "sha256": file_hash(args.j_groups), "count": j_count},
            "I_manifest_sha256": file_hash(provenance.manifest_path(args.i_groups)),
            "J_manifest_sha256": file_hash(provenance.manifest_path(args.j_groups)),
        },
        "I": {"numerator": str(i_value.numerator), "denominator": str(i_value.denominator)},
        "J": {"numerator": str(j_value.numerator), "denominator": str(j_value.denominator)},
        "kJ_minus_I": {"numerator": str(difference.numerator), "denominator": str(difference.denominator)},
        "quotient_kJ_over_I": {
            "numerator": str(quotient.numerator),
            "denominator": str(quotient.denominator),
            "decimal_70_digits": format(gmpy2.mpfr(quotient), ".70f"),
        },
        "deficit_one_minus_quotient": {
            "numerator": str(deficit.numerator),
            "denominator": str(deficit.denominator),
            "decimal_70_digits": format(gmpy2.mpfr(deficit), ".70f"),
        },
        "certified_strictly_above_one": bool(i_value > 0 and difference > 0),
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "k": args.k,
        "quotient": payload["quotient_kJ_over_I"]["decimal_70_digits"],
        "deficit": payload["deficit_one_minus_quotient"]["decimal_70_digits"],
        "certified": payload["certified_strictly_above_one"],
    }, indent=2))


if __name__ == "__main__":
    main()

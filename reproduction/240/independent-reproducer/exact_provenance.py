"""Canonical provenance records for checkpointed exact signature contractions."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import platform
import sys

import gmpy2


SCHEMA = "primegaps-stadlmann-exact-group-context-v1"
RECORD_SCHEMA = "primegaps-stadlmann-exact-group-v1"
SUPPORT = {
    "A": ["-3/400", "253/1000"],
    "B": ["3/20", "3/20", "17/100", "17/100", "17/100", "17/100"],
    "U": "521/2000",
    "R": "491/2000",
    "delta": "7/250",
    "epsilon": "3/400",
    "c1": 0,
    "c2": 0,
    "degree": 21,
    "basis_dimension": 846,
    "basis": "m_(2 lambda)(t)*(U-sum(t))^b, 2|lambda|+b<=21",
}


def canonical_bytes(value) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def value_hash(value) -> str:
    return sha256(canonical_bytes(value)).hexdigest()


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def terms_payload(terms) -> list[list[object]]:
    return [
        [list(term.signature), term.slack, str(term.coefficient.numerator),
         str(term.coefficient.denominator)]
        for term in terms
    ]


def signature_payload(signatures) -> list[list[int]]:
    return [list(signature) for signature in sorted(signatures)]


def build_context(*, kind, k, terms, signatures, candidate, binding, evaluator_files):
    candidate = Path(candidate)
    binding = Path(binding)
    candidate_digest = file_hash(candidate)
    binding_payload = json.loads(binding.read_text())
    if not binding_payload.get("exact_termwise_match"):
        raise ValueError("candidate binding does not record an exact termwise match")
    if binding_payload.get("k") != k or binding_payload.get("candidate_sha256") != candidate_digest:
        raise ValueError("candidate binding does not match this k/candidate")
    payload = {
        "schema": SCHEMA,
        "kind": kind,
        "k": k,
        "degree": 21,
        "candidate_sha256": candidate_digest,
        "candidate_binding_sha256": file_hash(binding),
        "terms_sha256": value_hash(terms_payload(terms)),
        "group_count": len(signatures),
        "expected_signatures_sha256": value_hash(signature_payload(signatures)),
        "support": SUPPORT,
        "support_sha256": value_hash(SUPPORT),
        "evaluator_files": {
            Path(path).name: file_hash(Path(path)) for path in evaluator_files
        },
        "environment": {
            "arithmetic": "gmpy2.mpq exact rational",
            "python_implementation": platform.python_implementation(),
            "python_version": platform.python_version(),
            "byteorder": sys.byteorder,
            "gmpy2_version": gmpy2.version(),
            "gmp_version": gmpy2.mp_version(),
            "mpfr_version": gmpy2.mpfr_version(),
            "mpc_version": gmpy2.mpc_version(),
        },
    }
    return payload, value_hash(payload)


def manifest_path(groups_path: Path) -> Path:
    return Path(str(groups_path) + ".manifest.json")

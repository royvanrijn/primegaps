#!/usr/bin/env python3
"""Render the finalized D27 arithmetic enclosure as a small Lean module."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import tempfile


I_UPPER_SHA256 = "6084ceb8438c7dca6c33d935d1c05ffa8723f6c0c9a474a9d1750c04568246ac"
UNRESTRICTED_J_SHA256 = "27bb19f7f1f0ab57159153e6499934e6deac69edc3763d62d5a6661c3b1aee9a"


def exact_fraction(value: object, field: str) -> Fraction:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{field} must be a [numerator, denominator] pair")
    numerator, denominator = (int(part) for part in value)
    if denominator <= 0:
        raise ValueError(f"{field} denominator must be positive")
    return Fraction(numerator, denominator)


def render(result: dict[str, object]) -> str:
    if result.get("schema") != "primegaps-boundary-J-arb-result-v1":
        raise ValueError("unexpected result schema")
    if not result.get("checkpoint_complete") or not result.get("certification_eligible"):
        raise ValueError("D27 checkpoint is not complete and certification-eligible")

    gate = result.get("exact_gate")
    if not isinstance(gate, dict):
        raise ValueError("result has no exact_gate")
    if gate.get("i_upper_sha256") != I_UPPER_SHA256:
        raise ValueError("exact I input does not match the formalized candidate")
    if gate.get("unrestricted_j_sha256") != UNRESTRICTED_J_SHA256:
        raise ValueError("exact unrestricted J input does not match the formalized candidate")
    if not gate.get("certified_strictly_above_one"):
        raise ValueError("the finalized enclosure does not certify the strict crossing")

    lower = exact_fraction(gate.get("normalized_correction_lower"), "normalized_correction_lower")
    required = exact_fraction(gate.get("required_normalized_correction"), "required_normalized_correction")
    legal_lower = exact_fraction(gate.get("normalized_legal_kJ_lower"), "normalized_legal_kJ_lower")
    i_upper = exact_fraction(gate.get("normalized_i_upper"), "normalized_i_upper")
    if not lower > required or not legal_lower > i_upper:
        raise ValueError("exact gate fields do not contain a strict crossing")

    input_hash = result.get("input_sha256", "unknown")
    manifest_hash = result.get("manifest_sha256", "unknown")
    return f'''import PrimeGaps236.Stadlmann.BoundaryCertificate

/-!
# Generated D27 boundary arithmetic certificate

Generated from a complete Arb enclosure. The source checkpoint SHA-256 is
`{input_hash}` and its manifest SHA-256 is `{manifest_hash}`.

This proves only the exact scalar comparison. The analytic theorem connecting
the Arb enclosure to the shaped-support integral remains a separate obligation.
-/

namespace Gaps236.Stadlmann

def d27BoundaryArithmeticCertificate : BoundaryArithmeticCertificate where
  lower := {lower.numerator} / {lower.denominator}
  passes := by norm_num [normalizedI, normalizedUnrestrictedKJ]

end Gaps236.Stadlmann
'''


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    result = json.loads(args.result.read_text())
    text = render(result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=args.output.parent, prefix=f".{args.output.name}.", delete=False
    ) as stream:
        temporary = Path(stream.name)
        stream.write(text)
    temporary.replace(args.output)


if __name__ == "__main__":
    main()

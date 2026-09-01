"""Rationalization and exact certificates for a strict Rayleigh inequality."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import json
from math import gcd
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from .exact_matrix import ExactSymmetricMatrix


CERTIFICATE_FORMAT = "primegaps-rayleigh-certificate-v1"


def _fraction_payload(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _parse_fraction(payload: object, field: str) -> Fraction:
    if not isinstance(payload, dict) or set(payload) != {"numerator", "denominator"}:
        raise ValueError(f"{field} must contain numerator and denominator")
    value = Fraction(int(payload["numerator"]), int(payload["denominator"]))
    return value


def _primitive(vector: Sequence[int]) -> tuple[int, ...]:
    values = tuple(int(value) for value in vector)
    divisor = 0
    for value in values:
        divisor = gcd(divisor, abs(value))
    if divisor == 0:
        raise ValueError("certificate vector cannot be zero")
    values = tuple(value // divisor for value in values)
    first_nonzero = next(value for value in values if value)
    if first_nonzero < 0:
        values = tuple(-value for value in values)
    return values


def _candidate_scales(max_scale: int) -> tuple[int, ...]:
    if max_scale < 1:
        raise ValueError("max_scale must be positive")
    candidates = {1, max_scale}
    power = 1
    while power <= max_scale:
        for multiplier in (1, 2, 5):
            value = multiplier * power
            if value <= max_scale:
                candidates.add(value)
        power *= 10
    return tuple(sorted(candidates))


def integer_direction(vector: Iterable[float], scale: int) -> tuple[int, ...]:
    """Round a normalized float direction at one common integer scale."""
    values = np.asarray(tuple(vector), dtype=float)
    if values.ndim != 1 or len(values) == 0 or not np.all(np.isfinite(values)):
        raise ValueError("candidate vector must be a finite non-empty sequence")
    maximum = float(np.max(np.abs(values)))
    if maximum == 0.0:
        raise ValueError("candidate vector cannot be zero")
    if scale < 1 or scale > 2**52:
        raise ValueError("scale must be between 1 and 2**52")
    rounded = np.rint(values * (float(scale) / maximum))
    return _primitive(tuple(int(value) for value in rounded))


@dataclass(frozen=True)
class RationalizedCandidate:
    vector: tuple[int, ...]
    scale: int
    quadratic_form_m1: Fraction
    quadratic_form_m2: Fraction
    attempted_scales: int

    @property
    def quotient(self) -> Fraction:
        return self.quadratic_form_m2 / self.quadratic_form_m1

    @property
    def proves_strict_inequality(self) -> bool:
        return self.quadratic_form_m1 > 0 and self.quadratic_form_m2 > self.quadratic_form_m1


def rationalize_candidate(
    m1: ExactSymmetricMatrix,
    m2: ExactSymmetricMatrix,
    numerical_vector: Iterable[float],
    *,
    max_scale: int = 1_000_000,
) -> RationalizedCandidate:
    """Choose the best exact quotient over logarithmically spaced scales.

    A single common scale avoids the enormous least-common-multiple created by
    rationalizing every coordinate independently. Every candidate is evaluated
    against the exact matrices; floating-point scores are never used to select
    the certificate.
    """
    if m1.dimension != m2.dimension:
        raise ValueError("matrix dimensions differ")
    if max_scale > 2**52:
        raise ValueError("max_scale must not exceed 2**52 at float64 precision")
    values = tuple(float(value) for value in numerical_vector)
    if len(values) != m1.dimension:
        raise ValueError("matrix/vector dimensions differ")
    best: RationalizedCandidate | None = None
    seen: set[tuple[int, ...]] = set()
    scales = _candidate_scales(max_scale)
    for scale in scales:
        vector = integer_direction(values, scale)
        if vector in seen:
            continue
        seen.add(vector)
        q1 = m1.quadratic_form(vector)
        if q1 <= 0:
            continue
        q2 = m2.quadratic_form(vector)
        candidate = RationalizedCandidate(vector, scale, q1, q2, len(scales))
        if best is None or candidate.quotient > best.quotient:
            best = candidate
    if best is None:
        raise ValueError("no rationalized candidate has a positive M1 quadratic form")
    return best


@dataclass(frozen=True)
class ExactCertificate:
    dimension: int
    m1_sha256: str
    m2_sha256: str
    vector: tuple[int, ...]
    quadratic_form_m1: Fraction
    quadratic_form_m2: Fraction
    difference: Fraction
    quotient: Fraction
    rationalization_scale: int

    def to_dict(self) -> dict[str, object]:
        return {
            "format": CERTIFICATE_FORMAT,
            "claim": "c^T M2 c > c^T M1 c",
            "dimension": self.dimension,
            "matrix_sha256": {"m1": self.m1_sha256, "m2": self.m2_sha256},
            "vector": [str(value) for value in self.vector],
            "quadratic_form_m1": _fraction_payload(self.quadratic_form_m1),
            "quadratic_form_m2": _fraction_payload(self.quadratic_form_m2),
            "difference_m2_minus_m1": _fraction_payload(self.difference),
            "rayleigh_quotient": _fraction_payload(self.quotient),
            "rationalization_scale": str(self.rationalization_scale),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> "ExactCertificate":
        if payload.get("format") != CERTIFICATE_FORMAT:
            raise ValueError("unsupported certificate format")
        if payload.get("claim") != "c^T M2 c > c^T M1 c":
            raise ValueError("unsupported certificate claim")
        hashes = payload.get("matrix_sha256")
        if not isinstance(hashes, dict) or set(hashes) != {"m1", "m2"}:
            raise ValueError("certificate matrix hashes are malformed")
        raw_vector = payload.get("vector")
        if not isinstance(raw_vector, list):
            raise ValueError("certificate vector must be a list")
        return cls(
            dimension=int(payload["dimension"]),
            m1_sha256=str(hashes["m1"]),
            m2_sha256=str(hashes["m2"]),
            vector=tuple(int(value) for value in raw_vector),
            quadratic_form_m1=_parse_fraction(payload.get("quadratic_form_m1"), "quadratic_form_m1"),
            quadratic_form_m2=_parse_fraction(payload.get("quadratic_form_m2"), "quadratic_form_m2"),
            difference=_parse_fraction(
                payload.get("difference_m2_minus_m1"), "difference_m2_minus_m1"
            ),
            quotient=_parse_fraction(payload.get("rayleigh_quotient"), "rayleigh_quotient"),
            rationalization_scale=int(payload["rationalization_scale"]),
        )


@dataclass(frozen=True)
class VerificationResult:
    valid: bool
    quotient: Fraction
    difference: Fraction
    message: str


def build_certificate(
    m1: ExactSymmetricMatrix,
    m2: ExactSymmetricMatrix,
    candidate: RationalizedCandidate,
) -> ExactCertificate:
    """Build a certificate only when the exact strict inequality holds."""
    if m1.dimension != m2.dimension or len(candidate.vector) != m1.dimension:
        raise ValueError("matrix/vector dimensions differ")
    original_q1 = m1.quadratic_form(candidate.vector)
    original_q2 = m2.quadratic_form(candidate.vector)
    if (
        original_q1 != candidate.quadratic_form_m1
        or original_q2 != candidate.quadratic_form_m2
    ):
        raise ValueError("candidate quadratic forms do not match supplied matrices")
    vector = _primitive(candidate.vector)
    q1 = m1.quadratic_form(vector)
    q2 = m2.quadratic_form(vector)
    if q1 <= 0:
        raise ValueError("certificate requires c^T M1 c > 0")
    if q2 <= q1:
        raise ValueError("candidate does not prove c^T M2 c > c^T M1 c")
    return ExactCertificate(
        dimension=m1.dimension,
        m1_sha256=m1.semantic_sha256(),
        m2_sha256=m2.semantic_sha256(),
        vector=vector,
        quadratic_form_m1=q1,
        quadratic_form_m2=q2,
        difference=q2 - q1,
        quotient=q2 / q1,
        rationalization_scale=candidate.scale,
    )


def verify_certificate(
    m1: ExactSymmetricMatrix,
    m2: ExactSymmetricMatrix,
    certificate: ExactCertificate,
) -> VerificationResult:
    """Recompute hashes and exact forms, rejecting any inconsistent field."""
    if certificate.dimension != m1.dimension or m1.dimension != m2.dimension:
        raise ValueError("certificate/matrix dimensions differ")
    if len(certificate.vector) != certificate.dimension:
        raise ValueError("certificate vector dimension is wrong")
    if certificate.vector != _primitive(certificate.vector):
        raise ValueError("certificate vector is not in canonical primitive form")
    if certificate.m1_sha256 != m1.semantic_sha256():
        raise ValueError("certificate M1 hash does not match supplied M1")
    if certificate.m2_sha256 != m2.semantic_sha256():
        raise ValueError("certificate M2 hash does not match supplied M2")
    q1 = m1.quadratic_form(certificate.vector)
    q2 = m2.quadratic_form(certificate.vector)
    difference = q2 - q1
    if certificate.quadratic_form_m1 != q1:
        raise ValueError("recorded M1 quadratic form is incorrect")
    if certificate.quadratic_form_m2 != q2:
        raise ValueError("recorded M2 quadratic form is incorrect")
    if certificate.difference != difference:
        raise ValueError("recorded quadratic-form difference is incorrect")
    if q1 == 0:
        raise ValueError("certificate Rayleigh quotient has zero denominator")
    quotient = q2 / q1
    if certificate.quotient != quotient:
        raise ValueError("recorded Rayleigh quotient is incorrect")
    if q1 <= 0:
        raise ValueError("certificate does not establish a positive denominator")
    if difference <= 0:
        raise ValueError("certificate does not establish the strict inequality")
    return VerificationResult(True, quotient, difference, "exact strict inequality verified")


def save_certificate(path: str | Path, certificate: ExactCertificate) -> None:
    Path(path).write_text(
        json.dumps(certificate.to_dict(), sort_keys=True, separators=(",", ":")) + "\n"
    )


def load_certificate(path: str | Path) -> ExactCertificate:
    payload = json.loads(Path(path).read_text())
    if not isinstance(payload, dict):
        raise ValueError("certificate JSON must be an object")
    return ExactCertificate.from_dict(payload)

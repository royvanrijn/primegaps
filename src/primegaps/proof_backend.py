"""End-to-end numerical search and exact certificate command line interface."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
from typing import Sequence

from .certificate import (
    ExactCertificate,
    RationalizedCandidate,
    build_certificate,
    load_certificate,
    rationalize_candidate,
    save_certificate,
    verify_certificate,
)
from .eigen import GeneralizedEigenResult, solve_generalized_eigenproblem
from .exact_matrix import ExactSymmetricMatrix, load_matrix_pair


@dataclass(frozen=True)
class ProofBackendResult:
    numerical: GeneralizedEigenResult
    rationalized: RationalizedCandidate
    certificate: ExactCertificate | None

    @property
    def proves_strict_inequality(self) -> bool:
        return self.certificate is not None


def solve_and_certify(
    m1: ExactSymmetricMatrix,
    m2: ExactSymmetricMatrix,
    *,
    method: str = "auto",
    dense_threshold: int = 1_200,
    tolerance: float = 1e-11,
    max_iterations: int = 180,
    max_scale: int = 1_000_000,
    seed: int = 0,
    prefer_sparse: bool = True,
) -> ProofBackendResult:
    """Find the best numerical quotient and attempt an exact certificate."""
    if m1.dimension != m2.dimension:
        raise ValueError("matrix dimensions differ")
    scipy_available = importlib.util.find_spec("scipy") is not None
    use_sparse = (
        prefer_sparse
        and scipy_available
        and m1.storage == "sparse-upper"
        and m2.storage == "sparse-upper"
        and method in {"auto", "sparse"}
    )
    if use_sparse:
        numerical_m1 = m1.to_scipy_sparse()
        numerical_m2 = m2.to_scipy_sparse()
    else:
        if method == "sparse":
            raise RuntimeError("sparse solving requires sparse exact inputs and SciPy")
        numerical_m1 = m1.to_dense_float()
        numerical_m2 = m2.to_dense_float()
    numerical = solve_generalized_eigenproblem(
        numerical_m1,
        numerical_m2,
        method=method,
        dense_threshold=dense_threshold,
        tolerance=tolerance,
        max_iterations=max_iterations,
        seed=seed,
    )
    rationalized = rationalize_candidate(
        m1, m2, numerical.vector, max_scale=max_scale
    )
    certificate = (
        build_certificate(m1, m2, rationalized)
        if rationalized.proves_strict_inequality
        else None
    )
    return ProofBackendResult(numerical, rationalized, certificate)


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _result_payload(result: ProofBackendResult) -> dict[str, object]:
    candidate = result.rationalized
    return {
        "numerical_quotient": result.numerical.quotient,
        "exact_candidate_quotient": _fraction_text(candidate.quotient),
        "exact_difference_m2_minus_m1": _fraction_text(
            candidate.quadratic_form_m2 - candidate.quadratic_form_m1
        ),
        "rationalization_scale": candidate.scale,
        "certificate_produced": result.proves_strict_inequality,
        "diagnostics": result.numerical.diagnostics.to_dict(),
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="primegaps-proof",
        description="Numerical generalized-eigenvalue search plus exact certificate replay",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    solve = subparsers.add_parser("solve", help="solve an exact matrix pair and attempt a certificate")
    solve.add_argument("matrices", type=Path, help="matrix-pair JSON or packed int64 NPZ")
    solve.add_argument("--certificate", type=Path, required=True, help="output certificate JSON")
    solve.add_argument("--result", type=Path, help="optional numerical diagnostics JSON")
    solve.add_argument("--method", choices=("auto", "dense", "iterative", "sparse"), default="auto")
    solve.add_argument("--dense-threshold", type=int, default=1_200)
    solve.add_argument("--tolerance", type=float, default=1e-11)
    solve.add_argument("--max-iterations", type=int, default=180)
    solve.add_argument("--max-scale", type=int, default=1_000_000)
    solve.add_argument("--seed", type=int, default=0)
    solve.add_argument(
        "--dense-input",
        action="store_true",
        help="materialize sparse exact matrices as dense numerical arrays",
    )

    verify = subparsers.add_parser("verify", help="replay a certificate using exact arithmetic")
    verify.add_argument("matrices", type=Path, help="matrix-pair JSON or packed int64 NPZ")
    verify.add_argument("certificate", type=Path, help="certificate JSON")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.command == "verify":
        m1, m2 = load_matrix_pair(args.matrices)
        certificate = load_certificate(args.certificate)
        verified = verify_certificate(m1, m2, certificate)
        print(
            json.dumps(
                {
                    "valid": verified.valid,
                    "message": verified.message,
                    "exact_quotient": _fraction_text(verified.quotient),
                    "exact_difference_m2_minus_m1": _fraction_text(verified.difference),
                },
                sort_keys=True,
            )
        )
        return 0

    m1, m2 = load_matrix_pair(args.matrices)
    result = solve_and_certify(
        m1,
        m2,
        method=args.method,
        dense_threshold=args.dense_threshold,
        tolerance=args.tolerance,
        max_iterations=args.max_iterations,
        max_scale=args.max_scale,
        seed=args.seed,
        prefer_sparse=not args.dense_input,
    )
    payload = _result_payload(result)
    if result.certificate is not None:
        save_certificate(args.certificate, result.certificate)
        payload["certificate"] = str(args.certificate)
        payload["certificate_m1_sha256"] = result.certificate.m1_sha256
        payload["certificate_m2_sha256"] = result.certificate.m2_sha256
    else:
        payload["failure"] = "best rationalized candidate does not prove the strict inequality"
    output = json.dumps(payload, sort_keys=True)
    print(output)
    if args.result:
        args.result.write_text(output + "\n")
    return 0 if result.certificate is not None else 2


if __name__ == "__main__":
    raise SystemExit(main())

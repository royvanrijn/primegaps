#!/usr/bin/env python3
"""Exploratory physical-face model for an ideal quarter-rough prime detector.

This is deliberately a discovery calculation, not a certificate.  It ports the
positive midpoint/Dickman cap model from the pinned PrimeGaps186 source to
float64, exposes the 77-by-77 quadratic forms, and permits the ambient
dimension to be changed from 40 to 39 or 38.  The published source masks and
trial basis are retained.  No source-loss or parity-correlation theorem is
assumed beyond the explicitly reported idealization.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from fractions import Fraction
from functools import cache
from hashlib import sha256
import json
from itertools import product
from math import comb, prod
from pathlib import Path
import tempfile

import numpy as np
from scipy import linalg
from scipy.signal import fftconvolve

from primegaps.parity import rough_factor_constants as compute_rough_factor_constants


def repository_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "pyproject.toml").is_file():
            return parent
    raise FileNotFoundError("could not locate the repository root")


ROOT = repository_root()
DEFAULT_INPUT = ROOT / "reproduction/186/physical-parity-input.json"


@cache
def _cap_block_partitions(signature: tuple[int, ...]):
    if not signature:
        return (((), 1),)
    q, answer = signature[-1], defaultdict(int)
    for blocks, multiplicity in _cap_block_partitions(signature[:-1]):
        answer[tuple(sorted((*blocks, q)))] += multiplicity
        for value, copies in Counter(blocks).items():
            joined = list(blocks)
            joined.remove(value)
            answer[tuple(sorted(joined + [value + q]))] += multiplicity * copies
    return tuple(sorted(answer.items(), key=lambda item: (len(item[0]), item[0])))


def moment_terms(count: int, signature: tuple[int, ...]):
    if count < 0 or any(q < 0 for q in signature):
        raise ValueError("negative moment dimension or exponent")
    return tuple(
        (blocks, multiplicity * prod(range(count - len(blocks) + 1, count + 1)))
        for blocks, multiplicity in _cap_block_partitions(tuple(sorted(signature)))
        if len(blocks) <= count
    )


@cache
def fiber_splits(signature: tuple[int, ...]):
    counts, answer = tuple(sorted(Counter(signature).items())), []
    for chosen in product(*(range(count + 1) for _, count in counts)):
        remaining, exponent, multiplicity = [], 0, 1
        for (power, count), selected in zip(counts, chosen, strict=True):
            remaining.extend([power] * (count - selected))
            exponent += power * selected
            multiplicity *= comb(count, selected)
        answer.append((tuple(remaining), exponent, multiplicity))
    return tuple(answer)


@dataclass(frozen=True)
class PhysicalInputs:
    CAP_SIGNATURES: tuple[tuple[int, ...], ...]
    CAP_COEFFICIENTS: tuple[tuple[int, ...], ...]
    CAP_SHELL_DATA: dict[str, tuple[tuple[Fraction, Fraction], ...]]
    DERIVED_INPUTS: dict[str, dict[str, Fraction]]
    RHO_STAR: Fraction
    OUTER_RADIUS: Fraction
    CENTER: Fraction
    PROFILE: tuple[Fraction, Fraction, Fraction]
    UPSTREAM: dict[str, str]

    @staticmethod
    def moment_terms(count: int, signature: tuple[int, ...]):
        return moment_terms(count, signature)

    @staticmethod
    def fiber_splits(signature: tuple[int, ...]):
        return fiber_splits(signature)


def load_inputs(path: Path = DEFAULT_INPUT) -> PhysicalInputs:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "primegaps.physical-parity-input.v1":
        raise ValueError("unexpected physical parity input schema")
    profile = payload["profile"]
    return PhysicalInputs(
        CAP_SIGNATURES=tuple(
            tuple(int(Fraction(value)) for value in signature)
            for signature in payload["signatures"]
        ),
        CAP_COEFFICIENTS=tuple(tuple(map(int, row)) for row in payload["coefficients"]),
        CAP_SHELL_DATA={
            label: tuple((Fraction(upper), Fraction(cap)) for upper, cap in rows)
            for label, rows in payload["cap_shell_data"].items()
        },
        DERIVED_INPUTS={
            "hybrid": {key: Fraction(value) for key, value in payload["hybrid"].items()}
        },
        RHO_STAR=Fraction(payload["rho_star"]),
        OUTER_RADIUS=Fraction(payload["outer_radius"]),
        CENTER=Fraction(payload["center"]),
        PROFILE=(
            Fraction(profile["weight"]),
            Fraction(profile["slow"]),
            Fraction(profile["fast"]),
        ),
        UPSTREAM=dict(payload["upstream"]),
    )


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, prefix=path.name, delete=False
    ) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        temporary = Path(stream.name)
    temporary.replace(path)


def truncated_convolution(left: np.ndarray, right: np.ndarray, length: int) -> np.ndarray:
    # Direct convolution is an important small-mesh oracle: high powers have
    # long exact-zero prefixes, for which an FFT's absolute roundoff can exceed
    # the genuine far-tail coefficients by many orders of magnitude.
    result = (
        np.convolve(left, right)[:length]
        if length <= 8192
        else fftconvolve(left, right)[:length]
    )
    # FFT roundoff is governed by the l1 convolution norm, including in the
    # exact leading-zero region where comparison with max(result) is misleading.
    scale = max(float(np.sum(np.abs(left)) * np.sum(np.abs(right))), np.finfo(float).tiny)
    if float(np.min(result)) < -2e-11 * scale:
        raise ArithmeticError("positive FFT convolution acquired a material negative part")
    return np.maximum(result, 0.0)


def rough_factor_constants(beta: float) -> dict[str, float]:
    """Serialize the reusable degree-two/three rough constants for this run.

    Prime powers and repeated prime factors are lower order at this scale.
    """
    constants = compute_rough_factor_constants(beta)
    return {
        **asdict(constants),
        "detector_degree": constants.detector_degree,
        "rough_carrier": constants.rough_carrier,
    }


@dataclass(frozen=True)
class Shell:
    lower: Fraction
    upper: Fraction
    ceiling: Fraction


class FloatPhysicalCapModel:
    def __init__(self, upstream, *, dimension: int, intervals: int):
        if not 4 <= dimension < intervals:
            raise ValueError("dimension must lie between 4 and the mesh size")
        self.upstream = upstream
        self.k = int(dimension)
        self.intervals = int(intervals)
        self.n = intervals - dimension
        self.S = upstream.OUTER_RADIUS
        self.hq = self.S / intervals
        self.h = float(self.hq)
        self.center = float(upstream.CENTER)
        self.signatures = tuple(map(tuple, upstream.CAP_SIGNATURES))
        self.descriptors = tuple((signature, degree) for signature in self.signatures for degree in range(7))
        self.basis_size = len(self.descriptors)
        self.midpoints = (np.arange(self.n, dtype=float) + 0.5) * self.h
        weight, slow, fast = map(float, upstream.PROFILE)
        self.root = weight / (1 + slow * self.midpoints) + (1 - weight) / (
            1 + fast * self.midpoints
        )
        self.Z = float(self.root @ self.root)
        self.coordinate_weight = self.root * self.root / self.Z
        self.radial = (np.arange(self.n, dtype=float) + dimension / 2) * self.h - self.center
        self.radial_powers = tuple(self.radial**degree for degree in range(13))
        self.shells = {
            label: self._aligned_shells(label) for label in ("outer", "base", "enlarged", "full")
        }
        finite_caps = sorted(
            {shell.ceiling for shells in self.shells.values() for shell in shells}
        )
        self.caps = tuple(finite_caps)
        self.outer_masks = tuple(self._shell_mask(shell, dimension) for shell in self.shells["outer"])
        self.inner_allowed = np.zeros((3, len(self.caps), self.n), dtype=bool)
        for role, label in enumerate(("base", "enlarged", "full")):
            for shell in self.shells[label]:
                for layer, cap in enumerate(self.caps):
                    if cap <= shell.ceiling:
                        self.inner_allowed[role, layer] |= self._shell_mask(shell, dimension - 1)
        self._survival: dict[Fraction, np.ndarray] = {}
        self._weighted: dict[tuple[Fraction, int], np.ndarray] = {}
        self._probability_power: dict[tuple[Fraction, int], np.ndarray] = {}
        self._block_product: dict[tuple[Fraction, tuple[int, ...]], np.ndarray] = {}
        self._moment: dict[tuple[Fraction, int, tuple[int, ...]], np.ndarray] = {}

    def _align(self, value) -> Fraction:
        value = Fraction(value)
        return (value // self.hq) * self.hq

    def _aligned_shells(self, label: str) -> tuple[Shell, ...]:
        result = []
        lower = Fraction(0)
        for upper, cap in self.upstream.CAP_SHELL_DATA[label]:
            aligned = self._align(cap)
            if result and result[-1].ceiling == aligned:
                lower = result.pop().lower
            result.append(Shell(lower, Fraction(upper), aligned))
            lower = Fraction(upper)
        return tuple(result)

    def _shell_mask(self, shell: Shell, count: int) -> np.ndarray:
        indices = np.arange(self.n, dtype=np.int64) + count
        return (indices > int(shell.lower // self.hq)) & (
            indices <= int(shell.upper // self.hq)
        )

    def survival(self, cap: Fraction) -> np.ndarray:
        cap = self._align(cap)
        if cap in self._survival:
            return self._survival[cap]
        cap_index = int(cap / self.hq)
        endpoints = np.ones(self.n + 1, dtype=float)
        for j in range(cap_index, self.n):
            decrement = 0.5 * (
                endpoints[j - cap_index] / j
                + endpoints[j + 1 - cap_index] / (j + 1)
            )
            endpoints[j + 1] = max(0.0, endpoints[j] - decrement)
        density = np.ones(self.n, dtype=float)
        if cap_index < self.n:
            density[cap_index:] = 0.5 * (
                endpoints[cap_index:self.n] + endpoints[cap_index + 1 : self.n + 1]
            )
        self._survival[cap] = density
        return density

    def weighted(self, cap: Fraction, exponent: int) -> np.ndarray:
        key = self._align(cap), int(exponent)
        if key not in self._weighted:
            self._weighted[key] = (
                self.coordinate_weight * self.survival(key[0]) * self.midpoints**key[1]
            )
        return self._weighted[key]

    def probability_power(self, cap: Fraction, count: int) -> np.ndarray:
        key = self._align(cap), int(count)
        if key not in self._probability_power:
            if count == 0:
                value = np.zeros(self.n, dtype=float)
                value[0] = 1.0
            elif count == 1:
                value = self.weighted(cap, 0)
            else:
                half = self.probability_power(cap, count // 2)
                value = truncated_convolution(half, half, self.n)
                if count & 1:
                    value = truncated_convolution(value, self.weighted(cap, 0), self.n)
            self._probability_power[key] = value
        return self._probability_power[key]

    def block_product(self, cap: Fraction, blocks: tuple[int, ...]) -> np.ndarray:
        key = self._align(cap), tuple(blocks)
        if key not in self._block_product:
            if not blocks:
                value = self.probability_power(cap, 0)
            else:
                value = self.weighted(cap, blocks[0])
                for exponent in blocks[1:]:
                    value = truncated_convolution(value, self.weighted(cap, exponent), self.n)
            self._block_product[key] = value
        return self._block_product[key]

    def moment(self, cap: Fraction, count: int, signature: tuple[int, ...]) -> np.ndarray:
        key = self._align(cap), int(count), tuple(sorted(signature))
        if key not in self._moment:
            result = np.zeros(self.n, dtype=float)
            for blocks, coefficient in self.upstream.moment_terms(count, key[2]):
                term = truncated_convolution(
                    self.probability_power(cap, count - len(blocks)),
                    self.block_product(cap, blocks),
                    self.n,
                )
                result += coefficient * term
            self._moment[key] = result
        return self._moment[key]

    def mass_matrix(self) -> np.ndarray:
        result = np.zeros((self.basis_size, self.basis_size), dtype=float)
        for shell, mask in zip(self.shells["outer"], self.outer_masks):
            moments = {}
            for right, (right_signature, right_degree) in enumerate(self.descriptors):
                for left in range(right + 1):
                    left_signature, left_degree = self.descriptors[left]
                    signature = tuple(sorted(left_signature + right_signature))
                    if signature not in moments:
                        moments[signature] = self.moment(shell.ceiling, self.k, signature)
                    value = float(
                        np.dot(
                            self.radial_powers[left_degree + right_degree][mask],
                            moments[signature][mask],
                        )
                    )
                    result[left, right] += value
                    if left != right:
                        result[right, left] += value
        return result

    @cache
    def _affine_shell_features(self, shell_index: int) -> dict[tuple[int, ...], np.ndarray]:
        shell = self.shells["outer"][shell_index]
        mask = self.outer_masks[shell_index]
        answer: dict[tuple[int, ...], np.ndarray] = {}
        for basis_index, (signature, degree) in enumerate(self.descriptors):
            for remaining, exponent, multiplicity in self.upstream.fiber_splits(signature):
                if remaining not in answer:
                    answer[remaining] = np.zeros((self.basis_size, self.n), dtype=float)
                radial = self.radial_powers[degree] * mask
                fiber = self.root * self.survival(shell.ceiling) * self.midpoints**exponent
                correlation = fftconvolve(radial, fiber[::-1])[self.n - 1 : 2 * self.n - 1]
                answer[remaining][basis_index] += multiplicity * correlation
        return answer

    def _layer_features(self, cap: Fraction) -> dict[tuple[int, ...], np.ndarray]:
        answer: dict[tuple[int, ...], np.ndarray] = {}
        for shell_index, shell in enumerate(self.shells["outer"]):
            if cap > shell.ceiling:
                continue
            for signature, values in self._affine_shell_features(shell_index).items():
                if signature not in answer:
                    answer[signature] = values.copy()
                else:
                    answer[signature] += values
        return answer

    def face_matrices(self) -> dict[str, np.ndarray]:
        result = {name: np.zeros((self.basis_size, self.basis_size), dtype=float) for name in ("J0", "Jplus", "Jtail")}
        previous: Fraction | None = None
        for layer, cap in enumerate(self.caps):
            features = self._layer_features(cap)
            if not features:
                previous = cap
                continue
            masks = (
                self.inner_allowed[0, layer],
                self.inner_allowed[1, layer] & ~self.inner_allowed[0, layer],
                self.inner_allowed[2, layer] & ~self.inner_allowed[1, layer],
            )
            moment_differences = {}
            signatures = tuple(sorted(features, key=lambda value: (len(value), value)))
            for right_index, right_signature in enumerate(signatures):
                right = features[right_signature]
                for left_index in range(right_index + 1):
                    left_signature = signatures[left_index]
                    left = features[left_signature]
                    joined = tuple(sorted(left_signature + right_signature))
                    if joined not in moment_differences:
                        values = self.moment(cap, self.k - 1, joined).copy()
                        if previous is not None:
                            values -= self.moment(previous, self.k - 1, joined)
                        scale = max(float(np.max(np.abs(values))), np.finfo(float).tiny)
                        if float(np.min(values)) < -2e-6 * scale:
                            raise ArithmeticError(
                                "a cap-layer moment difference became materially negative: "
                                f"layer={layer} signature={joined} min={np.min(values):.6e} "
                                f"maxabs={scale:.6e}"
                            )
                        moment_differences[joined] = np.maximum(values, 0.0)
                    moment = moment_differences[joined]
                    for name, mask in zip(result, masks):
                        if not np.any(mask):
                            continue
                        weighted_right = right[:, mask] * moment[mask]
                        block = left[:, mask] @ weighted_right.T
                        result[name] += block
                        if left_index != right_index:
                            result[name] += block.T
            previous = cap
        normalization = self.k * self.h / self.Z
        for name in result:
            result[name] = normalization * (result[name] + result[name].T) / 2
        result["Jfull"] = result["J0"] + result["Jplus"] + result["Jtail"]
        return result

    def matrices(self) -> dict[str, np.ndarray]:
        result = {"I": self.mass_matrix()}
        result.update(self.face_matrices())
        return result


def rayleigh(matrix: np.ndarray, metric: np.ndarray, vector: np.ndarray) -> float:
    return float(vector @ matrix @ vector) / float(vector @ metric @ vector)


def spectral_generalized_maximum(
    matrix: np.ndarray, metric: np.ndarray, cutoff: float
) -> dict[str, object]:
    diagonal = np.diag(metric)
    if np.any(diagonal <= 0):
        raise ArithmeticError("mass matrix has a nonpositive diagonal")
    scale = 1 / np.sqrt(diagonal)
    equilibrated = scale[:, None] * metric * scale[None, :]
    values, vectors = linalg.eigh((equilibrated + equilibrated.T) / 2)
    keep = values > cutoff * values[-1]
    whitening = scale[:, None] * (vectors[:, keep] / np.sqrt(values[keep]))
    projected = whitening.T @ matrix @ whitening
    objective_values, objective_vectors = linalg.eigh((projected + projected.T) / 2)
    vector = whitening @ objective_vectors[:, -1]
    vector /= np.max(np.abs(vector))
    quotient = rayleigh(matrix, metric, vector)
    residual = matrix @ vector - quotient * (metric @ vector)
    projected_residual = whitening.T @ residual
    relative_residual = float(
        np.linalg.norm(projected_residual)
        / max(np.linalg.norm(whitening.T @ (matrix @ vector)), np.finfo(float).tiny)
    )
    return {
        "cutoff": cutoff,
        "retained_dimension": int(np.count_nonzero(keep)),
        "mass_eigenvalue_min": float(values[0]),
        "mass_eigenvalue_max": float(values[-1]),
        "quotient": quotient,
        "projected_relative_residual": relative_residual,
        "vector": vector,
    }


def evaluate_dimension(upstream, dimension: int, intervals: int, beta: float) -> tuple[dict[str, object], dict[str, np.ndarray]]:
    model = FloatPhysicalCapModel(upstream, dimension=dimension, intervals=intervals)
    matrices = model.matrices()
    hybrid = upstream.DERIVED_INPUTS["hybrid"]
    a, b = Fraction(hybrid["a"]), Fraction(hybrid["b"])
    hybrid_matrix = matrices["J0"] + float(a + b) * matrices["Jplus"] + float(b) * matrices["Jtail"]
    published_vector = np.asarray(
        [Fraction(value, 10**10) for row in upstream.CAP_COEFFICIENTS for value in row],
        dtype=float,
    )
    rho = float(upstream.RHO_STAR)
    fixed = {
        "hybrid_score": rho * rayleigh(hybrid_matrix, matrices["I"], published_vector),
        "full_face_score": rho * rayleigh(matrices["Jfull"], matrices["I"], published_vector),
        "I": float(published_vector @ matrices["I"] @ published_vector),
        "J0": float(published_vector @ matrices["J0"] @ published_vector),
        "Jplus": float(published_vector @ matrices["Jplus"] @ published_vector),
        "Jtail": float(published_vector @ matrices["Jtail"] @ published_vector),
    }
    cutoff_runs = []
    for cutoff in (1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13):
        item = spectral_generalized_maximum(matrices["Jfull"], matrices["I"], cutoff)
        cutoff_runs.append({key: value for key, value in item.items() if key != "vector"})
    chosen = spectral_generalized_maximum(matrices["Jfull"], matrices["I"], 1e-11)
    score = rho * float(chosen["quotient"])
    constants = rough_factor_constants(beta)
    contributions = {
        "plus_omega": score * constants["omega_choose_1"],
        "minus_2_choose_2": -2.0 * score * constants["omega_choose_2"],
        "plus_3_choose_3": 3.0 * score * constants["omega_choose_3"],
    }
    gross = sum(abs(value) for value in contributions.values())
    eta = max(0.0, score - 1.0)
    record = {
        "dimension": dimension,
        "intervals": intervals,
        "convolution_length": model.n,
        "fixed_published_k40_coefficients": fixed,
        "full_face_optimization": {
            "selected_cutoff": 1e-11,
            "score": score,
            "quotient_before_rho_star": float(chosen["quotient"]),
            "retained_dimension": chosen["retained_dimension"],
            "projected_relative_residual": chosen["projected_relative_residual"],
            "cutoff_sensitivity": cutoff_runs,
        },
        "factorial_contributions_normalized_by_I": contributions,
        "signed_sum": sum(contributions.values()),
        "gross_absolute_contribution": gross,
        "parity_error_budget": {
            "eta_max_for_abs_error_le_eta_I": eta,
            "common_relative_error_budget_if_each_unsigned_term_has_same_relative_bound": eta / gross if gross else 0.0,
        },
    }
    matrices["hybrid"] = hybrid_matrix
    matrices["selected_vector"] = np.asarray(chosen["vector"])
    return record, matrices


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimensions", type=int, nargs="+", default=(40, 39, 38))
    parser.add_argument("--intervals", type=int, default=2048)
    parser.add_argument("--beta", type=float, default=0.250001)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--matrix-output", type=Path)
    args = parser.parse_args()
    upstream = load_inputs(args.input)
    records = []
    saved = {}
    for dimension in args.dimensions:
        print(f"BUILD dimension={dimension} intervals={args.intervals}", flush=True)
        record, matrices = evaluate_dimension(upstream, dimension, args.intervals, args.beta)
        records.append(record)
        print(json.dumps(record, sort_keys=True), flush=True)
        if args.matrix_output is not None:
            for name, matrix in matrices.items():
                saved[f"k{dimension}_{name}"] = matrix
    matrix_receipt = None
    if args.matrix_output is not None:
        args.matrix_output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.matrix_output, **saved)
        matrix_receipt = {
            "path": str(args.matrix_output),
            "sha256": file_hash(args.matrix_output),
        }
    payload = {
        "schema": "primegaps.physical-parity-viability.v1",
        "status": "exploratory-float64-idealized-rough-and-parity-asymptotics",
        "upstream": {
            **upstream.UPSTREAM,
            "input_path": str(args.input),
            "input_sha256": file_hash(args.input),
        },
        "model": {
            "beta": args.beta,
            "rho_star": str(upstream.RHO_STAR),
            "rough_factor_constants": rough_factor_constants(args.beta),
            "assumptions": [
                "perfect x^beta-rough detector",
                "perfect asymptotics for the signed factorial-count combination",
                "full PrimeGaps186 face support with no source losses",
                "PrimeGaps186 outer cap shells and 77-dimensional trial basis retained",
                "dimension is varied while continuous radii, caps, profile and rho_star are frozen",
            ],
        },
        "dimensions": records,
        "matrix_receipt": matrix_receipt,
    }
    atomic_json(args.output, payload)
    print(json.dumps(payload, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()

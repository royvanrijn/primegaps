#!/usr/bin/env python3
"""Reduced finite test of the residue-coloured parity architecture.

The experiment keeps the actual frozen k=39 physical trial vector, but reduces
the arithmetic to three non-target colours and four small prime atoms.  It
compiles compatible divisor pairs

    (d, e) -> (q, CRT colour word) -> c_i(q, a),

then pairs the coefficient with exactly enumerated prime and rough-semiprime
progression discrepancies on (X, 2X].  This is a viability diagnostic, not an
asymptotic model or a theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import floor, gcd, log, prod, sqrt
from pathlib import Path
import tempfile

import numpy as np


def repository_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "pyproject.toml").is_file():
            return parent
    raise FileNotFoundError("could not locate repository root")


ROOT = repository_root()
PHYSICAL_INPUT = ROOT / "reproduction/186/physical-parity-input.json"
VECTOR_INPUT = ROOT / "reproduction/186/physical-parity-k39-vectors.json"
TARGET_SHIFT = 0
COLOUR_SHIFTS = (2, 36, 48)
PRIME_ATOMS = (5, 7, 11, 13)
PHYSICAL_LOG_SIZES = (0.18, 0.22, 0.26, 0.30)


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, prefix=path.name, delete=False
    ) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        temporary = Path(stream.name)
    temporary.replace(path)


def sieve(limit: int) -> np.ndarray:
    values = np.ones(limit + 1, dtype=np.bool_)
    values[:2] = False
    for prime in range(2, floor(sqrt(limit)) + 1):
        if values[prime]:
            values[prime * prime : limit + 1 : prime] = False
    return values


def integer_cuberoot(number: int) -> int:
    root = int(round(number ** (1 / 3)))
    while (root + 1) ** 3 <= number:
        root += 1
    while root**3 > number:
        root -= 1
    return root


def sector_values(x: int) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Return primes and z-rough semiprimes in (x,2x].

    We choose z=floor((2x)^(1/3)).  Requiring factors strictly larger than z
    makes three such factors exceed 2x, so the finite carrier has degree two
    without an asymptotic endpoint convention.
    """
    primality = sieve(2 * x)
    all_primes = np.flatnonzero(primality).astype(np.int64)
    primes = all_primes[(all_primes > x) & (all_primes <= 2 * x)]
    z = integer_cuberoot(2 * x)
    small_factors = all_primes[
        (all_primes > z) & (all_primes <= floor(sqrt(2 * x)))
    ]
    semiprime_parts: list[np.ndarray] = []
    for left in small_factors:
        lower = max(int(left), x // int(left) + 1, z + 1)
        upper = (2 * x) // int(left)
        lo = int(np.searchsorted(all_primes, lower, side="left"))
        hi = int(np.searchsorted(all_primes, upper, side="right"))
        if lo < hi:
            semiprime_parts.append(left * all_primes[lo:hi])
    semiprimes = (
        np.concatenate(semiprime_parts).astype(np.int64, copy=False)
        if semiprime_parts
        else np.empty(0, dtype=np.int64)
    )
    if semiprimes.size and (
        int(semiprimes.min()) <= x or int(semiprimes.max()) > 2 * x
    ):
        raise AssertionError("semiprime enumeration escaped the interval")
    return primes, semiprimes, {
        "x": x,
        "interval": [x + 1, 2 * x],
        "roughness_cutoff_z": z,
        "effective_beta_log_x_z": log(z) / log(x),
        "degree_two_finite_carrier": (z + 1) ** 3 > 2 * x,
        "prime_count": int(primes.size),
        "rough_semiprime_count": int(semiprimes.size),
    }


def decode(code: int, base: int, width: int) -> tuple[int, ...]:
    answer = []
    for _ in range(width):
        answer.append(code % base)
        code //= base
    return tuple(answer)


def encode(values: tuple[int, ...], base: int) -> int:
    answer = 0
    multiplier = 1
    for value in values:
        answer += value * multiplier
        multiplier *= base
    return answer


def crt(primes: tuple[int, ...], residues: tuple[int, ...]) -> int:
    if not primes:
        return 0
    modulus = prod(primes)
    return sum(
        residue
        * (modulus // prime)
        * pow(modulus // prime, -1, prime)
        for prime, residue in zip(primes, residues, strict=True)
    ) % modulus


@dataclass(frozen=True)
class PhysicalTrial:
    signatures: tuple[tuple[int, ...], ...]
    shells: tuple[tuple[float, float, float], ...]
    center: float
    rho: float
    vector: np.ndarray

    def support_cap(self, coordinates: tuple[float, ...]) -> float | None:
        total = sum(coordinates)
        maximum = max(coordinates, default=0.0)
        for index, (lower, upper, cap) in enumerate(self.shells):
            inside_lower = total > lower or (index == 0 and total == 0.0)
            if inside_lower and total <= upper + 2e-15 and maximum <= cap + 2e-15:
                return cap
        return None

    def evaluate(self, coordinates: tuple[float, ...]) -> float:
        if self.support_cap(coordinates) is None:
            return 0.0
        total = sum(coordinates)
        power_sums = {
            power: sum(value**power for value in coordinates)
            for signature in self.signatures
            for power in signature
        }
        value = 0.0
        index = 0
        for signature in self.signatures:
            angular = prod(power_sums[power] for power in signature)
            for degree in range(7):
                value += self.vector[index] * angular * (total - self.center) ** degree
                index += 1
        if index != self.vector.size:
            raise AssertionError((index, self.vector.size))
        return value


def load_trials() -> tuple[dict[str, PhysicalTrial], dict[str, object]]:
    payload = json.loads(PHYSICAL_INPUT.read_text(encoding="utf-8"))
    vector_payload = json.loads(VECTOR_INPUT.read_text(encoding="utf-8"))
    if vector_payload.get("schema") != "primegaps.physical-parity-k39-vectors.v1":
        raise ValueError("unexpected k39 vector input schema")
    shells = []
    lower = 0.0
    for upper, cap in payload["cap_shell_data"]["outer"]:
        shells.append((lower, float(Fraction(upper)), float(Fraction(cap))))
        lower = float(Fraction(upper))
    common = {
        "signatures": tuple(
            tuple(int(Fraction(power)) for power in signature)
            for signature in payload["signatures"]
        ),
        "shells": tuple(shells),
        "center": float(Fraction(payload["center"])),
        "rho": float(Fraction(payload["rho_star"])),
    }
    trials = {}
    vectors = {}
    for name, values in vector_payload["vectors"].items():
        vector = np.asarray(values, dtype=float)
        trials[name] = PhysicalTrial(**common, vector=vector)
        vectors[name] = {"length": int(vector.size)}
    return trials, {
        "physical_input": {
            "path": str(PHYSICAL_INPUT.relative_to(ROOT)),
            "sha256": file_sha256(PHYSICAL_INPUT),
        },
        "vector_input": {
            "path": str(VECTOR_INPUT.relative_to(ROOT)),
            "sha256": file_sha256(VECTOR_INPUT),
            "vectors": vectors,
            "source": vector_payload["source"],
        },
    }


def configuration_data(
    x: int, primes: tuple[int, ...], colours: int, trial: PhysicalTrial
) -> tuple[np.ndarray, np.ndarray, list[tuple[int, ...]], dict[str, object]]:
    base = colours + 1
    assignments = [decode(code, base, len(primes)) for code in range(base ** len(primes))]
    logarithmic_sizes = PHYSICAL_LOG_SIZES
    if len(logarithmic_sizes) != len(primes):
        raise AssertionError("physical mesh and prime atoms have different lengths")
    raw_values = np.zeros(len(assignments), dtype=float)
    mobius_signs = np.ones(len(assignments), dtype=float)
    supported = np.zeros(len(assignments), dtype=np.bool_)
    for code, assignment in enumerate(assignments):
        coordinates = tuple(
            sum(size for size, assigned in zip(logarithmic_sizes, assignment, strict=True) if assigned == colour)
            for colour in range(1, colours + 1)
        )
        value = trial.evaluate(coordinates)
        if value != 0.0:
            supported[code] = True
            raw_values[code] = value
            mobius_signs[code] = -1.0 if sum(item != 0 for item in assignment) % 2 else 1.0
    signed = raw_values * mobius_signs
    norm = float(np.linalg.norm(signed))
    if not norm:
        raise ArithmeticError("reduced physical trial vanished")
    lambdas = signed / norm
    return lambdas, supported, assignments, {
        "prime_atoms": list(primes),
        "coarse_physical_fragment_sizes": list(logarithmic_sizes),
        "actual_log_x_prime_exponents_used_only_for_envelope_checks": [
            log(prime) / log(x) for prime in primes
        ],
        "coordinate_colours": colours,
        "ambient_trial_dimension": 39,
        "active_reduced_colours": colours,
        "configuration_count": len(assignments),
        "supported_configuration_count": int(np.count_nonzero(supported)),
        "lambda_normalization": "sum_d Lambda_d^2 = 1 on the reduced configuration set",
        "pre_normalization_l2": norm,
    }


def compatible_pairs(
    assignments: list[tuple[int, ...]], supported: np.ndarray, colours: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Enumerate the 1+3K compatible local states at every prime."""
    width = len(assignments[0])
    base = colours + 1
    left_codes: list[int] = []
    right_codes: list[int] = []
    union_codes: list[int] = []
    for states in product(range(1 + 3 * colours), repeat=width):
        left = []
        right = []
        union = []
        for state in states:
            if state == 0:
                left.append(0)
                right.append(0)
                union.append(0)
                continue
            offset = state - 1
            colour = offset // 3 + 1
            mode = offset % 3
            left.append(colour if mode in (0, 2) else 0)
            right.append(colour if mode in (1, 2) else 0)
            union.append(colour)
        left_code = encode(tuple(left), base)
        right_code = encode(tuple(right), base)
        if not (supported[left_code] and supported[right_code]):
            continue
        left_codes.append(left_code)
        right_codes.append(right_code)
        union_codes.append(encode(tuple(union), base))
    return (
        np.asarray(left_codes, dtype=np.int32),
        np.asarray(right_codes, dtype=np.int32),
        np.asarray(union_codes, dtype=np.int32),
    )


def state_arithmetic(
    assignments: list[tuple[int, ...]],
    primes: tuple[int, ...],
    colour_residues: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    moduli = np.ones(len(assignments), dtype=np.int64)
    residues = np.zeros(len(assignments), dtype=np.int64)
    masks = np.zeros(len(assignments), dtype=np.int64)
    totients = np.ones(len(assignments), dtype=np.int64)
    for code, assignment in enumerate(assignments):
        selected_primes = tuple(
            prime for prime, colour in zip(primes, assignment, strict=True) if colour
        )
        selected_residues = tuple(
            colour_residues[colour - 1] for colour in assignment if colour
        )
        mask = sum(1 << index for index, colour in enumerate(assignment) if colour)
        moduli[code] = prod(selected_primes)
        residues[code] = crt(selected_primes, selected_residues)
        masks[code] = mask
        totients[code] = prod(prime - 1 for prime in selected_primes)
        if gcd(int(moduli[code]), int(residues[code])) != 1 and moduli[code] != 1:
            raise AssertionError((code, assignment, moduli[code], residues[code]))
    return moduli, residues, masks, totients


def progression_errors(
    values: np.ndarray,
    moduli: np.ndarray,
    residues: np.ndarray,
    totients: np.ndarray,
) -> np.ndarray:
    errors = np.zeros(moduli.size, dtype=float)
    cached_counts: dict[int, np.ndarray] = {}
    total = int(values.size)
    for code, (modulus, residue, phi) in enumerate(zip(moduli, residues, totients, strict=True)):
        q = int(modulus)
        if q == 1:
            continue
        if q not in cached_counts:
            cached_counts[q] = np.bincount(values % q, minlength=q)
        errors[code] = float(cached_counts[q][int(residue)]) - total / int(phi)
    return errors


def safe_ratio(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator else None


def singular_summary(matrix: np.ndarray) -> dict[str, object]:
    singular = np.linalg.svd(matrix, compute_uv=False)
    squares = singular * singular
    energy = float(squares.sum())
    cumulative = np.cumsum(squares) / energy if energy else np.zeros_like(squares)
    captures = {}
    for rank in (1, 2, 4, 8, 16, 32):
        if rank <= singular.size:
            captures[str(rank)] = float(cumulative[rank - 1])
    top = float(singular[0]) if singular.size else 0.0
    return {
        "shape": list(matrix.shape),
        "frobenius_squared": energy,
        "singular_values": [float(value) for value in singular],
        "singular_values_relative_to_first": [
            float(value / top) if top else 0.0 for value in singular
        ],
        "cumulative_frobenius_energy": captures,
        "stable_rank": energy / (top * top) if top else 0.0,
        "relative_rank_at_1e-6": int(np.count_nonzero(singular > 1e-6 * top)) if top else 0,
    }


def global_colour_table(
    values: np.ndarray,
    assignments: list[tuple[int, ...]],
    masks: np.ndarray,
    colours: int,
) -> np.ndarray:
    """Isometric padding of ragged c(q,chi) fibres to global colour words."""
    width = len(assignments[0])
    base = colours + 1
    words = tuple(product(range(1, colours + 1), repeat=width))
    table = np.zeros((1 << width, len(words)), dtype=float)
    for column, word in enumerate(words):
        for mask in range(1 << width):
            restricted = tuple(
                word[index] if mask & (1 << index) else 0 for index in range(width)
            )
            code = encode(restricted, base)
            active = int(mask).bit_count()
            table[mask, column] = values[code] / sqrt(colours ** (width - active))
            if int(masks[code]) != mask:
                raise AssertionError((code, mask, masks[code]))
    return table


def quantiles(values: np.ndarray) -> dict[str, float]:
    if not values.size:
        return {}
    return {
        label: float(np.quantile(values, point))
        for label, point in (("min", 0), ("p25", 0.25), ("median", 0.5), ("p75", 0.75), ("max", 1))
    }


def analyse_trial(
    lambdas: np.ndarray,
    pair_left: np.ndarray,
    pair_right: np.ndarray,
    pair_union: np.ndarray,
    assignments: list[tuple[int, ...]],
    moduli: np.ndarray,
    residues: np.ndarray,
    masks: np.ndarray,
    prime_errors: np.ndarray,
    semiprime_errors: np.ndarray,
    colours: int,
) -> dict[str, object]:
    raw_weights = lambdas[pair_left] * lambdas[pair_right]
    state_count = len(assignments)
    coefficients = np.bincount(pair_union, weights=raw_weights, minlength=state_count)
    raw_l1_by_state = np.bincount(
        pair_union, weights=np.abs(raw_weights), minlength=state_count
    )
    collisions = np.bincount(pair_union, minlength=state_count)
    grouped_l1 = float(np.abs(coefficients).sum())
    raw_l1 = float(np.abs(raw_weights).sum())
    coefficient_l2_squared = float(coefficients @ coefficients)
    occupied = collisions > 0
    material = np.abs(coefficients) > max(1e-14 * float(np.max(np.abs(coefficients))), 1e-300)
    cancellation_per_state = raw_l1_by_state[material] / np.abs(coefficients[material])

    test_integers = sorted(
        set(range(4096))
        | {int(residues[code]) for code in np.flatnonzero(occupied)}
        | {
            int(residues[code] + moduli[code])
            for code in np.flatnonzero(occupied)
        }
    )
    compiler_max_error = 0.0
    active_lambda = np.flatnonzero(lambdas)
    active_coefficients = np.flatnonzero(occupied)
    for integer in test_integers:
        divisor_sum = float(
            lambdas[
                active_lambda[
                    integer % moduli[active_lambda] == residues[active_lambda]
                ]
            ].sum()
        )
        coefficient_sum = float(
            coefficients[
                active_coefficients[
                    integer % moduli[active_coefficients]
                    == residues[active_coefficients]
                ]
            ].sum()
        )
        compiler_max_error = max(
            compiler_max_error, abs(divisor_sum * divisor_sum - coefficient_sum)
        )
    if compiler_max_error > 2e-12:
        raise AssertionError(f"coefficient compiler error {compiler_max_error:.3e}")

    raw_prime = raw_weights * prime_errors[pair_union]
    raw_semiprime = raw_weights * semiprime_errors[pair_union]
    grouped_prime = coefficients * prime_errors
    grouped_semiprime = coefficients * semiprime_errors
    total_prime = float(grouped_prime.sum())
    total_semiprime = float(grouped_semiprime.sum())
    projected_total = total_prime - 2 * total_semiprime
    denominator = abs(projected_total)

    raw_gross = float(np.abs(raw_prime).sum() + 2 * np.abs(raw_semiprime).sum())
    raw_projected_l1 = float(np.abs(raw_prime - 2 * raw_semiprime).sum())
    grouped_gross = float(np.abs(grouped_prime).sum() + 2 * np.abs(grouped_semiprime).sum())
    grouped_projected_l1 = float(np.abs(grouped_prime - 2 * grouped_semiprime).sum())

    coefficient_table = global_colour_table(coefficients, assignments, masks, colours)
    projected_table = global_colour_table(
        coefficients * (prime_errors - 2 * semiprime_errors),
        assignments,
        masks,
        colours,
    )
    table_energy = float(np.square(coefficient_table).sum())
    if not np.isclose(table_energy, coefficient_l2_squared, rtol=2e-13, atol=1e-30):
        raise AssertionError((table_energy, coefficient_l2_squared))

    dyadic = []
    for exponent in sorted({int(floor(log(int(q), 2))) if q > 1 else 0 for q in moduli[occupied]}):
        lower = 1 << exponent
        upper = 1 << (exponent + 1)
        selected = occupied & (moduli >= lower) & (moduli < upper)
        if not np.any(selected):
            continue
        raw_bin = float(raw_l1_by_state[selected].sum())
        grouped_bin = float(np.abs(coefficients[selected]).sum())
        dyadic.append(
            {
                "q_interval": [lower, upper],
                "raw_pair_count": int(collisions[selected].sum()),
                "occupied_crt_states": int(np.count_nonzero(selected)),
                "sum_abs_raw_coefficients": raw_bin,
                "sum_abs_grouped_coefficients": grouped_bin,
                "sum_grouped_abs_squared": float(np.square(coefficients[selected]).sum()),
                "crt_collision_cancellation_ratio": safe_ratio(raw_bin, grouped_bin),
            }
        )

    exact_scalar_ratio = safe_ratio(abs(total_prime) + 2 * abs(total_semiprime), denominator)
    return {
        "coefficient_metrics": {
            "coefficient_compiler_identity": {
                "statement": "(sum_d Lambda_d 1_(d constraints))(m)^2 = sum_(q,a) c(q,a) 1_(m=a mod q)",
                "tested_integers": len(test_integers),
                "maximum_absolute_error": compiler_max_error,
            },
            "raw_compatible_pair_terms": int(raw_weights.size),
            "occupied_crt_states_before_coefficient_cancellation": int(np.count_nonzero(occupied)),
            "material_nonzero_crt_coefficients": int(np.count_nonzero(material)),
            "sum_abs_raw_pair_coefficients": raw_l1,
            "sum_abs_crt_aggregated_coefficients": grouped_l1,
            "sum_q_a_abs_c_squared": coefficient_l2_squared,
            "crt_collision_cancellation_ratio": safe_ratio(raw_l1, grouped_l1),
            "raw_pairs_per_occupied_crt_state": safe_ratio(int(raw_weights.size), int(np.count_nonzero(occupied))),
            "per_state_raw_l1_over_abs_c_quantiles": quantiles(cancellation_per_state),
            "dyadic_q": dyadic,
        },
        "parity_cancellation": {
            "definition": "For blocks b, R=(sum_b |P_b|+2|S_b|)/|sum_b(P_b-2S_b)|. Raw blocks are compatible (d,e); aggregated blocks are (q,a).",
            "E_prime": total_prime,
            "E_semiprime": total_semiprime,
            "E_prime_minus_2_E_semiprime": projected_total,
            "exact_scalar_ratio": exact_scalar_ratio,
            "scalar_ratio_invariant_under_exact_crt_aggregation": True,
            "before_crt_aggregation": {
                "gross": raw_gross,
                "projected_block_l1": raw_projected_l1,
                "cancellation_ratio": safe_ratio(raw_gross, denominator),
                "within_block_prime_semiprime_ratio": safe_ratio(raw_gross, raw_projected_l1),
                "across_block_ratio": safe_ratio(raw_projected_l1, denominator),
            },
            "after_crt_aggregation": {
                "gross": grouped_gross,
                "projected_block_l1": grouped_projected_l1,
                "cancellation_ratio": safe_ratio(grouped_gross, denominator),
                "within_block_prime_semiprime_ratio": safe_ratio(grouped_gross, grouped_projected_l1),
                "across_block_ratio": safe_ratio(grouped_projected_l1, denominator),
            },
            "gross_reduction_from_crt_aggregation": safe_ratio(raw_gross, grouped_gross),
        },
        "singular_value_decay": {
            "coefficient_table": singular_summary(coefficient_table),
            "projected_error_table": singular_summary(projected_table),
            "table_definition": "rows are physical prime-subset moduli; columns are full colour words; each restricted fibre is padded by K^(inactive/2), so Frobenius^2=sum_(q,a)|c(q,a)|^2",
        },
    }


def run_one_scale(
    x: int,
    trials: dict[str, PhysicalTrial],
    colours: int,
) -> dict[str, object]:
    prime_values, semiprime_values, sector = sector_values(x)
    first_trial = next(iter(trials.values()))
    primes = PRIME_ATOMS
    one_prime_envelope = first_trial.rho * float(Fraction(190370, 262499))
    pair_modulus_envelope = 2 * first_trial.rho * float(Fraction(2742997, 2624989))
    if max(log(prime) / log(x) for prime in primes) > one_prime_envelope:
        raise ValueError(f"x={x} violates the one-prime physical envelope")
    if log(prod(primes)) / log(x) > pair_modulus_envelope:
        raise ValueError(f"x={x} violates the generated pair-modulus envelope")

    # Support depends only on the shell, not on the trial vector.  Build the
    # combinatorial compiler once using the first vector, then reevaluate the
    # lambda values for every mesh vector.
    _, support, assignments, common_model = configuration_data(
        x, primes, colours, first_trial
    )
    pair_left, pair_right, pair_union = compatible_pairs(assignments, support, colours)
    colour_residues = tuple(TARGET_SHIFT - shift for shift in COLOUR_SHIFTS)
    moduli, residues, masks, totients = state_arithmetic(
        assignments, primes, colour_residues
    )
    prime_errors = progression_errors(
        prime_values, moduli, residues, totients
    )
    semiprime_errors = progression_errors(
        semiprime_values, moduli, residues, totients
    )

    results = {}
    trial_models = {}
    for name, trial in trials.items():
        lambdas, this_support, this_assignments, model = configuration_data(
            x, primes, colours, trial
        )
        if not np.array_equal(this_support, support) or this_assignments != assignments:
            raise AssertionError("trial vectors changed the physical support")
        trial_models[name] = model
        results[name] = analyse_trial(
            lambdas,
            pair_left,
            pair_right,
            pair_union,
            assignments,
            moduli,
            residues,
            masks,
            prime_errors,
            semiprime_errors,
            colours,
        )

    return {
        "sector_enumeration": sector,
        "common_reduced_model": {
            **common_model,
            "target_shift": TARGET_SHIFT,
            "colour_shifts": list(COLOUR_SHIFTS),
            "target_minus_colour_residue_integers": list(colour_residues),
            "local_residues_by_prime": {
                str(prime): [value % prime for value in colour_residues]
                for prime in primes
            },
            "compatible_pair_terms_after_support": int(pair_left.size),
            "maximum_generated_modulus": int(moduli[np.bincount(pair_union, minlength=len(assignments)) > 0].max()),
            "one_prime_exponent_envelope": one_prime_envelope,
            "pair_modulus_exponent_envelope": pair_modulus_envelope,
            "actual_maximum_prime_exponent": max(log(prime) / log(x) for prime in primes),
            "actual_full_atom_modulus_exponent": log(prod(primes)) / log(x),
        },
        "trial_normalizations": trial_models,
        "trials": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--x", type=int, nargs="+", default=(8_000_000, 16_000_000, 32_000_000, 64_000_000))
    parser.add_argument("--colours", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.colours != 3:
        raise ValueError("the recorded experiment fixes three reduced colours")
    trials, inputs = load_trials()
    payload = {
        "schema": "primegaps.parity-crt-finite.v1",
        "status": "exploratory-reduced-finite-model-not-an-asymptotic-or-theorem",
        "inputs": inputs,
        "assumptions": [
            "three non-target colours use the actual H39 shifts 2,36,48 relative to target shift 0",
            "four CRT prime atoms 5,7,11,13 and no presieving limit; the chosen three residues are nonzero and distinct at every atom",
            "coarse physical fragment sizes 0.18,0.22,0.26,0.30 are attached independently of the small CRT-prime labels",
            "k=39 physical trial restricted to three active coordinates; the other 36 coordinates are zero",
            "divisor pairs are compatible only when a shared prime has the same colour",
            "progression errors are count(q,a)-total/phi(q), centered separately in the prime and rough-semiprime sectors",
            "z=floor((2X)^(1/3)) and both semiprime factors exceed z",
            "Lambda is normalized to unit discrete l2 separately at each X and mesh vector",
        ],
        "scales": [run_one_scale(x, trials, args.colours) for x in args.x],
    }
    atomic_json(args.output, payload)
    print(json.dumps({
        "schema": payload["schema"],
        "output": str(args.output),
        "output_sha256": file_sha256(args.output),
        "scales": [item["sector_enumeration"]["x"] for item in payload["scales"]],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

"""Deterministic CRT and rational reconstruction for exact contractions."""

from __future__ import annotations

from math import gcd, isqrt


def crt_pair(a: int, modulus_a: int, b: int, modulus_b: int):
    """Combine two residues with coprime positive moduli."""
    if modulus_a <= 0 or modulus_b <= 0:
        raise ValueError("CRT moduli must be positive")
    if gcd(modulus_a, modulus_b) != 1:
        raise ValueError("CRT moduli must be coprime")
    adjustment = ((b - a) * pow(modulus_a, -1, modulus_b)) % modulus_b
    modulus = modulus_a * modulus_b
    return (a + modulus_a * adjustment) % modulus, modulus


def crt(residues, moduli):
    residues = tuple(int(value) for value in residues)
    moduli = tuple(int(value) for value in moduli)
    if not residues or len(residues) != len(moduli):
        raise ValueError("supply equally many non-empty residues and moduli")
    value, modulus = residues[0] % moduli[0], moduli[0]
    for residue, next_modulus in zip(residues[1:], moduli[1:]):
        value, modulus = crt_pair(value, modulus, residue, next_modulus)
    return value, modulus


def rational_residue(numerator: int, denominator: int, prime: int) -> int:
    """Map a rational to F_p, rejecting a non-invertible denominator."""
    denominator %= prime
    if not denominator:
        raise ZeroDivisionError("rational denominator vanishes modulo prime")
    return (numerator % prime) * pow(denominator, -1, prime) % prime


def rational_reconstruction(
    residue: int,
    modulus: int,
    *,
    numerator_bound: int | None = None,
    denominator_bound: int | None = None,
):
    """Recover n/d from residue modulo modulus under explicit uniqueness bounds.

    If bounds are omitted, both use floor(sqrt(modulus/2)). A successful result
    is checked for coprimality, positive denominator, bounds, and congruence.
    """
    if modulus <= 1:
        raise ValueError("modulus must exceed one")
    residue %= modulus
    default_bound = isqrt(modulus // 2)
    numerator_bound = (
        default_bound if numerator_bound is None else int(numerator_bound)
    )
    denominator_bound = (
        default_bound if denominator_bound is None else int(denominator_bound)
    )
    if numerator_bound < 0 or denominator_bound < 1:
        raise ValueError("invalid reconstruction bounds")
    if 2 * numerator_bound * denominator_bound >= modulus:
        raise ValueError("bounds do not guarantee unique reconstruction")

    old_r, r = modulus, residue
    old_t, t = 0, 1
    while r and abs(r) > numerator_bound:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_t, t = t, old_t - quotient * t
    numerator, denominator = r, t
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    common = gcd(abs(numerator), denominator)
    if common:
        numerator //= common
        denominator //= common
    if (
        abs(numerator) > numerator_bound
        or not 0 < denominator <= denominator_bound
        or gcd(abs(numerator), denominator) != 1
        or (numerator - residue * denominator) % modulus
    ):
        raise ValueError("no rational satisfies the reconstruction bounds")
    return numerator, denominator


def is_prime_64(value: int) -> bool:
    """Deterministic Miller-Rabin primality test for unsigned 64-bit integers."""
    value = int(value)
    if value < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for prime in small:
        if value == prime:
            return True
        if value % prime == 0:
            return False
    exponent = value - 1
    shifts = 0
    while exponent % 2 == 0:
        shifts += 1
        exponent //= 2
    for base in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if base % value == 0:
            continue
        witness = pow(base, exponent, value)
        if witness in (1, value - 1):
            continue
        for _ in range(shifts - 1):
            witness = witness * witness % value
            if witness == value - 1:
                break
        else:
            return False
    return True


def descending_primes(start: int = (1 << 61) - 1):
    """Yield deterministic large primes suitable for CRT."""
    candidate = int(start) | 1
    while candidate >= 3:
        if is_prime_64(candidate):
            yield candidate
        candidate -= 2

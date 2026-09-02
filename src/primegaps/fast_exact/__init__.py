"""Accelerated exact arithmetic for the reproduced 240 certificate.

The frozen verifier under :mod:`reproduction/240` remains the trust oracle.
This package supplies independently tested accelerators: closed-zero density
recurrences, pair-first J contraction, FLINT polynomial kernels, modular/CRT
arithmetic, and append-only candidate-independent moment caches.
"""

from .fast_i import (
    apply_zero_block,
    orbit_status_densities,
    positive_status_density,
    signature_moments,
    signature_value,
)
from .modular_exact import (
    crt,
    descending_primes,
    rational_reconstruction,
    rational_residue,
)
from .moment_cache import IMomentCache, JFunctionalCache

__all__ = [
    "IMomentCache",
    "JFunctionalCache",
    "apply_zero_block",
    "crt",
    "descending_primes",
    "orbit_status_densities",
    "positive_status_density",
    "rational_reconstruction",
    "rational_residue",
    "signature_moments",
    "signature_value",
]

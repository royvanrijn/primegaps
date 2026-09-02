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
from .j_block import (
    JBlockOperator,
    MarginalMap,
    accumulate_candidate_gram,
    accumulate_feature_gram_blocks,
    accumulate_gram_difference,
    candidate_feature_values,
    factorized_feature_values,
    load_block_operator,
    projected_feature_values,
    save_block_operator,
)

__all__ = [
    "IMomentCache",
    "JFunctionalCache",
    "JBlockOperator",
    "MarginalMap",
    "accumulate_candidate_gram",
    "accumulate_feature_gram_blocks",
    "accumulate_gram_difference",
    "candidate_feature_values",
    "factorized_feature_values",
    "load_block_operator",
    "save_block_operator",
    "apply_zero_block",
    "crt",
    "descending_primes",
    "orbit_status_densities",
    "positive_status_density",
    "projected_feature_values",
    "rational_reconstruction",
    "rational_residue",
    "signature_moments",
    "signature_value",
]

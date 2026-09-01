# primegaps

Exploratory code around Julia Stadlmann's 2026 bounded-prime-gaps paper
[arXiv:2608.31126](https://arxiv.org/abs/2608.31126), which improves the unconditional
bound from `H_1 <= 246` to `H_1 <= 240`.

The immediate goal is to reproduce the paper's `k=49`, degree-21 generalized
eigenvalue certificate and then probe `k=48,47,...`. This repository contains
independently checkable layers of that reconstruction:

- the published support `T_k(delta, A, B, epsilon)`;
- generation/counting of the symmetric-polynomial basis used by Polymath/Stadlmann;
- reconstructed exact rational `C_{m,i}`/`D_{m,i}` coefficient formulas;
- symmetry-compressed factorial moments and sparse simplex certificate replay,
  cross-pollinated from Axiom Math's PrimeGapsLib;
- a low-dimensional exact `I/J/K` matrix API for sparse polynomial bases;
- a geometry diagnostic for support sampling;
- a numerical generalized-eigenvalue and exact-certificate proof backend;
- an exact-rational distribution-region oracle for support-cell pairs;
- notes on what is and is not yet reproduced.

**Nothing in this repository currently proves a bound below 240.** The exact
integral kernel is now independently testable, but symmetry-compressed `k=49`
matrix assembly and the published certificate are not yet reproduced.

See [the Section 5 engine notes](docs/section5-engine.md) for the reconstructed
formulas, source errata, exact API, and verification boundary.

See [the PrimeGapsLib cross-pollination notes](docs/primegapslib-crosspollination.md)
for the monomial-signature representation, nilpotent factorial-moment ladder,
closed simplex/marginal formulas, sparse grouped contraction, and exact replay
of the released `M_{50,1/25}>4` source certificate.

## Proof backend

Once exact symmetric matrices are supplied, the backend finds a numerical
maximizer, rationalizes it to a primitive integer vector, and independently
checks

    c^T M2 c > c^T M1 c

with Python big-integer arithmetic:

    primegaps-proof solve matrices.json --certificate certificate.json
    primegaps-proof verify matrices.json certificate.json

The first command can use a complete dense solve, memory-reduced Lanczos, exact
block decomposition, or optional SciPy sparse solving. The second command never
runs an eigensolver: it checks semantic matrix hashes and recomputes both exact
quadratic forms, their difference, and their quotient. See
[docs/proof_backend.md](docs/proof_backend.md) for formats, diagnostics, and the
degree-21 through degree-30 memory envelope.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
python -m primegaps.scan
pytest
```

## Distribution-region oracle

`primegaps.is_certified(region_a, region_b, minorant)` checks the
Bombieri--Vinogradov range first and otherwise applies the hypotheses of
Propositions 2 and 3.  A positive result includes exact rational checks for all
global Type I/II/III inequalities and universal partition witnesses for the
continuous `Xi(B_a,B_b,m_a,m_b,delta)` cell; it does not rely on sampling.

```python
from primegaps import Minorant, RegionCell, is_certified

cell = RegionCell(a_upper="0.253", large_count=3,
                  large_sum_bound="0.17", delta="0.028")
certificate = is_certified(cell, cell, Minorant("0.38", "0.4", "0.4"))

assert certificate.certified
print(certificate.theorem)
print(certificate.as_dict())
```

The oracle is sound but conservative: a negative answer says that these
implemented theorem witnesses do not certify the pair, not that every possible
partition argument must fail.  Decimal strings are recommended for exact input.

The scan prints the published parameters, basis dimensions for degree 21 and 27,
and a Monte-Carlo support-geometry diagnostic. The latter is deliberately *not* a
sieve score; it is useful mainly for rejecting naive "maximize support volume"
strategies.

# primegaps

Exploratory code around Julia Stadlmann's 2026 bounded-prime-gaps paper
[arXiv:2608.31126](https://arxiv.org/abs/2608.31126), which improves the unconditional
bound from `H_1 <= 246` to `H_1 <= 240`.

The paper's `k=49`, degree-21 variational certificate has now been reproduced
with exact rational arithmetic.  For one fixed 846-term rational polynomial,
the independently assembled value is

    49 J(F) / I(F) = 1.0011632465949216560417861678682244509240906847502...

With the same support and degree at `k=48`, the exact baseline is

    48 J(F) / I(F) = 0.9969233513526357503888760066573328995217614432838...
    1 - 48 J(F) / I(F) = 0.0030766486473642496111239933426671004782385567162...

This repository contains independently checkable layers of the reconstruction:

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

This reproduces the variational certificate used for `H_1 <= 240`; it does not
claim a bound below 240. See [the exact 240 reproduction](docs/reproduction-240.md)
for conventions, commands, hashes, independent checks, and the `k=48` status.

See [the faster exact backend](docs/fast-exact-backend.md) for the independently
checked closed-zero recurrence, pair-first J assembly, FLINT and modular
polynomial kernels, reusable moment caches, benchmarks, and checkpointed run
commands. The original reproduction remains the trust oracle.

See [the boundary-only Arb J certifier](docs/boundary-j-arb-certifier.md) for
the outward-rounded legal-minus-unrestricted scalar verifier, its frozen D21
oracle calibration, and the detached D27 checkpoint workflow.

See [the Gap 236 formalization crosswalk](docs/formalization-gap236.md) for the
pinned AxiomMath/PrimeGapsLib theorem architecture, the explicit admissible
48-tuple of diameter 236, and the exact separation between the reusable
finite/endgame proof and the new shaped-support/distribution `DHL48` obligation.

See [the Section 5 engine notes](docs/section5-engine.md) for the reconstructed
formulas, source errata, exact API, and verification boundary.

See [the `P3.II.delta` frontier](docs/p3ii-delta-frontier.md) for the primary
analytic counterfactual target, the full numerical `A_max -> lambda_48` curve,
the estimated crossing of 1, and the next binding constraint.

See [the Type-IIc incomplete-rectangle theorem](docs/typeiic-incomplete-rectangles.md)
for the checked analytic saving on the originating rectangles, its safe
endpoint `A=2029/8000`, and the precise verification boundary. The production
distribution oracle remains unchanged pending a fully typeset human review.

See [the research status audit](docs/research-status.md) for a ledger-wide
summary of exact results, numerical screens, negative experiments, superseded
work, and the present blockers. In particular, no degree or legal-support
experiment currently supplies an exact `k=48` certificate.

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

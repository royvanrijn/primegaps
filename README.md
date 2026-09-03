# primegaps

Research code for reproducing and improving unconditional bounded-prime-gap results.
The repository started from Julia Stadlmann's 2026 `H_1 <= 240` construction and now
also tracks the much stronger physical-fragment approach released in
[`openai/PrimeGaps186`](https://github.com/openai/PrimeGaps186).

## Current frontier

`openai/PrimeGaps186` supplies a conditional Lean development of

```text
DHL[40,2] -> H_1 <= 186
```

and a directed-Arb Python certificate. This is the new primary research basis, but the
public package is **not yet treated here as a closed unconditional record**: its final Lean
theorems depend on two finite-field estimates and one large numerical-integral axiom. The
first two appear to match established Katz/Deligne and Fouvry--Kowalski--Michel bounds; the
149 physical source bounds plus three cap bounds still need an independent fresh replay and
a small verified certificate bridge.

The immediate verification target is therefore to close/replay 186. The immediate frontier
target is

```text
DHL[39,2] -> H_1 <= 182.
```

See [`docs/primegaps186-impact.md`](docs/primegaps186-impact.md) for the construction,
conditionality audit, exact margin/loss ledger, integration plan and proposed `k=39`
experiments.

The September 2026 Axiom Math draft proving `H_1 <= 212` remains an important independently
checkable baseline. Its exact parameters and the relationship with our older work are recorded
in [`docs/bgp212-impact.md`](docs/bgp212-impact.md).

## Established work in this repository

Stadlmann's `k=49`, degree-21 variational certificate has been reproduced with exact rational
arithmetic. For the fixed 846-term rational polynomial,

```text
49 J(F) / I(F)
  = 1.0011632465949216560417861678682244509240906847502...
```

The same published construction at `k=48` gives the exact deficit

```text
48 J(F) / I(F)
  = 0.9969233513526357503888760066573328995217614432838...
```

The repository includes:

- exact Stadlmann support and symmetric-basis reconstruction;
- exact rational `C_{m,i}` / `D_{m,i}` coefficient formulas;
- symmetry-compressed simplex and marginal moments;
- an exact reproduction of the released `H_1 <= 240` variational witness;
- fast exact, FLINT and reusable-moment backends;
- projected/matrix-free numerical optimization;
- outward-rounded boundary-only Arb certification;
- exact distribution-region and minorant feasibility tools;
- shadow-price, minimum-breakthrough and sparse-SOS experiments;
- exact admissible-tuple replay scripts.

The pre-186 `k=48`, degree-27 Arb run may still be completed as an independent validation of
that discovery/certification pipeline, but it is no longer a record attempt.

## PrimeGaps186 replay utilities

Exact finite/rational checks that do not require the expensive physical-integral computation:

```bash
python scripts/verify_h40_tuple.py
python scripts/check_primegaps186_margin.py
```

The first verifies both the published 40-tuple of diameter 186 and its admissible 39-element
prefix of diameter 182. The second replays Corollary 2.6's exact rational margin and identifies
the dominant source-loss classes. Neither substitutes for running the upstream Arb certificate.

## Key documents

- [`docs/primegaps186-impact.md`](docs/primegaps186-impact.md): audit and research pivot to 186/182.
- [`docs/bgp212-impact.md`](docs/bgp212-impact.md): the 212 construction and remaining additive ideas.
- [`docs/reproduction-240.md`](docs/reproduction-240.md): exact Stadlmann `H_1 <= 240` replay.
- [`docs/boundary-j-arb-certifier.md`](docs/boundary-j-arb-certifier.md): scalar rigorous boundary verifier.
- [`docs/fast-exact-backend.md`](docs/fast-exact-backend.md): accelerated exact assembly and caches.
- [`docs/typeiic-incomplete-rectangles.md`](docs/typeiic-incomplete-rectangles.md): experimental Type-IIc work.
- [`docs/minimum-breakthrough.md`](docs/minimum-breakthrough.md): coupled analytic slack diagnostics.
- [`docs/research-status.md`](docs/research-status.md): historical ledger through the pre-186 phase.
- [`FINDINGS.md`](FINDINGS.md): detailed reconstruction notes and negative experiments.

## Proof backend

Given exact symmetric matrices, the backend finds a numerical maximizer, rationalizes it and
independently verifies the strict quadratic-form inequality:

```bash
primegaps-proof solve matrices.json --certificate certificate.json
primegaps-proof verify matrices.json certificate.json
```

The verifier does not run an eigensolver; it checks semantic matrix hashes and recomputes the
exact quadratic forms using integer arithmetic.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
python -m primegaps.scan
pytest
```

## Historical distribution oracle

`primegaps.is_certified(region_a, region_b, minorant)` implements the Stadlmann Proposition 2/3
feasibility conditions using exact rational checks and continuous order-statistic witnesses. It
remains useful for the older support family, but it does not model PrimeGaps186's physical
fragment/source-ladder construction.

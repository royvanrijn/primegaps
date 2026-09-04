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

An idealized parity extension of the physical state has now passed its numerical viability gate.
The frozen full-face model extrapolates to `1.0192774814` for `k=39` and `1.0127855640` for
`k=38`. Sweeping the roughness threshold shows that `beta=1/3+epsilon` is the natural first
analytic target: the exact detector drops from degree three to degree two, its cancellation
condition falls from `7.1610` to `3.7726`, and the triprime bilinear block disappears. This is
not a theorem. A subsequent dyadic audit refutes the proposed Friedlander--Iwaniec `B_F`
route for the already rough sequence: its active inner variable is forced prime, so the outer
absolute value destroys the needed Mobius cancellation. The direct signed Liouville target
remains open and now requires a decomposition that preserves factor-count cancellation globally
after CRT-coloured coefficient aggregation. See
[`docs/physical-parity-viability.md`](docs/physical-parity-viability.md) and
[`docs/bf-dyadic-audit.md`](docs/bf-dyadic-audit.md).

Artificially cutting the full-face operator by its generated outer-plus-inner
modulus exponent puts the first scanned `k=39` crossing at about `theta=0.5176`:
`Lambda_39(0.517625)=0.9999536` and
`Lambda_39(0.517750)=1.0000968`.  Thus this idealized model needs more than
`1/2`, but not the full `0.5485994` envelope.

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
python scripts/check_physical_parity_viability.py
python scripts/check_physical_parity_modulus_reach.py
python scripts/sweep_physical_parity_beta.py
python scripts/check_bf_dyadic_audit.py
```

The first verifies both the published 40-tuple of diameter 186 and its admissible 39-element
prefix of diameter 182. The second replays Corollary 2.6's exact rational margin and identifies
the dominant source-loss classes. The final four replay the idealized physical-parity score,
locate its generated-modulus crossing, map its beta-dependent cancellation tradeoff, and check
the obstruction to the proposed asymptotic-sieve route. None substitutes for running the
upstream Arb certificate or proving a
valid replacement for the direct rough weighted Liouville target.

## Key documents

- [`docs/primegaps186-impact.md`](docs/primegaps186-impact.md): audit and research pivot to 186/182.
- [`docs/physical-parity-viability.md`](docs/physical-parity-viability.md): roughness-threshold
  sweep, `k=39/38` viability, the direct signed parity target and sector-coupled alternatives.
- [`docs/bf-dyadic-audit.md`](docs/bf-dyadic-audit.md): complete indexed block family and
  refutation of the proposed asymptotic-sieve `B_F` route.
- [`docs/bgp212-impact.md`](docs/bgp212-impact.md): the 212 construction and remaining additive ideas.
- [`docs/formalization-gap236.md`](docs/formalization-gap236.md): reusable historical Lean endgame architecture.
- [`paper/gap236/README.md`](paper/gap236/README.md): conditional historical Gap236 manuscript and proof boundaries.
- [`docs/reproduction-240.md`](docs/reproduction-240.md): exact Stadlmann `H_1 <= 240` replay.
- [`docs/boundary-j-arb-certifier.md`](docs/boundary-j-arb-certifier.md): scalar rigorous boundary verifier.
- [`docs/fast-exact-backend.md`](docs/fast-exact-backend.md): accelerated exact assembly and caches.
- [`docs/section5-engine.md`](docs/section5-engine.md): reconstructed distribution formulas and exact API.
- [`docs/p3ii-delta-frontier.md`](docs/p3ii-delta-frontier.md): historical support-extension numerical study.
- [`docs/typeiic-incomplete-rectangles.md`](docs/typeiic-incomplete-rectangles.md): experimental Type-IIc work.
- [`docs/minimum-breakthrough.md`](docs/minimum-breakthrough.md): coupled analytic slack diagnostics.
- [`docs/primegapslib-crosspollination.md`](docs/primegapslib-crosspollination.md): reusable sparse exact-integration ideas.
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

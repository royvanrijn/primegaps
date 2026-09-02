# Analytic constraint shadow prices

This experiment ranks theorem bottlenecks by a finite counterfactual gain. For
each stable analytic constraint ID `C_i`, it recomputes the best objective over
supports whose only failed check is `C_i`:

```text
Delta lambda_i = max(score | all checks except C_i) - max(score | all checks).
```

`Delta lambda_i` is a finite one-constraint relaxation value. It is not an
infinitesimal KKT multiplier. A relaxed optimum is never a proof certificate.

## Safety boundary

`primegaps.is_certified` has no relaxation option and remains the sound theorem
oracle. `primegaps.support_constraint_failures` reports failures under stable
IDs without changing certificate semantics. Local IDs `P3.local.A` through
`P3.local.E` refer to the implemented sufficient order-statistic witness
families; failure there is conservative and is not an impossibility theorem.

The registry includes every Proposition 2 parameter/domain inequality, all four
global Proposition 3 checks, and all five local Proposition 3 partition
conditions. Structural support validity and mismatched cell metadata are not
relaxable.

## Cheap replay

An expensive scoring engine first emits `primegaps.scored-support-grid.v1`:

```json
{
  "schema": "primegaps.scored-support-grid.v1",
  "minorant": {"xi1": "0.38", "xi2": "0.4", "xi3": "0.4"},
  "score_kind": "largest-generalized-eigenvalue",
  "uncertainty_kind": "replicate-standard-error",
  "candidates": [{
    "candidate_id": "published",
    "support": {"delta": 0.028, "epsilon": 0.0075, "A": [-0.0075, 0.253], "B": [[0.15, 0.15, 0.17]]},
    "score": 1.001,
    "score_standard_error": 0.0001
  }]
}
```

`B` must contain a complete row satisfying `SupportParameters.validate`; the
short row above is schematic. Replay is deterministic and does no integration:

```bash
python scripts/rank_constraint_relaxations.py scored-grid.json \
  --output shadow-prices.json
```

Each result has status `measured` only when the grid contains at least one
candidate newly admitted by that relaxation. `unprobed` is deliberately
different from a measured zero. Rankings include measured constraints only.
The replay caches no combined index and can be rebuilt directly from the scored
grid.

## 2026-09-02 degree-21 pilot

The initial targeted screen fixed the published direct-prime minorant
`(xi1,xi2,xi3)=(0.38,0.4,0.4)` and used the prior degree-21 four-vector-bank
QMC scoring method at `k=49`. Each radius received a newly optimized
unrestricted vector bank. Boundary corrections used four scrambled-Sobol
replicates with `2^10` points per exact-large-count stratum and seeds
`64901,64902,64903,64904`.

| rank | relaxed constraint | newly admitted | baseline score | relaxed score | `Delta lambda` | combined SE* |
|---:|---|---:|---:|---:|---:|---:|
| 1 | `P3.II.delta` | 3 | 1.001811157 | 1.004588066 | 0.002776910 | 0.000028119 |
| 2 | `P3.local.D` (Type IIc) | 3 | 1.001811157 | 1.001811157 | 0 | 0 |

`*` The displayed uncertainty combines replicate standard errors as if the two
score estimates were independent. The raw runs share seeds, but the current
artifact does not retain paired-difference covariance, so this is not a rigorous
error bound.

The best `P3.II.delta` counterfactual moved the terminal endpoint from `0.253`
to `0.2537`, still satisfying every other implemented analytic check. The three
larger outer-band `B` candidates admitted by waiving Type IIc all scored below
the baseline, so this grid assigns Type IIc a measured zero. The other 18 IDs
are `unprobed`, not zero-valued.

This pilot is a numerical screen over a narrow candidate family, not a global
optimum or exact matrix certificate. Its comparative ranking remains useful:
the global Type-II delta inequality dominates the nearby Type-IIc `B` moves.
Its numerical endpoint estimate is superseded by the higher-resolution,
translated-simplex [`P3.II.delta` frontier](p3ii-delta-frontier.md), which
removes the old unbounded importance weight and measures the complete interval
to the next surviving constraint.

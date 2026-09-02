# Research status audit

Status date: 2026-09-02.

This is a durable synthesis of the local append-only research store. It does
not turn numerical observations into theorems. The latest structural replay
read all 48 records in 24 ledger shards, including negative and refuted
records, and checked all 79
referenced immutable objects by SHA-256. The current arXiv entry remains
`2608.31126v1` at this date.

## Established results

| Area | Durable conclusion | Boundary |
|---|---|---|
| 240 reproduction | A fixed 846-term `k=49`, `D=21` rational polynomial has exact `49J/I=1.00116324659...`; `49J-I>0`. | Reproduces `H_1<=240`; it is not a new bound. |
| `k=48` baseline | The analogous fixed candidate has exact `48J/I=0.99692335135...`. | This candidate is below threshold by `0.00307664865...`. |
| Distribution oracle | The published support is certified by exact Proposition 2/3 and BV checks. | Rejection by a local witness is conservative. |
| Exact engines | The frozen and accelerated I/J contractions agree on all 2,714 groups at `k=48` and `k=49`. | This validates evaluation, not optimization. |
| Moment caches | I and J moments replay exactly; higher degree requests only missing moments. | Cache context changes with support or dimension. |
| PrimeGapsLib bridge | The external 1,295-term 246 certificate replays with exact positive difference. | It is a separate, older support geometry. |
| Proof backend | Numerical search plus exact integer quadratic-form replay is implemented through the tested dimension envelope. | Matrices must be supplied independently. |

## Numerical, negative, and incomplete experiments

| Experiment | Audited conclusion | Status |
|---|---|---|
| Support geometry | The reliable two-band degree-21 screen improves the published support slightly, but the best legal `k=48` screen at the analytic ceiling is `0.9982613325` (SE `3.0018e-6`). | Numerical; no crossing. |
| Narrow `delta=0.014` support | The earlier apparent gain is driven by a rare large-integrand event; seven of eight constant-weight replicates are below the `delta=0.028` baseline. | Refuted as a positive lead. |
| `P3.II.delta` relaxation | A translated-simplex screen crosses one near `A_max=0.2536077308`, inside the counterfactual interval ending at `0.253777778055...`. | Numerical and analytically illegal without a new lemma. |
| Degree `22--27` sweep | An unrestricted proxy first crosses one at `D=26`, but it omits the published `B` cutoffs. Randomized cutoff correction and float replay fail calibration. | Inconclusive on legal support. |
| Exact `D=24` candidate | Exact I is complete and positive; J stopped at `224/7338` groups. | Incomplete; no quotient. |
| Column generation through `D=27` | The unrestricted proxy stops at `0.9980139970`. | Negative discovery result. |
| Nominal `D=22` calibration | The emitted rational polynomial has exact `48J/I=0.2686576479...`; the builder used the wrong Jacobi conversion for its intended vector. | Exact value valid, calibration refuted. |
| Support-adapted summary bases | Stable low-degree spaces score around `0.84` at best in the checked `k=49` screen. | Negative; angular information is essential. |
| Harman/Buchstab catalog | Retained-mass estimates and conditional exceptional-region choices are available. | Numerical quadrature lacks interval enclosures. |
| Surgical Type-IIc deletion | The bad demand occurs in the core negative and positive direct Type-II Buchstab terms, not exceptional pieces A/B. A rigorously enclosed mandatory positive slice costs at least `0.07158848` mass. | Deletion-only variants are rejected; a genuinely new signed identity is not ruled out. |
| Modular exact J | Checked residues agree, but one prime is no faster than rational FLINT and many primes are required. | Correct but currently negative for performance. |

## Analytic target

The active `P3.II.delta` branch is

```text
xi2/4 + 11/16 - 3*A - 2*epsilon >= delta,
```

coming from `48*omega+16*delta_star-4*gamma<-1` in the Type-IIc `wts3`
estimate. At the numerical crossing, the missing normalized slack is
`0.0013231926`, equivalently an exponent saving `0.0211710816` in the parent
estimate. Epsilon padding and dyadic bookkeeping are far too small. Recorded
examples show that neither the complete two-variable sum nor the outer `H^4`
pair count admits the required uniform saving. Any successful lemma must use
the actual incomplete ranges and preserve cancellation before the existing
absolute-value and supremum steps. No such lemma is currently proved.

## Superseded work retained for audit

- Mixed-context and incompletely bound exact checkpoints were rejected and
  archived before the sealed 240 reproduction.
- Randomized QMC/control-variate attempts at certifying the 240 threshold were
  superseded by exact rational contraction.
- The original Dirichlet-tilted `P3.II.delta` estimator was archived after its
  unbounded denominator weight was found to shift the crossing.
- The original “no score engine exists” support audit was superseded once the
  exact and numerical engines landed.
- The initial D=22 calibration vector and the narrow-delta support lead remain
  recorded as refuted/negative evidence and must not be revived as candidates.

## Research-store integrity

All ledger JSON parses, dependency targets resolve, and every referenced
immutable object exists with the recorded content hash. Seven records
have a `record_sha256` that does not equal the canonical JSON body currently in
their shard: the distribution-region record, the three PrimeGapsLib records,
the first support-geometry audit, the legal-support-geometry record, and the
Type-IIc incomplete-average record. Their claims have separate source/test or
object evidence, but their ledger identity needs an append-only correction from
the owning agent; this audit does not rewrite another agent's history.

Replay this structural audit without rerunning any mathematical experiment:

```bash
python scripts/audit_research_state.py \
  --output .research/index/research-audit.json
```

Pre-existing work without a completed ledger record is excluded from every
conclusion above. In particular, scratch sensitivity and adaptive-enrichment
runs are not promoted merely because output files exist.

## Current gates

1. Build a stable full-support matrix constructor using the exact moment caches.
2. Resolve a legal higher-degree candidate before launching more expensive
   exact contractions.
3. Prove or refute the required incomplete-range Type-IIc saving. Do not spend
   on deletion-only minorants of the current Buchstab identity: their rigorous
   optimistic no-`K` gate already fails.
4. Preserve negative results and append corrections rather than overwriting
   historical ledgers.

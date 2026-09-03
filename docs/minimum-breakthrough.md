# Minimum-breakthrough optimization

The minimum-breakthrough experiment inverts the earlier one-at-a-time shadow
price screen. For every scored support it computes simultaneous slacks for all
registered Proposition 2 and Proposition 3 conditions, then solves the finite
replay problem

\[
  \min_S \sum_i w_i s_i(S)
  \quad\text{subject to}\quad
  \widehat\lambda_k(S)-z\,\operatorname{SE}(S)\ge 1.
\]

The supplied support set and numerical variational engine determine the search
space. This is not a proof that no cheaper support exists outside that set.

## Slack convention

All default slacks and weights use raw exponent units (`w_i=1`). A global
inequality receives its exact rational shortfall from the closure of the stated
inequality. For example, if a condition requires `L >= delta`, its slack is
`max(0, delta-L)`. Equality in a strict condition consequently has zero
numerical cost but is not a theorem certificate.

For a local condition `P3.local.A` through `P3.local.E`, one scalar slack is
added to every capacity displayed for that condition. The implementation
enumerates the same order-statistic partitions as the sound distribution
oracle and records the least common capacity enlargement that makes a witness
work. The support-wide slack is the maximum over every nonvacuous Cartesian
cell pair. Since the oracle implements a sufficient witness family, local
slack is a diagnostic for that family rather than an impossibility result for
all potential Proposition 3 arguments.

`primegaps.is_certified` is unchanged. Counterfactual slacks never produce a
distribution certificate.

## Numerical gate

Discovery uses a five-dimensional projected D21 space: four stable unrestricted
directions plus the independently recorded D21 candidate direction. Every target
`k` gets its own deterministic unrestricted D21 base. D27 is deliberately not
used in this benchmark.

Legal `J` is evaluated as an importance-sampled legal-minus-unrestricted
control variate. Legal `I` is bounded above by the unrestricted-simplex `I`, so
a projected quotient above one is conservative for the same fixed direction.
The uniform translated-strata estimator is excluded: its recorded D21
calibration missed exact `I` by `0.00108345` because of rare high-leverage
polynomial tails.

Search directions are frozen before validation. Independent importance-control
replicates report their mean, standard error, minimum, and
`mean - 2*standard_error`; only the last quantity is used for a claimed
numerical crossing. These remain numerical leads, not exact rational
certificates.

## Cheap replay

The tracked replay does no integration:

```bash
PYTHONPATH=src python scripts/minimize_breakthrough.py scored-grid.json \
  --score-z 2 --output minimum-breakthrough.json
```

An optional JSON object passed with `--weights` overrides selected weights;
unspecified registered constraints retain weight one. The output includes the
complete slack vector, worst cell pair for each nonzero local slack, and every
candidate diagnostic so alternative weightings can be replayed without an
analytic or variational recomputation.

## Idealized controls

Two controls accompany the finite support search:

1. the best zero-slack support found in the Proposition-3 search family, which
   measures what perfect numerical optimization inside the current analytic
   restrictions can do;
2. full-simplex D21 generalized-eigenvalue calculations at unrestricted
   distribution exponents `theta=0.525`, `0.60`, and `66/107`.

For the second control, `U=R=theta/2`; after normalization the quotient is
linear in `theta`, so one full D21 `theta=1` matrix per `k` supplies all three
values. These are degree-21 numerical optima at the reported spectral cutoffs,
not rigorous infinite-dimensional upper bounds.

## D21 screen (2026-09-03)

The reported screen fixes `delta=0.028`, `epsilon=0.0085`, and the
prime-indicator minorant `xi=(0.38,0.4,0.4)`. It searches a fine endpoint grid
and three structured outer-band profiles. Scores use a five-dimensional D21
Ritz space trained at `2^15` points. Each selected direction is then frozen and
checked with four independent `2^17`-point importance-control replicates. The
gate is `mean - 2*SE >= 1`. Unit cost means one per raw exponent unit.

| k | cheapest validated support in the grid | mean | SE | mean-2SE | unit cost |
|---:|---|---:|---:|---:|---:|
| 47 | `A=0.2548`, current B profile | 1.0003667901 | 0.0000098434 | 1.0003471032 | 0.048633324 |
| 46 | `A=0.2560`, current B profile | 1.0005370760 | 0.0001360629 | 1.0002649501 | 0.109833324 |
| 43 | `A=0.2596`, `B2+tail` | 1.0006735611 | 0.0001168064 | 1.0004399482 | 0.319379991 |
| 42 | `A=0.2610`, `B2+tail` | 1.0006792849 | 0.0001704540 | 1.0003383769 | 0.399333324 |

`B2+tail` changes the outer-band `B_2` from `0.15` to `0.16` and the
`m>=3` cap from `0.17` to `0.18`; the lower band is unchanged. It is genuinely
useful for `k=43,42`: the unchanged-B support at the same endpoint failed the
independent two-SE gate, while moving to the next endpoint with unchanged B
cost more.

The complete nonzero slack vectors are:

| k | P3.I | P3.II.range | P3.II.delta | P3.III | local C | local D |
|---:|---:|---:|---:|---:|---:|---:|
| 47 | 0.000533334 | 0.036799990 | 0.004900000 | 0 | 0 | 0.006400000 |
| 46 | 0.005333334 | 0.079999990 | 0.008500000 | 0 | 0 | 0.016000000 |
| 43 | 0.019733334 | 0.209599990 | 0.019300000 | 0.011600000 | 0.004480000 | 0.054666667 |
| 42 | 0.025333334 | 0.259999990 | 0.023500000 | 0.016500000 | 0.014000000 | 0.060000000 |

Removing every local support/B restriction while retaining the current
asymmetric endpoints `U=0.2615`, `R=0.2445` gives full-D21 optima
`0.993484649`, `0.989235171`, `0.975676975`, and `0.970845951` for
`k=47,46,43,42`. Thus perfect local Proposition-3 optimization at the current
global endpoint does not reach any target.

The unrestricted symmetric-exponent controls are:

| k | break-even theta | theta=0.525 | theta=0.60 | theta=66/107 |
|---:|---:|---:|---:|---:|
| 47 | 0.517107551 | 1.015262684 | 1.160300211 | 1.192831992 |
| 46 | 0.519491030 | 1.010604552 | 1.154976631 | 1.187359154 |
| 43 | 0.527195580 | 0.995835359 | 1.138097553 | 1.170006830 |
| 42 | 0.529982682 | 0.990598406 | 1.132112464 | 1.163853935 |

These are finite-grid, degree-21 numerical results, not exact certificates or
global impossibility theorems. In particular, “cheapest” is conditional on the
declared weights, fixed minorant/delta/epsilon, structured support family, and
projected D21 search space.

# `P3.II.delta` counterfactual frontier

This experiment makes the global Type II delta inequality in Proposition 3 the
primary analytic relaxation target and measures the degree-21 numerical
response

\[
  A_{\max}\longmapsto \lambda_{48}=48J/I.
\]

It is a numerical screening result, not an exact rational certificate and not a
global optimization over every support family.

## Fixed family

The screen holds fixed

- `delta = 0.028`, support `epsilon = 0.0085`;
- `A = (-epsilon, 0.23, A_max)`;
- inner-band `B_1=B_2=0.18`, `B_m=0.20` for `m>=3`;
- outer-band `B_1=B_2=0.15`, `B_m=0.17` for `m>=3`;
- `(xi_1,xi_2,xi_3)=(0.38,0.4,0.4)`, hence the prime indicator;
- the symmetric degree-21 basis and the four spectral cutoffs
  `3e-12,1e-12,3e-13,1e-13`.

Every sampled endpoint chose the `1e-13` vector-bank member.

## Analytic interval

With the theorem's `epsilon=10^-10`, the second branch of `P3.II.delta` binds
at

\[
 A_{\max}=\frac{1265833333}{5000000000}=0.2531666666.
\]

After relaxing only that inequality, the next condition is `P3.II.range`,
which binds at

\[
 A_{\max}=\frac{913600001}{3600000000}
 =0.253777778055\ldots.
\]

The exact distribution oracle finds no intervening local-witness failure for
the fixed two-band family. The binary64 endpoint used for the final numerical
score lies infinitesimally to one side of the exact equality; the constraint
transition above is computed with rational arithmetic.

## Curve and crossing

The full 26-point frontier is stored in immutable object
`5034f2c9a086f7ff61c130b3c29ac6fcefa43eedb8536daa82894298ddfad14b`.
Selected points are:

| `A_max` | `lambda_48` | replicate SE |
|---:|---:|---:|
| 0.2530000000 | 0.9976043632 | 2.94e-6 |
| 0.2531666666 | 0.9982613325 | 3.00e-6 |
| 0.2533000000 | 0.9987868943 | 3.05e-6 |
| 0.2534000000 | 0.9991811009 | 3.09e-6 |
| 0.2535000000 | 0.9995753320 | 3.12e-6 |
| 0.2536000000 | 0.9999694601 | 3.15e-6 |
| 0.2536250002 | 1.0000679920 | 3.16e-6 |
| 0.2537000000 | 1.0003636695 | 3.19e-6 |
| 0.2537777781 | 1.0006702180 | 3.22e-6 |

An independent five-point, 16-replicate local screen gives slope
`d lambda_48 / d A_max = 3.94309718` and crossing

\[
 A_{\max}^{(\lambda_{48}=1)}=0.2536077308,
\]

with common-seed randomized-QMC 95% interval
`[0.2536068027, 0.2536086590]`. Thus the estimated relaxation required is

- `0.0004410642` beyond the unrelaxed `P3.II.delta` ceiling; or
- `0.0006077308` beyond the paper's `A_max=0.253` support.

The crossing consumes about `72.17%` of the one-constraint interval and leaves
about `0.0001700472` in `A_max` before `P3.II.range` binds. At the previously
sampled `A_max=0.2537`, the corrected estimate is already
`lambda_48=1.0003636695`.

The local raw output and cheap-replay summary are immutable objects
`09af66511b12249abd9e9cb4da2e6b1fc48d4c68d221e1248565ffe1947341b5`
and `1f02e5d9575221baa2978ac67bad5fc130ea6d5434b961db6e2902df8a170244`.

## Estimator correction

The earlier single-Dirichlet importance sampler had an unbounded weight and a
long right tail dominated by the `m=2` denominator correction. Its finite runs
placed the crossing too far left. The promoted sampler partitions by the exact
number `m` of coordinates above `delta`, translates the nominated coordinates
by `delta`, and samples the residual simplex uniformly. Its constant normalized
volume multiplier is

\[
  {d\choose m}\left(\frac{R-m\delta}{R}\right)^d.
\]

This removes the heavy importance weights. A closed-form inclusion--exclusion
check of the exact-`m` probability agrees with the translated-simplex estimate.

## Rebuild and replay

Install the numerical dependencies with `pip install -e '.[sparse]'`. The
expensive calculations are explicit:

```bash
PYTHONPATH=src python scripts/sweep_p3ii_delta_frontier.py \
  --shifted-strata --interior-count 16 --log2-n 9 \
  --seeds 48901,48902,48903,48904,48905,48906,48907,48908 \
  --output .research/work/p3ii-delta-frontier/delta-frontier-root/curve-shifted-raw.json

PYTHONPATH=src python scripts/sweep_p3ii_delta_frontier.py \
  --shifted-strata --only 0.253602,0.253606,0.253610,0.253614,0.253618 \
  --log2-n 10 \
  --seeds 49001,49002,49003,49004,49005,49006,49007,49008,49009,49010,49011,49012,49013,49014,49015,49016 \
  --output .research/work/p3ii-delta-frontier/delta-frontier-root/crossing-shifted-raw.json
```

Cheap replay reads only those recorded outputs:

```bash
python scripts/analyze_p3ii_delta_frontier.py \
  --curve .research/work/p3ii-delta-frontier/delta-frontier-root/curve-shifted-raw.json \
  --crossing .research/work/p3ii-delta-frontier/delta-frontier-root/crossing-shifted-raw.json \
  --output .research/work/p3ii-delta-frontier/delta-frontier-root/curve-summary.json
```

The sampling interval does not include degree truncation, vector-bank
truncation, or optimization outside the fixed support family. An exact/rational
certificate remains necessary before claiming a bound below 240.


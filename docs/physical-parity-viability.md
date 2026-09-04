# Physical parity through quarter-rough almost-primes

## Result

The idealized experiment is viable.  With a perfect
`x^(1/4+epsilon)`-rough detector and perfect asymptotics for the signed
degree-three factorial-moment combination, the frozen PrimeGaps186 physical
model gives the following production-mesh predictions:

| dimension | ideal full-face score | margin above 1 | common relative error across the three unsigned terms |
|---:|---:|---:|---:|
| 39 | 1.0192774814 | 0.0192774814 | 0.264108% |
| 38 | 1.0127855640 | 0.0127855640 | 0.176290% |

Thus the answer to the requested stop/go test is **go**.  Perfect parity does
not give `1.04` in this frozen 77-dimensional basis, but it clears both targets.
The admissible absolute error in the final normalized parity combination is
about `1.93% I` for `k=39` and `1.28% I` for `k=38`.  Cancellation makes the
corresponding termwise relative requirements much smaller.

This is an exploratory binary64 mesh extrapolation, not a certificate and not
a prime-distribution theorem.  Its purpose is to identify the analytic theorem
worth attacking.

## The arithmetic detector

Fix `1/4 < beta < 1/3` and put `z=x^beta`.  If `m` is in a fixed dyadic
interval of size `x` and `P^-(m)>z`, then, for all sufficiently large `x`,
`N=Omega(m)` is at most three.  On `N=0,1,2,3`,

```text
1_{N=1} = N - 2 C(N,2) + 3 C(N,3).
```

Equivalently, writing `lambda_L(m)=(-1)^Omega(m)`,

```text
1_{N=1} = (1-lambda_L(m))/2 - C(N,3).
```

The second form isolates the genuinely missing information: a weighted
Liouville correlation.  The rough count and the third factorial moment are
nonnegative sieve quantities.

At leading `x/log x` scale, repeated prime factors are negligible.  The
semiprime and triprime constants are

```text
c2(beta) = log((1-beta)/beta),
c3(beta) = (1/6) integral integral du dv / (u v (1-u-v)),
```

where `u,v >= beta` and `u+v <= 1-beta`.  Hence

```text
E N          = 1 + 2 c2 + 3 c3,
E C(N,2)     = c2 + 3 c3,
E C(N,3)     = c3.
```

For `beta=0.250001`, the constants are

```text
c2             = 1.0986069553418876
c3             = 0.14721698020054072
E N            = 3.638864851285397
E C(N,2)       = 1.5402578959435096
E C(N,3)       = 0.14721698020054072
signed identity = 1.0000000000000000
```

The reusable implementation is in `primegaps.parity`.  It is intentionally
separate from the physical quadrature: the repository's existing
`FactorialMomentTable` integrates symmetric trial polynomials and is not a
table of the arithmetic quantities `C(Omega(m),r)`.

## How roughness enters the physical state

The experiment adds the factor-count state `N in {1,2,3}` to each erased target
face.  Its three components use the same physical full-face quadratic form,
multiplied by the three arithmetic constants above.  This is an idealization;
the current source theorem supplies no signed asymptotic with which to justify
that tensor-product rule.

There is nevertheless a useful structural compatibility.  The roughness
threshold `beta=0.250001` is `beta/rho*=0.952389` in physical-fragment units,
whereas the largest target fragment admitted by the full-face shell is
`190370/262499=0.725222`.  Thus roughness forces the target-coordinate sieve
divisors to be trivial throughout this model.  The other coordinates retain
their physical fragment state.  Under the ideal asymptotic this opens

```text
Jfull = J0 + Jplus + Jtail,
```

rather than the currently certified signed hybrid face with source losses.
The crude pair-support envelope is `2 rho* S = 0.5485994`; proving parity
cancellation under the resulting structured, high-modulus congruence weight is
the central obstruction.

## Numerical calculation

The calculation ports the pinned PrimeGaps186 positive midpoint/Dickman cap
model into explicit 77 by 77 matrices.  It freezes `rho*`, radii, caps, profile,
angular signatures and radial degrees, varies the dimension, and optimizes
`rho* Jfull/I` using a spectrally truncated generalized eigenproblem.  Direct
convolution was used through the finest discovery mesh because FFT roundoff is
unreliable in the long exact-zero prefixes of high convolution powers.

| mesh intervals | optimized k=39 | optimized k=38 |
|---:|---:|---:|
| 1,024 | 0.9994308621 | 0.9935784374 |
| 2,048 | 1.0095225559 | 1.0033423633 |
| 4,096 | 1.0145211148 | 1.0081804026 |
| 8,192 | 1.0170096787 | 1.0105895791 |
| 98,304, extrapolated | 1.0192774814 | 1.0127855640 |

The production prediction fits a quadratic in `1/N` to the 2,048--8,192
evaluations of the finest-mesh vector.  Fitting separately optimized values and
using all meshes changes the production prediction by at most `7.1e-6` in the
reported cases.  As a calibration, the identical extrapolation applied to the
published `k=40` hybrid vector predicts `1.0002004816`, only `-5.60e-6` below
the published rigorous lower-endpoint score `1.0002060794`.  This calibration
is evidence for the scale of the discretization error, not an interval bound.

The production factorial contributions are:

| dimension | `+N` | `-2 C(N,2)` | `+3 C(N,3)` | signed sum | absolute sum |
|---:|---:|---:|---:|---:|---:|
| 39 | 3.7090130008 | -3.1399003777 | 0.4501648584 | 1.0192774814 | 7.2990782369 |
| 38 | 3.6853897907 | -3.1199019237 | 0.4472976970 | 1.0127855640 | 7.2525894114 |

If the final signed asymptotic has normalized error at most `eta I`, then one
needs `eta < score-1`.  If instead every unsigned term is known only to a
common relative error `epsilon`, the triangle inequality requires
`epsilon < (score-1)/(absolute sum)`, producing the last column of the first
table.

## The exact missing correlation

Let the physical Selberg weight be written abstractly as

```text
w_F(n) = (sum_d Lambda_d product_j 1_{d_j | n+h_j})^2,
```

where `d` is a vector of squarefree divisors.  For target coordinate `i`, put

```text
A_{i,r} = sum_n w_F(n) 1_{P^-(n+h_i)>z} C(Omega(n+h_i),r),
L_i     = sum_n w_F(n) 1_{P^-(n+h_i)>z} lambda_L(n+h_i).
```

The prime-weighted sum is exactly `A_{i,1}-2A_{i,2}+3A_{i,3}`, or equivalently
`(A_{i,0}-L_i)/2-A_{i,3}`, once the rough support has `Omega<=3`.

Write

```text
r_beta   = 1 + c2 + c3,
ell_beta = -1 + c2 - c3 = -0.04861002485865315.
```

Assuming the positive rough and third-moment asymptotics, the direct theorem
needed from parity is the full-face weighted estimate

```text
(rho*/I_F) sum_i (L_i - ell_beta Jfull_i) = o(1).                 (P)
```

For the displayed finite viability margins it is enough to replace `o(1)` in
absolute value by `0.0192774814` for `k=39`, or by `0.0127855640` for `k=38`,
in the same normalization.  A genuine asymptotic proves both.

Expanding the square shows exactly what `(P)` asks us to cancel.  Since
roughness forces `d_i=e_i=1`,

```text
L_i = sum_{d,e: d_i=e_i=1} Lambda_d Lambda_e
        sum_m lambda_L(m) 1_{P^-(m)>z}
          product_{j != i} 1_{m = h_i-h_j (mod [d_j,e_j])}.       (E)
```

Thus the missing object is not an unweighted Chowla sum for a fixed pair.  It
is Liouville cancellation after summing a signed 77-dimensional family of
structured lcm congruences reaching polynomial-size moduli.

A Friedlander--Iwaniec asymptotic-sieve route would prove a weighted bilinear
axiom for the nonnegative sequence

```text
a_i(m) = 1_{P^-(m)>z} w_F(m-h_i).
```

With `gamma(v,C)=sum_{c|v,c<=C} mu(c)`, the required family has the form

```text
sum_u | sum_{V<v<=2V, uv in the target interval}
          gamma(v,C) mu(uv) a_i(uv) | <= epsilon_x sum_m a_i(m),  (B_F)
```

uniformly in the dyadic ranges used by the asymptotic-sieve decomposition,
with `epsilon_x -> 0`.  Expanding `a_i` in `(B_F)` gives exactly the lcm
congruences in `(E)`.  This coefficient-level estimate is the right target:
bounding every congruence class separately would discard cancellation among
the physical coefficients and demand far more than the viability calculation
uses.  To turn a non-asymptotic bound for each dyadic block into the numeric
budgets above, the asymptotic-sieve/Buchstab decomposition must also retain its
explicit block constants; the present experiment does not pretend that one
block may consume the whole aggregate budget.

Friedlander and Iwaniec's axiom `(B)` is precisely a sum of absolute values of
inner Mobius bilinear forms, with a truncated-divisor coefficient `gamma`; see
pp. 1043--1046 and Theorem 1 of [*Asymptotic sieve for
primes*](https://arxiv.org/abs/math/9811186).  Their paper explains that the
Mobius sign changes supply the information excluded by the classical parity
barrier.  Establishing `(B_F)` for the high-modulus, rough-conditioned physical
weight is the new theorem, not a consequence of the existing PrimeGaps186
source bounds.

Tao and Teräväinen's 2026 work supports the general moment strategy: their
Theorem 1.1 combines a growing-dimensional Maynard-type sieve with second and
fourth moment calculations, while Theorem 3.1 supplies log-power two-point
multiplicative cancellation.  But their correlation modulus and shifts are
only polylogarithmic and the conclusion excludes a logarithmically sparse set
of scales; see [*Quantitative correlations and some problems on prime factors
of consecutive integers*](https://arxiv.org/abs/2512.01739).  It therefore
does not imply `(P)` or `(B_F)` for the polynomial-modulus physical weight.

## Multi-candidate moment variant

One can instead study the number of prime-parity candidates among all rough
tuple coordinates.  A second/fourth-moment argument would require weighted
two- and four-point correlations of their Liouville signs under the same
physical weight.  This could avoid an asymptotic for any predetermined pair,
but the four-point input is stronger than current two-point technology and the
same polynomial-modulus conditioning remains.  It is a legitimate secondary
route, not a shortcut around `(P)`.

## Reproduction

The cheap replay recomputes all rough constants, signed contributions and
error budgets from the recorded mesh result:

```bash
python scripts/check_physical_parity_viability.py
```

The expensive discovery calculation is explicit and separate:

```bash
pip install -e '.[sparse]'
python experiments/physical_parity_viability.py \
  --dimensions 40 39 38 --intervals 8192 \
  --output result.json --matrix-output matrices.npz
```

Run the expensive command independently at 1,024, 2,048, 4,096 and 8,192
intervals, then pass the result and matrix files in increasing-mesh order to
`experiments/analyze_physical_parity_meshes.py`.  The pinned numerical inputs
are in `reproduction/186/physical-parity-input.json`; bulk matrices remain
content-addressed research objects rather than tracked files.

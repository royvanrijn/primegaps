# Physical parity through beta-rough almost-primes

## Corrected status

The numerical experiment remains a useful **ideal-geometry viability bound**:
with a perfect prime-sector projector on the full PrimeGaps186 physical face,
the frozen 77-dimensional model predicts

| dimension | ideal full-face score | margin above 1 |
|---:|---:|---:|
| 39 | 1.0192774814 | 0.0192774814 |
| 38 | 1.0127855640 | 0.0127855640 |

However, the Friedlander--Iwaniec bilinear axiom proposed in the previous
version of this note is structurally false for the quarter-rough sequence.  The
correct classification is therefore:

```text
ideal physical geometry:       GO
factorial prime identity:      exact
proposed F-I axiom B_F:        NO-GO (pointwise sign obstruction)
actual parity architecture:    open
```

The experiment shows that the physical geometry has enough latent score for
`k=39` and `k=38`.  It does **not** identify a valid route to the required
parity theorem.

This is an exploratory binary64 mesh extrapolation, not a certificate and not
a prime-distribution theorem.  Its purpose is to identify the analytic theorem
worth attacking.  The threshold sweep below makes `beta=1/3+epsilon` the first
target for that theorem; the quarter-rough calculation remains the baseline
that maximizes the rough carrier and the cancellation condition number.

## The arithmetic detector

Fix `1/4 < beta < 1/3`, put `z = X^beta`, and restrict a target integer
`m ~ X` by `P^-(m) > z`.  For sufficiently large `X`,
`N = Omega(m)` then lies in `{1,2,3}`.  On `N = 0,1,2,3`,

```text
1_{N=1} = N - 2 C(N,2) + 3 C(N,3).
```

Equivalently,

```text
1_{N=1} = (1-lambda_L(m))/2 - C(N,3).
```

For `beta = 0.250001`, the leading `X/log X` constants are

```text
prime                         1.0000000000000000
semiprime                     1.0986069553418876
triprime                      0.14721698020054072
E Omega                       3.638864851285397
E C(Omega,2)                  1.5402578959435096
E C(Omega,3)                  0.14721698020054072
signed factorial identity     1.0000000000000000
absolute factorial condition  7.161031583774038
```

The identity and these constants remain correct.  The failure is in the
suggested analytic route, not in the idealized finite-dimensional calculation.

## Numerical experiment retained

The calculation ports the pinned PrimeGaps186 midpoint/Dickman cap model into
explicit `77 x 77` matrices, freezes the profile, masks, signatures, radial
degrees and `rho*`, and optimizes

```text
rho* Jfull / I,
Jfull = J0 + Jplus + Jtail.
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

## Sweeping the roughness threshold

The quarter-rough endpoint is the worst-conditioned point of the usable
degree-three interval, not a distinguished optimum.  The ideal full-face score
is independent of `beta`, but the unsigned mass that must cancel is

```text
kappa_3(beta) = 1 + 4 c2(beta) + 12 c3(beta),   1/4 < beta <= 1/3.
```

For `beta>1/3`, a rough integer has at most two prime factors and the exact
identity drops to

```text
1_{N=1} = N - 2 C(N,2) = (1-lambda_L)/2,
kappa_2(beta) = 1 + 4 c2(beta).
```

The transition is continuous in the leading constants: `c3` tends to zero as
`(81/4)(1/3-beta)^2`.  It is structurally discontinuous in the proof, because
the nonnegative triprime correction disappears strictly above `1/3`.

| `beta` | detector degree | rough carrier `1+c2+c3` | `c3` | condition `kappa` | fraction of quarter-rough condition | k=39 common relative budget | k=38 common relative budget |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.250001 | 3 | 2.245824 | 0.147217 | 7.161032 | 1.000 | 0.2641% | 0.1763% |
| 0.275 | 3 | 2.039869 | 0.070468 | 5.723220 | 0.799 | 0.3305% | 0.2206% |
| 0.300 | 3 | 1.869964 | 0.022666 | 4.661185 | 0.651 | 0.4058% | 0.2708% |
| 0.325 | 3 | 1.732294 | 0.001407 | 3.940433 | 0.550 | 0.4800% | 0.3204% |
| `1/3 - 10^-6` | 3 | 1.693152 | `2.03e-11` | 3.772607 | 0.527 | 0.5013% | 0.3346% |
| `1/3 + 10^-6` | 2 | 1.693143 | 0 | 3.772571 | 0.527 | 0.5013% | 0.3346% |
| 0.350 | 2 | 1.619039 | 0 | 3.476157 | 0.485 | 0.5441% | 0.3632% |
| 0.400 | 2 | 1.405465 | 0 | 2.621860 | 0.366 | 0.7214% | 0.4815% |
| 0.450 | 2 | 1.200671 | 0 | 1.802683 | 0.252 | 1.0492% | 0.7003% |

The factorization geometry shrinks at the same time.  Semiprimes have their
smaller prime exponent in `[beta,1/2]`, a strip of width `1/2-beta`.  Below
`1/3`, a triprime split into one prime and a complementary semiprime has

```text
u in [beta, 1-2 beta],       1-u in [2 beta, 1-beta],
```

and the three-factor simplex has slack `1-3 beta`.  Above `1/3` that whole
prime-times-semiprime block is empty.  The remaining weighted Liouville
bilinear estimate is still essential; degree two simplifies its decomposition
but does not prove it.

The price is a strictly deeper roughness condition `P^-(m)>x^beta`.  In the
frozen physical geometry the largest target fragment exponent is only
`0.19036993`, so every swept value still forces target-coordinate Selberg
divisors to be trivial.  The crude physical pair-support exponent is
`0.5485994`, whose ratio to `beta` falls from `2.1944` at `0.250001`, to
`1.6458` at `1/3`, and to `1.3715` at `0.4`.  This ratio is only a scale
diagnostic, not a proved rough-sieve level: any proposed theorem must separately
show that it can impose the deeper sifting condition while retaining the
physical lcm weight.

Consequently the first analytic target should be `beta=1/3+epsilon`, not the
quarter-rough endpoint.  It captures almost all of the cheap structural gain:
the condition number falls by `47.3%`, both relative error budgets increase by
about `90%`, and the triprime correction vanishes, for an increase of only
about `1/12` in the roughness exponent.  Values such as `0.4` are useful second
targets if the rough-sieve theorem has enough depth.

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

Applying the same extrapolation to the published `k=40` hybrid vector misses
its rigorous lower-endpoint score by only `-5.60e-6`.  This is useful
calibration, not an interval proof.

### Minimal generated-modulus reach at k=39

The full face need not be treated as an indivisible endpoint.  For a retained
38-coordinate face state `Y` and its 39-coordinate outer extension
`Y +_i X`, define the artificially truncated erased-coordinate operator

```text
E_{i,theta} F(Y) = integral F(Y +_i X)
    1_{rho* (|Y| + |Y +_i X|) <= theta} dnu(X),
J_theta = sum_i E_{i,theta}^* E_{i,theta},
Lambda_39(theta) = rho* lambda_max(J_theta, I).
```

This is the generated outer-plus-inner modulus exponent used by the physical
source transfer.  On the mesh, the cut is made conservatively with both upper
cell endpoints.  It is imposed inside the erased-coordinate integral, before
squaring; replacing it by the stronger symmetric radius cut
`2 rho* |Y +_i X| <= theta` would throw away useful asymmetric edges.

Four direct-convolution meshes and the same production extrapolation give:

| theta | production `Lambda_39(theta)` |
|---:|---:|
| `0.500000` | `0.9787438843` |
| `0.510000` | `0.9913935182` |
| `0.515000` | `0.9971224045` |
| `0.517000` | `0.9993292535` |
| `0.517500` | `0.9999012600` |
| `0.517625` | `0.9999535605` |
| `0.517750` | `1.0000967522` |
| `0.518000` | `1.0003821396` |
| `0.520000` | `1.0025203787` |
| `0.525000` | `1.0073642488` |
| `0.530000` | `1.0115525121` |
| `0.5485994` | `1.0192774815` |

Thus the scanned crossing is in `[0.517625, 0.517750]`, with linear estimate
`theta_c = 0.51766554`.  Separately optimizing at every mesh gives
`0.51763978`; fitting all four fixed-vector meshes gives a value near `0.51760`.
The defensible conclusion is therefore `theta_c about 0.5176`, not a claim at
the fifth decimal place.  In particular, the model does not cross below the
Bombieri--Vinogradov boundary, but it also does not require the full
`0.5485994` envelope.  Only about 36% of the full exponent interval above
one-half is needed by this blunt low-pass truncation.

This strengthens the case for a hybrid theorem: parity cancellation restricted
to selected high-value source states could plausibly preserve a crossing while
avoiding the last `0.031` of modulus exponent.  The calculation remains an
exploratory binary64 operator experiment.  It neither proves monotonicity
between sampled theta values nor supplies the required correlation estimate.

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

## Why the proposed `B_F` is impossible

Let

```text
a_i(m) = 1_{P^-(m) > z} w_F(m-h_i) >= 0.
```

Friedlander and Iwaniec's axiom `(B)` requires, for every admissible inner
length `V` and **every** `1 <= C <= XD^{-1}`,

```text
sum_u | sum_{V < v <= 2V, uv <= X}
          gamma(v,C) mu(uv) a_i(uv) |
    <= A_i(X) (log X)^(-2^22),

gamma(v,C) = sum_{d|v, d<=C} mu(d).
```

Their technically reduced hypothesis `(B')` is weaker but still requires a
`(log X)^(-3)` saving.  The exponent in the displayed original axiom is
`-2^22`, not the decimal integer `-222`.

Now take the mandatory case `C=1`.  Then `gamma(v,1)=1`.  In every active block
of the proposed application,

```text
2V < sqrt(X) < X^(2 beta) = z^2.
```

If `a_i(uv) != 0`, then `uv` is `z`-rough, hence its divisor `v` is also
`z`-rough.  Since `1 < v < z^2`, `v` cannot contain two prime factors and is
therefore prime.  For fixed `u`, every nonzero Mobius term consequently has the
same sign:

```text
mu(uv) = -mu(u)       when u is squarefree and v does not divide u,
mu(uv) = 0            otherwise.
```

Thus

```text
| sum_v mu(uv) a_i(uv) |
  = sum_v mu(uv)^2 a_i(uv).
```

The inner absolute value is an exact positive mass; there is no Mobius
cancellation to estimate.  In the modeled decomposition a nonempty dyadic
prime block has size of order `A_i(X)/log X`, far larger than either
`A_i(X)(log X)^(-3)` or the original
`A_i(X)(log X)^(-2^22)` requirement.

Two apparent escapes do not help:

- blocks with `2V <= z` are empty, but the active blocks are not;
- choosing `C >= 2V` makes `gamma(v,C)=0`, but the axiom must also hold for
  `C=1`.

Therefore every nonempty active, recombined `B_F` block fails for structural
reasons.  This is not a missing estimate and cannot be repaired by stronger
trace-function or dispersion technology while retaining the same outer
absolute-value architecture.

## Why a large modulus exponent does not transfer

Expanding the physical Selberg square gives a residue-coloured coefficient
family of the schematic form

```text
c_i(q,a) = sum_{d,e producing (q,a)} Lambda_d Lambda_e,
```

with

```text
q <= X^0.5485994,
p | q  =>  p <= X^0.1903699274 < X^beta.
```

This is only a support envelope.  Maynard's and Lichtman's large-modulus
results require well-factorable or triply well-factorable coefficients in the
precise convolutional sense, for every requested factorization of their level.
The physical `c_i(q,a)` are residue-coloured CRT aggregates and have not yet
been defined with their discrete normalization, let alone proved to satisfy
those hypotheses.  Numerical compatibility of the exponent ranges is not a
coefficient transfer theorem.

Before invoking any averaged distribution theorem, the project needs the
actual dyadic norm and factorization data, beginning with

```text
sum_{Q < q <= 2Q} sum_{a mod q} |c_i(q,a)|^2.
```

The indexed `(U,V,q)` family, exact physical modulus envelopes, missing
coefficient norm, and block classifications are in the
[dyadic audit](bf-dyadic-audit.md).  The direct signed statement `(P)` remains
a meaningful target, but it needs a different decomposition which preserves
cancellation between factor-count sectors before taking absolute values.

The audit also records that the already rough sequence violates the paper's
dimension-one density axiom at every prime below `z`, and that no analogue of
its required remainder axiom with `D>X^(2/3)` has been supplied.

Tao and Teräväinen's 2026 work supports the general moment strategy: their
Theorem 1.1 combines a growing-dimensional Maynard-type sieve with second and
fourth moment calculations, while Theorem 3.1 supplies log-power two-point
multiplicative cancellation.  But their correlation modulus and shifts are
only polylogarithmic and the conclusion excludes a logarithmically sparse set
of scales; see [*Quantitative correlations and some problems on prime factors
of consecutive integers*](https://arxiv.org/abs/2512.01739).  It therefore
does not imply `(P)` for the polynomial-modulus physical weight.

## Viable replacement 1: a global factor-count projector

The prime, semiprime and triprime sectors must be combined **before** any outer
absolute value.  A well-conditioned way to do this is the generating function

```text
Z_i(omega) = sum_m a_i(m) omega^Omega(m),
omega in {1, i, -1, -i}.
```

Because `Omega(m) <= 3`, discrete Fourier inversion gives the exact pointwise
identity

```text
1_{Omega(m)=1}
  = (1/4) sum_{r=0}^3 i^(-r) (i^r)^Omega(m).
```

At `beta=0.250001`, the four leading sector transforms are

```text
Z(1)    =  2.245823935542428
Z(i)    = -1.0986069553418876 + 0.8527830197994593 i
Z(-1)   = -0.04861002485865315
Z(-i)   = conjugate(Z(i)).
```

The corresponding absolute projection condition is about

```text
( |Z(1)| + |Z(i)| + |Z(-1)| + |Z(-i)| ) / 4
  = 1.2689817929354914,
```

rather than `7.1610315838` for the factorial basis.  This does not prove any
correlation estimate, but it substantially improves the numerical conditioning
of a future sector-coupled theorem.

The analytic target would be a **single global residue-coloured estimate** for
these four completely multiplicative phases, or directly for their Fourier
projector, with the absolute value taken only after summing the sector,
modulus, residue and dyadic contributions.  Applying a Friedlander--Iwaniec
absolute value separately for each fixed outer variable would recreate the
same obstruction.

This architecture is related to Selberg--Delange generating functions and to a
Fourier decomposition of `Omega mod 4`; it is not an application of the
classical asymptotic-sieve axiom `(B)`.

## Viable replacement 2: a sector-coupled Buchstab operator

Keep the three arithmetic sectors as a vector

```text
P = (prime, semiprime, triprime)
```

through the Buchstab/Heath--Brown decomposition.  Extracting one rough prime
factor shifts the factor-count state.  One can encode this by a small transfer
matrix, or equivalently by the generating variable `omega` above.

For each dyadic geometry block, estimate only the projected combination

```text
E_prime - 2 E_semiprime + 3 E_triprime
```

(or its root-of-unity equivalent), and postpone absolute values until after the
three blocks share common outer variables and CRT-coloured coefficients.

This is the architecture explicitly required by the sign obstruction: the
prime-only contribution in one Friedlander--Iwaniec inner block must be allowed
to cancel against semiprime and triprime contributions elsewhere in the same
joint operator.  No such aligned decomposition has yet been derived.

## Viable replacement 3: global residue-coloured dispersion

Define the exact coefficient compiler

```text
(d,e) -> (q,a) -> c_i(q,a)
```

including:

- the discrete normalization of `Lambda_d`;
- CRT consistency and collisions;
- all signs and multiplicities;
- the physical masks and target-coordinate restriction;
- dyadic `q` and residue grouping.

Then write the projected parity error as

```text
sum_{q,a} c_i(q,a) E_projected(X;q,a).
```

A plausible proof would use an `l2`/spectral or dispersion estimate for this
**combined** coefficient family.  It must exploit cancellation among the
physical coefficients and among factor-count sectors.  Bounding each residue
class separately, or replacing `c_i(q,a)` by its support envelope, loses the
only cancellation the viability experiment is trying to use.

## Optional structural salvage of an F-I-style axiom

An F-I inner sum can only have Mobius sign changes if the active interval
contains rough composites.  A necessary condition is roughly

```text
V_min > z^2 = X^(2 beta).
```

Lowering `beta` enough could remove the prime-only obstruction, but then:

- `Omega(m)` may be as large as `floor(1/beta)`;
- the prime projector needs a correspondingly deeper factor-count basis;
- target-coordinate sieve divisors are no longer automatically trivial;
- the physical state and coefficient compiler become much larger.

This is worth a cheap exponent audit, but it is not presently the preferred
route.

## Immediate experiments

### 1. Discrete coefficient compiler

On a reduced prime universe and coarse physical mesh, construct the actual
`c_i(q,a)` and report by dyadic `Q`:

```text
l1 norm,
l2 norm,
number of occupied residues,
CRT collision/cancellation ratio,
best rank-r approximation,
well-factorable and triply-factorable residuals.
```

This is the prerequisite for any claim involving Maynard, Lichtman, Pascadi or
a spectral large sieve.

### 2. Root-of-unity projector sweep

For `beta` between `0.250001` and just below `1/3`, compute the transformed
sector constants and the resulting common-error condition number.  Reuse the
existing full-face score; no new physical matrix calculation is needed for
this first screen.

### 3. Sector-coupled finite model

At moderate finite `X`, expand the coarse physical coefficients and compare

```text
sum_blocks |prime block| + 2|semiprime block| + 3|triprime block|
```

with

```text
|sum_blocks (prime block - 2 semiprime block + 3 triprime block)|.
```

Do this before and after CRT aggregation.  A persistent large cancellation
ratio is evidence for a global operator theorem; a ratio near one is a quick
no-go.

### 4. F-I structural threshold audit

Compute the exponent of the smallest mandatory inner block and determine the
largest `beta` for which `V_min > X^(2 beta)`.  Record the implied maximum
factor count and projector dimension.  This decides whether lowering roughness
could ever make a local Mobius-bilinear architecture plausible.

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

The independent beta sweep reuses those frozen scores and does not rerun the
physical matrices:

```bash
python scripts/sweep_physical_parity_beta.py \
  --output experiments/physical_parity_beta_sweep.json
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

The modulus-reach scan and its cheap replay are separate:

```bash
python experiments/physical_parity_modulus_reach.py \
  --intervals 8192 \
  --theta 1/2 51/100 103/200 517/1000 207/400 4141/8000 \
          2071/4000 259/500 13/25 21/40 53/100 2742997/5000000 \
  --output reach-n8192.json --matrix-output reach-n8192.npz
python scripts/check_physical_parity_modulus_reach.py
```

Repeat the first command at all four meshes and use
`experiments/analyze_physical_parity_modulus_reach.py` to rebuild the tracked
summary.  This expensive computation is never run by tests or by the replay.

The recorded scores should be read as a bound on what a perfect sector-coupled
prime detector could extract from the frozen physical geometry, not as evidence
for the rejected `B_F` axiom.

## Primary references

- J. Friedlander and H. Iwaniec, *Asymptotic sieve for primes*, especially
  axioms `(B)`, `(B1)--(B3)` and the reduced hypothesis `(B')`.
- J. Maynard, *Primes in arithmetic progressions to large moduli II:
  well-factorable estimates*, Definitions 1--2.
- J. D. Lichtman, *Primes in arithmetic progressions to large moduli, and
  Goldbach beyond the square-root barrier*, Definitions 1.3--1.4.

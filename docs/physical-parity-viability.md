# Physical parity through quarter-rough almost-primes

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

The directly computed and extrapolated scores are

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

## Reproduction

The ideal numerical calculation and cheap arithmetic replay remain available:

```bash
pip install -e '.[sparse]'
python experiments/physical_parity_viability.py \
  --dimensions 40 39 38 --intervals 8192 \
  --output result.json --matrix-output matrices.npz
python scripts/check_physical_parity_viability.py
```

The recorded scores should now be read as a bound on what a perfect
sector-coupled prime detector could extract from the frozen physical geometry,
not as evidence for the rejected `B_F` axiom.

## Primary references

- J. Friedlander and H. Iwaniec, *Asymptotic sieve for primes*, especially
  axioms `(B)`, `(B1)--(B3)` and the reduced hypothesis `(B')`.
- J. Maynard, *Primes in arithmetic progressions to large moduli II:
  well-factorable estimates*, Definitions 1--2.
- J. D. Lichtman, *Primes in arithmetic progressions to large moduli, and
  Goldbach beyond the square-root barrier*, Definitions 1.3--1.4.

# Findings / reconstruction status

Date: 2026-09-01

Source: Julia Stadlmann, *Bounded gaps between primes*, arXiv:2608.31126v1.

## What is solidly established from the paper

- The paper proves `H_1 <= 240`, corresponding to the shortest admissible
  49-tuple. Polymath8b's `246` used a 50-tuple.
- Stadlmann obtains the new bound with symmetric polynomials satisfying
  `2a+b <= 21`; Polymath used degree `<= 27` for 246. The author explicitly
  states that time/memory prevented the degree-27 calculation and expects more
  resources to improve the bound further.
- The published final parameters are:
  - `epsilon = 0.0075`
  - `delta = 0.028`
  - `A = (-epsilon, 0.253)`
  - `B[1,1] = B[1,2] = 0.15`
  - `B[1,m] = 0.17` for `m >= 3`
  - `xi_1 = 0.38`, `xi_2 = xi_3 = 0.4`
- For those parameters Proposition 2 allows using the prime indicator directly
  (`c1=c2=0`), i.e. the general Harman-minorant machinery is not needed for 240.
- The final numerical certificate is a vector `c` with generalized Rayleigh
  quotient `c M2 c^T / (c M1 c^T) > 1` for `k=49`.

## Implemented so far

- Exact published support `T_k(delta,A,B,epsilon)` and validation.
- Symmetric basis enumeration. At `k=49` the basis grows from **846 functions at
  D=21** to **2526 at D=27**, so a dense matrix grows from about 0.72M to 6.38M
  entries (~8.8x) before considering exact-integral construction cost.
- A support-geometry Monte-Carlo diagnostic. Uniform volume is a bad proxy for
  sieve gain: at `k≈49` most uniform-simplex points do not probe the interesting
  `B` boundary. The actual objective must remain the `I/J/K` integral quotient.
- The Section 5.2.1 **exact k=1 recurrence base cases** are now implemented as
  rational coefficient vectors for
  `int_0^delta t^a(1-t)^b dt` (`C`) and
  `int_delta^1 t^a(1-t)^b dt` (`D`). Tests verify exactly that C+D equals the
  integer beta integral for multiple exponents, including the published
  `delta=7/250`.
- The omitted multidimensional formulas have now been reconstructed.
  `small_cube_coefficients` computes `C_{m,i}` by exact inclusion--exclusion and
  Dirichlet integration; `large_simplex_coefficients` computes `D_{m,i}` after
  shifting the large simplex. They reduce to every tested k=1 base case and
  agree exactly across reciprocal-delta chamber boundaries.
- A reference exact `I/J/K` assembler accepts sparse polynomial basis functions
  and rational support parameters and returns rational matrix entries. Tiny
  matrices agree with independent nested Gauss--Legendre quadrature, the
  full-simplex Dirichlet identity, and a case with an active `B` boundary.
- Axiom Math's PrimeGapsLib certificate layer has been independently audited at
  commit `1faa7b14e82ddebc2772dfb9153922f01b106477` and its reusable machinery
  ported without a Lean dependency. Canonical exponent signatures, erasure
  closure, adaptive overlap-profile moments, the nilpotent `(I+N)^k` moment
  ladder, exact monomial-symmetric simplex/marginal pairings, and degree-grouped
  sparse contractions are now implemented in `primegaps.symmetric`.
- Replaying PrimeGapsLib's 1,295-term `k50e25d25n1295.json` source certificate
  locally uses 272 erase-closed signatures, 138 mass groups, and 172 marginal
  groups. Exact integer arithmetic gives `50*J-4*I > 0` and quotient
  approximately `4.000000917784551`. This reproduces the external 246
  variational certificate, not Stadlmann's k=49 support matrix.
- The complete generalized-eigenvalue *end of the pipeline* is implemented:
  diagonally equilibrated symmetric generalized eigenvalue search in floating
  point, conditioning/residual diagnostics, candidate-vector rationalization,
  and a portable exact certificate. Dense, matrix-free Lanczos, simultaneous
  block, and optional sparse paths are available. The standalone verifier checks
  semantic matrix hashes and replays both quadratic forms using Python integers.
  This mirrors the
  paper's strategy: floating point locates `c`; exact arithmetic proves it.
- The `k=49`, `D=21` bottleneck has now been reproduced end-to-end. A fixed
  846-term decimal-rational candidate was contracted by exact rational,
  symmetry-compressed `I` and `J` evaluators over the published support. Cheap
  replay checks all source/candidate/support/signature hashes and gives
  `49J/I = 1.0011632465949216560417861678682244509240906847502...`, hence
  `49J-I > 0` exactly.
- An accelerated exact backend now preserves that frozen evaluator as an oracle
  while changing the implementation: it collapses zero-coordinate density
  transitions combinatorially, assembles J feature-pair-first, uses carry-free
  FLINT polynomial encoding over `QQ` or `GF(p)`, supports bounded CRT rational
  reconstruction, and persists candidate-independent I/J moments. Exact I
  assembly uses about one quarter of the original summed worker-task time at
  both `k=48` and `k=49`; every one of the 2,714 rows matches the frozen result.
  At `k=48`, all
  2,714 accelerated J rows also match; a lightly contended prototype used
  3.84x less summed worker-task time and 1.84x less wall time despite half as
  many workers. The manifest-bound `k=49` J run also matches all 2,714 rows,
  uses 3.47x less summed worker-task time, and reproduces the exact
  positive `49J-I`; the accelerated backend has therefore passed its oracle
  gate. The same manifest-bound code at `k=48` reproduces the exact published
  deficit. The modular backend is exact in checked
  residues but showed no single-prime speed advantage, so the FLINT rational
  path remains the default until geometry can be shared across batched primes.
- Proposition 2/3 distribution feasibility is now executable through
  `primegaps.is_certified`. It checks the Harman-minorant inequalities, the
  global Type I/II/III hypotheses, and continuous support-cell partition
  conditions using exact rational order-statistic witnesses. The published
  support's full Cartesian product of `(j,m)` cells is certified (empty cells
  are identified vacuously), including its non-BV `A_1+A_1=0.506` pairs.

## Distribution statement endpoint issue

The printed Type IIc condition in Proposition 3 quantifies
`omega_0 in [-epsilon, omega(j,j')]` while requiring
`sum_{i in I_4} y_i <= 8 omega_0`. For negative `omega_0`, even the empty set
has sum zero and cannot meet that bound. The proof of the final theorem only
uses the two positive-capacity bins, so it implicitly treats the modulus range
at or below `x^(1/2)` by Bombieri--Vinogradov. The executable certificate makes
that split explicit: BV covers `omega_0 <= 0`, and the Proposition 3 Type IIc
partition is checked uniformly on `0 <= omega_0 <= omega`.

## What the paper does and does not specify

Section 5 confirms that all required integrals are represented as polynomials in
`delta` with coefficients `C_{m,i}` and `D_{m,i}` and that complicated mixed
small/large-coordinate integrals reduce to matrix products of these vectors.
It only sketches the multidimensional coefficient recurrence and says that the
full implementation will be uploaded later. The present reconstruction supplies
exact closed forms, a low-dimensional reference assembler, and the
symmetry-compressed assembly used for the exact `k=49` reproduction.

## Coupled-optimization plan

Do **not** optimize the historical stages independently. Once exact `M1/M2`
construction works, use the final generalized eigenvalue as the objective and
search jointly over:

1. support parameters `(delta, epsilon, A, B)`;
2. Harman parameters / retained prime mass `(xi_1,xi_2,xi_3,c1,c2)`;
3. test-function representation;
4. analytically certified modulus regimes.

The outer optimizer should ask which currently-forbidden support interaction
would raise the final eigenvalue most. That turns analytic theorem improvement
into a targeted separation-oracle problem instead of blindly maximizing a
single distribution exponent.

## Primary analytic target: `P3.II.delta`

The first measured analytic relaxation is now the primary target.  In the
fixed degree-21 two-band family with `delta=0.028`, support `epsilon=0.0085`,
and the published prime-indicator minorant, `P3.II.delta` binds exactly at
`A_max=0.2531666666`.  Relaxing it alone leaves a finite interval until
`P3.II.range` binds at `A_max=0.253777778055...`; no local witness becomes
binding first.

A 26-point translated-simplex stratified-QMC frontier gives an essentially
linear `A_max -> lambda_48` curve.  An independent local screen estimates the
crossing at `A_max=0.2536077308`, with randomized-QMC 95% interval
`[0.2536068027,0.2536086590]` and slope `3.94309718`.  This is a relaxation of
about `0.0004410642` beyond the unrelaxed analytic ceiling (or `0.0006077308`
beyond the paper's `A_max=0.253`).  At `A_max=0.2537`,
`lambda_48=1.0003636695`; at the next constraint it is `1.0006702180`.

The former Dirichlet-tilted finite screen is superseded: an unbounded
importance weight in its `m=2` denominator correction created a long right tail
and placed the crossing too far left.  The replacement samples each exact
large-coordinate-count stratum as a translated residual simplex with a
constant analytic volume weight.  These remain numerical screening values, not
an exact certificate or a global optimum.  See
[the full frontier and reproducibility notes](docs/p3ii-delta-frontier.md).

## Current gate / next milestone

The exact `k=49`, `D=21` gate is complete. With the same published support,
degree, rationalization rule, and frozen evaluator, `k=48` gives
`48J/I = 0.9969233513526357503888760066573328995...`, an exact deficit
`1-48J/I = 0.0030766486473642496111239933426671005...`.
Next:

1. turn the `P3.II.delta` crossing support into an exact/rational `k=48`
   certificate candidate;
2. determine what analytic input could supply the measured `0.0004411`
   relaxation without violating `P3.II.range`;
3. run D=22,23,... until threshold or D=27;
4. jointly optimize the simple support family, then enable the full Proposition
   3 / Harman-minorant degrees of freedom.

No result below 240 should be claimed until an exact/rational certificate has
been independently verified.

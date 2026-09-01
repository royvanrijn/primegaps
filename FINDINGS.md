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
- The complete generalized-eigenvalue *end of the pipeline* is implemented:
  diagonally equilibrated symmetric generalized eigenvalue search in floating
  point, conditioning/residual diagnostics, candidate-vector rationalization,
  and a portable exact certificate. Dense, matrix-free Lanczos, simultaneous
  block, and optional sparse paths are available. The standalone verifier checks
  semantic matrix hashes and replays both quadratic forms using Python integers.
  This mirrors the
  paper's strategy: floating point locates `c`; exact arithmetic proves it.
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
exact closed forms and a low-dimensional reference assembler. The remaining
engineering problem is symmetry-compressed assembly at `k=49`, not the rational
integral kernel itself.

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

## Current gate / next milestone

This repository still does **not** independently reproduce `k=49`. The next
milestone is now narrower: replace the reference status enumeration with a
symmetry-compressed assembler using the verified C/D kernel, then reproduce the
paper's intermediate matrix data or certificate when it becomes available.
After that:

1. reproduce the paper's `k=49`, D=21 quotient > 1;
2. same support, `k=48`, D=21;
3. D=22,23,... until threshold or D=27;
4. jointly optimize the simple support family;
5. enable full Proposition 3 / Harman-minorant degrees of freedom;
6. expand to multi-regime support geometry only if useful.

No result below 240 should be claimed until an exact/rational certificate has
been independently verified.

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
- The complete generalized-eigenvalue *end of the pipeline* is implemented:
  symmetric generalized eigenvalue search in floating point plus exact rational
  Rayleigh-quotient verification for a candidate certificate. This mirrors the
  paper's strategy: floating point locates `c`; exact arithmetic proves it.

## What the paper does and does not specify

Section 5 confirms that all required integrals are represented as polynomials in
`delta` with coefficients `C_{m,i}` and `D_{m,i}` and that complicated mixed
small/large-coordinate integrals reduce to matrix products of these vectors.
However, it only sketches the multidimensional coefficient recurrence and says
that the full implementation will be uploaded later. So the remaining task is
not numerical integration: it is reconstructing that omitted multidimensional
recurrence/decomposition exactly.

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
milestone is now narrower: extend the verified one-dimensional C/D base cases to
`k>1`, then assemble `M1` for a tiny basis and cross-check it against brute-force
high-precision integration in low dimensions. After that:

1. reproduce the paper's `k=49`, D=21 quotient > 1;
2. same support, `k=48`, D=21;
3. D=22,23,... until threshold or D=27;
4. jointly optimize the simple support family;
5. enable full Proposition 3 / Harman-minorant degrees of freedom;
6. expand to multi-regime support geometry only if useful.

No result below 240 should be claimed until an exact/rational certificate has
been independently verified.

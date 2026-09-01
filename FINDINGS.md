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
  Stadlmann explicitly says the more general results were developed with better
  `H_1` bounds in mind.
- The final numerical certificate is a vector `c` with a generalized Rayleigh
  quotient `c M2 c^T / (c M1 c^T) > 1` for `k=49`.

## New computational observations in this repository

The symmetric basis size grows from **846 functions at D=21** to **2526 at
D=27** for `k=49`. A dense matrix therefore grows from about 0.72M to 6.38M
entries (~8.8x); this explains part of the resource jump even before accounting
for the much more expensive exact integral construction.

A uniform-volume Monte-Carlo sample inside the outer simplex almost never sees
Stadlmann's `B` constraint at `k≈49`: most coordinates are much smaller than
`delta=0.028`. Consequently a naive objective such as "maximize support volume"
is a poor surrogate for the actual sieve objective. The important mass is set by
the test function / `I,J,K` integrals, not Euclidean volume. This diagnostic is
implemented in `primegaps.scan` and intentionally carries no bound claim.

## Coupled-optimization plan

Do **not** optimize the historical stages independently. Once exact `M1/M2`
construction works, use the final generalized eigenvalue as the objective and
search jointly over:

1. support parameters `(delta, epsilon, A, B)`;
2. Harman parameters / retained prime mass `(xi_1,xi_2,xi_3,c1,c2)`;
3. test-function representation;
4. analytically certified modulus regimes.

The outer optimizer should be allowed to ask which currently-forbidden support
interaction would raise the final eigenvalue most. That turns analytic theorem
improvement into a targeted separation-oracle problem instead of blindly
maximizing a distribution exponent.

## Critical missing milestone

This repository does **not** yet independently reproduce the `k=49` result. The
paper only sketches the Section 5 coefficient recurrences for the exact integrals
and says the author's code will be published later. Reconstructing those
recurrences and obtaining the paper's `M1,M2` eigenvalue is the next gate.

After that, the first experiments should be, in order:

1. same support, `k=48`, D=21;
2. D=22,23,... until the threshold or D=27;
3. optimize the simple published support family;
4. enable the full Proposition 3 / Harman-minorant degrees of freedom;
5. only then expand to multi-regime support geometry.

No result below 240 should be claimed until an exact/rational certificate has
been independently verified.

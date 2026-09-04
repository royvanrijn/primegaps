# PrimeGaps186 audit and research pivot

Status date: 2026-09-03.

Upstream snapshot: [`openai/PrimeGaps186@61340d0b74163003b32756bb16e91d9209a5e330`](https://github.com/openai/PrimeGaps186/tree/61340d0b74163003b32756bb16e91d9209a5e330).

## Executive conclusion

The repository contains a serious and highly novel route to

```text
DHL[40,2] -> H_1 <= 186.
```

It is a much stronger architecture than the Stadlmann/BGP212 construction and should become the
primary research basis. It is **not yet a self-contained unconditional record in its public form**.
The final Lean declarations are conditional on three explicit axioms:

1. a normalized rank-three hyper-Kloosterman bound;
2. a rank-two Kloosterman correlation bound;
3. 149 physical source-integral inequalities and three cap bounds.

The first two are recognizable normalizations of published Deligne/Katz and
Fouvry--Kowalski--Michel estimates. The third is backed by a substantial directed-Arb Python
calculation, but that calculation has not been converted into a Lean proof and no passing JSON
receipt is committed. The formalization metadata describes the result as conditional and
self-assessed, with no independent human semantic review.

Our immediate job is therefore twofold:

1. independently replay and close the three input boundaries;
2. use the physical-fragment construction, rather than the now-superseded 212 support, as the basis
   for the next target `DHL[39,2] -> H_1 <= 182`.

## Static formalization audit

At the pinned revision:

- `PrimeGaps186.lean` is about 10 MB / 210,637 lines;
- it imports Mathlib and contains the full conditional sieve-to-gap chain;
- a static whole-file scan finds exactly three declarations beginning with `axiom`;
- it contains no `sorry` tokens, no `unsafe` declarations and no `opaque` declarations;
- the final declarations are:
  - `PrimeGap186.dhl_40_2`;
  - `PrimeGap186.infinite_two_prime_translates_admissibleTuple`;
  - `PrimeGap186.primeGapLiminf_le_186`.

Comparator, Nanoda and Lean's kernel reportedly accept those conditional proofs. This establishes
internal proof consistency relative to the three inputs. It does not by itself establish that the
formal definitions faithfully encode every intended analytic statement, nor that the numerical
input is true.

The human-readable source `Improved Gaps Between Primes`, repeatedly cited as the accompanying
main paper, is not present in the repository and was not found in the public search performed for
this audit. Until it appears, semantic review has to work from the Lean development and the
23-page numerical note.

## The construction

This is not a larger Stadlmann support or a higher-degree version of the same polynomial sieve.
It changes the state space.

### Physical prime-factor law

A sieve coordinate is represented by a finite multiset of logarithmic prime fragments. The
limiting fragment law is a scale-invariant Poisson process on `(0,c]` with intensity `dv/v`; its
total-size pushforward is expressed through the Dickman density. Erasing one coordinate means
integrating over this physical fragment law.

This retains factorization information that a single scalar coordinate and a worst-case support
polytope discard. Source failures can therefore be charged according to their actual fragment
configuration instead of excluding a whole support interaction.

### Fixed trial

The numerical trial has:

```text
k                         = 40
mesh cells per coordinate = 98,304
angular signatures        = 11
radial degrees             = 0..6
rational coefficients     = 77
rho*                       = 2624989 / 10000000 = 0.2624989
```

The angular signatures are

```text
empty, (2), (3), (4), (5), (6), (2,2), (2,3), (2,4), (3,3), (2,2,2).
```

Each is multiplied by a degree-at-most-six radial polynomial and by the positive two-pole profile

```text
g(t) = (21/200)/(1+t/100) + (179/200)/(1+(907/5)t).
```

The effective pair-modulus exponent is nearly

```text
2 rho* = 0.5249978.
```

That is qualitatively beyond the `0.514` generated level in the BGP212 datum, but it is achieved for
this structured physical source problem rather than as a generic level of distribution.

### Source ladders and signed restoration

The proof combines:

- an old/full-prime source ladder;
- a new minorant source ladder;
- a prime minorant with limiting mass deficit below `1/50000`;
- a signed hybrid restoration across base, enlarged and unrestricted face regions;
- 97 positive source-loss components grouped by role and source order.

The hybrid coefficients are exact rationals. The negative unrestricted-face coefficient is retained
throughout, while every support/source failure is covered by an explicit nonnegative loss term.
This is a concrete instance of the non-greedy strategy we were seeking: accept a tiny minorant
deficit and signed bookkeeping in order to unlock a substantially stronger distribution geometry.

## Numerical certificate boundary

The numerical note prints

```text
23685317816e-24 <= I_H <= 23685317890e-24
J_lambda,H >= 90248755123e-24
L+ / I_H^- <= 696075110e-12 < 697e-6
rho* J_lambda,H / I_H > 500103/500000
```

and concludes with the exact coarse margin

```text
500103/500000 - 1 - (2624989/10^7)(697/10^6)
  = 230382667/10^13
  > 1/50000.
```

Replay the printed rational arithmetic with

```bash
python scripts/check_primegaps186_margin.py
```

This does not recompute the integrals. The upstream certificate script does that from the embedded
exact inputs using:

```text
160-bit Arb arithmetic
224-bit cap convolutions
192-bit source convolutions
NumPy binary64 reductions with explicit directed-error accounting
python-flint 0.9.0
FLINT 3.6.0 plus the post-release signed-FFT correction
```

The correction is important: stock FLINT 3.6.0 predates the merged fix for a rare wrong-result bug
in signed `fmpz_poly_mul_SS`. The script contains a mandatory regression test and refuses to run on
the defective build.

## Where the current certificate loses margin

The six printed group totals sum to `696075110` units at scale `10^-12`:

| group | units | share |
|---|---:|---:|
| outer order `5/2` | 622,829,241 | 89.4773% |
| outer order `2` | 38,927,522 | 5.5924% |
| new-inner order `5/2` | 32,422,390 | 4.6579% |
| new-inner order `2` | 1,405,159 | 0.2019% |
| old-inner order `5/2` | 435,544 | 0.0626% |
| old-inner order `2` | 55,254 | 0.0079% |

By component type:

```text
rank-two: 602,422,937 units = 86.55%
low:       79,845,205 units = 11.47%
high:      13,806,968 units =  1.98%
```

The largest individual losses are almost all outer order-`5/2` rank-two components. This is the
right target for making the 186 certificate more robust and cheaper to certify.

It is **not**, by itself, the route to the next `k`. Eliminating every printed source loss while
freezing all other ingredients improves the quotient by less than

```text
rho* * 696075110e-12 = 0.000182719.
```

The frontier therefore has to move through the trial, physical masks, source reach or minorant
tradeoff, not merely through tighter numerical bounds on the existing 97 terms.

## Trust and replay ledger

### A. Finite tuple

The included 40-tuple is

```text
{0,2,6,12,20,26,30,32,36,42,48,50,56,60,68,72,78,86,90,92,
 98,102,110,116,120,126,132,138,140,146,152,156,158,162,168,
 170,176,180,182,186}.
```

It has cardinality 40, diameter 186 and is admissible. Its 39-element prefix is also admissible and
has diameter 182, giving the exact finite endpoint for the next target. Replay both with

```bash
python scripts/verify_h40_tuple.py
```

### B. Finite-field inputs

The two finite-field axioms look dischargeable from existing literature:

- Katz's Kloosterman sheaf has rank `n` and weight `n-1`; for `n=3`, the raw trace bound is `3p`,
  matching the repository's normalization by `p`.
- Fouvry--Kowalski--Michel prove `|S(alpha,beta)| <= 8 sqrt(p)` for normalized classical
  Kloosterman sums. Undoing the two `sqrt(p)` normalizations gives the repository's
  `8 p sqrt(p)` bound.

These normalizations still deserve a line-by-line formal audit, especially the excluded poles,
zero parameters and conventions for complex norm.

### C. Physical numerical input

This is the main computational gate. A credible replay should:

1. build a pinned FLINT revision containing the signed-FFT fix;
2. run the upstream script from a clean environment and retain the complete fresh JSON receipt;
3. check all 97 tasks, 149 raw forms, exact masks and final hashes;
4. reproduce critical convolutions with an independent backend, preferably fixed-denominator
   integers plus NTT/CRT or a second Arb implementation;
5. turn the receipt into a small deterministic verifier and eventually a Lean theorem replacing
   `physical_integral_bounds`.

### D. Semantic proof review

The Lean file is internally checked but the metadata explicitly records no independent semantic
review. The audit should focus first on the interfaces connecting:

```text
finite-field estimates
  -> source distribution bounds
  -> physical fragment realization
  -> signed restoration
  -> positive weighted first moment
  -> DHL[40,2].
```

Comparator establishes theorem-statement matching, not fidelity of every intermediate definition
to the intended paper mathematics.

## What remains reusable from this repository

The old Stadlmann/BGP212 functional and the new physical-fragment functional are not directly
interchangeable, so the D27 coefficient vector and Type-IIc endpoint do not simply plug in.

The reusable machinery is still substantial:

- manifest and hash binding;
- exact rational endpoint replay;
- outward-rounded Arb certification;
- checkpointed long-running computations;
- independent tuple verification;
- candidate-space / matrix-free optimization;
- shadow pricing of analytic and source constraints;
- sparse symmetric feature handling;
- independent backends for exact polynomial convolution.

The old D27 run should finish as a validation of this infrastructure, but further record-oriented
optimization should start from PrimeGaps186.

## Frontier plan: `k=39`, hence `H_1 <= 182`

### 1. Parameterize and measure `k=39` before proving anything new

The present code and Lean proof hardcode 40. Generalize the numerical engine to `k` and evaluate a
39-coordinate analogue using the same source theorem and a re-optimized trial. The 39-element
prefix of the displayed tuple already supplies diameter 182.

The first output must be a trustworthy deficit curve, not a guessed `40/39` scaling. In Maynard-type
variational problems the optimal function changes with dimension, so fixed-vector scaling can be
quite misleading.

### 2. Turn the entire restored certificate into a quadratic operator

The first decisive part of this experiment is complete. The cap forms have
been assembled in the actual 77-dimensional basis and combined with the exact
signed hybrid coefficients:

```text
H = J0 + (a+b) J+ + b Jt.
```

Inverse-mesh extrapolation gives `0.9943709810` for `k=39` and
`0.9883950587` for `k=38`, before source losses. Because all 97 loss forms
below are positive semidefinite and subtracted, `H` is already a decisive upper
screen: the full frozen-geometry operator cannot cross one. Accordingly, the
expensive assembly of the 97 individual loss matrices is not warranted for
this no-go decision. The exact calculation, calibration, and replay boundary
are in [the restoration and factorization note](physical-restoration-factorization.md).

For fixed profile, masks, source ladders, Young parameters and hybrid coefficients, every cap form
and every one of the 97 loss terms is quadratic in the 77 trial coefficients. Build matrices

```text
I, J0, J+, Jt, L1, ..., L97
```

and optimize the actual restored objective

```text
rho* (J_lambda,H - sum L_i) / I,
```

not the raw cap quotient. This is a 77-dimensional generalized eigenproblem before basis expansion,
so it should be far cheaper than the current from-scratch scalar certificate.

That full loss-matrix build remains useful only after a geometry or
hybrid/source change raises the signed cap screen above one; at the present
parameters it can only strengthen the negative result.

Then add angular signatures and radial degrees adaptively using a Davidson/column-generation loop.
The trial space is compact compared with the 846/2526-term Stadlmann spaces, leaving obvious
variational freedom to test.

### 3. Optimize the upstream trade rather than only the trial

Jointly vary:

- `rho*` and the old/new source ladders;
- the permitted minorant deficit;
- hybrid coefficients;
- radial shell and largest-fragment caps;
- the two profile poles and weights;
- Young parameters for the dominant outer order-`5/2` rank-two losses.

The present `1/50000` minorant deficit is exceptionally small. A deliberately weaker minorant may
be globally better if it allows a materially larger `rho*` or a less restrictive physical mask. The
objective should be the final restored quotient, producing a Pareto frontier

```text
minorant mass loss <-> source reach <-> physical support <-> final quotient.
```

This is the strongest concrete incarnation yet of the project's original “de-optimize an earlier
step to unlock a larger later gain” thesis.

## Target hierarchy

```text
close/replay DHL[40,2] -> H_1 <= 186
DHL[39,2]              -> H_1 <= 182
DHL[38,2]              -> H_1 <= 176
```

The immediate research target is 182. The immediate verification target is making 186 independent,
reproducible and no longer conditional on an opaque numerical axiom.

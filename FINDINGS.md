# Findings / reconstruction status

Date: 2026-09-02

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
- Exact J functionals can now be cached independently of the candidate and
  polynomial degree. A fully cached replay does not rebuild target densities or
  reintegrate support geometry; higher-degree runs request only absent exponent
  pairs. This removes a major repeated cost, but it is an acceleration result,
  not a new variational certificate.

## Negative and inconclusive `k=48` searches

The higher-degree and support searches completed after the original status
write-up do not establish a bound below 240:

- An unrestricted deterministic degree sweep (which omits the published `B`
  cutoffs) reaches `0.9982088741, 0.9985996972, 0.9992942888, 0.9996691406,
  1.0000734506, 1.0003862453` at `D=22,...,27` under its aggressive spectral
  cutoff. The apparent first crossing at `D=26` is therefore a discovery proxy,
  not a value on the legal support. The corresponding randomized cutoff
  correction and a floating-point exact replay were refuted by impossible
  quotients and catastrophic cancellation. Exact `I` for one `D=24` candidate
  is positive, but exact `J` is incomplete at `224/7338` groups.
- Column generation from the exact `D=21` candidate added 1,032 of the 1,680
  unused columns through `D=27`; its unrestricted proxy stopped below one at
  `0.9980139970`. An earlier nominal `D=22` candidate evaluates exactly to
  `48J/I=0.2686576479...`, but that vector was produced with an incorrect
  signature-dependent Jacobi conversion and is refuted as a calibration of the
  intended eigendirection.
- Compact support-adapted trial spaces depending only on total mass, large-count
  stratum, and large-coordinate excess were not competitive. Their best stable
  checked degree-4 screen was about `0.84`, far below the exact global
  degree-21 `k=49` quotient.
- A legal irregular-support search at degree 21 found no geometry that closes
  `k=48`. At the exact `P3.II.delta` ceiling the strongest recorded legal screen
  remains `lambda_48=0.9982613325` with replicate SE `3.0018e-6`. The former
  `delta=0.014` lead is now rejected as a heavy-tail estimator artifact.

## Distribution statement endpoint issue

The printed Type IIc condition in Proposition 3 quantifies
`omega_0 in [-epsilon, omega(j,j')]` while requiring
`sum_{i in I_4} y_i <= 8 omega_0`. For negative `omega_0`, even the empty set
has sum zero and cannot meet that bound. The proof of the final theorem only
uses the two positive-capacity bins, so it implicitly treats the modulus range
at or below `x^(1/2)` by Bombieri--Vinogradov. The executable certificate makes
that split explicit: BV covers `omega_0 <= 0`, and the Proposition 3 Type IIc
partition is checked uniformly on `0 <= omega_0 <= omega`.

## Rejected external claim: Shi (2025), `H_1 <= 234`

Yuhang Shi's July 2025 ResearchGate preprint *A Weighted Distribution of
Primes and a New Unconditional Bound on Gaps Between Primes* has been audited
and should **not** be treated as an input or revisited in its present form.
The negative conclusion is mathematical, not sociological:

- Lemma 4.3 is false as stated. At `sigma=1/2`, standard zero counting over
  its semiprime-modulus primitive-character family gives
  `N*(Q,1/2,T) >> Q^2 T/log Q`, contradicting its claimed
  `Q^2 T exp(-c sqrt(log Q))` upper bound.
- The Type-II character factorization drops the condition `mn<=x`; Appendix B
  gives an unnormalized zero-detector inequality, then drops its `L(s,chi)`
  factor, and asserts the required exponential saving only in prose. Lemma 4.4
  also changes a fixed `delta(C)>0` into `delta >> 1/(C log q)`.
- Its single rough-times-prime convolution is not a proof of standard
  well-factorability, and no lemma matches its sparse modulus weight to the
  least-common-multiple coefficients generated by the Maynard sieve.
- The cited Polymath8b source contains no rigorous `k=48, 3.998` certificate,
  and the paper uses mutually inconsistent normalizations of the Maynard
  threshold.
- An independent MILP optimization proves that an admissible subset of
  `[0,234]` has at most 47 elements. Diameter 236 admits 48 and reproduces the
  maintained MIT tuple. Thus the printed ellipsis cannot hide the claimed
  48-tuple of width 234.

The convolution, rough-number pre-sieving, Parseval step, and
diagonal/off-diagonal vocabulary are standard ingredients. Weighted levels
beyond `1/2` are already known at stronger exponents for weights that genuinely
match their sieve applications. Nothing specific to this preprint remains as a
promising research lead.

Reopen this line only for a substantially revised manuscript that supplies all
of: explicit `R,S` ranges; a correct zero-density theorem and proof; a
coefficient-level Maynard support transfer; a reproducible `M_k` certificate;
and a complete admissible tuple consistent with the claimed endpoint. See
[the detailed audit](docs/shi-2025-audit.md). The finite endpoint check is
reproducible with `scripts/audit_shi_tuple.py`.

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

The analytic loss behind the active branch has also been traced. It comes from
the final `m` term in the Type-IIc `wts3` estimate, through
`48*omega+16*delta_star-4*gamma<-1`; it is not removable by epsilon or dyadic
bookkeeping. Reaching the numerical crossing requires an exponent saving
`0.0211710816` in that term. Uniformly improving the complete exponential sum
or the outer `H^4` pair count is impossible by the recorded sharp examples, so
any useful lemma must retain cancellation in the actual incomplete ranges
before the current absolute-value/supremum steps. No such lemma is proved.

A surgical branch-deletion alternative has now also been rejected inside the
current Buchstab identity. The Case-IIc demand is already generated by a
negative direct Type-II semiprime correction and a positive direct Type-II
middle branch; it is not caused by either exceptional five-prime loss in the
minorant catalog. The negative term cannot be deleted without making the
candidate positive on rough semiprimes. Even granting a free reorganization of
that obstruction, an exactly identified high-sum subset of the positive term
has rigorously enclosed normalized mass
`[0.071588483537845070, 0.071771543742745046]`. Hence
`rho<0.928411516462154930`, and the optimistic no-`K` score at the measured
far endpoint is below `0.929033754571895163`. A survivor would need a raw score
above `1.077108569064980194`. No full `I/J/K` run is warranted for deletion-only
variants of this decomposition.

## Current gate / next milestone

The exact `k=49`, `D=21` gate is complete. With the same published support,
degree, rationalization rule, and frozen evaluator, `k=48` gives
`48J/I = 0.9969233513526357503888760066573328995...`, an exact deficit
`1-48J/I = 0.0030766486473642496111239933426671005...`.
Next:

1. build a numerically stable full-support matrix constructor, reusing the exact
   candidate-independent I/J moment caches, so higher-degree eigenspaces can be
   evaluated without the refuted monomial-basis float replay;
2. use that constructor to resolve the incomplete `D=24` candidate or generate
   a better legal-support candidate before spending on `D=25--27` exact jobs;
3. investigate the actual short-rectangle Type-IIc completion before the
   triangle inequality and either prove the required averaged saving or record
   a counterexample;
4. only after those gates, jointly optimize support and Harman-minorant degrees
   of freedom.

No result below 240 should be claimed until an exact/rational certificate has
been independently verified.

The ledger-wide classification, including superseded and negative experiments,
is recorded in [the research status audit](docs/research-status.md).

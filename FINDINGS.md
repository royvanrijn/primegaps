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

- A first sparse positive-semidefinite sum-of-squares screen used an
  oracle-derived induced four-cycle, 60 localized D7 components, and `2^16`
  scrambled-Sobol points per exact-count stratum. At both `k=47` and `k=46`
  the full SDP optimizer had numerical rank one and matched the best legal
  rank-one clique to within `5.1e-9` relative; it therefore showed no advantage,
  far below the proposed `2.2%` signal threshold. A signed four-cycle control
  does produce the expected rank-two `sqrt(2)` advantage, so the null result is
  not a solver degeneracy. This rejects only the first D7 bank, not all sparse
  support graphs; see [the sparse SOS screen](docs/sparse-sos.md).

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

Subsequent fixed-candidate calibration found that uniform sampling inside those
translated strata is itself not a reliable high-degree threshold estimator. At
the independently audited safe endpoint and degree 21, `2^15` samples per
stratum gave normalized `kJ=1.0000632111` instead of the independently
replicated importance-control value `0.9999667252`, and estimated
`I=0.9999951011` instead of the exact `0.9989116509`. The exact-count partition
is valid, but the rare high-leverage polynomial tail inside a boundary stratum
still needs importance sampling or deterministic integration. The published
translated-simplex frontier must therefore be treated as superseded as a
quantitative crossing estimate pending recalculation.

## Candidate-space J accumulation

The numerical J bottleneck is no longer target routing or persistent
signature-pair blocks. Each integration batch now forms the candidate feature
matrix `G=FM` and performs one candidate-space Gram update; a projected run
forms `GQ` before accumulation. At degree 27 this reduces the stored
upper-triangle from 31,546,512 marginal-feature entries to 3,191,601 candidate
entries. A real 256-row batch took 69 ms for feature evaluation, 11 ms to form
`G`, and 28 ms for the Gram update. Four complete D27, `2^15` legal-minus-
unrestricted corrections took 13.9--14.8 seconds each after a one-time analytic
unrestricted build.

Raw full D27 matrices are still too ill-conditioned to optimize directly.
Projecting each batch before accumulation and cross-validating nested
unrestricted eigenspaces gives a stable four-direction subspace: with the D21
direction included, the four-seed training value is `1.0004338031` and the
leave-one-out minimum is `1.0003994559`. Prefix eight already overfits. The
rationalized prefix-four candidate has four independent `2^17` estimates of
normalized `kJ` between `1.0002556783` and `1.0004605125` (mean
`1.0003688580`, standard error `4.40e-5`). Its exact normalized unrestricted-
simplex `I` is `0.9999955585079013`. Since the legal support is a subset of that
simplex and `F^2` is nonnegative, this is an exact upper bound for legal `I`;
the abandoned legal-`I` checkpoint is not on the certification critical path.
This is still a numerical `J` discovery, not a certificate: exact legal `J`
remains the final gate. The exact target-free unrestricted contraction gives
normalized `kJ_simplex=1.0005162874604419` in 442.5 seconds. Consequently the
only remaining exact variational calculation is the B-boundary correction,
which must be greater than `-0.0005207289525405`. The four independent
numerical corrections range from `-0.0002606091` to `-0.0000557749`.

The same `p^T H q` rewrite is exact on the certificate path, but it is not yet
a cache-fill speedup by itself. On the old 23-pair D21 benchmark, a Python/GMP
Hankel contraction did not finish in 390 seconds and a Sage/FLINT-matrix form
did not finish in 180 seconds. The measured cost is construction of a separate
density-weighted moment table for every pair and cell. The next exact change is
therefore to persist raw moments once per geometry cell and batch every density
correlation from that table. The direct target and per-pair fallbacks must not
be used for the full D27 correction.

The analytic loss behind the active branch has also been traced. It comes from
the final `m` term in the Type-IIc `wts3` estimate, through
`48*omega+16*delta_star-4*gamma<-1`; it is not removable by epsilon or dyadic
bookkeeping. Reaching the numerical crossing requires an exponent saving
`0.0211710816` in that term. Uniformly improving the complete exponential sum
or the outer `H^4` pair count is impossible by the recorded sharp examples, so
the successful argument retains the originating affine structure and
cancellation in the actual incomplete ranges before the old supremum step.

The resulting checked terminal condition is
`6-22*gamma+72*delta+216*omega<0`, equivalently
`delta<11*gamma/36+2/3-3*A`. Its asymptotic ceiling is
`856/3375=0.2536296296...`; the independently hostile-audited safe endpoint is
`2029/8000=0.253625`, with unpadded structural margin `0.001`. The proof uses
large-common-divisor affine sparsity, a rank-three Kloosterman shifted-
correlation estimate with Plancherel averaging, exact nonunit and `q0`
conductor bookkeeping, pre-lift Möbius removal of `w1`, and an explicit
`q3t=(q3,(w2,m))` charge for phase-dead masks. See
[the Type-IIc theorem note](docs/typeiic-incomplete-rectangles.md). This is an
analytic result; the `lambda_48>1` screen is still numerical and the production
oracle remains unchanged pending a fully typeset human review.

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

## Minimum-breakthrough inversion (D21)

A simultaneous-slack optimizer now evaluates every registered Proposition 2/3
condition rather than waiving one condition at a time. With unit raw-exponent
weights, fixed `delta=0.028`, `epsilon=0.0085`, and prime-indicator minorant
`(0.38,0.4,0.4)`, the independently validated D21 grid leaders are: `k=47`,
`A=0.2548` with the current B profile (cost `0.048633324`); `k=46`, `A=0.2560`
with the current B profile (cost `0.109833324`); `k=43`, `A=0.2596` with outer
`B_2=0.16`, `B_m=0.18` for `m>=3` (cost `0.319379991`); and `k=42`, the same B
relaxation at `A=0.2610` (cost `0.399333324`). Four independent `2^17`-point
importance-control replicates put every selected direction above one, and the
respective mean-minus-two-SE values are `1.0003471032`, `1.0002649501`,
`1.0004399482`, and `1.0003383769`.

The unchanged-B candidates at the first apparent crossings for `k=43,42`
failed independent validation, while the simultaneous B relaxation at the same
endpoint passed; this reverses their discovery-only ranking. Perfect removal
of local support restrictions at the current asymmetric endpoint remains below
one for all four k. At unrestricted symmetric exponent `theta=0.525`, full D21
clears only `k=47,46`; `theta=0.60` and `66/107` clear all four. These are
finite-grid numerical leads, not theorem certificates or global optima. D27 is
excluded from this benchmark. Full vectors and controls are documented in
[the minimum-breakthrough note](docs/minimum-breakthrough.md).

## Physical parity viability

A genuinely different prime detector passes the idealized full-face numerical
gate, but not the frozen signed-restoration gate.
For `beta=0.250001`, a beta-rough target has at most three prime factors, so
`1_{Omega=1}=Omega-2*C(Omega,2)+3*C(Omega,3)`.  The associated leading rough
constants are `3.6388648513`, `1.5402578959`, and `0.1472169802`; their signed
combination is exactly one, but their absolute signed condition number is
`7.1610315838`.

Adding those three arithmetic states to the frozen PrimeGaps186 full physical
face and reoptimizing the 77-dimensional trial gives production-mesh
predictions `1.0192774814` at `k=39` and `1.0127855640` at `k=38`.  The final
normalized parity errors may therefore be at most `0.0192774814 I` and
`0.0127855640 I`; if all three unsigned terms have the same relative error,
the limits shrink to `0.264108%` and `0.176290%`.  A four-mesh direct-
convolution calculation supports the extrapolation, but it is not a rigorous
enclosure. This idealized calculation deliberately omits the signed hybrid
restoration and is not the operative frozen-geometry score.

The missing signed operator has now been assembled in the same 77-dimensional
basis:

```text
H = J0 + (a+b) Jplus + b Jtail,
a = 2479900401/2500000000,
b = -843183/1000000000.
```

Its production-mesh extrapolations are `1.0001866542`, `0.9943709810`, and
`0.9883950587` for `k=40,39,38`, respectively. Separate optimization at every
mesh gives `1.0002083815`, `0.9943911126`, and `0.9884136585`. The `k=40`
cross-check is only `2.3021e-6` above the published rigorous fixed-vector lower
endpoint. Since every one of the 97 source-cover forms is nonnegative and is
subtracted after restoration, this loss-free matrix is an upper screen for the
full fixed-geometry operator. Therefore `k=39` and `k=38` are already below one
before source losses: trial-only optimization is a numerical NO-GO, and the
geometry or hybrid/source parameters must move.

The exact source-ladder oracle reconstructs 29 old and 43 new rows and compares
the union of exact order-three failures with grouped nonlargest `H_{5/2}` after
the separate largest-fragment and opposite-root guards. Across 373,857 exactly
classified critical and seeded cell configurations it records no false
negative, and it provides a strict-overcoverage witness for each of the outer,
old-inner, and new-inner groups. Thus `H_{5/2}` is not equivalent to the exact
row predicate; it is the safe majorant used by the source proof. The finite
census is not itself a theorem or a physical-law probability calculation. See
[the restoration and factorization note](docs/physical-restoration-factorization.md).

The generated-modulus reach has now been resolved numerically for `k=39`.
Writing `E_{i,theta}` for the erased-coordinate operator restricted by
`rho*(outer total + retained-face total) <= theta`, and
`Lambda_39(theta)=rho*lambda_max(sum_i E_{i,theta}^*E_{i,theta},I)`, the
production extrapolation gives `0.9787438843` at `theta=1/2`, `0.9913935182`
at `0.51`, and `1.0025203787` at `0.52`.  The primary crossing bracket is
`[0.517625,0.517750]`, with interpolation `0.51766554`; independent fitting
choices range from about `0.51760` to `0.51767`.  The answer is therefore the
middle regime: strictly above one-half, but far short of the full `0.5485994`
envelope.  This is a float64 mesh experiment, not a certificate or a theorem.

After expanding the squared physical Selberg weight, the direct signed target
requires Liouville cancellation across its structured, polynomial-size lcm
congruences.  The initially proposed Friedlander--Iwaniec `B_F` route for the
already rough nonnegative sequence is false.  In every nonempty active block,
the inner variable is forced prime; at the required value `C=1`, its Mobius
sign is constant for fixed outer variable and the outer absolute value gives a
positive mass rather than cancellation.  The proposal also omitted the
`(log x)^-222` saving, the physical coefficient norm, and the separate
`D>x^(2/3)` remainder axiom.  The complete indexed block family and
classification are in [the dyadic audit](docs/bf-dyadic-audit.md); the
numerical derivation and direct target `(P)` remain in
[the physical-parity note](docs/physical-parity-viability.md).

The roughness endpoint has now been swept rather than frozen.  The ideal score
stays unchanged, while the direct factorial-moment condition number decreases
monotonically from `7.1610315838` at `beta=0.250001` to `3.7725887222` at the
`1/3` transition.  Strictly above `1/3`, the exact detector becomes
`N-2*C(N,2)=(1-lambda_L)/2`: the triprime constant and its prime-times-semiprime
bilinear block disappear.  At `beta=1/3+epsilon` the k=39 and k=38 common
relative-error budgets are approximately `0.5013%` and `0.3346%`, about 1.9
times the quarter-rough budgets.  At `beta=0.4` they rise to `0.7214%` and
`0.4815%`, with condition number `2.6218604324`.  These gains cost the deeper
condition `P^-(m)>x^beta`; the recorded physical pair-support-to-beta ratio
falls from `2.1944` to `1.6458` at the transition and is only a scale
diagnostic, not a distribution theorem.  The first analytic target is therefore
just above `1/3`, with `0.4` a secondary target if the rough sieve has enough
depth.

A reduced nondegenerate discrete compiler now evaluates the frozen `k=39`
physical trial on 211 divisor configurations and maps 8,455 compatible
`(d,e)` pairs to 256 residue-coloured `(q,a)` states.  With unit discrete
`l2` normalization it gives `sum_(q,a)|c(q,a)|^2=0.4152968949`; CRT collision
aggregation reduces coefficient `l1` by a factor `3.638754`.  Across
`X=8,16,32,64` million, the requested blockwise parity cancellation ratio is
`54.42--115.47` before aggregation and `8.51--37.36` afterwards.  The latter
is persistent, but it comes almost entirely from cancellation across distinct
CRT blocks: the within-`(q,a)` prime/semiprime ratio is only `1.026--1.041`.
The literal ratio formed from the two fully summed scalar errors is exactly
one at every scale because their signs oppose and the projected minus sign
reinforces them.  The coefficient table has stable rank `2.205`; its top eight
singular directions contain `99.629%` of energy, while the projected error
table's top eight contain `86.09--94.33%`.  This is a qualified GO for studying
one global CRT-coloured dispersion operator and a NO-GO for local aligned
prime/semiprime cancellation.  See
[the finite experiment](docs/physical-parity-crt-finite.md).

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
3. typeset and human-review the checked short-rectangle Type-IIc theorem, then
   promote only the safe rational endpoint into the distribution oracle;
4. only after those gates, jointly optimize support and Harman-minorant degrees
   of freedom.
5. on the standard physical-source branch, move the frozen geometry or
   hybrid/source parameters until the signed-restoration screen crosses one;
   independently, on the full-face parity branch, pursue a theorem for the
   now-specified global CRT-coloured projected operator. Do not rely on local
   prime/semiprime cancellation, and test the observed across-CRT signal
   against larger compressed models and residue/sign controls.

No result below 240 should be claimed until an exact/rational certificate has
been independently verified.

The ledger-wide classification, including superseded and negative experiments,
is recorded in [the research status audit](docs/research-status.md).

# Impact of the September 2026 `H_1 <= 212` draft

Status date: 2026-09-03.

Primary source: F. Charton, L. Hong, K. Lau, K. Ono, G. Remy, H. C. Siu,
A. A. Swaminathan, J. Thorner, and Y. Xie, *A New Bound for Small Gaps
Between Primes*, preliminary draft dated 2026-09-03.

The fetched 46-page PDF has SHA-256
`2b307ae26046dffd2153fe7466e5a9a302976dc3ba5142151a8a8104a37b3d80`.
This repository does not mirror the manuscript; use the authors' published URL.

## What changed

The new draft proves `DHL[45,2]` and combines it with an explicit admissible
45-tuple of diameter 212. Its exact rescaled variational witness satisfies

```text
J_T(F) / I_T(F) = 4.00438409833460131937... > 4.
```

The paper still uses the direct prime indicator, a degree-21 symmetric rational
polynomial space of dimension 846, and a one-band Stadlmann support. Its main
advance is not a new generic distribution exponent. It preserves the exact
divisors consumed by the Type I/II/III proofs, verifies factor extraction over
the full continuous rough-factor domain, and explicitly bridges the transition
through the half-level.

The parameter datum is

```text
k = 45
omega = 7/1000
A = (-1/125, 257/1000)
support epsilon = 1/125
delta = 41/2500
xi = (19/50, 2/5, 2/5)
analytic epsilon = 10^-10
B_m = (777,794,875,917,953,983,1016,1042,1063,1081)/5000, m=1..10
B_m = 1081/5000, m>=11
```

Thus the total support cap is `A_1+epsilon=53/200=0.265`, the marginal cap is
`A_1-epsilon=249/1000=0.249`, and generated positive-level moduli reach
`x^(2*A_1)=x^0.514`.

## What is superseded

The in-progress `k=48` project can no longer produce a record at 236. It remains
valuable as an independent validation of the exact integral engine, projected
higher-degree search, and boundary-only Arb certifier, but 212 is now the
external target to beat.

The new draft also subsumes several research directions which had been open in
this repository:

- retaining unsimplified divisor intervals instead of uniform worst-case bounds;
- exact continuous packing rather than checking one extremal profile;
- a complete transition argument around `x^(1/2)`;
- optimized nonconstant rough caps `B_m`.

These should be imported as the new arithmetic baseline rather than rediscovered.

## What still composes with our work

### Higher degree and projected discovery

The draft uses degree 21. Our directly evaluated orthogonal-feature pipeline,
projected optimizer, and final scalar Arb certifier are designed for degrees up
to at least 27. They remain applicable after replacing the old support by the
new parameter datum.

Degree alone should not be assumed to yield `k=44`, but it is genuine unspent
variational freedom. The right experiment is to reproduce the authors' D21
witness first, then compare D21--D27 at `k=44` on exactly the new legal support.

### The incomplete-range Type-IIc estimate

The draft's Lemma 5.6 retains the exact `r,u,d_1` intervals but still uses the
three terminal walls

```text
8 omega + 4 delta + 2 gamma < 1
32 omega + 10 delta - gamma < 0
48 omega + 16 delta - 4 gamma < -1.
```

The last wall is the active one in its Table 3 datum. The experimental theorem
assembled in this repository attacks that wall itself by preserving cancellation
in the incomplete rectangles. It is therefore not subsumed by the 212 draft.

It cannot yet be inserted verbatim: our hostile audit was performed at
`delta=7/250`, whereas the new paper uses `delta=41/2500` and a different exact
factor-extraction family. The argument must be re-derived and independently
reviewed in the notation and ranges of Lemma 5.6.

If that transfer succeeds with the same structural inequality, exact arithmetic
shows that the old terminal wall would stop being the bottleneck. Holding
`delta=41/2500` and `xi_2=2/5` fixed, the paper's `A_1=257/1000` could move only
to

```text
A_1 = 2059/8000 = 0.257375,
```

because Proposition 7.8(II)'s *first* branch would then bind before the improved
terminal wall. This is a gain of `3/8000=0.000375` in `A_1`, raising the generated
modulus level from `0.514` to `0.51475`. Replay the exact comparison with

```bash
python scripts/check_bgp212_headroom.py
```

This is structural headroom, not a variational score and not a theorem until the
Type-IIc transfer is proved.

### Multi-band support

The new proof permits multiple bands in its general definitions but chooses one
band. Our multi-band support representation and joint support/function optimizer
remain relevant. A natural next construction is to keep the new low-total-mass
region generous while using stricter rough caps only in higher-total-mass bands.
The authors' exact continuum-packing certificate supplies the correct arithmetic
feasibility model for such a search.

## Immediate target hierarchy

The tuple consequences are

```text
DHL[45,2] -> H_1 <= 212
DHL[44,2] -> H_1 <= 210
DHL[43,2] -> H_1 <= 200
```

Accordingly:

1. `k=44` / 210 is the immediate target. Search the new legal support with the
   released D21 witness, D27 projection, and the transferred Type-IIc wall.
2. `k=43` / 200 is the first materially stronger target. It probably needs a
   further arithmetic gain: movement of the next Type-II wall, a successful
   multi-band construction, or distribution matched to the actual sieve
   coefficients rather than uniform coverage of every generated modulus.
3. The old `k=48` computation should finish as an independent engine validation,
   but no further optimization should be based on the superseded support.

## Reproduction blockers

The preliminary PDF states the exact quotient but does not contain or link the
846-entry rational coefficient vector. Its formal appendix also says that the
exact rational variational calculation is supplied to the Lean development as a
hypothesis. At this status date, the public `AxiomMath/PrimeGapsLib` repository
still contains the 246 formalization rather than the new 212 sources.

The exact Table 3 datum, every Appendix B/Table 6 row, all five unsimplified
modulus classes, and the quantified Proposition 7.8(A)--(E) packing problem are
now encoded in `primegaps.bgp212`. This import found three stale-edit
discrepancies in the draft; see `reproduction/212/README.md`. All three leave
the symbolic statements or relevant inequality signs unchanged.

The remaining external dependency is the release of the rational vector,
packing certificate trees, and formal source. The vector can alternatively be
reconstructed with our numerical pipeline and exactified independently. The
packing trees must be obtained or independently regenerated before the
arithmetic half is a replay rather than an encoded obligation.

## Recommended next actions

- Reconstruct the D21 vector on the exact new support and replay its quotient.
- Add the displayed H45 tuple and verify cardinality, diameter, and admissibility.
- Obtain or reconstruct the D21 rational witness and reproduce its exact quotient.
- Run the stable projected search at `k=44`, first at D21 and then D22--D27.
- Re-audit the incomplete-rectangle Type-IIc argument at `delta=41/2500`.
- Jointly optimize `A_1`, `delta`, the full `B_m` staircase, and the variational
  function; do not sweep any of them in isolation.
- Exactify only the best fixed candidate with the boundary-only Arb certifier.

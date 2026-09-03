# Boundary-only Arb J certifier

The final degree-27 variational gate is evaluated as an exact
legal-minus-unrestricted correction. The unrestricted `I` upper bound and
unrestricted `J` value are closed-form exact rationals; only cells changed by
the `B_m` restrictions are sent through geometry.

For each affected cell, the implementation constructs legal and unrestricted
marginal slices over exact rational affine geometry, factors
`L_legal R_legal - L_full R_full` before multiplication, routes the correction
into target-signature polynomials, constructs every raw monomial geometry
moment on first use, and reuses those moments for all target contractions on
that cell. Exact inputs are converted into a Sage `RealBallField(128)`, and
outward-rounded rational endpoints are recorded.

Density-weighted moments are never rebuilt per signature pair. Target
densities are compiled for one `(large, shifted)` status at a time and reused
over all cells in that status. Unaffected cells are absent from the hot loop.
Cells with the same exact slice operators share the routed candidate
polynomial. The published one-band D21 support has 284 affected cells in 30
statuses. The current two-band D27 support has 1,616 affected cells in 32
statuses and 515 distinct slice regimes; the counts are intentionally
different.

## Independent calibration

The complete 284-cell D21 run took 981.83 seconds on the development machine
under concurrent load. An independent closed-form unrestricted calculation
took 14.59 seconds. The exact frozen legal-oracle correction is contained in
the Arb interval:

```text
exact normalized correction  -0.0017438933429449071
Arb normalized lower          -0.0017438937392852984
Arb normalized upper          -0.0017438929482970595
```

The normalized interval width is about `7.91e-10`. This validates the complete
cell partition, density normalization, correction sign, polynomial routing,
geometry moments, and outward-rounded accumulation against the frozen rational
D21 result. A small exact rational unit test independently checks the pipeline.

## Detached D27 calculation

Start or resume from the repository root:

```bash
scripts/run_d27_boundary_certificate.sh start
```

The default is four workers to limit memory and interference. Override it only
at launch with `BOUNDARY_J_WORKERS=N`. Monitor without attaching:

```bash
scripts/run_d27_boundary_certificate.sh status
```

The append-only cell checkpoint, bound manifest, log, and PID live under
`.research/work/failed-experiment-revival/failed-ranker-20260902/`. Restarting
verifies the manifest and skips completed cells. A benchmark made with
`--limit-cells` or `--target-limit` is ineligible for certification.

After the process exits, produce the exact gate result:

```bash
scripts/run_d27_boundary_certificate.sh finalize
```

Finalization checks unique cell indices and completeness, sums rational Arb
endpoints, normalizes them exactly, binds the existing exact unrestricted `J`
and exact unrestricted `I` upper bound by candidate hash, and reports either
`certified_strictly_above_one: true` or an inconclusive result with the next
suggested precision. Checkpoint replay never repeats geometry.

## Scope

This is a scalar verifier for one fixed rational candidate, not a replacement
for the floating-point discovery operator. Arb proves containment of the
implemented finite exact integral. The surrounding analytic implication still
depends on the separately recorded Type-IIc endpoint assumptions and must be
reported with those assumptions.

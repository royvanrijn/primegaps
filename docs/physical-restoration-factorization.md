# Signed physical restoration and the exact source oracle

## Result

Two previously missing finite experiments have now been made executable.

First, the frozen PrimeGaps186 cap model is assembled as explicit `77 x 77`
quadratic forms and optimized with the **signed** hybrid restoration

```text
H = J0 + (a+b) Jplus + b Jtail,
a = 2479900401/2500000000,
b = -843183/1000000000.
```

This is materially different from the ideal full-face form
`Jfull = J0 + Jplus + Jtail`. Quadratic extrapolation in inverse mesh, using
the `2048`, `4096`, and `8192` direct-convolution meshes and evaluating the
finest-mesh optimizer at every mesh, gives:

| dimension | signed-restoration score | deficit from 1 | ideal full-face control |
|---:|---:|---:|---:|
| 40 | 1.0001866542 | -0.0001866542 | 1.0256101858 |
| 39 | 0.9943709810 | 0.0056290190 | 1.0192845867 |
| 38 | 0.9883950587 | 0.0116049413 | 1.0127922261 |

The independently optimized mesh values give the close cross-checks
`1.0002083815`, `0.9943911126`, and `0.9884136585`. At `k=40`, the first of
these is `0.0000023021` above the published rigorous fixed-vector lower
endpoint `1.0002060794`, which is the expected ordering within the observed
mesh error.

Every one of the 97 source-cover costs in the full fixed-geometry certificate
is a positive-semidefinite quadratic form which is subtracted from `H`.
Consequently the displayed signed-restoration eigenvalue is an **upper
screen** for the full loss-charged operator. The experiment does not pretend
to assemble those 97 expensive matrices: it makes that calculation unnecessary
for the frozen `k=39/38` decision, because even the upper screen is below one.
Trial-only reoptimization at the current geometry is therefore a numerical
no-go. A viable `k=39` attempt must move the geometry, hybrid coefficients, or
source construction before paying for a full directed certificate.

This is a binary64 mesh extrapolation, not a rigorous eigenvalue enclosure.
The complete machine-readable record is
[`experiments/physical_restored_operator.json`](../experiments/physical_restored_operator.json).

## Exact factorization oracle

The second experiment reconstructs all 29 old and 43 new source-ladder rows
with `fractions.Fraction`. It then compares, after the largest-fragment and
opposite-root guards used by PrimeGaps186 Lemma 1.4:

```text
union of the exact order-three row failures
versus
the grouped nonlargest H_{5/2} failure.
```

The order-three inventory is:

| group | source rows | exact configurations classified | false negatives |
|---|---:|---:|---:|
| outer | 19 | 134,607 | 0 |
| old inner | 4 | 79,208 | 0 |
| new inner | 15 | 160,042 | 0 |
| **total** | **38** | **373,857** | **0** |

Each total combines a deterministic critical-value census through seven
active fragments with a seeded uniform cell stress census through fourteen
fragments. The classification is exact rational arithmetic; the selection of
configurations is finite and therefore not a proof or a probability statement
about the physical fragment law.

The oracle also records and exactly replays one strict-overcoverage witness in
each group. Thus grouped `H_{5/2}` is not equivalent to the union of exact row
predicates. This is expected: the source argument uses it as a safe majorant
after separate guards. The test found no violation of that direction in its
finite census, while the published lemma—not the census—supplies the general
implication.

Endpoint conventions are explicit. The first sorted occurrence alone is
excluded as the largest witness; a tied second occurrence is nonlargest, and
the obstruction uses the full inclusive prefix at that value. Small inactive
fragments may supply residual radial mass without entering either obstruction.

The record is
[`experiments/physical_factorization_oracle.json`](../experiments/physical_factorization_oracle.json),
and the reusable exact predicates are in
[`src/primegaps/physical.py`](../src/primegaps/physical.py).

## Reproduction

Cheap replay validates the saved inputs, fits, group parameters, census
classification, and strict-overcoverage witnesses without repeating either
calculation:

```bash
python scripts/check_physical_restored_operator.py
python scripts/check_physical_factorization_oracle.py
```

The opt-in computations are separate:

```bash
python experiments/physical_restored_operator.py
python experiments/physical_factorization_oracle.py
```

The restoration calculation intentionally refuses meshes above `8192`: the
current float model switches from direct convolution to FFT there, and exact
leading-zero regions make that path unsuitable for this extrapolation.

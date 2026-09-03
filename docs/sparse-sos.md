# Sparse PSD sum-of-squares sieve screen

The rank-one sieve ansatz can be relaxed from one coefficient vector
`Q=cc^T` to a finite sum of squares `Q=V^T V`.  For preassembled bilinear
forms this gives

```text
maximize    <M_objective,Q>
subject to  <M_I,Q> = 1
            Q >= 0
            Q_ab = 0 for forbidden support interactions.
```

`primegaps.sos` implements the sparse SDP, its explicit SOS factorization, and
the correct rank-one comparator.  A feasible rank-one vector can use only one
clique of the analytic compatibility graph, so the comparator enumerates
maximal cliques and solves a generalized eigenproblem on each.  The SDP is
solved in dual form with CVXOPT; its PSD dual variable is the requested `Q`.

The test suite includes a signed induced four-cycle for which the best
rank-one value is `1` and the rank-two PSD value is `sqrt(2)`.  This guards
against an implementation that always collapses to rank one.  It also records
the useful structural prescreen that chordal masks cannot help: the sparse PSD
cone of a chordal graph decomposes into clique-supported PSD matrices, and a
linear normalized objective then selects one clique.

## First sieve-derived screen

The first bank uses four individually self-certified `(A_j,m,B_jm)` cells.
The production distribution oracle gives the induced cycle

```text
0 -- 1
|    |
3 -- 2
```

with `0--2` and `1--3` forbidden, both because the implemented Type-IIc local
witness fails.  Each region carries the 28-dimensional symmetric D7 polynomial
space.  Its 15 strongest locally whitened modes are retained, giving 60
localized components.  `I` and `kJ` are accumulated by exact-large-count
stratified scrambled Sobol integration with `2^16` points per stratum.  The
same matrices are supplied to the rank-one and PSD optimizers.

| k | best rank-one | sparse PSD | relative PSD advantage | PSD rank |
|---:|---:|---:|---:|---:|
| 47 | 0.166578785463 | 0.166578785308 | -9.34e-10 | 1 |
| 46 | 0.168342211365 | 0.168342210516 | -5.04e-9 | 1 |

The tiny negative differences are solver feasibility/normalization error, not
a rank-one improvement.  Both returned PSD matrices have numerical rank one;
their forbidden entries are below `3.4e-17`.  Thus this deliberately
non-chordal 60-component D7 bank shows no PSD advantage, far below the proposed
`2.2%` signal threshold.

The absolute quotients are weak and remain numerical QMC values.  This result
kills only this first localized bank.  It is not evidence that every richer
support graph has zero advantage, and it is not an exact sieve certificate.
Because the optimizer already collapsed to rank one, no second-stage
modulus/residue cancellation audit was triggered.  Any future positive mask
screen must perform that coefficient-level audit before being interpreted as a
valid sieve construction.

## Rebuild and replay

The expensive build is explicit and requires the `sdp` optional dependencies
(or the repository's Sage environment):

```bash
PYTHONPATH=src sage -python scripts/run_sparse_sos_experiment.py \
  --k 47 46 --degree 7 --modes 15 --log2-samples 16 --seed 260903 \
  --output-directory .research/work/sparse-psd-sos/root-20260903/run \
  --summary .research/work/sparse-psd-sos/root-20260903/result.json
```

Cheap replay reads the recorded matrices and factors only.  It does not run
QMC or an optimizer:

```bash
python scripts/replay_sparse_sos_experiment.py \
  .research/work/sparse-psd-sos/root-20260903/result.json
```

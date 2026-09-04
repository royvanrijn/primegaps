# BGP212 reproduction baseline

This directory pins the exact public data from the preliminary 2026-09-03 draft
*A New Bound for Small Gaps Between Primes* (`H_1 <= 212`).  It is the baseline
for all new searches in this repository.

## Replayed now

- the complete rational parameter datum from Table 3;
- support consistency and hereditary cap conditions;
- the Harman/direct-prime scalar inequalities;
- the global Type I/II/III walls in Proposition 7.8;
- exact symbolic definitions of all five Section 5 modulus classes, including
  the unsimplified Type-IIc intervals for `r`, `u`, and `d_1`;
- the complete quantified continuum-packing problem in Proposition 7.8(A)--(E),
  including the two Type-I partitions and the `(gamma, omega_0)` rectangle in
  Type IIc;
- every printed row of Appendix B/Table 6, alongside an independent
  recomputation from Table 3 and Proposition 7.8;
- the displayed H45 tuple: cardinality 45, diameter 212, and admissibility.

Run:

```bash
python scripts/verify_bgp212_parameters.py
```

The verifier emits a machine-readable description of the modulus classes and
packing problem. Pass `--output result.json` to persist it atomically.

## Source audit

The exact replay exposes three apparent stale-edit issues in the preliminary
draft rather than silently normalizing them:

- Appendix B prints the first Type-II range wall as
  `69599997/2000000000`. Direct substitution in Proposition 7.8(II) gives
  `347999991/10000000000`, larger by `3/5000000000`. Both are positive, so
  this does not change whether the displayed parameter point passes.
- Section 9.2 writes `A_1=513/2000` and `delta=179/10000` inside one expansion,
  while Table 3 gives `257/1000` and `41/2500`. The two pairs have the same
  combination `3*A_1+delta=3937/5000`, so the subsequent dominant-wall slack
  is unchanged.
- The proof of Lemma 9.1 substitutes the old rescaled caps
  `(31/50, 31/50, 17/25, ...)`; Table 3 instead starts
  `(777/1250, 397/625, 7/10, ...)`. The lemma statement itself remains symbolic
  in `C`, so this is stale explanatory prose rather than a changed statement.

These are transcription-level findings about the public draft, not objections
to the sign of the relevant analytic margins.

The support is also provided in the numerical-builder format as `support.json`.
Its cap list is expanded through `floor(1/delta)=60`; the paper's cap is constant
from `m=11` onward.

## Not yet independently replayed

The preliminary PDF reports the exact degree-21 quotient

```text
J_T(F*) / I_T(F*) = 4.00438409833460131937... > 4
```

but does not include the 846 rational coefficients of `F*`. It likewise reports
455 exact continuum-packing root certificates without publishing their
certificate trees. We now encode the exact statement and independently derive
the expected count `91 * 5 = 455`, but that count is not a replay of the trees.
The public PrimeGapsLib revision available at the status date still contains the
246 formalization. These are external inputs we must obtain or reconstruct.

A full reproduction has three independent gates:

1. replay the 455 arithmetic packing certificates and the half-level transition;
2. obtain/reconstruct `F*` and reproduce the exact rational quotient;
3. connect those two facts through the support-specific sieve criterion to
   `DHL[45,2]`, then verify the finite H45 endgame.

The old D27 `k=48` Arb job may finish as an engine validation, but it is not on
this reproduction path and must not be used as the current record baseline.

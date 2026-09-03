# BGP212 reproduction baseline

This directory pins the exact public data from the preliminary 2026-09-03 draft
*A New Bound for Small Gaps Between Primes* (`H_1 <= 212`).  It is the baseline
for all new searches in this repository.

## Replayed now

- the complete rational parameter datum from Table 3;
- support consistency and hereditary cap conditions;
- the Harman/direct-prime scalar inequalities;
- the global Type I/II/III walls in Proposition 7.8;
- the displayed H45 tuple: cardinality 45, diameter 212, and admissibility.

Run:

```bash
python scripts/verify_bgp212_parameters.py
```

The support is also provided in the numerical-builder format as `support.json`.
Its cap list is expanded through `floor(1/delta)=60`; the paper's cap is constant
from `m=11` onward.

## Not yet independently replayed

The preliminary PDF reports the exact degree-21 quotient

```text
J_T(F*) / I_T(F*) = 4.00438409833460131937... > 4
```

but does not include the 846 rational coefficients of `F*`.  It likewise reports
455 exact continuum-packing root certificates without publishing their certificate
tree.  The public PrimeGapsLib revision available at the status date still contains
the 246 formalization.  These are external inputs we must obtain or reconstruct.

A full reproduction has three independent gates:

1. replay the 455 arithmetic packing certificates and the half-level transition;
2. obtain/reconstruct `F*` and reproduce the exact rational quotient;
3. connect those two facts through the support-specific sieve criterion to
   `DHL[45,2]`, then verify the finite H45 endgame.

The old D27 `k=48` Arb job may finish as an engine validation, but it is not on
this reproduction path and must not be used as the current record baseline.

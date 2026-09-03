# Gap 236 formalization crosswalk

This note maps the active `k=48` research to the public formal proof architecture in
[AxiomMath/PrimeGapsLib](https://github.com/AxiomMath/PrimeGapsLib), pinned here at
commit `1faa7b14e82ddebc2772dfb9153922f01b106477`.

The interactive page at <https://primegaps.cicada71.net/> is useful as a navigation /
proof-blueprint aid, but it is not a pinned proof dependency in this repository.  We
prefer the primary source files and theorem names below so the crosswalk is auditable.

## The useful split

PrimeGapsLib's `Gap246` proof already separates the two logically different jobs we
need:

1. **analytic/sieve statement (`DHL[k,2]`)**: every admissible `k`-tuple has infinitely
   many translates containing at least two primes;
2. **finite tuple/endgame**: instantiate that statement with an explicit admissible
   tuple and convert two primes in a bounded interval into a bound on consecutive
   prime gaps.

At the pinned revision, the relevant primary files are:

- `PrimeGapsTheory/Gap246/Endgame/Main.lean`
- `PrimeGapsTheory/Gap246/Endgame/Witness.lean`
- `PrimeGapsTheory/Gap246/Tuple/H50.lean`
- `PrimeGapsTheory/Gap246/Sieve/Certificate/Hypothesis.lean`

`Gaps246.thm_dhl` has essentially the theorem shape we want at `k=48`:

```text
for every strictly increasing admissible h : Fin k -> Nat,
{n | at least two of n + h_i are prime} is infinite.
```

`Gaps246.thm_main` then uses only the explicit tuple, its admissibility/diameter,
and generic endgame machinery to turn that statement into `H_1 <= 246`.

## Our exact target interface

For `k=48`, isolate the following proposition:

```text
DHL48 :=
  for every strictly increasing admissible h : Fin 48 -> Nat,
  infinitely many translates n+h contain at least two primes.
```

Then the desired result decomposes as

```text
new Type-IIc / distribution theorem
        +
fixed support + fixed rational variational witness
        |
        v
      DHL48
        |
        + explicit admissible H48, diameter 236
        v
      H_1 <= 236.
```

The last arrow is not new number theory.  It is the same finite/endgame argument as
PrimeGapsLib's `Gaps246.thm_main`, with `48` and `236` substituted for `50` and
`246`.  A source-level adapter is tracked at `formalization/Gap236Endgame.lean`.

## Explicit H48 witness

The repository now uses the following independently generated representative:

```text
{0, 6, 8, 14, 18, 24, 26, 48, 50, 54, 56, 60,
 66, 68, 74, 78, 80, 84, 90, 96, 98, 104, 110, 116,
 120, 126, 134, 138, 144, 150, 158, 164, 168, 176, 180, 186,
 188, 194, 200, 204, 206, 210, 216, 224, 228, 230, 234, 236}
```

It has 48 elements, minimum 0, maximum 236, and omits at least one residue class
modulo every prime `p <= 47`; primes larger than 48 are automatic.  Run

```bash
python scripts/verify_h48_tuple.py
```

for the exact finite replay.  OEIS A008407 lists the minimum permitted diameter at
`k=48` as 236, and the MIT narrow-admissible-tuples database is the maintained
computational source.  We only need existence of an admissible diameter-236 tuple for
the gap theorem; the representative need not be identical to another database's
representative.

## What reuses PrimeGapsLib directly

The following pieces are already generic enough in PrimeGapsLib to reuse in a future
formalization:

- `Finset.Admissible` and finite tuple verification;
- `orderEmbOfFin` / strict ordering of a finite tuple;
- `ContainsAtLeastPrimes`;
- the conversion from infinitely many translates containing two primes to infinitely
  many consecutive prime gaps bounded by the tuple diameter;
- the final `Filter.atTop` statement used by `Gaps246.thm_main`.

Thus a formal `DHL48` theorem plus the explicit H48 is enough to make the last step
look almost exactly like the existing `246` theorem.

## What does *not* plug in unchanged

It would be incorrect to claim that the current result is merely an
`ExistsEpsCert 48` instance.

`PrimeGaps.ExistsEpsCert` and `Gaps246.prop_witness` are tailored to the Polymath8b
enlarged-simplex/Bombieri--Vinogradov setup.  Our candidate uses Stadlmann's shaped
support and, at the safe endpoint, the new incomplete-rectangle Type-IIc input.
Therefore the missing formal bridge is the **support/distribution -> DHL48** theorem,
not the tuple/endgame.

A clean future PrimeGapsLib-style module should have a new theorem interface such as

```text
Stadlmann48Certificate -> Stadlmann48Distribution -> DHL48
```

where `Stadlmann48Certificate` binds the exact rational candidate/support and the
rigorous `J`/`I` inequality, and `Stadlmann48Distribution` contains the precise
analytic hypotheses (including the reviewed Type-IIc replacement).

## Current proof-obligation ledger

| Layer | Gap246 analogue | Gap236 status |
|---|---|---|
| finite tuple | `H50`, `card_H50`, `diameter_H50`, `admissible_H50` | explicit H48; exact Python replay; Lean source adapter added |
| variational witness | `ExistsEpsCert 50` | D27 prefix-four candidate; unrestricted pieces exact; Arb boundary certificate currently running |
| distribution theorem | Bombieri--Vinogradov + split-support theorem | Stadlmann Proposition 2/3 plus checked incomplete-rectangle Type-IIc improvement; full human/typeset review still required |
| positivity / DHL | `Gaps246.prop_witness` + `Gaps246.thm_dhl` | mathematical bridge still to be written/formalized for shaped support |
| final gap endgame | `Gaps246.thm_main` | same argument after `DHL48`; source adapter tracked |

## Promotion rule

Do not state a proved unconditional `H_1 <= 236` merely because the Lean endgame
adapter is straightforward.  Promote the final result only after both independent
gates pass:

1. the fixed D27 rational candidate has a rigorous legal-support variational
   certificate (the current Arb calculation is designed to provide this); and
2. the new Type-IIc theorem is fully reviewed with all hypotheses matching the
   support used by that certificate.

Once those hold, expressing the result in the PrimeGapsLib/Cicada proof framing is a
strong presentation choice: the genuinely new work is isolated to one analytic DHL48
bridge, while the finite admissible-tuple and prime-gap endgame can be checked using
the same architecture already used for 246.

## External references

- Julia Stadlmann, *Bounded gaps between primes*, arXiv:2608.31126.
- AxiomMath/PrimeGapsLib, pinned commit
  `1faa7b14e82ddebc2772dfb9153922f01b106477`.
- A. V. Sutherland / MIT narrow admissible tuples: <https://math.mit.edu/~primegaps/>.
- OEIS A008407: <https://oeis.org/A008407> (`a(48)=236`).
- D. H. J. Polymath, *Variants of the Selberg sieve, and bounded intervals containing
  many primes*, Research in the Mathematical Sciences 1:12 (2014).

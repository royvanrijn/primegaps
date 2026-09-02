# PrimeGapsLib cross-pollination

This repository ports the reusable certificate machinery from Axiom Math's
[PrimeGapsLib](https://github.com/AxiomMath/PrimeGapsLib), pinned for this audit
at commit `1faa7b14e82ddebc2772dfb9153922f01b106477`.  The port is Python and exact
rational/integer arithmetic; it does not introduce Lean as a build dependency.

## What the external certificate represents

The sole source certificate
`PrimeGapsCert/Gap246/k50e25d25n1295.json` is a list of 1,295 triples

```text
[a, alpha, coefficient]
```

representing the sparse symmetric polynomial

```text
sum coefficient * (1 + epsilon - sum(t))^a * m_alpha(t).
```

Here `alpha` is the zero-free multiset of exponents of a monomial-symmetric
polynomial.  This avoids expanding an orbit in 50 variables.  Only 272
signatures are required after closing the source signatures under one-part
erasure.

For signatures `alpha` and `beta`, PrimeGapsLib reduces every symmetric simplex
pairing to the factorial moment

```text
F_k(alpha,beta) = sum_A sum_B product_i (A_i + B_i)!,
```

where `A` and `B` range over the distinct exponent-vector embeddings of their
signatures in `k` coordinates.  The exact mass pairing is

```text
R^(k+d) * (a1+a2)! / (k+d)! * F_k(alpha1,alpha2),
d = a1+a2+sum(alpha1)+sum(alpha2).
```

Integrating one coordinate first gives the exact marginal pairing as a sum over
the distinct erased parts `r1` and `r2`.  Its factors separate into:

- two beta factors `r! a! / (a+r+1)!`;
- a radial factor depending only on two total degrees;
- `F_(k-1)(erase(alpha1,r1), erase(alpha2,r2))`.

These formulas are implemented by `enlarged_simplex_pairing` and
`simplex_marginal_pairing` in `primegaps.symmetric`.

## Factorial-moment algorithms transferred

Two complementary algorithms are included.

1. `factorial_moment` groups equal exponents into multiplicity classes and
   enumerates overlap-count profiles.  It estimates both orientations and puts
   the smaller branching profile on the recursive side.
2. `FactorialMomentTable` writes the coordinate-splitting transition as
   `I + N`.  Erasing one part from one or both signatures makes `N` nilpotent.
   Its short ladder is computed once, and moments in any requested dimension
   are reconstructed as `sum_r choose(k,r) N^r`.  A level can be nonzero only
   when `r <= part_count <= 2r`, so inactive signature pairs are skipped.

The second method is the matrix-assembly path.  On the 246 source it computes
both dimensions 49 and 50 over the 272-signature closure in a fraction of a
second on the recorded environment.

## Sparse contraction transferred

`evaluate_sparse_symmetric_certificate` reproduces the structure of the packed
PrimeGapsLib checker while leaving Lean-specific numeral packing out of the
runtime representation:

- use one globally cleared positive denominator, so the final comparison is
  integer-only;
- group mass terms by `(slack degree, signature degree)`;
- aggregate marginal terms after each possible signature erasure into
  `(residual degree, radial degree, erased signature)` features;
- cache group-by-signature moment transforms;
- contract the smaller group against the cached transform of the larger group;
- evaluate only the diagonal and strict upper triangle, doubling the latter;
- cache the small scalar, beta-factor, and radial-factor tables by their degree
  arguments.

For Axiom Math's source certificate this reduces 1,295 squared pair terms to:

| layer | count |
|---|---:|
| erase-closed signatures | 272 |
| mass degree groups | 138 |
| mass transforms actually queried | 5,533 |
| nonzero marginal features | 1,504 |
| marginal degree groups | 172 |
| marginal transforms actually queried | 6,057 |

The local replay obtains an exact quotient of approximately
`4.000000917784551`, with a strictly positive integer difference
`50 * marginal - 4 * mass`.  This independently checks the numerical content
of the source JSON through the published closed formulas; it does not replace
the Lean proof of those formulas or prove the newer bound 240.

## Rebuild

Clone the exact external revision and run the tracked replay command explicitly:

```bash
git clone https://github.com/AxiomMath/PrimeGapsLib.git .research/PrimeGapsLib
git -C .research/PrimeGapsLib checkout 1faa7b14e82ddebc2772dfb9153922f01b106477
python scripts/verify_primegapslib_certificate.py \
  .research/PrimeGapsLib/PrimeGapsCert/Gap246/k50e25d25n1295.json
```

The default command is a replay of recorded source coefficients, not an
eigensolve and not a regeneration of the 1,295 coefficients.

## Relevance to the k=49 Stadlmann assembler

The old certificate uses a concentric enlarged/shrunken simplex, whereas
Stadlmann's support has small/large status cuts and `B_m` boundaries.  Its
radial scalar table therefore cannot be copied verbatim.  The reusable boundary
is nevertheless substantial:

- store each basis orbit as a canonical exponent signature, never as `k!`
  permutations;
- close only the signatures reached by the C/D and marginal erasure operators;
- precompute factorial moments for all required effective dimensions in one
  nilpotent ladder;
- factor matrix entries into degree-only rational kernels and signature-only
  moments;
- group assembly work by degree/status descriptors and cache transforms;
- assemble upper triangles directly into packed exact storage;
- keep expensive generation separate from a small hash-and-arithmetic replay.

That integration is complete for the fixed degree-21 Stadlmann candidates: the
frozen 240 reproducer replaces the full-simplex scalar by the verified C/D
status kernel, and the accelerated backend preserves the signature DAG and
grouped contractions while adding candidate-independent moment caches. The
remaining integration problem is a numerically stable full-support matrix
constructor for discovering higher-degree vectors; it must not fall back to the
refuted monomial-basis floating-point replay.

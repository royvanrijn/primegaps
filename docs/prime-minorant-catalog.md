# Prime-minorant catalog

This catalog covers only the Harman/Buchstab decomposition. It does not choose
the GPY test function, optimize `F`, or change the support geometry. Generate
the full optimizer-facing JSON with

```bash
python scripts/build_minorant_catalog.py --order 160 --comparison-order 96 \
  > prime-minorant-catalog.json
```

The JSON records retained mass, the pointwise negative bound `c2`, every base
Type I/II/III convolution regime, any additional exceptional regime that must
be distributed, and the modulus family over which distribution is required.

## Stadlmann's parameterized family

Write `t=xi2`. Proposition 2 of Stadlmann,
*Bounded gaps between primes*, arXiv:2608.31126v1, gives a minorant provided

```text
2 xi1 + 3 t < 2,       xi1 + 9 t < 4,
2 xi1 + t > 1,         17 t < 7,          t <= xi3.
```

For fixed `t`, the mass depends only on `t`. The admissible interval is

```text
(1-t)/2 < xi1 < min(1-3t/2, 4-9t),       t < 7/17.
```

Increasing `xi1` shrinks the required Type I class without changing mass.
Increasing `xi3` enlarges the Type III class without changing mass. Hence the
regime-demand frontier is

```text
xi3 = t,       xi1 -> min(1-3t/2,4-9t) from below.
```

The required base sequence classes are:

| class | convolution structure | exponent conditions |
|---|---|---|
| I | `alpha * beta_smooth` | smooth factor `N=x^gamma`, `gamma >= xi1` |
| II | `alpha_SW * beta_SW` | `t <= gamma <= 1-t` |
| III | `alpha * psi1_smooth * psi2_smooth * psi3_smooth` | `1-2xi3 <= gamma_i <= xi3`, `gamma_i+gamma_j >= 1-xi3` |

Every class must equidistribute over the support-supplied `Q*`: squarefree
moduli

```text
q = e e' product_i(f_i) product_i(f'_i),
```

where `ee'` is `x^delta`-smooth and the remaining factor sizes satisfy the
fixed support's `A/B` bounds. The minorant also imposes the compatibility

```text
max_j B[j,1] < beta = 1-2t.
```

This is an input constraint for the global checker, not a choice made here.

At `t<=0.4` both exceptional regions below are empty. In particular the
published `(xi1,xi2,xi3)=(0.38,0.4,0.4)` candidate is the prime indicator
itself, with retained mass `1` and `c2=0`.

## Four discard/retention choices

For `0.4<t<7/17`, the Buchstab identity isolates two nonnegative exceptional
five-prime pieces:

```text
A: beta<a4<a3<a2<a1<t, a1+a2<t,
   a2+a3+a4>1-t, a5>a4.

B: beta<a2,a3,a4,a5,a6<8t-3, a2>a3, a4>a3,
   a2+a4<t, a3+a4+a5>1-t, a6>a5.
```

Let their prime-mass integrals be `L_A(t)` and `L_B(t)`. Their pointwise
multiplicities are at most `4` and `20`. The published Proposition 2 minorant
discards both. The Buchstab identity also gives three conditional alternatives:

| candidate | retained mass | `c2` | extra distribution required over the same moduli |
|---|---:|---:|---|
| discard A and B | `1-L_A-L_B` | 24 | none |
| discard A, retain B | `1-L_A` | 4 | B |
| retain A, discard B | `1-L_B` | 20 | A |
| retain A and B | `1` | 0 | A and B (equivalently direct primes) |

The last three rows are conditional derived candidates, not claims that the
present distribution theorems cover the extra sequence. They let the global
optimizer price a precise theorem upgrade against mass and `c2`.

Selected baseline values (discard both) are:

| `t=xi2` | estimated retained mass | `beta=1-2t` |
|---:|---:|---:|
| 0.4005 | 0.999999998912 | 0.1990 |
| 0.4025 | 0.999999319764 | 0.1950 |
| 0.40481 | 0.999990666429 | 0.19038 |
| 0.4075 | 0.999944686940 | 0.1850 |
| 0.4090 | 0.999885082182 | 0.1820 |
| 0.4105 | 0.999786615598 | 0.1790 |
| 0.4115 | 0.999692431862 | 0.1770 |

At the literature anchor `t=0.40481`, the estimates are
`L_A=2.417331852e-7` and `L_B=9.091838091e-6`. Thus retaining B while
discarding A would improve the mass to `0.999999758267` and reduce `c2` from
`24` to `4`, but requires equidistribution of exceptional regime B over `Q*`.

## Baker--Irving alternative

Baker and Irving, *Bounded intervals containing many primes*,
arXiv:1505.01815v1, give a distinct `eta`-parameterized decomposition for

```text
0 < eta < 22/3295.
```

Its equivalent regime parameters and smooth-modulus exponent are

```text
xi1 = 199/600 + 119 eta/240,
xi2 = xi3 = 2/5 + eta,
beta = 1/5 - 2 eta,
theta = 1/2 + 7/300 + 17 eta/120.
```

Here distribution is proved only for squarefree `x^delta`-smooth
`q <= x^(theta-epsilon)`, not for the broader `Q*` family. If `I(E(eta))` is
the Lemma 2 polytope integral, the two exceptional multiplicities contribute
`I` and `5I`. The same four discard choices therefore have losses `6I`, `I`,
`5I`, and `0`, and `c2=24,4,20,0` respectively.

Selected published-baseline estimates (discard both) are:

| `eta` | `theta` | estimated retained mass |
|---:|---:|---:|
| 0.00100 | 0.523475000 | 0.999999997291 |
| 0.00300 | 0.523758333 | 0.999999780563 |
| 0.00481 | 0.524014750 | 0.999998549601 |
| 0.00600 | 0.524183333 | 0.999996487737 |
| 0.00667 | 0.524278250 | 0.999994635431 |

This older family is not dominated merely because its modulus family is
narrower: its second exceptional loss is much smaller. The global checker must
compare the actual theorem coverage, not merge it with the Stadlmann `Q*`
family.

## Numerical status

The parameter inequalities, convolution polytopes, mass formulas, and
pointwise multiplicities are analytic. The displayed integral values are
deterministic numerical estimates, not rigorous enclosures. The generator
evaluates the innermost integral in closed form and applies nested
Gauss--Legendre quadrature to the remaining three dimensions. It records the
change between orders 96 and 160 for every candidate. A proof-bearing optimizer
must replace or enclose these values with certified interval integration before
using a close mass comparison as a theorem.

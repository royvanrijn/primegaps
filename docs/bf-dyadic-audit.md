# Dyadic audit of the proposed physical parity bilinear axiom

## Verdict

The proposed `B_F` is not a difficult unproved lemma. In the stated
quarter-rough specialization it is false in every nonempty active
Friedlander--Iwaniec block. The obstruction occurs before any estimate for
the physical lcm moduli is invoked.

Let `X<m<=2X`, `z=X^beta`, `beta=0.250001`, and

```text
a_i(m) = 1_{P^-(m)>z} w_F(m-h_i).
```

The actual Friedlander--Iwaniec axiom requires, for every

```text
sqrt(D)/Delta < V < sqrt(X)/delta,
1 <= C <= X/D,
```

the power-logarithmic bound

```text
sum_u |sum_{V<v<=2V, X<uv<=2X}
          gamma(v,C) mu(uv) a_i(uv)|
    << A_i(X) (log X)^(-222).
```

It is not the weaker statement with an unspecified `epsilon_X -> 0`.

Take the required value `C=1`. Then `gamma(v,1)=1`. Rough support makes
every contributing `v` z-rough. The upper endpoint is less than `z^2`, since
`beta>1/4` and `V<sqrt(X)/delta`; hence every contributing `v` is prime. If
`mu(uv)` is nonzero, then

```text
mu(uv) = mu(u) mu(v) = -mu(u).
```

For fixed `u` this sign is constant throughout the inner sum, while `a_i` is
nonnegative. Thus the absolute value removes the sign exactly:

```text
|sum_v mu(uv) a_i(uv)| = sum_v mu(uv)^2 a_i(uv).
```

There is no Mobius cancellation left to estimate.

## Every block

Put

```text
Y       = sqrt(D)/Delta,
Z       = sqrt(X)/delta,
V_j     = 2^j Y,                         0 <= j < J, V_j < Z,
U_j     asy X/V_j,
Q_F     = X^(2742997/5000000) = X^0.5485994,
Q_ell   = (2^ell,2^(ell+1)] cap [1,Q_F].
```

The product cutoff gives the exact outer interval

```text
X/(2V_j) < u <= 2X/V_j,
```

which is at most three ordinary dyadic `U` intervals. The potential Cartesian
family is `(j,ell)`, including the singleton `q=1`. No values of `D`,
`Delta`, or `delta` were supplied in the proposal, so this indexed family is
the only canonical meaning of "every dyadic block". The replay script can
expand every row after those parameters are declared.

For every full, `q`-recombined `V_j` block the requested data are:

| support in `(V_j,2V_j]` | `U` | recombined `q`-range | coefficient norm | required saving | parent-block classification |
|---|---|---|---|---|---|
| `2V_j <= z` | `X/(2V_j)<u<=2X/V_j` | `1<=q<=Q_F` | zero on the active support | none | already controlled (empty) |
| `V_j<z<2V_j<z^2` | same | `1<=q<=Q_F` | at `C=1`, `||gamma||_infty=1` and `||gamma||_2^2=#` admissible primes | `A_i(log X)^-222` | impossible when nonempty |
| `z<=V_j<Z` | same | `1<=q<=Q_F` | same | `A_i(log X)^-222`; natural size about `A_i/log X` in the viability model | impossible |
| `V_j>=z^2` | same | `1<=q<=Q_F` | not reached asymptotically because `Z<z^2` | n/a | outside current ranges |

The physical expansion supplies only the envelope for `q`. If

```text
c_i(q,a) = sum_{d,e producing (q,a)} Lambda_d Lambda_e,
```

then a dispersion or trace estimate needs a declared norm such as

```text
C_i,2(Q)^2 = sum_{Q<q<=2Q} sum_(a mod q) |c_i(q,a)|^2.
```

That norm is not defined or computable from the repository's 77 limiting
polynomial coefficients: the discrete `Lambda_d` normalization, presieving,
and residue-coloured aggregation have not been specified. Splitting in `q`
and applying the triangle inequality is not harmless, because cancellation
between the signed lcm coefficients is what reconstructs the nonnegative
square `w_F`. Consequently the generated `(j,ell)` rows record the exact
`q` envelope and the missing norm, but do not pretend that an individual
`Q_ell` slice inherits the impossibility proof. That proof applies to each
full, `q`-recombined `V_j` block.

The full-face envelope itself is exact from the recorded rational parameters:

```text
q <= X^(2 rho* T) = X^0.5485994,
p|q => p <= X^0.1903699274... < X^beta.
```

Thus `(uv,q)=1` on rough support, but this does not restore a Mobius sign.

## Classification

- **Already controlled:** only empty rough-support blocks and the irrelevant
  choices `C>=2V`, for which `gamma(v,C)=0` on prime `v`. Axiom `(B)` also
  requires `C=1`, so the latter does not help.
- **Potentially reachable by dispersion or trace functions:** none of the
  proposed `B_F` blocks. Analytic cancellation cannot prove a quantity is
  small after it has become an exact positive mass.
- **Genuinely new parity blocks:** these belong to a replacement for `B_F`,
  not to `B_F` itself. One must retain cancellation between the prime,
  semiprime, and triprime sectors (or across `u`) before taking an absolute
  value, while preserving the residue-coloured lcm coefficients.
- **Impossible under current ranges:** every nonempty active block at `C=1`.
  Under the same tensor-product rough-factor model used for numerical
  viability, an interior dyadic `V` block has size `asymp A_i/log X`, not
  `O(A_i log^-222 X)`. Summing the `O(log X)` blocks recovers a positive
  proportion of the semiprime/triprime mass.

## Two further applicability failures

Even apart from the sign obstruction, the 1998 theorem cannot currently be
applied to this `a_i`.

1. Its density axiom has `sum_{p<=y} g(p)=log log y+c+o(1)`. For the
   already `z`-rough sequence, divisibility by every prime `p<z` has density
   zero.
2. It requires a remainder axiom to a level `D>X^(2/3)`. No such `R_F`
   estimate is stated. The pair-lcm support bound `q<=X^0.5485994` is not a
   replacement for that axiom.

Therefore the direct signed target `(P)` may remain meaningful, but the
displayed Friedlander--Iwaniec route through the already quarter-rough
sequence must be withdrawn. A viable replacement is an architectural change,
not a blockwise strengthening of known dispersion or trace-function bounds.

## Sources checked

- John Friedlander and Henryk Iwaniec, *Asymptotic sieve for primes*, Annals
  of Mathematics 148 (1998), equations `(R)`, `(B)`, `(B1)--(B3)`, Theorem 1,
  and Sections 7--10:
  https://annals.math.princeton.edu/articles/13036
- James Maynard, *Primes in arithmetic progressions to large moduli II:
  Well-factorable estimates*, definition of well-factorability and Theorem
  1.1: https://arxiv.org/abs/2006.07088
- Jared Duker Lichtman, *Primes in arithmetic progressions to large moduli,
  and Goldbach beyond the square-root barrier*, triply well-factorable
  hypotheses and residue uniformity: https://arxiv.org/abs/2309.08522
- Existing repository records `ef778a1254...` and `417b502a...`, which check
  that the physical lcm coefficients carry CRT residue colours and do not
  currently satisfy a published fixed-residue transfer theorem.

## Replay

```bash
python scripts/check_bf_dyadic_audit.py
python scripts/check_bf_dyadic_audit.py \
  --log2-x 64 --d-exponent 3/4 --capital-delta-exponent 1/8 \
  --log2-little-delta 4 --expand
```

The second command is an explicitly labelled example. It is not a choice of
the missing asymptotic-sieve parameters for the physical sequence.

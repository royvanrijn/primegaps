# Audit of Shi's claimed `H_1 <= 234`

## Status

**Rejected; do not use as an analytic input.** This audit concerns the revised
12-page version uploaded to ResearchGate on 30 July 2025. The paper does not
establish its weighted distribution theorem, its `k=48` corollary, or its
claimed diameter 234.

The conclusion rests on several independent failures. It does not depend on
whether the preprint has been peer reviewed or cited.

## Fatal analytic failures

### Lemma 4.3 is false over its stated range

The lemma counts zeros of primitive Dirichlet L-functions with moduli
`q=s_1s_2`, where the two factors are distinct primes near `S` and `Q` is near
`S^2`. There are asymptotically `Q/(log Q)^2` such moduli. Each has
`(s_1-2)(s_2-2)` primitive characters, so the family contains
`asymp Q^2/(log Q)^2` primitive characters.

For a primitive character of conductor near `Q`, Riemann--von Mangoldt zero
counting gives

```text
N(T,chi) = (T/pi) log(qT/(2*pi*e)) + O(log(q(T+2))).
```

For growing polylogarithmic `T`, functional-equation symmetry therefore gives

```text
N*(Q,1/2,T) >> Q^2 T/log Q.
```

Shi's Lemma 4.3, which states no restriction excluding `sigma=1/2`, instead
gives

```text
N*(Q,1/2,T) << Q^2 T exp(-c sqrt(log Q)).
```

The ratio of the lower scale to the asserted upper scale contains
`exp(c sqrt(L))/L`, with `L=log Q`, and tends to infinity for every fixed
`c>0`. Thus the lemma is false as written. Restricting it to `sigma` near one
would remove this literal contradiction but would still require a new proof.

### The displayed reductions are not valid proofs

- Equation (4.1) contains `mn<=x`, but its character transform is subsequently
  factored into independent `m`- and `n`-sums. This includes terms with `mn>x`
  and is false unless additional rectangular support or a Mellin separation is
  supplied. Neither appears.
- The fourth-moment variance is called equivalent to a new third-moment form
  without an inequality, coefficient norms, or parameter ranges.
- Failure of the asserted off-diagonal power saving is said to produce a
  positive-density family of large rough character sums. No large-values lemma
  rules out concentration in a few terms.
- Lemma 4.4 replaces a large unsmoothed sum by one arbitrary smoothing and says
  it remains large. Cancellation makes that implication invalid without a
  partition-of-unity argument. Its claimed `delta(C)>0` becomes
  `delta >> 1/(C log q)` in the last line and hence is not fixed.
- The finite Euler product called a short Dirichlet polynomial has
  `2^(pi(z))` terms and maximum index `prod_{p<=z}p=exp(vartheta(z))`. No
  truncation or fundamental-lemma error estimate is supplied.
- Appendix B's first zero-density inequality does not normalize or construct
  its auxiliary polynomial. Scaling that polynomial scales the right side but
  not the zero count. In the next display the factor `L(s,chi)` disappears.
  The calculation that follows derives only the standard large-sieve scale;
  the decisive exponential saving is asserted in prose.

The main theorem itself specifies only that `R,S` are chosen appropriately. It
does not give their individual ranges, coefficient hypotheses, uniformity, or
the inequalities from which `1/400` is obtained.

## No transfer to the Maynard sieve

The theorem weights a maximum absolute progression error by a nonnegative
sequence supported on one rough-times-prime modulus class. The Maynard sieve
requires control of moduli assembled from least common multiples `[d_i,e_i]`
throughout its multidimensional support. The preprint gives no domination,
identity, decomposition, or restricted-support theorem connecting those
coefficients.

The name "well-factorable" does not supply this. Standard well-factorability
requires a bounded convolution for every factorization of the level; Shi
displays only one fixed convolution. Published applications emphasize that
even related sieve weights require a detailed factorization theorem.

## `k=48`: conditional arithmetic versus missing evidence

For two primes the ordinary Maynard condition is

```text
M_k > 2/theta.
```

At `theta=1/2+1/400=201/400`, the threshold is
`800/201=3.980099502487562...`. Consequently, a rigorous ordinary level of
distribution at this theta plus a rigorous `M_48>=3.998` certificate would
indeed make the scalar comparison work. Neither premise is established.

Polymath8b publishes `M_54>4.00238` on the standard simplex and
`M_(50,1/25)>4.0043` for its enlarged support. Its concluding discussion says
that `k=49` might be reached with further computation. It does not contain the
preprint's attributed rigorous `k=48,3.998` calculation.

The preprint also gives incompatible normalizations: Section 5 asks for a
theta-dependent `M_k>2`, while Appendix A asks for `rho_k>4`. Its displayed
definition of `S_k` contains no theta, but the next lines introduce an unknown
`O(Delta theta)` variation and then incorrectly demote its first-order
contribution to `O(Delta theta^2)`.

## Diameter 234 is independently impossible for `k=48`

Admissibility modulo 2 forces all tuple entries to have one parity. After
translation, it is enough to select even numbers from `[0,D]`. For every odd
prime `p<=47`, the mixed-integer model in
[`scripts/audit_shi_tuple.py`](../scripts/audit_shi_tuple.py) selects at least
one omitted residue and forbids selected integers from it. Primes greater than
48 are automatic.

Using NumPy 1.26.4, SciPy 1.16.1, and its bundled HiGHS 1.8.0 solver with zero
MIP gap:

```text
python scripts/audit_shi_tuple.py --diameter 234
# maximum_cardinality: 47; dual bound: 47; 4244 nodes

python scripts/audit_shi_tuple.py --diameter 236
# maximum_cardinality: 48; dual bound: 48; 10462 nodes
```

The diameter-236 solution is exactly the tuple supplied by the maintained MIT
Prime Gaps database. The preprint's abbreviated diameter-234 list therefore
cannot be completed to 48 admissible entries.

## Salvage assessment

The rough-times-prime convolution and character-orthogonality rewrite are
standard dispersion-method patterns. Weighted levels beyond one half already
reach `4/7` in Bombieri--Friedlander--Iwaniec, `3/5` in Maynard, and
`66/107 approximately 0.617` in Lichtman for appropriately factorable weights.
Those theorems are not interchangeable with a maximum-absolute-error theorem,
but they show that the bare idea of a weighted level beyond one half is not new.

The only potentially new ingredients would be (1) the asserted
maximum-absolute-error estimate for this sparse family and (2) a coefficient
theorem matching it to Maynard's support. These are exactly the missing parts.
There is therefore no concrete inequality, weight, or numerical certificate
from this version worth carrying into the active `k=48` program.

## Sources

- Yuhang Shi, *A Weighted Distribution of Primes and a New Unconditional Bound
  on Gaps Between Primes*, revised ResearchGate preprint (2025), Definition
  3.1, Theorem 4.1, Proposition 4.2, Lemmas 4.3--4.4, Section 5, Appendices A--B:
  https://www.researchgate.net/publication/393888742_A_Weighted_Distribution_of_Primes_and_a_New_Unconditional_Bound_on_Gaps_Between_Primes
- James Maynard, *Small gaps between primes*, Proposition 4.2:
  https://annals.math.princeton.edu/2015/181-1/p07
- D.H.J. Polymath, *Variants of the Selberg sieve, and bounded intervals
  containing many primes*, Theorems 22, 23, 27 and Section 8.4:
  https://link.springer.com/article/10.1186/s40687-014-0012-7
- James Maynard, *Primes in arithmetic progressions to large moduli II:
  Well-factorable estimates*, Definition 1 and Theorems 1.1--1.2:
  https://arxiv.org/abs/2006.07088
- Jared Duker Lichtman, *Primes in arithmetic progressions to large moduli,
  and Goldbach beyond the square-root barrier*, Theorem 1.4:
  https://arxiv.org/abs/2309.08522
- Michael Bennett, Greg Martin, Kevin O'Bryant, Andrew Rechnitzer, *Counting
  zeros of Dirichlet L-functions*, Theorem 1.1:
  https://arxiv.org/abs/2005.02989
- MIT Prime Gaps database: https://math.mit.edu/~primegaps/

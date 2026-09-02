# Exact Section 5 integral engine

This module reconstructs the exact integral layer sketched in Section 5.2 of
Julia Stadlmann, *Bounded gaps between primes*, arXiv:2608.31126v1, pp. 31--33.
It does not search for parameters, eigenvalues, or certificates.

## Reconstructed C and D formulas

Put `m = floor(1/delta)` and

```text
T_s(k) = {0 <= t_i <= delta, sum(t_i) <= 1},
T_b(k) = {delta <= t_i <= 1, sum(t_i) <= 1}.
```

For a non-negative exponent vector `a` and integer `b >= 0`, Section 5.2.1
defines `C_m(k,a,b)` and `D_m(k,a,b)` as the coefficient vectors in `delta` of
the corresponding monomial integrals. The paper only describes the recurrence.
The following closed forms are the recurrence after coefficient extraction.

For `S` contained in `{1,...,k}`, shift `t_i = y_i + delta` for `i in S`, and
write `q_i` for the selected power of `y_i` after expanding the shifted
monomial. Inclusion--exclusion gives

```text
C_m(k,a,b) = sum_{S, |S|<=m} (-1)^|S|
             sum_{0<=q_i<=a_i, i in S}
             [prod_{i in S} binom(a_i,q_i) delta^(a_i-q_i)]
             [b! prod_i e_i! / (b+k+sum_i e_i)!]
             (1-|S| delta)^(b+k+sum_i e_i),
```

where `e_i=q_i` on `S` and `e_i=a_i` otherwise. Terms with `|S|>m` have empty
interior in the chamber `1/(m+1) < delta <= 1/m`.

For `D`, shift every coordinate by `delta`. If `k>m`, the region is empty;
otherwise

```text
D_m(k,a,b) = sum_{0<=q_i<=a_i}
             [prod_i binom(a_i,q_i) delta^(a_i-q_i)]
             [b! prod_i q_i! / (b+k+sum_i q_i)!]
             (1-k delta)^(b+k+sum_i q_i).
```

Expanding the last powers produces rational coefficient vectors. These are
implemented by `small_cube_coefficients` and `large_simplex_coefficients`.
They exactly reduce to the paper's one-dimensional base cases and agree across
every chamber boundary `delta=1/m`.

Both formulas use the Dirichlet identity

```text
integral prod_i y_i^e_i (R-sum_i y_i)^b dy
  = b! prod_i e_i! / (b+k+sum_i e_i)! * R^(b+k+sum_i e_i).
```

## I/J/K assembly

`primegaps.integrals` decomposes each coordinate at `delta`, shifts all large
coordinates, and removes upper bounds on small coordinates by
inclusion--exclusion. Coordinates with the same large/small status collapse to
their group sums using the same Dirichlet identity. For `J` and `K`, integrating
the distinguished last coordinate leaves a rational piecewise polynomial in
the common-coordinate total and common large-coordinate total. A rational line
arrangement partitions this two-dimensional domain, and Green's theorem
integrates every polynomial cell exactly. No numerical geometry is used.

Example:

```python
from fractions import Fraction

from primegaps.integrals import ExactSupportParameters, exact_ijk_matrices, monomial

support = ExactSupportParameters.from_values(
    delta=Fraction(1, 4),
    epsilon=Fraction(1, 10),
    A=(Fraction(-1, 10), Fraction(2, 5)),
    B=((Fraction(1, 2),) * 4,),
)
basis = [
    monomial((0, 0)),
    {(1, 0): Fraction(1), (0, 1): Fraction(1)},
]
matrices = exact_ijk_matrices(basis, support)
assert matrices.I[0][0] == Fraction(1, 8)
```

A basis function is a sparse mapping from an exponent tuple to a rational
coefficient. All returned entries are `Fraction` objects. This matrix assembler
remains the low-dimensional audit oracle and uses an explicit status
decomposition. The symmetry-compressed `k=49` contraction is now implemented
in the frozen 240 reproducer and the accelerated exact backend; both were
checked against this reference kernel before their full runs. They reuse the
same C/D arithmetic without making this deliberately simple API scalable.

## Source errata and conventions

The arXiv v1 source has several defects that affect reconstruction:

- The displayed chamber condition says `floor(1/delta) in [m,m+1]`; the
  meaningful condition, used here, is `m=floor(1/delta)` or equivalently
  `1/(m+1)<delta<=1/m` away from boundaries.
- Section 5.2.2 interchanges some `r`, `s`, `a`, and `b` indices. The engine
  follows the domains, not those stray indices.
- The printed `K` definition constrains `t'_k` but neither integrates it nor
  places it in the integrand. The API uses the dimensionally consistent
  existential reading of that constraint and documents it explicitly.
  Stadlmann's published `H_1<=240` parameters have `c2=0`, so `K` does not
  affect that certificate.

## Verification

The tests check:

- reduction to all tested one-dimensional `C/D` base cases;
- exact chamber-boundary continuity;
- independent Gauss--Legendre quadrature for nontrivial two-dimensional
  `C/D` examples;
- independent nested quadrature for every entry of tiny `I/J/K` matrices;
- the full-simplex Dirichlet identity and a support with an active `B` cut.

The numerical checks are tests only. Matrix construction itself is rational.

# Candidate-independent J block operator

Repeated `J` evaluation no longer needs target signatures or polynomial
products.  The reusable operator in `primegaps.fast_exact.j_block` implements

\[
c \longmapsto f \longmapsto Kf \longmapsto M^T KMc = Jc,
\]

where `M` is the sparse marginal map and each `K[sigma,tau]` is a dense, small
signature-pair block.  `JBlockOperator.matvec` exposes this operation directly;
`as_scipy_linear_operator` makes it usable by Lanczos, LOBPCG, or Davidson-style
solvers without constructing the dense candidate-basis matrix.

## Exact algebra

For cached moments `M[r,s] = L*(x^r z^s)`, the contraction of two polynomial
families is

\[
L^*(P_i Q_j) = A H B^T,
\qquad H_{uv}=M_{a_u+c_v,b_u+d_v}.
\]

`contract_polynomial_families` implements this identity.  It does not construct
`P_i Q_j`, call FLINT polynomial multiplication, route a product to target
candidate dictionaries, or convert product coefficients back to Python
rationals.  The offline compiler inverts all target routes once into a
signature-pair index.  `compile_signature_pair_block` combines the target
functionals for one indexed signature pair and then performs this Hankel
contraction.  Small exact tests compare the resulting quadratic form with the
frozen exact verifier.

## Stable numerical construction

For discovery, candidate-space accumulation is now preferred over persisting
all signature-pair blocks.  In each integration batch, assemble `G = F M` with
`candidate_feature_values` and accumulate `G.T @ W @ G` once.  For a projection
`Q`, `projected_feature_values` computes `G Q` directly from the marginal
feature blocks and the sparse map, so `projected_J(Q)` needs only
`(G Q).T @ W @ (G Q)`.  The production builder is:

```bash
OPENBLAS_NUM_THREADS=16 PYTHONPATH=src sage -python \
  scripts/build_numerical_j_candidate_matrix.py \
  --degree 27 --log2-n 15 --batch-log2 9 \
  --support-config .research/work/.../support.json \
  --unrestricted-base .research/work/.../unrestricted-d27.npz \
  --output .research/work/.../legal-j-d27.json
```

The unrestricted baseline must use the same `(k, U, R)` and basis convention.
With it, the builder uses common random points and accumulates the correction as

\[
\operatorname{sym}((G_L+G_U)^T W (G_L-G_U)).
\]

Rows where `G_L == G_U` are removed before the matrix multiplication.  Direct
importance integration without the unrestricted control variate is not stable
for the raw high-degree coefficient basis.

At D=27 there are 2,526 candidate columns and 7,912 marginal feature columns.
Candidate-space accumulation retains 3,191,601 upper-triangle entries instead
of 31,546,512.  On a real 256-row legal batch, feature evaluation took 69 ms,
assembling `G` took 11 ms, and the single Gram update took 28 ms, with about
236 MiB peak RSS.  At D=21, the candidate-space correction reproduced the old
signature-block value for a fixed certificate to `1.05e-9` while building in
6.56 seconds.

`scripts/build_projected_numerical_j.py` implements the smaller `(GQ)` path and
`scripts/analyze_projected_j_replicates.py` cross-validates nested subspaces.
Residual-based Davidson enrichment additionally needs a stable streamed `JQ`
action.  Applying a raw full D=27 matrix is not acceptable: in the current
calibration it produced replicate Ritz values ranging from negative values to
the tens.  `scripts/optimize_projected_j.py` is therefore an experimental
diagnostic, not the current production optimizer.

For the final scalar certificate, computing legal `I` exactly is optional.  The
legal support is contained in the full `U`-simplex, so positivity gives
`I_legal <= I_simplex`.  For the current rationalized D27 candidate, the exact
normalized upper bound is `I_simplex = 0.9999955585079013`.  Thus an exact legal
`kJ` above this number certifies the quotient directly.  This replaces an
estimated hour-plus legal-`I` run by a roughly five-minute closed-form rational
calculation and leaves the legal J boundary correction as the only exact gate.
The exact target-free unrestricted J contraction takes 442.5 seconds at D27
and gives normalized `kJ_simplex = 1.0005162874604419`. Therefore the exact
boundary correction only has to exceed `-0.0005207289525405`; four independent
importance-QMC estimates lie between `-0.0002606091` and `-0.0000557749`.

### Exact-m sampling caveat

Translated exact-large-count simplices are an exact partition, but uniform
sampling inside them is not sufficient for high-degree optimization.  At D=21
and `2^15` samples per stratum it returned normalized `kJ = 1.0000632111`, while
the calibrated importance control variate returned `0.9999667252`; it also
missed exact I by `0.00108345`.  The sparse boundary rows were found, but the
high-leverage polynomial tail within those rows was not.  Exact-m partitioning
should be combined with a within-stratum tail proposal or deterministic
boundary integration before it replaces the importance control variate.

## Signature-block construction

The older reusable-operator path evaluates the orthogonal basis directly:

1. Build `MarginalMap` once from the optimizer basis.
2. Evaluate common-variable symmetric factors and integrated Jacobi factors.
3. Use `factorized_feature_values` to form columns indexed by
   `(signature, power, radial_degree)`.
4. Stream integration batches through `accumulate_feature_gram_blocks`.
5. Persist and hash the completed operator with `save_block_operator`.

This builds signature blocks directly.  It remains useful when a matrix-free
operator must be loaded repeatedly, but is no longer the preferred way to
accumulate a numerical discovery matrix.  The tracked builder is:

```bash
PYTHONPATH=src sage -python scripts/build_numerical_j_block_operator.py \
  --degree 21 --log2-n 15 --batch-log2 9 \
  --output-directory .research/work/j-block-d21
```

By default it starts from the exact-Gauss unrestricted block form and samples
only `legal - unrestricted`, implementing the numerical control variate.  Pass
`--direct` to estimate all of legal `J`.  Blocks are stored in one hashed `.npy`
array and memory-mapped on load.

## Production-size runtime measurements

The benchmark uses the real D=21 and D=27 sparse marginal layouts with
deterministic synthetic block values, so it measures operator storage and hot
loop cost independently of an offline integration method.

| Degree | candidate dimension | marginal features | upper blocks | float64 storage | one `Jc` |
|---:|---:|---:|---:|---:|---:|
| 21 | 846 | 2,370 | 9,730 | 21.8 MiB | 0.0180 s |
| 27 | 2,526 | 7,912 | 69,751 | 240.7 MiB | 0.133 s |

The D=21 timing is the mean of 20 applications; D=27 is the mean of five.  These
are runtime/dataflow measurements, not accuracy claims about synthetic blocks.

A real D=21 control-variate operator with 32,768 QMC points built in 15.55
seconds with 181 MiB peak RSS, loaded and hash-checked in 0.059 seconds, and
applied in 0.0187 seconds.  On the published unrestricted source direction, two
independent seeds differed by about 6.0 parts per million.  That is a discovery
diagnostic, not a proof or a global operator-error bound; increase the sample
count and use independent replicates near the decision boundary.

## Exact-cache conditioning result

Directly converting the exact monomial target cache to floating blocks is not a
safe production path.  For the real D=21 empty-signature block, the raw block
condition estimate was about `8.3e53`.  Contracting the floating block with the
published candidate gave `-5.24e-68`, while the exact target-group result is
`5.6422e-90`.  Expanding the source Jacobi basis into monomials before the
contraction was also unstable because the change of basis itself loses roughly
30 decimal digits.  Long double does not provide enough headroom.

Accordingly, `scripts/build_j_block_operator.py` is explicitly guarded as an
experimental low-degree/conditioning tool.  Production numerical blocks must be
built from directly evaluated orthogonal features, and the final promising
candidate must still go through the exact rational scalar verifier.

## Geometry and control-variate experiments

An explicit per-cell geometry-moment table was exact but did not accelerate the
one-target cache-fill benchmark.  The current cached polygon-moment primitive
already reuses geometry within a process: direct contraction took 13.89 seconds,
the explicit table took 14.08 seconds, and FLINT cross-correlation took 15.65
seconds.  The extra layer is therefore not promoted.

The exact decomposition `J_legal = J_simplex - Delta J_B` remains promising.
For the published k=49 geometry, 284 of 737 cells (38.5%) differ from the
unrestricted slices; 453 can be supplied entirely by the closed unrestricted
form.  This corrects an earlier 242/683 census that omitted unrestricted-only
statuses with more large coordinates than the finite legal `B_m` table.
`compile_signature_pair_block(control_variate=True)` and the scalar exact path
include those statuses and implement the algebraic correction.  A cache used
for this mode must include moments needed by both legal and unrestricted slice
polynomials, through every unrestricted large-count status.  The direct exact
target fallback remains slow; production cache fill still needs geometry-cell
moments shared across signature pairs before the Hankel contractions.  Merely
replacing the product by an exact `p.T @ H @ q` at pair scope is insufficient:
the old 23-pair D21 benchmark did not finish in 390 seconds with Python/GMP
contractions or in 180 seconds with a Sage/FLINT matrix contraction.  In both
cases, rebuilding each pair's density-weighted moment table dominated.

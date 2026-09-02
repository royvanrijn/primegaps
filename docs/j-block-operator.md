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

The production discovery path should evaluate the orthogonal basis directly:

1. Build `MarginalMap` once from the optimizer basis.
2. Evaluate common-variable symmetric factors and integrated Jacobi factors.
3. Use `factorized_feature_values` to form columns indexed by
   `(signature, power, radial_degree)`.
4. Stream integration batches through `accumulate_feature_gram_blocks`.
5. Persist and hash the completed operator with `save_block_operator`.

This builds signature blocks directly.  It never materializes the dense `J`
matrix, and candidate vectors do not enter the geometry/integration stage.
The tracked builder is:

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
For the published k=49 geometry, only 242 of 683 cells (35.4%) differ from the
unrestricted slices; 441 can be supplied entirely by the closed unrestricted
form.  `compile_signature_pair_block(control_variate=True)` implements the
algebraic correction, but a cache used for this mode must include moments needed
by both the legal and unrestricted slice polynomials.  This path is structurally
validated, not yet production-timed.

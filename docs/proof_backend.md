# Generalized-eigenvalue proof backend

This backend begins with supplied matrices. It does not construct or derive any
integrals.

## Contract

For exact symmetric matrices \(M_1,M_2\), the backend returns:

1. the best numerical generalized Rayleigh quotient found and its vector;
2. symmetry, positive-definiteness, conditioning, block, convergence, and
   backward-residual diagnostics;
3. the best exact quotient among bounded common-scale integer approximations to
   that vector;
4. when possible, a JSON certificate for
   \(c^T M_2 c > c^T M_1 c\).

The numerical stage is only a search. The exact verifier accepts the claim only
if all of the following hold:

- the primitive integer vector has the declared dimension;
- semantic SHA-256 hashes of both supplied matrices match the certificate;
- the recorded quadratic forms, difference, and quotient equal freshly
  computed exact values;
- \(c^T M_1 c>0\) and \(c^T(M_2-M_1)c>0\).

Changing any matrix entry, vector coordinate, or claimed exact value invalidates
the certificate.

## Python API

    from primegaps.exact_matrix import ExactSymmetricMatrix
    from primegaps.proof_backend import solve_and_certify
    from primegaps.certificate import verify_certificate

    m1 = ExactSymmetricMatrix.from_dense([[2, 0], [0, 1]])
    m2 = ExactSymmetricMatrix.from_dense([[5, 0], [0, 1]])

    result = solve_and_certify(m1, m2)
    print(result.numerical.quotient)
    assert result.certificate is not None
    verify_certificate(m1, m2, result.certificate)

ExactSymmetricMatrix.from_sparse accepts upper- or lower-triangular
(row, column, value) triples and combines duplicates. Fraction values are
supported. Packed and sparse representations of the same exact matrix have the
same semantic hash.

## File formats and CLI

The JSON input stores a common denominator and either a row-major packed upper
triangle or sorted nonzero upper entries:

    {
      "format": "primegaps-exact-matrix-pair-v1",
      "m1": {
        "dimension": 2,
        "denominator": 1,
        "storage": "packed-upper",
        "numerators": [2, 0, 1]
      },
      "m2": {
        "dimension": 2,
        "denominator": 1,
        "storage": "sparse-upper",
        "entries": [[0, 0, 5], [1, 1, 1]]
      }
    }

For large dense matrices, NPZ avoids the memory overhead of millions of JSON
integer objects. It must contain four arrays:

- m1_numerators and m2_numerators: one-dimensional packed-upper integer arrays;
- m1_denominator and m2_denominator: scalar positive integer arrays.

NPZ numerators are limited to the NumPy integer dtype used in the file. Exact
quadratic-form accumulation converts every term to an unbounded Python integer,
so it cannot overflow.

    primegaps-proof solve matrices.npz \
      --certificate certificate.json \
      --method auto \
      --dense-threshold 1200 \
      --max-iterations 180 \
      --max-scale 1000000

    primegaps-proof verify matrices.npz certificate.json

Solve exits with status 2 and does not write a certificate when the best
rationalized candidate fails the strict inequality. Verify performs no
floating-point work and does not repeat the expensive eigenvalue search.

## Numerical regimes and diagnostics

The solver first finds connected components of the joint off-diagonal nonzero
pattern of \(M_1\) and \(M_2\). Exact simultaneous blocks are solved separately;
the largest block eigenvalue is the global answer. Every dense block is
diagonally equilibrated to unit diagonal in \(M_1\), then handled by:

- dense: Cholesky whitening and a complete symmetric eigensolve;
- iterative: matrix-free, fully reorthogonalized Lanczos, which does not form
  the dense whitened operator;
- auto: dense up to the configured threshold, iterative above it.

SciPy sparse inputs use scipy.sparse.linalg.eigsh. Install the sparse extra with
pip install -e '.[sparse]'.

Diagnostics report relative symmetry defects, diagonal range, estimated smallest
and largest eigenvalues and condition number of \(M_1\), estimated condition
after equilibration, Cholesky diagonal range, generalized backward residual,
iteration/convergence status, and simultaneous block sizes. Condition numbers
for the NumPy path are power-iteration estimates, not rigorous bounds. The exact
certificate remains decisive even when numerical conditioning is poor.

## Degree and memory envelope

For the repository's \(k=49\) symmetric basis:

| Degree | Dimension | One float64 dense matrix | One int64 packed triangle |
|---:|---:|---:|---:|
| 21 | 846 | 5.5 MiB | 2.7 MiB |
| 22 | 1,041 | 8.3 MiB | 4.1 MiB |
| 23 | 1,236 | 11.7 MiB | 5.8 MiB |
| 24 | 1,508 | 17.3 MiB | 8.7 MiB |
| 25 | 1,780 | 24.2 MiB | 12.1 MiB |
| 26 | 2,153 | 35.4 MiB | 17.7 MiB |
| 27 | 2,526 | 48.7 MiB | 24.3 MiB |
| 28 | 3,034 | 70.2 MiB | 35.1 MiB |
| 29 | 3,542 | 95.7 MiB | 47.9 MiB |
| 30 | 4,226 | 136.3 MiB | 68.1 MiB |

At degree 30, forming two inputs plus Cholesky and a dense whitened operator can
approach a gigabyte once library workspaces and copies are included.
Packed-int64 NPZ plus iterative solving avoids the extra whitened matrix;
sparse/block storage reduces both exact replay time and numerical memory when the
supplied matrices genuinely have that structure.

Candidate rationalization uses one common scale rather than independent
continued-fraction denominators. This keeps coefficient growth predictable and
tries scales 1, 2, and 5 times powers of ten through max_scale, selecting the
best quotient by exact arithmetic.

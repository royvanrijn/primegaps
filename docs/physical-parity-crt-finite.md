# Reduced CRT-coloured parity experiment

## Verdict

The first nondegenerate finite coefficient compilation gives a qualified
**GO** for a global residue-coloured operator and a **NO-GO** for local
prime/semiprime cancellation.

Across `X=8,16,32,64` million, the blockwise cancellation ratio after CRT
aggregation is between `8.5131` and `37.3645`.  It is therefore large and
persistent under the requested blockwise interpretation.  But within each
aligned `(q,a)` block the prime/semiprime ratio is only `1.0262--1.0407`.
Almost all of the surviving gain comes from cancellation **across** distinct
CRT-coloured blocks.

There is a necessary notational warning.  If `E_prime` and `E_semiprime` mean
the two fully summed scalar errors, then

```text
(|E_prime| + 2|E_semiprime|) / |E_prime - 2 E_semiprime|
```

is invariant under exact CRT aggregation.  In this experiment it is exactly
`1` at all four scales because the two scalar sector errors have opposite
signs, so the minus sign reinforces them.  The meaningful before/after
quantity is instead the blockwise ratio

```text
R(B) = sum_(b in B) (|P_b| + 2|S_b|)
       / |sum_(b in B) (P_b - 2S_b)|,
```

with raw blocks `b=(d,e)` before aggregation and `b=(q,a)` afterwards.

Thus this finite model does not support a theorem based on aligned
prime-versus-semiprime cancellation.  It does identify the genuinely global
operator

```text
sum_(q,a) c_i(q,a) (E_prime(q,a) - 2 E_semiprime(q,a))
```

as numerically viable enough to study further.

## Reduced model

The compiler uses the frozen `k=39` physical vectors from the `4096`- and
`8192`-interval viability meshes.  It restricts the 39-dimensional symmetric
trial to three active non-target coordinates and attaches four coarse physical
fragment sizes

```text
0.18, 0.22, 0.26, 0.30.
```

The arithmetic prime atoms are `5,7,11,13`.  The three colours use the actual
admissible `H39` shifts `2,36,48` relative to target shift `0`; their local
target-minus-shift residues are nonzero and distinct at every atom.  At the
smallest finite scale, every prime atom and their full product `5005` lie
inside the recorded one-prime and pair-modulus exponent envelopes.

For every squarefree divisor configuration, the model evaluates the frozen
physical trial and includes its Möbius sign.  The resulting `Lambda_d` are
normalized by

```text
sum_d Lambda_d^2 = 1.
```

Compatible pairs may share a prime only in the same colour.  Their union gives
the modulus and colour word; CRT gives `a`, and collision aggregation gives

```text
c_i(q,a) = sum_((d,e) -> (q,a)) Lambda_d Lambda_e.
```

The reduced model has 211 supported divisor configurations, 8,455 compatible
ordered pairs and 256 occupied CRT states.  On 4,096 test integers, the direct
square and compiled coefficient sum agree to at most `2.22e-16`.

For the arithmetic side, let `z=floor((2X)^(1/3))`.  The semiprime carrier
contains `m=pq` in `(X,2X]` with both primes greater than `z`; three such
factors cannot occur.  The sector errors are the exactly enumerated centered
progression discrepancies

```text
E_s(q,a) = # {m in sector s: m = a (mod q)} - #sector_s / phi(q).
```

The physical mesh and the small CRT labels are deliberately decoupled.  This
is a structural finite model, not a claim that the integers `5,7,11,13` have
the displayed physical logarithmic sizes.

## Main measurements

The primary column below uses the finer `k39_n8192` vector.

| `X` | `sum_(q,a)|c|^2` | raw/grouped coefficient `l1` | `R` before CRT | `R` after CRT | within-block sector ratio | across-CRT ratio | literal scalar ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8,000,000 | 0.415296895 | 3.638754 | 54.4183 | 18.8684 | 1.02899 | 18.3369 | 1.00000 |
| 16,000,000 | 0.415296895 | 3.638754 | 115.4735 | 37.3645 | 1.02625 | 36.4089 | 1.00000 |
| 32,000,000 | 0.415296895 | 3.638754 | 61.6648 | 19.5725 | 1.03504 | 18.9099 | 1.00000 |
| 64,000,000 | 0.415296895 | 3.638754 | 27.5422 | 8.51314 | 1.04067 | 8.18043 | 1.00000 |

CRT collision aggregation reduces the gross block mass by factors
`2.8841, 3.0904, 3.1506, 3.2352`, respectively, while leaving a large
blockwise ratio.  The two independently generated physical mesh vectors agree
closely: the coefficient `l2` values differ by `0.281%`, the CRT collision
ratio by `0.065%`, and the post-CRT blockwise ratio by at most `4.77%`.

## Singular-value decay

To preserve the ragged colour fibres, rows are prime-subset moduli and columns
are full colour words.  A value on a restricted word is repeated over its
global extensions with factor `3^(-inactive/2)`.  Hence the table is an
isometric padding and its squared Frobenius norm is exactly
`sum_(q,a)|c(q,a)|^2`.

For the coefficient table, the top four singular directions contain
`93.2788%` of Frobenius energy and the top eight contain `99.6287%`; the stable
rank is `2.2051`.  Multiplication by the finite projected progression error
makes the operator less compressible:

| `X` | projected stable rank | top-4 energy | top-8 energy | numerical rank at `10^-6` |
|---:|---:|---:|---:|---:|
| 8,000,000 | 2.2474 | 77.1477% | 94.3256% | 15 |
| 16,000,000 | 3.6674 | 67.3738% | 89.0370% | 15 |
| 32,000,000 | 5.3756 | 61.3190% | 86.0918% | 15 |
| 64,000,000 | 2.4722 | 76.3671% | 90.2271% | 15 |

The coefficient compiler therefore has strong low-rank decay, while the
arithmetic operator has useful but not spectacular decay.

## Interpretation and next gate

The user-specified viability gate is passed only in its global blockwise
sense.  It is not passed by the literal scalar sector quotient, and it is not
passed locally inside `(q,a)`.  Any theorem should therefore keep the complete
CRT-coloured sum intact; a blockwise triangle inequality would discard the
observed factor `8--37`.

The next finite experiment should enlarge both the prime-atom universe and the
number of active colours using a compressed state dynamic program, then add
residue-profile and sign/permutation controls.  The decisive question is
whether the across-CRT ratio remains separated from those null controls once
the reduced model is no longer a 256-state slice.

## Reproduction

The explicit computation is separate from replay and takes a few seconds on
the recorded machine:

```bash
.venv/bin/python experiments/physical_parity_crt_finite.py \
  --output experiments/physical_parity_crt_finite.json
```

Cheap replay checks hashes, the compiler identity and all recorded gates
without rerunning the sieve:

```bash
.venv/bin/python scripts/check_physical_parity_crt_finite.py
```

This remains a binary64 finite experiment.  It is neither a distribution
estimate nor evidence for an asymptotic power saving by itself.

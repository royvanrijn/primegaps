# Faster exact 240-certificate calculations

The frozen evaluator in `reproduction/240/independent-reproducer` remains the
gold oracle. The accelerated backend in `primegaps.fast_exact` changes loop
order and arithmetic representation, but not the support, candidate, or exact
integrals. It is intended to make repeated degree and support experiments
practical; it is not evidence for a new prime-gap bound by itself.

## Implemented improvements

1. **Closed-zero status recurrence.** A monomial signature has only a few
   positive exponents and roughly 40 zero-exponent coordinates. The positive
   coordinates are processed by the original exact recurrence. All `z` zero
   coordinates are then collapsed in one multinomial convolution:

   ```text
   choose(z,a) choose(z-a,b) (-1)^b,
   ```

   where `a` zero coordinates enter the large block and `b` enter the shifted
   inclusion-exclusion block. This is an exact regrouping, not an asymptotic
   approximation. The positive core is independent of `k`, so the same core is
   reused for nearby dimensions.

2. **Pair-first J assembly.** The frozen evaluator visits each target
   signature and repeats its feature-pair products. The accelerated evaluator
   forms each unordered feature-pair product once per support cell and routes
   it to every target signature in the current chunk. For the degree-21
   candidate there are 9,730 unordered feature pairs but 165,166 target-pair
   contributions, or 16.974 target outputs per pair on average.

3. **Compiled polynomial kernels.** A bivariate polynomial is encoded as a
   univariate FLINT polynomial by `x^a z^b -> q^(a*s+b)`, with stride `s` larger
   than every possible `z` exponent. Exponent addition then has no carries, so
   FLINT multiplication is exactly the original bivariate multiplication. Both
   rational (`QQ`) and prime-field (`GF(p)`) backends use this encoding.

4. **Modular arithmetic and CRT.** J contractions can run unchanged over an
   explicitly selected prime field. Deterministic 64-bit prime checking, CRT,
   and rational reconstruction are included. Reconstruction requires explicit
   numerator and denominator bounds and rejects a result unless
   `2*N*D < product(primes)`, which is the uniqueness condition. This path is
   useful only when enough primes and independently justified bounds are
   supplied.

5. **Candidate-independent exact caches.** `IMomentCache` stores exact
   `(signature, radial slack)` moments. `JFunctionalCache` stores the action of
   a fixed target density/support cell on candidate monomials. Both files are
   context-bound, append-only, detect conflicting values, and allow candidate
   coefficients to be changed without rebuilding existing moments.
   The J runner also keeps a bounded per-worker LRU of compiled slice
   polynomials so target chunks do not repeatedly encode the same slice.

6. **Reuse across degree and dimension.** Raising degree requests only new
   slack powers or monomial exponents from the two moment caches. Changing from
   `k=48` to `k=49` reuses the positive-signature status recurrence and applies
   a different closed zero block. Support changes deliberately invalidate the
   moment-cache context.

## Measurements

The historical exact jobs and accelerated jobs were run on the same host. Each
group records elapsed worker-task time, not an operating-system CPU counter.
Summing those values approximates CPU time when workers are fully scheduled but
overstates it under contention. Wall times are also affected by worker count
and concurrent load, so worker-task and wall comparisons are shown separately.

| calculation | frozen | accelerated | observed change |
| --- | ---: | ---: | ---: |
| `k=48` I, summed task time | 4.92 h | 1.285 h | 3.83x less task time |
| `k=48` I, wall | 4,427 s / 4 workers | 290 s / 16 workers | 15.3x faster wall |
| `k=49` I, summed task time | 5.08 h | 1.196 h | 4.25x less task time |
| `k=49` I, wall | 4,578 s / 4 workers | 269 s / 16 workers | 17.0x faster wall |
| first 32 `k=48` J groups, one worker | 2,410.6 s | 355.8 s | 6.78x faster |
| full `k=48` J prototype, summed task time | 27.60 h | 7.179 h | 3.84x less task time |
| full `k=48` J, wall | 4,143 s / 24 workers | 2,251 s / 12 workers | 1.84x faster wall |
| manifest-bound `k=48` J under contention | 27.60 h | 10.703 h | 2.58x less task time |
| same contention run, wall | 4,143 s / 24 workers | 3,296 s / 12 workers | 1.26x faster wall |
| full `k=49` J, summed task time | 28.84 h | 8.313 h | 3.47x less task time |
| full `k=49` J, wall | 4,191 s / 24 workers | 2,631 s / 12 workers | 1.59x faster wall |
| same 32 J groups, one 61-bit modular image | 355.8 s rational | 357.3 s modular | no speedup |
| one cached `k=48` I group, worker time | 0.844 s build | 0.00044 s replay | about 1,900x replay |

Every one of the 2,714 accelerated I rows matched the frozen rational numerator
and denominator at both `k=48` and `k=49`. All 2,714 pair-first/FLINT `k=48` J
rows also matched, as did all 2,714 `k=49` J rows. The accelerated finalizer
reproduced the exact published-certificate quotient and positive `49J-I`, then
the same manifest-bound code reproduced the unchanged k=48 deficit. The k=48
prototype ran with near-full worker scheduling; the later manifest-bound k=48
run was deliberately left unchanged while host load rose sharply, so its
summed-task and wall speedups are conservative contention measurements.

For representative slice polynomials with 276 and 231 terms, 100 polynomial
products took 5.117 seconds in the Python dictionary kernel and 0.157 seconds
in the carry-free FLINT encoding (32.6x for the isolated kernel). End-to-end J
speedup is smaller because density construction, slice construction, geometry,
and exact moment integration remain.

The modular run at `p=2305843009213693951` matched all 32 frozen rational
residues, but was not faster than the compiled rational run. Since the observed
group numerators and denominators require many 61-bit images before bounded
rational reconstruction is unique, this modular implementation is currently a
correct experimental backend, not the recommended fast path. It becomes
promising only if several primes are evaluated together while sharing geometry
and moment work.

## Commands

The exact runners require SageMath (for FLINT) and `gmpy2` (the latter is also
available through `pip install -e '.[exact]'`). Outputs are
checkpointed JSONL files; a sidecar manifest binds `k`, backend, candidate hash,
and frozen-verifier hash. Reusing a path with a different context is rejected.

```bash
PYTHONPATH=src sage -python scripts/run_fast_exact_i.py \
  --k 49 \
  --candidate reproduction/240/independent-reproducer/candidate-k49-d21.json \
  --workers 16 --moment-cache .research/work/k49-I-moments.jsonl \
  --output .research/work/fast-k49-I.jsonl

PYTHONPATH=src sage -python scripts/run_fast_exact_j.py \
  --k 49 \
  --candidate reproduction/240/independent-reproducer/candidate-k49-d21.json \
  --workers 12 --chunk-size 32 --compiled \
  --output .research/work/fast-k49-J.jsonl

PYTHONPATH=src python scripts/finalize_fast_exact.py \
  --k 49 --i-groups .research/work/fast-k49-I.jsonl \
  --j-groups .research/work/fast-k49-J.jsonl \
  --output .research/work/fast-k49-result.json
```

One modular image uses the same J runner:

```bash
PYTHONPATH=src sage -python scripts/run_fast_exact_j.py \
  --k 49 \
  --candidate reproduction/240/independent-reproducer/candidate-k49-d21.json \
  --workers 12 --chunk-size 32 --modular-prime 2305843009213693951 \
  --output .research/work/fast-k49-J-p2305843009213693951.jsonl
```

Repeat `--modular-prime` in one invocation to share density, slice, and support
geometry work across several CRT images; the output then carries parallel
`primes` and `residues` arrays.

After enough distinct primes have been run, exact rows can be reconstructed
with mandatory proof bounds:

```bash
PYTHONPATH=src python scripts/reconstruct_modular_exact.py \
  .research/work/fast-k49-J-p*.jsonl \
  --bounds .research/work/j-rational-bounds.jsonl \
  --output .research/work/fast-k49-J-reconstructed.jsonl
```

The accelerated result should always be checked row-for-row against the frozen
backend before its first use in a new support or basis regime. For a local
reproduction tree containing the oracle checkpoints, use
`scripts/compare_exact_checkpoints.py`; it rejects missing, foreign, duplicate,
or unequal rational rows unless an explicit prefix comparison is requested.

## Audited outputs

The manifest-bound full runs are preserved as immutable local research objects:

| artifact | SHA-256 |
| --- | --- |
| accelerated `k=49` J rows | `8a3fb5d9311257f9a872093f91d6e046e0a2379255b53201aabfee6143770ecb` |
| accelerated `k=49` result | `2387d85c0d14b58d67bb9786fc1794fc3408bfc42d14e7a8c4fe5f33c271929a` |
| accelerated `k=48` J rows | `f14eafac9660cb3d8d9e97d6168ca85f752eb9762dc9355eec85cfe6da4ea00e` |
| accelerated `k=48` result | `465d83562813c86f92e36963c1cfe67d039143af4a93a13c9e6fd5ea4f3241cb` |
| combined experiment summary | `98cfdebe8b98d3e17f504e66b60ac2ad111a144c18bda5b143f9ddd24cafbf43` |

The canonical positive research record is
`d577ca9c7f309207c00c453ff8f605b9da5cc5ec65cda3db1cc9686f4cf7f47e`.
The modular no-speedup finding is separately recorded as negative record
`d14689e3a833320e7652d0c10baef907a01a7541cfda37912e16b1bb9ab24b1f`.

## Reusing moments

The cache classes intentionally require an explicit context dictionary. Put
every mathematical input that changes a moment—`k`, `delta`, total cap, the
`B_m` caps, and evaluator schema/version—into that context. For I, request only
missing slacks with `IMomentCache.missing`, calculate them with
`fast_i.signature_moments`, and append them. A later candidate contracts its
atoms with `IMomentCache.evaluate_atoms` without support integration.

For J, `fast_j.density_weighted_moments` turns one fixed target density and
support cell into a functional on requested `(x_power,z_power)` monomials.
Store those values under a stable target/cell identifier with
`JFunctionalCache`; `evaluate` then contracts a new candidate polynomial.
When degree rises, `missing` returns only the exponents not already on disk.
The cache detects any attempt to append a different rational for an existing
key. Cache files are acceleration artifacts under `.research/`, never trusted
certificate inputs by themselves.

`run_fast_exact_i.py --moment-cache PATH` performs this missing-moment logic in
parallel workers and appends in the parent process, avoiding concurrent writes.
In a one-group exact check, cache replay returned the identical numerator and
denominator while reducing worker time from 0.844 seconds to 0.00044 seconds;
process startup then dominated the 0.033-second wall time.

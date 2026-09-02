# Frozen exact 240 reproduction

This directory contains the fixed-input evaluator used for the independent
`k=49`, `D=21` reproduction and the unchanged `k=48` baseline. It performs no
optimization. The source and candidate files are byte-for-byte copies of the
frozen versions bound by the sealed `k=49` and `k=48` manifests.

The calculation requires Sage's Python with `gmpy2`; it is intentionally an
explicit, expensive operation. From the repository root, set a disposable
output directory and run each `I`/`J` pair separately:

```bash
mkdir -p .research/work/reproduction-240/rebuild
PYTHONPATH=src sage -python -u reproduction/240/independent-reproducer/run_parallel_exact.py \
  --kind I --k 49 \
  --candidate reproduction/240/independent-reproducer/candidate-k49-d21.json \
  --binding reproduction/240/independent-reproducer/candidate-k49-binding.json \
  --workers 4 \
  --output .research/work/reproduction-240/rebuild/k49-I.jsonl
PYTHONPATH=src sage -python -u reproduction/240/independent-reproducer/run_parallel_exact.py \
  --kind J --k 49 \
  --candidate reproduction/240/independent-reproducer/candidate-k49-d21.json \
  --binding reproduction/240/independent-reproducer/candidate-k49-binding.json \
  --workers 24 \
  --output .research/work/reproduction-240/rebuild/k49-J.jsonl
```

Finalize without recomputing any group:

```bash
PYTHONPATH=src:reproduction/240/independent-reproducer sage -python \
  reproduction/240/independent-reproducer/finalize_exact.py \
  --k 49 \
  --i-groups .research/work/reproduction-240/rebuild/k49-I.jsonl \
  --j-groups .research/work/reproduction-240/rebuild/k49-J.jsonl \
  --candidate reproduction/240/independent-reproducer/candidate-k49-d21.json \
  --binding reproduction/240/independent-reproducer/candidate-k49-binding.json \
  --output .research/work/reproduction-240/rebuild/k49-result.json
```

For `k=48`, replace `49` by `48` in the arguments and output names. The worker
counts affect only runtime. Checkpoints are resumable. During computation the
runner rejects mixed contexts, duplicate/unexpected signatures, and invalid
record schemas; finalization additionally checks row `k`/kind and the complete
candidate-derived signature sets. Once sealed, manifest hashes reject later
alteration of the completed records or summary.

The independently recorded outputs are
[`exact-k49-d21-result.json`](exact-k49-d21-result.json) and
[`exact-k48-d21-result.json`](exact-k48-d21-result.json). See
[`docs/reproduction-240.md`](../../docs/reproduction-240.md)
for the mathematical conventions, trusted numbers, provenance hashes, and
independent checks.

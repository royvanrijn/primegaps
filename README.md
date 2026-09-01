# primegaps

Exploratory code around Julia Stadlmann's 2026 bounded-prime-gaps paper
[arXiv:2608.31126](https://arxiv.org/abs/2608.31126), which improves the unconditional
bound from `H_1 <= 246` to `H_1 <= 240`.

The immediate goal is to reproduce the paper's `k=49`, degree-21 generalized
eigenvalue certificate and then probe `k=48,47,...`. This repository currently
contains the first, independently checkable layer of that reconstruction:

- the published support `T_k(delta, A, B, epsilon)`;
- generation/counting of the symmetric-polynomial basis used by Polymath/Stadlmann;
- a geometry diagnostic for support sampling;
- notes on what is and is not yet reproduced.

**Nothing in this repository currently proves a bound below 240.** The exact
Section 5 integral recurrence / `M1,M2` construction is still the critical missing
piece.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
python -m primegaps.scan
pytest
```

The scan prints the published parameters, basis dimensions for degree 21 and 27,
and a Monte-Carlo support-geometry diagnostic. The latter is deliberately *not* a
sieve score; it is useful mainly for rejecting naive "maximize support volume"
strategies.

# Exploratory distance below 212

This is an intentionally non-rigorous first screen of the exact Table 3 support
from the preliminary BGP212 draft.  It asks how far the same degree-21 framework
appears to be from `k=44`, `k=43`, and `k=42` before any further support or
number-theory improvement.

The existing projected legal-minus-unrestricted engine was run in the top eight
unrestricted generalized eigendirections, with degree 21 and three independent
importance-control seeds (`2^16` samples each).  The paper reports the exact
physical `k=45` quotient

```text
4.00438409833460131937... / 4 = 1.0010960245836503...
```

so the mean `k=45` numerical error can be used as a crude additive calibration.
It is not a rigorous error model.

| k | raw mean | raw SD | calibrated heuristic | apparent deficit | tuple target |
|---:|---:|---:|---:|---:|---:|
| 45 | 1.00097118 | 0.00010105 | 1.00109602 (paper) | passes | 212 |
| 44 | 0.99691961 | 0.00043863 | 0.99704445 | 0.00296 | 210 |
| 43 | 0.99270035 | 0.00109424 | 0.99282520 | 0.00717 | 200 |
| 42 | 0.98755078 | 0.00026668 | 0.98767562 | 0.01232 | 196 |

These numbers strongly change the immediate research ranking:

- `k=44` looks like a roughly three-per-thousand variational problem, not a
  two-percent jump.  Higher degree plus a small legal support improvement is a
  credible route to 210.
- `k=43` appears below threshold by roughly seven-per-thousand.  It is not ruled
  out by this framework; it likely needs joint optimization of degree, `A`,
  `delta`, and the full `B_m` staircase, or another arithmetic relaxation.
- `k=42` is about one percent away in this small projected model and is a useful
  medium-term target, not an immediate exact run.

## Important limitations

The screen does not contain the authors' 846-entry exact vector, does not search
the full degree-21 space, and does not certify the legal support.  A 64-direction
projection was less stable and under-reproduced the known `k=45` value, so the
small projected space is being used only as a local distance indicator.  No
statement here can be promoted into a prime-gap result.

The immediate replay path remains:

1. obtain or reconstruct the BGP212 degree-21 witness and reproduce `k=45` exactly;
2. upgrade the exact distribution oracle to the paper's unsimplified divisor
   windows and 455 continuum-packing certificates;
3. use the stable projected engine at `k=44`, then exactify only a candidate with
   comfortable replicated margin;
4. test the repository's incomplete-range Type-IIc idea against the new
   `delta=41/2500` family, where the next Type-II wall allows at most
   `A_1=2059/8000` without another analytic improvement;
5. search multi-band and optimized cap staircases jointly rather than sweeping
   one parameter at a time.

Machine-readable values are in
`reproduction/212/exploratory-d21-frontier.json`.

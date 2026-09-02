# Numerical degree-21 screening engine

These files preserve the numerical vector-bank and support-boundary QMC engine
used for support screening. They are separate from the exact 240 certificate:
floating point locates promising support geometries, while exact arithmetic is
still required for a theorem.

`scripts/sweep_p3ii_delta_frontier.py` is the supported driver for the
`P3.II.delta` curve. Its `--shifted-strata` estimator should be used for recorded
results. The older Dirichlet-tilted estimator remains available only to replay
historical screens; its denominator correction can have a severe long tail.

The numerical implementation was promoted from the research shards associated
with the independently reproduced degree-21 calculation so the recorded
frontier can be rebuilt without relying on disposable `.research/` source.


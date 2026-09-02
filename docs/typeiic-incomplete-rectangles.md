# Type-IIc incomplete-rectangle saving

Status: research theorem assembled and independently checked on 2026-09-02.
The production distribution oracle has not yet been changed; the proof should
be typeset and human-reviewed first.

## Result

For the originating amended 2026 Type-IIc family, with
`delta=7/250` and `xi2=2/5`, the actual incomplete rectangles support the
replacement terminal condition

\[
 6-22\gamma+72\delta+216\omega<0.
\]

Since `omega=A-1/4`, this is

\[
 \delta < \frac{11\gamma}{36}+\frac23-3A.
\]

The limiting exponent calculation gives `A<856/3375=0.2536296296...`.
Because several analytic estimates retain an unspecified `x^O(epsilon)`
factor, the checked promotion point is the safer rational value

\[
 A=\frac{2029}{8000}=0.253625.
\]

It has unpadded structural margin `1/1000` and lies above the numerical
`lambda_48=1` crossing `0.2536077308`. This establishes the requested
Type-IIc saving analytically. It does not turn the randomized numerical
`k=48` screen into an exact variational certificate.

## Where the saving comes from

Write

\[
 G=(\lambda,\widetilde\lambda)=w_2g,
 \qquad \lambda=Ga,\quad \widetilde\lambda=Gb,
 \qquad (a,b)=1.
\]

The originating congruence forces `G|k`. With `k=Gj`, the lifted variables
satisfy

\[
 aY-bX=d\{cj+(a-b)B\},
\]

and the phase becomes

\[
 e_m\!\left(\frac{AG\{cj+(a-b)B\}}{XY}\right).
\]

This yields two complementary estimates.

1. For large `g`, affine-line completion and the exact divisor split give
   \[
   \Sigma_5\ll x^{O(\varepsilon)}L\left\{
     (w_2,m)T\left(\frac{ND}{m}+N+D\right)
     +\frac{mK}{g}+\frac{D^2}{q_0}
     +D\sqrt{\frac{mT}{q_0}}
   \right\}.
   \]
   The only critical term is `mK/g`.

2. For small `g`, a two-dimensional Fourier slice of the actual lifted phase
   is a rank-three Kloosterman pullback, minus a lower-weight deleted-line
   term. The translated correlations have square-root cancellation. A
   Plancherel average is essential: the zero-frequency transform has a real
   square-root spike, so a pointwise supremum would be false. Exact
   reduced-conductor completion handles the complementary large numerator-gcd
   range.

The geometric input is the connected `SL_3` monodromy and determinant of the
rank-three Kloosterman sheaf, together with trace-function
quasi-orthogonality; see [Katz](https://web.math.princeton.edu/~nmk/Katz-GKM.pdf)
and [Fouvry--Kowalski--Michel](https://arxiv.org/abs/1210.0851v5).

## Arithmetic bookkeeping that is part of the theorem

Several plausible shortcuts are false and are not used.

- The amended 2026 modulus is not assumed to divide `P(x^delta)`. The
  A-process factor is the named source divisor
  \[
  q_A=\frac{r_1}{(r_1,t s q_0)}.
  \]
- A power-sized auxiliary `w1` mask cannot be bounded by Fourier `L1` without
  a power loss. Möbius inversion is performed before the lift; the resulting
  sublattices only shorten the physical ranges and cost `x^o(1)`.
- Put `t=(w2,m)`, `s=(ell_j,m/t)`, `qt=(q0,t)`, `qs=(q0,s)`, and
  `qr=q0/(qt*qs)`. The true residual conductor and `j` period are
  \[
  R_0=\frac{m}{tsq_r},\qquad Q_j=[s,q_0],\qquad
  R_0Q_j=\frac{m q_t}{t}.
  \]
  The false clean-coefficient identity `R0*Qj=m` is not used.
- At phase-dead `t` primes, the post-lift graph mask can have normalized
  Fourier `L1` mass of order `p`. Its precise source factor is
  \[
  q_{3,t}=(q_3,t)\mid \frac{t}{(q_0,t)}.
  \]
  The residual part `q3/q3t` remains in the residual trace/class treatment.

After charging this sharp mask loss, if
`tau=log_x(t)` and `kappa_t=log_x((q0,t))`, the required A-factor exponent
changes by at most

\[
 \kappa_t-7\tau\le -6\tau.
\]

Removing `t` from `q_A` costs at most `tau`, so the structural margin improves
by at least `5 tau`. Thus `t=1` is the worst face. At `A=2029/8000`, the
remaining cube-capacity margin is at least `0.0539999984`.

## Verification boundary

The final hostile audit returned YES at `A=2029/8000` with strict analytic
epsilon. Its canonical record is
`f41320375247d751972ec1d6a959abb0914ea87831864c7c73acdb4855cd3d17`;
the independent `q3t` verification is
`0ed15026fd496a4efff6e56e40fe7397de864d25b948a13eb468e904c74660b7`.

All corrected component replays and the 15 focused distribution/frontier
tests pass. The audit does not endorse equality at `856/3375`, the exact
formally padded limiting decimal, or a generic improvement for arbitrary
Type-II sums.


# Direct source-barrier viability

## Conclusion: INCONCLUSIVE at the Stage-1 semantic gate

The requested probability/expectation comparison was not run. The stated bad
event includes a largest-fragment witness, while PrimeGaps186's outer
order-\(5/2\) low/rank-two/third-factorial cover is proved only for
**nonlargest** witnesses after a separate cap argument has handled the largest
fragment. This is not a floating-point discrepancy: the exact one-point
configuration \(\Pi=\{z\}\) has

\[
  \mu z-T=
  \frac{32664308770638763}{169342645250703360}>0,
  \qquad C(\{z\})=0.
\]

Thus \(C\not\geq 1_B\) for the \(B\) in the experiment specification, so the
proposed nonnegative excess decomposition and the ratio
\(\mathbb E[C]/\mathbb P(B)\) do not measure overcoverage. In accordance with
the instruction to stop on any Stage-1 majorization violation, no Monte Carlo,
DP comparison, stress test, or weighted root-side check was performed.

## Extracted source group and slice

The snapshot is
[`openai/PrimeGaps186@61340d0`](https://github.com/openai/PrimeGaps186/tree/61340d0b74163003b32756bb16e91d9209a5e330).
The selected group is `G2`, called `outer_h25` by the Python certificate and
index `g = 1` in Lean. It has one aligned cap piece, with radial cell-index sum
\(95599\leq r\leq98263\), so there is no cap ambiguity within the group.

| parameter | exact value | decimal |
|---|---:|---:|
| \(n\) | \(40\) | 40 |
| \(\xi\) | \(1038826867984921151804142858423732482601802307/30904730932085735018956267409864291787032494080000\) | \(0.00003361384605702541\) |
| \(z\) | \(31942200065/64511729664=46580h\) | \(0.4951378645614731\) |
| \(T\) | \(14400682015049/13781139750220=S+e/2\) | \(1.0449558074337872\) |
| \(\mu\) | \(5/2\) | 2.5 |
| split \(p_G\) | \(13481830255/64511729664=19660h\) | \(0.2089826195207935\) |
| mesh \(h\) | \(2742997/258046918656\) | \(0.00001062983822587963\) |
| radial lower endpoint | \(16430591763736936545249922448197799591/16161921199408696007503616565983946000\) | \(1.0166236774089747\) |
| radial upper endpoint | \(269644834091/258046918656=98303h\) | \(1.0449449871186451\) |

The Python derivation is in
[`prime_gap_186_certificate.py`, lines 200–261 and 342–425](https://github.com/openai/PrimeGaps186/blob/61340d0b74163003b32756bb16e91d9209a5e330/prime_gap_186_certificate.py#L200),
and its component schedule is at
[lines 264–339](https://github.com/openai/PrimeGaps186/blob/61340d0b74163003b32756bb16e91d9209a5e330/prime_gap_186_certificate.py#L264).
The independent Lean transcription is in
[`PrimeGaps186.lean`, lines 736–893](https://github.com/openai/PrimeGaps186/blob/61340d0b74163003b32756bb16e91d9209a5e330/PrimeGaps186.lean#L736),
with the radial/cap mask and cover at
[lines 927–997](https://github.com/openai/PrimeGaps186/blob/61340d0b74163003b32756bb16e91d9209a5e330/PrimeGaps186.lean#L927).
The numerical note gives the same data in equations (1.12), (1.23)–(1.26), and
(1.34)–(1.47), especially Lemmas 1.4–1.7.

Endpoint conventions are source-faithful: activation is strict
\(\xi<p\), the cap is inclusive \(p\leq z\), low and rank bins are
lower-open/upper-closed, and the high component counts \(p>p_G\).

## The stated bad event

The implementation uses the note's inclusive obstruction

\[
B=\left\{\max_{p\in\Pi,\ p>\xi}
  \left(\sum_{q\in\Pi:q\geq p}q+(\mu-1)p\right)>T\right\}.
\]

The inclusive prefix preserves the upstream convention when fragment sizes
coincide. A Poisson point process is simple almost surely, in which case this
agrees with \(S_{j-1}+\mu u_j>T\).

## Exact cover implemented

On the selected radial slice the radial mask is one, and the cap guard is
\(\Pi\subseteq(\xi,z]\). Let \(N(A)\) be the number of points in \(A\).
The script implements all 35 components of `physicalSourceCover`:

1. The 22 low bins \((\ell_j,u_j]\) use

   \[
   C_{L_j}=N((\ell_j,u_j])
   \exp\!\left(\theta_j\left[
      \sum_{q>\ell_j}q+(\mu-1)u_j-T
   \right]\right),
   \qquad \theta_j=\left\lceil\frac7{u_j}\right\rceil.
   \]

   The boundaries are exactly the `G2` list in §1.5.1 of the note:
   \((2^j\xi)_{j=0}^8\), followed by
   \(1/100,3/200,9/400,27/800\), then
   \(a_0,\ldots,a_7,(a_7+p_G)/2,p_G\), where
   \(a_j=(6/5)^j/20\).

2. Put \(q_0=T/(\mu+1)\). The 12 rank-two bins are the affine images of

   \[
   0,\frac1{16},\ldots,\frac8{16},\frac58,\frac34,\frac78,1
   \]

   from \(q_0\) to \(z\). For a bin \((a_j,b_j]\),

   \[
   C_{P_j}=\sum_{q\in\Pi\cap(a_j,b_j]}
     1_{N((q,\infty))=0}
     1_{N(((T-q)/\mu,q])\geq2}.
   \]

3. The third-factorial component is

   \[
   C_H={N((p_G,\infty))\choose3}.
   \]

The complete cover is \(C=\sum_j C_{L_j}+\sum_j C_{P_j}+C_H\).
These formulas are the point-configuration specialization of the Lean
definition cited above; the Python numerical code evaluates positive kernel
envelopes for their integrals rather than enumerating point configurations.

## What the source actually proves

PrimeGaps186 does not claim that this \(C\) covers the stated \(B\). On page 9
of the numerical note, Lemma 1.4 says that after the largest-fragment and
opposite-root caps have been imposed, every retained order-three failure is a
“nonlargest \(H_{5/2}\) failure.” Section 1.5.3 then begins with a
“second-largest witness,” and §1.5.4 handles rank at least three.

For the source-faithful event

\[
B_{\mathrm{nl}}=
\left\{\exists j\geq2:
S_{j-1}+\mu u_j>T\right\},
\]

the elementary cover proof is:

- If a violating witness \(p\leq p_G\), its unique low bin has a point count at
  least one and a strictly positive exponent, so its low component is at least
  one.
- If the violating witness above \(p_G\) is second-largest, call the largest
  point \(q\). Then \(q+\mu p>T\), hence
  \(p>(T-q)/\mu\), and \(q>T/(\mu+1)\); the rank bin containing \(q\) fires.
- If a violating witness above \(p_G\) has rank at least three, then at least
  three points exceed \(p_G\), so the factorial component is at least one.

This is the argument in equations (1.36)–(1.47). It proves
\(C\geq1_{B_{\mathrm{nl}}}\), not \(C\geq1_B\).

For the stated event, take \(\Pi=\{z\}\). Since \(\mu z>T\), \(B\) occurs at
index one. The point is above all low bins, rank two requires at least two
points, and \({1\choose3}=0\). Therefore \(C=0\).

## Boundary checks

The executable checks the empty state, points exactly at and just above
\(\xi\), the split and just above it, one point between \(T/\mu\) and \(z\), one
point exactly at \(z\), a barely violating two-point configuration, and three
points above the split. The two one-large configurations are the only
majorization failures in this list. The exact detailed values and nonzero
components are in the JSON result.

## Renewal recurrence and discretization

The script includes the requested recurrence

\[
Q_w(s,u)=\left(\frac\xi u\right)^n w(s)
+\int_\xi^{\min(u,(T-s)/\mu)}
 \left(\frac vu\right)^n Q_w(s+v,v)n\frac{dv}{v}.
\]

`renewal_dp` rewrites this as

\[
u^nQ_w(s,u)=\xi^n w(s)+
\int_\xi^{\min(u,(T-s)/\mu)}n v^{n-1}Q_w(s+v,v)\,dv.
\]

It uses linear grids in \(u\in[\xi,z]\) and \(s\in[0,T]\), descending dynamic
programming in \(s\), linear interpolation at \(s+v\), and trapezoidal
quadrature in \(v\). This is ready for 256-by-256 and 512-by-512 coarse runs,
including terminal weights \(1,s,s^2\), but running it would compute the
probability of the invalidly paired \(B\).

For \(B_{\mathrm{nl}}\), the first selected point must bypass the ordinary
barrier test. Its safe initial condition is instead

\[
Q_{\mathrm{nl},w}(0,z)=\left(\frac\xi z\right)^n w(0)
+\int_\xi^z\left(\frac vz\right)^n Q_w(v,v)n\frac{dv}{v}.
\]

That modified initial condition should be implemented and independently
checked only after confirming that \(B_{\mathrm{nl}}\) is the intended target.

There is a second semantic choice to make at that point: the full source cover
contains a radial mask, while the proposed PPP sampling law is unconditional.
Either the sampling law must be conditioned on the selected radial slice, or
both the bad event and cover must include the slice indicator. Treating the
mask as identically one while sampling configurations outside the slice does
not reproduce a source integral.

## MC/DP and stress-test status

| quantity | MC | DP-256 | DP-512 |
|---|---:|---:|---:|
| \(\mathbb P(B)\) | not run | not run | not run |
| \(\mathbb E[C]\) | not run | — | — |
| \(\mathbb E[C]/\mathbb P(B)\) | invalid | — | — |

| stress \(\eta\) | \(\Gamma(\eta)\) |
|---:|---:|
| 0 | invalid: Stage-1 majorization failure |
| 0.01 | not run |
| 0.02 | not run |
| 0.04 | not run |

No first-violation distribution was estimated. The deterministic failing
configuration has first violating index 1, fragment \(z\), inclusive prefix
mass \(z\), and barrier excess shown at the start of this document.

## Reproduction

Run:

```bash
python experiments/direct_source_barrier_viability.py
```

This regenerates
`experiments/direct_source_barrier_viability.json` using exact rational source
derivations and repeats all boundary assertions. The direction is neither a
GO nor a NO-GO on cover tightness. It is **INCONCLUSIVE** until the event and
radial sampling semantics are aligned with the source.

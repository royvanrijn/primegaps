# Exact reproduction of the 240 certificate

Status: `k=49`, `D=21` reproduced exactly; unchanged `k=48` baseline sealed.

This calculation independently reconstructs the variational certificate in
Section 5 of Julia Stadlmann's *Bounded gaps between primes*
([arXiv:2608.31126](https://arxiv.org/abs/2608.31126)). It uses the paper and
the repository's low-dimensional exact integrator as specifications. It does
not introduce or claim a new optimization result.

## Fixed conventions

All printed decimals are interpreted as exact rationals:

| parameter | exact value |
| --- | ---: |
| `epsilon` | `3/400` |
| `delta` | `7/250` |
| `A` | `(-3/400, 253/1000)` |
| total cap `U=A_1+epsilon` | `521/2000` |
| shared-coordinate cap `R=A_1-epsilon` | `491/2000` |
| `B_1,B_2` | `3/20` |
| `B_m`, `m>=3` | `17/100` |
| `c_1,c_2` | `0,0` |

Only six `B_m` entries need to be materialized: seven large coordinates have
sum strictly greater than `7 delta = 0.196 > 0.17`, so all statuses with seven
or more large coordinates are empty.

The 846-dimensional basis is represented as

    m_(2 lambda)(t) (U - sum(t))^b,  2|lambda| + b <= 21.

This is an invertible radial change of basis within the paper's degree-21
polynomial span. The even signature is `2 lambda`, not `lambda`. The paper's
introduction and theorem proof specify degree 21; the isolated `B_19` in the
detailed basis discussion is treated as a source typo.

Because `c_1=c_2=0`, the required Rayleigh quotient is exactly `k J(F)/I(F)`;
the general `K` form is irrelevant here.

## Exact result at k=49

For a fixed 846-term rational candidate, exact replay gives

    49 J(F) / I(F)
      = 1.0011632465949216560417861678682244509240906847502660897556997617934938...

The stored numerator of `49J-I` is positive, so this is a strict rational
inequality, not a floating-point confidence interval. The result artifact has
SHA-256

    2d6ff1239167c21ac79055deb7143d2bf52a24daafc82204d6529362ca0a81b5

and the fixed candidate has SHA-256

    c840f99232b6c821b1f63aa81e496d1e850a4f5b482e5822fbf537c06be90815.

The complete exact block files contain 2,714 product signatures each:

| artifact | SHA-256 |
| --- | --- |
| exact `I` groups | `f8af2b9831abbf3baa164b76bed91c14605fae39fb146bea209c28806c8c5be1` |
| `I` manifest | `3432d381a19b22b217b3c1243e7500d585d1189bf9d09909b951815ed0d0b039` |
| exact `J` groups | `3af2cdaf212b15d0d050cbf7f83afcb909446506f757271deca39214b228914f` |
| `J` manifest | `da9752e0fbcb8371e9226212bffa229211f882acbb8f077b22a800167a302963` |

## Independent checks

The high-dimensional calculation was admitted only after the following exact
checks passed:

- the status-density generating program matched explicit labeled-coordinate
  enumeration in four complete cases and the repository's reference support
  integrator in five final-integral cases;
- monomial-symmetric products, mixed radial terms, ordinary compressed `J`, and
  grouped `J` matched the independent low-dimensional slice engine in dimensions
  two and three, including repeated signatures and cancellation;
- candidate conversion from the numerical Jacobi basis matched all 846 stored
  rational terms exactly;
- finalization independently derives the expected `I` and `J` signature sets
  from the candidate and rejects missing, duplicate, foreign, or mixed-context
  rows;
- the current repository test suite passes (`65 passed, 2 skipped` in the
  non-Sage environment), including tracked hash, exact-identity, status-DP,
  low-dimensional `I`/`J`, grouped-`J`, cache, and frontier checks. The two
  skipped tests exercise optional Sage/FLINT kernels.

## Reproducibility boundary

The expensive calculation is explicit and checkpointed. Every group row is an
exact rational and carries a context hash. The launch manifest binds `k`, kind,
candidate and candidate-binding hashes, all support constants, the candidate's
canonical term hash, the independently derived 2,714-signature-set hash, Python
and GMP-family versions, and hashes of the verifier, status DP, runner,
provenance helper, and low-dimensional geometry kernel. A completed manifest
also binds the final records and summary hashes and cannot silently re-bless an
altered file.

The cheap finalizer never reruns the contraction. It verifies those bindings,
sums the recorded exact rationals, and checks the signs of `I` and `kJ-I`.
Research objects and canonical ledger records live under `.research/`. The
byte-identical rebuild entry points, both fixed candidates, and both result
artifacts are tracked under `reproduction/240/`.

## Exact baseline at k=48

The same support, `D=21` basis, 12-significant-decimal rationalization rule,
exact evaluator, and replay checks give

    48 J(F) / I(F)
      = 0.9969233513526357503888760066573328995217614432838426120218406063357052...

and therefore

    1 - 48 J(F) / I(F)
      = 0.0030766486473642496111239933426671004782385567161573879781593936642948...

The deficit numerator is the positive integer recorded in the result artifact;
equivalently, `48J-I` has its exact negative. Thus the unchanged setup does not
certify `k=48`, and no bound below 240 is claimed.

| artifact | SHA-256 |
| --- | --- |
| fixed candidate | `451301f90d6f5ded94f352a44cea935326fd9bb6dbac698fc45e829476f5479d` |
| exact `I` groups | `1a994c158ccc69bd9e70139da700f041e8fdddbd664bbc9c2a2272a781394a24` |
| `I` manifest | `97e0e92d1db64f056d7061b03653be26c6316fd3f2e4014c780b452e1ebfe917` |
| exact `J` groups | `e06d7cbd721c606adcd15ff101f399a9a05a857968a834b2b9dc08c03ede41a3` |
| `J` manifest | `0f472ff77effb143de54c855eb0d5c752dc148a115e63277158bd212e7b8e58b` |
| exact result | `fe186869bcbc51490bb7e2dbd499c80d7c5fa310fa2143fae0bb0e53548f02dc` |

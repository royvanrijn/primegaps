import PrimeGapsTheory.NumberTheory.Admissible

/-!
# The `DHL[k, 2]` interface needed for a bounded gap

This definition is independent of any particular sieve or distribution theorem.
It is the seam at which the future Stadlmann-style analytic argument connects to
the already-formalized finite/endgame argument.
-/

open Finset

namespace Gaps236

/-- `DHL[k, 2]`: every strictly increasing admissible `k`-tuple has infinitely
many translates containing at least two primes. -/
def DHL2 (k : ℕ) : Prop :=
  ∀ (h : Fin k → ℕ), StrictMono h →
    Finset.Admissible (Finset.image h Finset.univ) →
    {n : ℕ | 2 ≤ #{i : Fin k | (n + h i).Prime}}.Infinite

/-- The analytic statement sufficient for the explicit 48-tuple endgame. -/
abbrev DHL48 : Prop := DHL2 48

end Gaps236

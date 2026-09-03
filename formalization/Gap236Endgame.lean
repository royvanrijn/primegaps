import PrimeGaps236.Main

/-!
# Compatibility facade

The compiled formalization now lives under `PrimeGaps236`. This file preserves
the theorem name used by the original source-level adapter.
-/

namespace Gaps236

theorem thm_main_of_dhl48 (hDHL : DHL48) :
    ∃ᶠ m in Filter.atTop,
      Nat.nth Nat.Prime (m + 1) - Nat.nth Nat.Prime m ≤ 236 :=
  frequently_prime_gap_le_236_of_dhl48 hDHL

end Gaps236

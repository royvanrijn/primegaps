/-
Source-level Gap236 endgame adapter for AxiomMath/PrimeGapsLib.

Target upstream revision:
  1faa7b14e82ddebc2772dfb9153922f01b106477

This file is intentionally not claimed as compiled by this Python repository.
It isolates the finite/endgame theorem that should be checked inside the pinned
PrimeGapsLib environment.  The genuinely new analytic obligation is to prove
`DHL48`; this file does not assume that the current computational research has
already established it.
-/

import PrimeGapsTheory.Gap246.Endgame.Main

@[expose] public section

open Nat
open Finset MeasureTheory Filter
open scoped PrimeGaps BigOperators

namespace Gaps236

/-- An explicit admissible 48-tuple with diameter 236. -/
def H48 : Finset ℕ :=
  {0, 6, 8, 14, 18, 24, 26, 48, 50, 54, 56, 60,
   66, 68, 74, 78, 80, 84, 90, 96, 98, 104, 110, 116,
   120, 126, 134, 138, 144, 150, 158, 164, 168, 176, 180, 186,
   188, 194, 200, 204, 206, 210, 216, 224, 228, 230, 234, 236}

/-- `H48` has exactly 48 elements. -/
theorem card_H48 : #H48 = 48 := by
  set_option maxRecDepth 4000 in decide

/-- The diameter of `H48` is 236. -/
theorem diameter_H48 : H48.diameter = 236 := by
  set_option maxRecDepth 4000 in decide

/-- `H48` omits a residue class modulo every prime. -/
theorem admissible_H48 : H48.Admissible := by
  set_option maxRecDepth 4000 in decide

/-- Abstract `DHL[48,2]` interface, deliberately separated from the analytic proof. -/
def DHL48 : Prop :=
  ∀ (h : Fin 48 → ℕ), StrictMono h →
    Finset.Admissible (Finset.image h Finset.univ) →
    {n : ℕ | 2 ≤ #{i : Fin 48 | (n + h i).Prime}}.Infinite

/-- The finite/endgame implication: `DHL[48,2]` plus H48 gives infinitely many
consecutive prime gaps at most 236.  This is the `k=48` analogue of
`Gaps246.thm_main`; no distribution estimate occurs in this proof. -/
theorem thm_main_of_dhl48 (hDHL : DHL48) :
    ∃ᶠ m in Filter.atTop,
      Nat.nth Nat.Prime (m + 1) - Nat.nth Nat.Prime m ≤ 236 := by
  set h : Fin 48 → ℕ := ⇑(H48.orderEmbOfFin card_H48) with hh
  have himg : Finset.image h Finset.univ = H48 :=
    H48.image_orderEmbOfFin_univ card_H48
  have hmono : StrictMono h := (H48.orderEmbOfFin card_H48).strictMono
  have hinj : Function.Injective h := hmono.injective
  have hadm : Finset.Admissible (Finset.image h Finset.univ) := by
    rw [himg]
    exact admissible_H48
  have hInf : {n : ℕ | 2 ≤ #{i : Fin 48 | (n + h i).Prime}}.Infinite :=
    hDHL h hmono hadm
  have hfreq : ∃ᶠ n' in Filter.atTop,
      ContainsAtLeastPrimes n' (Finset.image h Finset.univ).diameter 2 := by
    haveI : Nonempty (Fin 48) := ⟨⟨0, by norm_num⟩⟩
    set H := Finset.image h Finset.univ with hHdef
    have hHne : H.Nonempty := Finset.univ_nonempty.image h
    set mn := H.min' hHne with hmndef
    set mx := H.max' hHne with hmxdef
    have hmnmx : mn ≤ mx := Finset.min'_le_max' H hHne
    have hdiam : H.diameter = mx - mn := Finset.diameter_eq_max_sub_min hHne
    rw [Filter.frequently_atTop]
    intro M
    obtain ⟨n, hnmem, hnM⟩ := hInf.exists_gt M
    have hrle : 2 ≤ #{i : Fin 48 | (n + h i).Prime} := hnmem
    refine ⟨n + mn, by omega, le_trans hrle ?_⟩
    apply Finset.card_le_card_of_injOn (fun i ↦ n + h i)
    · intro i hi
      rw [Finset.mem_coe, Finset.mem_filter] at hi
      have hiH : h i ∈ H := Finset.mem_image_of_mem h (Finset.mem_univ i)
      have h1 : mn ≤ h i := Finset.min'_le H (h i) hiH
      have h2 : h i ≤ mx := Finset.le_max' H (h i) hiH
      simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_Icc]
      refine ⟨⟨by omega, ?_⟩, hi.2⟩
      rw [hdiam]
      omega
    · intro i _ j _ hij
      have hij' : n + h i = n + h j := hij
      exact hinj (by omega)
  have hdiam236 : (Finset.image h Finset.univ).diameter = 236 := by
    rw [himg]
    exact diameter_H48
  rw [hdiam236] at hfreq
  have hfinal := frequently_prime_gap_le_of_frequently_interval 2 236 (by norm_num) hfreq
  simpa using hfinal

end Gaps236

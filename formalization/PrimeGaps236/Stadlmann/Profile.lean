import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

/-!
# Exact metadata for the active degree-27 candidate

This binds the formalization to the candidate file by exact rational data. It
does not claim that the candidate satisfies the variational inequality.
-/

namespace Gaps236.Stadlmann

/-- The support and basis metadata that identify a rational candidate. -/
structure RationalCandidateProfile where
  candidateSha256 : String
  k : ℕ
  degree : ℕ
  termCount : ℕ
  delta : ℚ
  epsilon : ℚ
  a : List ℚ
  b : List (List ℚ)
  deriving DecidableEq

/-- The candidate currently being certified by the D27 boundary calculation. -/
def d27Profile : RationalCandidateProfile where
  candidateSha256 := "d1cecf4db18f29da80cb2f5784ecfb9ee98ccc1221055b4a45ed3844427339ac"
  k := 48
  degree := 27
  termCount := 2526
  delta := 7 / 250
  epsilon := 17 / 2000
  a := [-17 / 2000, 23 / 100, 2029 / 8000]
  b :=
    [[9 / 50, 9 / 50, 1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5],
     [3 / 20, 3 / 20, 17 / 100, 17 / 100, 17 / 100, 17 / 100, 17 / 100]]

theorem d27Profile_dimensions :
    d27Profile.k = 48 ∧ d27Profile.degree = 27 ∧
      d27Profile.termCount = 2526 ∧ d27Profile.a.length = 3 ∧
      d27Profile.b.length = 2 ∧ ∀ row ∈ d27Profile.b, row.length = 7 := by
  norm_num [d27Profile]

theorem d27Profile_positive_scales :
    0 < d27Profile.delta ∧ 0 < d27Profile.epsilon := by
  norm_num [d27Profile]

/-- Outer simplex radius `A₂ + ε`, matching the exact-I computation. -/
theorem d27_outer_radius :
    (2029 / 8000 : ℚ) + d27Profile.epsilon = 2097 / 8000 := by
  norm_num [d27Profile]

/-- Shared-coordinate radius `A₂ - ε`, matching the exact-J computation. -/
theorem d27_shared_radius :
    (2029 / 8000 : ℚ) - d27Profile.epsilon = 1961 / 8000 := by
  norm_num [d27Profile]

end Gaps236.Stadlmann

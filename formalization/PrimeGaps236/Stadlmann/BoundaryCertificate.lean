import PrimeGaps236.Stadlmann.Profile

/-!
# Exact arithmetic interface for the pending D27 boundary calculation

The expensive calculation should emit a rational lower enclosure for the
boundary correction. Once its analytic soundness has been proved, the theorem
below reduces the variational crossing to exact ordered-field arithmetic.
-/

namespace Gaps236.Stadlmann

/-- Exact normalized unrestricted-simplex `I` for the active candidate. -/
def normalizedI : ℚ :=
  63945606120994044765287234573206684123131509281001157227985717635473047081518624668570734398828147978696297096807489435608181841239387053982750774851160251809680801952808048388513179357752902779545960764862965011252753123554242815834903615109 /
  63945890136159826402769406258137727436563468559511459803402157962419294299737887711948363721006480747397120000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

/-- Exact normalized `48 * J` on the unrestricted simplex. -/
def normalizedUnrestrictedKJ : ℚ :=
  30112843080690785377651165865151512458640386839403415031981979711577733819671231391455597548076799260393003525252224808137205175920288936486539368737445671867568192206800712905282910866890385556415130718752709689100415794067914129596601525923851415386271073025248650358732970950715278042939134098639394557847979290569875313619223375831628698915892938705296607860685340119026243722940010699586485567563670466974473 /
  30097304219928933581135646273585177918072770844534963073235536928598579212585939769205388411923254619329437108582510334382865581582016583537033576609678354811646284981296361946462006968959829617533592567126738787740860849068355411110733155951855622289839757022374108188167372800000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

/-- A small, independently checkable sanity property of the recorded exact
unrestricted values. -/
theorem known_unrestricted_crossing :
    normalizedI < normalizedUnrestrictedKJ := by
  norm_num [normalizedI, normalizedUnrestrictedKJ]

/-- The arithmetic part of a boundary certificate. `lower` is intended to be
a rigorous rational lower enclosure for the normalized legal-minus-
unrestricted `48 * J` correction. -/
structure BoundaryArithmeticCertificate where
  lower : ℚ
  passes : normalizedI < normalizedUnrestrictedKJ + lower

/-- If the computed rational bound is a sound lower enclosure for the actual
correction, and the known unrestricted quantities bound the actual legal
integrals in the stated directions, then the desired strict crossing follows.

This theorem cleanly separates cheap exact arithmetic from the future analytic
proof that the emitted enclosure really bounds the shaped-support integral. -/
theorem variational_crossing_of_boundary_certificate
    (cert : BoundaryArithmeticCertificate)
    {legalI legalKJ correction : ℝ}
    (hI : legalI ≤ (normalizedI : ℝ))
    (hCorrection : (cert.lower : ℝ) ≤ correction)
    (hJ : (normalizedUnrestrictedKJ : ℝ) + correction ≤ legalKJ) :
    legalI < legalKJ := by
  have hpass : (normalizedI : ℝ) <
      (normalizedUnrestrictedKJ : ℝ) + cert.lower := by
    exact_mod_cast cert.passes
  linarith

end Gaps236.Stadlmann

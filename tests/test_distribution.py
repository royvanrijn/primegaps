from fractions import Fraction

from primegaps.distribution import Minorant, RegionCell, cells_from_support, is_certified
from primegaps.support import stadlmann_240_parameters


PUBLISHED_MINORANT = Minorant("0.38", "0.4", "0.4")


def test_bv_pair_has_short_certificate():
    a = RegionCell("0.24", 2, "0.12", "0.028", support_max="0.253")
    certificate = is_certified(a, a, PUBLISHED_MINORANT)
    assert certificate
    assert certificate.theorem == "Bombieri–Vinogradov + Proposition 2"
    assert certificate.modulus_exponent_bound == Fraction(12, 25)
    assert not certificate.partitions


def test_published_non_bv_cell_pair_gets_all_type_witnesses():
    cell = RegionCell("0.253", 3, "0.17", "0.028")
    certificate = is_certified(cell, cell, PUBLISHED_MINORANT)
    assert certificate, certificate.as_dict()
    assert certificate.minorant_kind == "prime-indicator"
    assert [w.condition[0] for w in certificate.partitions] == ["A", "B", "C", "D", "E"]
    type_iic = certificate.partitions[3]
    assert type_iic.worst_bin_sums[0] <= type_iic.capacities[0]
    assert type_iic.worst_bin_sums[1] <= type_iic.capacities[1]


def test_all_published_cell_pairs_are_certified_or_vacuously_empty():
    cells = cells_from_support(stadlmann_240_parameters())
    for a in cells:
        for b in cells:
            certificate = is_certified(a, b, PUBLISHED_MINORANT)
            assert certificate, (a.label, b.label, certificate.as_dict())


def test_failed_global_hypothesis_names_the_failed_check():
    cell = RegionCell("0.28", 3, "0.17", "0.028")
    certificate = is_certified(cell, cell, PUBLISHED_MINORANT)
    assert not certificate
    assert "global" in certificate.reason
    assert any(not check.passed for check in certificate.checks if check.name.startswith("P3."))


def test_failed_partition_is_not_misstated_as_a_counterexample():
    cell = RegionCell("0.253", 2, "0.24", "0.028")
    certificate = is_certified(cell, cell, PUBLISHED_MINORANT)
    assert not certificate
    assert "No universal partition witness" in certificate.reason
    assert "conservative" in certificate.caveats[0]


def test_proposition_2_minorant_parameters_are_checked_first():
    cell = RegionCell("0.24", 0, 0, "0.028")
    certificate = is_certified(cell, cell, Minorant("0.1", "0.1", "0.1"))
    assert not certificate
    assert "Proposition 2" in certificate.reason


def test_decimal_inputs_produce_exact_boundary_decisions():
    left = RegionCell("0.25", 0, 0, "0.028", support_max="0.253")
    right = RegionCell("0.25", 0, 0, "0.028", support_max="0.253")
    certificate = is_certified(left, right, PUBLISHED_MINORANT)
    assert certificate
    assert certificate.modulus_exponent_bound == Fraction(1, 2)

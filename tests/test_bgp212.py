from fractions import Fraction

from primegaps.bgp212 import (
    modulus_classes,
    packing_problem,
    parameters,
    recomputed_table6_rows,
    reported_table6_rows,
    section_9_stale_datum_discrepancy,
    table6_source_discrepancies,
)


Q = Fraction


def test_table3_exact_datum_and_support_invariants():
    p = parameters()
    assert p.k == 45
    assert p.omega == Q(7, 1000)
    assert p.a0 == -p.support_epsilon == Q(-1, 125)
    assert p.a1 == Q(257, 1000)
    assert p.total_cap == Q(53, 200)
    assert p.marginal_cap == Q(249, 1000)
    assert p.delta == Q(41, 2500)
    assert len(p.rough_caps) == 60
    assert p.rough_caps[:3] == (Q(777, 5000), Q(794, 5000), Q(875, 5000))
    assert p.rough_caps[9:] == (Q(1081, 5000),) * 51
    assert p.maximum_nonempty_rough_count == 13

    for left, right in zip(p.rough_caps, p.rough_caps[1:]):
        assert left <= right <= left + p.delta


def test_five_modulus_classes_transcribe_factor_windows_exactly():
    p = parameters()
    classes = modulus_classes()
    assert tuple(item.identifier for item in classes) == (
        "D_I",
        "D_IIa",
        "D_IIb",
        "D_IIc",
        "D_III",
    )
    values = {
        "omega": p.omega,
        "gamma": Q(2, 5),
        "delta": p.delta,
        "epsilon": Q(1, 10**6),
        "theta": Q(257, 500),
        "rho": Q(2, 5),
    }
    for modulus_class in classes:
        assert modulus_class.modulus_upper_exponent.evaluate(values) == Q(257, 500)
        for case in modulus_class.cases:
            for interval in case.factor_intervals:
                assert interval.width(values) == p.delta

    type_iic = classes[3]
    intervals = type_iic.cases[0].factor_intervals
    assert [(item.factor, item.divides) for item in intervals] == [
        ("r", "d"),
        ("u", "d/r"),
        ("d1", "r"),
    ]
    assert all(wall.holds(values) for wall in type_iic.analytic_walls)


def test_continuous_packing_problem_has_all_455_quantified_roots():
    p = parameters()
    problem = packing_problem(p)
    assert problem.maximum_count == 13
    assert len(problem.ordered_positive_count_pairs) == 91
    assert tuple(condition.identifier for condition in problem.conditions) == (
        "A",
        "B",
        "C",
        "D",
        "E",
    )
    assert problem.expected_root_count == problem.reported_successful_roots == 455

    # Condition A contains two distinct existential partitions.  Condition D
    # is uniform over a rectangle and its partition may vary over that domain.
    assert len(problem.conditions[0].partition_requirements) == 2
    type_iic = problem.conditions[3]
    domains = {item.variable: (item.lower, item.upper) for item in type_iic.parameter_domains}
    assert domains["gamma"] == (
        p.xi2 - p.analytic_epsilon,
        Q(1, 3) + 8 * p.omega + Q(7, 3) * p.delta + 3 * p.analytic_epsilon,
    )
    assert domains["omega0"] == (Q(0), p.omega)
    capacities = type_iic.partition_requirements[0].capacities
    low_endpoint = {"gamma": domains["gamma"][0], "omega0": Q(0)}
    assert capacities[0].evaluate(low_endpoint) == p.xi2 - 2 * p.delta - 2 * p.analytic_epsilon
    assert capacities[1].evaluate(low_endpoint) == Q(1, 2) - p.xi2
    assert capacities[2].evaluate(low_endpoint) == p.delta - p.analytic_epsilon
    assert capacities[3].evaluate(low_endpoint) == 0

    assert problem.rough_profile(13, 13).left_cap == Q(1081, 5000)


def test_table6_replay_exposes_the_single_numeric_mismatch():
    reported = reported_table6_rows()
    recomputed = recomputed_table6_rows()
    assert len(reported) == len(recomputed) == 21
    assert all(row.slack >= 0 for row in reported)
    assert all(row.slack > 0 for row in reported if row.strict)
    assert all(row.slack >= 0 for row in recomputed)

    discrepancies = table6_source_discrepancies()
    assert discrepancies == (
        {
            "identifier": "global.II.range",
            "reported_left": "0",
            "reported_right": "69599997/2000000000",
            "recomputed_left": "0",
            "recomputed_right": "347999991/10000000000",
            "right_difference": "3/5000000000",
        },
    )


def test_section9_stale_components_cancel_in_the_combined_wall():
    discrepancy = section_9_stale_datum_discrepancy()
    assert discrepancy["table3_A1"] != discrepancy["section9_printed_A1"]
    assert discrepancy["table3_delta"] != discrepancy["section9_printed_delta"]
    assert discrepancy["table3_combination"] == discrepancy["section9_printed_combination"] == "3937/5000"

import importlib.util
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/sweep_p3ii_delta_frontier.py"
SPEC = importlib.util.spec_from_file_location("p3ii_delta_frontier", SCRIPT)
frontier = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(frontier)


def test_exact_constraint_transition_is_p3ii_delta_to_range():
    bounds = frontier.exact_constraint_bounds()
    assert bounds["P3.II.delta.branch2"] == Fraction(1265833333, 5000000000)
    assert bounds["P3.II.range"] == Fraction(913600001, 3600000000)
    later = [
        value for name, value in bounds.items()
        if not name.startswith("P3.II.delta")
    ]
    assert bounds["P3.II.range"] == min(later)


def test_frontier_grid_spans_both_exact_constraint_endpoints():
    points = frontier.curve_grid(16)
    assert points[0] == frontier.Decimal("0.253")
    assert frontier.Decimal("0.2537") in points
    assert float(points[-1]) == float(frontier.exact_constraint_bounds()["P3.II.range"])


def test_shifted_stratum_volume_weight_is_translated_simplex_ratio():
    dimension = 48
    large_count = 2
    radius = 0.262
    delta = 0.028
    expected = math.comb(dimension, large_count) * (
        (radius - large_count * delta) / radius
    ) ** dimension
    assert frontier.shifted_stratum_volume_weight(
        dimension, large_count, radius, delta
    ) == expected
    assert frontier.shifted_stratum_volume_weight(4, 5, radius, delta) == 0.0

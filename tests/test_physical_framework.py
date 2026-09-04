from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys

from primegaps.physical import (
    exact_group_failure,
    grouped_five_halves_failure,
    inclusive_obstruction,
    nonlargest_five_halves_obstruction,
    primegaps186_source_data,
)


ROOT = Path(__file__).resolve().parents[1]


def test_exact_source_inventory_and_inclusive_tie_convention():
    old, new, groups, mesh = primegaps186_source_data()
    assert (len(old), len(new)) == (29, 43)
    assert [len(group.rows) for group in groups] == [19, 4, 15]
    assert mesh == Fraction(2742997, 2624989) / 98304

    points = (Fraction(1, 10), Fraction(1, 10), Fraction(1, 20))
    obstruction = inclusive_obstruction(points, Fraction(1, 100), lambda p: p)
    assert obstruction == Fraction(3, 10)
    # The largest of two tied values is excluded, while its second occurrence
    # is a valid nonlargest witness with the full inclusive prefix.
    assert nonlargest_five_halves_obstruction(
        points[:2], Fraction(1, 100)
    ) == Fraction(7, 20)


def test_recorded_strict_h25_overcovers_are_exact():
    payload = json.loads(
        (ROOT / "experiments" / "physical_factorization_oracle.json").read_text()
    )
    _, _, groups, _ = primegaps186_source_data()
    by_name = {group.name: group for group in groups}
    for row in payload["groups"]:
        group = by_name[row["group"]]
        witness = row["strict_overcover_witness"]
        fragments = tuple(Fraction(value) for value in witness["fragments"])
        total = Fraction(witness["total_mass"])
        assert not exact_group_failure(group, fragments, total)
        assert grouped_five_halves_failure(group, fragments, total)


def test_physical_restored_operator_record_replays():
    completed = subprocess.run(
        [sys.executable, ROOT / "scripts" / "check_physical_restored_operator.py"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    replay = json.loads(completed.stdout)
    assert replay["status"] == "checked-exploratory-record"
    assert replay["dimensions"]["39"]["restored_score_upper_screen"] < 0.995
    assert replay["dimensions"]["39"]["full_face_score"] > 1.019


def test_physical_factorization_oracle_record_replays():
    completed = subprocess.run(
        [sys.executable, ROOT / "scripts" / "check_physical_factorization_oracle.py"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    replay = json.loads(completed.stdout)
    assert replay["status"] == "checked-exact-predicate-record"
    assert all(row["false_negatives"] == 0 for row in replay["groups"].values())

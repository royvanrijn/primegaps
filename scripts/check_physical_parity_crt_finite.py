#!/usr/bin/env python3
"""Cheap replay of the recorded reduced CRT-parity experiment."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULT = ROOT / "experiments/physical_parity_crt_finite.json"


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(right), 1e-300)


def replay(path: Path = DEFAULT_RESULT) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["schema"] == "primegaps.parity-crt-finite.v1"
    assert payload["status"] == "exploratory-reduced-finite-model-not-an-asymptotic-or-theorem"

    physical = payload["inputs"]["physical_input"]
    vector = payload["inputs"]["vector_input"]
    assert file_sha256(ROOT / physical["path"]) == physical["sha256"]
    assert file_sha256(ROOT / vector["path"]) == vector["sha256"]

    scales = payload["scales"]
    assert [item["sector_enumeration"]["x"] for item in scales] == [
        8_000_000,
        16_000_000,
        32_000_000,
        64_000_000,
    ]

    rows = []
    maximum_mesh_sensitivity = 0.0
    for item in scales:
        sector = item["sector_enumeration"]
        model = item["common_reduced_model"]
        assert sector["degree_two_finite_carrier"]
        assert model["actual_maximum_prime_exponent"] <= model["one_prime_exponent_envelope"]
        assert model["actual_full_atom_modulus_exponent"] <= model["pair_modulus_exponent_envelope"]
        assert model["compatible_pair_terms_after_support"] == 8455
        assert model["maximum_generated_modulus"] == 5005

        trials = item["trials"]
        primary = trials["k39_n8192"]
        control = trials["k39_n4096"]
        coefficients = primary["coefficient_metrics"]
        parity = primary["parity_cancellation"]
        before = parity["before_crt_aggregation"]
        after = parity["after_crt_aggregation"]
        coefficient_svd = primary["singular_value_decay"]["coefficient_table"]
        projected_svd = primary["singular_value_decay"]["projected_error_table"]

        assert coefficients["coefficient_compiler_identity"]["maximum_absolute_error"] <= 3e-16
        assert abs(coefficients["sum_q_a_abs_c_squared"] - 0.41529689489688126) <= 2e-15
        assert abs(coefficients["crt_collision_cancellation_ratio"] - 3.6387539702270675) <= 2e-14
        assert parity["E_prime"] * parity["E_semiprime"] < 0
        assert abs(parity["exact_scalar_ratio"] - 1.0) <= 2e-15
        assert after["cancellation_ratio"] >= 8.0
        assert before["cancellation_ratio"] > after["cancellation_ratio"]
        assert parity["gross_reduction_from_crt_aggregation"] >= 2.8
        assert after["within_block_prime_semiprime_ratio"] < 1.05
        assert after["across_block_ratio"] >= 8.0
        assert coefficient_svd["cumulative_frobenius_energy"]["8"] >= 0.996
        assert projected_svd["cumulative_frobenius_energy"]["8"] >= 0.86

        sensitivity = relative_difference(
            control["parity_cancellation"]["after_crt_aggregation"]["cancellation_ratio"],
            after["cancellation_ratio"],
        )
        maximum_mesh_sensitivity = max(maximum_mesh_sensitivity, sensitivity)
        rows.append(
            {
                "x": sector["x"],
                "sum_q_a_abs_c_squared": coefficients["sum_q_a_abs_c_squared"],
                "crt_collision_cancellation_ratio": coefficients["crt_collision_cancellation_ratio"],
                "literal_scalar_ratio": parity["exact_scalar_ratio"],
                "block_ratio_before_crt": before["cancellation_ratio"],
                "block_ratio_after_crt": after["cancellation_ratio"],
                "within_block_prime_semiprime_ratio": after["within_block_prime_semiprime_ratio"],
                "across_crt_block_ratio": after["across_block_ratio"],
                "projected_top8_energy": projected_svd["cumulative_frobenius_energy"]["8"],
            }
        )

    assert maximum_mesh_sensitivity < 0.05
    return {
        "schema": "primegaps.parity-crt-finite-replay.v1",
        "result_path": str(path),
        "result_sha256": file_sha256(path),
        "status": "checked",
        "interpretation": {
            "literal_scalar_gate": "NO-GO in this model: ratio equals one because sector errors reinforce",
            "blockwise_global_operator_gate": "GO for further study: post-CRT ratio is at least eight at every scale",
            "cancellation_location": "across CRT-coloured blocks, not within aligned prime/semiprime blocks",
        },
        "maximum_relative_mesh_vector_sensitivity": maximum_mesh_sensitivity,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    args = parser.parse_args()
    print(json.dumps(replay(args.result), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

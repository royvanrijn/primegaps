from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "export_d27_lean_certificate.py"
SPEC = importlib.util.spec_from_file_location("export_d27_lean_certificate", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def eligible_result() -> dict[str, object]:
    return {
        "schema": "primegaps-boundary-J-arb-result-v1",
        "checkpoint_complete": True,
        "certification_eligible": True,
        "input_sha256": "checkpoint-hash",
        "manifest_sha256": "manifest-hash",
        "exact_gate": {
            "i_upper_sha256": MODULE.I_UPPER_SHA256,
            "unrestricted_j_sha256": MODULE.UNRESTRICTED_J_SHA256,
            "certified_strictly_above_one": True,
            "normalized_correction_lower": ["0", "1"],
            "required_normalized_correction": ["-1", "2"],
            "normalized_legal_kJ_lower": ["5", "4"],
            "normalized_i_upper": ["1", "1"],
        },
    }


def test_render_emits_kernel_checked_certificate() -> None:
    rendered = MODULE.render(eligible_result())
    assert "lower := 0 / 1" in rendered
    assert "passes := by norm_num [normalizedI, normalizedUnrestrictedKJ]" in rendered
    assert "checkpoint-hash" in rendered


def test_render_rejects_non_certifying_result() -> None:
    result = eligible_result()
    result["exact_gate"]["certified_strictly_above_one"] = False
    with pytest.raises(ValueError, match="does not certify"):
        MODULE.render(result)

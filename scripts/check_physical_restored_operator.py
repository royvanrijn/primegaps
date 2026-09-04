#!/usr/bin/env python3
"""Cheap replay of the signed physical restoration operator experiment."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import isclose
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RECORD = ROOT / "experiments" / "physical_restored_operator.json"


def _file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _fit(meshes: list[int], values: list[float], start: int) -> float:
    coefficients = np.polyfit(
        1 / np.asarray(meshes[start:], dtype=float),
        np.asarray(values[start:], dtype=float),
        2,
    )
    return float(np.polyval(coefficients, 1 / 98_304))


def check(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "primegaps.physical-restored-operator.v1":
        raise ValueError("unexpected restored-operator schema")
    implementation = ROOT / "experiments" / "physical_restored_operator.py"
    if payload["implementation_sha256"] != _file_hash(implementation):
        raise ArithmeticError("restored-operator implementation hash changed")
    for path_key, hash_key in (
        ("physical_input", "physical_input_sha256"),
        ("matrix_builder", "matrix_builder_sha256"),
    ):
        source = ROOT / payload["inputs"][path_key]
        if payload["inputs"][hash_key] != _file_hash(source):
            raise ArithmeticError(f"restored-operator input hash changed: {source}")
    meshes = payload["meshes"]
    if meshes != [1024, 2048, 4096, 8192]:
        raise ArithmeticError("unexpected direct-convolution mesh ladder")
    operator = payload["operator"]
    if operator["basis_dimension"] != 77:
        raise ArithmeticError("restored operator does not have 77 dimensions")
    a, b = Fraction(operator["a"]), Fraction(operator["b"])
    if a + b != Fraction(operator["a_plus_b"]) or b >= 0:
        raise ArithmeticError("signed hybrid coefficients changed")

    checked = {}
    for row in payload["dimensions"]:
        dimension = int(row["dimension"])
        fixed = [
            float(row["finest_restored_vector_scores"][str(mesh)])
            for mesh in meshes
        ]
        prediction = _fit(meshes, fixed, 1)
        recorded = float(row["restored_quadratic_inverse_mesh_fit"]["production"])
        if not isclose(prediction, recorded, abs_tol=2e-14):
            raise ArithmeticError(f"k={dimension} restored fit does not replay")
        if not isclose(1 - prediction, row["production_deficit_from_one"], abs_tol=2e-14):
            raise ArithmeticError(f"k={dimension} deficit does not replay")
        checked[str(dimension)] = {
            "restored_score_upper_screen": prediction,
            "deficit": 1 - prediction,
            "full_face_score": row["full_face_optimized_fit"]["production"],
        }

    if not checked["39"]["restored_score_upper_screen"] < 0.995 < 1:
        raise ArithmeticError("k=39 signed-restoration classification changed")
    if not checked["39"]["full_face_score"] > 1.019:
        raise ArithmeticError("k=39 full-face control no longer crosses one")
    if not checked["38"]["restored_score_upper_screen"] < 0.99:
        raise ArithmeticError("k=38 signed-restoration classification changed")
    calibration = payload["calibration"]
    if abs(calibration["difference"]) > 1e-5:
        raise ArithmeticError("k=40 published calibration missed its tolerance")
    return {
        "schema": "primegaps.physical-restored-operator-replay.v1",
        "status": "checked-exploratory-record",
        "dimensions": checked,
        "calibration_difference": calibration["difference"],
        "classification": "k=39 and k=38 below one before source losses",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("record", type=Path, nargs="?", default=DEFAULT_RECORD)
    args = parser.parse_args()
    print(json.dumps(check(args.record), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Count exact J geometry pieces touched by the B-boundary correction."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = next(parent for parent in (HERE, *HERE.parents) if (parent / "pyproject.toml").is_file())
VERIFIER_PATH = (
    ROOT / "reproduction/240/independent-reproducer/exact_symmetric_verifier.py"
)
spec = importlib.util.spec_from_file_location(
    "control_variate_geometry_verifier", VERIFIER_PATH
)
verifier = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = verifier
spec.loader.exec_module(verifier)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=49)
    args = parser.parse_args()
    common_dimension = args.k - 1
    maximum_offset = int(verifier.R // verifier.DELTA)
    counts = {
        "status_orientations": 0,
        "unaffected_status_orientations": 0,
        "cells": 0,
        "boundary_cells": 0,
    }
    by_large = {}
    for large in range(min(common_dimension, len(verifier.B)) + 1):
        subtotal = {key: 0 for key in counts}
        for shifted in range(maximum_offset - large + 1):
            total_offset = (large + shifted) * verifier.DELTA
            large_offset = large * verifier.DELTA
            for left_large in (False, True):
                if large + int(left_large) > len(verifier.B):
                    continue
                left_limit = verifier._support_limit(large, left_large)
                for right_large in (False, True):
                    if large + int(right_large) > len(verifier.B):
                        continue
                    right_limit = verifier._support_limit(large, right_large)
                    geometry_specs = (
                        verifier.RadialSlice(
                            0, 0, left_large, support_limit=left_limit
                        ),
                        verifier.RadialSlice(
                            0, 0, right_large, support_limit=right_limit
                        ),
                    )
                    kind, cells = verifier._slice_geometry(
                        large > 0,
                        common_dimension > large,
                        total_offset,
                        large_offset,
                        geometry_specs,
                    )
                    iterable = cells if kind != "point" else (cells[0],)
                    orientation_affected = False
                    cell_count = 0
                    affected_count = 0
                    for cell in iterable:
                        if kind == "polygons":
                            _polygon, sample = cell
                        elif kind in ("xintervals", "zintervals"):
                            _start, _end, sample = cell
                        elif kind == "point":
                            sample = cell
                        else:
                            continue
                        cell_count += 1
                        legal_left = verifier._slice_polynomial(
                            geometry_specs[0], total_offset, large_offset, sample
                        )
                        legal_right = verifier._slice_polynomial(
                            geometry_specs[1], total_offset, large_offset, sample
                        )
                        full_left = verifier._slice_polynomial(
                            verifier.RadialSlice(0, 0, left_large),
                            total_offset,
                            large_offset,
                            sample,
                        )
                        full_right = verifier._slice_polynomial(
                            verifier.RadialSlice(0, 0, right_large),
                            total_offset,
                            large_offset,
                            sample,
                        )
                        affected = (
                            legal_left != full_left or legal_right != full_right
                        )
                        affected_count += int(affected)
                        orientation_affected |= affected
                    subtotal["status_orientations"] += 1
                    subtotal["unaffected_status_orientations"] += int(
                        not orientation_affected
                    )
                    subtotal["cells"] += cell_count
                    subtotal["boundary_cells"] += affected_count
        by_large[str(large)] = subtotal
        for key, value in subtotal.items():
            counts[key] += value
    payload = {
        "schema": "primegaps-J-control-variate-geometry-v1",
        "k": args.k,
        "support": {
            "delta": str(verifier.DELTA),
            "U": str(verifier.U),
            "R": str(verifier.R),
            "B": [str(value) for value in verifier.B],
        },
        "counts": counts,
        "boundary_cell_fraction": counts["boundary_cells"] / counts["cells"],
        "affected_orientation_fraction": 1 - (
            counts["unaffected_status_orientations"]
            / counts["status_orientations"]
        ),
        "by_large": by_large,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Cheap replay of a recorded sparse-SOS screen without QMC or optimization."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import tempfile

import numpy as np


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("summary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.summary.read_text())
    if payload.get("schema") != "primegaps-sparse-sos-screen-v1":
        raise ValueError("unknown summary schema")
    replay = []
    for result in payload["results"]:
        artifact = Path(result["artifact"])
        if file_hash(artifact) != result["artifact_sha256"]:
            raise ValueError(f"artifact hash mismatch: {artifact}")
        with np.load(artifact, allow_pickle=False) as arrays:
            matrix_i = arrays["matrix_i"]
            objective = arrays["matrix_objective"]
            vector = arrays["rank_one_vector"]
            q_matrix = arrays["psd_matrix"]
            factor = arrays["sos_factor"]
            groups = arrays["group_of_component"]
            allowed = arrays["allowed_mask"].astype(bool)
        rank_one_value = float(vector @ objective @ vector / (vector @ matrix_i @ vector))
        psd_value = float(np.sum(objective * q_matrix) / np.sum(matrix_i * q_matrix))
        forbidden = [
            abs(float(q_matrix[left, right]))
            for left in range(len(groups))
            for right in range(left + 1, len(groups))
            if not allowed[int(groups[left]), int(groups[right])]
        ]
        factor_error = float(
            np.linalg.norm(q_matrix - factor.T @ factor) / max(np.linalg.norm(q_matrix), 1e-300)
        )
        if not np.isclose(rank_one_value, result["rank_one"]["value"], rtol=2e-8):
            raise ValueError("rank-one replay mismatch")
        if not np.isclose(psd_value, result["psd"]["value"], rtol=2e-8):
            raise ValueError("PSD replay mismatch")
        if max(forbidden, default=0.0) > 1e-7:
            raise ValueError("forbidden PSD coefficient is nonzero")
        if factor_error > 1e-6:
            raise ValueError("stored SOS factor does not reconstruct Q")
        replay.append(
            {
                "k": result["k"],
                "rank_one_value": rank_one_value,
                "psd_value": psd_value,
                "relative_advantage": psd_value / rank_one_value - 1.0,
                "psd_rank": result["psd"]["rank"],
                "forbidden_max_abs": max(forbidden, default=0.0),
                "factor_relative_error": factor_error,
            }
        )
    output = {"schema": "primegaps-sparse-sos-replay-v1", "results": replay}
    rendered = json.dumps(output, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            "w", dir=args.output.parent, prefix=args.output.name, delete=False
        ) as stream:
            stream.write(rendered)
            temporary = Path(stream.name)
        temporary.replace(args.output)
    print(rendered, end="")


if __name__ == "__main__":
    main()

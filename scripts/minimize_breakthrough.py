#!/usr/bin/env python3
"""Cheap replay of a scored support set under simultaneous theorem slacks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile

from primegaps.breakthrough import minimum_breakthrough
from primegaps.distribution import Minorant
from primegaps.shadow_prices import ScoredSupport
from primegaps.support import SupportParameters


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        temporary = Path(stream.name)
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Minimize weighted simultaneous analytic slack over scored supports"
    )
    parser.add_argument("scores", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--weights", type=Path)
    parser.add_argument("--target", type=float, default=1.0)
    parser.add_argument("--score-z", type=float, default=0.0)
    args = parser.parse_args()

    payload = json.loads(args.scores.read_text())
    minorant_data = payload["minorant"]
    minorant = Minorant(
        minorant_data["xi1"], minorant_data["xi2"], minorant_data["xi3"]
    )
    candidates = []
    for row in payload["candidates"]:
        support_data = row["support"]
        candidates.append(
            ScoredSupport(
                row["candidate_id"],
                SupportParameters(
                    support_data["delta"],
                    support_data["epsilon"],
                    tuple(support_data["A"]),
                    tuple(tuple(values) for values in support_data["B"]),
                ),
                row["score"],
                row.get("score_standard_error", 0.0),
            )
        )
    weights = json.loads(args.weights.read_text()) if args.weights else None
    result = minimum_breakthrough(
        candidates,
        minorant,
        weights,
        target_score=args.target,
        score_standard_error_multiplier=args.score_z,
    ).as_dict()
    result["score_kind"] = payload.get("score_kind")
    result["uncertainty_kind"] = payload.get("uncertainty_kind")
    atomic_json(args.output, result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

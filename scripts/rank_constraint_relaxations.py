#!/usr/bin/env python3
"""Replay a scored support grid under one-at-a-time theorem relaxations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from primegaps.distribution import Minorant
from primegaps.shadow_prices import ScoredSupport, rank_constraint_relaxations
from primegaps.support import SupportParameters


def _support(payload: dict[str, object]) -> SupportParameters:
    return SupportParameters(
        delta=float(payload["delta"]),
        epsilon=float(payload["epsilon"]),
        A=tuple(float(value) for value in payload["A"]),
        B=tuple(tuple(float(value) for value in row) for row in payload["B"]),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rank finite score gains from relaxing analytic constraints one at a time"
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--constraints",
        help="comma-separated stable constraint IDs; default is every registered constraint",
    )
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    if payload.get("schema") != "primegaps.scored-support-grid.v1":
        raise ValueError("input schema must be primegaps.scored-support-grid.v1")
    raw_minorant = payload["minorant"]
    minorant = Minorant(raw_minorant["xi1"], raw_minorant["xi2"], raw_minorant["xi3"])
    candidates = tuple(
        ScoredSupport(
            candidate_id=item["candidate_id"],
            support=_support(item["support"]),
            score=float(item["score"]),
            score_standard_error=float(item.get("score_standard_error", 0.0)),
        )
        for item in payload["candidates"]
    )
    kwargs = {}
    if args.constraints:
        kwargs["constraint_ids"] = tuple(
            value.strip() for value in args.constraints.split(",") if value.strip()
        )
    report = rank_constraint_relaxations(candidates, minorant, **kwargs).as_dict()
    report["score_kind"] = payload.get("score_kind", "unspecified")
    report["uncertainty_kind"] = payload.get("uncertainty_kind", "unspecified")
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()

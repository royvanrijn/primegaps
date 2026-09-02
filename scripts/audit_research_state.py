#!/usr/bin/env python3
"""Cheap structural replay of the local append-only research store."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Iterator


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _record_hash(record: dict[str, Any]) -> str:
    body = {key: value for key, value in record.items() if key != "record_sha256"}
    canonical = json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    return _sha256_bytes(canonical)


def _walk(value: Any) -> Iterator[tuple[str, Any]]:
    if isinstance(value, dict):
        for key, child in value.items():
            yield key, child
            yield from _walk(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk(child)


def audit(root: Path) -> dict[str, Any]:
    research = root / ".research"
    ledger_files = sorted((research / "ledger").glob("*/*.jsonl"))
    records: list[tuple[Path, int, dict[str, Any]]] = []
    parse_issues: list[dict[str, Any]] = []

    for ledger in ledger_files:
        for number, line in enumerate(ledger.read_text().splitlines(), 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
                if not isinstance(record, dict):
                    raise TypeError("record is not a JSON object")
                records.append((ledger, number, record))
            except (json.JSONDecodeError, TypeError) as exc:
                parse_issues.append(
                    {
                        "location": f"{ledger.relative_to(root)}:{number}",
                        "issue": str(exc),
                    }
                )

    known_hashes = {
        record.get("record_sha256")
        for _, _, record in records
        if isinstance(record.get("record_sha256"), str)
    }
    counts = Counter(record.get("record_sha256") for _, _, record in records)
    record_hash_issues: list[dict[str, Any]] = []
    dependency_issues: list[dict[str, Any]] = []
    object_hashes: set[str] = set()

    for ledger, number, record in records:
        location = f"{ledger.relative_to(root)}:{number}"
        claimed = record.get("record_sha256")
        computed = _record_hash(record)
        if claimed != computed:
            record_hash_issues.append(
                {"location": location, "claimed": claimed, "computed": computed}
            )
        if counts[claimed] > 1:
            record_hash_issues.append(
                {"location": location, "claimed": claimed, "issue": "duplicate"}
            )
        for dependency in record.get("dependencies", []):
            if dependency not in known_hashes:
                dependency_issues.append(
                    {"location": location, "dependency": dependency}
                )
        supersedes = record.get("supersedes", [])
        if isinstance(supersedes, str):
            supersedes = [supersedes]
        for previous in supersedes:
            if previous not in known_hashes:
                dependency_issues.append(
                    {"location": location, "supersedes": previous}
                )
        for key, value in _walk(record):
            if key == "object_sha256" and isinstance(value, str):
                object_hashes.add(value)

    object_issues: list[dict[str, str]] = []
    for digest in sorted(object_hashes):
        path = research / "objects" / "sha256" / digest
        if not path.is_file():
            object_issues.append({"sha256": digest, "issue": "missing"})
        elif _sha256_file(path) != digest:
            object_issues.append({"sha256": digest, "issue": "hash mismatch"})

    work_dirs = sorted(
        agent
        for topic in (research / "work").iterdir()
        if topic.is_dir()
        for agent in topic.iterdir()
        if agent.is_dir()
    ) if (research / "work").is_dir() else []
    ledger_pairs = {(r.get("topic"), r.get("agent")) for _, _, r in records}

    return {
        "schema": "primegaps-research-audit-v1",
        "ledger_files": len(ledger_files),
        "records": len(records),
        "outcomes": dict(sorted(Counter(r.get("outcome") for _, _, r in records).items())),
        "verifications": dict(
            sorted(
                Counter(
                    r.get("verification_status", r.get("verification"))
                    for _, _, r in records
                ).items()
            )
        ),
        "parse_issues": parse_issues,
        "record_hash_issues": record_hash_issues,
        "dependency_issues": dependency_issues,
        "object_references": len(object_hashes),
        "object_issues": object_issues,
        "live_work_directories": len(work_dirs),
        "unledgered_work": [
            str(path.relative_to(root))
            for path in work_dirs
            if (path.parent.name, path.name) not in ledger_pairs
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(args.root.resolve())
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()

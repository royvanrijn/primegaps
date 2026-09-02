from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_research_state.py"
SPEC = importlib.util.spec_from_file_location("audit_research_state", SCRIPT)
audit_module = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(audit_module)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_store(root: Path) -> Path:
    objects = root / ".research" / "objects" / "sha256"
    ledger = root / ".research" / "ledger" / "topic"
    work = root / ".research" / "work" / "topic" / "agent"
    objects.mkdir(parents=True)
    ledger.mkdir(parents=True)
    work.mkdir(parents=True)
    payload = b"immutable output\n"
    object_hash = _sha256(payload)
    (objects / object_hash).write_bytes(payload)
    record = {
        "agent": "agent",
        "artifacts": {"outputs": [{"object_sha256": object_hash}]},
        "dependencies": [],
        "kind": "experiment",
        "outcome": "positive",
        "statement": "fixture",
        "topic": "topic",
        "verification_status": "checked",
    }
    record["record_sha256"] = audit_module._record_hash(record)
    (ledger / "agent.jsonl").write_text(json.dumps(record) + "\n")
    return objects / object_hash


def test_audit_accepts_canonical_record_and_object(tmp_path):
    _write_store(tmp_path)
    result = audit_module.audit(tmp_path)
    assert result["records"] == 1
    assert result["record_hash_issues"] == []
    assert result["dependency_issues"] == []
    assert result["object_issues"] == []
    assert result["unledgered_work"] == []


def test_audit_reports_object_and_record_tampering(tmp_path):
    object_path = _write_store(tmp_path)
    object_path.write_bytes(b"changed\n")
    ledger = tmp_path / ".research" / "ledger" / "topic" / "agent.jsonl"
    record = json.loads(ledger.read_text())
    record["statement"] = "changed after hashing"
    ledger.write_text(json.dumps(record) + "\n")
    result = audit_module.audit(tmp_path)
    assert len(result["record_hash_issues"]) == 1
    assert result["object_issues"] == [
        {"sha256": object_path.name, "issue": "hash mismatch"}
    ]

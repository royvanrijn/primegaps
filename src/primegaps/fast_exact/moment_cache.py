"""Append-only exact moment caches for repeated candidate evaluation."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path


HEADER_SCHEMA = "primegaps-fast-exact-moment-cache-v1"
RECORD_SCHEMA = "primegaps-fast-exact-I-moments-v1"
J_RECORD_SCHEMA = "primegaps-fast-exact-J-functionals-v1"
J_STATUS_RECORD_SCHEMA = "primegaps-fast-exact-J-density-statuses-v1"


def canonical_bytes(value) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode()


def value_hash(value) -> str:
    return sha256(canonical_bytes(value)).hexdigest()


def rational_payload(value) -> list[str]:
    return [str(value.numerator), str(value.denominator)]


def j_density_status_record(context_hash, target, statuses):
    return {
        "schema": J_STATUS_RECORD_SCHEMA,
        "context_sha256": context_hash,
        "target": list(target),
        "statuses": [list(status) for status in statuses],
    }


def j_functional_record(context_hash, functional_id, moments):
    return {
        "schema": J_RECORD_SCHEMA,
        "context_sha256": context_hash,
        "functional_id": str(functional_id),
        "moments": {
            f"{exponent[0]},{exponent[1]}": rational_payload(value)
            for exponent, value in moments.items()
        },
    }


class IMomentCache:
    """Candidate-independent exact I(signature, slack) values.

    A signature may appear in several append records as a higher-degree run
    adds new slack powers. Existing values must agree exactly.
    """

    def __init__(self, path: str | Path, *, context: dict, rational):
        self.path = Path(path)
        self.context = context
        self.context_hash = value_hash(context)
        self.rational = rational
        self.values: dict[tuple[int, ...], dict[int, object]] = {}
        if self.path.exists():
            self._load()
        else:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            header = {
                "schema": HEADER_SCHEMA,
                "context": context,
                "context_sha256": self.context_hash,
            }
            self.path.write_text(
                json.dumps(header, sort_keys=True) + "\n"
            )

    def _load(self):
        lines = self.path.read_text().splitlines()
        if not lines:
            raise ValueError("empty moment cache")
        header = json.loads(lines[0])
        if (
            header.get("schema") != HEADER_SCHEMA
            or header.get("context") != self.context
            or header.get("context_sha256") != self.context_hash
        ):
            raise ValueError("moment-cache context mismatch")
        for line_number, line in enumerate(lines[1:], 2):
            record = json.loads(line)
            if (
                record.get("schema") != RECORD_SCHEMA
                or record.get("context_sha256") != self.context_hash
            ):
                raise ValueError(
                    f"invalid moment-cache record on line {line_number}"
                )
            signature = tuple(record["signature"])
            bucket = self.values.setdefault(signature, {})
            for raw_slack, (numerator, denominator) in record["moments"].items():
                slack = int(raw_slack)
                value = self.rational(int(numerator), int(denominator))
                if slack in bucket and bucket[slack] != value:
                    raise ValueError(
                        f"conflicting cached moment {signature}, slack={slack}"
                    )
                bucket[slack] = value

    def missing(self, signature, slacks):
        bucket = self.values.get(tuple(signature), {})
        return tuple(sorted(set(slacks) - set(bucket)))

    def append(self, signature, moments):
        signature = tuple(signature)
        bucket = self.values.setdefault(signature, {})
        fresh = {}
        for slack, value in moments.items():
            slack = int(slack)
            if slack in bucket:
                if bucket[slack] != value:
                    raise ValueError(
                        f"conflicting moment {signature}, slack={slack}"
                    )
                continue
            bucket[slack] = value
            fresh[str(slack)] = rational_payload(value)
        if not fresh:
            return False
        record = {
            "schema": RECORD_SCHEMA,
            "context_sha256": self.context_hash,
            "signature": list(signature),
            "moments": fresh,
        }
        with self.path.open("a") as stream:
            stream.write(json.dumps(record, sort_keys=True) + "\n")
        return True

    def evaluate_atoms(self, atoms):
        total = self.rational(0)
        for (signature, slack), coefficient in atoms.items():
            try:
                moment = self.values[tuple(signature)][int(slack)]
            except KeyError as error:
                raise KeyError(
                    f"missing I moment {signature}, slack={slack}"
                ) from error
            total += coefficient * moment
        return total


class JFunctionalCache:
    """Candidate-independent density/geometry functionals for J.

    Each functional maps a candidate monomial ``x^a z^b`` to its exact
    integral after multiplying by one fixed target density on one fixed support
    cell.  Increasing the polynomial degree merely appends newly requested
    exponents; all old values remain reusable.
    """

    def __init__(self, path: str | Path, *, context: dict, rational):
        self.path = Path(path)
        self.context = context
        self.context_hash = value_hash(context)
        self.rational = rational
        self.values: dict[str, dict[tuple[int, int], object]] = {}
        self.density_statuses: dict[tuple[int, ...], tuple[tuple[int, int], ...]] = {}
        if self.path.exists():
            self._load()
        else:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            header = {
                "schema": HEADER_SCHEMA,
                "context": context,
                "context_sha256": self.context_hash,
            }
            self.path.write_text(json.dumps(header, sort_keys=True) + "\n")

    def _load(self):
        lines = self.path.read_text().splitlines()
        if not lines:
            raise ValueError("empty moment cache")
        header = json.loads(lines[0])
        if (
            header.get("schema") != HEADER_SCHEMA
            or header.get("context") != self.context
            or header.get("context_sha256") != self.context_hash
        ):
            raise ValueError("moment-cache context mismatch")
        for line_number, line in enumerate(lines[1:], 2):
            self._load_record(json.loads(line), line_number)

    def _load_record(self, record, line_number, *, retain=True):
        if record.get("context_sha256") != self.context_hash:
            raise ValueError(
                f"invalid J-functional record on line {line_number}"
            )
        if record.get("schema") == J_STATUS_RECORD_SCHEMA:
            target = tuple(int(value) for value in record["target"])
            statuses = tuple(
                sorted(
                    tuple(int(value) for value in status)
                    for status in record["statuses"]
                )
            )
            if any(
                len(status) != 2 or min(status) < 0 for status in statuses
            ):
                raise ValueError(
                    f"invalid J density status on line {line_number}"
                )
            existing = self.density_statuses.get(target)
            if existing is not None and existing != statuses:
                raise ValueError(
                    f"conflicting cached J density statuses for {target}"
                )
            if retain:
                self.density_statuses[target] = statuses
            return
        if record.get("schema") != J_RECORD_SCHEMA:
            raise ValueError(
                f"invalid J-functional record on line {line_number}"
            )
        functional_id = record["functional_id"]
        bucket = self.values.get(functional_id, {})
        for raw_exponent, payload in record["moments"].items():
            exponent = tuple(int(value) for value in raw_exponent.split(","))
            if len(exponent) != 2:
                raise ValueError(f"invalid exponent on line {line_number}")
            value = self.rational(int(payload[0]), int(payload[1]))
            if exponent in bucket and bucket[exponent] != value:
                raise ValueError(
                    f"conflicting cached J functional {functional_id}"
                )
            if retain:
                if functional_id not in self.values:
                    bucket = self.values.setdefault(functional_id, {})
                bucket[exponent] = value

    def ingest(self, path, *, retain=True):
        """Validate and append a worker-produced record shard."""
        with Path(path).open() as source, self.path.open("a") as destination:
            for line_number, line in enumerate(source, 1):
                record = json.loads(line)
                self._load_record(record, line_number, retain=retain)
                destination.write(json.dumps(record, sort_keys=True) + "\n")

    def append_density_statuses(self, target, statuses):
        """Record the complete nonzero density-status set for one target.

        This small index lets a replay enumerate active support cells without
        reconstructing the target density merely to discover where it is
        nonzero.  It is candidate-independent and must agree exactly if seen
        again.
        """
        target = tuple(int(value) for value in target)
        statuses = tuple(
            sorted(
                set(tuple(int(value) for value in status) for status in statuses)
            )
        )
        if any(len(status) != 2 or min(status) < 0 for status in statuses):
            raise ValueError("J density statuses must be non-negative pairs")
        existing = self.density_statuses.get(target)
        if existing is not None:
            if existing != statuses:
                raise ValueError(
                    f"conflicting J density statuses for {target}"
                )
            return False
        self.density_statuses[target] = statuses
        record = j_density_status_record(self.context_hash, target, statuses)
        with self.path.open("a") as stream:
            stream.write(json.dumps(record, sort_keys=True) + "\n")
        return True

    def missing(self, functional_id, exponents):
        bucket = self.values.get(str(functional_id), {})
        return tuple(sorted(set(exponents) - set(bucket)))

    def append(self, functional_id, moments):
        functional_id = str(functional_id)
        bucket = self.values.setdefault(functional_id, {})
        fresh = {}
        for raw_exponent, value in moments.items():
            exponent = tuple(int(power) for power in raw_exponent)
            if len(exponent) != 2 or min(exponent) < 0:
                raise ValueError("J exponents must be non-negative pairs")
            if exponent in bucket:
                if bucket[exponent] != value:
                    raise ValueError(
                        f"conflicting J functional {functional_id}, {exponent}"
                    )
                continue
            bucket[exponent] = value
            fresh[f"{exponent[0]},{exponent[1]}"] = rational_payload(value)
        if not fresh:
            return False
        record = {
            "schema": J_RECORD_SCHEMA,
            "context_sha256": self.context_hash,
            "functional_id": functional_id,
            "moments": fresh,
        }
        with self.path.open("a") as stream:
            stream.write(json.dumps(record, sort_keys=True) + "\n")
        return True

    def evaluate(self, functional_id, polynomial):
        bucket = self.values.get(str(functional_id), {})
        total = self.rational(0)
        for exponent, coefficient in polynomial.items():
            try:
                moment = bucket[tuple(exponent)]
            except KeyError as error:
                raise KeyError(
                    f"missing J moment {functional_id}, exponent={exponent}"
                ) from error
            total += coefficient * moment
        return total

# Research protocol

Treat this repository as reproducible mathematical research, not as a scratchpad.
Generated research state belongs under `.research/`; it is local, disposable, and
ignored by Git. Tracked source, tests, and documentation must be sufficient to
rebuild any result that the repository relies on.

## Before starting

1. Read `README.md`, `FINDINGS.md`, the relevant source and tests, and the current
   `.research/` ledgers (when present).
2. Search both positive and negative records before repeating work. Read cited
   literature and distinguish published facts from our computations and guesses.
3. Choose one narrow topic and a unique agent name. Work only in
   `.research/work/<topic>/<agent>/` until a result is ready to promote.

## Local research store

Use this layout:

```text
.research/
  ledger/<topic>/<agent>.jsonl   # one append-only shard per agent
  objects/sha256/<digest>        # immutable outputs, sources, and snapshots
  work/<topic>/<agent>/          # active scratch work
  archive/<topic>/<record-hash>/ # negative, superseded, or unfruitful work
  index/                          # disposable indexes rebuilt from ledger shards
```

Never edit another agent's shard or mutable work directory. Objects are immutable;
write through a temporary file and rename only after verifying its SHA-256. Any
combined index is a cache, never a source of truth. Coordinate before modifying
shared tracked code, and keep promotion commits small and topic-specific.

## Ledger contract

Append one canonical JSON object per completed literature check, claim, theorem,
derivation, or experiment. Each record must contain:

- `record_sha256`: SHA-256 of canonical JSON with this field omitted;
- `kind`: `literature`, `claim`, `theorem`, `derivation`, or `experiment`;
- `outcome`: `positive`, `negative`, or `inconclusive`;
- a precise statement and the assumptions/parameter range;
- source citations or URLs, with page/theorem/equation locators where possible;
- hashes of every input, script, snapshot, and output needed to audit the result;
- the command, code revision (plus dirty-diff hash), environment, and random seed;
- dependencies by record hash, the topic, agent, and UTC timestamp;
- verification status: `reported`, `reproduced`, `checked`, or `refuted`.

Do not overwrite history. Corrections append a new record that names the superseded
record. Record failures and counterexamples with the same care as successes. A
numerical observation is not a theorem; a repository claim is not established
until its assumptions and independent verification are recorded.

## Calculate once; replay cheaply

Separate expensive calculation from replay:

1. **Compute** explicitly, never as an implicit dependency of tests or replay.
2. Store outputs by SHA-256 and append the ledger record immediately.
3. Move failed, negative, abandoned, and superseded scripts, notes, and outputs to
   the matching `archive/<topic>/<record-hash>/`; do not silently delete them.
4. **Replay** reads records and recorded outputs only. It verifies schema, hashes,
   provenance, dependency links, and claimed comparisons; it must not rerun the
   underlying numerical or symbolic calculation.

Keep the default pipeline fast: validate ledgers and report recorded conclusions
in seconds. Expensive rebuilds must be separate, opt-in commands. Generated data
may be removed at any time, so every promoted result needs a documented rebuild
command using tracked code and declared inputs. Promote only the minimal durable
code, tests, and conclusions; leave bulk data and abandoned avenues in
`.research/`.

## Handoff

Before leaving a topic, append a record even if the outcome is negative or
inconclusive, archive inactive work, and state the best next step. A new agent
should be able to learn what was tried, what worked, what failed, and why without
rerunning a calculation.

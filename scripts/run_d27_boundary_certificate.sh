#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
STATE="$ROOT/.research/work/failed-experiment-revival/failed-ranker-20260902"
OUTPUT="$STATE/arb-d27-boundary-full-v4.jsonl"
LOG="$STATE/arb-d27-boundary-full-v4.log"
PID_FILE="$STATE/arb-d27-boundary-full-v4.pid"
RESULT="$STATE/arb-d27-boundary-full-v4-result.json"
TMP="$STATE/tmp"

is_running() {
    [[ -s "$PID_FILE" ]] && kill -0 "$(<"$PID_FILE")" 2>/dev/null
}

case "${1:-status}" in
    start)
        if is_running; then
            echo "already running: PID $(<"$PID_FILE")"
            exit 0
        fi
        mkdir -p "$STATE" "$TMP"
        cd "$ROOT"
        nohup setsid env \
            TMPDIR="$TMP" \
            PYTHONPATH=src \
            sage -python scripts/certify_boundary_j_arb.py run \
            --k 48 \
            --candidate "$STATE/candidate-k48-d27-prefix4.json" \
            --verifier "$STATE/multiband_exact.py" \
            --output "$OUTPUT" \
            --workers "${BOUNDARY_J_WORKERS:-4}" \
            --precision "${BOUNDARY_J_PRECISION:-128}" \
            --stride 192 \
            >>"$LOG" 2>&1 </dev/null &
        pid=$!
        temporary="$PID_FILE.tmp"
        printf '%s\n' "$pid" >"$temporary"
        mv "$temporary" "$PID_FILE"
        echo "started PID $pid"
        echo "log: $LOG"
        echo "checkpoint: $OUTPUT"
        ;;
    status)
        if is_running; then
            state="running"
        else
            state="not running"
        fi
        rows=0
        [[ ! -f "$OUTPUT" ]] || rows=$(wc -l <"$OUTPUT")
        echo "$state; checkpointed cells: $rows/1616"
        [[ ! -f "$LOG" ]] || tail -n 12 "$LOG"
        ;;
    finalize)
        if is_running; then
            echo "calculation is still running: PID $(<"$PID_FILE")" >&2
            exit 1
        fi
        cd "$ROOT"
        TMPDIR="$TMP" PYTHONPATH=src sage -python scripts/certify_boundary_j_arb.py finalize \
            --input "$OUTPUT" \
            --output "$RESULT" \
            --unrestricted-j "$STATE/exact-prefix4-D27-unrestricted-J.json" \
            --i-upper "$STATE/exact-prefix4-D27-unrestricted-I.json"
        echo "result: $RESULT"
        ;;
    *)
        echo "usage: $0 {start|status|finalize}" >&2
        exit 2
        ;;
esac

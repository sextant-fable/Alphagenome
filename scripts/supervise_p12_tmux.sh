#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="/home/zelinli6/Alphagenome"
PYTHON="/home/zelinli6/miniconda3/envs/alphagenome/bin/python"
LOG_DIR="$REPO_ROOT/logs/v2_p12_replicate_holdout_retry_20260819"
LOG_PATH="$LOG_DIR/p12_tmux_supervisor.log"
OLD_PID="${1:?existing P12 runner PID is required}"

mkdir -p "$LOG_DIR"
exec >>"$LOG_PATH" 2>&1
echo "[$(date -Is)] supervisor started; watching pid=$OLD_PID"

while kill -0 "$OLD_PID" 2>/dev/null; do
  sleep 30
done

echo "[$(date -Is)] watched runner exited; checking execution and controller state"
sleep 5

cd "$REPO_ROOT"
if "$PYTHON" - <<'PY'
import json
from pathlib import Path
p = Path("alphagenome_custom/metadata/v2/replicate_holdout_v1/p12_replicate_holdout_execution.json")
if p.is_file() and json.loads(p.read_text()).get("status") == "completed":
    raise SystemExit(0)
raise SystemExit(1)
PY
then
  echo "[$(date -Is)] P12 execution is complete; supervisor exiting"
  exit 0
fi

status_json="$($PYTHON -m scripts.v2_phase_controller status)"
status="$(printf '%s' "$status_json" | "$PYTHON" -c 'import json,sys; print(json.load(sys.stdin)["status"])')"
phase="$(printf '%s' "$status_json" | "$PYTHON" -c 'import json,sys; print(json.load(sys.stdin)["current_phase"])')"
if [[ "$phase" != "P12" || "$status" == "COMPLETE" ]]; then
  echo "[$(date -Is)] controller phase=$phase status=$status; supervisor exiting"
  exit 0
fi

echo "[$(date -Is)] reopening P12 for resume"
"$PYTHON" -m scripts.v2_phase_controller reopen --phase P12 \
  --reason "Detached tmux continuation after the prior P12 runner exited; reuse only hash-validated completed records." \
  >>"$LOG_PATH" 2>&1 || true
"$PYTHON" -m scripts.v2_phase_controller approve --gate G4 \
  --scope p12_replicate_holdout_matrix \
  --note "Approved detached tmux continuation of P12-R1; reuse hash-validated completed records, physical GPUs 2 and 3, final-test prohibited." \
  >>"$LOG_PATH" 2>&1
exec "$PYTHON" -m scripts.v2_phase_controller run --phase P12

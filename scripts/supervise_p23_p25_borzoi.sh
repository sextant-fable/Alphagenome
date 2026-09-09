#!/usr/bin/env bash
set -euo pipefail

# Persistent hand-off for the user-authorized P23 -> P24 -> P25 chain. The
# script stops on any failed review; it never bypasses the controller or reads
# chromosome X/locked-test data.
cd /home/zelinli6/Alphagenome
PYTHON=python
LOG=logs/supervise_p23_p25_borzoi.log
mkdir -p logs
exec >>"$LOG" 2>&1
echo "[$(date -Is)] supervisor started"

while true; do
  phase=$($PYTHON - <<'PY'
import json
print(json.load(open('alphagenome_custom/metadata/v2/execution_state.json'))['current_phase'])
PY
)
  status=$($PYTHON - <<'PY'
import json
print(json.load(open('alphagenome_custom/metadata/v2/execution_state.json'))['status'])
PY
)
  echo "[$(date -Is)] phase=$phase status=$status"
  if [[ "$phase" == P23 && "$status" == COMPLETE ]]; then
    $PYTHON -m scripts.v2_phase_controller start-p24-borzoi-adapter
    $PYTHON -m scripts.v2_phase_controller approve --gate G4 --scope p24_borzoi_adapter --note "User-authorized official Borzoi source feasibility audit; native human/mouse heads excluded; physical GPUs 2/3 only; no X or locked-test access."
    CUDA_VISIBLE_DEVICES=3 $PYTHON -m scripts.v2_phase_controller run --phase P24
    continue
  fi
  if [[ "$phase" == P24 && "$status" == COMPLETE ]]; then
    $PYTHON -m scripts.v2_phase_controller start-p25-borzoi-adapter-pilot
    $PYTHON -m scripts.v2_phase_controller approve --gate G4 --scope p25_borzoi_adapter_pilot --note "User-authorized task-faithful Borzoi C. elegans adapter pilot; physical GPU 3; I-V only; no native human/mouse head leaderboard, X or locked-test access."
    CUDA_VISIBLE_DEVICES=3 conda run -n borzoi_py310 python -m scripts.v2_phase_controller run --phase P25
    continue
  fi
  if [[ "$phase" == P25 && "$status" == COMPLETE ]]; then
    echo "[$(date -Is)] P25 completed; supervisor exiting"
    exit 0
  fi
  if [[ "$status" == FAIL || "$status" == APPROVAL_REQUIRED ]]; then
    echo "[$(date -Is)] stopping on controller boundary phase=$phase status=$status"
    exit 1
  fi
  sleep 30
done

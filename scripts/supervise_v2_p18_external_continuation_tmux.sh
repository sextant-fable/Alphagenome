#!/usr/bin/env bash
# Persistent P18 continuation: verified source transfer and I-V-only reprocessing.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
source /home/zelinli6/miniconda3/etc/profile.d/conda.sh
conda activate alphagenome

run_stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log_dir="logs/v2_p18_external_continuation_${run_stamp}"
mkdir -p "$log_dir"

finish() {
  local exit_code="$1"
  printf '%s\tp18_continuation_exit_%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$exit_code" >>"$log_dir/status.tsv"
  python -m scripts.v2_phase_controller status >"$log_dir/final_controller_status.json" 2>&1 || true
}
trap 'code=$?; finish "$code"; exit "$code"' EXIT

printf '%s\tp18_continuation_started\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >"$log_dir/status.tsv"
python -m scripts.v2_phase_controller run --phase P18 2>&1 | tee "$log_dir/controller.log"

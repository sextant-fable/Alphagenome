#!/usr/bin/env bash
# Persistent P18 external-source transfer sidecar. The runner verifies G1 scope.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
source /home/zelinli6/miniconda3/etc/profile.d/conda.sh
conda activate alphagenome

run_stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log_dir="logs/v2_p18_external_download_${run_stamp}"
mkdir -p "$log_dir"

finish() {
  local exit_code="$1"
  printf '%s\tcontroller_exit_%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$exit_code" >>"$log_dir/status.tsv"
  python -m scripts.v2_phase_controller status >"$log_dir/final_controller_status.json" 2>&1 || true
}
trap 'code=$?; finish "$code"; exit "$code"' EXIT

printf '%s\tdownload_sidecar_started\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >"$log_dir/status.tsv"
python -m scripts.run_v2_p18_external_rna_download \
  --transport-manifest results/v2_p18_external_source_transport_20260822T185343Z/p18_external_source_transport.tsv \
  --work-dir shared/source_reads/v2/p18_external_fastqs \
  --audit-dir results/v2_p18_external_download

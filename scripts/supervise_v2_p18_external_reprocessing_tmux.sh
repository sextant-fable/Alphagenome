#!/usr/bin/env bash
# Wait for the registered P18 transfer, then run I--V-only external reprocessing.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
source /home/zelinli6/miniconda3/etc/profile.d/conda.sh
conda activate alphagenome

run_stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log_dir="logs/v2_p18_external_reprocessing_${run_stamp}"
mkdir -p "$log_dir"

finish() {
  local exit_code="$1"
  printf '%s\treprocessing_exit_%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$exit_code" >>"$log_dir/status.tsv"
  python -m scripts.v2_phase_controller status >"$log_dir/final_controller_status.json" 2>&1 || true
}
trap 'code=$?; finish "$code"; exit "$code"' EXIT

printf '%s\twaiting_for_checksum_verified_download\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >"$log_dir/status.tsv"
while [[ ! -f results/v2_p18_external_download/p18_download_summary.json ]]; do
  sleep 60
done
printf '%s\treprocessing_started\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"$log_dir/status.tsv"
python -m scripts.run_v2_p18_external_iv_reprocessing \
  --transport-manifest results/v2_p18_external_source_transport_20260822T185343Z/p18_external_source_transport.tsv \
  --download-audit-dir results/v2_p18_external_download \
  --fastq-root shared/source_reads/v2/p18_external_fastqs \
  --work-dir shared/source_reads/v2/p18_external_iv_work \
  --output-dir alphagenome_custom/tracks/rna_seq_v2_external_p18 \
  --threads 8

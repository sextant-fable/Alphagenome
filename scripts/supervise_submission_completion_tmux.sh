#!/usr/bin/env bash
# Persistent supervisor for the first manuscript-completion lanes.
# It launches only registered P13 work plus metadata-only P14 discovery.
set -euo pipefail

mode="${1:-full}"
if [[ "$mode" != "full" && "$mode" != "resume-p13" ]]; then
  printf 'Usage: %s [full|resume-p13]\n' "$0" >&2
  exit 2
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
source /home/zelinli6/miniconda3/etc/profile.d/conda.sh
conda activate alphagenome

run_stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log_dir="logs/submission_completion_${run_stamp}"
metadata_dir="results/v2_p14_metadata_collection_${run_stamp}"
collision_dir="results/v2_p14_collision_screen_${run_stamp}"
mkdir -p "$log_dir"

write_status() {
  local state="$1"
  printf '%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$state" | tee -a "$log_dir/supervisor_status.tsv"
}

metadata_lane() {
  write_status "lane0_metadata_collection_started"
  python -m scripts.collect_v2_p14_candidate_metadata \
    --candidate 'fig5:developmental_timecourse:PRJNA231838' \
    --candidate 'fig5:early_l4:PRJNA903192' \
    --candidate 'fig5:wt_mutant_series_controls:PRJNA420450' \
    --candidate 'fig6:wild_strain_eqtl:PRJNA669810' \
    --output-dir "$metadata_dir" \
    >"$log_dir/lane0_metadata_collection.log" 2>&1
  python -m scripts.audit_v2_p14_external_collision \
    --external-registry "$metadata_dir/p14_external_candidate_registry.tsv" \
    --output-dir "$collision_dir" \
    >"$log_dir/lane0_collision_audit.log" 2>&1
  write_status "lane0_metadata_collection_completed_manual_gate0_review_required"
}

finish() {
  local exit_code="$1"
  write_status "supervisor_finished_exit_${exit_code}"
  python -m scripts.v2_phase_controller status >"$log_dir/final_controller_status.json" 2>&1 || true
  if [[ "$exit_code" -ne 0 ]]; then
    printf 'Supervisor stopped with exit code %s. Review %s and reattach this tmux session.\n' "$exit_code" "$log_dir" >&2
  fi
}

trap 'code=$?; finish "$code"; exit "$code"' EXIT

write_status "supervisor_started"
metadata_pid=""
if [[ "$mode" == "full" ]]; then
  metadata_lane &
  metadata_pid=$!
else
  write_status "lane0_metadata_collection_reused_from_prior_completed_screen"
fi

write_status "p13_iv_replicate_agreement_started"
python -m scripts.v2_phase_controller run --phase P13 \
  >"$log_dir/p13_controller.log" 2>&1
write_status "p13_iv_replicate_agreement_completed"

if [[ -n "$metadata_pid" ]]; then
  wait "$metadata_pid"
  write_status "lane0_joined"
fi

printf 'P13 and metadata-only P14 preparation completed. No raw external RNA, GPU training, or chromosome-X access was started.\n' | tee -a "$log_dir/supervisor_status.tsv"

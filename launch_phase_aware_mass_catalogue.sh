#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "usage: $0 MANIFEST BASE BASELINE_PROFILE_DIR OUTPUT_DIR WORKERS" >&2
    exit 2
fi

manifest=$1
base=$2
baseline_profile_dir=$3
output_dir=$4
workers=$5
repo=$(cd "$(dirname "$0")" && pwd)
python_bin=${PANSY_PROFILE_PYTHON:-python3}

mkdir -p "$output_dir/process_logs" "$output_dir/worker_logs" "$output_dir/profiles" "$output_dir/scratch"

export repo python_bin manifest base baseline_profile_dir output_dir workers
seq 0 $((workers - 1)) | xargs -P "$workers" -I{} bash -c '
    worker=$1
    "$python_bin" "$repo/run_phase_aware_mass_catalogue.py" \
        --manifest "$manifest" \
        --base "$base" \
        --baseline-profile-dir "$baseline_profile_dir" \
        --output-dir "$output_dir" \
        --worker-index "$worker" \
        --worker-count "$workers" \
        > "$output_dir/process_logs/worker_$(printf "%03d" "$worker").log" 2>&1
' _ {}

date -u +%FT%TZ > "$output_dir/COMPLETE"

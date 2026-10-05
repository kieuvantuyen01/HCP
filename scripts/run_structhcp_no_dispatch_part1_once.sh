#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

export HCP_ABLATION_SHARD_COUNT=2
export HCP_ABLATION_SHARD_INDEX=0
export HCP_ABLATION_CORES="${HCP_ABLATION_CORES:-1-3}"
export HCP_ABLATION_RUN_DIR="${HCP_ABLATION_RUN_DIR:-$repo_root/benchmark-runs/structhcp-no-dispatch-1800s-part1}"

exec "$script_dir/run_ablation_condition_once.sh" no-dispatch-no-two-cut

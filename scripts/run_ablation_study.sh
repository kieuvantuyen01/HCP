#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tier="${1:-screening}"
conditions="${2:-full,no-decomposition,no-dispatch,no-series,no-two-cut,no-repair,undirected-cubic,one-alternating-worker}"
graphs_dir="${HCP_GRAPHS_DIR:-$repo_root/FHCPCS-col}"

case "$tier" in
  smoke)
    cutoff=10
    sample_args=(--sample 5)
    ;;
  screening)
    cutoff=60
    sample_args=()
    ;;
  confirmation)
    cutoff=300
    sample_args=()
    ;;
  final)
    cutoff=1800
    sample_args=()
    ;;
  *)
    echo "Usage: $0 {smoke|screening|confirmation|final} [comma-separated-conditions]" >&2
    exit 2
    ;;
esac

if [[ -n "${HCP_ABLATION_BINARY:-}" ]]; then
  binary="$HCP_ABLATION_BINARY"
  if [[ ! -x "$binary" ]]; then
    echo "HCP_ABLATION_BINARY is not executable: $binary" >&2
    exit 2
  fi
else
  cargo build --release --manifest-path "$repo_root/src/cegar-fix/Cargo.toml"
  binary="$repo_root/src/cegar-fix/target/release/cegar-fix"
fi

python3 "$repo_root/tools/benchmark_runner_ablation.py" \
  --binary "$binary" \
  --graphs-dir "$graphs_dir" \
  --conditions "$conditions" \
  --seeds "1,42,137,777,2026" \
  --timeout "$cutoff" \
  --output "$repo_root/logs/ablation_${tier}_${cutoff}s.jsonl" \
  --tour-dir "$repo_root/scratch/ablation_${tier}_${cutoff}s_tours" \
  "${sample_args[@]}"

python3 "$repo_root/tools/analyze_ablation_results.py" \
  --input "$repo_root/logs/ablation_${tier}_${cutoff}s.jsonl" \
  --timeout "$cutoff" \
  --cactus-out "$repo_root/logs/ablation_${tier}_${cutoff}s_cactus.csv" \
  > "$repo_root/logs/ablation_${tier}_${cutoff}s_report.md"

echo "Report: $repo_root/logs/ablation_${tier}_${cutoff}s_report.md"

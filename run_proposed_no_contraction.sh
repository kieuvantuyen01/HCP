#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINARY="$REPO_DIR/src/cegar-fix/target/release/cegar-fix"
OUTPUT_DIR="$REPO_DIR/benchmark-runs/proposed-no-contraction"

if [[ "${SKIP_BUILD:-0}" != "1" ]]; then
    cargo build --release --manifest-path "$REPO_DIR/src/cegar-fix/Cargo.toml"
fi

mkdir -p "$OUTPUT_DIR"

exec python3 "$REPO_DIR/tools/benchmark_runner_ablation.py" \
    --binary "$BINARY" \
    --graphs-dir "$REPO_DIR/FHCPCS-col" \
    --profiles proposed-no-contraction \
    --profile-config "$REPO_DIR/benchmarks/proposed_variants.json" \
    --seeds 1 \
    --timeout 1800 \
    --expected-status SATISFIABLE \
    --cores 0 \
    --output "$OUTPUT_DIR/results.jsonl" \
    --tour-dir "$OUTPUT_DIR/tours" \
    "$@"

#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
run_dir="${HCP_ABLATION_RUN_DIR:-$repo_root/benchmark-runs/ablation-final-once-1001}"
binary="$repo_root/src/cegar-fix/target/release/cegar-fix"
graphs_dir="${HCP_GRAPHS_DIR:-$repo_root/FHCPCS-col}"
all_conditions="full,no-dispatch,no-dispatch-no-two-cut,no-decomposition,no-two-cut"
seed=1
timeout=1800
memory_limit_mb=7168
shard_count="${HCP_ABLATION_SHARD_COUNT:-1}"
shard_index="${HCP_ABLATION_SHARD_INDEX:-0}"
cores="${HCP_ABLATION_CORES:-1-3}"

if [[ "$#" -ne 1 ]]; then
  echo "Usage: $0 {full|no-dispatch|no-dispatch-no-two-cut|no-decomposition|no-two-cut}" >&2
  exit 2
fi

condition="$1"
case "$condition" in
  full|no-dispatch|no-dispatch-no-two-cut|no-decomposition|no-two-cut) ;;
  *)
    echo "Unsupported final condition: $condition" >&2
    exit 2
    ;;
esac

if ! [[ "$shard_count" =~ ^[1-9][0-9]*$ ]]; then
  echo "HCP_ABLATION_SHARD_COUNT must be a positive integer." >&2
  exit 2
fi
if ! [[ "$shard_index" =~ ^[0-9]+$ ]] || (( shard_index >= shard_count )); then
  echo "HCP_ABLATION_SHARD_INDEX must be in [0, HCP_ABLATION_SHARD_COUNT)." >&2
  exit 2
fi

selected_graph_numbers=()
selected_graphs=()
for ((graph_number = 1; graph_number <= 1001; graph_number++)); do
  if (( (graph_number - 1) % shard_count == shard_index )); then
    selected_graph_numbers+=("$graph_number")
    selected_graphs+=("graph${graph_number}")
  fi
done
selected_graph_numbers_csv="$(IFS=,; echo "${selected_graph_numbers[*]}")"
selected_graphs_csv="$(IFS=,; echo "${selected_graphs[*]}")"

mkdir -p "$run_dir/raw-logs" "$run_dir/runlim" "$run_dir/tours"

runlim_bin="${RUNLIM_BIN:-}"
if [[ -z "$runlim_bin" && -x "$repo_root/bin/runlim" ]]; then
  runlim_bin="$repo_root/bin/runlim"
fi
if [[ -z "$runlim_bin" ]]; then
  runlim_bin="$(command -v runlim || true)"
fi
if [[ -z "$runlim_bin" && -x /usr/local/bin/runlim ]]; then
  runlim_bin=/usr/local/bin/runlim
fi
if [[ -z "$runlim_bin" || ! -x "$runlim_bin" ]]; then
  echo "runlim is required to enforce the ${memory_limit_mb}-MB memory limit." >&2
  echo "Run ./setup_env.sh or set RUNLIM_BIN=/absolute/path/to/runlim." >&2
  exit 2
fi
runlim_version="$("$runlim_bin" --version 2>&1 | head -n 1)"
if [[ "$runlim_version" != "2.0.0rc14" ]]; then
  echo "Expected runlim 2.0.0rc14; found: ${runlim_version:-unknown}" >&2
  echo "Run ./setup_env.sh to install the pinned benchmark version." >&2
  exit 2
fi

lock_dir="$run_dir/.active-condition-run"
if ! mkdir "$lock_dir" 2>/dev/null; then
  echo "Another final ablation script is active in $run_dir." >&2
  echo "Run the four condition scripts one at a time on this machine." >&2
  exit 2
fi
trap 'rmdir "$lock_dir" 2>/dev/null || true' EXIT

graph_count="$(find "$graphs_dir" -maxdepth 1 -type f -name 'graph*.col' | wc -l | tr -d ' ')"
if [[ "$graph_count" != "1001" ]]; then
  echo "Expected exactly 1001 graph*.col files in $graphs_dir; found $graph_count" >&2
  exit 2
fi

if ! git -C "$repo_root" diff --quiet -- || \
   ! git -C "$repo_root" diff --cached --quiet --; then
  echo "Tracked source files have uncommitted changes. Commit or stash them before the final run." >&2
  exit 2
fi

cargo build --release --locked --manifest-path "$repo_root/src/cegar-fix/Cargo.toml"

REPO_ROOT="$repo_root" BINARY="$binary" GRAPHS_DIR="$graphs_dir" RUN_DIR="$run_dir" \
CONDITIONS="$all_conditions" CONDITION="$condition" SEED="$seed" TIMEOUT="$timeout" \
MEMORY_LIMIT_MB="$memory_limit_mb" RUNLIM_BIN="$runlim_bin" \
SHARD_COUNT="$shard_count" SHARD_INDEX="$shard_index" \
SELECTED_GRAPH_NUMBERS="$selected_graph_numbers_csv" CORES="$cores" python3 - <<'PY'
import hashlib
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

repo = Path(os.environ["REPO_ROOT"])
binary = Path(os.environ["BINARY"])
run_dir = Path(os.environ["RUN_DIR"])
graphs_dir = Path(os.environ["GRAPHS_DIR"]).resolve()
metadata_path = run_dir / "metadata.json"
results_path = run_dir / "results.jsonl"

def command_output(command):
    try:
        return subprocess.check_output(command, cwd=repo, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None

def physical_memory_bytes():
    try:
        return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (OSError, ValueError):
        value = command_output(["sysctl", "-n", "hw.memsize"])
        return int(value) if value else None

graph_paths = list(graphs_dir.glob("graph*.col"))
expected_names = {f"graph{i}.col" for i in range(1, 1002)}
actual_names = {path.name for path in graph_paths}
if actual_names != expected_names:
    raise SystemExit(
        "Dataset mismatch: "
        f"missing={sorted(expected_names - actual_names)}, "
        f"unexpected={sorted(actual_names - expected_names)}"
    )
graph_paths.sort(key=lambda path: int(path.stem[5:]))

shard_count = int(os.environ["SHARD_COUNT"])
shard_index = int(os.environ["SHARD_INDEX"])
selected_graph_numbers = [
    int(value) for value in os.environ["SELECTED_GRAPH_NUMBERS"].split(",")
]
expected_graph_numbers = [
    number
    for number in range(1, 1002)
    if (number - 1) % shard_count == shard_index
]
if selected_graph_numbers != expected_graph_numbers:
    raise SystemExit("Internal error: graph shard selection is inconsistent")

manifest_lines = []
for path in graph_paths:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_lines.append(f"{digest}  {path.name}")
manifest_text = "\n".join(manifest_lines) + "\n"

now = datetime.now(timezone.utc).isoformat()
condition = os.environ["CONDITION"]
is_deadline_clean_dispatch = condition == "no-dispatch-no-two-cut"
metadata = {
    "started_at_utc": now,
    "protocol": (
        "deadline-single-seed-clean-dispatch"
        if is_deadline_clean_dispatch else "final-single-run"
    ),
    "protocol_document": (
        "docs/research/experiments/2026-10-05-structhcp-no-dispatch-single-seed-protocol.md"
        if is_deadline_clean_dispatch
        else "docs/research/experiments/2026-09-28-final-single-run-ablation-protocol.md"
    ),
    "instances": 1001,
    "conditions": os.environ["CONDITIONS"].split(","),
    "seed": int(os.environ["SEED"]),
    "timeout_seconds": int(os.environ["TIMEOUT"]),
    "memory_limit_mb": int(os.environ["MEMORY_LIMIT_MB"]),
    "matrix_runs_expected": 1001 * len(os.environ["CONDITIONS"].split(",")),
    "runs_expected": len(selected_graph_numbers),
    "shard_count": shard_count,
    "shard_index": shard_index,
    "shard_number": shard_index + 1,
    "shard_selection": "round-robin by natural graph number",
    "selected_graph_numbers": selected_graph_numbers,
    "cpu_affinity": os.environ["CORES"],
    "execution": "one nonoverlapping graph shard; one solver process at a time; up to three solver-internal threads",
    "git_commit": command_output(["git", "rev-parse", "HEAD"]),
    "git_branch": command_output(["git", "branch", "--show-current"]),
    "git_status_porcelain": command_output(["git", "status", "--porcelain", "--untracked-files=no"]),
    "binary": str(binary.resolve()),
    "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
    "graphs_dir": str(graphs_dir),
    "dataset_manifest": str((run_dir / "instances.sha256").resolve()),
    "dataset_manifest_sha256": hashlib.sha256(manifest_text.encode()).hexdigest(),
    "platform": platform.platform(),
    "machine": platform.machine(),
    "processor": platform.processor() or command_output(["sysctl", "-n", "machdep.cpu.brand_string"]),
    "logical_cpus": os.cpu_count(),
    "memory_bytes": physical_memory_bytes(),
    "runlim_binary": str(Path(os.environ["RUNLIM_BIN"]).resolve()),
    "runlim_binary_sha256": hashlib.sha256(
        Path(os.environ["RUNLIM_BIN"]).read_bytes()
    ).hexdigest(),
    "runlim_version": command_output([os.environ["RUNLIM_BIN"], "--version"]),
    "python": platform.python_version(),
    "rustc": command_output(["rustc", "--version"]),
    "cargo": command_output(["cargo", "--version"]),
}

if results_path.exists() and results_path.stat().st_size:
    if not metadata_path.exists():
        raise SystemExit("Cannot resume: results.jsonl exists without metadata.json")
    previous = json.loads(metadata_path.read_text())
    frozen_keys = (
        "conditions", "seed", "timeout_seconds", "memory_limit_mb",
        "git_commit", "binary_sha256",
        "graphs_dir", "dataset_manifest_sha256", "platform", "machine",
        "processor", "logical_cpus", "memory_bytes", "runlim_binary_sha256",
        "runlim_version", "rustc", "cargo", "shard_count", "shard_index",
        "selected_graph_numbers", "cpu_affinity",
    )
    mismatches = [key for key in frozen_keys if previous.get(key) != metadata.get(key)]
    if mismatches:
        raise SystemExit("Cannot resume: frozen metadata differs for " + ", ".join(mismatches))
    metadata["started_at_utc"] = previous["started_at_utc"]
    metadata["invocations"] = previous.get("invocations", [])
    metadata["completed_conditions"] = previous.get("completed_conditions", [])
    metadata["condition_completed_at_utc"] = previous.get(
        "condition_completed_at_utc", {}
    )
    if "completed_at_utc" in previous:
        metadata["completed_at_utc"] = previous["completed_at_utc"]
    if "shard_completed_at_utc" in previous:
        metadata["shard_completed_at_utc"] = previous["shard_completed_at_utc"]

metadata.setdefault("invocations", []).append({
    "condition": os.environ["CONDITION"],
    "shard_number": shard_index + 1,
    "shard_count": shard_count,
    "started_at_utc": now,
})

(run_dir / "instances.sha256").write_text(manifest_text)
metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
PY

export LC_ALL=C
export OMP_NUM_THREADS=1
export RAYON_NUM_THREADS=3
export PYTHONUNBUFFERED=1

runner=(
  python3 "$repo_root/tools/benchmark_runner_ablation.py"
  --binary "$binary"
  --graphs-dir "$graphs_dir"
  --graphs "$selected_graphs_csv"
  --conditions "$condition"
  --seeds "$seed"
  --timeout "$timeout"
  --cores "$cores"
  --memory-limit-mb "$memory_limit_mb"
  --runlim "$runlim_bin"
  --runlim-log-dir "$run_dir/runlim"
  --output "$run_dir/results.jsonl"
  --tour-dir "$run_dir/tours"
  --raw-log-dir "$run_dir/raw-logs"
)

if command -v caffeinate >/dev/null 2>&1; then
  caffeinate -ims "${runner[@]}" 2>&1 | tee -a "$run_dir/progress-$condition.log"
else
  "${runner[@]}" 2>&1 | tee -a "$run_dir/progress-$condition.log"
fi

RESULTS="$run_dir/results.jsonl" METADATA="$run_dir/metadata.json" \
CONDITION="$condition" SHARD_COUNT="$shard_count" SHARD_INDEX="$shard_index" \
SELECTED_GRAPH_NUMBERS="$selected_graph_numbers_csv" python3 - <<'PY'
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

results_path = Path(os.environ["RESULTS"])
metadata_path = Path(os.environ["METADATA"])
records = [json.loads(line) for line in results_path.read_text().splitlines() if line.strip()]
conditions = (
    "full", "no-dispatch", "no-dispatch-no-two-cut",
    "no-decomposition", "no-two-cut",
)
current_condition = os.environ["CONDITION"]
shard_count = int(os.environ["SHARD_COUNT"])
shard_index = int(os.environ["SHARD_INDEX"])
selected_graph_numbers = {
    int(value) for value in os.environ["SELECTED_GRAPH_NUMBERS"].split(",")
}
expected_features = {
    "full": "name=full seed=1 dispatch=true series=true two_cut=true repair=true directed_cubic=true alternating_workers=3",
    "no-dispatch": "name=no-dispatch seed=1 dispatch=false series=true two_cut=true repair=true directed_cubic=false alternating_workers=3",
    "no-dispatch-no-two-cut": "name=no-dispatch-no-two-cut seed=1 dispatch=false series=true two_cut=false repair=true directed_cubic=false alternating_workers=3",
    "no-decomposition": "name=no-decomposition seed=1 dispatch=false series=false two_cut=false repair=true directed_cubic=false alternating_workers=3",
    "no-two-cut": "name=no-two-cut seed=1 dispatch=true series=true two_cut=false repair=true directed_cubic=true alternating_workers=3",
}
allowed_keys = {(f"graph{i}", condition, 1) for i in range(1, 1002) for condition in conditions}
expected_current_keys = {
    (f"graph{i}", current_condition, 1) for i in selected_graph_numbers
}
actual_keys = {(r["graph"], r["condition"], int(r["seed"])) for r in records}
actual_current_keys = {key for key in actual_keys if key[1] == current_condition}

errors = []
if len(actual_keys) != len(records):
    errors.append(f"duplicate keys: records={len(records)}, unique={len(actual_keys)}")
if not actual_keys <= allowed_keys:
    errors.append(
        f"unexpected graph-condition-seed keys={len(actual_keys - allowed_keys)}"
    )
if actual_current_keys != expected_current_keys:
    errors.append(
        f"{current_condition} key mismatch: "
        f"missing={len(expected_current_keys - actual_current_keys)}, "
        f"unexpected={len(actual_current_keys - expected_current_keys)}"
    )
if any(float(r["timeout_limit"]) != 1800.0 for r in records):
    errors.append("at least one record has a non-1800-second cutoff")
if any(int(r.get("memory_limit_mb", -1)) != 7168 for r in records):
    errors.append("at least one record has a non-7168-MB memory limit")
bad_features = [
    r for r in records
    if r.get("feature_vector") != expected_features.get(r.get("condition"))
]
if bad_features:
    errors.append(f"{len(bad_features)} records have a missing or unexpected feature vector")
if any(r.get("status") == "SATISFIABLE" and r.get("verified") is not True for r in records):
    errors.append("at least one SAT result lacks an independently verified tour")
if any(r.get("status") == "SATISFIABLE" and not Path(r.get("tour_path", "")).is_file() for r in records):
    errors.append("at least one verified SAT record is missing its archived tour")
if any(
    not Path(r.get(field, "")).is_file()
    for r in records
    for field in ("stdout_path", "stderr_path", "runlim_log_path")
):
    errors.append("at least one record is missing a raw stdout/stderr/runlim log")
if any(r.get("peak_memory_mb") is None for r in records):
    errors.append("at least one record is missing peak-memory telemetry")
if any(not r.get("runlim_status") for r in records):
    errors.append("at least one record is missing its final runlim status")
if any(float(r.get("cpu_total_time", -1)) < 0 for r in records):
    errors.append("at least one record is missing valid aggregate CPU time")
unsat_runs = [r for r in records if r.get("status") == "UNSATISFIABLE"]
if unsat_runs:
    errors.append(f"{len(unsat_runs)} runs reported UNSAT on known-Hamiltonian instances")
bad_runs = [r for r in records if r.get("status") in ("ERROR", "UNKNOWN")]
if bad_runs:
    errors.append(f"{len(bad_runs)} runs ended in ERROR/UNKNOWN")

metadata = json.loads(metadata_path.read_text())
if any(r.get("binary_sha256") != metadata["binary_sha256"] for r in records):
    errors.append("at least one record has a different binary hash")
if any(r.get("git_commit") != metadata["git_commit"] for r in records):
    errors.append("at least one record has a different Git revision")
now = datetime.now(timezone.utc).isoformat()
completed_conditions = []
for condition in conditions:
    expected = {(f"graph{i}", condition, 1) for i in range(1, 1002)}
    present = {key for key in actual_keys if key[1] == condition}
    if present == expected:
        completed_conditions.append(condition)
metadata["records"] = len(records)
metadata["status_counts"] = dict(Counter(r["status"] for r in records))
metadata["postflight_errors"] = errors
metadata["completed_conditions"] = completed_conditions
condition_times = metadata.get("condition_completed_at_utc", {})
if not errors and shard_count == 1:
    condition_times[current_condition] = now
if not errors:
    metadata["shard_completed_at_utc"] = now
metadata["condition_completed_at_utc"] = condition_times
if set(completed_conditions) == set(conditions):
    if len(records) != 5005:
        errors.append(f"all conditions are present but expected 5005 records, found {len(records)}")
    elif not errors:
        metadata["completed_at_utc"] = now
metadata["shard_complete"] = not errors
metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")

if errors:
    raise SystemExit("Postflight validation failed: " + "; ".join(errors))
print(
    f"Postflight validation passed for {current_condition}: "
    f"{len(expected_current_keys)} unique, independently checked runs "
    f"(shard {shard_index + 1}/{shard_count})."
)
print("Completed conditions: " + ", ".join(completed_conditions))
PY

python3 "$repo_root/tools/analyze_ablation_results.py" \
  --input "$run_dir/results.jsonl" \
  --timeout "$timeout" \
  --expected-seeds 1 \
  --cactus-out "$run_dir/cactus.csv" \
  | tee "$run_dir/report.md"

echo "Condition completed: $condition"
echo "Current report: $run_dir/report.md"
echo "Raw results: $run_dir/results.jsonl"

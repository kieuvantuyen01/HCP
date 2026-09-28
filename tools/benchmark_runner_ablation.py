#!/usr/bin/env python3
"""Run reproducible, independently verified ablations for the unified HCP solver."""

import argparse
import glob
import hashlib
import json
import os
import re
import resource
import shlex
import shutil
import signal
import subprocess
import sys
import time
from typing import Dict, List, Optional, Set, Tuple


CONDITIONS = (
    "full",
    "no-decomposition",
    "no-dispatch",
    "no-series",
    "no-two-cut",
    "no-repair",
    "undirected-cubic",
    "one-alternating-worker",
)

ALIASES = {
    "proposed": "full",
    "proposed-nodecomposition": "no-decomposition",
    "proposed-norepair": "no-repair",
    "proposed-nocontraction": "no-series",
    "no-contraction": "no-series",
}


def canonical_condition(value: str) -> str:
    normalized = value.strip().lower().replace("_", "-")
    normalized = ALIASES.get(normalized, normalized)
    if normalized not in CONDITIONS:
        raise ValueError(
            f"Unknown condition '{value}'. Expected one of: {', '.join(CONDITIONS)}"
        )
    return normalized


def natural_sort_key(path: str):
    base = os.path.basename(path)
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r"(\d+)", base)]


def parse_col_graph(col_path: str) -> Tuple[int, int, Dict[int, Set[int]]]:
    num_vertices = 0
    num_edges = 0
    adj: Dict[int, Set[int]] = {}
    with open(col_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("c"):
                continue
            if line.startswith("p"):
                parts = line.split()
                if len(parts) == 3:
                    num_vertices, num_edges = int(parts[1]), int(parts[2])
                elif len(parts) >= 4:
                    num_vertices, num_edges = int(parts[2]), int(parts[3])
                for vertex in range(1, num_vertices + 1):
                    adj[vertex] = set()
            elif line.startswith("e "):
                _, u_text, v_text, *_ = line.split()
                u, v = int(u_text), int(v_text)
                adj.setdefault(u, set()).add(v)
                adj.setdefault(v, set()).add(u)
    return num_vertices, num_edges, adj


def parse_hcp_tour(tour_path: str) -> List[int]:
    tour: List[int] = []
    in_tour = False
    with open(tour_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line == "TOUR_SECTION":
                in_tour = True
                continue
            if line in ("-1", "EOF"):
                break
            if in_tour:
                try:
                    tour.append(int(line))
                except ValueError:
                    pass
    return tour


def verify_tour(
    tour: List[int], num_vertices: int, adj: Dict[int, Set[int]]
) -> Tuple[bool, str]:
    if len(tour) != num_vertices:
        return False, f"Tour length {len(tour)} != expected {num_vertices}"
    if len(set(tour)) != len(tour):
        return False, "Tour contains a duplicate vertex"
    for vertex in tour:
        if vertex not in adj:
            return False, f"Vertex {vertex} is absent from the original graph"
    for index, u in enumerate(tour):
        v = tour[(index + 1) % len(tour)]
        if v not in adj[u]:
            return False, f"Edge ({u}, {v}) is absent from the original graph"
    return True, "Valid Hamiltonian cycle"


def select_graphs(
    graphs_arg: Optional[str], sample_arg: Optional[int], graphs_dir: str
) -> List[str]:
    if graphs_arg:
        selected = []
        for raw_item in graphs_arg.split(","):
            item = raw_item.strip()
            candidates = (item, os.path.join(graphs_dir, item),
                          os.path.join(graphs_dir, f"{item}.col"))
            match = next((path for path in candidates if os.path.exists(path)), None)
            if match is None:
                raise FileNotFoundError(f"Graph file not found: {item}")
            selected.append(os.path.abspath(match))
        return selected

    all_graphs = sorted(
        glob.glob(os.path.join(graphs_dir, "*.col")), key=natural_sort_key
    )
    if not all_graphs:
        raise FileNotFoundError(f"No .col graphs found in {graphs_dir}")
    if sample_arg is None or sample_arg >= len(all_graphs):
        return [os.path.abspath(path) for path in all_graphs]
    if sample_arg <= 0:
        raise ValueError("--sample must be positive")
    if sample_arg == 1:
        return [os.path.abspath(all_graphs[0])]
    indices = [
        round(i * (len(all_graphs) - 1) / (sample_arg - 1))
        for i in range(sample_arg)
    ]
    return [os.path.abspath(all_graphs[index]) for index in indices]


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit(repo_root: str) -> Optional[str]:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repo_root, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def load_completed_keys(path: str) -> Set[Tuple[str, str, int, float, str]]:
    keys: Set[Tuple[str, str, int, float, str]] = set()
    if not os.path.exists(path):
        return keys
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            try:
                record = json.loads(line)
                keys.add((
                    record["graph_path"],
                    canonical_condition(str(record["condition"])),
                    int(record["seed"]),
                    float(record["timeout_limit"]),
                    str(record["binary_sha256"]),
                ))
            except (json.JSONDecodeError, KeyError, TypeError, ValueError):
                continue
    return keys


def extract_int_max(pattern: str, text: str) -> int:
    values = [int(match.group(1)) for match in re.finditer(pattern, text)]
    return max(values) if values else 0


def run_benchmark_instance(
    binary_path: str,
    binary_sha256: str,
    commit: Optional[str],
    graph_path: str,
    condition: str,
    seed: int,
    timeout: float,
    tour_dir: str,
    raw_log_dir: Optional[str],
    cores: Optional[str],
    graph_cache: Dict[str, Tuple[int, int, Dict[int, Set[int]]]],
) -> Dict:
    graph_base = os.path.splitext(os.path.basename(graph_path))[0]
    os.makedirs(tour_dir, exist_ok=True)
    tour_path = os.path.abspath(
        os.path.join(tour_dir, f"{graph_base}_{condition}_s{seed}.hcp")
    )
    if os.path.exists(tour_path):
        os.remove(tour_path)

    solver_cmd = [
        os.path.abspath(binary_path),
        "--input", os.path.abspath(graph_path),
        "--ablation", condition,
        "--seed", str(seed),
        "--timeout", str(timeout),
        "--output-tour", tour_path,
    ]
    taskset = shutil.which("taskset")
    command = [taskset, "-c", cores] + solver_cmd if cores and taskset else solver_cmd

    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()
    timed_out = False
    stdout = ""
    stderr = ""
    returncode = 0
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
        # The process-level limit is the published end-to-end wall-clock cutoff.
        # The solver receives the same deadline, but the harness remains the
        # authority when a blocking backend call does not poll its terminator.
        stdout, stderr = process.communicate(timeout=timeout)
        returncode = process.returncode
        wall_time = time.perf_counter() - start
    except subprocess.TimeoutExpired:
        timed_out = True
        wall_time = timeout
        try:
            os.killpg(os.getpgid(process.pid), signal.SIGKILL)
            stdout, stderr = process.communicate(timeout=2.0)
        except (OSError, subprocess.TimeoutExpired):
            pass
        returncode = -1
    except KeyboardInterrupt:
        # The solver runs in its own session so an interrupted benchmark must
        # explicitly reap it before the checkpointed driver exits.
        try:
            os.killpg(os.getpgid(process.pid), signal.SIGKILL)
            process.communicate(timeout=2.0)
        except (OSError, subprocess.TimeoutExpired):
            pass
        raise

    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu_user_time = usage_after.ru_utime - usage_before.ru_utime
    cpu_system_time = usage_after.ru_stime - usage_before.ru_stime

    stdout_path = None
    stderr_path = None
    if raw_log_dir:
        os.makedirs(raw_log_dir, exist_ok=True)
        run_stem = f"{graph_base}_{condition}_s{seed}"
        stdout_path = os.path.abspath(os.path.join(raw_log_dir, f"{run_stem}.stdout.log"))
        stderr_path = os.path.abspath(os.path.join(raw_log_dir, f"{run_stem}.stderr.log"))
        with open(stdout_path, "w", encoding="utf-8") as handle:
            handle.write(stdout)
        with open(stderr_path, "w", encoding="utf-8") as handle:
            handle.write(stderr)

    if "s SATISFIABLE" in stdout:
        status = "SATISFIABLE"
    elif "s UNSATISFIABLE" in stdout:
        status = "UNSATISFIABLE"
    elif timed_out or "TIMEOUT" in stdout or "TIMEOUT" in stderr:
        status = "TIMEOUT"
    elif returncode != 0:
        status = "ERROR"
    else:
        status = "UNKNOWN"

    solve_match = re.search(r"solve_time_sec:\s*([0-9.eE+-]+)", stdout)
    verify_match = re.search(r"verify_time_sec:\s*([0-9.eE+-]+)", stdout)
    solve_time = float(solve_match.group(1)) if solve_match else wall_time
    verify_time = float(verify_match.group(1)) if verify_match else None
    cegar_iterations = extract_int_max(r"\[cegar_engine\]\s+Iter\s+(\d+)", stdout)
    alternating_rounds = extract_int_max(r"\[alternating_pairs\]\s+Round\s+(\d+)", stdout)
    feature_match = re.search(r"^\[ablation\]\s+(.+)$", stdout, re.MULTILINE)
    feature_vector = feature_match.group(1).strip() if feature_match else None

    verified = None
    verification_msg = None
    if status == "SATISFIABLE":
        if not os.path.exists(tour_path):
            verified, verification_msg = False, "Solver did not create a tour file"
            status = "ERROR"
        else:
            if graph_path not in graph_cache:
                graph_cache[graph_path] = parse_col_graph(graph_path)
            num_vertices, _, adjacency = graph_cache[graph_path]
            verified, verification_msg = verify_tour(
                parse_hcp_tour(tour_path), num_vertices, adjacency
            )
            if not verified:
                status = "ERROR"

    return {
        "graph": graph_base,
        "graph_path": os.path.abspath(graph_path),
        "condition": condition,
        "condition_name": condition,
        "seed": seed,
        "status": status,
        "wall_time": round(wall_time, 6),
        "cpu_user_time": round(cpu_user_time, 6),
        "cpu_system_time": round(cpu_system_time, 6),
        "cpu_total_time": round(cpu_user_time + cpu_system_time, 6),
        "solve_time": round(solve_time, 6),
        "verify_time": round(verify_time, 6) if verify_time is not None else None,
        "cegar_iterations": cegar_iterations,
        "alternating_rounds": alternating_rounds,
        "feature_vector": feature_vector,
        "verified": verified,
        "verification_msg": verification_msg,
        "tour_path": tour_path if os.path.exists(tour_path) else None,
        "timeout_limit": timeout,
        "returncode": returncode,
        "error_message": stderr.strip() or None,
        "stdout_path": stdout_path,
        "stderr_path": stderr_path,
        "command": shlex.join(command),
        "binary_sha256": binary_sha256,
        "git_commit": commit,
    }


def main() -> None:
    # Reset both handlers explicitly: shells may launch background jobs with
    # SIGINT ignored, and an interrupted final run must still reap its solver.
    def interrupt_as_keyboard(_signum, _frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGINT, interrupt_as_keyboard)
    signal.signal(signal.SIGTERM, interrupt_as_keyboard)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="src/cegar-fix/target/release/cegar-fix")
    parser.add_argument("--graphs", help="Comma-separated graph names or paths")
    parser.add_argument("--graphs-dir", default="FHCPCS-col")
    parser.add_argument("--sample", type=int)
    parser.add_argument("--conditions", default=",".join(CONDITIONS))
    parser.add_argument("--seeds", default="1,42,137,777,2026")
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--output", default="logs/ablation_results.jsonl")
    parser.add_argument("--cores", default=None,
                        help="Linux taskset core list; omitted on other platforms")
    parser.add_argument("--tour-dir", default="scratch/ablation_tours")
    parser.add_argument("--raw-log-dir",
                        help="Directory for per-run stdout/stderr audit logs")
    parser.add_argument("--balanced-condition-order", action="store_true",
                        help="Rotate the first condition across consecutive graphs")
    parser.add_argument("--no-resume", action="store_true",
                        help="Start a fresh run and replace the JSONL output")
    parser.add_argument("--list-conditions", action="store_true")
    args = parser.parse_args()

    if args.list_conditions:
        print("\n".join(CONDITIONS))
        return
    if args.timeout <= 0:
        parser.error("--timeout must be positive")

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    def resolve(path: str) -> str:
        return path if os.path.isabs(path) else os.path.join(repo_root, path)

    binary_path = resolve(args.binary)
    if not os.path.isfile(binary_path):
        parser.error(f"binary does not exist: {binary_path}")

    try:
        conditions = list(dict.fromkeys(
            canonical_condition(value) for value in args.conditions.split(",")
        ))
        seeds = list(dict.fromkeys(
            int(value.strip()) for value in args.seeds.split(",") if value.strip()
        ))
        if not seeds or any(seed <= 0 for seed in seeds):
            raise ValueError("all seeds must be positive integers")
        selected_graphs = select_graphs(args.graphs, args.sample, resolve(args.graphs_dir))
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))

    output_path = resolve(args.output)
    tour_dir = resolve(args.tour_dir)
    raw_log_dir = resolve(args.raw_log_dir) if args.raw_log_dir else None
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    binary_hash = sha256_file(binary_path)
    completed = set() if args.no_resume else load_completed_keys(output_path)
    existing_hashes = {key[4] for key in completed}
    if existing_hashes and existing_hashes != {binary_hash}:
        parser.error(
            "output contains a different solver binary; choose a new --output "
            "or use --no-resume to replace it"
        )
    commit = git_commit(repo_root)

    scheduled = []
    for graph_index, graph in enumerate(selected_graphs):
        graph_conditions = conditions
        if args.balanced_condition_order and len(conditions) > 1:
            offset = graph_index % len(conditions)
            graph_conditions = conditions[offset:] + conditions[:offset]
        scheduled.extend(
            (graph, condition, seed)
            for condition in graph_conditions
            for seed in seeds
        )
    pending = [
        (graph, condition, seed)
        for graph, condition, seed in scheduled
        if (os.path.abspath(graph), condition, seed, float(args.timeout), binary_hash)
        not in completed
    ]
    print(f"Graphs: {len(selected_graphs)} | Conditions: {conditions}")
    print(f"Seeds: {seeds} | Timeout: {args.timeout}s | Pending runs: {len(pending)}")
    print(f"Output: {output_path}")

    graph_cache: Dict[str, Tuple[int, int, Dict[int, Set[int]]]] = {}
    output_mode = "w" if args.no_resume else "a"
    with open(output_path, output_mode, encoding="utf-8") as output:
        for index, (graph_path, condition, seed) in enumerate(pending, 1):
            record = run_benchmark_instance(
                binary_path, binary_hash, commit, graph_path, condition, seed,
                args.timeout, tour_dir, raw_log_dir, args.cores, graph_cache,
            )
            output.write(json.dumps(record, sort_keys=True) + "\n")
            output.flush()
            print(
                f"[{index}/{len(pending)}] {record['graph']} | {condition} | "
                f"seed={seed} | {record['status']} | {record['wall_time']:.3f}s | "
                f"verified={record['verified']}"
            )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Benchmark interrupted; completed JSONL records are safe to resume.",
              file=sys.stderr)
        raise SystemExit(130)

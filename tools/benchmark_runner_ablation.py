#!/usr/bin/env python3
"""
Automated Multi-Seed Ablation Benchmark Runner & Telemetry Collector for HCP Solvers.

Supports:
- Ablation conditions: C0 (Baseline), C1 (Baseline + H2-Refined), C2 (Baseline + H1), C3 (Full Hybrid H1+H2)
- Named paired profiles and repeated runs
- CPU affinity binding via taskset
- Subprocess timeout containment with process-group kill
- End-to-end 100% sound tour verification against raw .col graph
- Real-time JSONL telemetry logging
"""

import argparse
import glob
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from typing import Dict, List, Optional, Set, Tuple


DEFAULT_PROFILE_CONFIG = "benchmarks/proposed_variants.json"


def natural_sort_key(path: str):
    """Natural sorting key (graph1.col < graph2.col < ... < graph10.col)."""
    base = os.path.basename(path)
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', base)]


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as in_f:
        for chunk in iter(lambda: in_f.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def arguments_to_map(arguments: List[str], slug: str) -> Dict[str, str]:
    """Convert a profile's flag/value list to a map and reject ambiguous profiles."""
    if len(arguments) % 2 != 0:
        raise ValueError(f"profile {slug!r} must contain flag/value pairs")
    result: Dict[str, str] = {}
    for flag, value in zip(arguments[::2], arguments[1::2]):
        if not flag.startswith("--"):
            raise ValueError(f"profile {slug!r} contains a non-flag argument: {flag!r}")
        if flag in result:
            raise ValueError(f"profile {slug!r} repeats flag {flag!r}")
        result[flag] = value
    return result


def validate_single_factor_profiles(profile_config: Dict, profile_arguments: Dict[str, List[str]]) -> None:
    """Enforce the declared one-factor contrasts before any benchmark is run."""
    reference_slug = str(profile_config.get("reference_profile", "")).lower()
    paired_factors = profile_config.get("paired_factors")
    if not reference_slug or not isinstance(paired_factors, dict):
        return
    if reference_slug not in profile_arguments:
        raise ValueError(f"reference profile {reference_slug!r} is missing")

    reference = arguments_to_map(profile_arguments[reference_slug], reference_slug)
    for raw_slug, raw_factor in paired_factors.items():
        slug = str(raw_slug).lower()
        factor = str(raw_factor)
        if slug not in profile_arguments:
            raise ValueError(f"declared paired profile {slug!r} is missing")
        candidate = arguments_to_map(profile_arguments[slug], slug)
        if set(candidate) != set(reference):
            raise ValueError(f"profile {slug!r} does not use the same flags as {reference_slug!r}")
        changed = [flag for flag in reference if reference[flag] != candidate[flag]]
        expected_flag = factor if factor.startswith("--") else "--" + factor.replace("_", "-")
        if changed != [expected_flag]:
            raise ValueError(
                f"profile {slug!r} must differ from {reference_slug!r} only at "
                f"{expected_flag!r}; changed flags: {changed}"
            )


def parse_col_graph(col_path: str) -> Tuple[int, int, Dict[int, Set[int]]]:
    """Parse a DIMACS-like .col graph file into (|V|, |E|, adjacency map)."""
    num_vertices = 0
    num_edges = 0
    adj: Dict[int, Set[int]] = {}
    with open(col_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('c'):
                continue
            if line.startswith('p'):
                parts = line.split()
                if len(parts) == 3:
                    num_vertices = int(parts[1])
                    num_edges = int(parts[2])
                elif len(parts) >= 4:
                    num_vertices = int(parts[2])
                    num_edges = int(parts[3])
                for v in range(1, num_vertices + 1):
                    adj[v] = set()
            elif line.startswith('e '):
                parts = line.split()
                u, v = int(parts[1]), int(parts[2])
                if u not in adj:
                    adj[u] = set()
                if v not in adj:
                    adj[v] = set()
                adj[u].add(v)
                adj[v].add(u)
    return num_vertices, num_edges, adj


def parse_hcp_tour(tour_path: str) -> List[int]:
    """Parse TSPLIB .hcp tour file."""
    tour: List[int] = []
    in_tour = False
    with open(tour_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('c') or line.startswith('NAME') or line.startswith('TYPE') or line.startswith('DIMENSION'):
                continue
            if line == 'TOUR_SECTION':
                in_tour = True
                continue
            if line == '-1' or line == 'EOF':
                break
            if in_tour:
                try:
                    tour.append(int(line))
                except ValueError:
                    continue
    return tour


def verify_tour(tour: List[int], num_vertices: int, adj: Dict[int, Set[int]]) -> Tuple[bool, str]:
    """
    Validates 100% exact-2 degree, uniqueness, length, and edge membership on raw uncontracted graph.
    """
    if len(tour) != num_vertices:
        return False, f"Tour length {len(tour)} != expected {num_vertices}"

    seen = set()
    for v in tour:
        if v in seen:
            return False, f"Duplicate vertex {v} in tour"
        if v < 1 or v > num_vertices or v not in adj:
            return False, f"Vertex {v} outside valid range 1..{num_vertices}"
        seen.add(v)

    for i in range(len(tour)):
        u = tour[i]
        v = tour[(i + 1) % len(tour)]
        if v not in adj[u]:
            return False, f"Edge ({u}, {v}) does not exist in graph"

    return True, "Valid Hamiltonian cycle"


def parse_duration_to_seconds(val_str: str, unit: str) -> float:
    """Convert duration string with unit (s, ms, µs/us, ns) to float seconds."""
    val = float(val_str)
    if unit == 's':
        return val
    elif unit == 'ms':
        return val / 1_000.0
    elif unit in ('µs', 'us'):
        return val / 1_000_000.0
    elif unit == 'ns':
        return val / 1_000_000_000.0
    return val


def select_graphs(graphs_arg: Optional[str], sample_arg: Optional[int], graphs_dir: str) -> List[str]:
    """Select graph paths based on explicit list or equispaced sampling."""
    if graphs_arg:
        selected = []
        for item in graphs_arg.split(','):
            item = item.strip()
            if not item:
                continue
            if os.path.exists(item):
                selected.append(os.path.abspath(item))
            elif os.path.exists(os.path.join(graphs_dir, item)):
                selected.append(os.path.abspath(os.path.join(graphs_dir, item)))
            elif os.path.exists(os.path.join(graphs_dir, f"{item}.col")):
                selected.append(os.path.abspath(os.path.join(graphs_dir, f"{item}.col")))
            else:
                raise FileNotFoundError(f"Graph file not found: {item}")
        return selected

    all_graphs = sorted(glob.glob(os.path.join(graphs_dir, "*.col")), key=natural_sort_key)
    if not all_graphs:
        raise FileNotFoundError(f"No .col graphs found in {graphs_dir}")

    if sample_arg is not None:
        n = int(sample_arg)
        if n <= 0:
            raise ValueError("--sample must be >= 1")
        if n >= len(all_graphs):
            return [os.path.abspath(p) for p in all_graphs]
        if n == 1:
            return [os.path.abspath(all_graphs[0])]
        indices = [round(i * (len(all_graphs) - 1) / (n - 1)) for i in range(n)]
        return [os.path.abspath(all_graphs[i]) for i in indices]

    return [os.path.abspath(p) for p in all_graphs]


def run_benchmark_instance(
    binary_path: str,
    graph_path: str,
    condition: Optional[int],
    profile: Optional[str],
    profile_arguments: Optional[List[str]],
    seed: int,
    timeout: float,
    tour_dir: str,
    raw_dir: str,
    expected_status: str,
    cores: Optional[str] = None,
    graph_cache: Optional[Dict[str, Tuple[int, int, Dict[int, Set[int]]]]] = None,
) -> Dict:
    """Execute a single benchmark run with telemetry collection and soundness verification."""
    graph_base = os.path.splitext(os.path.basename(graph_path))[0]
    os.makedirs(tour_dir, exist_ok=True)
    os.makedirs(raw_dir, exist_ok=True)
    variant_name = profile if profile is not None else f"c{condition}"
    safe_variant = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant_name)
    tour_path = os.path.abspath(os.path.join(tour_dir, f"{graph_base}_{safe_variant}_s{seed}.hcp"))

    if os.path.exists(tour_path):
        try:
            os.remove(tour_path)
        except OSError:
            pass

    solver_cmd = [
        os.path.abspath(binary_path),
        "-i", os.path.abspath(graph_path),
        "--timeout", str(timeout),
        "--output-tour", tour_path,
    ]
    if profile is not None:
        solver_cmd.extend(profile_arguments or [])
    else:
        solver_cmd.extend(["--ablation", str(condition)])

    taskset_bin = shutil.which("taskset")
    if cores and taskset_bin:
        cmd = [taskset_bin, "-c", cores] + solver_cmd
    else:
        cmd = solver_cmd

    start_time = time.perf_counter()
    timed_out = False
    stdout = ""
    stderr = ""
    returncode = 0

    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            preexec_fn=os.setsid,
        )
        # Give a small buffer beyond solver's internal timeout before hard SIGKILL
        py_timeout = timeout + max(5.0, timeout * 0.1)
        stdout, stderr = proc.communicate(timeout=py_timeout)
        wall_time = time.perf_counter() - start_time
        returncode = proc.returncode
    except subprocess.TimeoutExpired:
        timed_out = True
        wall_time = timeout
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except OSError:
            pass
        try:
            stdout, stderr = proc.communicate(timeout=2.0)
        except Exception:
            stdout, stderr = "", ""
        returncode = -1

    # Classify run status
    if "s SATISFIABLE" in stdout:
        status = "SATISFIABLE"
    elif timed_out or "TIMEOUT" in stdout:
        status = "TIMEOUT"
    elif "s UNKNOWN" in stdout:
        status = "UNKNOWN"
    elif "s UNSATISFIABLE" in stdout:
        status = "UNSATISFIABLE"
    elif returncode != 0:
        status = "ERROR"
    else:
        status = "UNKNOWN"

    if expected_status == "SATISFIABLE" and status == "UNSATISFIABLE":
        status = "INVALID_UNSATISFIABLE"

    raw_base = f"{graph_base}_{safe_variant}_s{seed}"
    stdout_path = os.path.abspath(os.path.join(raw_dir, raw_base + ".stdout.txt"))
    stderr_path = os.path.abspath(os.path.join(raw_dir, raw_base + ".stderr.txt"))
    with open(stdout_path, 'w', encoding='utf-8') as stdout_f:
        stdout_f.write(stdout)
    with open(stderr_path, 'w', encoding='utf-8') as stderr_f:
        stderr_f.write(stderr)

    # Extract solver reported overall time
    total_time = wall_time
    m_time = re.search(r'overall time\s*=\s*([0-9\.]+)\s*(s|ms|µs|us|ns)', stdout)
    if m_time:
        total_time = parse_duration_to_seconds(m_time.group(1), m_time.group(2))

    # Extract SAT solving time
    sat_times = []
    for m in re.finditer(r'sat solving time\s*=\s*([0-9\.]+)\s*(s|ms|µs|us|ns)', stdout):
        sat_times.append(parse_duration_to_seconds(m.group(1), m.group(2)))
    sat_time = sum(sat_times) if sat_times else total_time

    # Extract CEGAR increments count
    m_inc = re.search(r'overall incremented number\s*=\s*(\d+)', stdout)
    if m_inc:
        increments = int(m_inc.group(1))
    else:
        incs = [int(m.group(1)) for m in re.finditer(r'incremented number\s*=\s*(\d+)', stdout)]
        increments = max(incs) if incs else 0

    # Extract block clauses added
    m_block = re.search(r'overall number of added block clauses\s*=\s*(\d+)', stdout)
    if m_block:
        block_clause_count = int(m_block.group(1))
    else:
        blocks = [int(m.group(1)) for m in re.finditer(r'number of added block clauses\s*=\s*(\d+)', stdout)]
        block_clause_count = sum(blocks) if blocks else 0

    # Extract subcycle count from latest round
    subcycles = [int(m.group(1)) for m in re.finditer(r'number of subcycles found\s*=\s*(\d+)', stdout)]
    num_cycles = subcycles[-1] if subcycles else 0

    # Verify tour if reported SATISFIABLE
    verified = None
    verification_msg = ""
    if status == "SATISFIABLE":
        if os.path.exists(tour_path):
            if graph_cache is not None and graph_path in graph_cache:
                num_v, num_e, adj = graph_cache[graph_path]
            else:
                num_v, num_e, adj = parse_col_graph(graph_path)
                if graph_cache is not None:
                    graph_cache[graph_path] = (num_v, num_e, adj)

            tour = parse_hcp_tour(tour_path)
            valid, msg = verify_tour(tour, num_v, adj)
            verified = valid
            verification_msg = msg
            if not valid:
                status = "INVALID_SATISFIABLE"
        else:
            verified = False
            verification_msg = "Tour file was not created by solver"
            status = "INVALID_SATISFIABLE"

    condition_labels = {
        0: "C0 (Baseline)",
        1: "C1 (Baseline + H2-Refined)",
        2: "C2 (Baseline + H1)",
        3: "C3 (Full Hybrid H1+H2)",
    }

    record = {
        "graph": graph_base,
        "graph_path": graph_path,
        "condition": condition,
        "condition_name": (
            profile if profile is not None
            else condition_labels.get(condition, f"C{condition}")
        ),
        "variant": variant_name,
        "solver_arguments": solver_cmd[1:],
        "seed": seed,
        "solver_seed_applied": False,
        "expected_status": expected_status,
        "status": status,
        "wall_time": round(wall_time, 6),
        "total_time": round(total_time, 6),
        "sat_time": round(sat_time, 6),
        "increments": increments,
        "increment_count": increments,
        "block_clause_count": block_clause_count,
        "num_cycles": num_cycles,
        "verified": verified,
        "verification_msg": verification_msg if verified is not None else None,
        "tour_path": tour_path if os.path.exists(tour_path) else None,
        "stdout_path": stdout_path,
        "stderr_path": stderr_path,
        "timeout_limit": timeout,
        "error_message": stderr.strip() if returncode != 0 and stderr else None,
    }

    return record


def main():
    parser = argparse.ArgumentParser(
        description="Automated multi-seed ablation benchmark runner and telemetry collector."
    )
    parser.add_argument(
        "--binary",
        default="src/cegar-fix/target/release/cegar-fix",
        help="Path to cegar-fix binary (default: src/cegar-fix/target/release/cegar-fix)",
    )
    parser.add_argument(
        "--graphs",
        default=None,
        help="Comma-separated graph names or paths (e.g. graph1,graph10)",
    )
    parser.add_argument(
        "--graphs-dir",
        default="FHCPCS-col",
        help="Directory containing benchmark .col graphs (default: FHCPCS-col)",
    )
    parser.add_argument(
        "--sample",
        type=int,
        default=None,
        help="Sample N equispaced graphs from graphs-dir",
    )
    parser.add_argument(
        "--conditions",
        default=None,
        help="Comma-separated ablation indices. Defaults to 0,1,2,3 when --profiles is absent.",
    )
    parser.add_argument(
        "--profiles",
        default=None,
        help=(
            "Comma-separated named profiles: proposed,proposed-no-decomposition,"
            "proposed-no-contraction,proposed-no-repair. "
            "Cannot be combined with --conditions."
        ),
    )
    parser.add_argument(
        "--profile-config",
        default=DEFAULT_PROFILE_CONFIG,
        help=f"Named-profile JSON file (default: {DEFAULT_PROFILE_CONFIG})",
    )
    parser.add_argument(
        "--seeds",
        default="1,42,137,777,2026",
        help="Comma-separated repetition labels; the current binary does not expose a SAT seed.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=1800.0,
        help="Timeout in seconds per instance (default: 1800.0)",
    )
    parser.add_argument(
        "--expected-status",
        choices=("SATISFIABLE", "UNKNOWN"),
        default="SATISFIABLE",
        help="Expected benchmark status; FHCP Challenge instances are known Hamiltonian.",
    )
    parser.add_argument(
        "--output",
        default="logs/ablation_results.jsonl",
        help="Output JSONL telemetry file (default: logs/ablation_results.jsonl)",
    )
    parser.add_argument(
        "--cores",
        default="0,1",
        help="CPU core affinity string for taskset (default: 0,1)",
    )
    parser.add_argument(
        "--tour-dir",
        default="scratch/ablation_tours",
        help="Directory to store intermediate generated tours (default: scratch/ablation_tours)",
    )
    parser.add_argument(
        "--raw-dir",
        default=None,
        help="Directory for per-run stdout/stderr (default: raw beside the JSONL output).",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Append to an existing JSONL file and skip completed graph/variant/seed tuples.",
    )

    args = parser.parse_args()

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    binary_path = os.path.join(repo_root, args.binary) if not os.path.isabs(args.binary) else args.binary
    if not os.path.exists(binary_path):
        print(f"Error: binary {binary_path} does not exist. Please compile it first.", file=sys.stderr)
        sys.exit(1)

    graphs_dir = os.path.join(repo_root, args.graphs_dir) if not os.path.isabs(args.graphs_dir) else args.graphs_dir
    selected_graphs = select_graphs(args.graphs, args.sample, graphs_dir)

    if args.conditions and args.profiles:
        parser.error("--conditions and --profiles are mutually exclusive")

    profile_config_path: Optional[str] = None
    profile_arguments: Dict[str, List[str]] = {}
    profile_labels: Dict[str, str] = {}
    profiles: List[Optional[str]]
    conditions: List[Optional[int]]
    if args.profiles:
        profile_config_path = (
            args.profile_config if os.path.isabs(args.profile_config)
            else os.path.join(repo_root, args.profile_config)
        )
        with open(profile_config_path, 'r', encoding='utf-8') as profile_f:
            profile_config = json.load(profile_f)
        for variant in profile_config.get("variants", []):
            slug = str(variant["slug"]).lower()
            profile_arguments[slug] = [str(arg) for arg in variant.get("arguments", [])]
            profile_labels[slug] = str(variant.get("name", slug))
        try:
            validate_single_factor_profiles(profile_config, profile_arguments)
        except ValueError as exc:
            parser.error(str(exc))
        requested_profiles = [p.strip().lower() for p in args.profiles.split(',') if p.strip()]
        unknown_profiles = [p for p in requested_profiles if p not in profile_arguments]
        if unknown_profiles:
            parser.error(f"unknown profile(s): {', '.join(unknown_profiles)}")
        profiles = requested_profiles
        conditions = [None]
    else:
        conditions = [
            int(c.strip())
            for c in (args.conditions or "0,1,2,3").split(',')
            if c.strip()
        ]
        profiles = [None]
    seeds = [int(s.strip()) for s in args.seeds.split(',') if s.strip()]

    output_path = os.path.join(repo_root, args.output) if not os.path.isabs(args.output) else args.output
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    if os.path.exists(output_path) and not args.resume:
        parser.error(f"output already exists: {output_path}; choose a new path or pass --resume")

    completed: Set[Tuple[str, str, int]] = set()
    if args.resume and os.path.exists(output_path):
        with open(output_path, 'r', encoding='utf-8') as existing_f:
            for line in existing_f:
                try:
                    old = json.loads(line)
                    completed.add((old["graph"], old.get("variant", f"c{old.get('condition')}"), int(old["seed"])))
                except (json.JSONDecodeError, KeyError, TypeError, ValueError):
                    continue

    tour_dir = os.path.join(repo_root, args.tour_dir) if not os.path.isabs(args.tour_dir) else args.tour_dir
    os.makedirs(tour_dir, exist_ok=True)
    raw_dir = (
        os.path.abspath(args.raw_dir)
        if args.raw_dir
        else os.path.join(os.path.dirname(output_path), "raw")
    )
    os.makedirs(raw_dir, exist_ok=True)

    try:
        git_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repo_root, text=True
        ).strip()
        git_dirty = bool(subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=repo_root,
            text=True,
        ).strip())
    except (OSError, subprocess.CalledProcessError):
        git_commit = None
        git_dirty = None

    variants = profiles if args.profiles else [f"c{c}" for c in conditions]
    manifest_path = output_path + ".manifest.json"
    manifest = {
        "schema_version": 1,
        "command": sys.argv,
        "git_commit": git_commit,
        "git_dirty_tracked_files": git_dirty,
        "binary": os.path.abspath(binary_path),
        "binary_sha256": sha256_file(binary_path),
        "profile_config": profile_config_path,
        "profile_config_sha256": (
            sha256_file(profile_config_path) if profile_config_path else None
        ),
        "variants": variants,
        "timeout_sec": args.timeout,
        "expected_status": args.expected_status,
        "repetition_labels": seeds,
        "solver_seed_applied": False,
        "graphs": selected_graphs,
        "cpu_affinity": args.cores,
        "raw_dir": raw_dir,
        "tour_dir": os.path.abspath(tour_dir),
    }
    if args.resume and os.path.exists(manifest_path):
        with open(manifest_path, 'r', encoding='utf-8') as manifest_f:
            old_manifest = json.load(manifest_f)
        comparable_keys = (
            "git_commit", "binary_sha256", "profile_config_sha256", "variants",
            "timeout_sec", "expected_status", "repetition_labels", "graphs",
            "cpu_affinity"
        )
        mismatches = [
            key for key in comparable_keys if old_manifest.get(key) != manifest.get(key)
        ]
        if mismatches:
            parser.error(
                "cannot resume with changed provenance fields: " + ", ".join(mismatches)
            )
    elif not args.resume or not os.path.exists(manifest_path):
        with open(manifest_path, 'w', encoding='utf-8') as manifest_f:
            json.dump(manifest, manifest_f, indent=2)
            manifest_f.write("\n")

    total_runs = len(selected_graphs) * len(variants) * len(seeds)
    print("=" * 70)
    print("Ablation Benchmark Runner Initialized")
    print(f"Binary:     {binary_path}")
    print(f"Graphs:     {len(selected_graphs)} instances selected")
    print(f"Variants:   {variants}")
    print(f"Seeds:      {seeds}")
    print(f"Timeout:    {args.timeout}s")
    print(f"Cores:      {args.cores}")
    print(f"Output:     {output_path}")
    print(f"Total Runs: {total_runs}")
    print("=" * 70)

    graph_cache: Dict[str, Tuple[int, int, Dict[int, Set[int]]]] = {}
    run_idx = len(completed)

    with open(output_path, 'a', encoding='utf-8') as out_f:
        for graph_path in selected_graphs:
            graph_base = os.path.splitext(os.path.basename(graph_path))[0]
            jobs = (
                [(None, profile) for profile in profiles]
                if args.profiles
                else [(condition, None) for condition in conditions]
            )
            for condition, profile in jobs:
                for seed in seeds:
                    variant_name = profile if profile is not None else f"c{condition}"
                    if (graph_base, variant_name, seed) in completed:
                        continue
                    run_idx += 1
                    record = run_benchmark_instance(
                        binary_path=binary_path,
                        graph_path=graph_path,
                        condition=condition,
                        profile=profile,
                        profile_arguments=(profile_arguments.get(profile) if profile else None),
                        seed=seed,
                        timeout=args.timeout,
                        tour_dir=tour_dir,
                        raw_dir=raw_dir,
                        expected_status=args.expected_status,
                        cores=args.cores,
                        graph_cache=graph_cache,
                    )
                    if profile is not None:
                        record["condition_name"] = profile_labels[profile]

                    out_f.write(json.dumps(record) + "\n")
                    out_f.flush()

                    status_str = record["status"]
                    v_str = f"Verified: {record['verified']}" if record["verified"] is not None else ""
                    print(
                        f"[{run_idx:4d}/{total_runs:4d}] {graph_base:12s} | "
                        f"Variant: {profile_labels.get(variant_name, variant_name):28s} | Seed: {seed:4d} | "
                        f"Status: {status_str:12s} | "
                        f"Time: {record['wall_time']:8.3f}s | "
                        f"Inc: {record['increments']:3d} | "
                        f"{v_str}"
                    )

    print("\nBenchmark run completed successfully.")
    print(f"Results saved to: {output_path}")


if __name__ == '__main__':
    main()

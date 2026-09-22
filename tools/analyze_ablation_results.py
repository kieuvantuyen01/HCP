#!/usr/bin/env python3
"""
SAT Benchmark Telemetry Analyzer & Statistical Significance Validator.

Performs:
1. SAT Competition PAR-2 score calculation (2x timeout penalty)
2. Mean & Median runtime analysis on solved instances
3. Paired Wilcoxon Signed-Rank Test (C3 vs C0) with SciPy and robust fallback
4. Cumulative Cactus plot CSV export (c0_time, c1_time, c2_time, c3_time)
5. Comprehensive Markdown report formatting
"""

import argparse
import csv
import json
import math
import os
import statistics
import sys
from typing import Dict, List, Optional, Tuple


CONDITION_NAMES = {
    0: "Baseline (C0)",
    1: "Baseline + H2-Refined (C1)",
    2: "Baseline + H1 (C2)",
    3: "Full Hybrid H1+H2 (C3)",
}


def record_group(record: Dict):
    """Use a named profile when present, otherwise the legacy condition id."""
    return record.get("variant", record.get("condition"))


def load_records(jsonl_path: str) -> List[Dict]:
    """Load JSONL benchmark records from file."""
    records = []
    if not os.path.exists(jsonl_path):
        raise FileNotFoundError(f"Results file not found: {jsonl_path}")
    with open(jsonl_path, 'r', encoding='utf-8') as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"Warning: Failed to parse JSON on line {line_no}: {e}", file=sys.stderr)
    return records


def compute_par2(records: List[Dict], default_timeout: float) -> Tuple[float, List[float]]:
    """
    Compute PAR-2 score:
    PAR-2 = 1/N * sum(t_i if solved else 2 * T_timeout)
    """
    if not records:
        return 0.0, []

    costs = []
    for r in records:
        timeout_limit = float(r.get("timeout_limit", default_timeout))
        # Instance is strictly considered solved only if SATISFIABLE and certified verified
        is_solved = (r.get("status") == "SATISFIABLE" and r.get("verified") is True)
        if is_solved:
            # Use wall_time or total_time
            t = float(r.get("wall_time", r.get("total_time", timeout_limit)))
            costs.append(min(t, timeout_limit))
        else:
            costs.append(2.0 * timeout_limit)

    return sum(costs) / len(costs), costs


def manual_wilcoxon_signed_rank(x: List[float], y: List[float]) -> Tuple[float, float, str]:
    """
    Fallback paired Wilcoxon signed-rank test calculation.
    Compares x (e.g. C3) vs y (e.g. C0).
    """
    diffs = [b - a for a, b in zip(x, y)]
    non_zero = [(abs(d), 1 if d > 0 else -1) for d in diffs if abs(d) > 1e-9]
    if not non_zero:
        return 0.0, 1.0, "Identical performance across all paired instances (no differences)"

    sorted_diffs = sorted(non_zero, key=lambda item: item[0])
    ranks = []
    i = 0
    n = len(sorted_diffs)
    while i < n:
        j = i
        while j < n and abs(sorted_diffs[j][0] - sorted_diffs[i][0]) < 1e-9:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks.append((avg_rank, sorted_diffs[k][1]))
        i = j

    w_pos = sum(r for r, s in ranks if s > 0)
    w_neg = sum(r for r, s in ranks if s < 0)
    w_stat = min(w_pos, w_neg)

    n_pts = len(ranks)
    mean_w = n_pts * (n_pts + 1) / 4.0
    var_w = n_pts * (n_pts + 1) * (2 * n_pts + 1) / 24.0

    if var_w > 0:
        z = (w_stat - mean_w) / math.sqrt(var_w)
        p_val = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / math.sqrt(2.0))))
        p_val = max(0.0, min(1.0, p_val))
    else:
        p_val = 1.0

    conclusion = "Statistically Significant (p < 0.05)" if p_val < 0.05 else "Not Statistically Significant (p >= 0.05)"
    return w_stat, p_val, conclusion


def perform_wilcoxon_test(
    records_a: List[Dict],
    records_b: List[Dict],
    name_a: str,
    name_b: str,
    default_timeout: float,
) -> Dict:
    """Perform paired Wilcoxon signed-rank test on matched instances (graph, seed)."""
    # Key by (graph, seed)
    map_a = {(r["graph"], r["seed"]): r for r in records_a}
    map_b = {(r["graph"], r["seed"]): r for r in records_b}

    common_keys = sorted(list(set(map_a.keys()) & set(map_b.keys())))
    if not common_keys:
        return {
            "name_a": name_a,
            "name_b": name_b,
            "pairs_count": 0,
            "w_stat": None,
            "p_value": None,
            "conclusion": "No common instances found for paired test",
            "wins_a": 0,
            "wins_b": 0,
            "ties": 0,
        }

    costs_a = []
    costs_b = []
    wins_a = 0
    wins_b = 0
    ties = 0

    for k in common_keys:
        ra = map_a[k]
        rb = map_b[k]

        timeout_a = float(ra.get("timeout_limit", default_timeout))
        timeout_b = float(rb.get("timeout_limit", default_timeout))

        cost_a = float(ra.get("wall_time", ra.get("total_time", timeout_a))) if (ra.get("status") == "SATISFIABLE" and ra.get("verified") is True) else 2.0 * timeout_a
        cost_b = float(rb.get("wall_time", rb.get("total_time", timeout_b))) if (rb.get("status") == "SATISFIABLE" and rb.get("verified") is True) else 2.0 * timeout_b

        costs_a.append(cost_a)
        costs_b.append(cost_b)

        if cost_a < cost_b - 1e-6:
            wins_a += 1
        elif cost_b < cost_a - 1e-6:
            wins_b += 1
        else:
            ties += 1

    diffs = [ca - cb for ca, cb in zip(costs_a, costs_b)]
    non_zero_diffs = [d for d in diffs if abs(d) > 1e-9]

    w_stat: Optional[float] = None
    p_value: Optional[float] = None
    conclusion = ""

    # Attempt scipy first
    use_scipy = False
    try:
        from scipy.stats import wilcoxon
        if len(non_zero_diffs) > 0:
            res = wilcoxon(costs_a, costs_b, zero_method='pratt')
            w_stat = float(res.statistic)
            p_value = float(res.pvalue)
            use_scipy = True
            conclusion = "Statistically Significant (p < 0.05)" if p_value < 0.05 else "Not Statistically Significant (p >= 0.05)"
        else:
            w_stat = 0.0
            p_value = 1.0
            conclusion = "Identical performance across all paired instances (no differences)"
    except Exception:
        use_scipy = False

    if not use_scipy and w_stat is None:
        w_stat, p_value, conclusion = manual_wilcoxon_signed_rank(costs_a, costs_b)

    return {
        "name_a": name_a,
        "name_b": name_b,
        "pairs_count": len(common_keys),
        "w_stat": w_stat,
        "p_value": p_value,
        "conclusion": conclusion,
        "wins_a": wins_a,
        "wins_b": wins_b,
        "ties": ties,
    }


def export_cactus_csv(records: List[Dict], csv_path: str, conditions: List):
    """Export sorted solved runtimes for each condition into cactus plot CSV format."""
    solved_times: Dict[int, List[float]] = {}
    for c in conditions:
        c_records = [
            r for r in records
            if record_group(r) == c and r.get("status") == "SATISFIABLE" and r.get("verified") is True
        ]
        times = sorted([float(r.get("wall_time", r.get("total_time", 0.0))) for r in c_records])
        solved_times[c] = times

    max_solved = max([len(t) for t in solved_times.values()]) if solved_times else 0

    os.makedirs(os.path.dirname(os.path.abspath(csv_path)), exist_ok=True)
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        header = ["solved_instance_idx"] + [
            f"c{c}_time" if isinstance(c, int) else f"{c}_time"
            for c in conditions
        ]
        writer.writerow(header)

        for i in range(max_solved):
            idx = i + 1
            row = [idx]
            for c in conditions:
                times = solved_times.get(c, [])
                if i < len(times):
                    row.append(f"{times[i]:.6f}")
                else:
                    row.append("")
            writer.writerow(row)


def analyze_records(records: List[Dict], default_timeout: float) -> Tuple[List[Dict], Dict]:
    """Calculate summary statistics per condition."""
    conditions = sorted(set(record_group(r) for r in records), key=str)
    summary = []

    for c in conditions:
        cond_records = [r for r in records if record_group(r) == c]
        total_runs = len(cond_records)

        solved_records = [
            r for r in cond_records
            if r.get("status") == "SATISFIABLE" and r.get("verified") is True
        ]
        solved_count = len(solved_records)
        timeout_count = sum(1 for r in cond_records if r.get("status") == "TIMEOUT")
        unsat_count = sum(1 for r in cond_records if r.get("status") == "UNSATISFIABLE")
        error_count = sum(
            1 for r in cond_records
            if str(r.get("status", "")).startswith("INVALID")
            or r.get("status") == "ERROR"
            or (r.get("status") == "SATISFIABLE" and r.get("verified") is not True)
        )

        solve_rate = (solved_count / total_runs * 100.0) if total_runs > 0 else 0.0
        par2, _ = compute_par2(cond_records, default_timeout)

        solved_times = [float(r.get("wall_time", r.get("total_time", 0.0))) for r in solved_records]
        mean_time = statistics.mean(solved_times) if solved_times else float('nan')
        median_time = statistics.median(solved_times) if solved_times else float('nan')
        std_time = statistics.stdev(solved_times) if len(solved_times) > 1 else 0.0

        summary.append({
            "condition": c,
            "name": cond_records[0].get("condition_name", CONDITION_NAMES.get(c, str(c))),
            "runs": total_runs,
            "solved": solved_count,
            "timeouts": timeout_count,
            "unsat": unsat_count,
            "errors": error_count,
            "solve_rate": solve_rate,
            "par2": par2,
            "mean_solved_time": mean_time,
            "median_solved_time": median_time,
            "std_solved_time": std_time,
        })

    return summary, {c: [r for r in records if record_group(r) == c] for c in conditions}


def print_markdown_report(
    summary: List[Dict],
    wilcoxon_res: Optional[Dict],
    cactus_out: str,
):
    """Print standard SAT Competition formatted markdown report."""
    print("\n# SAT Benchmark & Multi-Seed Ablation Report")
    print("## 1. Summary Performance Matrix (PAR-2 & Solved Rates)\n")

    print("| Condition/Variant | Description | Runs | Solved | Timeouts | Unsat/Err | Solve Rate (%) | Mean Solved (s) | Median Solved (s) | PAR-2 Score (s) |")
    print("|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|")
    for s in summary:
        mean_str = f"{s['mean_solved_time']:.3f}" if not math.isnan(s['mean_solved_time']) else "N/A"
        median_str = f"{s['median_solved_time']:.3f}" if not math.isnan(s['median_solved_time']) else "N/A"
        group_label = f"C{s['condition']}" if isinstance(s['condition'], int) else str(s['condition'])
        print(
            f"| {group_label} | {s['name']} | "
            f"{s['runs']} | {s['solved']} | {s['timeouts']} | {s['unsat'] + s['errors']} | "
            f"{s['solve_rate']:.1f}% | {mean_str} | {median_str} | "
            f"**{s['par2']:.3f}** |"
        )

    print("\n## 2. Paired Statistical Significance Testing (Wilcoxon Signed-Rank Test)\n")
    if wilcoxon_res and wilcoxon_res.get("pairs_count", 0) > 0:
        w_str = f"{wilcoxon_res['w_stat']:.2f}" if wilcoxon_res['w_stat'] is not None else "N/A"
        p_str = f"{wilcoxon_res['p_value']:.4e}" if wilcoxon_res['p_value'] is not None else "N/A"
        print(f"- **Comparison**: {wilcoxon_res['name_a']} vs {wilcoxon_res['name_b']}")
        print(f"- **Paired Instances**: {wilcoxon_res['pairs_count']}")
        print(f"- **Wins / Ties / Losses**: {wilcoxon_res['name_a']} won {wilcoxon_res['wins_a']} | Ties {wilcoxon_res['ties']} | {wilcoxon_res['name_b']} won {wilcoxon_res['wins_b']}")
        print(f"- **Wilcoxon W-statistic**: {w_str}")
        print(f"- **p-value**: {p_str}")
        print(f"- **Outcome (alpha = 0.05)**: **{wilcoxon_res['conclusion']}**")
    else:
        print("Paired Wilcoxon comparison could not be performed (insufficient overlapping instances).")

    print("\n## 3. Cactus Plot Data Export\n")
    print(f"- Cactus plot coordinates successfully exported to: `{cactus_out}`")
    print("- Ready for visualization with standard SAT Competition plotting scripts.")


def main():
    parser = argparse.ArgumentParser(
        description="SAT Benchmark Telemetry Analyzer & Statistical Validator."
    )
    parser.add_argument(
        "--input",
        default="logs/ablation_results.jsonl",
        help="Input JSONL telemetry file (default: logs/ablation_results.jsonl)",
    )
    parser.add_argument(
        "--cactus-out",
        default="logs/cactus_ablation.csv",
        help="Output CSV for Cactus plot (default: logs/cactus_ablation.csv)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=1800.0,
        help="Default timeout threshold for PAR-2 penalty (default: 1800.0)",
    )
    parser.add_argument("--candidate", help="Candidate condition id or named profile.")
    parser.add_argument("--baseline", help="Baseline condition id or named profile.")

    args = parser.parse_args()

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_path = os.path.join(repo_root, args.input) if not os.path.isabs(args.input) else args.input
    cactus_path = os.path.join(repo_root, args.cactus_out) if not os.path.isabs(args.cactus_out) else args.cactus_out

    records = load_records(input_path)
    if not records:
        print(f"Error: No valid records found in {input_path}", file=sys.stderr)
        sys.exit(1)

    summary, grouped_records = analyze_records(records, args.timeout)

    conditions = [s["condition"] for s in summary]
    export_cactus_csv(records, cactus_path, conditions)

    wilcoxon_res = None
    candidate = args.candidate
    baseline = args.baseline
    if candidate is None and baseline is None:
        if "proposed" in grouped_records and "proposed-no-decomposition" in grouped_records:
            candidate = "proposed"
            baseline = "proposed-no-decomposition"
        elif 3 in grouped_records and 0 in grouped_records:
            candidate = 3
            baseline = 0
    elif candidate is None or baseline is None:
        parser.error("--candidate and --baseline must be provided together")

    if isinstance(candidate, str) and candidate.isdigit() and int(candidate) in grouped_records:
        candidate = int(candidate)
    if isinstance(baseline, str) and baseline.isdigit() and int(baseline) in grouped_records:
        baseline = int(baseline)

    if candidate in grouped_records and baseline in grouped_records:
        wilcoxon_res = perform_wilcoxon_test(
            grouped_records[candidate],
            grouped_records[baseline],
            name_a=str(candidate),
            name_b=str(baseline),
            default_timeout=args.timeout,
        )

    print_markdown_report(summary, wilcoxon_res, cactus_path)


if __name__ == '__main__':
    main()

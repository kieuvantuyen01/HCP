#!/usr/bin/env python3
"""Analyze named HCP ablations with PAR-2, cactus data, and paired tests."""

import argparse
import csv
import json
import math
import os
import statistics
import sys
from typing import Dict, List, Optional, Tuple


REFERENCE = "full"


def load_records(path: str) -> List[Dict]:
    records = []
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
                record["condition"] = str(record["condition"])
                records.append(record)
            except (json.JSONDecodeError, KeyError) as error:
                print(f"Skipping line {line_number}: {error}", file=sys.stderr)
    return records


def solved(record: Dict) -> bool:
    return record.get("status") == "SATISFIABLE" and record.get("verified") is True


def par2_cost(record: Dict, default_timeout: float) -> float:
    timeout = float(record.get("timeout_limit", default_timeout))
    return min(float(record.get("wall_time", timeout)), timeout) if solved(record) else 2 * timeout


def compute_par2(records: List[Dict], default_timeout: float) -> float:
    return statistics.mean(par2_cost(record, default_timeout) for record in records)


def normal_approx_wilcoxon(x: List[float], y: List[float]) -> Tuple[float, float]:
    differences = [a - b for a, b in zip(x, y) if abs(a - b) > 1e-12]
    if not differences:
        return 0.0, 1.0
    ordered = sorted(enumerate(differences), key=lambda pair: abs(pair[1]))
    ranks = [0.0] * len(ordered)
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and abs(abs(ordered[end][1]) - abs(ordered[index][1])) < 1e-12:
            end += 1
        average_rank = (index + 1 + end) / 2
        for position in range(index, end):
            ranks[position] = average_rank
        index = end
    positive = sum(rank for rank, (_, diff) in zip(ranks, ordered) if diff > 0)
    negative = sum(rank for rank, (_, diff) in zip(ranks, ordered) if diff < 0)
    statistic = min(positive, negative)
    n = len(ordered)
    mean = n * (n + 1) / 4
    variance = n * (n + 1) * (2 * n + 1) / 24
    z_value = (statistic - mean) / math.sqrt(variance) if variance else 0.0
    p_value = math.erfc(abs(z_value) / math.sqrt(2))
    return statistic, min(1.0, max(0.0, p_value))


def paired_test(
    reference_records: List[Dict],
    ablation_records: List[Dict],
    default_timeout: float,
) -> Dict:
    ref = {(record["graph"], int(record["seed"])): record for record in reference_records}
    abl = {(record["graph"], int(record["seed"])): record for record in ablation_records}
    common = sorted(set(ref) & set(abl))
    ref_costs = [par2_cost(ref[key], default_timeout) for key in common]
    abl_costs = [par2_cost(abl[key], default_timeout) for key in common]
    if not common:
        return {"pairs": 0, "statistic": None, "p_raw": None,
                "reference_wins": 0, "ablation_wins": 0, "ties": 0}

    try:
        from scipy.stats import wilcoxon
        if all(abs(a - b) <= 1e-12 for a, b in zip(ref_costs, abl_costs)):
            statistic, p_value = 0.0, 1.0
        else:
            result = wilcoxon(ref_costs, abl_costs, zero_method="pratt",
                              alternative="two-sided")
            statistic, p_value = float(result.statistic), float(result.pvalue)
    except (ImportError, ValueError):
        statistic, p_value = normal_approx_wilcoxon(ref_costs, abl_costs)

    return {
        "pairs": len(common),
        "statistic": statistic,
        "p_raw": p_value,
        "reference_wins": sum(a < b - 1e-9 for a, b in zip(ref_costs, abl_costs)),
        "ablation_wins": sum(b < a - 1e-9 for a, b in zip(ref_costs, abl_costs)),
        "ties": sum(abs(a - b) <= 1e-9 for a, b in zip(ref_costs, abl_costs)),
    }


def holm_adjust(comparisons: Dict[str, Dict]) -> None:
    ranked = sorted(
        ((name, result["p_raw"]) for name, result in comparisons.items()
         if result["p_raw"] is not None),
        key=lambda pair: pair[1],
    )
    running_max = 0.0
    count = len(ranked)
    for rank, (name, p_value) in enumerate(ranked):
        adjusted = min(1.0, p_value * (count - rank))
        running_max = max(running_max, adjusted)
        comparisons[name]["p_holm"] = running_max
    for result in comparisons.values():
        result.setdefault("p_holm", None)


def summarize(records: List[Dict], timeout: float) -> Tuple[List[Dict], Dict[str, List[Dict]]]:
    grouped: Dict[str, List[Dict]] = {}
    for record in records:
        grouped.setdefault(record["condition"], []).append(record)
    rows = []
    for condition in sorted(grouped, key=lambda name: (name != REFERENCE, name)):
        subset = grouped[condition]
        solved_records = [record for record in subset if solved(record)]
        solved_times = [float(record["wall_time"]) for record in solved_records]
        rows.append({
            "condition": condition,
            "runs": len(subset),
            "graphs": len({record["graph"] for record in subset}),
            "seeds": len({int(record["seed"]) for record in subset}),
            "solved": len(solved_records),
            "timeouts": sum(record.get("status") == "TIMEOUT" for record in subset),
            "errors": sum(record.get("status") in ("ERROR", "UNKNOWN") for record in subset),
            "median": statistics.median(solved_times) if solved_times else math.nan,
            "par2": compute_par2(subset, timeout),
        })
    return rows, grouped


def export_cactus(records: List[Dict], path: str, conditions: List[str]) -> None:
    times = {
        condition: sorted(float(record["wall_time"]) for record in records
                          if record["condition"] == condition and solved(record))
        for condition in conditions
    }
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["solved_index"] + [f"{condition}_time" for condition in conditions])
        for index in range(max((len(values) for values in times.values()), default=0)):
            writer.writerow([index + 1] + [
                f"{times[condition][index]:.6f}" if index < len(times[condition]) else ""
                for condition in conditions
            ])


def print_report(rows: List[Dict], comparisons: Dict[str, Dict], cactus_path: str) -> None:
    print("# Named Ablation Study Report\n")
    print("| Condition | Graphs | Seeds | Runs | Solved | Timeouts | Errors | Median solved (s) | PAR-2 (s) |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        median = "N/A" if math.isnan(row["median"]) else f"{row['median']:.3f}"
        print(f"| {row['condition']} | {row['graphs']} | {row['seeds']} | {row['runs']} | "
              f"{row['solved']} | {row['timeouts']} | {row['errors']} | {median} | "
              f"{row['par2']:.3f} |")

    print("\n## Paired comparison against `full`\n")
    print("| Removed/changed condition | Pairs | Full wins | Ties | Ablation wins | Raw p | Holm p |")
    print("|---|---:|---:|---:|---:|---:|---:|")
    for condition in sorted(comparisons):
        result = comparisons[condition]
        raw = "N/A" if result["p_raw"] is None else f"{result['p_raw']:.4g}"
        holm = "N/A" if result["p_holm"] is None else f"{result['p_holm']:.4g}"
        print(f"| {condition} | {result['pairs']} | {result['reference_wins']} | "
              f"{result['ties']} | {result['ablation_wins']} | {raw} | {holm} |")
    print(f"\nCactus data: `{cactus_path}`")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="logs/ablation_results.jsonl")
    parser.add_argument("--cactus-out", default="logs/cactus_ablation.csv")
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--expected-seeds", type=int, default=5)
    args = parser.parse_args()

    if args.expected_seeds <= 0:
        parser.error("--expected-seeds must be positive")

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    resolve = lambda path: path if os.path.isabs(path) else os.path.join(repo_root, path)
    input_path = resolve(args.input)
    if not os.path.isfile(input_path):
        parser.error(f"input does not exist: {input_path}")
    records = load_records(input_path)
    if not records:
        parser.error("input contains no valid records")
    binary_hashes = {
        record.get("binary_sha256") for record in records if record.get("binary_sha256")
    }
    if len(binary_hashes) > 1:
        parser.error("input mixes multiple solver binaries; split the JSONL by binary hash")

    rows, grouped = summarize(records, args.timeout)
    incomplete = [
        row["condition"] for row in rows
        if row["seeds"] < args.expected_seeds
    ]
    if incomplete:
        print(
            f"Warning: fewer than {args.expected_seeds} seeds for: "
            + ", ".join(incomplete),
            file=sys.stderr,
        )
    conditions = [row["condition"] for row in rows]
    cactus_path = resolve(args.cactus_out)
    export_cactus(records, cactus_path, conditions)

    comparisons: Dict[str, Dict] = {}
    if REFERENCE in grouped:
        for condition, subset in grouped.items():
            if condition != REFERENCE:
                comparisons[condition] = paired_test(
                    grouped[REFERENCE], subset, args.timeout
                )
    holm_adjust(comparisons)
    print_report(rows, comparisons, cactus_path)


if __name__ == "__main__":
    main()

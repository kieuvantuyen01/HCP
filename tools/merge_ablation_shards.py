#!/usr/bin/env python3
"""Merge complete single-condition ablation shards into one audited matrix."""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple


CONDITIONS = ("full", "no-dispatch", "no-decomposition", "no-two-cut")
EXPECTED_FEATURES = {
    "full": (
        "name=full seed=1 dispatch=true series=true two_cut=true repair=true "
        "directed_cubic=true alternating_workers=3"
    ),
    "no-dispatch": (
        "name=no-dispatch seed=1 dispatch=false series=true two_cut=true "
        "repair=true directed_cubic=false alternating_workers=3"
    ),
    "no-decomposition": (
        "name=no-decomposition seed=1 dispatch=false series=false two_cut=false "
        "repair=true directed_cubic=false alternating_workers=3"
    ),
    "no-two-cut": (
        "name=no-two-cut seed=1 dispatch=true series=true two_cut=false "
        "repair=true directed_cubic=true alternating_workers=3"
    ),
}
FROZEN_KEYS = (
    "seed",
    "timeout_seconds",
    "memory_limit_mb",
    "git_commit",
    "binary_sha256",
    "dataset_manifest_sha256",
    "runlim_binary_sha256",
    "runlim_version",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_jsonl(path: Path) -> List[Dict]:
    records = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as error:
                raise ValueError(f"{path}:{line_number}: invalid JSON: {error}") from error
    return records


def graph_number(record: Dict) -> int:
    name = str(record.get("graph", ""))
    if not name.startswith("graph") or not name[5:].isdigit():
        raise ValueError(f"invalid graph name: {name!r}")
    return int(name[5:])


def locate_artifact(shard: Path, stored_path: object, subdirectory: str) -> Path:
    if not stored_path:
        raise ValueError(f"missing {subdirectory} artifact path in {shard}")
    candidate = shard / subdirectory / Path(str(stored_path)).name
    if not candidate.is_file():
        raise ValueError(f"missing artifact: {candidate}")
    return candidate


def validate_shard(shard: Path) -> Tuple[Dict, List[Dict]]:
    metadata_path = shard / "metadata.json"
    results_path = shard / "results.jsonl"
    manifest_path = shard / "instances.sha256"
    for path in (metadata_path, results_path, manifest_path):
        if not path.is_file():
            raise ValueError(f"required shard file is missing: {path}")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    for field in ("memory_limit_mb", "runlim_binary_sha256"):
        if metadata.get(field) is None:
            raise ValueError(f"shard metadata is missing {field}: {shard}")
    records = read_jsonl(results_path)
    if not records:
        raise ValueError(f"shard has no records: {shard}")
    if metadata.get("postflight_errors"):
        raise ValueError(f"shard reports postflight errors: {shard}")
    if sha256_file(manifest_path) != metadata.get("dataset_manifest_sha256"):
        raise ValueError(f"dataset manifest hash mismatch: {shard}")

    selected_graph_numbers = metadata.get("selected_graph_numbers")
    if selected_graph_numbers is None:
        selected_graph_numbers = list(range(1, 1002))
    if (
        not isinstance(selected_graph_numbers, list)
        or not selected_graph_numbers
        or any(not isinstance(number, int) for number in selected_graph_numbers)
        or len(selected_graph_numbers) != len(set(selected_graph_numbers))
        or any(number < 1 or number > 1001 for number in selected_graph_numbers)
    ):
        raise ValueError(f"invalid selected_graph_numbers in {shard}")
    selected_graph_numbers = set(selected_graph_numbers)

    keys = []
    condition_counts = Counter()
    for record in records:
        number = graph_number(record)
        condition = str(record.get("condition"))
        seed = int(record.get("seed", -1))
        key = (number, condition, seed)
        keys.append(key)
        condition_counts[condition] += 1
        if number < 1 or number > 1001:
            raise ValueError(f"graph outside graph1..graph1001 in {shard}: {key}")
        if condition not in CONDITIONS:
            raise ValueError(f"unexpected condition in {shard}: {condition}")
        if seed != metadata.get("seed"):
            raise ValueError(f"seed mismatch in {shard}: {key}")
        if float(record.get("timeout_limit", -1)) != float(
            metadata.get("timeout_seconds", -2)
        ):
            raise ValueError(f"timeout mismatch in {shard}: {key}")
        if int(record.get("memory_limit_mb", -1)) != int(
            metadata.get("memory_limit_mb", -2)
        ):
            raise ValueError(f"memory-limit mismatch in {shard}: {key}")
        if record.get("peak_memory_mb") is None:
            raise ValueError(f"missing peak-memory telemetry in {shard}: {key}")
        if not record.get("runlim_status"):
            raise ValueError(f"missing final runlim status in {shard}: {key}")
        if record.get("git_commit") != metadata.get("git_commit"):
            raise ValueError(f"Git revision mismatch in {shard}: {key}")
        if record.get("binary_sha256") != metadata.get("binary_sha256"):
            raise ValueError(f"binary hash mismatch in {shard}: {key}")
        if record.get("feature_vector") != EXPECTED_FEATURES[condition]:
            raise ValueError(f"feature-vector mismatch in {shard}: {key}")
        if record.get("status") in ("ERROR", "UNKNOWN", "UNSATISFIABLE"):
            raise ValueError(f"invalid final status in {shard}: {key}")
        if record.get("status") == "SATISFIABLE" and record.get("verified") is not True:
            raise ValueError(f"unverified SAT result in {shard}: {key}")
        locate_artifact(shard, record.get("stdout_path"), "raw-logs")
        locate_artifact(shard, record.get("stderr_path"), "raw-logs")
        locate_artifact(shard, record.get("runlim_log_path"), "runlim")
        if record.get("status") == "SATISFIABLE":
            locate_artifact(shard, record.get("tour_path"), "tours")

    if len(keys) != len(set(keys)):
        raise ValueError(f"duplicate graph-condition-seed key within {shard}")
    for condition, count in condition_counts.items():
        expected = {
            (number, condition, metadata["seed"])
            for number in selected_graph_numbers
        }
        actual = {key for key in keys if key[1] == condition}
        if count != len(selected_graph_numbers) or actual != expected:
            raise ValueError(
                f"condition {condition} does not match its declared graph shard "
                f"in {shard}: {count}/{len(selected_graph_numbers)} records"
            )
    return metadata, records


def copy_record_artifacts(record: Dict, shard: Path, stage: Path, final: Path) -> Dict:
    merged = dict(record)
    for field, subdirectory in (
        ("stdout_path", "raw-logs"),
        ("stderr_path", "raw-logs"),
        ("runlim_log_path", "runlim"),
        ("tour_path", "tours"),
    ):
        if not record.get(field):
            merged[field] = None
            continue
        source = locate_artifact(shard, record[field], subdirectory)
        relative = Path(subdirectory) / source.name
        destination = stage / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if sha256_file(destination) != sha256_file(source):
                raise ValueError(f"artifact-name collision with different content: {relative}")
        else:
            shutil.copy2(source, destination)
        merged[field] = str((final / relative).resolve())
    merged["source_shard"] = str(shard.resolve())
    return merged


def validate_condition_coverage(loaded: List[Tuple[Path, Dict, List[Dict]]]) -> None:
    """Require each condition represented by the supplied shards to be complete."""
    present_conditions = {
        str(record["condition"])
        for _, _, records in loaded
        for record in records
    }
    expected_numbers = set(range(1, 1002))
    for condition in present_conditions:
        actual_numbers = {
            graph_number(record)
            for _, _, records in loaded
            for record in records
            if record["condition"] == condition
        }
        if actual_numbers != expected_numbers:
            raise ValueError(
                f"condition {condition} is incomplete across shards: "
                f"missing={sorted(expected_numbers - actual_numbers)}, "
                f"unexpected={sorted(actual_numbers - expected_numbers)}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shards", nargs="+", help="Downloaded shard directories")
    parser.add_argument("--output-dir", required=True, help="New merged output directory")
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="Merge complete condition shards before all four conditions are available",
    )
    args = parser.parse_args()

    shards = [Path(path).expanduser().resolve() for path in args.shards]
    output = Path(args.output_dir).expanduser().resolve()
    if output.exists():
        parser.error(f"output already exists: {output}")
    if len(set(shards)) != len(shards):
        parser.error("the same shard directory was supplied more than once")

    try:
        loaded = [(shard, *validate_shard(shard)) for shard in shards]
        reference = loaded[0][1]
        for shard, metadata, _ in loaded[1:]:
            mismatches = [
                key for key in FROZEN_KEYS if metadata.get(key) != reference.get(key)
            ]
            if mismatches:
                raise ValueError(
                    f"incompatible shard {shard}; differing frozen fields: "
                    + ", ".join(mismatches)
                )

        all_keys = []
        present_conditions = set()
        for _, _, records in loaded:
            for record in records:
                key = (record["graph"], record["condition"], int(record["seed"]))
                all_keys.append(key)
                present_conditions.add(record["condition"])
        if len(all_keys) != len(set(all_keys)):
            raise ValueError("duplicate graph-condition-seed key across shards")
        validate_condition_coverage(loaded)
        if not args.allow_incomplete and present_conditions != set(CONDITIONS):
            missing = sorted(set(CONDITIONS) - present_conditions)
            raise ValueError("final matrix is missing conditions: " + ", ".join(missing))

        output.parent.mkdir(parents=True, exist_ok=True)
        stage = Path(tempfile.mkdtemp(prefix=f".{output.name}-", dir=output.parent))
        try:
            merged_records = []
            sources = []
            for shard, metadata, records in loaded:
                sources.append({
                    "path": str(shard),
                    "machine": metadata.get("machine"),
                    "platform": metadata.get("platform"),
                    "processor": metadata.get("processor"),
                    "logical_cpus": metadata.get("logical_cpus"),
                    "memory_bytes": metadata.get("memory_bytes"),
                    "memory_limit_mb": metadata.get("memory_limit_mb"),
                    "runlim_binary_sha256": metadata.get("runlim_binary_sha256"),
                    "runlim_version": metadata.get("runlim_version"),
                    "conditions": sorted({record["condition"] for record in records}),
                    "shard_count": metadata.get("shard_count", 1),
                    "shard_index": metadata.get("shard_index", 0),
                    "selected_graph_numbers": metadata.get(
                        "selected_graph_numbers", list(range(1, 1002))
                    ),
                    "records": len(records),
                })
                merged_records.extend(
                    copy_record_artifacts(record, shard, stage, output)
                    for record in records
                )

            condition_order = {name: index for index, name in enumerate(CONDITIONS)}
            merged_records.sort(
                key=lambda record: (
                    condition_order[record["condition"]],
                    graph_number(record),
                    int(record["seed"]),
                )
            )
            with (stage / "results.jsonl").open("w", encoding="utf-8") as handle:
                for record in merged_records:
                    handle.write(json.dumps(record, sort_keys=True) + "\n")
            shutil.copy2(shards[0] / "instances.sha256", stage / "instances.sha256")

            metadata = {
                "protocol": "final-single-run-equivalent-machines",
                "protocol_document": reference.get("protocol_document"),
                "merged_at_utc": datetime.now(timezone.utc).isoformat(),
                "instances": 1001,
                "conditions": [name for name in CONDITIONS if name in present_conditions],
                "completed_conditions": [
                    name for name in CONDITIONS if name in present_conditions
                ],
                "records": len(merged_records),
                "runs_expected": 4004,
                "status_counts": dict(Counter(r["status"] for r in merged_records)),
                "sources": sources,
            }
            for key in FROZEN_KEYS:
                metadata[key] = reference.get(key)
            (stage / "metadata.json").write_text(
                json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
            )
            os.replace(stage, output)
        except Exception:
            shutil.rmtree(stage, ignore_errors=True)
            raise

        analyzer = Path(__file__).resolve().with_name("analyze_ablation_results.py")
        report = subprocess.run(
            [
                sys.executable,
                str(analyzer),
                "--input",
                str(output / "results.jsonl"),
                "--timeout",
                str(reference["timeout_seconds"]),
                "--expected-seeds",
                "1",
                "--cactus-out",
                str(output / "cactus.csv"),
            ],
            check=True,
            text=True,
            capture_output=True,
        )
        (output / "report.md").write_text(report.stdout, encoding="utf-8")
        if report.stderr:
            print(report.stderr, file=sys.stderr, end="")
        print(
            f"Merged {len(merged_records)} records from {len(shards)} shard(s) "
            f"into {output}"
        )
    except (KeyError, OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()

# Final single-run ablation protocol on equivalent machines

Date frozen: 2026-09-28

## Scope

This is the resource-constrained final ablation matrix for the AAMAS 2027
manuscript. It intentionally replaces the five-seed plan only for this run. Each
of the 1,001 FHCP instances is executed once under each frozen condition:

| Condition | Interpretation |
|---|---|
| `full` | Complete proposed solver |
| `no-dispatch` | Disable the hub-cluster, alternating-pair, and directed-cubic front ends |
| `no-decomposition` | After dispatch is disabled, also disable forced-chain contraction and 2-cut decomposition |
| `no-two-cut` | Disable only separation-pair extraction and subpath stitching |

The first three conditions form the nested chain `full -> no-dispatch ->
no-decomposition`. The `full` versus `no-two-cut` contrast measures the 2-cut
stage directly. Other implemented ablations remain exploratory and are not part
of the final matrix.

The 1,001 `full` records already exist as two nonoverlapping partitions:
`benchmark_machine1` covers graphs 1--700 and `benchmark_machine2` covers
graphs 701--1001. The experiment owner confirms that they used the same solver
configuration and equivalent VM environments. Their merge is therefore frozen
as the single `full` condition and is not rerun. It serves both as the reference
for every ablation contrast and as the proposed solver's coverage row in the
related-work comparison. Archived elapsed times from external solvers are not
combined with this matrix for runtime or PAR-2 claims.

## Frozen controls

- Dataset: all `graph1.col` through `graph1001.col` in `FHCPCS-col`.
- Cutoff: 1,800 seconds end-to-end wall time for every run.
- Memory cap: 7,168 MB of process-tree space per run, enforced by
  `runlim` 2.0.0rc14, matching the accepted `full` artifact. Physical VM RAM is not the
  experimental limit. A run terminated by this cap is recorded as `MEMOUT` and
  remains unsolved for solved-count and PAR-2 calculations.
- Seed: the three new ablation jobs use CaDiCaL seed 1. The accepted legacy
  `full` artifact records the default invocation but not a seed field.
- Replication: one run per instance-condition pair, 4,004 runs total.
- Scheduling: `no-dispatch` and `no-decomposition` each run as one complete
  1,001-instance job. To reduce elapsed time, `no-two-cut` may be split across
  two equivalent VMs using a frozen round-robin partition: shard 1 contains
  odd-numbered graphs and shard 2 contains even-numbered graphs. The shards
  are disjoint and their union is exactly `graph1` through `graph1001`. Only
  one condition shard runs on a given VM at a time.
- Solver parallelism: at most three internal workers, with
  `RAYON_NUM_THREADS=3` and `OMP_NUM_THREADS=1`.
- Artifact: for each new ablation job, use one release binary from one clean Git
  revision and record the Git hash, binary SHA-256, dataset manifest, VM
  metadata, compiler versions, command line, raw stdout/stderr, wall time,
  aggregate child CPU time, peak memory, `runlim` binary hash and version,
  raw `runlim` output, and resolved
  feature vector. The accepted `full`
  CSV predates this metadata schema and does not record the Git revision or
  binary SHA-256; retain that fact as a provenance limitation. Machine identity
  is provenance, not a merge key.
- Validation: independently verify every reported Hamiltonian cycle against the
  original `.col` graph. An invalid or missing tour is an error, not a solve.
- Resume: retain completed records only when graph path, condition, seed, cutoff,
  memory cap, and binary hash all match.

## Measures and interpretation

Primary measures are verified solved count and PAR-2, with timeouts charged at
3,600 seconds. Secondary outputs are a cactus curve, solved-instance median,
and paired win/tie/loss counts. Pairwise tests use matched instances and Holm
correction, but one run cannot estimate run-to-run or seed variability. Claims
therefore apply only to the frozen seed-1 experiment. Prefer solved-count
differences and large PAR-2 effects over marginal runtime differences.

The coarse `no-decomposition` condition cannot identify the contribution of an
individual decomposition mechanism. Interpret `no-dispatch` versus
`no-decomposition` as the combined contribution of series contraction and
2-cut substitution in the general pipeline. Use `full` versus `no-two-cut` for
the direct 2-cut claim.

## Validity checks

1. Stop interpretation if a condition does not emit its expected resolved
   feature vector, any reported tour fails independent verification, or any run
   ends as `ERROR`/`UNKNOWN`. Treat an `UNSATISFIABLE` result as an error because
   every FHCP Challenge instance in this matrix is known to be Hamiltonian.
   Also stop if any new ablation record lacks the 7,168-MB cap, peak-memory
   telemetry, or its raw `runlim` record.
2. Require exactly one unique record for every graph-condition pair and exactly
   4,004 records before producing the final table.
3. Condition shards from different equivalent VMs may be merged. A condition
   may consist of multiple disjoint graph shards only when their audited union
   has exactly one record for every graph from 1 through 1001. For the three
   new ablation jobs, require the same binary SHA-256, Git revision,
   dataset-manifest SHA-256, cutoff, and seed. The already accepted `full` merge
   is the documented legacy exception; verify its exact 1--700/701--1001
   partition, configuration, limits, outcomes, and tours instead of inventing
   missing hashes.
4. Report timeouts as censored outcomes through solved count and PAR-2. Do not
   treat a timeout as evidence of unsatisfiability.
5. Report the single-run limitation in the manuscript and make no claim of seed
   robustness.

## Execution

The `full` condition is already complete. Do **not** run
`run_ablation_full_once.sh` again. From the repository root, run only the three
remaining entrypoints, one on each assigned VM. A VM may also run multiple
entrypoints sequentially; do not overlap two jobs on the same VM.

The entrypoints require `runlim`; they abort before building or resuming if it
is unavailable. Resume keys include the 7,168-MB cap, so results produced by an
older unlimited-memory harness cannot be reused:

```bash
./scripts/run_ablation_no_dispatch_once.sh
./scripts/run_ablation_no_decomposition_once.sh
./scripts/run_ablation_no_two_cut_part1_once.sh  # VM 1: odd graphs
./scripts/run_ablation_no_two_cut_part2_once.sh  # VM 2: even graphs
```

Re-running an interrupted condition safely resumes it when its frozen artifact
and metadata match. A lock prevents two scripts from sharing CPU and writing the
same shard concurrently. Outputs are written to
`benchmark-runs/ablation-final-once-1001/` on each VM.

After downloading each VM's complete output directory under a distinct local
name, merge and audit the three new-format shards:

```bash
python3 tools/merge_ablation_shards.py \
  --allow-incomplete \
  --output-dir benchmark-runs/ablation-final-once-1001-merged \
  /path/to/no-dispatch /path/to/no-decomposition \
  /path/to/no-two-cut-part1 /path/to/no-two-cut-part2
```

The merger rejects an incomplete supplied condition, duplicate
graph--condition keys, and any mismatch in the frozen fields for these new
jobs. It copies tours, raw logs, and `runlim` records, rewrites their paths, and regenerates
`report.md` and `cactus.csv`. The accepted legacy `full` CSV is joined with this
three-condition result at the analysis layer; it is not a new-format shard and
must not be passed directly to this merger.

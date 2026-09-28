# Final single-run ablation protocol

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

## Frozen controls

- Dataset: all `graph1.col` through `graph1001.col` in `FHCPCS-col`.
- Cutoff: 1,800 seconds end-to-end wall time for every run.
- Seed: CaDiCaL seed 1.
- Replication: one run per instance-condition pair, 4,004 runs total.
- Scheduling: four condition jobs on the same machine, run sequentially without
  overlap. Each job traverses all 1,001 instances in natural numeric order.
- Solver parallelism: at most three internal workers, with
  `RAYON_NUM_THREADS=3` and `OMP_NUM_THREADS=1`.
- Artifact: one release binary from one clean Git revision. Record the Git hash,
  binary SHA-256, dataset manifest, platform, compiler versions, command line,
  raw stdout/stderr, wall time, aggregate child CPU time, and resolved feature
  vector.
- Validation: independently verify every reported Hamiltonian cycle against the
  original `.col` graph. An invalid or missing tour is an error, not a solve.
- Resume: retain completed records only when graph path, condition, seed, cutoff,
  and binary hash all match.

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
2. Require exactly one unique record for every graph-condition pair and exactly
   4,004 records before producing the final table.
3. Do not merge results from different binary hashes, Git revisions, cutoffs,
   seeds, or machines into this matrix.
4. Report timeouts as censored outcomes through solved count and PAR-2. Do not
   treat a timeout as evidence of unsatisfiability.
5. Report the single-run limitation in the manuscript and make no claim of seed
   robustness.

## Execution

From the repository root, run these four entrypoints one at a time. They append
to one binary-checked result matrix and regenerate the report after every
completed condition:

```bash
./scripts/run_ablation_full_once.sh
./scripts/run_ablation_no_dispatch_once.sh
./scripts/run_ablation_no_decomposition_once.sh
./scripts/run_ablation_no_two_cut_once.sh
```

Re-running an interrupted condition safely resumes it when its frozen artifact
and metadata match. A lock prevents two condition scripts from sharing CPU and
writing the result matrix concurrently. Outputs are written to
`benchmark-runs/ablation-final-once-1001/`.

# Proposed–NoDecomposition benchmark profile

Date: 2026-09-23

## Purpose

This is the matched monolithic baseline for measuring the contribution of
structural decomposition. `Proposed` and `Proposed–NoDecomposition` use the
same `cegar-fix` binary and shared global CEGAR components: topology features,
global encoding selection, CaDiCaL backend, static cuts, CEGAR loop, cycle-repair
operators, external cutoff, and independent tour verifier. The sole configured
switch is `--structural-decomposition`; local/decomposed formulations are part
of the treatment and are consequently absent from the monolithic profile.

With the factor disabled, the solver:

- skips the B1 two-tier route and uses its global CEGAR fallback;
- retains the same degree-2 chain contraction as `Proposed`;
- skips inverse-3SAT, hub-hierarchical, ring-quotient, module-interface,
  modular-macro, and hub-partition decomposition paths.

Graph parsing, safe degree-2 triangle pruning, degree-2 contraction,
topology-based encoding selection, subtour cuts, and cycle repair remain
enabled in both profiles.

The exact arguments are registered in
`benchmarks/proposed_variants.json`. Each raw log prints either
`Structural decomposition: enabled` or `Structural decomposition: disabled`.

## Commands

Run the no-decomposition profile over all 1,001 graphs with a 1,800-second
cutoff:

```bash
./run_proposed_no_decomposition.sh
```

Resume the same run without rebuilding:

```bash
SKIP_BUILD=1 ./run_proposed_no_decomposition.sh --resume
```

Run the paired profiles together, which is preferred when neither sweep has
started:

```bash
python3 tools/benchmark_runner_ablation.py \
  --binary src/cegar-fix/target/release/cegar-fix \
  --graphs-dir FHCPCS-col \
  --profiles proposed,proposed-no-decomposition \
  --profile-config benchmarks/proposed_variants.json \
  --seeds 1 \
  --timeout 1800 \
  --expected-status SATISFIABLE \
  --cores 0 \
  --output benchmark-runs/proposed-pair/results.jsonl \
  --tour-dir benchmark-runs/proposed-pair/tours \
  --raw-dir benchmark-runs/proposed-pair/raw
```

The runner writes a provenance manifest beside `results.jsonl`, retains stdout
and stderr for every run, applies a hard process timeout, and verifies every SAT
tour against the original `.col` graph. Because every FHCP Challenge instance is
known Hamiltonian, an `UNSATISFIABLE` line is recorded as
`INVALID_UNSATISFIABLE`.

Generate the PAR-2 summary, paired comparison, and cactus data with:

```bash
python3 tools/analyze_ablation_results.py \
  --input benchmark-runs/proposed-pair/results.jsonl \
  --candidate proposed \
  --baseline proposed-no-decomposition \
  --timeout 1800 \
  --cactus-out benchmark-runs/proposed-pair/cactus.csv \
  > benchmark-runs/proposed-pair/analysis.md
```

The current binary does not expose a global SAT seed. Therefore `--seeds` is a
repetition label, not a claim of independent random seeds. Use one label for the
full sweep and treat repeated labels only as runtime-noise measurements until a
seed is wired through all SAT and randomized repair components.

## Pre-registered comparison

Primary metrics are verified solved count on 1,001 instances, verified solved
count on the residual 29, and PAR-2 at 1,800 seconds. Secondary metrics are
paired wins/losses, cactus curves, CEGAR increments, blocking clauses, and unique
verified solves.

The decomposition contribution is not supported if `Proposed` solves fewer
instances than `Proposed–NoDecomposition`, or if its full-set PAR-2 is worse.
A general efficiency claim is also rejected if `Proposed` regresses by more than
5% on the pre-registered easy/medium subset. Any SAT result without a verified
tour, and any UNSAT result on this known-Hamiltonian suite, is invalid rather
than a solve.

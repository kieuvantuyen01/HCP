# Proposed component-ablation profiles

## Experimental contract

The named profiles in `benchmarks/proposed_variants.json` compare the complete
`Proposed` configuration with three single-factor ablations. Every profile uses
the same solver binary, automatic topology classifier, SAT encoding policy,
CaDiCaL backend, refinement cuts, resource limit, and independent tour
verification.

| Profile | Structural decomposition | Degree-2 contraction | Cycle repair |
|---|---:|---:|---:|
| `proposed` | on | on | on |
| `proposed-no-decomposition` | off | on | on |
| `proposed-no-contraction` | on | off | on |
| `proposed-no-repair` | on | on | off |

`Proposed–NoContraction` passes the uncontracted graph and an identity
contractor to the same downstream solver. It keeps decomposition and repair
enabled. This isolates the effect of degree-2 chain contraction on reduced
graph size, CNF size, SAT conflicts, and runtime.

`Proposed–NoRepair` disables all cycle-cover-changing operators in the global
CEGAR path, including alternating-port and corridor repair, patchers,
absorbers, 2-opt/3-opt, macro crossover, and localized SAT repair. It also
disables two-tier splicer patching. The solver still rejects disconnected
cycle covers by adding valid refinement clauses. This isolates repair from the
SAT–CEGAR search rather than disabling connectivity refinement.

## Full runs

Run each new profile over all 1,001 FHCP instances with a 1,800-second limit:

```bash
./run_proposed_no_contraction.sh
./run_proposed_no_repair.sh
```

Resume an interrupted run without duplicating completed rows:

```bash
SKIP_BUILD=1 ./run_proposed_no_contraction.sh --resume
SKIP_BUILD=1 ./run_proposed_no_repair.sh --resume
```

For a strictly paired run using one binary and one manifest, run all four
profiles together:

```bash
python3 tools/benchmark_runner_ablation.py \
  --binary src/cegar-fix/target/release/cegar-fix \
  --graphs-dir FHCPCS-col \
  --profiles proposed,proposed-no-decomposition,proposed-no-contraction,proposed-no-repair \
  --profile-config benchmarks/proposed_variants.json \
  --seeds 1 \
  --timeout 1800 \
  --expected-status SATISFIABLE \
  --cores 0 \
  --output benchmark-runs/proposed-component-ablations/results.jsonl \
  --tour-dir benchmark-runs/proposed-component-ablations/tours
```

The current executable does not expose a globally controlled random seed.
Consequently, `--seeds` records repetition labels; it must not be described as
controlled random seeding. Use one cold-start sweep for the primary comparison
and report repeated runs separately as process-level robustness checks.

## Analysis

Analyze each paired contrast separately:

```bash
python3 tools/analyze_ablation_results.py \
  --input benchmark-runs/proposed-component-ablations/results.jsonl \
  --candidate proposed \
  --baseline proposed-no-contraction \
  --timeout 1800 \
  --cactus-out benchmark-runs/proposed-component-ablations/cactus-no-contraction.csv

python3 tools/analyze_ablation_results.py \
  --input benchmark-runs/proposed-component-ablations/results.jsonl \
  --candidate proposed \
  --baseline proposed-no-repair \
  --timeout 1800 \
  --cactus-out benchmark-runs/proposed-component-ablations/cactus-no-repair.csv
```

Count only independently verified tours as solved. An invalid SAT result or an
UNSAT result on a known Hamiltonian FHCP instance is an error, never a solve.
Report solved count, PAR-2, cactus curves, paired wins/losses, and timeout-aware
paired statistics. For mechanism evidence, also report contracted vertex/CNF
size for `NoContraction` and the number of covers and repair successes for
`NoRepair` once those telemetry counters are available.

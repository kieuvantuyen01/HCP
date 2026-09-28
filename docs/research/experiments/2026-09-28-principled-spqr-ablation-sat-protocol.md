# Principled SPQR SAT ablation protocol

Date frozen: 2026-09-28

## Question and independent variable

The experiment measures which mechanisms in the proposed solver account for its
performance on the 1,001-instance FHCP suite.  `full` is the reference condition.
Every named ablation changes one mechanism, except `no-decomposition`, which is a
deliberately coarse compatibility condition corresponding to the earlier
Proposed--NoDecomposition label.

| Condition | Change relative to `full` |
|---|---|
| `full` | All implemented mechanisms enabled |
| `no-dispatch` | Disable hub-cluster, alternating-pair, and directed-cubic dispatch |
| `no-series` | Disable forced degree-2 chain contraction; the dependent alternating-pair frontend is also unavailable |
| `no-two-cut` | Disable separation-pair extraction and subpath stitching |
| `no-repair` | Disable all implemented 2-opt/3-opt cycle/path repair steps |
| `undirected-cubic` | Route cubic graphs through the undirected engine instead of directed CEGAR |
| `one-alternating-worker` | Use one alternating-pair SAT worker instead of three |
| `no-decomposition` | Disable tailored dispatch, series contraction, and 2-cut decomposition |

The aliases `proposed-nodecomposition`, `proposed-norepair`, and
`proposed-nocontraction` map to `no-decomposition`, `no-repair`, and `no-series`.

## Fixed controls

- Seeds: `1,42,137,777,2026`; the runner passes each seed to CaDiCaL.
- Cutoff tiers: 60 s screening, 300 s confirmation, and 1,800 s final.
- Hardware, compiler profile, CPU affinity, solver binary, and Git commit must be
  identical within each paired comparison.
- Every reported SAT result must include a tour independently checked against the
  original, uncontracted graph.
- Primary metrics: solved count and PAR-2. Secondary metrics: cactus curves,
  median solved runtime, CEGAR iterations, and pairwise win/tie/loss counts.
- Statistical test: paired two-sided Wilcoxon signed-rank test over PAR-2 costs,
  with Holm correction across ablations and alpha 0.05. Report effect sizes and
  raw counts even when significance is not reached.

## Baselines

Run the same five seeds and cutoffs for: (1) the original Takehide Soh SAT-CEGAR
implementation, (2) the best CEGAR 2025 configuration, (3) the best configuration
from *An UPgrade for HCP Solving the Hamiltonian Cycle Problem via User
Propagators*, (4) a plain edge-degree CNF with CaDiCaL 1.9.5, and (5) the same
plain CNF with Kissat. Record exact solver versions, revisions, seeds, and command
lines. The two plain-CNF runs isolate CDCL-engine effects from the proposed
decomposition and repair mechanisms.

## Pre-registered failure criteria

These criteria are fixed before executing the smoke tests or benchmark runs:

1. Any invalid reported tour, crash, or condition that silently executes the same
   feature set as `full` is an implementation failure and invalidates that run.
2. A mechanism is not claimed to help unless removing it worsens either solved
   count or PAR-2 at 1,800 s, the direction is consistent in at least four of five
   seeds, and the paired effect survives Holm correction. Otherwise the result is
   reported as inconclusive or neutral.
3. If `full` solves over 5% fewer instances than the frozen pre-ablation full
   binary at any confirmation/final tier, stop and investigate the regression
   before interpreting ablations.
4. If a baseline or ablation is missing any of the five seeds, compare it only in
   an explicitly labeled incomplete analysis; do not use it for the main claim.
5. A timeout is charged at twice the cutoff for PAR-2. Solver-reported UNSAT on a
   challenge instance is checked independently before being treated as proven.

## Execution order

First run a small topology-stratified smoke set for all conditions and verify that
the log reports the intended feature vector. Then run the 60 s tier on all 1,001
instances, promote informative and unresolved cases to 300 s, and run the final
1,800 s comparison on the preregistered full suite or a clearly frozen residual
set. Preserve JSONL logs, tours, binary hash, Git commit, and machine metadata.

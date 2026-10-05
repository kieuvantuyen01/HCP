# StructHCP Minus Dispatch: Deadline Single-Seed Protocol

**Date frozen:** 2026-10-05

## Hypothesis

Disabling graph-family dispatch while retaining series contraction, repair, and
the generic SAT--CEGAR route reduces verified full-suite coverage relative to
the retained StructHCP result.

## Conditions and controls

- Reference: retained `no-two-cut` StructHCP artifact, 826/1,001 verified tours.
- New condition: `no-dispatch-no-two-cut` with `dispatch=false`, `series=true`,
  `two_cut=false`, `repair=true`, and `directed_cubic=false`.
- Dataset: all 1,001 FHCP Challenge instances.
- Seed: 1. No seed-robustness claim will be made.
- Cutoff: 1,800 seconds end-to-end wall time.
- Process-tree memory limit: 7,168 MB under runlim 2.0.0rc14.
- CPU affinity: logical CPUs 1--3; LKH remains isolated on CPU 0.
- Sharding: odd graph numbers on shard 1 and even graph numbers on shard 2.
  The two shards must share the binary hash, Git revision, dataset manifest,
  seed, cutoff, memory limit, and feature vector.

## Measures

Verified solved count and PAR-2 are primary. Also report timeout and memory-out
counts, cactus data, peak memory, CEGAR iterations, the paired solved-set
contingency, exact McNemar, and paired per-instance PAR-2 statistics. Audit
retention of the fixed 34 complementary certificates and the finalized 12.

## Pre-registered interpretation and failure thresholds

1. Any invalid or missing tour, wrong feature vector, `ERROR`, `UNKNOWN`, shard
   overlap/gap, or provenance mismatch invalidates the new condition.
2. Dispatch receives a coverage-benefit interpretation only if the new
   condition solves fewer instances than the retained StructHCP artifact and
   loses at least one of the fixed 34 complementary certificates. Otherwise
   the paper reports no supported aggregate coverage benefit from dispatch.
3. A difference of fewer than five solved instances is treated as practically
   marginal unless supported by the paired solved-set test and a material
   PAR-2 change.
4. Because only seed 1 is run, all conclusions are explicitly restricted to
   this frozen experiment; no robustness-across-seeds claim is permitted.
5. HCup 2026 and cut-set CEGAR 2025 remain the external coverage baselines. The
   existing `full` configuration is the already completed `StructHCP+2Cut`
   comparison and is not rerun.

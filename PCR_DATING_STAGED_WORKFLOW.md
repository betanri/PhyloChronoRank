# PCR Dating Staged Workflow

This is the default workflow for dating runs that will feed PCR / postfit comparisons.

## Runner

```bash
python3 scripts/run_dating_grid_staged.py \
  --phylogram=PATH/TO/phylogram.tre \
  --calibrations-csv=PATH/TO/calibrations.csv \
  --outdir=PATH/TO/dating_out \
  --out-prefix=my_dataset \
  --treepl-numsites=ALIGNMENT_LENGTH
```

Use `--reference-time-tree=...` instead of `--calibrations-csv=...` for congruification. Use `--start-stage=N` to resume at stage `N`.

The runner calls `scripts/run_dating_grid.R` once per stage. All stages write to the same `--outdir`.

## Stage order

1. `RelTime_MEGA` (MEGA-CC): point estimate, Tao analytical CI, and bootstrap CI (`megacc` rerun per replicate).
2. `treePL`: prime, cross-validation over `--treepl-smoothing`, then one `thorough` dating run with the CV-selected smoothing.
3. `chronos`: grid winners for all four models (`clock`, `correlated`, `discrete`, `relaxed`). Starting ages are seeded from the stage-2 treePL tree (the RelTime tree if treePL failed).
4. `treePL` bootstrap CI for the stage-2 tree.
5. `chronos` bootstrap CI for `clock`, `correlated`, `discrete` (`--chronos-ci-models`). `relaxed` stays plain.

Stages 4 and 5 reuse the stage-2 and stage-3 trees from their checkpoints. They add CIs to those exact trees and do not refit them.

## MAIN_OUTPUT_TREES policy

- `MAIN_OUTPUT_TREES` is refreshed at the end of every stage.
- Each selected candidate appears in its best available form: `<candidate>_with_CI.tre` when its CI exists, else `<candidate>.tre`.
- After stage 3 the folder holds a usable PCR candidate set:
  - `RelTime_MEGA_with_bootstrap_CI.tre` and `RelTime_MEGA_with_tao_CI.tre`
  - plain `treePL_smooth<value>.tre`
  - plain `chronos_<model>_lambda<value>[_k<k>].tre` winners
- Stage 4 replaces the plain treePL tree with `treePL_smooth<value>_with_CI.tre`.
- Stage 5 replaces the plain `chronos_clock`, `chronos_correlated`, and `chronos_discrete` trees with their `_with_CI.tre` versions; `chronos_relaxed` stays plain.
- For large trees, subset-tuning trees never go into `MAIN_OUTPUT_TREES`; only full-tree chronos reruns and full-tree treePL/RelTime trees do.
- `MAIN_OUTPUT_TREES/ci/` holds the matching CI tables.
- `candidates.csv` in `--outdir` points to the `MAIN_OUTPUT_TREES` files and is ready for `scripts/run_pcr.R`.

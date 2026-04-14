# PCR Dating Staged Workflow

This is the default workflow for future dating runs that will feed PCR / postfit comparisons.

## Canonical stage order

1. `RelTime MEGA` first, with both CI products.
2. `treePL` second, plain selected winner only.
3. `chronos` third, plain selected winners only for all four models:
   - `clock`
   - `correlated`
   - `discrete`
   - `relaxed`
4. `treePL` CI pass after the plain winner is already visible.
5. `chronos` CI pass after that, for:
   - `clock`
   - `correlated`
   - `discrete`
   `relaxed` stays plain.

## MAIN_OUTPUT_TREES policy

- Put trees into `MAIN_OUTPUT_TREES` as soon as each stage finishes.
- Do not wait for the full pipeline before exposing trees.
- Early state should show:
  - `RelTime MEGA` with CIs
  - `treePL` plain
  - `chronos` plain winners
- Later CI stages should replace the plain files in `MAIN_OUTPUT_TREES`:
  - replace plain `treePL` with `treePL_*_with_CI.tre`
  - replace plain `chronos_clock`, `chronos_correlated`, and `chronos_discrete` with `_with_CI.tre`
  - keep `chronos_relaxed` plain

## Runner

Use:

`/Users/ricardobetancur/Desktop/Proxy_Misplaced/PhyloChronoRank/scripts/run_dating_grid_staged.py`

instead of calling `run_dating_grid.R` directly when the goal is a PCR-ready run folder with incrementally updated `MAIN_OUTPUT_TREES`.

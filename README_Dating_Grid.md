# Dating Grid Guide

This page documents `scripts/run_dating_grid.R`, the repo-local helper for generating candidate chronograms before running PCR.

The script is designed for one specific use case:

- start from an unconstrained phylogram
- define one shared calibration source
- run the same effective calibration set through `chronos`, `treePL`, and `RelTime`
- collect all successful dated trees into one `candidates.csv` that can be passed directly into `scripts/run_pcr.R`

## What It Runs

`scripts/run_dating_grid.R` can run:

- `chronos` across a lambda grid and the four supported clock models: `clock`, `correlated`, `relaxed`, `discrete`
- `treePL` across a smoothing grid
- `RelTime` with the same merged node bounds used for the other methods

The main point is calibration consistency. The script first resolves one shared calibration table, maps pairwise calibrations onto MRCA nodes on the target phylogram, merges duplicate-node rows by interval intersection, drops empty intersections, and then passes that same resolved node-bound set to all three methods.

## Inputs

You must provide:

- one rooted phylogram with positive branch lengths
- one output directory
- exactly one calibration source:
  - a pairwise calibration CSV
  - or a reference backbone time tree for congruification

Optional:

- an exact `root_age`
- a reduced method list
- custom lambda and smoothing grids
- a method tag filter if your calibration CSV contains a `candidate` column

## Calibration Modes

### Option 1: Pairwise Calibration CSV

Required columns:

- `taxonA`
- `taxonB`
- `age_min`

Optional columns:

- `age_max`
- `candidate`

Notes:

- if `age_max` is missing, the script treats that row as minimum-only with an open upper bound
- if a `candidate` column exists and contains multiple method-specific tags, pass `--calibration-tag=...` so the script knows which subset to keep before the shared calibration merge
- if you pass `--root-age=...`, that exact root age is appended to the shared calibration set for all methods

### Option 2: Reference Backbone Time Tree

Instead of a calibration CSV, you can provide `--reference-time-tree=...`.

In that mode the script:

- congruifies the reference time tree onto the target phylogram
- converts the congruified pairwise calibrations into one shared pairwise table
- maps them onto MRCA nodes on the target tree
- merges duplicate-node bounds exactly once before method-specific runs

## Basic Usage

```bash
Rscript scripts/run_dating_grid.R \
  --phylogram=PATH/TO/phylogram.tre \
  --calibrations-csv=PATH/TO/calibrations.csv \
  --outdir=PATH/TO/dating_out
```

Or with congruification:

```bash
Rscript scripts/run_dating_grid.R \
  --phylogram=PATH/TO/phylogram.tre \
  --reference-time-tree=PATH/TO/reference_time_tree.tre \
  --outdir=PATH/TO/dating_out
```

## Common Options

```text
--methods=chronos,treepl,reltime
--chronos-lambdas=0.01,0.1,1,10,100
--chronos-models=clock,correlated,relaxed,discrete
--chronos-discrete-k=5
--chronos-retries=2
--treepl-smoothing=0.01,0.1,1,10,100
--treepl-bin=/absolute/path/to/treePL
--treepl-numsites=1000
--treepl-threads=1
--reltime-sites=1000
--root-age=123.4
--out-prefix=my_dataset
```

## Terap Example

Using the bundled Terap phylogram and calibration CSV:

```bash
Rscript scripts/run_dating_grid.R \
  --phylogram=examples/terapontoid/Terapontoid_ML_MAIN_phylogram_used.tree \
  --calibrations-csv=examples/terapontoid/Terapontoid_ML_MAIN_calibrations_used.csv \
  --outdir=examples/terapontoid/dating_grid_out \
  --chronos-lambdas=0.01,0.1,1,10,100 \
  --treepl-smoothing=0.01,0.1,1,10,100 \
  --chronos-discrete-k=5
```

If you want a fast smoke test first:

```bash
Rscript scripts/run_dating_grid.R \
  --phylogram=examples/terapontoid/Terapontoid_ML_MAIN_phylogram_used.tree \
  --calibrations-csv=examples/terapontoid/Terapontoid_ML_MAIN_calibrations_used.csv \
  --outdir=/tmp/pcr_dating_grid_smoke \
  --chronos-lambdas=0.1 \
  --treepl-smoothing=0.1 \
  --chronos-discrete-k=5
```

## Outputs

The script writes:

- `candidates.csv`
  - one row per successful dated tree
  - ready for `scripts/run_pcr.R`
- `<prefix>_all_runs_summary.csv`
  - combined run summary across all requested methods
- `<prefix>_calibration_pairs_used.csv`
  - the pairwise calibration table actually used
- `<prefix>_calibration_pairs_mapped.csv`
  - those pairwise rows after MRCA mapping onto the target tree
- `<prefix>_calibration_bounds_used.csv`
  - the final merged node-bound table shared by all methods
- `<prefix>_calibration_bounds_dropped.csv`
  - duplicate-node conflicts that collapsed to empty intervals and were dropped
- `<prefix>_run_metadata.txt`
  - a compact provenance summary

Method-specific outputs:

- `chronos/trees/*.tre`
- `chronos/<prefix>_chronos_runs.csv`
- `treepl/configs/*.cfg`
- `treepl/logs/*.log`
- `treepl/trees/*.tre`
- `treepl/<prefix>_treepl_runs.csv`
- `reltime/<prefix>_RelTime.tre`
- `reltime/<prefix>_RelTime_ci.csv`
- `reltime/<prefix>_RelTime_bounds_used.csv`
- `reltime/<prefix>_RelTime_run.csv`

## Feeding The Results Into PCR

Once the dating grid finishes, run PCR on the generated candidates:

```bash
Rscript scripts/run_pcr.R \
  --ref-tree=PATH/TO/phylogram.tre \
  --candidates-csv=PATH/TO/dating_out/candidates.csv \
  --calibrations-csv=PATH/TO/dating_out/<prefix>_calibration_pairs_used.csv \
  --outdir=PATH/TO/pcr_out
```

If you are working from primary fossil calibrations, that shared calibration CSV can be used directly in PCR’s gap layer. If your dates were generated from congruified or secondary ages, treat that gap layer as calibration slack rather than as independent fossil-fit evidence.

## Requirements

R packages:

- `ape`
- `geiger` for congruification mode
- `phytools` if the reference time tree needs ultrametric extension before congruification
- `quadprog` for full-bound `RelTime` projection

External binary:

- `treePL` is required only if `treePL` is included in `--methods`

The script looks for `treePL` in this order:

- `--treepl-bin=...`
- `TREEPL_BIN` environment variable
- `treePL` on `PATH`
- `../treePL` relative to the repo root

## Important Behavior

- all three methods use the same resolved calibration set after MRCA mapping and duplicate-node merging
- duplicate-node conflicts are merged by interval intersection
- empty intersections are dropped for everyone, not just for one method
- `RelTime` is run with the repo-local helper in `scripts/reltime_helpers.R`
- `RelTime` CI files are produced, but those widths can become numerically unstable when hard bounds create very short internal branches
- `chronos` failures and `treePL` failures are retained in the run summary even when they do not appear in `candidates.csv`

## Help

The script also prints its own usage summary:

```bash
Rscript scripts/run_dating_grid.R --help
```

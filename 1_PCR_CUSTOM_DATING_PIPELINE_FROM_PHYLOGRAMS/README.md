# PCR Custom Dating Pipeline From Phylograms

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

For `treePL`, the default repo path is the recommended two-step run:

- a `prime` pass first
- then the real optimized run using the optimizer hints reported by the prime pass
- with `thorough` enabled by default
- and with no explicit `opt` override unless you pass `--treepl-opt=...`

The main point is calibration consistency. The script first resolves one shared calibration table, maps pairwise calibrations onto MRCA nodes on the target phylogram, merges duplicate-node rows by interval intersection, drops empty intersections, and then passes that same resolved node-bound set to all three methods.

## Method Provenance

The three method paths in this repo are tied directly to the core method papers:

- `chronos`: [Paradis 2013, Molecular dating of phylogenies by likelihood methods: A comparison of models and a new information criterion](https://www.sciencedirect.com/science/article/abs/pii/S1055790313000651)
- `treePL`: [Smith and O'Meara 2012, treePL: divergence time estimation using penalized likelihood for large phylogenies](https://academic.oup.com/bioinformatics/article/28/20/2689/203074)
- `RelTime`: [Tamura et al. 2012, Estimating divergence times in large molecular phylogenies](https://pubmed.ncbi.nlm.nih.gov/23129628/) and [Tamura, Tao, and Kumar 2018, Theoretical Foundation of the RelTime Method for Estimating Divergence Times from Variable Evolutionary Rates](https://pubmed.ncbi.nlm.nih.gov/29893954/)
- `RelTime` confidence intervals: [Tao, Tamura, Mello, and Kumar 2020, Reliable confidence intervals for RelTime estimates of evolutionary divergence times](https://academic.oup.com/mbe/article/37/1/280/5602325)

Operationally, all three methods are exposed here through one R-driven workflow:

- `chronos` is run directly through `ape::chronos`
- `treePL` is driven from R by writing the control files and calling the external `treePL` binary from the script
- `RelTime` is implemented in repo-local R code derived from the relative-rate framework papers above, so this workflow does not require `MEGA`

In other words, the shipped `treePL` path here is not a minimal one-shot wrapper. It runs the recommended `prime + thorough` workflow from R, then validates the dated output tree before adding it to `candidates.csv`.

That means this repo is not just storing finished chronograms. It can also generate a comparable multi-method candidate set in one place, using one shared calibration resolution step before PCR scoring.

## Three Layers

This workflow is easiest to read if you keep three layers separate:

- dating / fit layer
  - `run_dating_grid.R` generates dated trees and method-specific run summaries
- tuning layer
  - lambda, smoothing, and `nb.rate.cat` grids are there to expose sensitivity, not to guarantee that one single setting must always be declared the winner by fit alone
- post-fit layer
  - PCR, run afterward with `scripts/run_pcr.R`, scores the finished candidate chronograms on a common biological diagnostic set

That separation matters. A tree can look best by a fit criterion inside one method family and still lose the broader PCR comparison after all methods are placed on the same post-fit scale.

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

## Large Trees And Subset Tuning

A large-tree subset strategy can be useful: tune on a smaller calibration-preserving subset, then rerun the selected settings on the full tree.

`run_dating_grid.R` does not automate that subset workflow. If you need it, the practical pattern is:

- build a reduced tree that preserves all calibration taxa
- keep deep or tempo-extreme parts of the tree rather than using a purely random subset
- tune `chronos` lambda values and `treePL` smoothing values on that subset
- rerun the selected settings on the full phylogram with the same shared calibration file

For modest tree sizes, the direct full-tree grid is usually simpler. For very large trees, subset tuning can be the difference between a usable screening pass and an impractically slow one.

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
--treepl-thorough=TRUE
--treepl-prime=TRUE
--reltime-sites=1000
--root-age=123.4
--out-prefix=my_dataset
```

## Terapontoidei Example

Using the bundled Terapontoidei phylogram and calibration CSV:

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

## What To Open First

If you are reviewing a new run, the most useful order is:

1. `<prefix>_run_metadata.txt`
2. `<prefix>_all_runs_summary.csv`
3. `chronos/<prefix>_chronos_runs.csv`
4. `treepl/<prefix>_treepl_runs.csv`
5. `candidates.csv`

That tells you:

- which calibrations survived the shared merge step
- which runs actually succeeded
- which candidate trees were promoted into the PCR-ready set

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

## How To Choose Among Grid Results

The script writes the full grid. It does not enforce one global selector for you.

The practical review pattern is:

- compare fit within a method family first
  - for `chronos`, inspect `chronos/<prefix>_chronos_runs.csv`
  - for `treePL`, inspect `treepl/<prefix>_treepl_runs.csv`
- if neighboring settings are near-tied or biologically different, keep more than one candidate rather than pretending the fit surface is sharper than it is
- then run PCR on the retained candidate set and compare methods on the same post-fit scale

In other words:

- use the dating grid to expose sensitivity
- use PCR to arbitrate among the finished candidate chronograms
- if fit and post-fit disagree, report both instead of collapsing them into one claim

## How To Read Disagreements

Disagreement across layers is normal and often biologically informative rather than a sign that one step failed.

Common patterns are:

- one setting looks best by fit, but a neighboring setting looks better by PCR
  - that usually means the likelihood-style fit surface is shallow relative to the biological shape differences PCR is detecting
- one method family dominates fit inside its own grid, but loses to another method family in PCR
  - that means the method-specific optimization target and the post-fit biological diagnostics are rewarding different properties
- a candidate satisfies calibrations tightly, but looks poor by pulse preservation or rate regularity
  - that is a real tradeoff worth stating explicitly

The practical rule is simple: do not force those layers into one synthetic winner if they are telling different stories. Report the fit preference, report the PCR preference, and explain the biological consequence of the difference.

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

So if you want the repo fallback to work without passing `--treepl-bin` or setting `TREEPL_BIN`, place the executable one directory above the repo root and name the file exactly `treePL`. In placeholder form, that fallback looks like `/PATH/TO/PARENT_FOLDER/treePL`.

## Important Behavior

- all three methods use the same resolved calibration set after MRCA mapping and duplicate-node merging
- duplicate-node conflicts are merged by interval intersection
- empty intersections are dropped for everyone, not just for one method
- `RelTime` is run with the repo-local helper in `scripts/reltime_helpers.R`
- `treePL` defaults to `thorough = TRUE` and `prime = TRUE`, with a real post-prime optimization pass rather than stopping after the priming step
- `RelTime` CI files are produced, but those widths can become numerically unstable when hard bounds create very short internal branches
- `chronos` failures and `treePL` failures are retained in the run summary even when they do not appear in `candidates.csv`

## Common Issues

- no `treePL` binary found
  - pass `--treepl-bin=/PATH/TO/treePL`, set `TREEPL_BIN`, or use the documented fallback location
- no successful dated trees written
  - inspect `<prefix>_all_runs_summary.csv` and the method-specific log files first
- dropped calibration bounds
  - check `<prefix>_calibration_bounds_dropped.csv`; these are duplicate-node conflicts that collapsed to empty intersections
- `chronos` starting-date failures
  - increase `--chronos-retries`, simplify the grid, or use subset tuning for very large trees
- suspicious congruification results
  - inspect `<prefix>_calibration_pairs_mapped.csv` and verify taxon overlap and reference-tree quality before blaming the dating method itself

## Help

The script also prints its own usage summary:

```bash
Rscript scripts/run_dating_grid.R --help
```

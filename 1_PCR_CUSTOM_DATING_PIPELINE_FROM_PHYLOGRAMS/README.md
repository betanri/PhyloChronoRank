# PCR Custom Dating Pipeline From Phylograms

This page documents `scripts/run_dating_grid.R`, the repo-local helper for generating candidate chronograms before running PCR.

> **For PCR-targeted runs, the default runner is `scripts/run_dating_grid_staged.py`, not `run_dating_grid.R` directly.** The staged runner drives the same `run_dating_grid.R` engine in five ordered stages and incrementally populates a `MAIN_OUTPUT_TREES/` folder so PCR-ready trees appear as each stage finishes rather than only at the end of the pipeline. See [`PCR_DATING_STAGED_WORKFLOW.md`](../PCR_DATING_STAGED_WORKFLOW.md) for the canonical spec. The stage order is:
>
> 1. `RelTime MEGA` (with both analytical and bootstrap CIs)
> 2. `treePL` plain winner
> 3. `chronos` plain winners for all four clock models (`clock`, `correlated`, `discrete`, `relaxed`), seeded by the treePL winner
> 4. `treePL` CI pass — replaces plain treePL in `MAIN_OUTPUT_TREES`
> 5. `chronos` CI pass for `clock`, `correlated`, `discrete` — replaces their plain counterparts in `MAIN_OUTPUT_TREES`; `chronos_relaxed` stays plain
>
> Early-stage `MAIN_OUTPUT_TREES` already contains a usable PCR candidate set (RelTime with CIs, plain treePL, plain chronos winners) before the expensive CI passes finish. The sections below describing `run_dating_grid.R` still apply — the staged runner just orchestrates calls to that same R script.
>
> ```bash
> python3 scripts/run_dating_grid_staged.py \
>   --phylogram=PATH/TO/phylogram.tre \
>   --calibrations-csv=PATH/TO/calibrations.csv \
>   --outdir=PATH/TO/dating_out \
>   --out-prefix=my_dataset
> ```

## What It Runs

`scripts/run_dating_grid.R` can run:

- `chronos` across a lambda grid and the four supported clock models: `clock`, `correlated`, `relaxed`, `discrete`
- `treePL` across a smoothing grid
- `RelTime` via MEGA-CC (`megacc` binary required; freely available from [megasoftware.net](https://www.megasoftware.net))
- optional uncertainty summaries: all three methods (`chronos`, `treePL`, `RelTime`) compute bootstrap CIs via Poisson branch-length resampling ([Paradis et al. 2023](https://doi.org/10.1016/j.ympev.2022.107652); `chronos` helper adapted from [josephwb/chronos](https://github.com/josephwb/chronos)). `RelTime` additionally reports analytical delta-method CIs parsed from MEGA's native output ([Tao et al. 2020](https://academic.oup.com/mbe/article/37/1/280/5602325))

The main point is calibration consistency: all three methods receive exactly the same age constraints. The script reads one shared calibration table, finds the tree node that corresponds to each calibration (MRCAs of the two taxa listed), combines any rows that map to the same node by keeping only the overlap of their age ranges, drops any calibrations that cannot be placed (taxa not found in the tree) or that produce contradictory constraints on the same node (e.g., one row says ≥10 Ma while another says ≤5 Ma), and then passes that single cleaned set of node-age bounds identically to chronos, treePL, and RelTime. The script reports any dropped calibrations at runtime and writes them to `*_calibration_bounds_dropped.csv` (see [Outputs](#outputs)). This ensures that any differences among the resulting chronograms come from the dating methods themselves, not from differences in calibration input.

## What You Get

In the standard full-method comparison, the pipeline returns six candidate chronograms:

- `chronos_clock`
- `chronos_discrete`
- `chronos_correlated`
- `chronos_relaxed`
- `treePL`
- `RelTime_MEGA` — produced by MEGA-CC's RelTime implementation

The `chronos` side is kept as four separate trees, one per clock model, rather than collapsed to a single selected `chronos` result. Variation among `chronos` clock models is often large, and in practice can exceed the difference between some `chronos` models and `treePL` or `RelTime`. Keeping all four therefore avoids hiding biologically meaningful variation before PCR evaluates it.

**Example output (terapontoid dataset, 105 tips, 19 calibrations from congruification):**

<img src="https://raw.githubusercontent.com/betanri/PhyloChronoRank/main/figures/terapontoid_pipeline_output.png" alt="Terapontoid pipeline output" width="100%">

## Method Provenance

The three method paths in this repo are tied directly to the core method papers:

- `chronos`: [Paradis 2013, Molecular dating of phylogenies by likelihood methods: A comparison of models and a new information criterion](https://www.sciencedirect.com/science/article/abs/pii/S1055790313000651)
- `chronos` bootstrap confidence intervals: [Paradis et al. 2023, Confidence intervals in molecular dating by maximum likelihood](https://doi.org/10.1016/j.ympev.2022.107652)
- `treePL`: [Smith and O'Meara 2012, treePL: divergence time estimation using penalized likelihood for large phylogenies](https://academic.oup.com/bioinformatics/article/28/20/2689/203074)
- `RelTime`: [Tamura et al. 2012, Estimating divergence times in large molecular phylogenies](https://pubmed.ncbi.nlm.nih.gov/23129628/) and [Tamura, Tao, and Kumar 2018, Theoretical Foundation of the RelTime Method for Estimating Divergence Times from Variable Evolutionary Rates](https://pubmed.ncbi.nlm.nih.gov/29893954/)
- `RelTime` confidence intervals: [Tao, Tamura, Mello, and Kumar 2020, Reliable confidence intervals for RelTime estimates of evolutionary divergence times](https://academic.oup.com/mbe/article/37/1/280/5602325)

Operationally, all methods are exposed here through one R-driven workflow:

- `chronos` is run directly through `ape::chronos`
- `treePL` is driven from R by writing the control files and calling the external `treePL` binary from the script
- `RelTime` requires the `megacc` binary (MEGA-CC, freely available from [megasoftware.net](https://www.megasoftware.net) for all platforms). The pipeline auto-generates all MEGA input files (.mao settings, calibration file, outgroup file), grafts a temporary mock outgroup at the root (required because MEGA's RelTime cannot calibrate the root node directly), runs `megacc`, strips the outgroup from the output, and cleans tip labels
- all three methods compute bootstrap CIs via Poisson branch-length resampling ([Paradis et al. 2023](https://doi.org/10.1016/j.ympev.2022.107652); `chronos` helper adapted from [josephwb/chronos](https://github.com/josephwb/chronos)). `RelTime` additionally reports analytical delta-method CIs parsed from MEGA's native output ([Tao et al. 2020](https://academic.oup.com/mbe/article/37/1/280/5602325))
- all CI-annotated trees are written as FigTree-ready NEXUS files with `height_95%_HPD` node-bar annotations

This repo can generate a comparable multi-method candidate set in one place, using one shared calibration resolution step before PCR scoring.

## Four Layers

This workflow is easiest to read if you keep four layers separate:

- dating / fit layer
  - `run_dating_grid.R` generates dated trees and method-specific run summaries for each candidate (`chronos`, `treePL`, `RelTime_MEGA`)
- model-selection layer (`chronos` only)
  - `chronos` is the only method in this comparison that supports multiple clock models (`clock`, `correlated`, `relaxed`, `discrete`). The best-fitting model is selected by `PHIIC` ([Paradis 2013](https://doi.org/10.1016/j.ympev.2013.02.008)) before candidates enter post-fit scoring. `treePL` and `RelTime_MEGA` each produce a single candidate — no model selection applies
- tuning layer
  - `chronos`: the `--chronos-lambdas` grid (default `0.01,0.1,1,10,100`) controls the smoothing penalty across all four clock models. For the `discrete` model only, `--chronos-discrete-k` (default `2,3,5`) additionally searches the number of rate categories. The best lambda (and k, for discrete) per model is selected by `PHIIC`
  - `treePL`: the `--treepl-smoothing` grid (default `0.01,0.1,1,10,100`) controls the smoothing parameter. The best smoothing value is selected by treePL's internal cross-validation objective
  - `RelTime_MEGA`: no tuning parameters — RelTime produces a single tree from the input phylogram and calibrations
- post-fit layer (see [2_PCR_POSTFIT_METRICS](../2_PCR_POSTFIT_METRICS/README.md))
  - PCR, run afterward with `scripts/run_pcr.R`, scores the finished candidate chronograms on a common biological diagnostic set

That separation matters. A tree can look best by a fit criterion inside one method family and still lose the broader PCR comparison after all methods are placed on the same post-fit scale.

## Inputs

You must provide:

- one rooted phylogram with positive branch lengths measured in substitutions per site (not in coalescent time units)
- one output directory
- exactly one calibration source:
  - a pairwise calibration CSV
  - or a reference backbone time tree for congruification

Optional:

- an exact `root_age`
- a reduced method list
- custom lambda and smoothing grids
- a method tag filter if your calibration CSV contains a `candidate` tag column

## Calibration Modes

### Option 1: Pairwise Calibration CSV

Required columns:

- `taxonA`
- `taxonB`
- `age_min`

Optional columns:

- `age_max`
- `candidate` (an optional tag column used only to filter rows when one CSV contains multiple calibration subsets or method-specific rows)

Notes:

- if `age_max` is missing, the script treats that row as minimum-only with an open upper bound
- if a `candidate` column exists, its values are just row labels for filtering; they do not name candidate trees
- if the same CSV mixes multiple tags, pass `--calibration-tag=...` so the script keeps the desired tagged rows before the shared calibration merge
- blank `candidate` cells can be used for calibration rows that should apply to all methods
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
Tree size also affects method behavior itself: [Paradis et al. 2023](https://doi.org/10.1016/j.ympev.2022.107652) showed that coverage and CI width vary with `n`, and our companion benchmark (`P1` and `P2`) tests this directly by comparing `MAE` across `n = 50`, `150`, and `300` tips for `treePL`, `chronos`, and `RelTime` under the same calibration design.

## Basic Usage

For PCR-targeted runs, prefer the staged wrapper so `MAIN_OUTPUT_TREES/` is populated incrementally:

```bash
python3 scripts/run_dating_grid_staged.py \
  --phylogram=PATH/TO/phylogram.tre \
  --calibrations-csv=PATH/TO/calibrations.csv \
  --outdir=PATH/TO/dating_out \
  --out-prefix=my_dataset
```

Or, to drive the R engine directly (single-pass, no incremental `MAIN_OUTPUT_TREES`):

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
--methods=chronos,treepl,reltime_mega
--megacc-bin=/absolute/path/to/megacc    # required for RelTime
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
--ci-sites=1000
--treepl-bootstrap-reps=100
--treepl-bootstrap-jobs=4
--megacc-bootstrap-reps=100
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
  - the representative PCR-ready candidate set
  - in the standard full-method comparison, this is the reduced 6-tree table: one retained `chronos` tree per clock model, one retained `treePL` tree, and `RelTime_MEGA`
  - ready for `scripts/run_pcr.R`
- `full_grid/candidates.csv`
  - one row per successful dated tree across the full grid
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
  - calibrations that were dropped during the shared merge step; the script also prints a summary of all drops at runtime (missing taxa, unmappable MRCAs, and conflicting-node intervals) so you can catch problems before the dating runs start
- `<prefix>_run_metadata.txt`
  - a compact provenance summary

<details>
<summary><strong>Method-specific outputs</strong></summary>

- `chronos/trees/*.tre`
- `chronos/ci/*.csv`
- `chronos/<prefix>_chronos_runs.csv`
- `treepl/configs/*.cfg`
- `treepl/logs/*.log`
- `treepl/trees/*.tre`
- `treepl/ci/*.csv`
- `treepl/<prefix>_treepl_runs.csv`
- `reltime_mega/<prefix>_RelTime_MEGA.tre`
- `reltime_mega/<prefix>_RelTime_MEGA_with_tao_CI.tre` — FigTree NEXUS with Tao analytical CIs
- `reltime_mega/<prefix>_RelTime_MEGA_tao_ci.csv`
- `reltime_mega/<prefix>_RelTime_MEGA_run.csv`
- `reltime_mega/<prefix>_mega_reltime.log`
- `uncertainty_summary_long.csv`

</details>

## What To Open First

If you are reviewing a new run, the most useful order is:

1. `<prefix>_run_metadata.txt`
2. `<prefix>_all_runs_summary.csv`
3. `chronos/<prefix>_chronos_runs.csv`
4. `treepl/<prefix>_treepl_runs.csv`
5. `uncertainty_summary_long.csv`
6. `candidates.csv`
7. `full_grid/candidates.csv`

That tells you:

- which calibrations survived the shared merge step
- which runs actually succeeded
- which CI summaries were written for the successful dated trees
- which candidate trees were promoted into the reduced PCR-ready set
- which additional successful trees remain available in the full grid

## Feeding The Results Into PCR

Once the dating grid finishes, run PCR on the reduced representative candidates:

```bash
Rscript scripts/run_pcr.R \
  --ref-tree=PATH/TO/phylogram.tre \
  --candidates-csv=PATH/TO/dating_out/candidates.csv \
  --calibrations-csv=PATH/TO/dating_out/<prefix>_calibration_pairs_used.csv \
  --uncertainty-csv=PATH/TO/dating_out/uncertainty_summary_long.csv \
  --outdir=PATH/TO/pcr_out
```

If you are working from primary fossil calibrations, that shared calibration CSV can be used directly in PCR’s gap layer. If your dates were generated from congruified or secondary ages, treat that gap layer as calibration slack rather than as independent fossil-fit evidence.

## How To Choose Among Grid Results

The script writes both the full grid and a reduced representative candidate set.

The practical review pattern is:

- inspect the reduced `candidates.csv` first if you want the standard 6-tree PCR comparison
- inspect `full_grid/candidates.csv` and the run summaries if you want to review the full sensitivity surface
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

## RelTime via MEGA-CC

The pipeline uses MEGA-CC's RelTime implementation ([Tamura et al. 2012](https://pubmed.ncbi.nlm.nih.gov/23129628/); [Tamura, Tao, and Kumar 2018](https://pubmed.ncbi.nlm.nih.gov/29893954/)) as its RelTime engine. The `megacc` binary is required and is freely available from [megasoftware.net](https://www.megasoftware.net) for all platforms.

The pipeline auto-generates all MEGA input files (.mao settings, calibration file, outgroup file) from the same shared calibration set used by `chronos` and `treePL`. Because MEGA's RelTime cannot constrain the root node directly — it requires a designated outgroup lineage external to all calibrated nodes — the pipeline excludes the root calibration and grafts a temporary single-lineage mock outgroup at the base of the tree before invoking `megacc`. Without this step, the root age is left unconstrained and deep-node confidence intervals can span the entire tree depth. After the run, the mock outgroup is stripped from the output tree, and tip labels are cleaned (MEGA adds quotes and replaces underscores with spaces).

Confidence intervals are computed two ways: (1) analytical delta-method CIs parsed from MEGA's native NEXUS output ([Tao et al. 2020](https://academic.oup.com/mbe/article/37/1/280/5602325)), and (2) Poisson branch-length bootstrap CIs ([Paradis et al. 2023](https://doi.org/10.1016/j.ympev.2022.107652)), where each replicate perturbs the phylogram's branch lengths, runs `megacc` on the perturbed tree, and collects node ages for percentile CIs. Both CI types produce FigTree-ready NEXUS trees with `height_95%_HPD` node-bar annotations.

### Usage

```bash
# megacc must be on PATH or specified with --megacc-bin:
Rscript scripts/run_dating_grid.R \
  --phylogram=tree.tre \
  --calibrations-csv=cals.csv \
  --outdir=out \
  --methods=chronos,treepl,reltime_mega

# If megacc is not on PATH:
Rscript scripts/run_dating_grid.R \
  --phylogram=tree.tre \
  --calibrations-csv=cals.csv \
  --outdir=out \
  --methods=chronos,treepl,reltime_mega \
  --megacc-bin=/path/to/megacc
```

## Requirements

R packages:

- `ape`
- `geiger` for congruification mode
- `phytools` if the reference time tree needs ultrametric extension before congruification
- `quadprog` for full-bound `RelTime` projection

External binaries:

- `treePL` is required only if `treepl` is included in `--methods`
- `megacc` (MEGA-CC) is required only if `reltime_mega` is included in `--methods`

The script looks for `treePL` in this order:

- `--treepl-bin=...`
- `TREEPL_BIN` environment variable
- `treePL` on `PATH`
- `../treePL` relative to the repo root

The script looks for `megacc` in this order:

- `--megacc-bin=...`
- `MEGACC_BIN` environment variable
- `megacc` on `PATH`

MEGA-CC can be downloaded from [megasoftware.net](https://www.megasoftware.net/). The command-line version (`megacc`) is what this pipeline uses — the GUI is not needed.

So if you want the repo fallback to work without passing `--treepl-bin` or setting `TREEPL_BIN`, place the executable one directory above the repo root and name the file exactly `treePL`. In placeholder form, that fallback looks like `/PATH/TO/PARENT_FOLDER/treePL`.

## Important Behavior

- all three methods use the same resolved calibration set after MRCA mapping and duplicate-node merging
- duplicate-node conflicts are merged by interval intersection
- empty intersections are dropped for everyone, not just for one method
- `RelTime` is run via MEGA-CC (`megacc`); the pipeline auto-generates all input files, grafts a mock outgroup, runs `megacc`, strips the outgroup, and cleans tip labels
- all three methods (`chronos`, `treePL`, `RelTime`) compute bootstrap CIs via Poisson branch-length resampling ([Paradis et al. 2023](https://doi.org/10.1016/j.ympev.2022.107652)); `RelTime` additionally reports Tao analytical CIs from MEGA's native output
- `uncertainty_summary_long.csv` combines uncertainty summaries across all methods on the same comparison scale
- the Tao-style analytical `RelTime` CI is reported separately only; on hard-bounded empirical trees it can live in a completely different numerical regime from the bootstrap widths because the analytical variance term can explode after bound projection compresses internal durations
- `treePL` defaults to `thorough = TRUE` and `prime = TRUE`, with a real post-prime optimization pass rather than stopping after the priming step
- `uncertainty_summary_long.csv` is written in PCR-ready format for `--uncertainty-csv`
- in the bracketed benchmarks, `RelTime` is run only under the constraints it can actually consume, which means internal max bounds are not used on the `RelTime` side; those runs therefore compare method-specific usable calibration information rather than a literal all-methods-share-every-bound setup
- `chronos` failures and `treePL` failures are retained in the run summary even when they do not appear in `candidates.csv`

## Common Issues

- no `treePL` binary found
  - pass `--treepl-bin=/PATH/TO/treePL`, set `TREEPL_BIN`, or use the documented fallback location
- no `megacc` binary found (only relevant if `reltime_mega` is in `--methods`)
  - pass `--megacc-bin=/PATH/TO/megacc`, set `MEGACC_BIN`, or place `megacc` on your `PATH`
  - download MEGA-CC from [megasoftware.net](https://www.megasoftware.net/)
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

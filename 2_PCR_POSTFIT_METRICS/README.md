# PCR Postfit Metrics

`PCR Postfit Metrics` is a post-fit evaluation framework for phylogeneticists who already have a set of competing chronograms and need to decide which one is the most biologically defensible.

In simple terms, the core idea is that divergence-time estimation is hard: the resulting chronogram can shift substantially with clock-model choice, tree priors, calibration priors, and other analytical decisions ([Lepage et al. 2007](https://doi.org/10.1093/molbev/msm193); [Warnock et al. 2015](https://doi.org/10.1098/rspb.2014.1013); [dos Reis et al. 2016](https://doi.org/10.1038/nrg.2015.8)). A practical response is to compare a defensible set of alternative chronograms after they have been estimated, rather than betting everything on a single, computationally intensive "gold-standard" analysis and then treating that one tree as settled.

It is method-agnostic. The candidates can come from `BEAST`, `MCMCTree`, `MrBayes`, `chronos`, `treePL`, `RelTime`, or any other dating workflow.

PCR starts from finished chronograms. If you need to generate a set of chronograms from three alternative methods (`chronos`, `treePL`, `RelTime`) using an unconstrained phylogram and a calibration set, you can do that here: [PCR Custom Dating Pipeline From Phylograms](../1_PCR_CUSTOM_DATING_PIPELINE_FROM_PHYLOGRAMS/README.md).

## What it evaluates

`PhyloChronoRank (PCR)` uses three core metric families. These are implementation-level diagnostics rather than named published indices; the citations below support the underlying ideas each family is trying to capture.

- `pulse preservation`: asks whether a dated tree keeps the same branching rhythm seen in the source phylogram. In practice, this means preserving clustered speciation bursts and quiet intervals instead of smearing them into evenly spaced splits. In this workflow, the pulse family is reported three ways: `burst loss` is the standalone burst-flattening submetric, `pulse preservation (burst)` is the burst-priority composite selector, and `pulse preservation (overall)` is the balanced composite selector. This follows the literature on extracting diversification tempo from phylogenies and on distinguishing burst-like from unusually regular branching patterns ([Nee et al. 1992](https://doi.org/10.1073/pnas.89.17.8322); [Pybus and Harvey 2000](https://doi.org/10.1098/rspb.2000.1278); [Ford et al. 2009](https://doi.org/10.1093/sysbio/syp018)).

- `gap burden`: asks how much extra unseen lineage history the dated tree implies relative to the calibration evidence. This is the same general idea as ghost-lineage and stratigraphic-congruence measures ([Huelsenbeck 1994](https://doi.org/10.1017/S009483730001294X); [Wills 1999](https://doi.org/10.1080/106351599260148); [O'Connor and Wills 2016](https://doi.org/10.1093/sysbio/syw039)). Lower is usually better, but it should be interpreted carefully: fossils usually provide minimum ages, not true lineage origins, so a tree that minimizes this too aggressively can simply be too young overall ([Parham et al. 2012](https://doi.org/10.1093/sysbio/syr107)). In PCR, this should be treated as a core metric only when the calibrations are primary evidence. With secondary or congruified ages, the same calculation becomes calibration slack against inherited ages rather than an independent biological diagnostic.

- `rate irregularity`: for each branch, divides the phylogram branch length (substitutions) by the chronogram branch duration (time) to get an implied evolutionary rate. The score rises when those implied rates are too dispersed, jump sharply from parent to child branch, produce too many outlier branches, or lose the positive autocorrelation expected among closely related lineages. This follows the penalized-likelihood and relaxed-clock literature on among-lineage rate variation and autocorrelation ([Sanderson 2002](https://doi.org/10.1093/oxfordjournals.molbev.a003974); [Drummond et al. 2006](https://doi.org/10.1371/journal.pbio.0040088); [Lepage et al. 2007](https://doi.org/10.1093/molbev/msm193); [Ho 2009](https://doi.org/10.1098/rsbl.2008.0729); [Tao et al. 2019](https://doi.org/10.1093/molbev/msz014)).

- `uncertainty width` (optional precision layer): measures how wide the confidence or credible intervals are around estimated node ages. Narrower intervals indicate greater precision, though precision and accuracy are distinct properties. Frequentist confidence intervals and Bayesian credibility intervals arise from different statistical frameworks and are not directly commensurable ([Paradis et al. 2023](https://doi.org/10.1016/j.ympev.2022.107652); [Drummond et al. 2006](https://doi.org/10.1371/journal.pbio.0040088); [Bromham 2019](https://doi.org/10.1016/j.tree.2019.01.017); [Tao et al. 2020a](https://academic.oup.com/mbe/article/37/1/280/5602325); [Tao et al. 2020b](https://academic.oup.com/mbe/article/37/6/1819/5771370)), so they are usually reported side by side rather than collapsed into one score ([Costa et al. 2022](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-09030-5); [Beavan et al. 2020](https://academic.oup.com/gbe/article/12/7/1087/5842139)). PhyloChronoRank scores this layer only when intervals derive from comparable methodology across candidates. In the bundled examples, Syngnatharia uses extracted HPD widths, while Terapontoidei and the bundled vertebrate release use matched bootstrap summaries for `chronos`, `treePL`, and `RelTime`. The analytical delta-method `RelTime` CI of [Tao et al. 2020a](https://academic.oup.com/mbe/article/37/1/280/5602325) is kept as a supplemental diagnostic because it can diverge sharply from bootstrap widths after bound projection.

<details>
<summary><strong>Compact formulas used in the current implementation</strong></summary>

`burst loss`

```text
burst_loss_clade = max(0, (burst_ref - burst_est) / (burst_ref + 1e-12))
mean_burst_loss = weighted mean across matched clades
weight_clade = log(1 + n_tips) * sqrt(n_events)
```

What it means: how much burstiness was flattened away in each matched clade, with larger and more event-rich clades given more weight.

`pulse preservation (overall)`

```text
local_error = 0.35 * mean_emd + 0.55 * mean_burst_loss + 0.10 * mean_centroid_shift
global_error = 0.35 * global_emd + 0.65 * global_burst_loss
pulse_overall = 0.80 * local_error + 0.20 * global_error + 0.20 * (1 - coverage)
```

What it means: a balanced pulse composite combining local clade rhythm, whole-tree rhythm, and how much of the reference pulse panel was actually matched. Here `emd` is the Earth Mover's Distance between relative event-time distributions, and `coverage = matched_clades / panel_clades`.

`pulse preservation (burst)`

```text
local_error_burst = 0.20 * mean_emd + 0.75 * mean_burst_loss + 0.05 * mean_centroid_shift
pulse_burst = 0.80 * local_error_burst + 0.20 * global_error + 0.20 * (1 - coverage)
```

What it means: the same pulse family, but with extra weight placed on keeping burst structure.

`gap layer`

```text
relative_gap_i = (node_age_i - age_min_i) / age_min_i
mean_relative_gap = mean(relative_gap_i)
```

What it means: the average amount of extra inferred lineage history beyond the calibration minima, scaled by the minimum ages. Use it as a core metric only when the calibration ages are primary evidence. With secondary or congruified ages, report it separately as calibration slack or omit it from the core rank.

`rate irregularity`

```text
branch_rate = phylogram_branch_length / dated_branch_duration
rate_irregularity = sd(log_rate) + mean_parent_child_jump + 2 * extreme_rate_frac + autocorr_penalty
autocorr_penalty = 1 - max(rate_autocorr_spearman, 0)
```

What it means: the score rises when branchwise rates are more dispersed, jump more sharply from parent to child, produce more extreme outlier branches, or lose positive autocorrelation.

`uncertainty width` (optional precision layer)

```text
mean_interval_width = mean(width_i)
median_interval_width = median(width_i)
```

What it means: lower values indicate narrower confidence or credible intervals and therefore greater precision. In the Syngnatharia example, these widths come from extracted HPD bars in the published figure rather than from interval metadata embedded in the Newick trees, but PCR can also summarize widths directly from annotated Newick trees when those intervals are present in embedded metadata.

`overall family-balanced rank`

```text
pulse_family_rank = mean(rank(burst_loss), rank(pulse_burst), rank(pulse_overall))
overall_mean_rank = mean(pulse_family_rank, available_nonpulse_family_ranks)
```

What it means: when pulse, gap, and rate are all used, each contributes `1/3` of the final rank. When gap is omitted because it is not independently informative, pulse and rate each contribute `1/2`.

</details>

## Run PCR on your own data

`PhyloChronoRank (PCR)` includes a standalone runner:

- `scripts/run_pcr.R`

At minimum, it expects:

- one reference phylogram
- one candidate manifest with `candidate,tree_file`

Optional:

- `--calibrations-csv` when you want PCR to score the calibration-fit / gap layer
- `--uncertainty-csv` when you want to report the optional precision layer from a separate interval-width table

If you do not provide `--uncertainty-csv`, PCR will still try to extract interval widths directly from annotated Newick trees.

Example commands:

```bash
Rscript scripts/run_pcr.R \
  --ref-tree=examples/terapontoid/Terapontoid_ML_MAIN_phylogram_used.tree \
  --candidates-csv=examples/terapontoid/candidates.csv \
  --outdir=out/terapontoid
```

```bash
Rscript scripts/run_pcr.R \
  --ref-tree=examples/syngnatharia/backbone_Raxml_besttree_matrix75.tre \
  --candidates-csv=examples/syngnatharia/candidates.csv \
  --calibrations-csv=examples/syngnatharia/calibrations_by_candidate.csv \
  --uncertainty-csv=examples/syngnatharia/uncertainty_summary_long.csv \
  --outdir=out/syngnatharia
```

To validate the bundled examples and the displayed README tables, run:

```bash
Rscript scripts/validate_examples.R
```

## Example 1: Empirical dataset with two competing chronograms (Syngnatharia)

### Visual choice before metrics

This example is different. It does not start from a `chronos` fit search. It starts from an earlier visual comparison among the `RAxML` phylogram, `MCMCTree`, and `RelTime`. In that original comparison, the practical choice was to favor `RelTime` because it visually preserved the diversification bursts in the phylogram better than `MCMCTree`, as discussed in [Santaquiteria et al. 2024](https://www.journals.uchicago.edu/doi/10.1086/733931).

That is the key point of this first example: the choice to prefer `RelTime` came first as a visual judgment. The post-fit metrics quantify that visual rationale.

### Quick takeaway

- `RelTime` is the core PCR winner in this comparison
- `RelTime` wins all three pulse summaries and also wins `rate irregularity`
- `MCMCTree` wins the simple calibration-fit layer through lower `mean relative gap`
- `MCMCTree` also has narrower extracted HPD bars, so it wins the optional precision layer on interval width
- `Figure A` is the original visual rationale from the paper; `Figure B` is the quantitative post-fit follow-up

### Figure A: Original visual rationale from the paper

![Syngnatharia paper figure showing the original visual choice](../examples/syngnatharia/Fig_S5_Burst_preservation.png)

This is the original paper figure from [Santaquiteria et al. 2024](https://www.journals.uchicago.edu/doi/10.1086/733931) that motivated the visual preference for `RelTime`. The RAxML phylogram on the left shows clustered branching bursts in several parts of the tree. In the middle panel, `MCMCTree` spreads many of those events out more evenly through time. In the right panel, `RelTime` better tracks the burst structure seen in the phylogram. That was the original rationale back then. This panel addresses the pulse issue only. It does not show the calibration-fit layer or the rate layer.

### Ranked post-fit results (lower is better)

The core PCR rank shown below is family-balanced across `pulse`, `mean relative gap`, and `rate irregularity`, so pulse contributes one-third of the final score. Here the gap layer is informative because the comparison uses primary calibration information summarized from Table S2 rather than secondary congruified ages. The uncertainty-width layer is shown separately as an additional precision consideration, not folded into the core winner, because interval width reflects method-dependent precision rather than the same chronogram-behavior axis captured by pulse, gap, and rate. In this example, that separate reporting appears in two places: the `uncertainty width` column in the table and the `Uncertainty width` panel in Figure B. The uncertainty widths are based on extracted HPD-width spreadsheets because the supplied Newick trees do not themselves contain embedded interval metadata.

| candidate | burst loss | pulse preservation (burst) | pulse preservation (overall) | mean relative gap | rate irregularity | uncertainty width (mean HPD width, Ma) | core overall mean rank (pulse = 1/3) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `RelTime` | `0.1620` | `0.1874` | `0.2077` | `0.2689` | `1.1514` | `11.41` | `1.33` |
| `MCMCTree` | `0.2396` | `0.2560` | `0.2600` | `0.2184` | `1.5531` | `6.24` | `1.67` |

### Figure B: Post-fit comparison across metric families

![Syngnatharia post-fit evaluation metric families](../figures/syngnatharia_postfit_metric_family_values.png)

Figure B shows the core PCR comparison and also displays `uncertainty width` as a separate optional precision layer. The three pulse panels are still averaged into one pulse-family contribution, and in the core PCR rank `pulse`, `mean relative gap`, and `rate irregularity` each contribute one-third.

That separation is deliberate. The RelTime literature does not support a simple story that wider intervals are automatically worse or are just a generic consequence of not using MCMC. Instead, broader RelTime intervals are often discussed as a precision-versus-coverage tradeoff: the analytical confidence-interval procedure explicitly propagates branch-length uncertainty and rate heterogeneity, which can yield wider intervals than other fast methods and, in some scenarios, wider intervals than Bayesian HPDs. Those broader intervals can improve coverage while reducing precision ([Tao et al. 2020](https://academic.oup.com/mbe/article/37/1/280/5602325); [Costa et al. 2022](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-09030-5); [Beavan et al. 2020](https://academic.oup.com/gbe/article/12/7/1087/5842139)). That is why PCR reports interval width here as an additional consideration rather than treating it as a fourth co-equal family in the core rank.

### Interpretation for this example

- `RelTime` is the core PCR winner because it leads all three pulse summaries and also leads `rate irregularity`, while losing only the simple calibration-fit layer
- `MCMCTree` has the lower `mean relative gap`, so it stays closer to the calibration minima on average in this scoring
- `MCMCTree` also has the narrower extracted HPD bars, so it wins the optional uncertainty-width layer on precision
- the post-fit metrics therefore support, rather than reverse, the original visual rationale from Figure A: `RelTime` better preserves the branching bursts seen in the RAxML phylogram
- this is exactly the kind of case where a visual choice made before these metrics existed can now be quantified explicitly instead of being left as impression only

### Practical decision rule

1. If you want the core PCR winner focused on chronogram behavior, choose `RelTime`.
2. If you care most about preserving diversification bursts and smoother implied rate behavior, choose `RelTime`.
3. If you care primarily about calibration fit and narrower interval estimates, `MCMCTree` wins `mean relative gap` and the optional uncertainty-width layer.
4. Report the tradeoff explicitly: here the core chronogram-behavior layer favors `RelTime`, while the calibration-plus-precision side favors `MCMCTree`.
5. The original visual choice to favor `RelTime` is supported quantitatively by the current core post-fit layer.

<details>
<summary><strong>Files behind this example</strong></summary>

- `examples/syngnatharia/candidates.csv`
- `examples/syngnatharia/calibrations_by_candidate.csv`
- `examples/syngnatharia/backbone_Raxml_besttree_matrix75.tre`
- `examples/syngnatharia/CalibratedTree_backbone_MCMCTree_matrix75_RAXML.tre`
- `examples/syngnatharia/CalibratedTree_backbone_RelTime_matrix75_RAXML.tre`
- `examples/syngnatharia/Fig_S5_Burst_preservation.png`
- `examples/syngnatharia/syngnatharia_HPD_summary.csv`
- `examples/syngnatharia/syngnatharia_HPD_widths_extracted.csv`
- `examples/syngnatharia/uncertainty_summary_long.csv`
- `examples/syngnatharia/postfit_metrics/syngnatharia_postfit_metrics.csv`
- `examples/syngnatharia/postfit_metrics/syngnatharia_fossil_gap_side_by_side.csv`
- `examples/syngnatharia/postfit_metrics/syngnatharia_tableS2_method_audit.csv`
- `figures/syngnatharia_postfit_metric_family_values.png`
- `scripts/run_pcr.R`
- `scripts/make_syngnatharia_postfit_figures.R`

</details>

## Example 2: Empirical dataset with six competing chronograms (Terapontoidei)

### Optional upstream fit context

This tab does not do model fitting; that workflow lives in tab 1. In this example, the upstream fit statistics and the post-fit results still point to the same two `chronos` models. `clock` has the best `PHIIC` in the fit summary. `discrete` has the best penalized log-likelihood. Under the core PCR comparison, `clock` is the strongest balanced `chronos` tree and `discrete` is the close runner-up.

### Quick takeaway

- the core PCR comparison in this example uses `pulse` and `rate irregularity`, not `gap burden`, because the dates come from secondary / congruified calibration ages rather than primary calibration evidence
- `chronos_clock` is the clear winner under the two-family core rank
- `chronos_discrete` is the close runner-up
- `RelTime` is the strongest pulse-preservation candidate, but it pays a large `rate irregularity` penalty
- `chronos_clock` is the best tree for `rate irregularity`
- the optional uncertainty layer is available here for all six trees on one shared bootstrap scale
- `treePL` is a middle-tier candidate rather than a leading tree in this example

In this bundled six-tree comparison, the selected `treePL` candidate uses `smooth = 0.01`; it is labeled simply `treePL` below because only one `treePL` tree is carried forward into the example.

### Figure A: Pulse-layer tree-shape comparison among bundled chronos trees

![Pulse preservation tree panel](../figures/branching_tempo_tree_panel_clean_v3.png)

This figure shows the pulse layer directly on alternative `chronos` trees (estimated with different clock models); `treePL` and `RelTime` are not shown here. It helps explain the upstream `chronos` pulse tradeoffs, but the quantitative post-fit comparison below also includes `treePL` and `RelTime`. This figure is only to illustrate the pulse issue; `gap burden` is not computed in this example, and the panel does not show the `rate irregularity` part of the broader PCR toolkit.

### Ranked post-fit results (lower is better)

`gap burden` is not computed in this example because the trees were dated using congruified / secondary calibration ages rather than primary calibration evidence. In that setting, a gap score would just measure how closely each method reproduced inherited ages rather than providing an independent biological diagnostic.

The overall mean rank below is therefore family-balanced across pulse and rate only. The three pulse summaries are shown separately for transparency, but they are first collapsed into one pulse-family contribution. So pulse as a whole contributes one-half of the final overall rank, and `rate irregularity` contributes the other half. The optional uncertainty-width layer is reported separately as a precision check and is not folded into the core rank. The last column reports that mean-rank value itself (`rank_mean_core`), not the separate ordinal finish position (`rank_mean_core_rank`).

The uncertainty-width contrasts here should be read cautiously. The shared table uses one bootstrap comparison scale: the four `chronos` rows come from the vendored parametric bootstrap helper of [Paradis et al. 2023](https://doi.org/10.1016/j.ympev.2022.107652), while `treePL` and `RelTime` come from repo-local bootstrap reruns under the same branch-length resampling design. The [Tao et al. 2020](https://academic.oup.com/mbe/article/37/1/280/5602325) analytical `RelTime` CI is supplemental and not reported in the table because, on these hard-bounded empirical trees, it lives in a completely different numerical universe from the bootstrap widths: after full-bound projection compresses internal durations, the analytical variance term can explode by orders of magnitude.

| candidate | burst loss | pulse preservation (burst) | pulse preservation (overall) | rate irregularity | uncertainty width (mean CI width, Ma) | overall mean rank (pulse = 1/2) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `chronos_clock` | `0.1348` | `0.1464` | `0.1604` | `2.5897` | `3.39` | `1.50` |
| `chronos_discrete` | `0.1349` | `0.1464` | `0.1605` | `2.5946` | `3.26` | `2.50` |
| `treePL` | `0.2492` | `0.2083` | `0.1995` | `3.1294` | `8.01` | `3.50` |
| `RelTime` | `0.1145` | `0.1328` | `0.1501` | `6.6066` | `8.70` | `3.50` |
| `chronos_relaxed` | `0.2797` | `0.2275` | `0.2140` | `4.6401` | `2.56` | `5.00` |
| `chronos_correlated` | `0.2797` | `0.2275` | `0.2140` | `4.6401` | `1.81` | `5.00` |

In short: `RelTime` minimizes all three pulse summaries, `chronos_clock` leads `rate irregularity`, and `chronos_clock` has a clean edge over `chronos_discrete` in the two-family core rank. `treePL` lands in the same mean-rank tier as `RelTime`, but for the opposite reason: better rate behavior with much weaker pulse preservation. On the optional precision layer, the four `chronos` trees are the narrowest on the shared bootstrap scale, `treePL` is broader than all four `chronos` trees but slightly narrower than `RelTime`, and `chronos_correlated` is the narrowest overall.

### Figure B: Post-fit comparison across metric families

![Post-fit evaluation metric families](../figures/postfit_metric_family_values.png)

Figure B uses the same family-balanced rule as the table. Even though three pulse panels are shown, they do not count as three separate halves. They are averaged into one pulse-family contribution, and that pulse family contributes one-half of the overall rank. `gap burden` is intentionally absent here because the calibration ages are secondary / congruified rather than primary evidence. `Uncertainty width` is shown separately as an optional precision layer.

### Interpretation for this example

- `chronos_clock` is the core PCR winner
- `chronos_discrete` is the runner-up
- `RelTime` is the strongest tree on all three pulse summaries
- `chronos_clock` is the best tree on `rate irregularity`
- the optional uncertainty layer is led by the non-clock `chronos` trees on width alone, with `chronos_correlated` narrowest and `chronos_relaxed` next
- `treePL` beats both non-clock `chronos` trees on the core comparison, but it still trails `chronos_clock` and `chronos_discrete`
- `treePL` is broader than all four `chronos` trees on the shared bootstrap layer, but slightly narrower than `RelTime`
- `RelTime` has broader bootstrap widths than any of the `chronos` candidates in this example, so it loses the optional precision layer
- the supplemental Tao-style `RelTime` CI file is not mixed into the table because on this hard-bounded point tree it lives in a completely different numerical regime from the bootstrap widths
- `chronos_correlated` and `chronos_relaxed` are the weakest candidates in the set under the post-fit layer

### Practical decision rule

1. If you want the strongest pulse-preservation candidate, choose `RelTime`.
2. If you want the smoothest implied rate behavior, choose `chronos_clock`.
3. If you want the narrowest bundled bootstrap intervals, `chronos_correlated` and `chronos_relaxed` are narrower on width alone, but they perform poorly on the core post-fit comparison.
4. If you want one concise core-PCR statement, report `chronos_clock` as the winner under the two-family comparison, with `chronos_discrete` as the close runner-up.
5. If an upstream fit-based selector and PCR point to different trees, report both explicitly rather than collapsing them into one claim.
6. In this example, `RelTime` behaves like a pulse specialist, while `treePL` acts as a mid-ranking rate-friendlier alternative rather than a leading tree. On the shared bootstrap uncertainty layer, `treePL` sits between the narrower `chronos` trees and the broader `RelTime` tree. The supplemental Tao-style `RelTime` CI file is discussed separately because it is not on the same numerical scale as the shared bootstrap summaries.

<details>
<summary><strong>Files behind this example</strong></summary>

- `examples/terapontoid/summary_terap_empirical_model_fits.csv`
- `examples/terapontoid/summary_terap_empirical_postfit_metrics.csv`
- `examples/terapontoid/uncertainty_summary_long.csv`
- `examples/terapontoid/candidates.csv`
- `examples/terapontoid/Terapontoid_ML_MAIN_calibrations_used.csv` (upstream congruified calibration ages used in dating; not scored in the core PCR rank)
- the six trees in `examples/terapontoid/`, including `Terapontoid_ML_MAIN_treePL_congruify.tre` and `Terapontoid_ML_MAIN_RelTime_full_bounds.tre`
- `examples/terapontoid/Terapontoid_ML_MAIN_chronos_dated_modelclock_ci.csv`
- `examples/terapontoid/Terapontoid_ML_MAIN_chronos_dated_modeldiscrete_ci.csv`
- `examples/terapontoid/Terapontoid_ML_MAIN_chronos_dated_modelcorrelated_ci.csv`
- `examples/terapontoid/Terapontoid_ML_MAIN_chronos_dated_modelrelaxed_ci.csv`
- `examples/terapontoid/Terapontoid_ML_MAIN_RelTime_bounds_used.csv`
- `examples/terapontoid/Terapontoid_ML_MAIN_treePL_congruify_bootstrap_ci.csv`
- `examples/terapontoid/Terapontoid_ML_MAIN_RelTime_full_bounds_bootstrap_ci.csv`
- `examples/terapontoid/Terapontoid_ML_MAIN_RelTime_full_bounds_ci.csv`
- `figures/branching_tempo_tree_panel_clean_v3.png`
- `figures/postfit_metric_family_values.png`
- `scripts/run_pcr.R`
- `scripts/chronos_ci_helpers.R`
- `scripts/reltime_helpers.R`
- `scripts/build_uncertainty_examples.R`
- `scripts/make_terapontoid_postfit_figures.R`
- `scripts/make_terapontoid_pulse_tree_panel.R`

</details>

## Example 3: Unpublished vertebrate dataset (derived outputs only)

### Quick takeaway

- `chronos_clock` is the core PCR winner in this comparison
- `treePL` and `RelTime` land in the next mean-rank tier for opposite reasons
- `RelTime` dominates the pulse and gap layers but is the worst tree on `rate irregularity`
- the optional uncertainty layer is available here for all six trees on one shared bootstrap scale
- `chronos_correlated`, `chronos_relaxed`, and `chronos_discrete` are all much worse on pulse, gap, and rate
- the raw trees and calibration table are not distributed here because this dataset is unpublished

### Selected candidates

This example uses `57` calibrations and compares six selected chronograms:

- `chronos_clock` with `lambda = 1`
- `chronos_correlated` with `lambda = 0.1`
- `chronos_relaxed` with `lambda = 0.1`
- `chronos_discrete` with `lambda = 0.1` and `nb_rate_cat = 5`
- `treePL` with best `smooth = 100`
- `RelTime`, projected onto the same merged full calibration bounds used by the local `chronos` and `treePL` pipelines

Only derived outputs are shown in this repository. The raw input trees and calibration table are withheld because the dataset is unpublished.

### Figure A: Relative tree-shape comparison across the selected chronograms

![Unpublished vertebrate tree panel](../figures/unpublished_vertebrate_tree_panel.png)

This panel compares the reference phylogram and the five originally selected chronograms after scaling each tree to its own maximum root-to-tip depth. The point is to show relative branching tempo, not absolute branch-length units. Tip labels are hidden, and the raw trees are not distributed. The added `RelTime` benchmark appears in the quantitative comparison below rather than in this legacy tree panel.

### Ranked post-fit results (lower is better)

The core PCR rank is family-balanced across `pulse`, `mean relative gap`, and `rate irregularity`, so pulse contributes one-third of the final score. The optional uncertainty-width layer is reported separately as a precision check and is not folded into the core rank.

The shared uncertainty table here uses one bootstrap comparison scale. The four `chronos` rows come from the vendored parametric bootstrap helper of [Paradis et al. 2023](https://doi.org/10.1016/j.ympev.2022.107652), while `treePL` and `RelTime` come from repo-local bootstrap reruns under the same branch-length resampling design. The [Tao et al. 2020](https://academic.oup.com/mbe/article/37/1/280/5602325) analytical `RelTime` CI is supplemental and not reported in the table because, on these hard-bounded empirical trees, it can live in a completely different numerical universe from the bootstrap widths after bound projection compresses internal durations.

| candidate | burst loss | pulse preservation (burst) | pulse preservation (overall) | mean relative gap | rate irregularity | uncertainty width (mean CI width, Ma) | core overall mean rank (pulse = 1/3) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `chronos_clock` | `0.1606` | `0.2106` | `0.2233` | `0.1235` | `1.5915` | `30.04` | `1.67` |
| `treePL` | `0.1646` | `0.2209` | `0.2332` | `0.1548` | `1.7128` | `11.14` | `2.67` |
| `RelTime` | `0.0449` | `0.1164` | `0.1539` | `0.0072` | `6.3142` | `6.17` | `2.67` |
| `chronos_discrete` | `0.3468` | `0.3172` | `0.2909` | `0.4129` | `3.5759` | `23.33` | `4.00` |
| `chronos_correlated` | `0.3547` | `0.3268` | `0.2989` | `0.3951` | `3.6171` | `18.68` | `4.33` |
| `chronos_relaxed` | `0.3635` | `0.3302` | `0.3009` | `0.4377` | `3.6971` | `22.43` | `5.67` |

### Figure B: Post-fit comparison across metric families

![Unpublished vertebrate post-fit metric families](../figures/unpublished_vertebrate_postfit_metric_family_values.png)

Figure B uses the same family-balanced rule as the table. The three pulse panels are shown separately for transparency, but together they count as one pulse family. `Uncertainty width` is shown separately as an optional precision layer.

### Interpretation for this example

- `chronos_clock` is the core PCR winner in this comparison
- `treePL` and `RelTime` share the next mean-rank tier, but for opposite reasons
- `RelTime` is the strongest tree on all three pulse summaries and also has the smallest `mean relative gap`
- `treePL` has the cleanest `rate irregularity` in the set while staying close to `chronos_clock` on the pulse layer
- `chronos_discrete` is a weak candidate, and `chronos_correlated` and `chronos_relaxed` are no better
- `chronos_correlated` and `chronos_relaxed` both score poorly across the full post-fit layer, especially on `rate irregularity` and `mean relative gap`
- the bundled uncertainty layer is available here for all six trees, with `treePL` intermediate on the shared bootstrap precision layer
- `RelTime` is the pulse-plus-gap specialist and also the narrowest candidate on the shared bootstrap precision layer
- the supplemental Tao-style `RelTime` CI file is not mixed into the table because on this hard-bounded point tree it lives in a completely different numerical regime from the bootstrap widths
- in this dataset, the post-fit layer supports `chronos_clock` as the best balanced tree, `treePL` as the rate-friendlier close alternative, and `RelTime` as the pulse-plus-gap specialist

### Practical decision rule

1. If you want one core PCR winner in this comparison, choose `chronos_clock`.
2. If you want the closest balanced alternative under the current bundled outputs, `treePL` is the next-best point-estimate candidate.
3. If your priority is preserving branching tempo and staying closest to the calibration minima, consider `RelTime`, but report its poor `rate irregularity` explicitly.
4. If you discuss uncertainty in the bundled release, note that the shared table compares `chronos`, `treePL`, and `RelTime` bootstrap widths, while the Tao-style `RelTime` CI is kept supplemental only.
5. If you report multiple candidate chronograms, the main contrast is `chronos_clock` as the balanced winner, `treePL` as the close rate-friendlier alternative, and `RelTime` as the pulse-plus-gap alternative.

<details>
<summary><strong>Files behind this example</strong></summary>

- `examples/unpublished_vertebrate/postfit_metrics/summary_unpublished_vertebrate_postfit_metrics.csv`
- `examples/unpublished_vertebrate/postfit_metrics/uncertainty_summary_long.csv`
- `figures/unpublished_vertebrate_tree_panel.png`
- `figures/unpublished_vertebrate_postfit_metric_family_values.png`
- `scripts/chronos_ci_helpers.R`
- `scripts/reltime_helpers.R`
- `scripts/build_uncertainty_examples.R`
- `scripts/make_unpublished_vertebrate_postfit_figure.R`

</details>

<details>
<summary><strong>Scope notes</strong></summary>

- The pulse-family weights are user-chosen defaults. A small fixed robustness check across five perturbation sets is included in `examples/weight_sensitivity/`; neither the pulse-family winner nor the core PCR winner changed in either bundled example. A broader sensitivity analysis across additional datasets has not yet been done.
- The pulse family treats the source phylogram as the reference for branching rhythm. That is useful when the question is whether a dated tree preserves the tempo structure visible in the starting phylogram, but it should not be read as proof that the phylogram itself is the true diversification history.
- `rate irregularity` is useful for comparing candidates within the same dataset, but the absolute values are not meant to be compared across unrelated datasets.
- Under point calibrations, the gap layer behaves as symmetric calibration slack: older and younger deviations from the calibration point are penalized equally.
- `gap burden` should be treated as a core family only when the calibration ages are primary evidence. With secondary or congruified ages, the same calculation becomes circular calibration slack and is better reported separately or omitted from the core rank.
- PCR reports raw scores and ranks. It does not yet attach bootstrap or permutation p-values to score differences.
- The framework evaluates point chronograms. It does not yet propagate posterior tree uncertainty through the post-fit scores.
- The optional `uncertainty width` layer is reported separately from the core PCR rank. In the bundled examples it comes from extracted HPD widths in Syngnatharia and from matched `chronos`, `treePL`, and `RelTime` bootstrap summaries in Terapontoidei and the bundled vertebrate release. The Tao-style analytical `RelTime` CI is kept as a supplemental diagnostic and discussed in prose when it diverges sharply from the bootstrap scale. This layer speaks to precision, not accuracy.

</details>

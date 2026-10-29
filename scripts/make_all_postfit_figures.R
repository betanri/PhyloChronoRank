## Generate the post-fit 3-family figures for the README examples (Part 2).
## Inputs are the pcr_rebuild / pcr_sparse6 / pcr_dense19 outputs; datasets
## whose outputs are not in the repo (unpublished) are skipped.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
if (length(file_arg)) {
  script_dir <- dirname(normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE))
} else {
  script_dir <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
}
base_dir <- normalizePath(file.path(script_dir, '..'), winslash = '/', mustWork = TRUE)
out_fig <- file.path(base_dir, 'figures')
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

## --- Helper ---
plot_panel <- function(vals, labels, cols, ttl, ylab, fam_label = '', higher_better = FALSE) {
  ylim_top <- max(vals, na.rm = TRUE) * 1.30
  if (all(is.na(vals) | vals == 0)) ylim_top <- 1
  bp <- barplot(vals, names.arg = labels, col = cols, las = 2,
                ylab = ylab, main = ttl,
                ylim = c(0, ylim_top))
  text(bp, vals, labels = sprintf('%.3f', vals), pos = 3, cex = 0.95)
  orient <- if (higher_better) 'Higher is better' else 'Lower is better'
  mtext(paste0(orient, '   ', fam_label), side = 1, line = 6.2, cex = 0.85)
}

## Sources that represent a method's NATIVE confidence intervals (i.e., what
## the dating method itself emits — not a synthetic post-hoc bootstrap added
## by PCR for chronos/treePL, which produces near-zero widths that would be
## misleading next to genuine CIs). Candidates outside this whitelist are
## excluded from the uncertainty-width panel.
.native_ci_sources <- function(sources) {
  s <- ifelse(is.na(sources), '', as.character(sources))
  s == 'tree_embedded_intervals' |
    grepl('^RelTime', s, ignore.case = TRUE)
}

make_figure <- function(csv_path, fig_path, main_title, short_labels, col_map, unit = 'Ma') {
  if (!file.exists(csv_path)) { message('Skip: ', csv_path); return(invisible(NULL)) }
  d_raw <- read.csv(csv_path, stringsAsFactors = FALSE)

  ## Core-metric panels (Families 1-3, family scores, overall rank) deduplicate
  ## rows that share identical core metrics: multiple RelTime CI variants score
  ## identically on core PCR and would produce identical bars in those panels.
  fp <- paste(d_raw$burst_loss, d_raw$internode_concordance,
              d_raw$rate_irregularity, sep = '|')
  d <- d_raw[!duplicated(fp), , drop = FALSE]
  ord <- d$candidate[order(d$rank_mean_core_rank, d$rank_mean_core)]
  d <- d[match(ord, d$candidate), ]
  labels <- vapply(ord, function(x) if (x %in% names(short_labels)) short_labels[[x]] else x, '')
  cols <- vapply(ord, function(x) if (x %in% names(col_map)) col_map[[x]] else 'grey60', '')

  ## Uncertainty-width panel uses a SEPARATE subset built from the raw CSV so
  ## that multiple CI variants for the same base method (e.g. RelTime Tao vs
  ## RelTime bootstrap) are both plotted side-by-side. Only candidates whose
  ## uncertainty comes from a method-native CI source are kept.
  ## Older summary CSVs (e.g. syngnatharia) predate the uncertainty_source
  ## column and only record tree-embedded CIs, so every finite width is
  ## effectively a method-native CI in that era.
  if ('uncertainty_source' %in% names(d_raw)) {
    d_unc <- d_raw[.native_ci_sources(d_raw$uncertainty_source), , drop = FALSE]
  } else {
    d_unc <- d_raw
  }
  d_unc <- d_unc[is.finite(d_unc$uncertainty_mean_width_ma) &
                   d_unc$uncertainty_mean_width_ma > 0, , drop = FALSE]
  ## Drop candidates whose CI width is <1% of the widest CI in the retained
  ## set: such bars render as visually invisible at the panel scale and
  ## misleadingly suggest "zero uncertainty" when in fact the CI procedure
  ## failed or collapsed (typical for strict-clock bootstraps on deep trees).
  if (nrow(d_unc)) {
    uw_max <- max(d_unc$uncertainty_mean_width_ma, na.rm = TRUE)
    d_unc <- d_unc[d_unc$uncertainty_mean_width_ma >= 0.01 * uw_max, , drop = FALSE]
  }
  ord_unc <- d_unc$candidate
  lab_unc <- vapply(ord_unc, function(x) if (x %in% names(short_labels)) short_labels[[x]] else x, '')
  col_unc <- vapply(ord_unc, function(x) if (x %in% names(col_map)) col_map[[x]] else 'grey60', '')

  png(fig_path, width = 3400, height = 2000, res = 170)
  layout(matrix(c(1,2,3,4, 5,6,7,8, 9,9,10,10), nrow = 3, byrow = TRUE),
         heights = c(1, 1, 1))
  par(mar = c(8, 5, 4, 1), oma = c(2.6, 0.2, 2.2, 0.2))

  ## Row 1: Family 1 (2 metrics) + Family 2 (2 metrics)
  plot_panel(d$burst_loss, labels, cols, 'Burst loss', 'Burst loss', 'Family 1: Radiation-zone fidelity')
  plot_panel(d$internode_concordance, labels, cols, 'Internode concordance', 'Concordance', 'Family 1: Radiation-zone fidelity', higher_better = TRUE)
  plot_panel(d$compression_score, labels, cols, 'Compression score', 'Compression', 'Family 2: Global fidelity')
  plot_panel(d$depth_r2, labels, cols, 'Depth R\u00B2', 'R\u00B2', 'Family 2: Global fidelity', higher_better = TRUE)
  ## Row 2: Family 3 (2 metrics) + uncertainty + empty
  plot_panel(d$rate_irregularity, labels, cols, 'Rate irregularity', 'Rate irregularity', 'Family 3: Rate & Calibration')
  vals_gap <- d$mean_relative_gap
  if (all(is.na(vals_gap))) {
    plot.new(); text(0.5, 0.5, 'Mean relative gap\nnot scored\n(all calibrations are fixed ages)', cex = 1.2)
  } else {
    plot_panel(vals_gap, labels, cols, 'Mean relative gap', 'Relative gap', 'Family 3: Rate & Calibration')
  }
  if (!nrow(d_unc)) {
    plot.new()
    text(0.5, 0.5,
         'Uncertainty width\nnot available\n(no method-native CIs)',
         cex = 1.2)
  } else {
    uw_ci   <- d_unc$uncertainty_mean_width_ma
    ylim_top <- max(uw_ci, na.rm = TRUE) * 1.30
    if (!is.finite(ylim_top) || ylim_top == 0) ylim_top <- 1
    bp <- barplot(uw_ci, names.arg = lab_unc, col = col_unc, las = 2,
                  ylab = paste0('Mean 95% CI width (', unit, ')'), main = 'Uncertainty width',
                  ylim = c(0, ylim_top))
    text(bp, uw_ci, labels = sprintf('%.3f', uw_ci), pos = 3, cex = 0.95)
    mtext('Lower is better   Method-native CIs only', side = 1, line = 6.2, cex = 0.85)
  }
  plot.new()  ## empty panel to fill row
  ## Row 3
  par(mar = c(8, 5, 4, 1))
  fam_vals <- rbind(d$family_radiation_score, d$family_global_score, d$family_ratecal_score)
  colnames(fam_vals) <- labels
  rownames(fam_vals) <- c('Fam 1: Radiation', 'Fam 2: Global', 'Fam 3: Rate+Cal')
  fam_cols <- c('#e41a1c', '#377eb8', '#4daf4a')
  bp <- barplot(fam_vals, beside = TRUE, col = fam_cols, las = 2,
                ylab = 'Family rank (lower is better)', main = 'Family scores',
                ylim = c(0, max(fam_vals, na.rm = TRUE) * 1.35),
                legend.text = rownames(fam_vals),
                args.legend = list(x = 'topright', cex = 0.85, bty = 'n'))
  text(bp, fam_vals, labels = sprintf('%.1f', fam_vals), pos = 3, cex = 0.8)
  mtext('Lower is better', side = 1, line = 6.2, cex = 0.85)

  bp2 <- barplot(d$rank_mean_core, names.arg = labels, col = cols, las = 2,
                 ylab = 'Mean rank across 3 families', main = 'Core overall rank (3-family balanced)',
                 ylim = c(0, max(d$rank_mean_core, na.rm = TRUE) * 1.4))
  text(bp2, d$rank_mean_core,
       labels = paste0(sprintf('%.2f', d$rank_mean_core), ' (rank ', d$rank_mean_core_rank, ')'),
       pos = 3, cex = 0.9)
  mtext('Lower is better', side = 1, line = 6.0, cex = 0.85)

  mtext(main_title, side = 3, outer = TRUE, line = 0.5, cex = 1.5, font = 2)
  mtext('Core PCR rank balances radiation-zone fidelity, global chronogram fidelity, and rate+calibration. Mean relative gap NA when all calibrations are fixed ages. Uncertainty width is not part of the core rank.',
        side = 1, outer = TRUE, line = 0.5, cex = 0.95)
  dev.off()
  message('Wrote: ', normalizePath(fig_path, winslash = '/'))
}

## ---- Labels and colours by method ----
## Candidate names carry their settings (e.g. chronos_clock_lambda1,
## treePL_smooth1000), so labels and colours follow the method prefix.
method_label <- function(x) {
  rules <- c('^chronos_clock' = 'Clock', '^chronos_correlated' = 'Correlated',
             '^chronos_relaxed' = 'Relaxed', '^chronos_discrete' = 'Discrete',
             '^treePL' = 'treePL', '^RelTime' = 'RelTime', '^MDCat' = 'MD-Cat',
             '^MCMCTree' = 'MCMCTree')
  for (k in names(rules)) if (grepl(k, x)) return(rules[[k]])
  x
}
method_colours <- c(Clock = '#1b9e77', Correlated = '#d95f0e', Relaxed = '#7570b3',
                    Discrete = '#2c7fb8', treePL = '#6baed6', RelTime = '#d7301f',
                    'MD-Cat' = '#e7298a', MCMCTree = '#525252')

figure_for <- function(csv_rel, png_name, title, unit = 'Ma') {
  csv <- if (grepl('^/', csv_rel)) csv_rel else file.path(base_dir, csv_rel)
  if (!file.exists(csv)) { message('Skip (not in repo): ', csv_rel); return(invisible(NULL)) }
  cands <- read.csv(csv, stringsAsFactors = FALSE)$candidate
  short <- vapply(cands, method_label, '')
  cols <- unname(method_colours[short])
  cols[is.na(cols)] <- 'grey60'
  make_figure(csv, file.path(out_fig, png_name), title,
              setNames(short, cands), setNames(cols, cands), unit = unit)
}

## ---- Datasets (README Examples 1-6) ----
figure_for('examples/syngnatharia/pcr_rebuild/summary_pcr_metrics.csv',
           'syngnatharia_postfit_metric_family_values.png',
           'Syngnatharia (fossil calibrations): post-fit evaluation, 3 families')
figure_for('examples/terapontoid/pcr_sparse6/summary_pcr_metrics.csv',
           'terapontoid_postfit_FEW.png',
           'Terapontoidei (6 congruified ages): post-fit evaluation, 3 families')
figure_for('examples/terapontoid/pcr_dense19/summary_pcr_metrics.csv',
           'terapontoid_postfit_MANY.png',
           'Terapontoidei (19 congruified ages): post-fit evaluation, 3 families')
figure_for('5_DATASETS_FEW_CALS/vertebrate/pcr_rebuild/summary_pcr_metrics.csv',
           'vertebrate_postfit_FEW.png',
           'Vertebrate (18 interval calibrations): post-fit evaluation, 3 families')
figure_for('5_DATASETS_MANY_CALS/vertebrate/pcr_rebuild/summary_pcr_metrics.csv',
           'vertebrate_postfit_MANY.png',
           'Vertebrate (57 interval calibrations): post-fit evaluation, 3 families')
figure_for('GOBIES/GOBY_TREES/McCraney_minus_kurtiform_cals_plus_eosphaeramia_apogoninae/pcr_rebuild/summary_pcr_metrics.csv',
           'gobiaria_postfit.png',
           'Gobiaria (75 congruified calibrations): post-fit evaluation, 3 families')
figure_for('5_DATASETS_FEW_CALS/tetraodontiformes/pcr_rebuild/summary_pcr_metrics.csv',
           'tetraodontiformes_postfit_FEW.png',
           'Tetraodontiformes (38 congruified ages): post-fit evaluation, 3 families')
figure_for('5_DATASETS_MANY_CALS/tetraodontiformes/pcr_rebuild/summary_pcr_metrics.csv',
           'tetraodontiformes_postfit_MANY.png',
           'Tetraodontiformes (120 congruified ages): post-fit evaluation, 3 families')
figure_for(file.path(base_dir, '..', 'Ostario', 'run1_plusAfro_plusCall', 'pcr_rebuild', 'summary_pcr_metrics.csv'),
           'ostariophysi_postfit_v2.png',
           'Otophysi (56 fossil calibrations): post-fit evaluation, 3 families')

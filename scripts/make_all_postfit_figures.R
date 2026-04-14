## Generate postfit 3-family figures for all datasets (vertebrate, gobiaria, tetraodontiformes)
## Each dataset produces a FEW and MANY figure.

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

make_figure <- function(csv_path, fig_path, main_title, short_labels, col_map) {
  if (!file.exists(csv_path)) { message('Skip: ', csv_path); return(invisible(NULL)) }
  d <- read.csv(csv_path, stringsAsFactors = FALSE)
  # Deduplicate by identical metric fingerprints, but keep rows that differ in

  # uncertainty width (so Tao-CI and bootstrap-CI variants both appear).
  fp_cols <- c('burst_loss', 'internode_concordance', 'rate_irregularity')
  uw <- d$uncertainty_mean_width_ma
  uw[is.na(uw)] <- -999          # treat NA as its own level
  fp <- paste(d$burst_loss, d$internode_concordance, d$rate_irregularity, round(uw, 4), sep = '|')
  d <- d[!duplicated(fp), ]
  ord <- d$candidate[order(d$rank_mean_core_rank, d$rank_mean_core)]
  d <- d[match(ord, d$candidate), ]
  labels <- vapply(ord, function(x) if (x %in% names(short_labels)) short_labels[[x]] else x, '')
  cols <- vapply(ord, function(x) if (x %in% names(col_map)) col_map[[x]] else 'grey60', '')

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
    plot.new(); text(0.5, 0.5, 'Mean relative gap\nnot scored\n(secondary calibrations)', cex = 1.2)
  } else {
    plot_panel(vals_gap, labels, cols, 'Mean relative gap', 'Relative gap', 'Family 3: Rate & Calibration')
  }
  uw <- d$uncertainty_mean_width_ma
  if (all(is.na(uw))) {
    plot.new(); text(0.5, 0.5, 'Uncertainty width\nnot available', cex = 1.2)
  } else {
    ## Only plot candidates that have CI data (non-NA uncertainty).
    ## Candidates without CI-embedded trees should not appear here.
    has_ci <- !is.na(uw)
    if (sum(has_ci) == 0) {
      plot.new(); text(0.5, 0.5, 'Uncertainty width\nnot available', cex = 1.2)
    } else {
      uw_ci   <- uw[has_ci]
      lab_ci  <- labels[has_ci]
      col_ci  <- cols[has_ci]
      ylim_top <- max(uw_ci, na.rm = TRUE) * 1.30
      if (ylim_top == 0) ylim_top <- 1
      bp <- barplot(uw_ci, names.arg = lab_ci, col = col_ci, las = 2,
                    ylab = 'Mean CI width (Ma)', main = 'Uncertainty width',
                    ylim = c(0, ylim_top))
      text(bp, uw_ci, labels = sprintf('%.3f', uw_ci), pos = 3, cex = 0.95)
      mtext('Lower is better   Optional precision layer', side = 1, line = 6.2, cex = 0.85)
    }
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
  mtext('Core PCR rank balances radiation-zone fidelity, global chronogram fidelity, and rate+calibration. mean_relative_gap NA when secondary cals. Uncertainty width separate.',
        side = 1, outer = TRUE, line = 0.5, cex = 0.95)
  dev.off()
  message('Wrote: ', normalizePath(fig_path, winslash = '/'))
}

## ---- Short labels and color maps ----
## Common building blocks
short_vert <- c(
  RelTime_MEGA = 'RelTime', RelTime_MEGA_with_bootstrap_CI = 'RelTime(Boot)',
  RelTime_MEGA_with_tao_CI = 'RelTime(Tao)',
  chronos_clock_lambda100 = 'Clock', chronos_clock_lambda0p1 = 'Clock',
  chronos_discrete_lambda100_k2 = 'Discrete', chronos_discrete_lambda0p01_k2 = 'Discrete',
  chronos_correlated_lambda0p01 = 'Correlated',
  chronos_relaxed_lambda0p01 = 'Relaxed', chronos_relaxed_lambda0p1 = 'Relaxed',
  treePL_smooth0p01 = 'treePL'
)
col_vert <- c(
  RelTime_MEGA = '#d7301f', RelTime_MEGA_with_bootstrap_CI = '#fc8d59',
  RelTime_MEGA_with_tao_CI = '#d7301f',
  chronos_clock_lambda100 = '#1b9e77', chronos_clock_lambda0p1 = '#1b9e77',
  chronos_discrete_lambda100_k2 = '#2c7fb8', chronos_discrete_lambda0p01_k2 = '#2c7fb8',
  chronos_correlated_lambda0p01 = '#d95f0e',
  chronos_relaxed_lambda0p01 = '#7570b3', chronos_relaxed_lambda0p1 = '#7570b3',
  treePL_smooth0p01 = '#6baed6'
)

short_gobi <- c(
  chronos_clock = 'Clock', chronos_correlated = 'Correlated',
  chronos_discrete = 'Discrete', treePL = 'treePL',
  RelTime_with_tao_CI = 'RelTime(Tao)', RelTime_with_bootstrap_CI = 'RelTime(Boot)'
)
col_gobi <- c(
  chronos_clock = '#1b9e77', chronos_correlated = '#d95f0e',
  chronos_discrete = '#2c7fb8', treePL = '#6baed6',
  RelTime_with_tao_CI = '#d7301f', RelTime_with_bootstrap_CI = '#fc8d59'
)

short_tforms <- c(
  RelTime_MEGA = 'RelTime', RelTime_MEGA_with_bootstrap_CI = 'RelTime(Boot)',
  RelTime_MEGA_with_tao_CI = 'RelTime(Tao)',
  chronos_clock_lambda1 = 'Clock',
  chronos_discrete_lambda1_k2 = 'Discrete',
  chronos_correlated_lambda0p01 = 'Correlated',
  chronos_relaxed_lambda0p01 = 'Relaxed',
  treePL_smooth0p01 = 'treePL'
)
col_tforms <- col_vert  # same colors, extended
col_tforms['chronos_clock_lambda1'] <- '#1b9e77'
col_tforms['chronos_discrete_lambda1_k2'] <- '#2c7fb8'

## ---- Syngnatharia ----
short_syn <- c(RelTime = 'RelTime', MCMCTree = 'MCMCTree')
col_syn <- c(RelTime = '#d7301f', MCMCTree = '#2c7fb8')
make_figure(
  file.path(base_dir, 'examples', 'syngnatharia', 'pcr_rerun_3fam', 'summary_pcr_metrics.csv'),
  file.path(out_fig, 'syngnatharia_postfit_metric_family_values.png'),
  'Syngnatharia: post-fit evaluation, 3 families',
  short_syn, col_syn
)

## ---- Terapontoidei ----
short_terap_few <- c(
  chronos_discrete = 'Discrete', chronos_clock = 'Clock',
  chronos_correlated = 'Correlated', treePL = 'treePL',
  chronos_relaxed = 'Relaxed', RelTime = 'RelTime',
  RelTime_with_tao_CI = 'RelTime(Tao)', RelTime_with_bootstrap_CI = 'RelTime(Boot)'
)
col_terap_few <- c(
  chronos_discrete = '#2c7fb8', chronos_clock = '#1b9e77',
  chronos_correlated = '#d95f0e', treePL = '#6baed6',
  chronos_relaxed = '#7570b3', RelTime = '#d7301f',
  RelTime_with_tao_CI = '#d7301f', RelTime_with_bootstrap_CI = '#fc8d59'
)
make_figure(
  file.path(base_dir, 'examples', 'terapontoid', 'pcr_rerun_3fam', 'summary_pcr_metrics.csv'),
  file.path(out_fig, 'terapontoid_postfit_FEW.png'),
  'Terapontoidei (6 calibrations): post-fit evaluation, 3 families',
  short_terap_few, col_terap_few
)

short_terap_many <- c(
  RelTime_MEGA = 'RelTime', RelTime_MEGA_with_bootstrap_CI = 'RelTime(Boot)',
  RelTime_MEGA_with_tao_CI = 'RelTime(Tao)',
  chronos_clock_lambda1 = 'Clock',
  chronos_discrete_lambda1_k2 = 'Discrete', treePL_smooth0p01 = 'treePL',
  chronos_correlated_lambda0p01 = 'Correlated', chronos_relaxed_lambda0p01 = 'Relaxed'
)
col_terap_many <- c(
  RelTime_MEGA = '#d7301f', RelTime_MEGA_with_bootstrap_CI = '#fc8d59',
  RelTime_MEGA_with_tao_CI = '#d7301f',
  chronos_clock_lambda1 = '#1b9e77',
  chronos_discrete_lambda1_k2 = '#2c7fb8', treePL_smooth0p01 = '#6baed6',
  chronos_correlated_lambda0p01 = '#d95f0e', chronos_relaxed_lambda0p01 = '#7570b3'
)
make_figure(
  file.path(base_dir, '5_DATASETS_MANY_CALS', 'terapontoidei', 'pcr_rerun_3fam_restored', 'summary_pcr_metrics.csv'),
  file.path(out_fig, 'terapontoid_postfit_MANY.png'),
  'Terapontoidei (17 calibrations): post-fit evaluation, 3 families',
  short_terap_many, col_terap_many
)

## ---- Gobiaria (single scheme, McCraney minus kurtiform) ----
make_figure(
  file.path(base_dir, 'GOBIES', 'GOBY_TREES',
            'McCraney_minus_kurtiform_cals_plus_eosphaeramia_apogoninae',
            'pcr_output', 'summary_pcr_metrics.csv'),
  file.path(out_fig, 'gobiaria_postfit.png'),
  'Gobiaria (75 congruified calibrations): post-fit evaluation, 3 families',
  short_gobi, col_gobi
)

## ---- Other 2 datasets (vertebrate, tetraodontiformes) ----
datasets <- list(
  list(name = 'vertebrate', few_ncal = 18, many_ncal = 57,
       short = short_vert, cols = col_vert),
  list(name = 'tetraodontiformes', few_ncal = 38, many_ncal = 120,
       short = short_tforms, cols = col_tforms)
)

for (ds in datasets) {
  nm <- ds$name
  make_figure(
    file.path(base_dir, '5_DATASETS_FEW_CALS', nm, 'pcr_rerun_3fam', 'summary_pcr_metrics.csv'),
    file.path(out_fig, paste0(nm, '_postfit_FEW.png')),
    paste0(tools::toTitleCase(nm), ' (', ds$few_ncal, ' calibrations): post-fit evaluation, 3 families'),
    ds$short, ds$cols
  )
  make_figure(
    file.path(base_dir, '5_DATASETS_MANY_CALS', nm, 'pcr_rerun_3fam', 'summary_pcr_metrics.csv'),
    file.path(out_fig, paste0(nm, '_postfit_MANY.png')),
    paste0(tools::toTitleCase(nm), ' (', ds$many_ncal, ' calibrations): post-fit evaluation, 3 families'),
    ds$short, ds$cols
  )
}

## ---- Ostariophysi ----
## 5 candidates: two RelTime CI variants (deduped for core metrics), chronos_clock,
## chronos_correlated, treePL.  All CI-embedded trees except chronos_correlated.
short_ostario <- c(
  RelTime_tao_CI = 'RelTime(Tao)', RelTime_bootstrap_CI = 'RelTime(Boot)',
  chronos_clock = 'Clock', treePL = 'treePL'
)
col_ostario <- c(
  RelTime_tao_CI = '#d7301f', RelTime_bootstrap_CI = '#fc8d59',
  chronos_clock = '#1b9e77', treePL = '#6baed6'
)
ostario_csv <- file.path(base_dir, '..', 'Ostario', 'run1_plusAfro_plusCall', 'pcr_output', 'summary_pcr_metrics.csv')
make_figure(
  ostario_csv,
  file.path(out_fig, 'ostariophysi_postfit_v2.png'),
  'Ostariophysi (56 calibrations): post-fit evaluation, 3 families',
  short_ostario, col_ostario
)

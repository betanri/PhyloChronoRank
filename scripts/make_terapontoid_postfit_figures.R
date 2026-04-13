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
plot_panel <- function(vals, labels, cols, ttl, ylab, note, fam_label = '', higher_better = FALSE) {
  ylim_top <- max(vals, na.rm = TRUE) * 1.30
  if (all(is.na(vals) | vals == 0)) ylim_top <- 1
  bp <- barplot(vals, names.arg = labels, col = cols, las = 2,
                ylab = ylab, main = ttl,
                ylim = c(0, ylim_top))
  text(bp, vals, labels = sprintf('%.3f', vals), pos = 3, cex = 0.95)
  orient <- if (higher_better) 'Higher is better' else 'Lower is better'
  mtext(paste0(orient, '   ', fam_label), side = 1, line = 6.2, cex = 0.85)
}

## ============================================================
## PANEL A: FEW calibrations (6 cal, 6 candidates)
## ============================================================
few_dir <- normalizePath(file.path(base_dir, 'examples', 'terapontoid'),
                         winslash = '/', mustWork = TRUE)
few_file <- file.path(few_dir, 'pcr_rerun_3fam', 'summary_pcr_metrics.csv')
if (!file.exists(few_file)) stop('Missing input file: ', few_file)
d_few <- read.csv(few_file, stringsAsFactors = FALSE)
ord_few <- d_few$candidate[order(d_few$rank_mean_core_rank, d_few$rank_mean_core)]
d_few <- d_few[match(ord_few, d_few$candidate), ]

label_map_few <- c(
  chronos_discrete = 'Discrete', chronos_clock = 'Clock',
  chronos_correlated = 'Correlated', treePL = 'treePL',
  chronos_relaxed = 'Relaxed', RelTime = 'RelTime'
)
labels_few <- unname(label_map_few[ord_few])
col_pool_few <- c(chronos_discrete = '#2c7fb8', chronos_clock = '#1b9e77',
                  chronos_correlated = '#d95f0e', treePL = '#6baed6',
                  chronos_relaxed = '#7570b3', RelTime = '#d7301f')
cols_few <- unname(col_pool_few[ord_few])

png(file.path(out_fig, 'terapontoid_postfit_FEW.png'), width = 3400, height = 2000, res = 170)
layout(matrix(c(1,2,3,4, 5,6,7,8, 9,9,10,10), nrow = 3, byrow = TRUE),
       heights = c(1, 1, 1))
par(mar = c(8, 5, 4, 1), oma = c(2.6, 0.2, 2.2, 0.2))

## Row 1
plot_panel(d_few$burst_loss, labels_few, cols_few, 'Burst loss', 'Burst loss', '', 'Family 1: Radiation-zone fidelity')
plot_panel(d_few$internode_concordance, labels_few, cols_few, 'Internode concordance', 'Concordance', '', 'Family 1: Radiation-zone fidelity', higher_better = TRUE)
plot_panel(d_few$compression_score, labels_few, cols_few, 'Compression score', 'Compression', '', 'Family 2: Global fidelity')
plot_panel(d_few$tempo_redistribution, labels_few, cols_few, 'Tempo redistribution', 'EMD', '', 'Family 2: Global fidelity')
## Row 2
plot_panel(d_few$depth_r2, labels_few, cols_few, 'Depth R\u00B2', 'R\u00B2', '', 'Family 2: Global fidelity', higher_better = TRUE)
plot_panel(d_few$rate_irregularity, labels_few, cols_few, 'Rate irregularity', 'Rate irregularity', '', 'Family 3: Rate & Calibration')
# mean_relative_gap is NA for this dataset
vals_gap_few <- d_few$mean_relative_gap
if (all(is.na(vals_gap_few))) {
  plot.new()
  text(0.5, 0.5, 'Mean relative gap\nnot scored\n(secondary calibrations)', cex = 1.2)
} else {
  plot_panel(vals_gap_few, labels_few, cols_few, 'Mean relative gap', 'Relative gap', '', 'Family 3: Rate & Calibration')
}
plot_panel(d_few$uncertainty_mean_width_ma, labels_few, cols_few, 'Uncertainty width', 'Mean CI width (Ma)', '', 'Optional precision layer')
## Row 3
par(mar = c(8, 5, 4, 1))
fam_vals_few <- rbind(d_few$family_radiation_score, d_few$family_global_score, d_few$family_ratecal_score)
colnames(fam_vals_few) <- labels_few
rownames(fam_vals_few) <- c('Fam 1: Radiation', 'Fam 2: Global', 'Fam 3: Rate+Cal')
fam_cols <- c('#e41a1c', '#377eb8', '#4daf4a')
bp <- barplot(fam_vals_few, beside = TRUE, col = fam_cols, las = 2,
              ylab = 'Family rank (lower is better)',
              main = 'Family scores',
              ylim = c(0, max(fam_vals_few, na.rm = TRUE) * 1.35),
              legend.text = rownames(fam_vals_few),
              args.legend = list(x = 'topright', cex = 0.85, bty = 'n'))
text(bp, fam_vals_few, labels = sprintf('%.1f', fam_vals_few), pos = 3, cex = 0.8)
mtext('Lower is better', side = 1, line = 6.2, cex = 0.85)

bp2 <- barplot(d_few$rank_mean_core, names.arg = labels_few, col = cols_few, las = 2,
               ylab = 'Mean rank across 3 families',
               main = 'Core overall rank (3-family balanced)',
               ylim = c(0, max(d_few$rank_mean_core, na.rm = TRUE) * 1.4))
text(bp2, d_few$rank_mean_core,
     labels = paste0(sprintf('%.2f', d_few$rank_mean_core), ' (rank ', d_few$rank_mean_core_rank, ')'),
     pos = 3, cex = 0.9)
mtext('Lower is better', side = 1, line = 6.0, cex = 0.85)

mtext('Terapontoidei (6 calibrations): post-fit evaluation across 7 metrics, 3 families', side = 3, outer = TRUE, line = 0.5, cex = 1.5, font = 2)
mtext('Core PCR rank balances radiation-zone fidelity, global chronogram fidelity, and rate+calibration. mean_relative_gap NA (secondary calibrations). Uncertainty width shown separately.', side = 1, outer = TRUE, line = 0.5, cex = 0.95)
dev.off()
message('Done. Wrote: ', normalizePath(file.path(out_fig, 'terapontoid_postfit_FEW.png'), winslash = '/'))

## ============================================================
## PANEL B: MANY calibrations (17 cal, 8 candidates → deduplicated to 6 unique)
## ============================================================
many_dir <- normalizePath(file.path(base_dir, '5_DATASETS_MANY_CALS', 'terapontoidei'),
                          winslash = '/', mustWork = FALSE)
many_file <- file.path(many_dir, 'pcr_rerun_3fam_restored', 'summary_pcr_metrics.csv')
if (file.exists(many_file)) {
  d_many <- read.csv(many_file, stringsAsFactors = FALSE)
  # Deduplicate: keep only one row per unique core rank (3 RelTime entries are identical)
  d_many <- d_many[!duplicated(d_many[, c('burst_loss', 'internode_concordance', 'rate_irregularity')]), ]
  ord_many <- d_many$candidate[order(d_many$rank_mean_core_rank, d_many$rank_mean_core)]
  d_many <- d_many[match(ord_many, d_many$candidate), ]

  label_map_many <- c(
    RelTime_MEGA = 'RelTime', chronos_clock_lambda1 = 'Clock',
    chronos_discrete_lambda1_k2 = 'Discrete', treePL_smooth0p01 = 'treePL',
    chronos_correlated_lambda0p01 = 'Correlated', chronos_relaxed_lambda0p01 = 'Relaxed'
  )
  labels_many <- unname(label_map_many[ord_many])
  labels_many[is.na(labels_many)] <- ord_many[is.na(labels_many)]
  col_pool_many <- c(RelTime_MEGA = '#d7301f', chronos_clock_lambda1 = '#1b9e77',
                     chronos_discrete_lambda1_k2 = '#2c7fb8', treePL_smooth0p01 = '#6baed6',
                     chronos_correlated_lambda0p01 = '#d95f0e', chronos_relaxed_lambda0p01 = '#7570b3')
  cols_many <- unname(col_pool_many[ord_many])
  cols_many[is.na(cols_many)] <- 'grey60'

  png(file.path(out_fig, 'terapontoid_postfit_MANY.png'), width = 3400, height = 2000, res = 170)
  layout(matrix(c(1,2,3,4, 5,6,7,8, 9,9,10,10), nrow = 3, byrow = TRUE),
         heights = c(1, 1, 1))
  par(mar = c(8, 5, 4, 1), oma = c(2.6, 0.2, 2.2, 0.2))

  ## Row 1
  plot_panel(d_many$burst_loss, labels_many, cols_many, 'Burst loss', 'Burst loss', '', 'Family 1: Radiation-zone fidelity')
  plot_panel(d_many$internode_concordance, labels_many, cols_many, 'Internode concordance', 'Concordance', '', 'Family 1: Radiation-zone fidelity', higher_better = TRUE)
  plot_panel(d_many$compression_score, labels_many, cols_many, 'Compression score', 'Compression', '', 'Family 2: Global fidelity')
  plot_panel(d_many$tempo_redistribution, labels_many, cols_many, 'Tempo redistribution', 'EMD', '', 'Family 2: Global fidelity')
  ## Row 2
  plot_panel(d_many$depth_r2, labels_many, cols_many, 'Depth R\u00B2', 'R\u00B2', '', 'Family 2: Global fidelity', higher_better = TRUE)
  plot_panel(d_many$rate_irregularity, labels_many, cols_many, 'Rate irregularity', 'Rate irregularity', '', 'Family 3: Rate & Calibration')
  vals_gap_many <- d_many$mean_relative_gap
  if (all(is.na(vals_gap_many))) {
    plot.new()
    text(0.5, 0.5, 'Mean relative gap\nnot scored\n(secondary calibrations)', cex = 1.2)
  } else {
    plot_panel(vals_gap_many, labels_many, cols_many, 'Mean relative gap', 'Relative gap', '', 'Family 3: Rate & Calibration')
  }
  uw <- d_many$uncertainty_mean_width_ma
  if (all(is.na(uw))) {
    plot.new()
    text(0.5, 0.5, 'Uncertainty width\nnot available\nfor all candidates', cex = 1.2)
  } else {
    uw[is.na(uw)] <- 0
    plot_panel(uw, labels_many, cols_many, 'Uncertainty width', 'Mean CI width (Ma)', '', 'Optional precision layer')
  }
  ## Row 3
  par(mar = c(8, 5, 4, 1))
  fam_vals_many <- rbind(d_many$family_radiation_score, d_many$family_global_score, d_many$family_ratecal_score)
  colnames(fam_vals_many) <- labels_many
  rownames(fam_vals_many) <- c('Fam 1: Radiation', 'Fam 2: Global', 'Fam 3: Rate+Cal')
  bp <- barplot(fam_vals_many, beside = TRUE, col = fam_cols, las = 2,
                ylab = 'Family rank (lower is better)',
                main = 'Family scores',
                ylim = c(0, max(fam_vals_many, na.rm = TRUE) * 1.35),
                legend.text = rownames(fam_vals_many),
                args.legend = list(x = 'topright', cex = 0.85, bty = 'n'))
  text(bp, fam_vals_many, labels = sprintf('%.1f', fam_vals_many), pos = 3, cex = 0.8)
  mtext('Lower is better', side = 1, line = 6.2, cex = 0.85)

  bp2 <- barplot(d_many$rank_mean_core, names.arg = labels_many, col = cols_many, las = 2,
                 ylab = 'Mean rank across 3 families',
                 main = 'Core overall rank (3-family balanced)',
                 ylim = c(0, max(d_many$rank_mean_core, na.rm = TRUE) * 1.4))
  text(bp2, d_many$rank_mean_core,
       labels = paste0(sprintf('%.2f', d_many$rank_mean_core), ' (rank ', d_many$rank_mean_core_rank, ')'),
       pos = 3, cex = 0.9)
  mtext('Lower is better', side = 1, line = 6.0, cex = 0.85)

  mtext('Terapontoidei (17 calibrations): post-fit evaluation across 7 metrics, 3 families', side = 3, outer = TRUE, line = 0.5, cex = 1.5, font = 2)
  mtext('Core PCR rank balances radiation-zone fidelity, global chronogram fidelity, and rate+calibration. mean_relative_gap NA (secondary calibrations). Uncertainty width shown separately.', side = 1, outer = TRUE, line = 0.5, cex = 0.95)
  dev.off()
  message('Done. Wrote: ', normalizePath(file.path(out_fig, 'terapontoid_postfit_MANY.png'), winslash = '/'))
} else {
  message('Skipping MANY figure: ', many_file, ' not found')
}

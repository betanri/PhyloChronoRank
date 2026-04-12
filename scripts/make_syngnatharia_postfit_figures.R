args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
if (length(file_arg)) {
  script_dir <- dirname(normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE))
} else {
  script_dir <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
}
base_dir <- normalizePath(file.path(script_dir, '..'), winslash = '/', mustWork = TRUE)
out_fig <- file.path(base_dir, 'figures')
example_dir <- normalizePath(file.path(base_dir, 'examples', 'syngnatharia'),
                             winslash = '/', mustWork = TRUE)
infile <- file.path(example_dir, 'pcr_rerun_3fam', 'summary_pcr_metrics.csv')
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(infile)) stop('Missing input file: ', infile)
d <- read.csv(infile, stringsAsFactors = FALSE)
ord <- d$candidate[order(d$rank_mean_core_rank, d$rank_mean_core)]
d <- d[match(ord, d$candidate), ]
labels <- ord
cols <- c(RelTime = '#1b9e77', MCMCTree = '#2c7fb8')
cols <- unname(cols[ord])

png(file.path(out_fig, 'syngnatharia_postfit_metric_family_values.png'), width = 3400, height = 2000, res = 170)
layout(matrix(c(1,2,3,4, 5,6,7,8, 9,9,10,10), nrow = 3, byrow = TRUE),
       heights = c(1, 1, 1))
par(mar = c(8, 5, 4, 1), oma = c(2.6, 0.2, 2.2, 0.2))

plot_panel <- function(vals, ttl, ylab, note, fam_label = '', higher_better = FALSE) {
  ylim_top <- max(vals, na.rm = TRUE) * 1.30
  if (all(is.na(vals) | vals == 0)) ylim_top <- 1
  bp <- barplot(vals, names.arg = labels, col = cols, las = 2,
                ylab = ylab, main = ttl,
                ylim = c(0, ylim_top))
  text(bp, vals, labels = sprintf('%.3f', vals), pos = 3, cex = 0.95)
  orient <- if (higher_better) 'Higher is better' else 'Lower is better'
  mtext(paste0(orient, '   ', fam_label), side = 1, line = 6.2, cex = 0.85)
}

## Row 1: Family 1 (radiation zone) + Family 2 first metric
plot_panel(d$burst_loss, 'Burst loss', 'Burst loss', '', 'Family 1: Radiation-zone fidelity')
plot_panel(d$internode_concordance, 'Internode concordance', 'Concordance', '', 'Family 1: Radiation-zone fidelity', higher_better = TRUE)
plot_panel(d$compression_score, 'Compression score', 'Compression', '', 'Family 2: Global fidelity')
plot_panel(d$tempo_redistribution, 'Tempo redistribution', 'EMD', '', 'Family 2: Global fidelity')

## Row 2: Family 2 (cont.) + Family 3 + Uncertainty
plot_panel(d$depth_r2, 'Depth R\u00B2', 'R\u00B2', '', 'Family 2: Global fidelity', higher_better = TRUE)
plot_panel(d$rate_irregularity, 'Rate irregularity', 'Rate irregularity', '', 'Family 3: Rate & Calibration')
plot_panel(d$mean_relative_gap, 'Mean relative gap', 'Relative gap', '', 'Family 3: Rate & Calibration')
plot_panel(d$uncertainty_mean_width_ma, 'Uncertainty width', 'Mean HPD width (Ma)', '', 'Optional precision layer')

## Row 3: Family scores + overall
par(mar = c(8, 5, 4, 1))
fam_vals <- rbind(d$family_radiation_score, d$family_global_score, d$family_ratecal_score)
colnames(fam_vals) <- labels
rownames(fam_vals) <- c('Fam 1: Radiation', 'Fam 2: Global', 'Fam 3: Rate+Cal')
fam_cols <- c('#e41a1c', '#377eb8', '#4daf4a')
bp <- barplot(fam_vals, beside = TRUE, col = fam_cols, las = 2,
              ylab = 'Family rank (lower is better)',
              main = 'Family scores',
              ylim = c(0, max(fam_vals, na.rm = TRUE) * 1.35),
              legend.text = rownames(fam_vals),
              args.legend = list(x = 'topright', cex = 0.85, bty = 'n'))
text(bp, fam_vals, labels = sprintf('%.1f', fam_vals), pos = 3, cex = 0.8)
mtext('Lower is better', side = 1, line = 6.2, cex = 0.85)

bp2 <- barplot(d$rank_mean_core, names.arg = labels, col = cols, las = 2,
               ylab = 'Mean rank across 3 families',
               main = 'Core overall rank (3-family balanced)',
               ylim = c(0, max(d$rank_mean_core, na.rm = TRUE) * 1.4))
text(bp2, d$rank_mean_core,
     labels = paste0(sprintf('%.2f', d$rank_mean_core), ' (rank ', d$rank_mean_core_rank, ')'),
     pos = 3, cex = 0.9)
mtext('Lower is better', side = 1, line = 6.0, cex = 0.85)

mtext('Syngnatharia: post-fit evaluation across 7 metrics, 3 families', side = 3, outer = TRUE, line = 0.5, cex = 1.5, font = 2)
mtext('Core PCR rank balances radiation-zone fidelity, global chronogram fidelity, and rate+calibration. Uncertainty width shown separately.', side = 1, outer = TRUE, line = 0.5, cex = 0.95)
dev.off()
message('Done. Wrote: ', normalizePath(file.path(out_fig, 'syngnatharia_postfit_metric_family_values.png'), winslash = '/'))

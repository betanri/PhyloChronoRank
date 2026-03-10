suppressPackageStartupMessages({library(grDevices)})
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
if (length(file_arg)) {
  script_dir <- dirname(normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE))
} else {
  script_dir <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
}
base_dir <- normalizePath(file.path(script_dir, '..'), winslash = '/', mustWork = TRUE)
infile <- file.path(base_dir, 'examples', 'unpublished_vertebrate', 'postfit_metrics', 'summary_unpublished_vertebrate_postfit_metrics.csv')
outfile <- file.path(base_dir, 'figures', 'unpublished_vertebrate_postfit_metric_family_values.png')
d <- read.csv(infile, stringsAsFactors = FALSE)
d <- d[order(d$rank_mean_core_rank, d$rank_mean_core), ]
label_map <- c(
  chronos_clock='Chronos clock',
  chronos_correlated='Chronos correlated',
  chronos_relaxed='Chronos relaxed',
  chronos_discrete='Chronos discrete',
  `treepl_best-smooth-100`='treePL (smooth=100)'
)
labels <- unname(label_map[d$candidate])
cols <- c('#1b9e77','#2c7fb8','#7570b3','#d95f0e','#6baed6')[match(d$candidate, c('chronos_clock','chronos_correlated','chronos_relaxed','chronos_discrete','treepl_best-smooth-100'))]
png(outfile, width = 2600, height = 1500, res = 170)
layout(matrix(1:6, nrow = 2, byrow = TRUE))
par(mar = c(10, 5, 4, 1), oma = c(2.5, 0.2, 2.3, 0.2))
plot_panel <- function(vals, ttl, ylab) {
  bp <- barplot(vals, names.arg = labels, col = cols, las = 2, ylab = ylab, main = ttl,
                ylim = c(0, max(vals, na.rm = TRUE) * 1.22), cex.names = 0.9)
  text(bp, vals, labels = sprintf('%.3f', vals), pos = 3, cex = 0.85)
  mtext('Lower is better', side = 1, line = 7.3, cex = 0.9)
}
plot_panel(d$burst_loss, 'Burst loss', 'Burst loss')
plot_panel(d$pulse_burst_selector_error, 'Pulse preservation (burst)', 'Selector error')
plot_panel(d$pulse_default_selector_error, 'Pulse preservation (overall)', 'Selector error')
plot_panel(d$mean_relative_gap, 'Mean relative gap', 'Mean relative gap')
plot_panel(d$rate_irregularity, 'Rate plausibility', 'Rate irregularity')
bp <- barplot(d$rank_mean_core, names.arg = labels, col = cols, las = 2,
              ylab = 'Mean rank across 3 families',
              main = 'Core overall rank (family-balanced; pulse = 1/3)',
              ylim = c(0, max(d$rank_mean_core, na.rm = TRUE) * 1.3), cex.names = 0.9)
text(bp, d$rank_mean_core, labels = paste0(sprintf('%.2f', d$rank_mean_core), ' (rank ', d$rank_mean_core_rank, ')'), pos = 3, cex = 0.82)
mtext('Lower is better', side = 1, line = 7.0, cex = 0.9)
mtext('Unpublished vertebrate example: post-fit evaluation metrics across selected chronograms', side = 3, outer = TRUE, line = 0.5, cex = 1.5, font = 2)
mtext('Core PCR rank is family-balanced: the 3 pulse panels count together as one pulse family (1/3), alongside mean relative gap (1/3) and rate plausibility (1/3).', side = 1, outer = TRUE, line = 0.7, cex = 1.0)
dev.off()
message('Wrote: ', normalizePath(outfile, winslash = '/'))

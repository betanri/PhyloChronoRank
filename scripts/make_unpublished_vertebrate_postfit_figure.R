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
labels <- ifelse(
  grepl('^treepl_best-smooth-', d$candidate),
  paste0('treePL (smooth=', sub('^treepl_best-smooth-', '', d$candidate), ')'),
  unname(c(
    chronos_clock='Chronos clock',
    chronos_correlated='Chronos correlated',
    chronos_relaxed='Chronos relaxed',
    chronos_discrete='Chronos discrete',
    RelTime='RelTime'
  )[d$candidate])
)
color_key <- ifelse(grepl('^treepl_best-smooth-', d$candidate), 'treePL', d$candidate)
cols <- unname(c(
  chronos_clock = '#1b9e77',
  chronos_correlated = '#2c7fb8',
  chronos_relaxed = '#7570b3',
  chronos_discrete = '#d95f0e',
  treePL = '#6baed6',
  RelTime = '#d7301f'
)[color_key])
has_uncertainty <- 'uncertainty_mean_width_ma' %in% names(d) && any(is.finite(d$uncertainty_mean_width_ma))
png(outfile, width = 3000, height = 1500, res = 170)
if (has_uncertainty) {
  layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 7), nrow = 2, byrow = TRUE))
} else {
  layout(matrix(1:6, nrow = 2, byrow = TRUE))
}
par(mar = c(10, 5, 4, 1), oma = c(2.5, 0.2, 2.3, 0.2))
plot_panel <- function(vals, ttl, ylab) {
  vals_num <- as.numeric(vals)
  vals_plot <- vals_num
  vals_plot[!is.finite(vals_plot)] <- 0
  ylim_max <- max(vals_num[is.finite(vals_num)], na.rm = TRUE)
  if (!is.finite(ylim_max) || ylim_max <= 0) ylim_max <- 1
  bp <- barplot(vals_plot, names.arg = labels, col = cols, las = 2, ylab = ylab, main = ttl,
                ylim = c(0, ylim_max * 1.22), cex.names = 0.88)
  ok <- is.finite(vals_num)
  if (any(ok)) {
    text(bp[ok], vals_num[ok], labels = sprintf('%.3f', vals_num[ok]), pos = 3, cex = 0.85)
  }
  if (any(!ok)) {
    text(bp[!ok], 0, labels = 'not scored', pos = 3, cex = 0.78, xpd = TRUE)
  }
  mtext('Lower is better', side = 1, line = 7.3, cex = 0.9)
}
plot_panel(d$burst_loss, 'Burst loss', 'Burst loss')
plot_panel(d$pulse_burst_selector_error, 'Pulse preservation (burst)', 'Selector error')
plot_panel(d$pulse_default_selector_error, 'Pulse preservation (overall)', 'Selector error')
plot_panel(d$mean_relative_gap, 'Mean relative gap', 'Mean relative gap')
plot_panel(d$rate_irregularity, 'Rate irregularity', 'Rate irregularity')
if (has_uncertainty) {
  plot_panel(d$uncertainty_mean_width_ma, 'Uncertainty width', 'Mean CI width (Ma)')
}
bp <- barplot(d$rank_mean_core, names.arg = labels, col = cols, las = 2,
              ylab = 'Mean rank across 3 families',
              main = 'Core overall rank (family-balanced; pulse = 1/3)',
              ylim = c(0, max(d$rank_mean_core, na.rm = TRUE) * 1.3), cex.names = 0.9)
text(bp, d$rank_mean_core, labels = paste0(sprintf('%.2f', d$rank_mean_core), ' (rank ', d$rank_mean_core_rank, ')'), pos = 3, cex = 0.82)
mtext('Lower is better', side = 1, line = 7.0, cex = 0.9)
mtext('Unpublished vertebrate example: post-fit evaluation metrics across selected chronograms', side = 3, outer = TRUE, line = 0.5, cex = 1.5, font = 2)
mtext(
  if (has_uncertainty) {
    'Core PCR rank is family-balanced: the 3 pulse panels count together as one pulse family (1/3), alongside mean relative gap (1/3) and rate irregularity (1/3). Uncertainty width is shown separately as an optional precision layer; blank bars mean not scored.'
  } else {
    'Core PCR rank is family-balanced: the 3 pulse panels count together as one pulse family (1/3), alongside mean relative gap (1/3) and rate irregularity (1/3).'
  },
  side = 1, outer = TRUE, line = 0.7, cex = 1.0
)
dev.off()
message('Wrote: ', normalizePath(outfile, winslash = '/'))

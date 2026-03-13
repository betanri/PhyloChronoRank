args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
if (length(file_arg)) {
  script_dir <- dirname(normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE))
} else {
  script_dir <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
}
base_dir <- normalizePath(file.path(script_dir, '..'), winslash = '/', mustWork = TRUE)
out_fig <- file.path(base_dir, 'figures')
pipeline_demo_dir <- normalizePath(file.path(base_dir, 'examples', 'terapontoid'),
                                   winslash = '/', mustWork = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

infile <- file.path(pipeline_demo_dir, 'summary_terap_empirical_postfit_metrics.csv')
if (!file.exists(infile)) stop('Missing input file: ', infile)
d <- read.csv(infile, stringsAsFactors = FALSE)
if (!('rank_burst_loss' %in% names(d))) d$rank_burst_loss <- rank(d$burst_loss, ties.method = 'min')
if (!('rank_pulse_family_mean' %in% names(d))) d$rank_pulse_family_mean <- rowMeans(cbind(d$rank_pulse_overall, d$rank_burst_loss, d$rank_pulse_burst))
if (!('rank_pulse_family' %in% names(d))) d$rank_pulse_family <- rank(d$rank_pulse_family_mean, ties.method = 'min')
if (!('rank_mean_core' %in% names(d))) d$rank_mean_core <- rowMeans(cbind(d$rank_pulse_family, d$rank_rate_irregularity), na.rm = TRUE)
if (!('rank_mean_core_rank' %in% names(d))) d$rank_mean_core_rank <- rank(d$rank_mean_core, ties.method = 'min')
ord <- d$candidate[order(d$rank_mean_core_rank, d$rank_mean_core, d$rank_pulse_family, d$rank_rate_irregularity)]
d <- d[match(ord, d$candidate), ]
label_map <- c(
  chronos_discrete = 'Discrete',
  chronos_clock = 'Clock',
  chronos_correlated = 'Correlated',
  treePL = 'treePL',
  chronos_relaxed = 'Relaxed',
  RelTime = 'RelTime'
)
labels <- unname(label_map[ord])
cols <- c('#1b9e77', '#2c7fb8', '#d95f0e', '#6baed6', '#7570b3', '#d7301f')
cols <- cols[match(ord, c('chronos_discrete', 'chronos_clock', 'chronos_correlated', 'treePL', 'chronos_relaxed', 'RelTime'))]

png(file.path(out_fig, 'postfit_metric_family_values.png'), width = 2800, height = 1500, res = 170)
layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 2, byrow = TRUE))
par(mar = c(8, 5, 4, 1), oma = c(2.4, 0.2, 2.2, 0.2))
plot_panel <- function(vals, ttl, ylab, note) {
  bp <- barplot(vals, names.arg = labels, col = cols, las = 2,
                ylab = ylab, main = ttl,
                ylim = c(0, max(vals, na.rm = TRUE) * 1.22), cex.names = 0.95)
  text(bp, vals, labels = sprintf('%.3f', vals), pos = 3, cex = 0.85)
  mtext(note, side = 1, line = 6.4, cex = 0.9)
}
plot_panel(d$burst_loss, 'Burst loss', 'Burst loss', 'Lower is better')
plot_panel(d$pulse_burst_selector_error, 'Pulse preservation (burst)', 'Selector error', 'Lower is better')
plot_panel(d$pulse_default_selector_error, 'Pulse preservation (overall)', 'Selector error', 'Lower is better')
plot_panel(d$rate_irregularity, 'Rate irregularity', 'Rate irregularity', 'Lower is better')

par(mar = c(8, 5, 4, 1))
bp <- barplot(d$rank_mean_core, names.arg = labels, col = cols, las = 2,
              ylab = 'Mean rank across 2 families',
              main = 'Core overall rank (family-balanced; pulse = 1/2)',
              ylim = c(0, max(d$rank_mean_core, na.rm = TRUE) * 1.25), cex.names = 0.95)
text(bp, d$rank_mean_core,
     labels = paste0(sprintf('%.2f', d$rank_mean_core), ' (rank ', d$rank_mean_core_rank, ')'),
     pos = 3, cex = 0.85)
mtext('Lower is better', side = 1, line = 6.0, cex = 0.9)
mtext('Terapontoid example: post-fit evaluation metrics across chronos models, treePL, and RelTime', side = 3, outer = TRUE, line = 0.5, cex = 1.55, font = 2)
mtext('Core PCR rank is family-balanced here: the 3 pulse panels count together as one pulse family (1/2), alongside rate irregularity (1/2). Gap burden is not computed because these trees were dated with congruified / secondary calibrations.', side = 1, outer = TRUE, line = 0.5, cex = 0.95)
dev.off()
message('Done. Wrote: ', normalizePath(file.path(out_fig, 'postfit_metric_family_values.png'), winslash = '/'))

args <- commandArgs(trailingOnly = TRUE)
kv <- list()
for (a in args) {
  if (!grepl('^--', a)) next
  parts <- strsplit(sub('^--', '', a), '=', fixed = TRUE)[[1]]
  key <- parts[1]
  val <- if (length(parts) > 1) paste(parts[-1], collapse = '=') else 'TRUE'
  kv[[key]] <- val
}
need <- c('ref-tree', 'candidates-csv', 'outdir')
miss <- need[!need %in% names(kv)]
if (length(miss)) {
  cat('Usage: Rscript scripts/run_pcr.R --ref-tree=REF.tre --candidates-csv=CANDIDATES.csv --outdir=OUT [--calibrations-csv=CALS.csv] [--uncertainty-csv=UNC.csv]\n', file = stderr())
  stop('Missing required args: ', paste(miss, collapse = ', '))
}
if (!requireNamespace('ape', quietly = TRUE)) stop('Package ape is required.')
script_dir <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = '/', mustWork = FALSE))
if (!nzchar(script_dir) || !dir.exists(script_dir)) script_dir <- getwd()
source(file.path(script_dir, 'pcr_helpers.R'))

candidates_csv <- normalizePath(kv[['candidates-csv']], winslash = '/', mustWork = TRUE)
candidates_base <- dirname(candidates_csv)
outdir <- kv[['outdir']]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir, winslash = '/', mustWork = TRUE)
ref_tree <- ape::read.tree(normalizePath(kv[['ref-tree']], winslash = '/', mustWork = TRUE))
min_tips <- if ('min-tips' %in% names(kv)) as.integer(kv[['min-tips']]) else 8L
min_events <- if ('min-events' %in% names(kv)) as.integer(kv[['min-events']]) else 4L
panel <- pcr_build_pulse_panel(ref_tree, min_tips = min_tips, min_events = min_events)
if (!length(panel)) stop('Could not build pulse panel from reference tree.')

cand <- read.csv(candidates_csv, stringsAsFactors = FALSE)
if (!all(c('candidate', 'tree_file') %in% names(cand))) stop('Candidates CSV must contain candidate,tree_file')
cand$tree_path <- vapply(cand$tree_file, pcr_resolve_path, FUN.VALUE = character(1), base_dir = candidates_base)

calibrations <- NULL
if ('calibrations-csv' %in% names(kv)) {
  calibrations <- read.csv(normalizePath(kv[['calibrations-csv']], winslash = '/', mustWork = TRUE), stringsAsFactors = FALSE, na.strings = c('', 'NA', 'Inf'))
  req <- c('taxonA','taxonB','age_min')
  if (!all(req %in% names(calibrations))) stop('Calibrations CSV must contain ', paste(req, collapse = ','))
  if (!('age_max' %in% names(calibrations))) calibrations$age_max <- NA_real_
}

uncertainty <- NULL
if ('uncertainty-csv' %in% names(kv)) {
  uncertainty <- read.csv(normalizePath(kv[['uncertainty-csv']], winslash = '/', mustWork = TRUE), stringsAsFactors = FALSE)
  if (!('candidate' %in% names(uncertainty))) stop('Uncertainty CSV must contain candidate')
}

rows <- list(); pulse_default_detail <- list(); pulse_burst_detail <- list(); gap_detail <- list(); rate_detail <- list()
for (i in seq_len(nrow(cand))) {
  tr <- ape::read.tree(cand$tree_path[i])
  p_overall <- pcr_score_pulse_panel(ref_tree, tr, panel, w_emd = 0.35, w_burst_loss = 0.55, w_centroid = 0.10)
  p_burst <- pcr_score_pulse_panel(ref_tree, tr, panel, w_emd = 0.20, w_burst_loss = 0.75, w_centroid = 0.05)
  rp <- pcr_rate_metrics(ref_tree, tr)
  row <- data.frame(candidate = cand$candidate[i], tree_file = cand$tree_file[i], stringsAsFactors = FALSE)
  row$pulse_default_selector_error <- p_overall$summary$selector_error
  row$burst_loss <- p_overall$summary$mean_burst_loss
  row$pulse_burst_selector_error <- p_burst$summary$selector_error
  row$rate_irregularity <- rp$summary$rate_irregularity
  row$gap_mode <- NA_character_
  row$mean_relative_gap <- NA_real_
  row$ghost_mean_ma <- NA_real_
  row$ghost_median_ma <- NA_real_
  row$min_violation_count <- NA_real_
  row$min_violation_sum_ma <- NA_real_
  row$max_violation_count <- NA_real_
  row$max_violation_sum_ma <- NA_real_
  if (!is.null(calibrations)) {
    cal_i <- pcr_subset_calibrations(calibrations, cand$candidate[i])
    if (nrow(cal_i)) {
      gp <- pcr_gap_metrics(tr, cal_i)
      row$gap_mode <- gp$summary$gap_mode
      row$mean_relative_gap <- gp$summary$mean_relative_gap
      row$ghost_mean_ma <- gp$summary$ghost_mean_ma
      row$ghost_median_ma <- gp$summary$ghost_median_ma
      row$min_violation_count <- gp$summary$min_violation_count
      row$min_violation_sum_ma <- gp$summary$min_violation_sum_ma
      row$max_violation_count <- gp$summary$max_violation_count
      row$max_violation_sum_ma <- gp$summary$max_violation_sum_ma
      gd <- gp$detail; gd$candidate <- cand$candidate[i]; gap_detail[[length(gap_detail)+1L]] <- gd
    }
  }
  if (!is.null(uncertainty)) {
    u <- uncertainty[uncertainty$candidate == cand$candidate[i], , drop = FALSE]
    if (nrow(u)) {
      for (nm in setdiff(names(u), 'candidate')) row[[nm]] <- u[[nm]][1]
    }
  }
  rows[[length(rows)+1L]] <- row
  pd <- p_overall$detail; pd$candidate <- cand$candidate[i]; pulse_default_detail[[length(pulse_default_detail)+1L]] <- pd
  pb <- p_burst$detail; pb$candidate <- cand$candidate[i]; pulse_burst_detail[[length(pulse_burst_detail)+1L]] <- pb
  rd <- rp$detail; rd$candidate <- cand$candidate[i]; rate_detail[[length(rate_detail)+1L]] <- rd
}
summary_df <- do.call(rbind, rows)
summary_df$rank_pulse_overall <- pcr_rank_low(summary_df$pulse_default_selector_error)
summary_df$rank_burst_loss <- pcr_rank_low(summary_df$burst_loss)
summary_df$rank_pulse_burst <- pcr_rank_low(summary_df$pulse_burst_selector_error)
summary_df$rank_pulse_family_mean <- rowMeans(cbind(summary_df$rank_burst_loss, summary_df$rank_pulse_burst, summary_df$rank_pulse_overall), na.rm = TRUE)
summary_df$rank_pulse_family <- pcr_rank_low(summary_df$rank_pulse_family_mean)
summary_df$rank_mean_relative_gap <- pcr_rank_low(summary_df$mean_relative_gap)
summary_df$rank_rate_irregularity <- pcr_rank_low(summary_df$rate_irregularity)
core_mat <- cbind(summary_df$rank_pulse_family, summary_df$rank_mean_relative_gap, summary_df$rank_rate_irregularity)
summary_df$core_n_families_ranked <- rowSums(is.finite(core_mat))
summary_df$rank_mean_core <- rowMeans(core_mat, na.rm = TRUE)
summary_df$rank_mean_core_rank <- pcr_rank_low(summary_df$rank_mean_core)
if ('uncertainty_mean_width_ma' %in% names(summary_df)) {
  summary_df$rank_uncertainty_mean_width <- pcr_rank_low(summary_df$uncertainty_mean_width_ma)
  ext_mat <- cbind(summary_df$rank_pulse_family, summary_df$rank_mean_relative_gap, summary_df$rank_rate_irregularity, summary_df$rank_uncertainty_mean_width)
  summary_df$extended_n_families_ranked <- rowSums(is.finite(ext_mat))
  summary_df$rank_mean_extended <- rowMeans(ext_mat, na.rm = TRUE)
  summary_df$rank_mean_extended_rank <- pcr_rank_low(summary_df$rank_mean_extended)
}
write.csv(summary_df, file.path(outdir, 'summary_pcr_metrics.csv'), row.names = FALSE)
if (length(pulse_default_detail)) write.csv(do.call(rbind, pulse_default_detail), file.path(outdir, 'pulse_default_details.csv'), row.names = FALSE)
if (length(pulse_burst_detail)) write.csv(do.call(rbind, pulse_burst_detail), file.path(outdir, 'pulse_burst_details.csv'), row.names = FALSE)
if (length(gap_detail)) write.csv(do.call(rbind, gap_detail), file.path(outdir, 'gap_details.csv'), row.names = FALSE)
if (length(rate_detail)) write.csv(do.call(rbind, rate_detail), file.path(outdir, 'rate_details.csv'), row.names = FALSE)
lines <- c(
  'PhyloChronoRank summary',
  paste0('Reference tree: ', normalizePath(kv[['ref-tree']], winslash = '/', mustWork = TRUE)),
  paste0('Candidates: ', nrow(summary_df)),
  paste0('Pulse panel clades: ', length(panel), ' (min_tips=', min_tips, ', min_events=', min_events, ')'),
  '',
  paste0('Best pulse preservation (overall): ', summary_df$candidate[which.min(summary_df$pulse_default_selector_error)]),
  paste0('Best burst loss: ', summary_df$candidate[which.min(summary_df$burst_loss)]),
  paste0('Best pulse preservation (burst): ', summary_df$candidate[which.min(summary_df$pulse_burst_selector_error)]),
  if (all(!is.na(summary_df$mean_relative_gap))) paste0('Best mean relative gap: ', summary_df$candidate[which.min(summary_df$mean_relative_gap)]) else 'Mean relative gap: not scored',
  paste0('Best rate plausibility: ', summary_df$candidate[which.min(summary_df$rate_irregularity)]),
  paste0('Core overall winner: ', summary_df$candidate[which.min(summary_df$rank_mean_core)]),
  if ('rank_mean_extended' %in% names(summary_df)) paste0('Extended overall winner: ', paste(summary_df$candidate[summary_df$rank_mean_extended == min(summary_df$rank_mean_extended, na.rm=TRUE)], collapse = ', ')) else 'Extended overall winner: not scored'
)
writeLines(lines, file.path(outdir, 'interpretation.txt'))
message('Wrote PCR outputs to: ', outdir)

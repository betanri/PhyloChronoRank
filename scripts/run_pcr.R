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
file_arg <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = FALSE)) else NA_character_
if (is.na(script_dir) || !nzchar(script_dir) || !dir.exists(script_dir)) script_dir <- getwd()
source(file.path(script_dir, 'pcr_helpers.R'))

pcr_read_tree_or_stop <- function(path, label) {
  tr <- try(ape::read.tree(path), silent = TRUE)
  if (inherits(tr, 'try-error') || !inherits(tr, 'phylo')) {
    stop('Failed to read ', label, ' tree: ', path)
  }
  tr
}

pcr_get_metric <- function(df, col) {
  if (is.null(df) || !nrow(df) || !(col %in% names(df))) return(NA_real_)
  df[[col]][1]
}

pcr_best_candidate_label <- function(values, candidates, prefix, fallback) {
  ok <- is.finite(values)
  if (!any(ok)) return(fallback)
  idx <- which.min(values[ok])[1]
  paste0(prefix, candidates[ok][idx])
}

candidates_csv <- normalizePath(kv[['candidates-csv']], winslash = '/', mustWork = TRUE)
candidates_base <- dirname(candidates_csv)
outdir <- kv[['outdir']]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir, winslash = '/', mustWork = TRUE)
ref_tree <- pcr_read_tree_or_stop(normalizePath(kv[['ref-tree']], winslash = '/', mustWork = TRUE), 'reference')
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
  tr <- pcr_read_tree_or_stop(cand$tree_path[i], paste0('candidate ', cand$candidate[i]))
  p_overall <- pcr_score_pulse_panel(ref_tree, tr, panel, w_emd = 0.35, w_burst_loss = 0.55, w_centroid = 0.10)
  p_burst <- pcr_score_pulse_panel(ref_tree, tr, panel, w_emd = 0.20, w_burst_loss = 0.75, w_centroid = 0.05)
  rp <- pcr_rate_metrics(ref_tree, tr)
  row <- data.frame(candidate = cand$candidate[i], tree_file = cand$tree_file[i], stringsAsFactors = FALSE)
  row$pulse_default_selector_error <- pcr_get_metric(p_overall$summary, 'selector_error')
  row$burst_loss <- pcr_get_metric(p_overall$summary, 'mean_burst_loss')
  row$pulse_burst_selector_error <- pcr_get_metric(p_burst$summary, 'selector_error')
  row$rate_irregularity <- pcr_get_metric(rp$summary, 'rate_irregularity')
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
  } else {
    u_tree <- pcr_extract_uncertainty_from_newick(cand$tree_path[i])
    if (!is.null(u_tree)) {
      for (nm in names(u_tree)) row[[nm]] <- u_tree[[nm]][1]
    }
  }
  rows[[length(rows)+1L]] <- row
  if (is.data.frame(p_overall$detail) && nrow(p_overall$detail)) {
    pd <- p_overall$detail
    pd$candidate <- cand$candidate[i]
    pulse_default_detail[[length(pulse_default_detail)+1L]] <- pd
  }
  if (is.data.frame(p_burst$detail) && nrow(p_burst$detail)) {
    pb <- p_burst$detail
    pb$candidate <- cand$candidate[i]
    pulse_burst_detail[[length(pulse_burst_detail)+1L]] <- pb
  }
  if (is.data.frame(rp$detail) && nrow(rp$detail)) {
    rd <- rp$detail
    rd$candidate <- cand$candidate[i]
    rate_detail[[length(rate_detail)+1L]] <- rd
  }
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
  pcr_best_candidate_label(summary_df$pulse_default_selector_error, summary_df$candidate, 'Best pulse preservation (overall): ', 'Pulse preservation (overall): not scored'),
  pcr_best_candidate_label(summary_df$burst_loss, summary_df$candidate, 'Best burst loss: ', 'Burst loss: not scored'),
  pcr_best_candidate_label(summary_df$pulse_burst_selector_error, summary_df$candidate, 'Best pulse preservation (burst): ', 'Pulse preservation (burst): not scored'),
  pcr_best_candidate_label(summary_df$mean_relative_gap, summary_df$candidate, 'Best mean relative gap: ', 'Mean relative gap: not scored'),
  pcr_best_candidate_label(summary_df$rate_irregularity, summary_df$candidate, 'Best rate irregularity: ', 'Rate irregularity: not scored'),
  if ('rank_uncertainty_mean_width' %in% names(summary_df)) pcr_best_candidate_label(summary_df$uncertainty_mean_width_ma, summary_df$candidate, 'Most precise by uncertainty width: ', 'Uncertainty width: not scored') else 'Uncertainty width: not scored',
  pcr_best_candidate_label(summary_df$rank_mean_core, summary_df$candidate, 'Core overall winner: ', 'Core overall winner: not scored')
)
writeLines(lines, file.path(outdir, 'interpretation.txt'))
message('Wrote PCR outputs to: ', outdir)

suppressPackageStartupMessages({
  if (!requireNamespace('ape', quietly = TRUE)) stop('Package ape is required.')
})
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE)) else getwd()
base_dir <- normalizePath(file.path(script_dir, '..'), winslash = '/', mustWork = TRUE)
source(file.path(script_dir, 'pcr_helpers.R'))

weight_sets <- data.frame(
  set = c('default','emd_up','burst_up','centroid_up','centroid_zero'),
  overall_emd = c(0.35,0.45,0.25,0.30,0.40),
  overall_burst = c(0.55,0.45,0.65,0.50,0.60),
  overall_centroid = c(0.10,0.10,0.10,0.20,0.00),
  burst_emd = c(0.20,0.30,0.10,0.15,0.25),
  burst_burst = c(0.75,0.65,0.85,0.65,0.75),
  burst_centroid = c(0.05,0.05,0.05,0.20,0.00),
  stringsAsFactors = FALSE
)

run_example <- function(name, ref_tree_file, candidates_csv, calibrations_csv = NULL) {
  ref <- ape::read.tree(ref_tree_file)
  cand <- read.csv(candidates_csv, stringsAsFactors = FALSE)
  cand_base <- dirname(candidates_csv)
  cand$tree_path <- vapply(cand$tree_file, pcr_resolve_path, FUN.VALUE = character(1), base_dir = cand_base)
  panel <- pcr_build_pulse_panel(ref, min_tips = 8L, min_events = 4L)
  cal <- NULL
  if (!is.null(calibrations_csv) && file.exists(calibrations_csv)) {
    cal <- read.csv(calibrations_csv, stringsAsFactors = FALSE, na.strings = c('', 'NA', 'Inf'))
    if (!('age_max' %in% names(cal))) cal$age_max <- NA_real_
  }
  rows <- list()
  for (i in seq_len(nrow(weight_sets))) {
    ws <- weight_sets[i, ]
    out <- list()
    for (j in seq_len(nrow(cand))) {
      tr <- ape::read.tree(cand$tree_path[j])
      p1 <- pcr_score_pulse_panel(ref, tr, panel, w_emd = ws$overall_emd, w_burst_loss = ws$overall_burst, w_centroid = ws$overall_centroid)
      p2 <- pcr_score_pulse_panel(ref, tr, panel, w_emd = ws$burst_emd, w_burst_loss = ws$burst_burst, w_centroid = ws$burst_centroid)
      rp <- pcr_rate_metrics(ref, tr)
      gap <- NA_real_
      if (!is.null(cal)) {
        gp <- pcr_gap_metrics(tr, pcr_subset_calibrations(cal, cand$candidate[j]))
        gap <- gp$summary$mean_relative_gap
      }
      out[[j]] <- data.frame(candidate = cand$candidate[j], pulse_overall = p1$summary$selector_error, burst_loss = p1$summary$mean_burst_loss, pulse_burst = p2$summary$selector_error, mean_relative_gap = gap, rate_irregularity = rp$summary$rate_irregularity, stringsAsFactors = FALSE)
    }
    d <- do.call(rbind, out)
    d$rank_pulse_family_mean <- rowMeans(cbind(pcr_rank_low(d$burst_loss), pcr_rank_low(d$pulse_burst), pcr_rank_low(d$pulse_overall)), na.rm = TRUE)
    d$rank_pulse_family <- pcr_rank_low(d$rank_pulse_family_mean)
    d$rank_gap <- pcr_rank_low(d$mean_relative_gap)
    d$rank_rate <- pcr_rank_low(d$rate_irregularity)
    d$rank_mean_core <- rowMeans(cbind(d$rank_pulse_family, d$rank_gap, d$rank_rate), na.rm = TRUE)
    rows[[i]] <- data.frame(example = name, set = ws$set,
                            pulse_family_winner = d$candidate[which.min(d$rank_pulse_family_mean)],
                            core_winner = d$candidate[which.min(d$rank_mean_core)],
                            stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

terap <- run_example(
  'terapontoid',
  file.path(base_dir, 'examples', 'terapontoid', 'Terapontoid_ML_MAIN_phylogram_used.tree'),
  file.path(base_dir, 'examples', 'terapontoid', 'candidates.csv'),
  file.path(base_dir, 'examples', 'terapontoid', 'Terapontoid_ML_MAIN_calibrations_used.csv')
)
syng <- run_example(
  'syngnatharia',
  file.path(base_dir, 'examples', 'syngnatharia', 'backbone_Raxml_besttree_matrix75.tre'),
  file.path(base_dir, 'examples', 'syngnatharia', 'candidates.csv'),
  file.path(base_dir, 'examples', 'syngnatharia', 'calibrations_by_candidate.csv')
)
all_res <- rbind(terap, syng)
out_dir <- file.path(base_dir, 'examples', 'weight_sensitivity')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(all_res, file.path(out_dir, 'pulse_weight_sensitivity_sets.csv'), row.names = FALSE)
summary_df <- do.call(rbind, lapply(split(all_res, all_res$example), function(d) {
  data.frame(
    example = d$example[1],
    n_sets = nrow(d),
    default_pulse_family_winner = d$pulse_family_winner[d$set == 'default'][1],
    default_core_winner = d$core_winner[d$set == 'default'][1],
    pulse_family_winner_stable = length(unique(d$pulse_family_winner)) == 1,
    core_winner_stable = length(unique(d$core_winner)) == 1,
    stringsAsFactors = FALSE
  )
}))
write.csv(summary_df, file.path(out_dir, 'pulse_weight_sensitivity_summary.csv'), row.names = FALSE)
message('Wrote weight sensitivity outputs to: ', out_dir)

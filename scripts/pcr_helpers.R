pcr_rank_low <- function(x) {
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  out[ok] <- rank(x[ok], ties.method = 'min')
  out
}

pcr_safe_get_mrca <- function(tr, tips) {
  if (!all(tips %in% tr$tip.label)) return(NA_integer_)
  node <- ape::getMRCA(tr, tips)
  if (is.null(node) || !length(node) || !is.finite(node)) return(NA_integer_)
  as.integer(node)
}

pcr_mean_if_any <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

pcr_median_if_any <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

pcr_metric_node_heights_norm <- function(tr) {
  n <- ape::Ntip(tr)
  d <- ape::node.depth.edgelength(tr)
  maxd <- max(d[seq_len(n)])
  if (!is.finite(maxd) || maxd <= 0) {
    h <- rep(0, length(d))
  } else {
    h <- 1 - (d / maxd)
  }
  names(h) <- as.character(seq_along(h))
  h
}

pcr_metric_node_ages <- function(tr) {
  n <- ape::Ntip(tr)
  d <- ape::node.depth.edgelength(tr)
  maxd <- max(d[seq_len(n)])
  age <- maxd - d
  names(age) <- as.character(seq_along(age))
  age
}

pcr_metric_wasserstein_1d <- function(x, y) {
  x <- sort(as.numeric(x))
  y <- sort(as.numeric(y))
  if (!length(x) || !length(y)) return(NA_real_)
  m <- max(length(x), length(y), 32L)
  p <- seq(0, 1, length.out = m)
  qx <- as.numeric(stats::quantile(x, p, type = 8, names = FALSE))
  qy <- as.numeric(stats::quantile(y, p, type = 8, names = FALSE))
  mean(abs(qx - qy))
}

pcr_metric_event_times_relative <- function(tr) {
  n <- ape::Ntip(tr)
  if (tr$Nnode < 2L) return(numeric(0))
  d <- ape::node.depth.edgelength(tr)
  max_tip <- max(d[seq_len(n)])
  if (!is.finite(max_tip) || max_tip <= 0) return(numeric(0))
  h <- 1 - (d / max_tip)
  hi <- h[(n + 1L):(n + tr$Nnode)]
  crown <- max(hi, na.rm = TRUE)
  if (!is.finite(crown) || crown <= 0) return(numeric(0))
  ev <- hi[hi < (crown - 1e-12)] / crown
  ev <- ev[is.finite(ev)]
  sort(pmax(0, pmin(1, ev)))
}

pcr_metric_burstiness_from_events <- function(ev) {
  if (!length(ev)) return(NA_real_)
  waits <- diff(c(0, ev, 1))
  if (length(waits) < 2L) return(NA_real_)
  m <- mean(waits)
  if (!is.finite(m) || m <= 0) return(NA_real_)
  stats::sd(waits) / m
}

pcr_build_pulse_panel <- function(ref, min_tips = 8L, min_events = 4L) {
  ref_nodes <- (ape::Ntip(ref) + 1L):(ape::Ntip(ref) + ref$Nnode)
  panel <- list()
  for (nd in ref_nodes) {
    sub <- try(ape::extract.clade(ref, nd), silent = TRUE)
    if (inherits(sub, 'try-error') || !inherits(sub, 'phylo')) next
    nt <- ape::Ntip(sub)
    if (nt < min_tips || nt >= ape::Ntip(ref)) next
    ev <- pcr_metric_event_times_relative(sub)
    if (length(ev) < min_events) next
    b <- pcr_metric_burstiness_from_events(ev)
    if (!is.finite(b)) next
    panel[[length(panel) + 1L]] <- list(
      ref_node = nd,
      tips = sort(sub$tip.label),
      n_tips = nt,
      n_events = length(ev),
      ev_ref = ev,
      burst_ref = b,
      centroid_ref = mean(ev)
    )
  }
  panel
}

pcr_empty_pulse_detail <- function() {
  data.frame(
    ref_node = integer(0),
    n_tips = integer(0),
    n_events = integer(0),
    emd = numeric(0),
    burst_ref = numeric(0),
    burst_est = numeric(0),
    burst_loss = numeric(0),
    centroid_shift = numeric(0),
    weight = numeric(0),
    stringsAsFactors = FALSE
  )
}

pcr_empty_pulse_summary <- function(panel_size, global_emd = NA_real_,
                                    global_burst_loss = NA_real_,
                                    global_error = NA_real_) {
  data.frame(
    matched_clades = 0L,
    panel_clades = panel_size,
    coverage = 0,
    mean_emd = NA_real_,
    mean_burst_loss = NA_real_,
    mean_centroid_shift = NA_real_,
    local_error = NA_real_,
    global_emd = global_emd,
    global_burst_loss = global_burst_loss,
    global_error = global_error,
    selector_error = NA_real_,
    selector_score = NA_real_,
    stringsAsFactors = FALSE
  )
}

pcr_score_pulse_panel <- function(ref, tr, panel,
                                  w_emd = 0.35,
                                  w_burst_loss = 0.55,
                                  w_centroid = 0.10,
                                  w_global = 0.20,
                                  coverage_penalty = 0.20) {
  global_ref <- pcr_metric_event_times_relative(ref)
  global_burst_ref <- pcr_metric_burstiness_from_events(global_ref)
  w_sum <- 0
  emd_sum <- 0
  burst_loss_sum <- 0
  centroid_sum <- 0
  matched <- 0L
  detail_rows <- list()
  for (k in seq_along(panel)) {
    cl <- panel[[k]]
    node_c <- pcr_safe_get_mrca(tr, cl$tips)
    if (!is.finite(node_c)) next
    sub_c <- try(ape::extract.clade(tr, node_c), silent = TRUE)
    if (inherits(sub_c, 'try-error') || !inherits(sub_c, 'phylo')) next
    if (!setequal(sub_c$tip.label, cl$tips)) next
    ev_c <- pcr_metric_event_times_relative(sub_c)
    if (length(ev_c) < 2L) next
    emd <- pcr_metric_wasserstein_1d(cl$ev_ref, ev_c)
    burst_c <- pcr_metric_burstiness_from_events(ev_c)
    if (!is.finite(emd) || !is.finite(burst_c)) next
    burst_loss <- max(0, (cl$burst_ref - burst_c) / (cl$burst_ref + 1e-12))
    centroid_shift <- abs(cl$centroid_ref - mean(ev_c))
    w <- log1p(cl$n_tips) * sqrt(cl$n_events)
    w_sum <- w_sum + w
    emd_sum <- emd_sum + (w * emd)
    burst_loss_sum <- burst_loss_sum + (w * burst_loss)
    centroid_sum <- centroid_sum + (w * centroid_shift)
    matched <- matched + 1L
    detail_rows[[length(detail_rows) + 1L]] <- data.frame(
      ref_node = cl$ref_node,
      n_tips = cl$n_tips,
      n_events = cl$n_events,
      emd = emd,
      burst_ref = cl$burst_ref,
      burst_est = burst_c,
      burst_loss = burst_loss,
      centroid_shift = centroid_shift,
      weight = w,
      stringsAsFactors = FALSE
    )
  }
  global_est <- pcr_metric_event_times_relative(tr)
  global_emd <- pcr_metric_wasserstein_1d(global_ref, global_est)
  global_burst_est <- pcr_metric_burstiness_from_events(global_est)
  global_burst_loss <- max(0, (global_burst_ref - global_burst_est) / (global_burst_ref + 1e-12))
  global_error <- (0.35 * global_emd) + (0.65 * global_burst_loss)
  if (w_sum <= 0 || matched == 0L) {
    return(list(
      summary = pcr_empty_pulse_summary(length(panel), global_emd, global_burst_loss, global_error),
      detail = pcr_empty_pulse_detail()
    ))
  }
  mean_emd <- emd_sum / w_sum
  mean_burst_loss <- burst_loss_sum / w_sum
  mean_centroid <- centroid_sum / w_sum
  pulse_error <- (w_emd * mean_emd) + (w_burst_loss * mean_burst_loss) + (w_centroid * mean_centroid)
  coverage <- matched / length(panel)
  selector_error <- ((1 - w_global) * pulse_error) + (w_global * global_error) + (coverage_penalty * (1 - coverage))
  selector_score <- 1 / (1 + selector_error)
  summary <- data.frame(
    matched_clades = matched,
    panel_clades = length(panel),
    coverage = coverage,
    mean_emd = mean_emd,
    mean_burst_loss = mean_burst_loss,
    mean_centroid_shift = mean_centroid,
    local_error = pulse_error,
    global_emd = global_emd,
    global_burst_loss = global_burst_loss,
    global_error = global_error,
    selector_error = selector_error,
    selector_score = selector_score,
    stringsAsFactors = FALSE
  )
  detail <- do.call(rbind, detail_rows)
  list(summary = summary, detail = detail)
}

pcr_metric_edge_signature <- function(tr) {
  n <- ape::Ntip(tr)
  all_tips <- ape::prop.part(tr)
  sig <- character(nrow(tr$edge))
  for (i in seq_len(nrow(tr$edge))) {
    ch <- tr$edge[i, 2]
    if (ch <= n) {
      sig[i] <- paste0('tip:', tr$tip.label[ch])
    } else {
      idx <- ch - n
      tips <- sort(tr$tip.label[all_tips[[idx]]])
      sig[i] <- paste0('clade:', paste(tips, collapse = '|'))
    }
  }
  sig
}

pcr_rate_metrics <- function(ref_tree, dated_tree) {
  ref_df <- data.frame(sig = pcr_metric_edge_signature(ref_tree), ref_branch = ref_tree$edge.length, stringsAsFactors = FALSE)
  est_sig <- pcr_metric_edge_signature(dated_tree)
  est_df <- data.frame(sig = est_sig, time_branch = dated_tree$edge.length, parent_node = dated_tree$edge[,1], child_node = dated_tree$edge[,2], stringsAsFactors = FALSE)
  merged <- merge(est_df, ref_df, by = 'sig', all = FALSE)
  merged <- merged[is.finite(merged$time_branch) & merged$time_branch > 0 & is.finite(merged$ref_branch) & merged$ref_branch > 0, , drop = FALSE]
  if (!nrow(merged)) {
    return(list(summary = data.frame(rate_n_edges = 0L, rate_mean = NA_real_, rate_median = NA_real_, rate_cv = NA_real_, log_rate_sd = NA_real_, log_rate_iqr = NA_real_, extreme_rate_frac = NA_real_, parent_child_pairs = 0L, parent_child_jump_mean = NA_real_, parent_child_jump_q95 = NA_real_, rate_autocorr_spearman = NA_real_, rate_irregularity = NA_real_, stringsAsFactors = FALSE), detail = NULL))
  }
  merged$rate <- merged$ref_branch / merged$time_branch
  merged$log_rate <- log(merged$rate)
  node_to_sig <- stats::setNames(est_sig, dated_tree$edge[,2])
  parent_sig <- unname(node_to_sig[as.character(merged$parent_node)])
  rate_map <- stats::setNames(merged$rate, merged$sig)
  parent_rate <- unname(rate_map[parent_sig])
  keep_pairs <- is.finite(parent_rate) & is.finite(merged$rate) & parent_rate > 0 & merged$rate > 0
  jump <- abs(log(merged$rate[keep_pairs]) - log(parent_rate[keep_pairs]))
  q <- stats::quantile(merged$log_rate, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  iqr_lr <- q[2] - q[1]
  if (is.finite(iqr_lr) && iqr_lr > 0) {
    lo <- q[1] - (1.5 * iqr_lr)
    hi <- q[2] + (1.5 * iqr_lr)
    extreme_frac <- mean(merged$log_rate < lo | merged$log_rate > hi)
  } else {
    extreme_frac <- 0
  }
  if (length(jump) >= 5) {
    autocorr <- suppressWarnings(stats::cor(log(parent_rate[keep_pairs]), log(merged$rate[keep_pairs]), method = 'spearman'))
  } else {
    autocorr <- NA_real_
  }
  autocorr_penalty <- if (is.finite(autocorr)) 1 - max(autocorr, 0) else if (stats::sd(merged$log_rate) < 1e-10) 0 else 1
  summary <- data.frame(
    rate_n_edges = nrow(merged),
    rate_mean = mean(merged$rate),
    rate_median = stats::median(merged$rate),
    rate_cv = stats::sd(merged$rate) / mean(merged$rate),
    log_rate_sd = stats::sd(merged$log_rate),
    log_rate_iqr = stats::IQR(merged$log_rate),
    extreme_rate_frac = extreme_frac,
    parent_child_pairs = sum(keep_pairs),
    parent_child_jump_mean = if (length(jump)) mean(jump) else NA_real_,
    parent_child_jump_q95 = if (length(jump)) as.numeric(stats::quantile(jump, 0.95, names = FALSE)) else NA_real_,
    rate_autocorr_spearman = autocorr,
    rate_irregularity = stats::sd(merged$log_rate) + ifelse(length(jump), mean(jump), 0) + (2 * extreme_frac) + autocorr_penalty,
    stringsAsFactors = FALSE
  )
  list(summary = summary, detail = merged)
}

pcr_gap_metrics <- function(dated_tree, calibrations, tol = 1e-4) {
  age_by_node <- pcr_metric_node_ages(dated_tree)
  detail <- calibrations
  detail$tree_mrca <- vapply(
    seq_len(nrow(detail)),
    function(i) pcr_safe_get_mrca(dated_tree, c(detail$taxonA[i], detail$taxonB[i])),
    integer(1)
  )
  detail$node_age <- as.numeric(age_by_node[as.character(detail$tree_mrca)])
  detail$ghost_gap_ma_raw <- detail$node_age - detail$age_min
  detail$ghost_gap_ma_raw[abs(detail$ghost_gap_ma_raw) < tol] <- 0
  is_point <- is.finite(detail$age_max) & abs(detail$age_max - detail$age_min) < tol
  detail$gap_mode_row <- ifelse(is_point, 'point', 'window')
  detail$ghost_gap_ma <- ifelse(is_point, abs(detail$ghost_gap_ma_raw), pmax(0, detail$ghost_gap_ma_raw))
  detail$ghost_relmin <- ifelse(is.finite(detail$age_min) & detail$age_min > 0,
                                detail$ghost_gap_ma / detail$age_min,
                                NA_real_)
  detail$window_width <- detail$age_max - detail$age_min
  detail$window_position <- ifelse(detail$window_width > 0, (detail$node_age - detail$age_min) / detail$window_width, NA_real_)
  detail$min_violation_ma <- pmax(0, detail$age_min - detail$node_age - tol)
  detail$max_violation_ma <- pmax(0, detail$node_age - detail$age_max - tol)
  detail$within_bounds <- (detail$min_violation_ma == 0) & (detail$max_violation_ma == 0)
  all_point <- all(is_point)
  any_point <- any(is_point)
  gap_mode <- if (all_point) 'point_calibration_slack' else if (any_point) 'mixed_point_and_window' else 'minimum_window_gap'
  summary <- data.frame(
    fossil_n_calibrations = nrow(detail),
    fossil_n_missing_node_age = sum(!is.finite(detail$node_age)),
    ghost_sum_ma = sum(detail$ghost_gap_ma, na.rm = TRUE),
    ghost_mean_ma = pcr_mean_if_any(detail$ghost_gap_ma),
    ghost_median_ma = pcr_median_if_any(detail$ghost_gap_ma),
    mean_relative_gap = pcr_mean_if_any(detail$ghost_relmin),
    relative_gap_median = pcr_median_if_any(detail$ghost_relmin),
    window_position_mean = pcr_mean_if_any(detail$window_position),
    window_position_median = pcr_median_if_any(detail$window_position),
    min_violation_count = sum(detail$min_violation_ma > 0, na.rm = TRUE),
    min_violation_sum_ma = sum(detail$min_violation_ma, na.rm = TRUE),
    max_violation_count = sum(detail$max_violation_ma > 0, na.rm = TRUE),
    max_violation_sum_ma = sum(detail$max_violation_ma, na.rm = TRUE),
    gap_mode = gap_mode,
    stringsAsFactors = FALSE
  )
  list(summary = summary, detail = detail)
}

pcr_subset_calibrations <- function(calibrations, candidate) {
  if (!('candidate' %in% names(calibrations))) return(calibrations)
  keep <- is.na(calibrations$candidate) | calibrations$candidate == '' | calibrations$candidate == 'all' | calibrations$candidate == candidate
  out <- calibrations[keep, , drop = FALSE]
  out$candidate <- NULL
  out
}

pcr_resolve_path <- function(path_value, base_dir) {
  if (grepl('^/', path_value)) return(normalizePath(path_value, winslash = '/', mustWork = TRUE))
  normalizePath(file.path(base_dir, path_value), winslash = '/', mustWork = TRUE)
}

pcr_extract_uncertainty_from_newick <- function(tree_path) {
  txt <- paste(readLines(tree_path, warn = FALSE), collapse = " ")
  if (!nzchar(txt)) return(NULL)

  # Common annotated-Newick forms:
  #   [&height_95%_HPD={1.2,3.4}]
  #   [&age_95%_HPD={1.2,3.4}]
  #   [&HPD={1.2,3.4}]
  #   [&CI={1.2,3.4}]
  pattern <- gregexpr(
    "(?i)([A-Za-z0-9_%.:-]+)\\s*=\\s*[\\{\\[]\\s*([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\\s*,\\s*([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\\s*[\\}\\]]",
    txt,
    perl = TRUE
  )
  hits <- regmatches(txt, pattern)[[1]]
  if (!length(hits)) return(NULL)

  keep_widths <- numeric(0)
  keep_keys <- character(0)

  for (h in hits) {
    m <- regexec(
      "(?i)([A-Za-z0-9_%.:-]+)\\s*=\\s*[\\{\\[]\\s*([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\\s*,\\s*([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\\s*[\\}\\]]",
      h,
      perl = TRUE
    )
    mm <- regmatches(h, m)[[1]]
    if (length(mm) != 4) next
    key <- mm[2]
    lo <- suppressWarnings(as.numeric(mm[3]))
    hi <- suppressWarnings(as.numeric(mm[4]))
    if (!is.finite(lo) || !is.finite(hi) || hi < lo) next

    key_l <- tolower(key)
    include <- grepl("hpd|ci|cred|interval|age|height|time", key_l)
    exclude <- grepl("length|rate|ucld|branch|clock|posterior|prob", key_l)
    if (!include || exclude) next

    keep_widths <- c(keep_widths, hi - lo)
    keep_keys <- c(keep_keys, key)
  }

  if (!length(keep_widths)) return(NULL)

  data.frame(
    uncertainty_source = "tree_embedded_intervals",
    uncertainty_n_bars = length(keep_widths),
    uncertainty_mean_width_ma = mean(keep_widths),
    uncertainty_median_width_ma = stats::median(keep_widths),
    uncertainty_sd_width_ma = if (length(keep_widths) > 1) stats::sd(keep_widths) else 0,
    uncertainty_min_width_ma = min(keep_widths),
    uncertainty_max_width_ma = max(keep_widths),
    uncertainty_q25_width_ma = as.numeric(stats::quantile(keep_widths, 0.25, names = FALSE, type = 8)),
    uncertainty_q75_width_ma = as.numeric(stats::quantile(keep_widths, 0.75, names = FALSE, type = 8)),
    uncertainty_keys_detected = paste(unique(keep_keys), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

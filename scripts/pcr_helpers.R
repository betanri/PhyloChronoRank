pcr_rank_low <- function(x) {
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  out[ok] <- rank(x[ok], ties.method = 'average')
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

## Internode clustering diagnostic (Actions B-C: Manceau et al. 2020 RCS model)
## Identifies radiation zones and computes:
## - internode CV concordance (phylogram vs chronogram)
## - spike-to-background rate ratio within radiation zones
pcr_metric_internode_cv <- function(tr) {
  ## CV of internode distances (branch lengths on internal backbone)
  n <- ape::Ntip(tr)
  if (tr$Nnode < 3L) return(NA_real_)
  int_edges <- which(tr$edge[,2] > n)  ## edges leading to internal nodes
  if (length(int_edges) < 3L) return(NA_real_)
  int_bl <- tr$edge.length[int_edges]
  int_bl <- int_bl[is.finite(int_bl) & int_bl > 0]
  if (length(int_bl) < 3L) return(NA_real_)
  stats::sd(int_bl) / mean(int_bl)
}

pcr_identify_radiation_zones <- function(ref, min_tips = 8L, min_internal = 4L,
                                         internode_ratio_threshold = 0.5) {
  ## Radiation zone: a clade where internal branch lengths are short

  ## relative to the tree-wide mean (internodes clustered together).
  n <- ape::Ntip(ref)
  if (n < min_tips) return(list())
  tree_mean_bl <- mean(ref$edge.length[ref$edge.length > 0], na.rm = TRUE)
  ref_nodes <- (n + 1L):(n + ref$Nnode)
  zones <- list()
  for (nd in ref_nodes) {
    sub <- try(ape::extract.clade(ref, nd), silent = TRUE)
    if (inherits(sub, "try-error") || !inherits(sub, "phylo")) next
    nt <- ape::Ntip(sub)
    if (nt < min_tips || nt >= n) next
    ns <- sub$Nnode
    int_edges <- which(sub$edge[,2] > nt)
    if (length(int_edges) < min_internal) next
    int_bl <- sub$edge.length[int_edges]
    int_bl <- int_bl[is.finite(int_bl) & int_bl > 0]
    if (length(int_bl) < min_internal) next
    mean_int <- mean(int_bl)
    if (mean_int / tree_mean_bl > internode_ratio_threshold) next  ## not a radiation
    zones[[length(zones) + 1L]] <- list(
      node = nd,
      tips = sort(sub$tip.label),
      n_tips = nt,
      n_internal = length(int_bl),
      internode_cv = stats::sd(int_bl) / mean_int,
      mean_internode = mean_int
    )
  }
  zones
}

pcr_score_internode_concordance <- function(ref, tr, zones) {
  ## For each radiation zone, compare internode CV between ref and candidate.
  ## Also compute spike-to-background rate ratio.
  if (!length(zones)) return(list(
    concordance_mean = NA_real_,
    spike_ratio_mean = NA_real_,
    n_zones = 0L,
    detail = data.frame()
  ))
  rows <- list()
  for (z in zones) {
    node_c <- pcr_safe_get_mrca(tr, z$tips)
    if (!is.finite(node_c)) next
    sub_c <- try(ape::extract.clade(tr, node_c), silent = TRUE)
    if (inherits(sub_c, "try-error") || !inherits(sub_c, "phylo")) next
    if (!setequal(sub_c$tip.label, z$tips)) next
    cv_c <- pcr_metric_internode_cv(sub_c)
    if (!is.finite(cv_c)) next
    ## Concordance: ratio of CVs, ideally near 1
    concordance <- if (z$internode_cv > 0) min(cv_c, z$internode_cv) / max(cv_c, z$internode_cv) else NA_real_

    ## Spike ratio: mean rate on internal branches / mean rate across clade
    ## Rate = ref_branch_length / chronogram_branch_length
    sub_r <- try(ape::extract.clade(ref, pcr_safe_get_mrca(ref, z$tips)), silent = TRUE)
    spike_ratio <- NA_real_
    if (inherits(sub_r, "phylo") && inherits(sub_c, "phylo")) {
      nt_c <- ape::Ntip(sub_c)
      int_edges_c <- which(sub_c$edge[,2] > nt_c)
      all_rates <- sub_r$edge.length / pmax(sub_c$edge.length, 1e-12)
      all_rates <- all_rates[is.finite(all_rates)]
      int_rates <- all_rates[int_edges_c[int_edges_c <= length(all_rates)]]
      if (length(int_rates) > 0 && length(all_rates) > 0 && mean(all_rates) > 0) {
        spike_ratio <- mean(int_rates) / mean(all_rates)
      }
    }
    rows[[length(rows) + 1L]] <- data.frame(
      node = z$node, n_tips = z$n_tips,
      cv_ref = z$internode_cv, cv_est = cv_c,
      concordance = concordance, spike_ratio = spike_ratio,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(list(
    concordance_mean = NA_real_, spike_ratio_mean = NA_real_,
    n_zones = 0L, detail = data.frame()
  ))
  detail <- do.call(rbind, rows)
  list(
    concordance_mean = mean(detail$concordance, na.rm = TRUE),
    spike_ratio_mean = mean(detail$spike_ratio, na.rm = TRUE),
    n_zones = nrow(detail),
    detail = detail
  )
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

## Node-depth concordance: R² of relative node depths between ref and candidate.
## Directly captures global distortion (e.g., tipward compression in relaxed models).
## Also computes depth_emd: EMD between the two depth distributions.
pcr_depth_concordance <- function(ref_tree, dated_tree) {
  n_ref <- ape::Ntip(ref_tree)
  n_est <- ape::Ntip(dated_tree)
  if (n_ref != n_est) return(list(depth_r2 = NA_real_, depth_slope = NA_real_, depth_emd = NA_real_))

  ## Get relative node depths (fraction of root-to-tip) for internal nodes
  d_ref <- ape::node.depth.edgelength(ref_tree)
  d_est <- ape::node.depth.edgelength(dated_tree)
  max_ref <- max(d_ref[seq_len(n_ref)])
  max_est <- max(d_est[seq_len(n_est)])
  if (max_ref <= 0 || max_est <= 0) return(list(depth_r2 = NA_real_, depth_slope = NA_real_, depth_emd = NA_real_))

  ## Match internal nodes by descendant tip sets
  ref_sig <- pcr_metric_edge_signature(ref_tree)
  est_sig <- pcr_metric_edge_signature(dated_tree)
  ## Node depths for internal nodes (parent of each edge)
  ref_node_depth <- stats::setNames(d_ref[ref_tree$edge[,2]] / max_ref, ref_sig)
  est_node_depth <- stats::setNames(d_est[dated_tree$edge[,2]] / max_est, est_sig)
  common <- intersect(names(ref_node_depth), names(est_node_depth))
  if (length(common) < 10) return(list(depth_r2 = NA_real_, depth_slope = NA_real_, depth_emd = NA_real_))

  x <- ref_node_depth[common]
  y <- est_node_depth[common]
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 10) return(list(depth_r2 = NA_real_, depth_slope = NA_real_, depth_emd = NA_real_))

  fit <- stats::lm(y ~ x)
  r2 <- summary(fit)$r.squared
  slope <- stats::coef(fit)[2]
  emd <- pcr_metric_wasserstein_1d(sort(x), sort(y))

  list(depth_r2 = r2, depth_slope = unname(slope), depth_emd = emd)
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
    rate_irregularity = stats::sd(merged$log_rate) + (2 * extreme_frac),
    stringsAsFactors = FALSE
  )
  list(summary = summary, detail = merged)
}

## Calibration-induced node compression:
## Detect near-polytomies in the chronogram that are NOT supported by the phylogram.
## Returns fraction of near-zero chronogram internodes not backed by short phylogram internodes.
pcr_compression_score <- function(ref_tree, dated_tree, quantile_threshold = 0.05) {
  n_ref <- ape::Ntip(ref_tree)
  n_est <- ape::Ntip(dated_tree)
  if (n_ref != n_est || n_ref < 10) return(list(compression_score = NA_real_, compression_locality = NA_real_, n_compressed = 0L, n_total_internal = 0L))

  ## Get internal branch lengths
  int_edges_ref <- which(ref_tree$edge[,2] > n_ref)
  int_edges_est <- which(dated_tree$edge[,2] > n_est)
  ref_int_bl <- ref_tree$edge.length[int_edges_ref]
  est_int_bl <- dated_tree$edge.length[int_edges_est]
  ref_int_bl <- ref_int_bl[is.finite(ref_int_bl)]
  est_int_bl <- est_int_bl[is.finite(est_int_bl)]
  if (length(ref_int_bl) < 5 || length(est_int_bl) < 5) return(list(compression_score = NA_real_, compression_locality = NA_real_, n_compressed = 0L, n_total_internal = 0L))

  ## Define thresholds
  est_total_depth <- max(ape::node.depth.edgelength(dated_tree))
  est_near_zero <- est_total_depth * 0.001  ## 0.1% of tree depth
  ref_short_threshold <- as.numeric(stats::quantile(ref_int_bl, quantile_threshold, names = FALSE))

  ## Match internal edges by tip-set signature
  ref_sig <- pcr_metric_edge_signature(ref_tree)
  est_sig <- pcr_metric_edge_signature(dated_tree)
  ref_int_df <- data.frame(sig = ref_sig[int_edges_ref], bl = ref_tree$edge.length[int_edges_ref], stringsAsFactors = FALSE)
  est_int_df <- data.frame(sig = est_sig[int_edges_est], bl = dated_tree$edge.length[int_edges_est], stringsAsFactors = FALSE)
  matched <- merge(est_int_df, ref_int_df, by = 'sig', suffixes = c('_est', '_ref'))
  matched <- matched[is.finite(matched$bl_est) & is.finite(matched$bl_ref), , drop = FALSE]
  if (nrow(matched) < 5) return(list(compression_score = NA_real_, compression_locality = NA_real_, n_compressed = 0L, n_total_internal = nrow(matched)))

  ## Near-zero in chronogram but NOT short in phylogram
  compressed <- matched$bl_est <= est_near_zero & matched$bl_ref > ref_short_threshold
  n_compressed <- sum(compressed)
  score <- n_compressed / nrow(matched)

  list(compression_score = score, compression_locality = NA_real_, n_compressed = n_compressed, n_total_internal = nrow(matched))
}

## Tempo redistribution on non-burst nodes:
## EMD between phylogram and chronogram event-time distributions,
## computed ONLY on nodes outside radiation zones.
pcr_tempo_redistribution_nonburst <- function(ref_tree, dated_tree, radiation_zones) {
  n <- ape::Ntip(ref_tree)
  if (n != ape::Ntip(dated_tree) || n < 10) return(NA_real_)

  ## Collect all tips inside any radiation zone
  zone_tips <- character(0)
  for (z in radiation_zones) zone_tips <- union(zone_tips, z$tips)

  ## Get relative event times for ALL internal nodes
  ref_ev_all <- pcr_metric_event_times_relative(ref_tree)
  est_ev_all <- pcr_metric_event_times_relative(dated_tree)

  if (!length(zone_tips)) {
    return(pcr_metric_wasserstein_1d(ref_ev_all, est_ev_all))
  }

  ## Identify which internal nodes are INSIDE radiation zones
  d_ref <- ape::node.depth.edgelength(ref_tree)
  d_est <- ape::node.depth.edgelength(dated_tree)
  max_ref <- max(d_ref[seq_len(n)])
  max_est <- max(d_est[seq_len(n)])
  if (max_ref <= 0 || max_est <= 0) return(NA_real_)

  ## Mark radiation zone nodes in ref
  zone_nodes_ref <- integer(0)
  for (z in radiation_zones) {
    nd <- z$node
    if (is.finite(nd)) {
      queue <- nd
      while (length(queue)) {
        current <- queue[1]; queue <- queue[-1]
        children <- ref_tree$edge[ref_tree$edge[,1] == current, 2]
        children <- children[children > n]
        zone_nodes_ref <- c(zone_nodes_ref, children)
        queue <- c(queue, children)
      }
      zone_nodes_ref <- c(zone_nodes_ref, nd)
    }
  }
  zone_nodes_ref <- unique(zone_nodes_ref)

  ## Non-burst internal nodes in ref
  all_int_ref <- (n + 1L):(n + ref_tree$Nnode)
  nonburst_ref <- setdiff(all_int_ref, zone_nodes_ref)
  if (length(nonburst_ref) < 5) return(NA_real_)

  ## Get normalized depths for non-burst ref nodes
  ref_depths <- d_ref[nonburst_ref] / max_ref

  ## Match non-burst ref nodes to est nodes by edge signature
  ref_sig <- pcr_metric_edge_signature(ref_tree)
  est_sig <- pcr_metric_edge_signature(dated_tree)
  ## Build node->signature map: for each internal node, find the edge where it is the child
  ref_child_to_edge <- stats::setNames(seq_len(nrow(ref_tree$edge)), ref_tree$edge[,2])
  est_child_to_edge <- stats::setNames(seq_len(nrow(dated_tree$edge)), dated_tree$edge[,2])

  ref_node_sig <- character(length(nonburst_ref))
  for (j in seq_along(nonburst_ref)) {
    nd <- nonburst_ref[j]
    eidx <- ref_child_to_edge[as.character(nd)]
    ref_node_sig[j] <- if (!is.na(eidx)) ref_sig[eidx] else NA_character_
  }

  ## Find matching est depths
  est_sig_map <- stats::setNames(d_est[dated_tree$edge[,2]] / max_est, est_sig)
  matched_est_depths <- est_sig_map[ref_node_sig]
  keep <- !is.na(ref_node_sig) & is.finite(ref_depths) & is.finite(matched_est_depths)
  if (sum(keep) < 5) return(NA_real_)

  pcr_metric_wasserstein_1d(sort(ref_depths[keep]), sort(matched_est_depths[keep]))
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
  ## Non-fixed calibrations: only window constraints are informative for gap scoring
  nonfixed <- !is_point
  mean_rel_gap_nonfixed <- pcr_mean_if_any(detail$ghost_relmin[nonfixed])
  summary <- data.frame(
    fossil_n_calibrations = nrow(detail),
    fossil_n_nonfixed = sum(nonfixed),
    fossil_n_missing_node_age = sum(!is.finite(detail$node_age)),
    ghost_sum_ma = sum(detail$ghost_gap_ma, na.rm = TRUE),
    ghost_mean_ma = pcr_mean_if_any(detail$ghost_gap_ma),
    ghost_median_ma = pcr_median_if_any(detail$ghost_gap_ma),
    mean_relative_gap = mean_rel_gap_nonfixed,
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

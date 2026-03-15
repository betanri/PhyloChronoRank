##############################################################
## chronos_ci_helpers.R
## Delta-method confidence intervals for chronos-dated trees
##
## Extends the Tao, Tamura, Mello & Kumar (2020) delta method
## (originally derived for RelTime) to any penalized-likelihood
## dated tree produced by ape::chronos().
##
## Key idea:
##   chronos returns node ages {T[c]} that were obtained by
##   optimising a penalized log-likelihood.  Although these ages
##   do NOT satisfy the RelTime formula T[c] = T[p]*H[c]/(b+H[c]),
##   we can compute "effective H" values that would produce the
##   same ages under the RelTime formula:
##
##     H_eff[c] = T[c] * b / (T[p] - T[c])
##
##   Plugging H_eff into the Tao et al. variance propagation
##   yields a first-order approximation of the node-age variance
##   that accounts for (i) site-sampling variance and (ii) rate
##   heterogeneity across branches, regardless of whether the
##   ages came from RelTime, chronos, or treePL.
##
## The resulting alpha and beta simplify to:
##   alpha = T[c] / T[p]           (fraction of parent's age)
##   beta  = (T[p] - T[c]) / T[p]  (complement)
##
## Branch-level variance (v_b) follows Tao et al. exactly:
##   v_S(b_j) = b_j / n_sites            [site-sampling]
##   v_R(b_j) = (b_j^2 / sum(b^2)) * RV  [rate heterogeneity]
##   v_b(j)   = v_S + v_R
##
## Dependencies: requires reltime_helpers.R to be sourced first
##   (for .compute_v_b, .compute_vH, .propagate_vt, ci_coverage,
##    ci_mean_width)
##
## Preserves all chronos attributes: PHIIC, ploglik, model,
##   lambda, call, etc.  The CI data.frame is returned alongside
##   the original tree, not embedded in tree attributes.
##############################################################

# ---- compute effective H from any dated tree's node ages ---------
# Given a phylogram (phy) and node ages from any dating method,
# compute H_eff[c] = T[c] * b / (T[p] - T[c]) for each child c.
# For tips (T[c] = 0), H_eff = 0 by definition.
.compute_H_eff <- function(phy, node_age) {
  n_tip <- Ntip(phy)
  total <- n_tip + phy$Nnode
  H_eff <- numeric(total)
  e  <- phy$edge
  el <- phy$edge.length

  for (k in seq_len(nrow(e))) {
    p  <- e[k, 1L]
    cc <- e[k, 2L]
    T_p <- node_age[p]
    T_c <- if (cc <= n_tip) 0 else node_age[cc]
    b   <- el[k]
    dur <- T_p - T_c
    if (dur > 1e-12 && T_c > 1e-12 && b > 1e-16) {
      H_eff[cc] <- T_c * b / dur
    }
  }
  H_eff
}

# ---- extract node ages from a dated phylo object ----------------
.extract_node_ages <- function(dated_tree) {
  n_tip  <- Ntip(dated_tree)
  n_node <- dated_tree$Nnode
  total  <- n_tip + n_node

  depth <- ape::node.depth.edgelength(dated_tree)
  if (length(depth) != total || any(!is.finite(depth))) {
    stop("Could not derive finite root-to-node depths from dated tree.")
  }

  max_tip_depth <- max(depth[seq_len(n_tip)], na.rm = TRUE)
  node_age <- pmax(max_tip_depth - depth, 0)
  node_age[seq_len(n_tip)] <- 0
  node_age
}

##############################################################
## Public API
##############################################################

#' Compute delta-method CIs for a dated tree using effective-H
#'
#' This is the core CI function.  It takes any phylogram and a
#' vector of node ages (from chronos, treePL, or any method) and
#' returns CI bounds for every internal node.
#'
#' @param phy       Original phylogram (branch lengths in substitutions).
#' @param node_age  Numeric vector of length (n_tip + n_node) giving
#'                  the age of each node; tips should be 0 for
#'                  ultrametric trees.
#' @param n_sites   Alignment length for site-sampling variance.
#' @param root_var  Variance of the root age (0 = exact root).
#' @param cal_min   Named list: minimum age bounds for CI truncation.
#' @param cal_max   Named list: maximum age bounds for CI truncation.
#' @return          data.frame(node, age, ci_lo, ci_hi, se)
chronos_ci <- function(phy, dated_tree, n_sites = 1000L,
                       root_var = 0, cal_min = NULL, cal_max = NULL) {
  n_tip     <- Ntip(phy)
  n_node    <- phy$Nnode
  root_node <- n_tip + 1L
  total     <- n_tip + n_node

  ## Extract node ages from the dated tree
  node_age <- .extract_node_ages(dated_tree)

  ## Compute ntip_s (tip counts per node)
  ntip_s <- integer(total)
  ntip_s[seq_len(n_tip)] <- 1L
  e <- phy$edge
  for (k in postorder(phy)) {
    ntip_s[e[k, 1L]] <- ntip_s[e[k, 1L]] + ntip_s[e[k, 2L]]
  }

  ## Branch-level variance (uses phylogram branch lengths + node ages)
  vb_res <- .compute_v_b(phy, node_age, n_sites)
  v_b    <- vb_res$v_b

  ## Variance of H (bottom-up; topology- and v_b-dependent only)
  vH <- .compute_vH(phy, v_b, ntip_s)

  ## Effective H from the dating method's node ages
  H_eff <- .compute_H_eff(phy, node_age)

  ## Propagate variance (top-down) using effective H
  v_t <- .propagate_vt(phy, H_eff, node_age, vH, v_b, v_root = root_var)

  ## Build CIs
  se    <- sqrt(pmax(v_t, 0))
  ci_lo <- pmax(node_age - 1.96 * se, 0)
  ci_hi <- node_age + 1.96 * se

  ## Truncate at calibration bounds
  if (!is.null(cal_min)) {
    for (nm in names(cal_min)) {
      nd <- suppressWarnings(as.integer(nm))
      if (is.finite(nd) && nd > n_tip && nd <= total)
        ci_lo[nd] <- max(ci_lo[nd], cal_min[[nm]])
    }
  }
  if (!is.null(cal_max)) {
    for (nm in names(cal_max)) {
      nd <- suppressWarnings(as.integer(nm))
      if (is.finite(nd) && nd > n_tip && nd <= total)
        ci_hi[nd] <- min(ci_hi[nd], cal_max[[nm]])
    }
  }
  ci_hi <- pmax(ci_hi, ci_lo)

  ## Fix root if exact calibration
  if (!is.finite(root_var) || root_var <= 0) {
    ci_lo[root_node] <- node_age[root_node]
    ci_hi[root_node] <- node_age[root_node]
  }

  int_nodes <- seq.int(n_tip + 1L, total)
  data.frame(
    node  = int_nodes,
    age   = node_age[int_nodes],
    ci_lo = ci_lo[int_nodes],
    ci_hi = ci_hi[int_nodes],
    se    = se[int_nodes]
  )
}

#' Summarize CI widths into PCR's optional uncertainty-table schema
#'
#' @param candidate Candidate label carried into PCR outputs.
#' @param ci_df     data.frame with ci_lo and ci_hi columns.
#' @param source    Provenance string for uncertainty_source.
#' @return          One-row data.frame ready for rbind across candidates.
ci_width_summary <- function(candidate, ci_df, source = "delta_method_ci") {
  need <- c("ci_lo", "ci_hi")
  if (!all(need %in% names(ci_df))) {
    stop("CI data frame must contain ci_lo and ci_hi columns.")
  }

  widths <- as.numeric(ci_df$ci_hi) - as.numeric(ci_df$ci_lo)
  widths <- widths[is.finite(widths) & widths >= 0]
  if (!length(widths)) {
    return(data.frame(
      candidate = candidate,
      uncertainty_source = source,
      uncertainty_n_bars = 0L,
      uncertainty_mean_width_ma = NA_real_,
      uncertainty_median_width_ma = NA_real_,
      uncertainty_sd_width_ma = NA_real_,
      uncertainty_min_width_ma = NA_real_,
      uncertainty_max_width_ma = NA_real_,
      uncertainty_q25_width_ma = NA_real_,
      uncertainty_q75_width_ma = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    candidate = candidate,
    uncertainty_source = source,
    uncertainty_n_bars = length(widths),
    uncertainty_mean_width_ma = mean(widths),
    uncertainty_median_width_ma = stats::median(widths),
    uncertainty_sd_width_ma = if (length(widths) > 1L) stats::sd(widths) else 0,
    uncertainty_min_width_ma = min(widths),
    uncertainty_max_width_ma = max(widths),
    uncertainty_q25_width_ma = as.numeric(stats::quantile(widths, 0.25, names = FALSE, type = 8)),
    uncertainty_q75_width_ma = as.numeric(stats::quantile(widths, 0.75, names = FALSE, type = 8)),
    stringsAsFactors = FALSE
  )
}

#' Run chronos model-selection and compute delta-method CIs
#'
#' Wrapper that calls run_chronos_modelselect (must already exist in
#' the session, from the benchmark runner), then adds CIs.
#'
#' @param phy           Phylogram.
#' @param calib_bundle  Calibration bundle from build_calibration_bundle.
#' @param n_sites       Alignment length for CI variance.
#' @param ...           Passed to run_chronos_modelselect.
#' @return              The cfit list from run_chronos_modelselect,
#'                      augmented with $ci (data.frame) and
#'                      $ci_coverage, $ci_width (if bt_true supplied).
run_chronos_with_ci <- function(phy, calib_bundle, n_sites = 1000L,
                                bt_true = NULL, ...) {
  if (!exists("run_chronos_modelselect", mode = "function")) {
    stop("run_chronos_modelselect not found; source the benchmark runner first.")
  }
  cfit <- run_chronos_modelselect(phy, calib_bundle = calib_bundle, ...)
  if (is.null(cfit$tree)) return(cfit)

  n_tip     <- Ntip(phy)
  root_node <- n_tip + 1L

  ## Derive root_var from calibration
  calib <- calib_bundle$chronos_calib
  root_calib <- calib[calib$node == root_node, , drop = FALSE]
  root_var <- 0
  if (nrow(root_calib) &&
      is.finite(root_calib$age.min[1]) && is.finite(root_calib$age.max[1]) &&
      abs(root_calib$age.max[1] - root_calib$age.min[1]) > 1e-10) {
    root_var <- ((root_calib$age.max[1] - root_calib$age.min[1]) / (2 * 1.96))^2
  }

  ## Build cal_min / cal_max for CI truncation
  cal_min <- cal_max <- NULL
  for (i in seq_len(nrow(calib))) {
    kk <- as.character(calib$node[i])
    cal_min[[kk]] <- calib$age.min[i]
    if (is.finite(calib$age.max[i])) cal_max[[kk]] <- calib$age.max[i]
  }

  ci_df <- chronos_ci(phy, cfit$tree,
                       n_sites = n_sites,
                       root_var = root_var,
                       cal_min = cal_min, cal_max = cal_max)
  cfit$ci <- ci_df

  ## Optional coverage/width if ground-truth provided
  if (!is.null(bt_true)) {
    cfit$ci_coverage <- ci_coverage(ci_df, bt_true)
    cfit$ci_width    <- ci_mean_width(ci_df)
  }
  cfit$ci_method <- "delta_method_effective_H"
  cfit
}

#' Apply chronos_ci to a treePL-dated tree
#'
#' Convenience wrapper for computing CIs on a treePL output.
#' treePL, like chronos, produces a dated tree from penalized
#' likelihood; the same effective-H delta method applies.
#'
#' @param phy          Original phylogram.
#' @param treepl_tree  treePL output (dated phylo).
#' @param calib_bundle Calibration bundle (for bounds).
#' @param n_sites      Alignment length.
#' @return             data.frame(node, age, ci_lo, ci_hi, se)
treepl_ci <- function(phy, treepl_tree, calib_bundle, n_sites = 1000L) {
  if (is.null(treepl_tree)) return(NULL)

  n_tip     <- Ntip(phy)
  root_node <- n_tip + 1L
  calib <- calib_bundle$chronos_calib
  root_calib <- calib[calib$node == root_node, , drop = FALSE]
  root_var <- 0
  if (nrow(root_calib) &&
      is.finite(root_calib$age.min[1]) && is.finite(root_calib$age.max[1]) &&
      abs(root_calib$age.max[1] - root_calib$age.min[1]) > 1e-10) {
    root_var <- ((root_calib$age.max[1] - root_calib$age.min[1]) / (2 * 1.96))^2
  }

  cal_min <- cal_max <- NULL
  for (i in seq_len(nrow(calib))) {
    kk <- as.character(calib$node[i])
    cal_min[[kk]] <- calib$age.min[i]
    if (is.finite(calib$age.max[i])) cal_max[[kk]] <- calib$age.max[i]
  }

  chronos_ci(phy, treepl_tree,
             n_sites = n_sites,
             root_var = root_var,
             cal_min = cal_min, cal_max = cal_max)
}

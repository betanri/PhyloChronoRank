##############################################################
## reltime_helpers.R
## RelTime point-estimate + confidence interval implementation
##
## Point estimates: Tamura et al. (2012) Mol. Biol. Evol. 29:1893-1907
##                  Tamura et al. (2018) Mol. Biol. Evol. 35:1770-1782
## Confidence intervals: Tao, Tamura, Mello & Kumar (2020)
##                       Mol. Biol. Evol. 37:280-290
##
## Core algorithm (point estimates)
## ---------------------------------
## For each node n with direct parent p and edge branch length b:
##   H[n] = weighted-mean substitution depth from n to descendant tips
##         (tip-count-weighted, so H[tip] = 0 by definition)
##   T[n] = T[p] * H[n] / (b + H[n])
##
## Derivation of the formula:
##   Under the local-clock approximation, the rate on branch (p→n)
##   equals the average rate within n's subtree: r_n = H[n] / T[n].
##   Branch length: b = r_n * (T[p] - T[n]) = (H[n]/T[n]) * (T[p]-T[n]).
##   Solving for T[n]: T[n] = T[p] * H[n] / (b + H[n]).
##
##   For a strict clock all H values are proportional to true ages and
##   the formula recovers exact node ages.
##
## Confidence intervals (Tao et al. 2020, delta method)
## -----------------------------------------------------
##   v(b_j) = v_S(b_j) + v_R(b_j)
##     v_S(b_j) = b_j / n_sites          [site-sampling variance]
##     v_R(b_j) = (b_j^2 / sum(b^2)) * RV(R)  [rate-heterogeneity component]
##
##   Rate variance:
##     r_j     = b_j / dur_j             [implied rate per branch]
##     V_obs   = mean((r_j - r_bar)^2)   [observed rate variance]
##     sv(r_j) = v_S(b_j) / dur_j^2     [sampling variance of rate j]
##     SV(R)   = mean(sv(r_j))
##     RV(R)   = max(0, V_obs - SV(R))  [true rate variance]
##
##   Variance propagation (pre-order):
##     h_c  = b_c + H[c]
##     a    = H[c] / h_c   (alpha)
##     be   = b_c / h_c    (beta)
##     v(T[c]) = a^2 * v(T[p])
##             + (T[p]*a / h_c)^2 * v_b[c]
##             + (T[p]*be / h_c)^2 * v(H[c])
##
##   v(H[c]) propagated bottom-up (post-order):
##     v(H[p]) = sum_c  (n_c/n_p)^2 * (v_b[c] + v(H[c]))
##
##   CI = T[n] +/- 1.96 * sqrt(v(T[n])), truncated at calibration bounds.
##   Output is reported for internal nodes only; tip ages are fixed at 0.
##############################################################

# ---- helper: build children lookup ----------------------------------------
.reltime_children <- function(phy) {
  total <- Ntip(phy) + phy$Nnode
  ch <- vector("list", total)
  e  <- phy$edge
  for (k in seq_len(nrow(e))) ch[[e[k, 1L]]] <- c(ch[[e[k, 1L]]], e[k, 2L])
  ch
}

# ---- core: compute H[n] (post-order) ---------------------------------------
# H[n] = tip-count-weighted MEAN substitution depth from n to its descendant tips.
# Key: when computing the contribution of child c to parent p, we need the MEAN
# path from c to its tips (H_sum[c] / ntip_s[c]), not the raw accumulated sum.
.compute_H <- function(phy) {
  n_tip <- Ntip(phy)
  total <- n_tip + phy$Nnode
  H_sum  <- numeric(total)   # raw accumulated sum (divided later)
  ntip_s <- integer(total)
  ntip_s[seq_len(n_tip)] <- 1L
  e  <- phy$edge
  el <- phy$edge.length
  for (k in postorder(phy)) {
    p <- e[k, 1L]; c <- e[k, 2L]
    # At this point ntip_s[c] is complete (all of c's subtree has been processed).
    # Use the MEAN path from c to its tips as the contribution to p.
    H_c_mean <- if (c <= n_tip || ntip_s[c] == 0L) 0 else H_sum[c] / ntip_s[c]
    H_sum[p] <- H_sum[p] + (el[k] + H_c_mean) * ntip_s[c]
    ntip_s[p] <- ntip_s[p] + ntip_s[c]
  }
  H <- numeric(total)
  for (nd in seq.int(n_tip + 1L, total)) {
    if (ntip_s[nd] > 0L) H[nd] <- H_sum[nd] / ntip_s[nd]
  }
  list(H = H, ntip_s = ntip_s)
}

# ---- core: compute v(H[n]) given v_b per edge (post-order) -----------------
# v(H[p]) = sum_c (n_c/n_p)^2 * (v_b[c] + v(H[c]))
# vH stores the VARIANCE OF THE MEAN, consistent with H storing the mean.
.compute_vH <- function(phy, v_b, ntip_s) {
  n_tip <- Ntip(phy)
  total <- n_tip + phy$Nnode
  vH_sum <- numeric(total)   # raw accumulated variance sum
  e  <- phy$edge
  for (k in postorder(phy)) {
    p <- e[k, 1L]; c <- e[k, 2L]
    np <- ntip_s[p]; nc <- ntip_s[c]
    if (np == 0L || nc == 0L) next
    # vH[c] is the variance of H[c] (mean); v_b[k] is variance of branch c's length
    vH_c <- if (c <= n_tip) 0 else vH_sum[c] / (nc * nc)
    vH_sum[p] <- vH_sum[p] + nc * nc * (v_b[k] + vH_c)
  }
  # Convert to variance of mean H[p]
  vH <- numeric(total)
  for (nd in seq.int(n_tip + 1L, total)) {
    n <- ntip_s[nd]
    if (n > 0L) vH[nd] <- vH_sum[nd] / (n * n)
  }
  vH
}

# ---- core: assign node ages (pre-order) ------------------------------------
# Uses rev(postorder()) which is a valid topological order (parents before
# children) for any rooted binary tree; no dependency on ape::preorder().
.assign_node_ages <- function(phy, H, root_age) {
  n_tip     <- Ntip(phy)
  total     <- n_tip + phy$Nnode
  root_node <- n_tip + 1L
  e  <- phy$edge
  el <- phy$edge.length
  node_age  <- numeric(total)
  node_age[root_node] <- root_age
  for (k in rev(postorder(phy))) {
    p <- e[k, 1L]; c <- e[k, 2L]
    b   <- el[k]
    h_c <- b + H[c]
    node_age[c] <- if (h_c > 1e-12) node_age[p] * H[c] / h_c else 0
    node_age[c] <- max(node_age[c], 0)
  }
  node_age
}

# ---- core: propagate variance (pre-order) ----------------------------------
.propagate_vt <- function(phy, H, node_age, vH, v_b, v_root = 0) {
  n_tip     <- Ntip(phy)
  total     <- n_tip + phy$Nnode
  root_node <- n_tip + 1L
  e  <- phy$edge
  el <- phy$edge.length
  v_t              <- numeric(total)
  v_t[root_node]   <- v_root
  for (k in rev(postorder(phy))) {
    p   <- e[k, 1L]; c <- e[k, 2L]
    b   <- el[k]
    h_c <- b + H[c]
    if (h_c < 1e-12) next
    T_p  <- node_age[p]
    a    <- H[c] / h_c          # alpha = dTc/dTp
    be   <- b   / h_c           # beta  = dTc/dH[c] * h_c/T_p
    v_t[c] <- a * a * v_t[p] +
              (T_p * a / h_c)^2 * v_b[k] +
              (T_p * be / h_c)^2 * vH[c]
  }
  v_t
}

# ---- apply internal minimum-bound calibrations -----------------------------
.apply_min_cals <- function(node_age, internal_cal, children_of, n_tip) {
  if (is.null(internal_cal) || nrow(internal_cal) == 0L) return(node_age)
  ic <- internal_cal[order(internal_cal$min_age, decreasing = TRUE), , drop = FALSE]
  for (i in seq_len(nrow(ic))) {
    nd    <- ic$node[i]
    min_a <- ic$min_age[i]
    if (!is.finite(min_a) || !is.finite(nd)) next
    if (is.na(node_age[nd]) || node_age[nd] >= min_a - 1e-8) next
    sf    <- min_a / max(node_age[nd], 1e-12)
    stack <- nd
    while (length(stack) > 0L) {
      cur   <- stack[1L]; stack <- stack[-1L]
      node_age[cur] <- node_age[cur] * sf
      for (cc in children_of[[cur]]) if (!is.na(cc) && cc > n_tip) stack <- c(stack, cc)
    }
  }
  node_age
}

# ---- compute branch-level variance (v_b) -----------------------------------
#  n_sites : alignment length (integer); controls site-sampling variance.
#            In simulations without a real alignment, use a representative
#            value (default 1000 bp; larger values → narrower CIs).
.compute_v_b <- function(phy, node_age, n_sites) {
  n_tip <- Ntip(phy)
  e  <- phy$edge
  el <- phy$edge.length
  m  <- nrow(e)

  # branch durations
  dur <- numeric(m)
  for (k in seq_len(m)) {
    p <- e[k, 1L]; c <- e[k, 2L]
    T_p <- node_age[p]
    T_c <- if (c <= n_tip) 0 else node_age[c]
    dur[k] <- max(T_p - T_c, 1e-12)
  }

  # implied rates and rate variance
  rates <- el / dur
  r_bar <- mean(rates, na.rm = TRUE)
  V_obs <- mean((rates - r_bar)^2, na.rm = TRUE)

  # site-sampling variance of each branch
  v_s   <- el / n_sites

  # sampling variance of each implied rate: sv(r_j) = v_S(b_j) / dur_j^2
  sv_r  <- v_s / dur^2
  SV_R  <- mean(sv_r, na.rm = TRUE)
  RV_R  <- max(0, V_obs - SV_R)

  # rate-heterogeneity component
  sum_b2 <- sum(el^2, na.rm = TRUE)
  v_r    <- if (sum_b2 > 1e-12) (el^2 / sum_b2) * RV_R else numeric(m)

  list(v_b = v_s + v_r, rates = rates, RV_R = RV_R, V_obs = V_obs)
}

##############################################################
## Public API
##############################################################

#' RelTime point-estimate dated tree
#'
#' @param phy          Rooted phylogram (branch lengths in substitutions).
#' @param root_age     Calibrated root age (time units, e.g. Ma).
#' @param internal_cal Optional data.frame(node = integer, min_age = numeric)
#'                     for minimum-bound internal calibrations (Benchmark C).
#' @return             Dated phylo object (branch lengths in time units),
#'                     or NULL on failure.
run_reltime <- function(phy, root_age, internal_cal = NULL) {
  n_tip     <- Ntip(phy)
  n_node    <- phy$Nnode
  root_node <- n_tip + 1L
  e  <- phy$edge
  el <- phy$edge.length

  if (!is.finite(root_age) || root_age <= 0) return(NULL)
  if (is.null(el) || any(!is.finite(el)) || any(el <= 0)) return(NULL)

  res      <- .compute_H(phy)
  H        <- res$H
  H_root   <- H[root_node]
  if (!is.finite(H_root) || H_root < 1e-12) return(NULL)

  node_age <- .assign_node_ages(phy, H, root_age)

  if (!is.null(internal_cal) && nrow(internal_cal) > 0L) {
    ch       <- .reltime_children(phy)
    node_age <- .apply_min_cals(node_age, internal_cal, ch, n_tip)
  }

  dated <- phy
  for (k in seq_len(nrow(e))) {
    p <- e[k, 1L]; c <- e[k, 2L]
    T_p <- node_age[p]
    T_c <- if (c <= n_tip) 0 else node_age[c]
    dated$edge.length[k] <- max(T_p - T_c, 1e-12)
  }

  bt <- try(branching.times(dated), silent = TRUE)
  if (inherits(bt, "try-error")) return(NULL)
  if (any(!is.finite(bt)) || any(bt < -1e-8) || any(bt > root_age + 1e-6)) return(NULL)
  dated
}

#' RelTime confidence intervals (Tao, Tamura, Mello & Kumar 2020)
#'
#' Must be called AFTER run_reltime() so that node_age values are available.
#'
#' @param phy          The same rooted phylogram used in run_reltime().
#' @param node_age     Vector of node ages (length n_tip + n_node), returned
#'                     by .assign_node_ages() or extracted from the dated tree.
#' @param n_sites      Alignment length used to estimate v_S(b_j) = b_j/n.
#'                     Set to Inf to suppress site-sampling variance entirely.
#' @param root_var     Variance of the root age (0 for exact calibration,
#'                     ((root_max - root_min)/(2*1.96))^2 for a root window).
#' @param cal_min      Named numeric: minimum age bound per calibration node
#'                     (used to truncate CI lower bounds).
#' @param cal_max      Named numeric: maximum age bound per calibration node.
#' @return             data.frame with columns: node, age, ci_lo, ci_hi, se
#'                     for internal nodes only; tips are omitted because their
#'                     age is fixed at 0 by convention.
reltime_ci <- function(phy, node_age, n_sites = 1000L,
                       root_var = 0, cal_min = NULL, cal_max = NULL) {
  n_tip     <- Ntip(phy)
  n_node    <- phy$Nnode
  root_node <- n_tip + 1L

  if (!is.finite(n_sites) || n_sites <= 0) n_sites <- Inf

  res  <- .compute_H(phy)
  H    <- res$H
  ntip_s <- res$ntip_s

  vb_res <- .compute_v_b(phy, node_age, n_sites)
  v_b    <- vb_res$v_b

  vH  <- .compute_vH(phy, v_b, ntip_s)
  v_t <- .propagate_vt(phy, H, node_age, vH, v_b, v_root = root_var)

  se      <- sqrt(pmax(v_t, 0))
  ci_lo   <- pmax(node_age - 1.96 * se, 0)
  ci_hi   <- node_age + 1.96 * se

  # truncate at calibration bounds
  if (!is.null(cal_min)) {
    for (nm in names(cal_min)) {
      nd <- suppressWarnings(as.integer(nm))
      if (is.finite(nd) && nd > n_tip)
        ci_lo[nd] <- max(ci_lo[nd], cal_min[[nm]])
    }
  }
  if (!is.null(cal_max)) {
    for (nm in names(cal_max)) {
      nd <- suppressWarnings(as.integer(nm))
      if (is.finite(nd) && nd > n_tip)
        ci_hi[nd] <- min(ci_hi[nd], cal_max[[nm]])
    }
  }
  # Keep the root fixed only when there is no root-age uncertainty.
  if (!is.finite(root_var) || root_var <= 0) {
    ci_lo[root_node] <- node_age[root_node]
    ci_hi[root_node] <- node_age[root_node]
  }

  int_nodes <- seq.int(n_tip + 1L, n_tip + n_node)
  data.frame(
    node  = int_nodes,
    age   = node_age[int_nodes],
    ci_lo = ci_lo[int_nodes],
    ci_hi = ci_hi[int_nodes],
    se    = se[int_nodes]
  )
}

#' Run RelTime and compute CIs in one call
#'
#' @param phy          Rooted phylogram.
#' @param root_age     Root calibration age (exact if root_window is NULL).
#' @param internal_cal Optional data.frame(node, min_age) for Benchmark C.
#' @param n_sites      Alignment length for CI variance (default 1000).
#' @param root_window  Optional c(min, max) for root uncertainty. If provided,
#'                     root_age is set to midpoint and root variance is derived
#'                     from window width.
#' @return             List with $tree (dated phylo), $ci (internal-node CI
#'                     data.frame), $rates (per-branch implied rates),
#'                     $RV_R (rate variance).
run_reltime_with_ci <- function(phy, root_age, internal_cal = NULL,
                                n_sites = 1000L, root_window = NULL) {
  n_tip     <- Ntip(phy)
  root_node <- n_tip + 1L
  e  <- phy$edge
  el <- phy$edge.length

  # root uncertainty
  root_var <- 0
  if (!is.null(root_window) && length(root_window) == 2L) {
    root_age <- mean(root_window)
    root_var <- ((root_window[2L] - root_window[1L]) / (2 * 1.96))^2
  }

  if (!is.finite(root_age) || root_age <= 0) return(NULL)
  if (is.null(el) || any(!is.finite(el)) || any(el <= 0)) return(NULL)

  res      <- .compute_H(phy)
  H        <- res$H
  ntip_s   <- res$ntip_s
  H_root   <- H[root_node]
  if (!is.finite(H_root) || H_root < 1e-12) return(NULL)

  node_age <- .assign_node_ages(phy, H, root_age)

  ch <- .reltime_children(phy)
  if (!is.null(internal_cal) && nrow(internal_cal) > 0L) {
    node_age <- .apply_min_cals(node_age, internal_cal, ch, n_tip)
  }

  # Build calibration bound lookups for CI truncation
  cal_min <- cal_max <- NULL
  if (!is.null(root_window) && length(root_window) == 2L && all(is.finite(root_window))) {
    root_window <- sort(root_window)
    cal_min[[as.character(root_node)]] <- root_window[1L]
    cal_max[[as.character(root_node)]] <- root_window[2L]
  } else {
    cal_min[[as.character(root_node)]] <- root_age
    cal_max[[as.character(root_node)]] <- root_age
  }
  if (!is.null(internal_cal) && nrow(internal_cal) > 0L) {
    for (i in seq_len(nrow(internal_cal))) {
      key <- as.character(internal_cal$node[i])
      cal_min[[key]] <- internal_cal$min_age[i]
    }
  }

  # CI computation
  vb_res <- .compute_v_b(phy, node_age, n_sites)
  v_b    <- vb_res$v_b
  vH     <- .compute_vH(phy, v_b, ntip_s)
  v_t    <- .propagate_vt(phy, H, node_age, vH, v_b, v_root = root_var)

  se    <- sqrt(pmax(v_t, 0))
  ci_lo <- pmax(node_age - 1.96 * se, 0)
  ci_hi <- node_age + 1.96 * se

  # Truncate at calibration bounds
  if (!is.null(cal_min)) {
    for (nm in names(cal_min)) {
      nd <- suppressWarnings(as.integer(nm))
      if (is.finite(nd)) ci_lo[nd] <- max(ci_lo[nd], cal_min[[nm]])
    }
  }
  if (!is.null(cal_max)) {
    for (nm in names(cal_max)) {
      nd <- suppressWarnings(as.integer(nm))
      if (is.finite(nd)) ci_hi[nd] <- min(ci_hi[nd], cal_max[[nm]])
    }
  }
  if (!is.finite(root_var) || root_var <= 0) {
    ci_lo[root_node] <- node_age[root_node]
    ci_hi[root_node] <- node_age[root_node]
  }

  # Build dated tree
  dated <- phy
  for (k in seq_len(nrow(e))) {
    p <- e[k, 1L]; c <- e[k, 2L]
    T_p <- node_age[p]
    T_c <- if (c <= n_tip) 0 else node_age[c]
    dated$edge.length[k] <- max(T_p - T_c, 1e-12)
  }

  bt <- try(branching.times(dated), silent = TRUE)
  if (inherits(bt, "try-error")) return(NULL)
  if (any(!is.finite(bt)) || any(bt < -1e-8) || any(bt > root_age + 1e-6)) return(NULL)

  int_nodes <- seq.int(n_tip + 1L, n_tip + phy$Nnode)
  ci_df <- data.frame(
    node  = int_nodes,
    age   = node_age[int_nodes],
    ci_lo = ci_lo[int_nodes],
    ci_hi = ci_hi[int_nodes],
    se    = se[int_nodes]
  )

  list(
    tree  = dated,
    ci    = ci_df,
    node_age = node_age,
    rates = vb_res$rates,
    RV_R  = vb_res$RV_R,
    V_obs = vb_res$V_obs
  )
}

#' CI coverage: fraction of internal nodes whose true age falls in [ci_lo, ci_hi]
ci_coverage <- function(ci_df, bt_true) {
  common <- intersect(as.character(ci_df$node), names(bt_true))
  if (!length(common)) return(NA_real_)
  idx <- match(common, as.character(ci_df$node))
  true_ages <- bt_true[common]
  in_ci <- true_ages >= ci_df$ci_lo[idx] - 1e-8 & true_ages <= ci_df$ci_hi[idx] + 1e-8
  mean(in_ci, na.rm = TRUE)
}

#' Mean CI width for internal nodes
ci_mean_width <- function(ci_df) {
  mean(ci_df$ci_hi - ci_df$ci_lo, na.rm = TRUE)
}

#' Merge pairwise calibration rows onto tree nodes using the same duplicate-node
#' intersection rule used in the local chronos/treePL pipelines.
#'
#' @param phy            Rooted phylogram.
#' @param calibration_df data.frame with taxonA,taxonB,age_min,age_max columns.
#' @return               data.frame(node, age_min, age_max, n_merged)
reltime_merge_calibration_bounds <- function(phy, calibration_df) {
  need <- c("taxonA", "taxonB", "age_min")
  if (!all(need %in% names(calibration_df))) {
    stop("Calibration data frame must contain taxonA, taxonB, age_min.")
  }
  cal <- calibration_df
  if (!("age_max" %in% names(cal))) cal$age_max <- Inf
  cal$taxonA <- as.character(cal$taxonA)
  cal$taxonB <- as.character(cal$taxonB)
  cal$age_min <- as.numeric(cal$age_min)
  cal$age_max <- as.numeric(cal$age_max)
  cal$age_max[is.na(cal$age_max)] <- Inf
  cal <- cal[is.finite(cal$age_min), , drop = FALSE]
  if (!nrow(cal)) {
    return(data.frame(
      node = integer(0),
      age_min = numeric(0),
      age_max = numeric(0),
      n_merged = integer(0)
    ))
  }

  mapped <- list()
  for (i in seq_len(nrow(cal))) {
    a <- cal$taxonA[i]
    b <- cal$taxonB[i]
    if (!(a %in% phy$tip.label) || !(b %in% phy$tip.label)) next
    nd <- ape::getMRCA(phy, c(a, b))
    if (is.null(nd) || !is.finite(nd)) next
    mapped[[length(mapped) + 1L]] <- data.frame(
      node = as.integer(nd),
      age_min = cal$age_min[i],
      age_max = cal$age_max[i],
      stringsAsFactors = FALSE
    )
  }
  if (!length(mapped)) {
    return(data.frame(
      node = integer(0),
      age_min = numeric(0),
      age_max = numeric(0),
      n_merged = integer(0)
    ))
  }

  mapped_df <- do.call(rbind, mapped)
  merged <- list()
  for (nd in sort(unique(mapped_df$node))) {
    idx <- mapped_df$node == nd
    mn <- max(mapped_df$age_min[idx], na.rm = TRUE)
    mx <- min(mapped_df$age_max[idx], na.rm = TRUE)
    if (!is.finite(mn) || is.na(mx) || mn > mx) next
    merged[[length(merged) + 1L]] <- data.frame(
      node = as.integer(nd),
      age_min = mn,
      age_max = mx,
      n_merged = sum(idx),
      stringsAsFactors = FALSE
    )
  }
  if (!length(merged)) {
    return(data.frame(
      node = integer(0),
      age_min = numeric(0),
      age_max = numeric(0),
      n_merged = integer(0)
    ))
  }
  do.call(rbind, merged)
}

# ---- project node ages onto full merged calibration bounds -----------------
# Solves a quadratic program:
#   minimize   sum_i (x_i - t0_i)^2
#   subject to lower_i <= x_i <= upper_i
#              x_parent - x_child >= eps
# for internal nodes i only, with tip ages fixed at 0.
.reltime_project_node_ages_qp <- function(phy, node_age_init, lower, upper, eps = 1e-6) {
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("Package quadprog is required for full-bound RelTime projection.")
  }
  n_tip <- Ntip(phy)
  int_nodes <- seq.int(n_tip + 1L, n_tip + phy$Nnode)
  p <- length(int_nodes)
  idx_of <- stats::setNames(seq_along(int_nodes), int_nodes)

  eq_cols <- list()
  eq_b <- numeric(0)
  ge_cols <- list()
  ge_b <- numeric(0)

  add_eq <- function(v, b) {
    eq_cols[[length(eq_cols) + 1L]] <<- v
    eq_b <<- c(eq_b, b)
  }
  add_ge <- function(v, b) {
    ge_cols[[length(ge_cols) + 1L]] <<- v
    ge_b <<- c(ge_b, b)
  }

  for (nd in int_nodes) {
    v <- numeric(p)
    v[idx_of[[as.character(nd)]]] <- 1
    lo <- lower[nd]
    hi <- upper[nd]
    if (is.finite(lo) && is.finite(hi) && abs(lo - hi) < 1e-10) {
      add_eq(v, lo)
    } else {
      if (is.finite(lo)) add_ge(v, lo)
      if (is.finite(hi)) add_ge(-v, -hi)
    }
  }

  for (k in seq_len(nrow(phy$edge))) {
    pnd <- phy$edge[k, 1L]
    cnd <- phy$edge[k, 2L]
    v <- numeric(p)
    v[idx_of[[as.character(pnd)]]] <- 1
    if (cnd > n_tip) v[idx_of[[as.character(cnd)]]] <- -1
    add_ge(v, eps)
  }

  Amat <- do.call(cbind, c(eq_cols, ge_cols))
  bvec <- c(eq_b, ge_b)
  Dmat <- diag(2, p)
  dvec <- 2 * node_age_init[int_nodes]

  sol <- quadprog::solve.QP(
    Dmat = Dmat,
    dvec = dvec,
    Amat = Amat,
    bvec = bvec,
    meq = length(eq_b)
  )

  out <- node_age_init
  out[int_nodes] <- sol$solution
  out
}

#' Run RelTime, then project the point estimates onto a full set of merged
#' node bounds before computing analytical CIs.
#'
#' @param phy            Rooted phylogram.
#' @param calibration_df Pairwise calibration table with taxonA,taxonB,age_min,age_max.
#' @param root_age       Optional root point estimate used for the initial
#'                       unconstrained RelTime profile. Defaults to the midpoint
#'                       of the merged root interval, or to the exact root age.
#' @param n_sites        Alignment length for CI variance calculations.
#' @param eps            Minimum internal branch duration enforced in the QP.
#' @return               List with tree, ci, node_age, bounds, and initial ages.
run_reltime_with_bounds_ci <- function(phy, calibration_df, root_age = NULL,
                                       n_sites = 1000L, eps = 1e-6) {
  bounds <- reltime_merge_calibration_bounds(phy, calibration_df)
  n_tip <- Ntip(phy)
  root_node <- n_tip + 1L

  root_row <- bounds[bounds$node == root_node, , drop = FALSE]
  if (!nrow(root_row) && (is.null(root_age) || !is.finite(root_age))) {
    stop("A root calibration or explicit root_age is required.")
  }

  root_var <- 0
  if (nrow(root_row)) {
    if (is.null(root_age) || !is.finite(root_age)) {
      if (is.finite(root_row$age_max[1])) {
        root_age <- mean(c(root_row$age_min[1], root_row$age_max[1]))
      } else {
        root_age <- root_row$age_min[1]
      }
    }
    if (is.finite(root_row$age_max[1]) &&
        abs(root_row$age_max[1] - root_row$age_min[1]) > 1e-10) {
      root_var <- ((root_row$age_max[1] - root_row$age_min[1]) / (2 * 1.96))^2
    }
  }
  if (!is.finite(root_age) || root_age <= 0) {
    stop("root_age must be positive after root calibration handling.")
  }

  base <- .compute_H(phy)
  node_age_init <- .assign_node_ages(phy, base$H, root_age)

  lower <- numeric(length(node_age_init))
  upper <- rep(Inf, length(node_age_init))
  if (nrow(bounds)) {
    lower[bounds$node] <- bounds$age_min
    upper[bounds$node] <- bounds$age_max
  } else {
    lower[root_node] <- root_age
    upper[root_node] <- root_age
  }

  node_age <- .reltime_project_node_ages_qp(
    phy = phy,
    node_age_init = node_age_init,
    lower = lower,
    upper = upper,
    eps = eps
  )

  dated <- phy
  for (k in seq_len(nrow(phy$edge))) {
    p <- phy$edge[k, 1L]
    c <- phy$edge[k, 2L]
    T_p <- node_age[p]
    T_c <- if (c <= n_tip) 0 else node_age[c]
    dated$edge.length[k] <- max(T_p - T_c, eps)
  }

  cal_min <- as.list(stats::setNames(lower[bounds$node], bounds$node))
  cal_max <- as.list(stats::setNames(upper[bounds$node], bounds$node))
  ci <- reltime_ci(
    phy = phy,
    node_age = node_age,
    n_sites = n_sites,
    root_var = root_var,
    cal_min = cal_min,
    cal_max = cal_max
  )

  list(
    tree = dated,
    ci = ci,
    node_age = node_age,
    bounds = bounds,
    initial_node_age = node_age_init
  )
}

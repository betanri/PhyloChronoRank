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

# ---- null-coalesce operator (safe for R < 4.1) ----------------------------
if (!exists("%||%", mode = "function"))
  `%||%` <- function(a, b) if (!is.null(a)) a else b

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
# ---- Arithmetic-mean H (used for CI variance computation only) -------------
.compute_H_arithmetic <- function(phy) {
  n_tip <- Ntip(phy)
  total <- n_tip + phy$Nnode
  H_sum  <- numeric(total)
  ntip_s <- integer(total)
  ntip_s[seq_len(n_tip)] <- 1L
  e  <- phy$edge
  el <- phy$edge.length
  for (k in postorder(phy)) {
    p <- e[k, 1L]; c <- e[k, 2L]
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

# Keep .compute_H as alias for backward compat (CI code calls it)
.compute_H <- .compute_H_arithmetic

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

# ---- core: MEGA-style node height computation and age assignment -----------
# Reimplements MEGA's AnchorClockTree algorithm (Tamura et al. 2012):
#   1. SetHa: compute ha0 (anchor height) for each node post-order
#   2. ComputeHeightOfNode: geometric-mean height formula with rate-ratio cap
#   3. RecomputeNodeHeightDown: propagate height adjustments to descendants
#   4. Global time factor + calibration clamping
#
# This replaces the old arithmetic-mean .assign_node_ages which did not match
# MEGA's output for trees with many calibrations.

.assign_node_ages_mega <- function(phy, lower, upper, root_age,
                                   max_rate_ratio = 1e10,
                                   skip_recompute_down = FALSE) {
  FP_CUTOFF <- 1e-12
  n_tip     <- Ntip(phy)
  total     <- n_tip + phy$Nnode
  root_node <- n_tip + 1L
  e  <- phy$edge
  el <- phy$edge.length

  # Build children list and parent lookup
  ch <- vector("list", total)
  parent_of <- integer(total)
  for (k in seq_len(nrow(e))) {
    p <- e[k, 1L]; c <- e[k, 2L]
    ch[[p]] <- c(ch[[p]], c)
    parent_of[c] <- p
  }

  # Branch length lookup: blen[c] = branch length from parent to c
  blen <- numeric(total)
  for (k in seq_len(nrow(e))) blen[e[k, 2L]] <- el[k]

  # Identify anchored (calibrated) nodes
  anchored <- logical(total)
  for (nd in seq.int(n_tip + 1L, total)) {
    if (lower[nd] > FP_CUTOFF || upper[nd] > FP_CUTOFF)
      anchored[nd] <- TRUE
  }

  # Compute minh/maxh (will be set after time factor)
  minh <- numeric(total)
  maxh <- rep(Inf, total)

  # ---- Step 1: SetHa (post-order) ------------------------------------------
  # ha0[n] = max(anchored_child_heights) from below
  # For tips: ha0 = 0 (h0)
  # For internal: ha0 = max(ha0_or_minh of children)
  ha0 <- numeric(total)  # tips = 0
  for (k in postorder(phy)) {
    p <- e[k, 1L]; c <- e[k, 2L]
    if (c <= n_tip) {
      child_h <- 0  # tip
    } else if (anchored[c]) {
      child_h <- minh[c]  # will be 0 initially, updated after time factor
    } else {
      child_h <- ha0[c]
    }
    ha0[p] <- max(ha0[p], child_h)
  }

  # ---- Step 2: Initial height computation (post-order, no calibration) ------
  # Simple RelTime: height[p] from branch lengths
  # First pass: compute raw heights using T[c] = T[p] * H[c] / (b + H[c])
  # Actually MEGA computes height bottom-up using the geometric mean formula.
  # height[nd] = distance from nd to tips (in relative time units)

  height <- numeric(total)
  rate <- rep(1, total)

  # Post-order: compute initial heights
  # For tips: height = 0
  # For internal nodes with children c1, c2:
  #   h01 = ha0[c1] (or minh[c1] if anchored)
  #   h02 = ha0[c2] (or minh[c2] if anchored)
  #   h1 = height[c1] + blen[c1]/rate[c1] - h01
  #   h2 = height[c2] + blen[c2]/rate[c2] - h02
  #   delta = |h01 - h02|
  #   h = (-delta + sqrt(delta^2 + 4*h1*h2)) / 2 + max(h01, h02)

  # BFS post-order
  post_visit <- integer(0)
  queue <- root_node
  while (length(queue)) {
    nd <- queue[1L]; queue <- queue[-1L]
    post_visit <- c(post_visit, nd)
    kids <- ch[[nd]]
    for (kid in kids) if (kid > n_tip) queue <- c(queue, kid)
  }
  post_visit <- rev(post_visit)  # post-order

  for (nd in post_visit) {
    kids <- ch[[nd]]
    if (is.null(kids) || length(kids) == 0) next
    c1 <- kids[1]; c2 <- kids[2]

    h01 <- if (c1 > n_tip && anchored[c1]) minh[c1] else if (c1 > n_tip) ha0[c1] else 0
    h02 <- if (c2 > n_tip && anchored[c2]) minh[c2] else if (c2 > n_tip) ha0[c2] else 0

    if (c1 <= n_tip) {
      h1 <- blen[c1]  # tip: height=0, rate=1
    } else if (rate[c1] < FP_CUTOFF) {
      h1 <- height[c1] - h01
    } else {
      h1 <- height[c1] + blen[c1] / rate[c1] - h01
    }

    if (c2 <= n_tip) {
      h2 <- blen[c2]
    } else if (rate[c2] < FP_CUTOFF) {
      h2 <- height[c2] - h02
    } else {
      h2 <- height[c2] + blen[c2] / rate[c2] - h02
    }

    h1 <- max(h1, 0)
    h2 <- max(h2, 0)

    if (h1 > FP_CUTOFF && h2 > FP_CUTOFF &&
        max(h1, h2) / min(h1, h2) < max_rate_ratio) {
      # Geometric mean formula
      delta <- abs(h01 - h02)
      h <- (-delta + sqrt(delta * delta + 4 * h1 * h2)) / 2 + max(h01, h02)
    } else {
      # Rate-ratio exceeded: constrained calculation
      if (h1 > h2) {
        h <- max(h1 / max_rate_ratio + h01, h02)
      } else if (h1 < h2) {
        h <- max(h2 / max_rate_ratio + h02, h01)
      } else {
        h <- max(h1 / max_rate_ratio + h01, h2 / max_rate_ratio + h02)
      }
    }

    height[nd] <- h

    # Compute rate as geometric mean of children's rates
    r1 <- if (c1 <= n_tip) 1 else rate[c1]
    r2 <- if (c2 <= n_tip) 1 else rate[c2]
    rate[nd] <- sqrt(r1 * r2)
  }

  # ---- Step 3: Global time factor ------------------------------------------
  # minRate[i] = lower[i] / height[i], maxRate[i] = upper[i] / height[i]
  # FTimeFactor = midpoint of feasible range across all anchored nodes
  r0 <- 0; r1 <- Inf
  for (nd in seq.int(n_tip + 1L, total)) {
    if (!anchored[nd] || height[nd] < FP_CUTOFF) next
    if (lower[nd] > FP_CUTOFF)
      r0 <- max(r0, lower[nd] / height[nd])
    if (is.finite(upper[nd]) && upper[nd] > FP_CUTOFF)
      r1 <- min(r1, upper[nd] / height[nd])
  }

  if (!is.finite(r1) || r1 < r0) {
    # Constraints conflict — use root age
    time_factor <- root_age / height[root_node]
  } else {
    time_factor <- (r0 + r1) / 2
  }

  # Set minh/maxh as calibration / time_factor
  for (nd in seq.int(n_tip + 1L, total)) {
    if (lower[nd] > FP_CUTOFF)
      minh[nd] <- lower[nd] / time_factor
    if (is.finite(upper[nd]) && upper[nd] > FP_CUTOFF)
      maxh[nd] <- upper[nd] / time_factor
  }

  # ---- Step 4: Recompute with calibration clamping -------------------------
  # Re-run SetHa with updated minh
  ha0 <- numeric(total)
  for (k in postorder(phy)) {
    p <- e[k, 1L]; c <- e[k, 2L]
    if (c <= n_tip) {
      child_h <- 0
    } else if (anchored[c]) {
      child_h <- minh[c]
    } else {
      child_h <- ha0[c]
    }
    ha0[p] <- max(ha0[p], child_h)
  }

  # Re-run ComputeHeightOfNode with anchored heights
  # Initialize from Step 2 heights (MEGA carries forward, doesn't start from zero)
  height2 <- height
  rate2 <- rate

  .recompute_down <- function(nd, dh) {
    if (abs(dh) < FP_CUTOFF) return()
    p <- parent_of[nd]
    if (p == 0) return()
    denom <- height2[p] - dh - ha0[nd]
    if (abs(denom) < FP_CUTOFF) return()
    dh_new <- (height2[nd] - ha0[nd]) / denom * dh
    # Clamp to bounds
    if (maxh[nd] < Inf && height2[nd] + dh_new > maxh[nd])
      dh_new <- maxh[nd] - height2[nd]
    if (minh[nd] > FP_CUTOFF && height2[nd] + dh_new < minh[nd])
      dh_new <- minh[nd] - height2[nd]
    if (abs(dh_new) < FP_CUTOFF) return()
    old_h <- height2[nd]
    height2[nd] <- height2[nd] + dh_new
    if (old_h > FP_CUTOFF) rate2[nd] <- rate2[nd] * old_h / height2[nd]
    kids <- ch[[nd]]
    if (!is.null(kids)) {
      for (kid in kids) if (kid > n_tip) .recompute_down(kid, dh_new)
    }
  }

  for (nd in post_visit) {
    kids <- ch[[nd]]
    if (is.null(kids) || length(kids) == 0) next
    c1 <- kids[1]; c2 <- kids[2]

    h01 <- if (c1 > n_tip && anchored[c1]) minh[c1] else if (c1 > n_tip) ha0[c1] else 0
    h02 <- if (c2 > n_tip && anchored[c2]) minh[c2] else if (c2 > n_tip) ha0[c2] else 0

    if (c1 <= n_tip) {
      h1 <- blen[c1]
    } else if (rate2[c1] < FP_CUTOFF) {
      h1 <- height2[c1] - h01
    } else {
      h1 <- height2[c1] + blen[c1] / rate2[c1] - h01
    }

    if (c2 <= n_tip) {
      h2 <- blen[c2]
    } else if (rate2[c2] < FP_CUTOFF) {
      h2 <- height2[c2] - h02
    } else {
      h2 <- height2[c2] + blen[c2] / rate2[c2] - h02
    }

    h1 <- max(h1, 0); h2 <- max(h2, 0)

    if (h1 > FP_CUTOFF && h2 > FP_CUTOFF &&
        max(h1, h2) / min(h1, h2) < max_rate_ratio) {
      delta <- abs(h01 - h02)
      h <- (-delta + sqrt(delta * delta + 4 * h1 * h2)) / 2 + max(h01, h02)
    } else {
      if (h1 > h2) h <- max(h1 / max_rate_ratio + h01, h02)
      else if (h1 < h2) h <- max(h2 / max_rate_ratio + h02, h01)
      else h <- max(h1 / max_rate_ratio + h01, h2 / max_rate_ratio + h02)
    }

    h0_old <- height2[nd]
    height2[nd] <- h

    r1c <- if (c1 <= n_tip) 1 else rate2[c1]
    r2c <- if (c2 <= n_tip) 1 else rate2[c2]
    rate2[nd] <- sqrt(r1c * r2c)

    if (!skip_recompute_down && abs(height2[nd] - h0_old) > FP_CUTOFF) {
      if (c1 > n_tip) .recompute_down(c1, height2[nd] - (h1 + h01))
      if (c2 > n_tip) .recompute_down(c2, height2[nd] - (h2 + h02))
    }

    # Clamp to calibration bounds
    if (maxh[nd] < Inf && height2[nd] > maxh[nd]) {
      old_h <- height2[nd]
      height2[nd] <- maxh[nd]
      if (!skip_recompute_down) {
        dh <- height2[nd] - old_h
        if (c1 > n_tip) .recompute_down(c1, dh)
        if (c2 > n_tip) .recompute_down(c2, dh)
      }
      if (old_h > FP_CUTOFF) rate2[nd] <- rate2[nd] * old_h / height2[nd]
    } else if (minh[nd] > FP_CUTOFF && height2[nd] < minh[nd]) {
      old_h <- height2[nd]
      height2[nd] <- minh[nd]
      if (!skip_recompute_down) {
        dh <- height2[nd] - old_h
        if (c1 > n_tip) .recompute_down(c1, dh)
        if (c2 > n_tip) .recompute_down(c2, dh)
      }
      if (old_h > FP_CUTOFF) rate2[nd] <- rate2[nd] * old_h / height2[nd]
    }
  }

  # ---- Step 5: Convert heights to ages using time factor --------------------
  node_age <- numeric(total)
  for (nd in seq.int(n_tip + 1L, total)) {
    node_age[nd] <- height2[nd] * time_factor
  }
  # Root age = root_age (or calibrated)
  node_age[root_node] <- root_age

  node_age
}

# Backward-compatible wrapper
.assign_node_ages <- function(phy, H, root_age, cal_ages = NULL) {
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

.reltime_prepare_bounds_context <- function(phy, calibration_df, root_age = NULL,
                                             use_densities = FALSE) {
  bounds <- reltime_merge_calibration_bounds(phy, calibration_df)
  n_tip <- Ntip(phy)
  total <- n_tip + phy$Nnode
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

  lower <- numeric(total)
  upper <- rep(Inf, total)
  if (nrow(bounds)) {
    lower[bounds$node] <- bounds$age_min
    upper[bounds$node] <- bounds$age_max
  } else {
    lower[root_node] <- root_age
    upper[root_node] <- root_age
  }

  # ---- Build density objects when requested --------------------------------
  densities <- NULL
  if (isTRUE(use_densities) && nrow(bounds)) {
    densities <- vector("list", total)
    for (i in seq_len(nrow(bounds))) {
      nd <- bounds$node[i]
      mn <- bounds$age_min[i]
      mx <- bounds$age_max[i]
      # Check if calibration_df has distribution info
      has_dist <- "dist_type" %in% names(calibration_df)
      if (has_dist) {
        # Find matching calibration row(s) for this node
        # Use the first match with a dist_type
        cal_rows <- which(
          calibration_df$taxonA %in% phy$tip.label &
          calibration_df$taxonB %in% phy$tip.label
        )
        matched <- FALSE
        for (ci in cal_rows) {
          mrca <- ape::getMRCA(phy, c(
            as.character(calibration_df$taxonA[ci]),
            as.character(calibration_df$taxonB[ci])
          ))
          if (!is.null(mrca) && mrca == nd &&
              !is.na(calibration_df$dist_type[ci]) &&
              nchar(as.character(calibration_df$dist_type[ci])) > 0) {
            dtype <- tolower(as.character(calibration_df$dist_type[ci]))
            if (dtype == "normal") {
              densities[[nd]] <- list(
                type = "normal",
                mean = as.numeric(calibration_df$dist_mean[ci]),
                stddev = as.numeric(calibration_df$dist_stddev[ci])
              )
            } else if (dtype == "lognormal") {
              densities[[nd]] <- list(
                type = "lognormal",
                offset = as.numeric(calibration_df$dist_offset[ci]),
                mean = as.numeric(calibration_df$dist_mean[ci]),
                stddev = as.numeric(calibration_df$dist_stddev[ci])
              )
            } else if (dtype == "exponential") {
              densities[[nd]] <- list(
                type = "exponential",
                time = as.numeric(calibration_df$dist_time[ci]),
                lambda = as.numeric(calibration_df$dist_lambda[ci])
              )
            } else {
              # uniform or unknown → use bounds
              densities[[nd]] <- list(type = "uniform", min = mn, max = mx)
            }
            matched <- TRUE
            break
          }
        }
        if (!matched) {
          # Default: convert hard bounds to normal density
          # Mode = midpoint, stddev so that bounds ≈ ±2σ
          densities[[nd]] <- .bounds_to_normal(mn, mx)
        }
      } else {
        # No dist_type column: auto-convert min/max to normal densities
        densities[[nd]] <- .bounds_to_normal(mn, mx)
      }
    }
    # Root always gets a tight density
    if (is.null(densities[[root_node]])) {
      densities[[root_node]] <- list(
        type = "normal",
        mean = root_age,
        stddev = max(abs(root_age) * 0.01, 1e-6)
      )
    }
  }

  list(
    bounds = bounds,
    root_age = root_age,
    root_var = root_var,
    lower = lower,
    upper = upper,
    densities = densities
  )
}

# Convert hard min/max bounds to a calibration density.
# For fixed calibrations (min == max): use normal with stddev = 5% of age.
#   This allows the density approach to treat fixed points as "soft targets"
#   rather than hard constraints, giving backbone branches room to breathe.
# For interval calibrations: normal with mode = midpoint, stddev = (max-min)/4
#   so that the 95% interval ≈ [min, max].
# For min-only calibrations: exponential with rate lambda = 1/min, offset = min.
.bounds_to_normal <- function(mn, mx) {
  if (!is.finite(mx)) {
    # Min-only calibration → exponential: peaks at offset, decays
    list(type = "exponential", time = mn, lambda = 1 / max(mn, 1e-6))
  } else if (abs(mx - mn) < 1e-10) {
    # Fixed calibration point → normal with 5% stddev (soft target)
    mid <- mean(c(mn, mx))
    list(type = "normal", mean = mid, stddev = max(mid * 0.05, 1e-6))
  } else {
    mid <- mean(c(mn, mx))
    sd <- (mx - mn) / 4
    list(type = "normal", mean = mid, stddev = sd)
  }
}

.run_reltime_with_bounds_context <- function(phy, ctx, eps = 1e-6,
                                              use_densities = FALSE,
                                              smooth_backbone = FALSE) {
  n_tip <- Ntip(phy)

  base <- .compute_H(phy)
  node_age_init <- .assign_node_ages(phy, base$H, ctx$root_age)

  if (isTRUE(use_densities) && !is.null(ctx$densities)) {
    node_age <- .reltime_project_node_ages_density(
      phy = phy,
      node_age_init = node_age_init,
      densities = ctx$densities,
      eps = eps
    )
  } else {
    node_age <- .reltime_project_node_ages_local(
      phy = phy,
      node_age_init = node_age_init,
      lower = ctx$lower,
      upper = ctx$upper,
      eps = eps
    )
  }

  # Optional: smooth near-zero backbone branches
  if (isTRUE(smooth_backbone)) {
    node_age <- .smooth_near_zero_branches(
      node_age = node_age,
      phy = phy,
      lower = ctx$lower,
      upper = ctx$upper,
      eps = eps
    )
  }

  dated <- phy
  for (k in seq_len(nrow(phy$edge))) {
    p <- phy$edge[k, 1L]
    c <- phy$edge[k, 2L]
    T_p <- node_age[p]
    T_c <- if (c <= n_tip) 0 else node_age[c]
    dated$edge.length[k] <- max(T_p - T_c, eps)
  }

  list(
    tree = dated,
    node_age = node_age,
    bounds = ctx$bounds,
    initial_node_age = node_age_init,
    lower = ctx$lower,
    upper = ctx$upper,
    root_age = ctx$root_age,
    root_var = ctx$root_var
  )
}

#' Run bounded RelTime without confidence intervals
#'
#' @param phy            Rooted phylogram.
#' @param calibration_df Pairwise calibration table with taxonA,taxonB,age_min,age_max.
#' @param root_age       Optional root point estimate used for the initial profile.
#' @param eps            Minimum internal branch duration enforced in the QP.
#' @return               List with tree, node_age, bounds, and projection context.
run_reltime_with_bounds <- function(phy, calibration_df, root_age = NULL,
                                    eps = 1e-6, use_densities = FALSE,
                                    smooth_backbone = FALSE) {
  ctx <- .reltime_prepare_bounds_context(phy, calibration_df,
                                          root_age = root_age,
                                          use_densities = use_densities)
  .run_reltime_with_bounds_context(phy, ctx, eps = eps,
                                    use_densities = use_densities,
                                    smooth_backbone = smooth_backbone)
}

#' Compute Tao-style analytical CIs for a bounded RelTime run
#'
#' @param phy       Rooted phylogram used for the bounded RelTime run.
#' @param rel_run   Output list from [run_reltime_with_bounds()] or
#'                  [reltime_bootstrap_ci()].
#' @param n_sites   Alignment length for CI variance calculations.
#' @return          data.frame(node, age, ci_lo, ci_hi, se)
#' Generic Tao-style analytical CIs for any dated tree
#'
#' Applies the Tao et al. (2020) delta-method variance to any phylogram +
#' chronogram pair, regardless of the dating method used.  Works for chronos,
#' treePL, RelTime, MCMCTree, etc.
#'
#' @param phy          Rooted phylogram (substitution branch lengths).
#' @param dated_tree   Corresponding dated tree (time branch lengths).
#' @param n_sites      Alignment length for Poisson branch-length variance.
#' @param node_bounds  Optional data.frame with columns node, age_min, age_max
#'                     for calibration-bound truncation of CIs.
#' @return             data.frame(node, age, ci_lo, ci_hi, se)
tao_analytical_ci <- function(phy, dated_tree, n_sites = 1000L,
                              node_bounds = NULL) {
  n_tip  <- Ntip(dated_tree)
  n_node <- dated_tree$Nnode

  ## Build node_age vector indexed by node number
  node_age <- numeric(n_tip + n_node)
  bt <- ape::branching.times(dated_tree)
  for (nm in names(bt)) node_age[as.integer(nm)] <- bt[nm]

  ## Build calibration bound lists
  cal_min <- cal_max <- NULL
  if (!is.null(node_bounds) && nrow(node_bounds) > 0L) {
    cal_min <- as.list(stats::setNames(node_bounds$age_min, node_bounds$node))
    cal_max <- as.list(stats::setNames(node_bounds$age_max, node_bounds$node))
  }

  reltime_ci(
    phy      = phy,
    node_age = node_age,
    n_sites  = n_sites,
    root_var = 0,
    cal_min  = cal_min,
    cal_max  = cal_max
  )
}

reltime_tao_ci_from_bounds_run <- function(phy, rel_run, n_sites = 1000L) {
  bounds <- rel_run$bounds
  cal_min <- if (nrow(bounds)) {
    as.list(stats::setNames(rel_run$lower[bounds$node], bounds$node))
  } else {
    list()
  }
  cal_max <- if (nrow(bounds)) {
    as.list(stats::setNames(rel_run$upper[bounds$node], bounds$node))
  } else {
    list()
  }
  reltime_ci(
    phy = phy,
    node_age = rel_run$node_age,
    n_sites = n_sites,
    root_var = rel_run$root_var,
    cal_min = cal_min,
    cal_max = cal_max
  )
}

.reltime_poisson_perturb_phy <- function(phy, n_sites, min_edge = 1e-12) {
  if (!is.finite(n_sites) || n_sites <= 0) {
    stop("n_sites must be a positive finite number for RelTime bootstrap.")
  }
  out <- phy
  lambda <- pmax(as.numeric(phy$edge.length) * n_sites, 0)
  out$edge.length <- stats::rpois(length(lambda), lambda) / n_sites
  out$edge.length <- pmax(out$edge.length, min_edge)
  out
}

#' Bootstrap confidence intervals for bounded RelTime point trees
#'
#' This is the shared empirical uncertainty path used in PCR comparisons.
#' Branch lengths are perturbed with a Poisson parametric bootstrap, the same
#' bounded RelTime point-dating procedure is rerun on each replicate, and node
#' ages are summarized with empirical quantiles. Tao-style analytical CIs remain
#' available as a separate supplemental RelTime-specific diagnostic.
#'
#' @param phy            Rooted phylogram.
#' @param calibration_df Pairwise calibration table with taxonA,taxonB,age_min,age_max.
#' @param root_age       Optional root point estimate used for the initial profile.
#' @param B              Number of bootstrap replicates.
#' @param n_sites        Alignment length used for Poisson branch-length resampling.
#' @param eps            Minimum internal branch duration enforced in the QP.
#' @param quiet          Suppress progress messages.
#' @param trees          Return successful bootstrap trees as well.
#' @param min_edge       Small positive floor for zero-count bootstrap branches.
#' @return               List with tree, ci, node_age, bounds, and bootstrap metadata.
reltime_bootstrap_ci <- function(phy, calibration_df, root_age = NULL,
                                 B = 100L, n_sites = 1000L, eps = 1e-6,
                                 quiet = FALSE, trees = FALSE,
                                 min_edge = 1e-12, use_densities = FALSE,
                                 smooth_backbone = TRUE) {
  B <- as.integer(B)
  if (!is.finite(B) || B < 1L) stop("B must be a positive integer.")
  if (!is.finite(n_sites) || n_sites <= 0) {
    stop("n_sites must be a positive finite number for RelTime bootstrap.")
  }

  ctx <- .reltime_prepare_bounds_context(phy, calibration_df, root_age = root_age,
                                          use_densities = use_densities)
  base_run <- .run_reltime_with_bounds_context(phy, ctx, eps = eps,
                                                use_densities = use_densities,
                                                smooth_backbone = smooth_backbone)

  n_tip <- Ntip(phy)
  int_nodes <- seq.int(n_tip + 1L, n_tip + phy$Nnode)
  boot_age <- matrix(NA_real_, nrow = length(int_nodes), ncol = B)
  boot_trees <- if (isTRUE(trees)) vector("list", B) else NULL
  ok <- logical(B)

  for (i in seq_len(B)) {
    if (!quiet) cat("\rRunning RelTime bootstrap:", i, "/", B)
    phy_boot <- try(.reltime_poisson_perturb_phy(phy, n_sites, min_edge = min_edge), silent = TRUE)
    if (inherits(phy_boot, "try-error")) next
    run_i <- try(.run_reltime_with_bounds_context(phy_boot, ctx, eps = eps,
                                                    use_densities = use_densities,
                                                    smooth_backbone = smooth_backbone), silent = TRUE)
    if (inherits(run_i, "try-error") || is.null(run_i)) next
    ages_i <- run_i$node_age[int_nodes]
    if (!all(is.finite(ages_i))) next
    boot_age[, i] <- ages_i
    ok[i] <- TRUE
    if (isTRUE(trees)) boot_trees[[i]] <- run_i$tree
  }
  if (!quiet) cat("\n")
  if (!any(ok)) stop("RelTime bootstrap produced no successful replicates")

  boot_ok <- boot_age[, ok, drop = FALSE]
  qfun <- function(probs) {
    apply(
      boot_ok,
      1L,
      stats::quantile,
      probs = probs,
      names = FALSE,
      na.rm = TRUE,
      type = 8
    )
  }

  ci_df <- data.frame(
    node = int_nodes,
    age = base_run$node_age[int_nodes],
    ci_lo = qfun(0.025),
    ci_hi = qfun(0.975),
    q50_lo = qfun(0.25),
    q50_hi = qfun(0.75),
    se = NA_real_,
    stringsAsFactors = FALSE
  )

  base_run$ci <- ci_df
  base_run$bootstrap_ok <- ok
  base_run$n_bootstrap_ok <- sum(ok)
  if (isTRUE(trees)) base_run$trees <- boot_trees[ok]
  base_run
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

# ---- calibration density helpers -------------------------------------------
# Supported distributions: uniform, normal, lognormal, exponential.
# Each function returns log-density (unnormalised is fine for optimisation).

.cal_density_log <- function(t, dist) {
  type <- tolower(dist$type)
  if (type == "uniform") {
    lo <- dist$min %||% 0
    hi <- dist$max %||% Inf
    if (t >= lo && t <= hi) return(0)  # log(1/(hi-lo)) is constant
    return(-Inf)
  }
  if (type == "normal") {
    mu <- dist$mean
    sd <- dist$stddev
    return(-0.5 * ((t - mu) / sd)^2)
  }
  if (type == "lognormal") {
    # offset: hard minimum; distribution is on (t - offset)
    off <- dist$offset %||% 0
    if (t <= off) return(-Inf)
    x <- t - off
    mu <- dist$mean
    sd <- dist$stddev
    return(-log(x) - 0.5 * ((log(x) - mu) / sd)^2)
  }
  if (type == "exponential") {
    # time = hard minimum (offset), lambda = rate
    off <- dist$time %||% 0
    lam <- dist$lambda
    if (t < off) return(-Inf)
    return(-lam * (t - off))
  }
  # fallback: uniform-like (no penalty)
  0
}

# Mode of each density (used as the "target" for initial attraction)
.cal_density_mode <- function(dist) {
  type <- tolower(dist$type)
  if (type == "uniform") {
    lo <- dist$min %||% 0
    hi <- dist$max %||% lo
    return(mean(c(lo, hi)))
  }
  if (type == "normal") return(dist$mean)
  if (type == "lognormal") {
    off <- dist$offset %||% 0
    mu <- dist$mean
    sd <- dist$stddev
    return(off + exp(mu - sd^2))
  }
  if (type == "exponential") return(dist$time %||% 0)
  NA_real_
}

# Effective min/max bounds from density (where density > exp(-9) of mode)
.cal_density_bounds <- function(dist) {
  type <- tolower(dist$type)
  if (type == "uniform") {
    return(c(dist$min %||% 0, dist$max %||% Inf))
  }
  if (type == "normal") {
    # ~4.2 sigma covers exp(-9) of mode
    return(c(max(0, dist$mean - 4.2 * dist$stddev),
             dist$mean + 4.2 * dist$stddev))
  }
  if (type == "lognormal") {
    off <- dist$offset %||% 0
    mu <- dist$mean
    sd <- dist$stddev
    # Approximate bounds
    lo <- off + exp(mu - 4 * sd)
    hi <- off + exp(mu + 4 * sd)
    return(c(lo, hi))
  }
  if (type == "exponential") {
    off <- dist$time %||% 0
    lam <- dist$lambda
    return(c(off, off + 9 / lam))
  }
  c(0, Inf)
}

# ---- project node ages with calibration densities --------------------------
# Extension of the MEGA-style local rescaling that replaces hard clamping
# with density-weighted attraction.
#
# For calibrated nodes, instead of snapping to min/max, the algorithm finds
# the age within the feasible parent–child interval that maximises the
# calibration density. Uncalibrated nodes are rescaled proportionally
# (preserving RelTime branch structure) when a calibrated ancestor shifts.
#
# Steps:
#   1. Propagate hard bounds from densities for hierarchical consistency
#   2. Global scaling factor from density modes
#   3. Pre-order density-guided placement + proportional subtree rescale
#   4. Iterative refinement (coordinate descent on log-density)
#   5. Enforce parent > child ordering
.reltime_project_node_ages_density <- function(phy, node_age_init, densities,
                                                eps = 1e-6, n_iter = 5L) {
  n_tip <- Ntip(phy)
  total <- n_tip + phy$Nnode
  root_node <- n_tip + 1L
  ch <- .reltime_children(phy)

  # Build parent lookup
  par <- integer(total)
  e <- phy$edge
  for (k in seq_len(nrow(e))) par[e[k, 2L]] <- e[k, 1L]

  # ---- BFS pre-order of internal nodes ------------------------------------
  visit <- integer(0)
  queue <- root_node
  while (length(queue)) {
    nd <- queue[1L]; queue <- queue[-1L]
    visit <- c(visit, nd)
    kids <- ch[[nd]]
    for (kid in kids) if (kid > n_tip) queue <- c(queue, kid)
  }

  # ---- Derive hard bounds from densities ----------------------------------
  lower <- numeric(total)
  upper <- rep(Inf, total)
  for (nd in visit) {
    if (!is.null(densities[[nd]])) {
      bds <- .cal_density_bounds(densities[[nd]])
      lower[nd] <- bds[1]
      upper[nd] <- bds[2]
    }
  }

  # ---- Step 1: propagate bounds for hierarchical consistency ---------------
  for (nd in visit) {
    if (!is.finite(upper[nd])) next
    for (kid in ch[[nd]])
      if (kid > n_tip) upper[kid] <- min(upper[kid], upper[nd])
  }
  for (nd in rev(visit)) {
    for (kid in ch[[nd]])
      if (kid > n_tip && is.finite(lower[kid]) && lower[kid] > 0)
        lower[nd] <- max(lower[nd], lower[kid])
  }

  # ---- Step 2: global scaling from density modes --------------------------
  node_age <- node_age_init
  r0 <- 0; r1 <- Inf
  for (nd in visit) {
    if (node_age_init[nd] < 1e-12 || is.null(densities[[nd]])) next
    bds <- .cal_density_bounds(densities[[nd]])
    if (bds[1] > 0) r0 <- max(r0, bds[1] / node_age_init[nd])
    if (is.finite(bds[2])) r1 <- min(r1, bds[2] / node_age_init[nd])
  }
  if (is.finite(r1) && r1 >= r0 && r0 > 0) {
    mode_ratios <- numeric(0)
    for (nd in visit) {
      if (node_age_init[nd] < 1e-12 || is.null(densities[[nd]])) next
      mr <- .cal_density_mode(densities[[nd]]) / node_age_init[nd]
      if (is.finite(mr) && mr >= r0 && mr <= r1)
        mode_ratios <- c(mode_ratios, mr)
    }
    f <- if (length(mode_ratios)) median(mode_ratios) else (r0 + r1) / 2
    f <- max(r0, min(r1, f))
    for (nd in visit) node_age[nd] <- node_age_init[nd] * f
  }

  # Helper: find feasible range for node nd
  .feasible_range <- function(nd) {
    parent_age <- if (par[nd] > 0) node_age[par[nd]] else Inf
    max_child_age <- 0
    for (kid in ch[[nd]])
      if (kid > n_tip) max_child_age <- max(max_child_age, node_age[kid])
    lo <- max(lower[nd], max_child_age + eps)
    hi <- min(upper[nd], parent_age - eps)
    c(lo, hi)
  }

  # Helper: rescale subtree below nd
  .rescale_subtree <- function(nd, sf, skip_calibrated = FALSE) {
    stack <- integer(0)
    for (kid in ch[[nd]]) if (kid > n_tip) stack <- c(stack, kid)
    while (length(stack)) {
      cur <- stack[1L]; stack <- stack[-1L]
      if (skip_calibrated && !is.null(densities[[cur]])) next
      node_age[cur] <<- node_age[cur] * sf
      for (kid in ch[[cur]]) if (kid > n_tip) stack <- c(stack, kid)
    }
  }

  # ---- Step 3: pre-order density-guided placement -------------------------
  for (nd in visit) {
    if (is.null(densities[[nd]])) next
    old_age <- node_age[nd]
    rng <- .feasible_range(nd)
    lo <- rng[1]; hi <- rng[2]
    if (lo > hi) { lo <- lower[nd]; hi <- upper[nd] }
    if (lo > hi) next

    new_age <- .golden_section_max(
      function(t) .cal_density_log(t, densities[[nd]]),
      lo, hi
    )

    if (abs(new_age - old_age) > 1e-12 && old_age > 1e-12) {
      sf <- new_age / old_age
      node_age[nd] <- new_age
      .rescale_subtree(nd, sf)
    }
  }

  # ---- Step 4: iterative refinement (coordinate descent) ------------------
  for (iter in seq_len(n_iter)) {
    moved <- FALSE
    for (nd in visit) {
      if (is.null(densities[[nd]])) next
      old_age <- node_age[nd]
      rng <- .feasible_range(nd)
      if (rng[1] > rng[2]) next

      new_age <- .golden_section_max(
        function(t) .cal_density_log(t, densities[[nd]]),
        rng[1], rng[2]
      )

      if (abs(new_age - old_age) > eps * 0.1 && old_age > 1e-12) {
        sf <- new_age / old_age
        node_age[nd] <- new_age
        .rescale_subtree(nd, sf, skip_calibrated = TRUE)
        moved <- TRUE
      }
    }
    if (!moved) break
  }

  # ---- Step 5: enforce parent > child + eps --------------------------------
  for (nd in visit) {
    for (kid in ch[[nd]]) {
      if (kid > n_tip && node_age[kid] >= node_age[nd] - eps)
        node_age[kid] <- node_age[nd] - eps
    }
  }

  node_age
}

# Golden-section search for maximum of f on [a, b]
.golden_section_max <- function(f, a, b, tol = 1e-8, max_iter = 100L) {
  gr <- (sqrt(5) + 1) / 2
  c <- b - (b - a) / gr
  d <- a + (b - a) / gr
  for (i in seq_len(max_iter)) {
    if (abs(b - a) < tol) break
    if (f(c) < f(d)) {
      a <- c
    } else {
      b <- d
    }
    c <- b - (b - a) / gr
    d <- a + (b - a) / gr
  }
  (a + b) / 2
}

# ---- post-processing: smooth near-zero backbone branches ------------------
# After calibration projection, some backbone branches collapse because an
# uncalibrated parent node gets pushed right on top of its calibrated child.
# This happens when RelTime's initial age for the parent is BELOW the child's
# calibrated age — the parent is forced upward to satisfy parent > child, but
# lands barely above the child.
#
# Fix: for each near-zero internal branch where the parent is uncalibrated,
# place the parent at the proportional midpoint between its calibrated child
# and its nearest calibrated ancestor.  This preserves the tree's temporal
# structure while giving backbone branches real duration.
.smooth_near_zero_branches <- function(node_age, phy, lower, upper,
                                        eps = 1e-6, threshold = 1e-5) {
  n_tip <- Ntip(phy)
  total <- n_tip + phy$Nnode
  root_node <- n_tip + 1L
  ch <- .reltime_children(phy)
  e <- phy$edge

  # Build parent lookup
  par_of <- integer(total)
  for (k in seq_len(nrow(e))) par_of[e[k, 2L]] <- e[k, 1L]

  # Identify constrained nodes. Interval-bounded nodes also need protection
  # here; otherwise the smoothing step can undo a valid lower/upper clamp
  # applied during the main projection.
  is_calibrated <- logical(total)
  for (nd in seq_len(total)) {
    if (nd > n_tip &&
        ((is.finite(lower[nd]) && lower[nd] > 0) ||
         (is.finite(upper[nd]) && upper[nd] < Inf)))
      is_calibrated[nd] <- TRUE
  }
  # Root is always anchored
  is_calibrated[root_node] <- TRUE

  # Multiple passes to handle cascading near-zero chains
  for (pass in 1:3) {
    any_fixed <- FALSE

    for (k in seq_len(nrow(e))) {
      child <- e[k, 2L]
      if (child <= n_tip) next
      parent <- e[k, 1L]
      branch_dur <- node_age[parent] - node_age[child]
      if (branch_dur > threshold) next      # not near-zero
      if (is_calibrated[parent]) next        # can't move calibrated parent

      # Find nearest calibrated ancestor above parent
      anc <- par_of[parent]
      while (anc > 0 && !is_calibrated[anc]) anc <- par_of[anc]
      anc_age <- if (anc > 0) node_age[anc] else node_age[root_node]

      # Find the calibrated child age (floor)
      cal_child_age <- node_age[child]

      # Place parent proportionally between ancestor and child
      # Use the initial RelTime proportion if available
      room <- anc_age - cal_child_age
      if (room < 2 * eps) next

      # Position at the midpoint between ancestor and child, but
      # respect other children's ages (parent must be > all children)
      max_child_age <- 0
      for (kid in ch[[parent]])
        if (kid > n_tip) max_child_age <- max(max_child_age, node_age[kid])

      # New parent age: midpoint of available range
      new_lo <- max_child_age + eps
      new_hi <- anc_age - eps
      if (is.finite(upper[parent])) new_hi <- min(new_hi, upper[parent])
      if (is.finite(lower[parent]) && lower[parent] > 0)
        new_lo <- max(new_lo, lower[parent])

      if (new_lo >= new_hi) next

      new_parent_age <- (new_lo + new_hi) / 2

      # Ensure the branch to child opens up meaningfully
      if (new_parent_age - cal_child_age < threshold) next

      node_age[parent] <- new_parent_age
      any_fixed <- TRUE
    }

    if (!any_fixed) break
  }

  # Final pass: enforce parent > child + eps (might be violated by moves above)
  # BFS pre-order
  visit <- integer(0)
  queue <- root_node
  while (length(queue)) {
    nd <- queue[1L]; queue <- queue[-1L]
    visit <- c(visit, nd)
    kids <- ch[[nd]]
    for (kid in kids) if (kid > n_tip) queue <- c(queue, kid)
  }
  for (nd in visit) {
    for (kid in ch[[nd]]) {
      if (kid > n_tip && node_age[kid] >= node_age[nd] - eps)
        node_age[kid] <- node_age[nd] - eps
    }
  }

  node_age
}

# ---- project node ages onto full merged calibration bounds -----------------
# MEGA-style local rescaling (Tamura et al. 2013, 2018; Tao et al. 2020):
#   1. Propagate bounds hierarchically for consistency
#   2. Try a single global scaling factor f (midpoint of feasible range)
#   3. Pre-order clamp-and-proportional-rescale for remaining violations
#   4. Enforce parent > child ordering
#
# Unlike the previous QP formulation, each calibration adjustment only
# rescales the subtree below the clamped node, preserving relative branch
# proportions and avoiding the artificial near-zero branches that QP
# introduced at backbone nodes.
.reltime_project_node_ages_local <- function(phy, node_age_init, lower, upper, eps = 1e-6) {
  n_tip <- Ntip(phy)
  total <- n_tip + phy$Nnode
  root_node <- n_tip + 1L
  ch <- .reltime_children(phy)
  parent_of <- integer(total)
  for (k in seq_len(nrow(phy$edge))) parent_of[phy$edge[k, 2L]] <- phy$edge[k, 1L]

  # ---- BFS pre-order of internal nodes ------------------------------------
  visit <- integer(0)
  queue <- root_node
  while (length(queue)) {
    nd <- queue[1L]; queue <- queue[-1L]
    visit <- c(visit, nd)
    kids <- ch[[nd]]
    for (kid in kids) if (kid > n_tip) queue <- c(queue, kid)
  }

  # ---- MEGA-style local rescaling (Tamura et al. 2012) --------------------
  # 1. Start with RelTime relative ages scaled so root = root_age.
  # 2. For each calibrated node (pre-order, root first): set its age to the
  #    calibrated value, then proportionally rescale all descendant ages
  #    that lie between this node and the next calibrated descendant(s).
  #    The rescaling maps the interval [0, old_age] -> [0, new_age] for
  #    descendants, preserving relative proportions within each segment.
  # 3. Iterate until stable (calibrations can interact when nested).

  node_age <- node_age_init

  # Identify fixed-point calibrated nodes
  is_cal <- logical(total)
  cal_age <- rep(NA_real_, total)
  for (nd in visit) {
    if (is.finite(lower[nd]) && is.finite(upper[nd]) &&
        lower[nd] > 0 && abs(upper[nd] - lower[nd]) < 1e-10) {
      is_cal[nd] <- TRUE
      cal_age[nd] <- lower[nd]
    }
  }

  max_iter <- 50L
  for (iter in seq_len(max_iter)) {
    max_change <- 0

    # Pre-order: at each calibrated node, fix age and rescale below
    for (nd in visit) {
      if (!is_cal[nd]) next

      target <- cal_age[nd]
      old_age <- node_age[nd]
      if (abs(old_age) < 1e-15) {
        node_age[nd] <- target
        max_change <- max(max_change, abs(target))
        next
      }

      sf <- target / old_age
      if (abs(sf - 1) < 1e-12) next

      max_change <- max(max_change, abs(target - old_age))
      node_age[nd] <- target

      # Rescale all descendants until hitting another calibrated node
      stack <- integer(0)
      for (kid in ch[[nd]]) if (kid > n_tip) stack <- c(stack, kid)
      while (length(stack)) {
        cur <- stack[1L]; stack <- stack[-1L]
        if (is_cal[cur]) next  # stop at next calibrated node
        node_age[cur] <- node_age[cur] * sf
        for (kid in ch[[cur]]) if (kid > n_tip) stack <- c(stack, kid)
      }
    }

    if (max_change < 1e-10) break
  }

  # ---- Feasibility clamping for non-fixed nodes ----------------------------
  # After fixed-point calibration enforcement, nodes with interval bounds still
  # need to be clamped explicitly. The previous implementation only respected
  # surrounding fixed calibrations, which meant min/max interval bounds could
  # silently be ignored. Here we intersect the local feasible interval implied
  # by calibrated ancestors/descendants with the node's own [lower, upper]
  # bounds and clamp into that range.

  for (nd in visit) {
    if (is_cal[nd]) next
    if (nd == root_node) next

    # Find nearest calibrated ancestor (walk up)
    cal_anc_age <- node_age[root_node]  # default: root
    p <- parent_of[nd]
    while (p > 0) {
      if (is_cal[p]) { cal_anc_age <- cal_age[p]; break }
      p <- parent_of[p]
    }

    # Find nearest calibrated descendant (walk down) — take the max age
    # among all calibrated descendants reachable without crossing another cal
    cal_desc_age <- 0  # default: tips at age 0
    stack <- integer(0)
    for (kid in ch[[nd]]) if (kid > n_tip) stack <- c(stack, kid)
    while (length(stack)) {
      cur <- stack[1L]; stack <- stack[-1L]
      if (is_cal[cur]) {
        cal_desc_age <- max(cal_desc_age, cal_age[cur])
      } else {
        for (kid in ch[[cur]]) if (kid > n_tip) stack <- c(stack, kid)
      }
    }

    # Intersect ancestor/descendant feasibility with the node's own bounds.
    lo <- cal_desc_age + eps
    hi <- cal_anc_age - eps
    if (is.finite(lower[nd]) && lower[nd] > 0) lo <- max(lo, lower[nd])
    if (is.finite(upper[nd])) hi <- min(hi, upper[nd])

    # If the interval collapses, honor the bound target as closely as possible
    # and let the parent>child pass below clean up any tiny ordering issues.
    if (lo > hi) {
      if (is.finite(lower[nd]) && lower[nd] > 0 && lower[nd] >= cal_desc_age) {
        node_age[nd] <- lower[nd]
      } else if (is.finite(upper[nd]) && upper[nd] <= cal_anc_age) {
        node_age[nd] <- upper[nd]
      } else {
        node_age[nd] <- max(min(node_age[nd], cal_anc_age - eps), cal_desc_age + eps)
      }
      next
    }

    if (node_age[nd] < lo) node_age[nd] <- lo
    if (node_age[nd] > hi) node_age[nd] <- hi
  }

  # ---- Enforce parent > child + eps without breaking lower bounds ----------
  # First push ancestors upward so constrained descendants cannot be dragged
  # below their own minima by a too-young parent.
  for (nd in rev(visit)) {
    if (nd == root_node) next
    p <- parent_of[nd]
    if (p <= 0) next
    needed_parent_age <- node_age[nd] + eps
    if (is.finite(lower[p]) && lower[p] > 0) {
      needed_parent_age <- max(needed_parent_age, lower[p])
    }
    if (node_age[p] < needed_parent_age) {
      node_age[p] <- needed_parent_age
    }
  }

  # Final cleanup: only lower a child if doing so does not violate the child's
  # own lower bound; otherwise raise the parent just enough to open the branch.
  for (nd in visit) {
    for (kid in ch[[nd]]) {
      if (kid <= n_tip) next
      if (node_age[kid] < node_age[nd] - eps) next
      kid_lower <- if (is.finite(lower[kid]) && lower[kid] > 0) lower[kid] else -Inf
      target_child <- node_age[nd] - eps
      if (target_child >= kid_lower) {
        node_age[kid] <- target_child
      } else {
        node_age[nd] <- kid_lower + eps
      }
    }
  }

  node_age
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
                                       n_sites = 1000L, eps = 1e-6,
                                       use_densities = FALSE,
                                       smooth_backbone = TRUE) {
  rel_run <- run_reltime_with_bounds(
    phy = phy,
    calibration_df = calibration_df,
    root_age = root_age,
    eps = eps,
    use_densities = use_densities,
    smooth_backbone = smooth_backbone
  )
  rel_run$ci <- reltime_tao_ci_from_bounds_run(
    phy = phy,
    rel_run = rel_run,
    n_sites = n_sites
  )
  rel_run
}

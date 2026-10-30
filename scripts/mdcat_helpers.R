##############################################################
## mdcat_helpers.R
##
## MD-Cat (Mai, Charvel & Mirarab, "Expectation-Maximization enables
## Phylogenetic Dating under a Categorical Rate Model";
## github.com/uym2/MD-Cat) dates a phylogram by maximum likelihood.
## Branch rates are drawn from k rate categories (default 50) whose
## probabilities are estimated by EM, so the rate distribution is not
## forced to be lognormal, autocorrelated or clock-like.
##
## MD-Cat takes calibrations only as fixed node ages (hard equality
## constraints), not as bounds. The shared node bounds are fitted with an
## active set:
##   1. bounds with min == max (e.g. congruified ages) are fixed as given;
##      if there are none, the interval with the smallest relative width
##      is fixed at its midpoint, so the tree has a time scale
##   2. MD-Cat runs; every other calibrated node whose age falls outside
##      its bounds is fixed at the bound it crossed
##   3. repeat 2 until all calibrated nodes are inside their bounds (or
##      `rounds` runs are used); then the start midpoint, if any, is
##      released and the fit repeated, so no arbitrary midpoint stays fixed
## At the optimum a crossed bound is binding, so this is the usual way to
## add inequality constraints to an equality-only fitter. Nodes that never
## cross a bound stay free.
##
## Fixed ages must be in order (ancestor older than descendant), with room
## for MD-Cat's minimum branch duration (1e-3 time units); they are put in
## order inside their bounds with a gap of 2e-3 per branch before each run.
## Point ages may move by at most 1e-4 of the root age for this.
##
## MD-Cat itself reports no intervals. mdcat_bootstrap_ci adds them with
## the Poisson branch-length bootstrap used for chronos, treePL and RelTime
## in this pipeline (see the comment above mdcat_bootstrap_rep).
##
## Install: python3 -m venv venv; venv/bin/pip install git+https://github.com/uym2/MD-Cat
## and pass venv/bin/md_cat.py as --mdcat-bin (or set MDCAT_BIN).
##############################################################

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

mdcat_resolve_bin <- function(user_bin = NULL) {
  candidates <- c(user_bin %||% NA_character_, Sys.getenv("MDCAT_BIN", unset = ""),
                  Sys.which("md_cat.py"))
  candidates <- unique(trimws(candidates[!is.na(candidates)]))
  for (cand in candidates[nzchar(candidates)]) {
    if (file.exists(cand)) return(normalizePath(cand, winslash = "/", mustWork = TRUE))
  }
  ""
}

## Put fixed ages in order (ancestor older than descendant by at least one
## gap per branch between them), each inside [lo, hi].
.mdcat_order_ages <- function(phy, nodes, age, lo, hi, gap) {
  parent <- integer(ape::Ntip(phy) + phy$Nnode)
  parent[phy$edge[, 2L]] <- phy$edge[, 1L]
  anc <- rep(NA_integer_, length(nodes))
  steps <- integer(length(nodes))
  depth <- integer(length(nodes))
  for (i in seq_along(nodes)) {
    a <- parent[nodes[i]]
    k <- 1L
    while (a > 0L) {
      j <- match(a, nodes)
      if (!is.na(j) && is.na(anc[i])) { anc[i] <- j; steps[i] <- k }
      depth[i] <- depth[i] + 1L
      a <- parent[a]
      k <- k + 1L
    }
  }
  for (i in order(-depth)) {
    j <- anc[i]
    if (!is.na(j) && age[j] < age[i] + gap * steps[i]) age[j] <- min(age[i] + gap * steps[i], hi[j])
  }
  for (i in order(depth)) {
    j <- anc[i]
    if (!is.na(j) && age[i] > age[j] - gap * steps[i]) age[i] <- max(age[j] - gap * steps[i], lo[i])
  }
  bad <- !is.na(anc) & age[anc] - age < gap * steps * (1 - 1e-9)
  if (any(bad)) stop("MD-Cat: could not order the fixed calibration ages inside their bounds at node(s) ",
                     paste(nodes[bad], collapse = ", "))
  age
}

## MD-Cat writes branch lengths with 4 decimals, so tips are not exactly
## level. Each terminal branch is reset so every tip sits at age 0.
.mdcat_level_tips <- function(tr) {
  n_tip <- ape::Ntip(tr)
  depth <- ape::node.depth.edgelength(tr)
  height <- max(depth[seq_len(n_tip)])
  term <- tr$edge[, 2L] <= n_tip
  tr$edge.length[term] <- height - depth[tr$edge[term, 1L]]
  tr
}

.mdcat_node_ages <- function(tr) {
  d <- ape::node.depth.edgelength(tr)
  max(d[seq_len(ape::Ntip(tr))]) - d
}

## One md_cat.py call. Returns the levelled tree, log-likelihood, solver and
## the largest calibration drift MD-Cat reported (see mdcat_run).
.mdcat_call <- function(mdcat_bin, in_tree, cal_file, stem, reps, ncat, seq_len, seed, threads) {
  args <- c("-i", in_tree, "-t", cal_file, "-b", "-o", paste0(stem, "_raw.tre"),
            "-p", as.integer(reps), "-k", as.integer(ncat), "-l", as.integer(seq_len),
            "--randSeed", as.integer(seed), "--threads", as.integer(threads))
  log_file <- paste0(stem, ".log")
  status <- suppressWarnings(system2(mdcat_bin, args, stdout = log_file, stderr = log_file))
  log_txt <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character(0)
  if (!file.exists(paste0(stem, "_raw.tre"))) {
    return(list(ok = FALSE, error = paste0("exit ", status, ": ", paste(tail(log_txt, 3L), collapse = " | "))))
  }
  txt <- paste(readLines(paste0(stem, "_raw.tre"), warn = FALSE), collapse = "")
  tr <- ape::read.tree(text = gsub("\\[[^]]*\\]", "", txt))
  tr$node.label <- NULL
  llh <- suppressWarnings(as.numeric(sub("^Best log-likelihood:\\s*", "",
                                         grep("^Best log-likelihood:", log_txt, value = TRUE))))
  solver <- unique(sub("^Solver \\[optimization\\]: (\\S+).*", "\\1",
                       grep("^Solver \\[optimization\\]: .* is in use", log_txt, value = TRUE)))
  drift <- suppressWarnings(as.numeric(sub(".*Violate by\\s*", "",
                                           grep("^Inconsistent divergence time", log_txt, value = TRUE))))
  list(ok = TRUE, tree = .mdcat_level_tips(tr), llh = if (length(llh)) llh[1L] else NA_real_,
       solver = paste(solver, collapse = ";"), max_drift = if (length(drift)) max(drift) else 0,
       log_file = log_file)
}

## Run MD-Cat on `phy` with the shared node bounds (active set, see top).
## `jobs` > 1 runs the `reps` random starts as parallel single-start
## processes (seeds seed + 1..reps) and keeps the best log-likelihood, which
## is what MD-Cat's own -p does in one process.
## `start_calibrations` (a _calibrations_used.csv table of an earlier run)
## starts the active set from the ages that run ended with fixed.
mdcat_run <- function(mdcat_bin, phy, node_bounds, work_dir, prefix,
                      reps = 10L, ncat = 50L, seq_len = 1000L, seed = 1L,
                      threads = 1L, jobs = 1L, rounds = 5L, quiet = FALSE,
                      start_calibrations = NULL) {
  bad <- grepl("[[:space:]+=():;,\\[\\]]", phy$tip.label)
  if (any(bad)) stop("MD-Cat cannot read tip labels with spaces or + = ( ) : ; , [ ]: ",
                     paste(head(phy$tip.label[bad], 3L), collapse = ", "))
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  f <- function(suffix) file.path(work_dir, paste0(prefix, "_MDCat", suffix))

  ## Every internal node gets a unique label. Without a MOSEK license MD-Cat
  ## solves with OSQP, whose tolerance can exceed MD-Cat's own consistency
  ## check (5e-4); MD-Cat then prints a warning that crashes on unlabelled
  ## nodes. With labels the warning is harmless, and the calibrated ages are
  ## checked below.
  in_tree <- phy
  in_tree$node.label <- paste0("mdcatN", seq_len(in_tree$Nnode))
  ape::write.tree(in_tree, f("_input.tre"))

  cal <- data.frame(node = node_bounds$node, taxonA = node_bounds$taxonA, taxonB = node_bounds$taxonB,
                    age_min = node_bounds$age_min, age_max = node_bounds$age_max, stringsAsFactors = FALSE)
  cal$age_max[is.na(cal$age_max)] <- Inf
  ## Nodes are found again from the taxa, so a tree read back from file
  ## (where ape may number nodes differently) still matches.
  cal$node <- unname(mapply(function(a, b) ape::getMRCA(phy, c(a, b)), cal$taxonA, cal$taxonB))
  scale <- max(cal$age_min, cal$age_max[is.finite(cal$age_max)])
  ## MD-Cat writes branch lengths with 4 decimals; summed over a root-to-tip
  ## path this moves node ages by a few 1e-4 time units, so the bounds check
  ## allows at least 1e-3 (it matters for trees in units of 100 Myr).
  tol <- max(1e-4 * scale, 1e-3)
  point <- abs(cal$age_max - cal$age_min) <= 1e-9
  cal$fixed <- ifelse(point, "point", "")
  cal$fixed_age <- ifelse(point, cal$age_min, NA_real_)
  if (!any(point)) {
    width <- ifelse(is.finite(cal$age_max), (cal$age_max - cal$age_min) / ((cal$age_max + cal$age_min) / 2), Inf)
    if (!any(is.finite(width))) stop("MD-Cat needs at least one calibration with a finite age (point or interval).")
    i <- which.min(width)
    cal$fixed[i] <- "start_midpoint"
    cal$fixed_age[i] <- (cal$age_min[i] + cal$age_max[i]) / 2
  }
  if (!is.null(start_calibrations)) {
    k <- match(paste(cal$taxonA, cal$taxonB), paste(start_calibrations$taxonA, start_calibrations$taxonB))
    prev <- ifelse(is.na(start_calibrations$fixed[k]), "", start_calibrations$fixed[k])
    if (!anyNA(k) && any(nzchar(prev))) {
      cal$fixed <- prev
      cal$fixed_age <- ifelse(nzchar(prev), start_calibrations$fixed_age[k], NA_real_)
    }
  }

  started <- Sys.time()
  fit <- NULL
  for (r in seq_len(rounds)) {
    use <- which(nzchar(cal$fixed))
    lo <- ifelse(cal$fixed[use] == "point", cal$fixed_age[use] - tol, cal$age_min[use])
    hi <- ifelse(cal$fixed[use] == "point", cal$fixed_age[use] + tol, cal$age_max[use])
    ## MD-Cat forces every branch to last at least 1e-3 time units, so fixed
    ## ages need at least that much room per branch between them.
    cal$fixed_age[use] <- .mdcat_order_ages(phy, cal$node[use], cal$fixed_age[use], lo, hi,
                                            max(2e-3, 1e-6 * scale))
    cal_file <- f(sprintf("_round%d_calibrations.txt", r))
    writeLines(sprintf("%s+%s %.10g", cal$taxonA[use], cal$taxonB[use], cal$fixed_age[use]), cal_file)
    stems <- f(if (jobs > 1L) sprintf("_round%d_start%02d", r, seq_len(reps)) else sprintf("_round%d", r))
    one <- function(k) .mdcat_call(mdcat_bin, f("_input.tre"), cal_file, stems[k],
                                   reps = if (jobs > 1L) 1L else reps, ncat = ncat, seq_len = seq_len,
                                   seed = if (jobs > 1L) seed + k else seed, threads = threads)
    res <- if (jobs > 1L && .Platform$OS.type == "unix") {
      parallel::mclapply(seq_along(stems), one, mc.cores = min(jobs, length(stems)), mc.preschedule = FALSE)
    } else {
      lapply(seq_along(stems), one)
    }
    ok <- vapply(res, function(x) is.list(x) && isTRUE(x$ok), logical(1))
    if (!any(ok)) stop("MD-Cat wrote no tree in round ", r, ": ", res[[1L]]$error %||% "worker crashed")
    llh <- vapply(res, function(x) if (is.list(x) && isTRUE(x$ok)) x$llh else -Inf, numeric(1))
    fit <- res[[which.max(llh)]]
    fit$n_starts_ok <- sum(ok)

    age <- .mdcat_node_ages(fit$tree)
    nd <- mapply(function(a, b) ape::getMRCA(fit$tree, c(a, b)), cal$taxonA, cal$taxonB)
    cal$mdcat_age <- age[nd]
    free <- !nzchar(cal$fixed)
    low <- free & cal$mdcat_age < cal$age_min - tol
    high <- free & cal$mdcat_age > cal$age_max + tol
    if (!quiet) cat(sprintf("    MD-Cat round %d: %d fixed, log-likelihood %.3f, %d free node(s) outside bounds\n",
                            r, length(use), fit$llh, sum(low | high)))
    if (!any(low | high)) {
      ## The start midpoint only gave the first run a time scale. Once other
      ## ages are fixed it is released; if it then crosses a bound it is
      ## fixed there in the next round.
      start <- which(cal$fixed == "start_midpoint")
      if (length(start) && sum(nzchar(cal$fixed)) > 1L && r < rounds) {
        cal$fixed[start] <- ""
        cal$fixed_age[start] <- NA_real_
        next
      }
      break
    }
    cal$fixed[low] <- sprintf("min_round%d", r)
    cal$fixed_age[low] <- cal$age_min[low]
    cal$fixed[high] <- sprintf("max_round%d", r)
    cal$fixed_age[high] <- cal$age_max[high]
  }
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  cal$within_bounds <- cal$mdcat_age >= cal$age_min - tol & cal$mdcat_age <= cal$age_max + tol
  write.csv(cal, f("_calibrations_used.csv"), row.names = FALSE)
  ape::write.tree(fit$tree, f(".tre"))
  list(tree = fit$tree, tree_file = f(".tre"), calibrations = cal, log_file = fit$log_file,
       seconds = seconds, log_likelihood = fit$llh, solver = fit$solver, rounds = r,
       n_starts_ok = fit$n_starts_ok, max_drift = fit$max_drift)
}

## ---- Bootstrap CIs ----
## Each replicate resamples every branch length as Poisson(n_sites * length)
## / n_sites, the scheme used for the chronos, treePL and RelTime intervals
## in this pipeline, and dates the result with one random start. The active
## set starts from the calibrations the point run ended with (`calibrations`,
## its _calibrations_used.csv table); bounds a replicate crosses are added as
## in mdcat_run, and a replicate still outside its bounds is replaced. Nodes
## the point run fixed at a bound stay fixed, so their intervals have zero
## width, like point calibrations. Replicate trees are kept in
## work_dir/rep_NNN, so a stopped run resumes where it ended.
mdcat_bootstrap_rep <- function(i, mdcat_bin, phy, calibrations, work_dir, prefix,
                                n_sites = 1000L, ncat = 50L, seq_len = 1000L, seed = 1L,
                                threads = 1L, rounds = 5L) {
  rep_dir <- file.path(work_dir, sprintf("rep_%03d", i))
  out <- file.path(rep_dir, paste0(prefix, "_MDCat.tre"))
  if (file.exists(out)) {
    tr <- try(ape::read.tree(out), silent = TRUE)
    if (inherits(tr, "phylo") && ape::Ntip(tr) == ape::Ntip(phy)) return(out)
  }
  set.seed((as.integer(seed) * 1000L + as.integer(i)) %% .Machine$integer.max)
  phy_b <- phy
  phy_b$edge.length <- pmax(stats::rpois(length(phy$edge.length), pmax(phy$edge.length, 0) * n_sites) / n_sites,
                            1e-12)
  res <- try(mdcat_run(mdcat_bin, phy_b, calibrations, rep_dir, prefix, reps = 1L, ncat = ncat,
                       seq_len = seq_len, seed = seed + i, threads = threads, jobs = 1L, rounds = rounds,
                       quiet = TRUE, start_calibrations = calibrations), silent = TRUE)
  if (inherits(res, "try-error") || any(!res$calibrations$within_bounds)) {
    unlink(out)
    return(NA_character_)
  }
  unlink(list.files(rep_dir, pattern = "(_input|_raw)\\.tre$", full.names = TRUE))
  out
}

## Node ages of `tr` for the clades in `sig` (sorted tip labels per node).
.mdcat_clade_signatures <- function(tr) {
  parts <- ape::prop.part(tr)
  labs <- attr(parts, "labels")
  vapply(parts, function(p) paste(sort(labs[p]), collapse = "\r"), character(1))
}

.mdcat_ages_by_clade <- function(tr, sig) {
  age <- .mdcat_node_ages(tr)
  own <- .mdcat_clade_signatures(tr)
  age[ape::Ntip(tr) + match(sig, own)]
}

## Percentile intervals (2.5-97.5 and 25-75) over B replicates, on the
## nodes of `dated_tree` (the MD-Cat point estimate).
mdcat_bootstrap_ci <- function(mdcat_bin, phy, calibrations, dated_tree, work_dir, prefix,
                               B = 100L, n_sites = 1000L, ncat = 50L, seq_len = 1000L, seed = 1L,
                               threads = 1L, jobs = 1L, rounds = 5L) {
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  one <- function(i) mdcat_bootstrap_rep(i, mdcat_bin, phy, calibrations, work_dir, prefix,
                                         n_sites = n_sites, ncat = ncat, seq_len = seq_len,
                                         seed = seed, threads = threads, rounds = rounds)
  run <- function(ids) {
    res <- if (jobs > 1L && .Platform$OS.type == "unix") {
      parallel::mclapply(ids, one, mc.cores = jobs, mc.preschedule = FALSE)
    } else {
      lapply(ids, one)
    }
    vapply(res, function(x) if (is.character(x) && length(x) == 1L) x else NA_character_, "")
  }
  trees <- run(seq_len(B))
  ## failed replicates are replaced, up to B extra attempts
  tried <- B
  while (sum(!is.na(trees)) < B && tried < 2L * B) {
    ids <- seq.int(tried + 1L, min(2L * B, tried + B - sum(!is.na(trees))))
    trees <- c(trees, run(ids))
    tried <- max(ids)
  }
  trees <- trees[!is.na(trees)]
  if (length(trees) < B) stop("MD-Cat bootstrap: only ", length(trees), " of ", tried, " replicates finished.")
  trees <- trees[seq_len(B)]

  sig <- .mdcat_clade_signatures(dated_tree)
  nodes <- ape::Ntip(dated_tree) + seq_along(sig)
  boot_age <- vapply(trees, function(f) .mdcat_ages_by_clade(ape::read.tree(f), sig), numeric(length(sig)))
  if (anyNA(boot_age)) stop("MD-Cat bootstrap: replicate trees do not match the point-estimate topology.")
  q <- function(p) apply(boot_age, 1L, stats::quantile, probs = p, names = FALSE, type = 8)
  ci <- data.frame(node = nodes, age = .mdcat_node_ages(dated_tree)[nodes],
                   ci_lo = q(0.025), ci_hi = q(0.975), q50_lo = q(0.25), q50_hi = q(0.75),
                   se = NA_real_, stringsAsFactors = FALSE)
  list(ci = ci, n_bootstrap_ok = B, n_bootstrap_attempts = tried)
}

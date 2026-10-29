##############################################################
## treepl_helpers.R
##
## Repo-local treePL support (Smith & O'Meara 2012). The dating
## engine runs treePL in the standard order:
##
##   1. prime  — treePL tests optimizers and prints the best
##               opt/optad/optcvad/moredetail lines, then exits.
##               Prime is stochastic, so it is repeated with
##               different seeds and the most frequent result is
##               used (the repeated-prime idea is from
##               github.com/tongjial/treepl_wrapper).
##   2. cv     — cross-validation over a geometric smoothing series
##               (cvstart -> cvstop by cvmultstep); the smoothing
##               with the lowest CV chi-square is selected. If the
##               best value is at the edge of the series, the
##               series is extended outward.
##   3. final  — thorough runs with the selected smoothing under
##               several optimizer settings (each distinct prime
##               result, treePL defaults, extra seeds). The run with
##               the lowest penalized-likelihood objective is kept,
##               and the node-age spread across settings is reported.
##
## treePL's penalized-likelihood objective is not comparable across
## smoothing values (the penalty grows with smoothing), so it is
## never used to pick smoothing. At one fixed smoothing it is
## comparable, so it is used to pick among optimizer settings.
##
## treePL seeds its random number generator from the clock unless the
## config has a `seed` line, so every config written here has one.
##
## Replicated random CV (median of accepted replicates, with outlier
## rejection), the full 144-setting optimizer screen with a deep seed search,
## and the bounds check of every dated tree follow A. Pyron's treePL benchmark
## for the VGP time tree (github.com/VGP/vgp-dating,
## dating/treePL/treepl_benchmark.py).
##
## Bootstrap CIs: branch lengths are perturbed with the same Poisson
## resampling scheme used for chronos and RelTime, treePL is rerun
## on each perturbed phylogram with the selected smoothing and primed
## lines, and node ages are summarized with empirical quantiles after
## matching internal nodes by descendant tip sets.
##############################################################

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

.treepl_fmt_num <- function(x) {
  out <- formatC(as.numeric(x), format = "f", digits = 12)
  sub("\\.?0+$", "", out)
}

treepl_resolve_bin <- function(user_bin, repo_dir) {
  candidates <- c(
    user_bin %||% NA_character_,
    Sys.getenv("TREEPL_BIN", unset = ""),
    Sys.which("treePL"),
    file.path(dirname(repo_dir), "treePL")
  )
  candidates <- unique(trimws(candidates))
  candidates <- candidates[nzchar(candidates)]
  for (cand in candidates) {
    if (file.exists(cand)) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }
  ""
}

treepl_write_cfg <- function(cfg_path, tree_path, out_tree_path, smooth, numsites,
                             node_bounds, thorough = TRUE, prime = FALSE,
                             opt = NULL, plsimaniter = NULL, pliter = NULL,
                             extra_lines = NULL, seed = NULL) {
  lines <- c(
    paste0("treefile = ", tree_path),
    paste0("outfile = ", out_tree_path),
    paste0("numsites = ", as.integer(numsites)),
    paste0("smooth = ", .treepl_fmt_num(smooth))
  )
  if (!is.null(seed) && is.finite(seed)) lines <- c(lines, paste0("seed = ", as.integer(seed)))
  if (isTRUE(thorough)) lines <- c(lines, "thorough")
  if (isTRUE(prime)) lines <- c(lines, "prime")
  if (!is.null(opt) && nzchar(opt)) {
    lines <- c(lines, paste0("opt = ", as.integer(opt)))
  }
  if (!is.null(plsimaniter) && nzchar(plsimaniter)) {
    lines <- c(lines, paste0("plsimaniter = ", as.integer(plsimaniter)))
  }
  if (!is.null(pliter) && nzchar(pliter)) {
    lines <- c(lines, paste0("pliter = ", as.integer(pliter)))
  }
  if (length(extra_lines)) {
    extra_lines <- trimws(as.character(extra_lines))
    extra_lines <- extra_lines[nzchar(extra_lines)]
    lines <- c(lines, extra_lines)
  }
  for (i in seq_len(nrow(node_bounds))) {
    tag <- paste0("N", node_bounds$node[i])
    lines <- c(
      lines,
      paste0("mrca = ", tag, " ", node_bounds$taxonA[i], " ", node_bounds$taxonB[i]),
      paste0("min = ", tag, " ", .treepl_fmt_num(node_bounds$age_min[i]))
    )
    if (is.finite(node_bounds$age_max[i])) {
      lines <- c(lines, paste0("max = ", tag, " ", .treepl_fmt_num(node_bounds$age_max[i])))
    }
  }
  writeLines(lines, cfg_path)
}

treepl_run_one <- function(treepl_bin, cfg_path, log_path, omp_threads = 1L) {
  if (file.exists(log_path)) unlink(log_path)
  suppressWarnings(system2(
    treepl_bin,
    args = cfg_path,
    stdout = log_path,
    stderr = log_path,
    env = paste0("OMP_NUM_THREADS=", as.integer(omp_threads))
  ))
}

treepl_parse_prime_lines <- function(log_path) {
  if (!file.exists(log_path)) return(character(0))
  txt <- readLines(log_path, warn = FALSE)
  marker <- grep("PLACE THE LINES BELOW IN THE CONFIGURATION FILE", txt, fixed = TRUE)
  if (!length(marker)) return(character(0))
  out <- trimws(txt[(marker[1L] + 1L):length(txt)])
  out <- out[
    grepl("^(opt(ad|cvad)?\\s*=|moredetail(ad|cvad)?\\b)", out) &
      !grepl("^[-]+$", out)
  ]
  unique(out[nzchar(out)])
}

## Step 1. Prime on the tree that will be dated. treePL exits after priming
## and writes an empty outfile, so the outfile is a throwaway temp path.
treepl_run_prime <- function(treepl_bin, tree_path, node_bounds, smooth, numsites,
                             cfg_path, log_path, thorough = TRUE,
                             opt = NULL, plsimaniter = NULL, pliter = NULL,
                             omp_threads = 1L, seed = NULL) {
  junk_out <- tempfile(pattern = "treepl_prime_", fileext = ".tre")
  on.exit(unlink(c(junk_out, paste0(junk_out, ".r8s")), force = TRUE), add = TRUE)
  treepl_write_cfg(
    cfg_path, tree_path, junk_out, smooth, numsites, node_bounds,
    thorough = thorough, prime = TRUE,
    opt = opt, plsimaniter = plsimaniter, pliter = pliter, seed = seed
  )
  treepl_run_one(treepl_bin, cfg_path, log_path, omp_threads = omp_threads)
  lines <- treepl_parse_prime_lines(log_path)
  if (!length(lines)) {
    stop("treePL prime step did not report optimizer settings; see ", log_path)
  }
  lines
}

## Repeated prime. Replicate i uses seed + i; cfg_path/log_path are functions
## of i. Returns the distinct outcomes ordered by frequency (ties keep the
## order of first appearance); outcomes$lines[[1]] is the modal result.
treepl_prime_consensus <- function(treepl_bin, tree_path, node_bounds, smooth, numsites,
                                   reps, seed, cfg_path, log_path, ...) {
  reps <- max(1L, as.integer(reps))
  runs <- vector("list", reps)
  for (i in seq_len(reps)) {
    res <- try(treepl_run_prime(treepl_bin, tree_path, node_bounds, smooth, numsites,
                                cfg_path(i), log_path(i), seed = seed + i, ...),
               silent = TRUE)
    runs[[i]] <- if (inherits(res, "try-error")) character(0) else res
  }
  ok <- lengths(runs) > 0L
  if (!any(ok)) stop("all ", reps, " treePL prime replicates failed; see ", log_path(1L))
  keys <- vapply(runs, paste, character(1), collapse = "; ")
  uniq <- unique(keys[ok])
  counts <- vapply(uniq, function(k) sum(keys == k), integer(1))
  ord <- order(-counts, match(uniq, keys))
  list(
    outcomes = data.frame(
      lines_key = uniq[ord], n = unname(counts[ord]),
      stringsAsFactors = FALSE
    ),
    lines = runs[match(uniq[ord], keys)],
    reps = data.frame(rep = seq_len(reps), seed = seed + seq_len(reps),
                      ok = ok, lines_key = ifelse(ok, keys, NA_character_),
                      stringsAsFactors = FALSE)
  )
}

## treePL CV walks a geometric series, so the smoothing grid must be one.
## A single-value series (cvstart == cvstop) returns junk scores in treePL,
## so at least two values are required. cvstop is set half a step below the
## smallest value so floating-point rounding cannot drop it.
treepl_cv_settings <- function(smoothing) {
  s <- sort(unique(as.numeric(smoothing)), decreasing = TRUE)
  s <- s[is.finite(s) & s > 0]
  if (length(s) < 2L) {
    stop("--treepl-smoothing needs at least two positive values for treePL cross-validation.")
  }
  ratio <- s[-1L] / s[-length(s)]
  if (max(abs(log(ratio) - log(ratio[1L]))) > 1e-6) {
    stop(
      "--treepl-smoothing must be a geometric series (constant ratio between values), ",
      "because treePL cross-validation steps from cvstart to cvstop by cvmultstep. ",
      "Got: ", paste(s, collapse = ",")
    )
  }
  list(
    grid = s,
    cvstart = s[1L],
    cvstop = s[length(s)] * sqrt(ratio[1L]),
    cvmultstep = ratio[1L]
  )
}

treepl_cv_lines <- function(cv_settings, cv_type = c("loo", "random"), cvoutfile) {
  cv_type <- match.arg(cv_type)
  c(
    "cv",
    if (identical(cv_type, "random")) "randomcv" else NULL,
    paste0("cvstart = ", .treepl_fmt_num(cv_settings$cvstart)),
    paste0("cvstop = ", .treepl_fmt_num(cv_settings$cvstop)),
    paste0("cvmultstep = ", .treepl_fmt_num(cv_settings$cvmultstep)),
    paste0("cvoutfile = ", cvoutfile)
  )
}

## treePL writes one line per smoothing value: "chisq: (<smooth>) <score>".
## Printed smoothing values are rounded, so each is snapped to the nearest
## grid value on the log scale.
treepl_parse_cv_scores <- function(cvoutfile, grid) {
  if (!file.exists(cvoutfile)) return(data.frame(smooth = numeric(0), cv_chisq = numeric(0)))
  txt <- readLines(cvoutfile, warn = FALSE)
  hits <- regmatches(txt, regexec("^chisq:\\s*\\(([^)]+)\\)\\s*(\\S+)", txt))
  hits <- Filter(function(h) length(h) == 3L, hits)
  if (!length(hits)) return(data.frame(smooth = numeric(0), cv_chisq = numeric(0)))
  printed <- suppressWarnings(as.numeric(vapply(hits, `[`, "", 2L)))
  chisq <- suppressWarnings(as.numeric(vapply(hits, `[`, "", 3L)))
  snapped <- vapply(printed, function(x) grid[which.min(abs(log(grid) - log(x)))], numeric(1))
  data.frame(smooth = snapped, cv_chisq = chisq, stringsAsFactors = FALSE)
}

## Per-fold chi-square terms from a CV log. treePL prints a running total per
## fold: "(chisq): <total>" for random CV, "<tip> rate .. chisq <total>" for
## LOOCV. Returns one numeric vector of fold terms per smoothing value, in run
## order. Needs treePL's nthreads = 1 (the default), so the totals are serial.
treepl_parse_cv_folds <- function(log_path, cv_type = c("loo", "random")) {
  cv_type <- match.arg(cv_type)
  if (!file.exists(log_path)) return(list())
  txt <- readLines(log_path, warn = FALSE)
  sec <- grep("^smoothing:", txt)
  if (!length(sec)) return(list())
  pat <- if (identical(cv_type, "random")) "^\\(chisq\\):\\s*(\\S+)\\s*$" else
    "^\\S+ rate \\S+ expe \\S+ dur \\S+ len \\S+ sq \\S+ chisq (\\S+)\\s*$"
  bounds <- c(sec, length(txt) + 1L)
  out <- lapply(seq_along(sec), function(i) {
    z <- txt[seq.int(bounds[i], bounds[i + 1L] - 1L)]
    hit <- regmatches(z, regexec(pat, z))
    tot <- suppressWarnings(as.numeric(vapply(Filter(function(h) length(h) == 2L, hit), `[`, "", 2L)))
    diff(c(0, tot))
  })
  ## the run ends with a plain dating pass that has no folds
  Filter(length, out)
}

## A single fold can blow up when treePL drives a rate to ~0 and the expected
## length of a removed tip becomes tiny (chi-square terms of 1e+39 and more).
## treePL prints only a running total, so every fold after a blow-up is lost
## in rounding. The robust score therefore uses the folds before the first
## blow-up seen at any smoothing value of the run: the same folds for every
## value, so the scores stay comparable. A term blows up if it is not finite,
## is above 1e12, or is more than `fold_cap` times the median normal term.
treepl_robust_cv_scores <- function(fold_terms, fold_cap = 1e4) {
  n <- unique(lengths(fold_terms))
  if (length(n) != 1L || n < 1L) return(NULL)
  m <- do.call(rbind, fold_terms)
  normal <- m[is.finite(m) & m > 0 & m < 1e12]
  limit <- if (length(normal)) min(1e12, fold_cap * stats::median(normal)) else 1e12
  bad <- apply(m, 2L, function(col) any(!is.finite(col) | col > limit))
  n_used <- if (any(bad)) which(bad)[1L] - 1L else n
  list(scores = if (n_used < 1L) rep(NA_real_, nrow(m)) else rowSums(m[, seq_len(n_used), drop = FALSE]),
       n_folds = n, n_used = n_used)
}

## Step 2. Cross-validate the smoothing grid. Returns one row per grid value.
treepl_cross_validate <- function(treepl_bin, tree_path, node_bounds, cv_settings,
                                  cv_type, numsites, prime_lines,
                                  cfg_path, log_path, cvoutfile, thorough = TRUE,
                                  opt = NULL, plsimaniter = NULL, pliter = NULL,
                                  omp_threads = 1L, seed = NULL) {
  junk_out <- tempfile(pattern = "treepl_cv_", fileext = ".tre")
  on.exit(unlink(c(junk_out, paste0(junk_out, ".r8s")), force = TRUE), add = TRUE)
  if (file.exists(cvoutfile)) unlink(cvoutfile)
  treepl_write_cfg(
    cfg_path, tree_path, junk_out, cv_settings$cvstart, numsites, node_bounds,
    thorough = thorough, prime = FALSE,
    opt = opt, plsimaniter = plsimaniter, pliter = pliter,
    extra_lines = c(prime_lines, treepl_cv_lines(cv_settings, cv_type, cvoutfile)),
    seed = seed
  )
  treepl_run_one(treepl_bin, cfg_path, log_path, omp_threads = omp_threads)
  scores <- treepl_parse_cv_scores(cvoutfile, cv_settings$grid)
  out <- data.frame(smooth = cv_settings$grid, stringsAsFactors = FALSE)
  out$cv_chisq_treepl <- scores$cv_chisq[match(out$smooth, scores$smooth)]
  out$cv_chisq <- out$cv_chisq_treepl
  out$cv_folds <- NA_integer_
  out$cv_folds_used <- NA_integer_
  robust <- treepl_robust_cv_scores(treepl_parse_cv_folds(log_path, cv_type))
  if (!is.null(robust) && length(robust$scores) == nrow(out)) {
    out$cv_chisq <- robust$scores
    out$cv_folds <- robust$n_folds
    out$cv_folds_used <- robust$n_used
  }
  out
}

## Random CV with replicates (after A. Pyron's VGP treePL benchmark). One
## random-CV run draws 10 groups of 10% of the tips, so a single run is noisy
## and a single blown-up fold ruins it. The series is run `reps` times with
## seeds seed + 1..reps.
##
## Selection (cv_score): within each replicate every smoothing value uses the
## same tip groups, so the values are compared inside the replicate, on the
## folds before the first blow-up (treepl_robust_cv_scores), and ranked. The
## smoothing with the lowest mean rank across replicates wins. A replicate
## whose first fold already blows up carries no information and is skipped.
## On some trees most replicates contain a blow-up, so whole-replicate
## rejection would leave too few scores.
##
## Diagnostic (cv_chisq): Pyron's rule. A replicate's full score is rejected
## if it is not finite or above 1e12 (until three are accepted) or above
## `outlier_factor` x the median accepted score; the value's cv_chisq is the
## median of the accepted scores when at least half are accepted.
## cfg_path, log_path and cvoutfile are functions of the replicate index.
treepl_cv_replicated <- function(treepl_bin, tree_path, node_bounds, cv_settings, reps, seed,
                                 cfg_path, log_path, cvoutfile, ...,
                                 outlier_factor = 100, fallback_ceiling = 1e12) {
  long <- vector("list", reps)
  for (r in seq_len(reps)) {
    started <- Sys.time()
    sc <- treepl_cross_validate(treepl_bin, tree_path, node_bounds, cv_settings, "random",
                                cfg_path = cfg_path(r), log_path = log_path(r),
                                cvoutfile = cvoutfile(r), seed = seed + r, ...)
    long[[r]] <- data.frame(smooth = sc$smooth, rep = r, seed = seed + r,
                            score = sc$cv_chisq_treepl, score_folds_before_blowup = sc$cv_chisq,
                            folds_used = sc$cv_folds_used, stringsAsFactors = FALSE)
    ## Builds without the `seed` option seed random CV from the clock (in
    ## seconds), so keep replicate starts at least one second apart.
    if (!any(grepl("random number seed", readLines(log_path(r), warn = FALSE), fixed = TRUE))) {
      wait <- 1.1 - as.numeric(difftime(Sys.time(), started, units = "secs"))
      if (wait > 0) Sys.sleep(wait)
    }
  }
  long <- do.call(rbind, long)
  long$accepted <- FALSE
  for (sm in unique(long$smooth)) {
    i <- which(long$smooth == sm)
    ok <- is.finite(long$score[i]) & long$score[i] >= 0 & long$score[i] <= fallback_ceiling
    repeat {
      if (sum(ok) < 3L) break
      keep <- ok & long$score[i] <= outlier_factor * stats::median(long$score[i][ok])
      if (identical(keep, ok)) break
      ok <- keep
    }
    long$accepted[i] <- ok
  }
  ## paired ranks inside each informative replicate
  long$rank_in_rep <- NA_real_
  for (r in unique(long$rep)) {
    i <- which(long$rep == r)
    if (all(is.finite(long$score_folds_before_blowup[i]))) {
      long$rank_in_rep[i] <- rank(signif(long$score_folds_before_blowup[i], 6), ties.method = "average")
    }
  }
  need <- ceiling(reps / 2)
  out <- data.frame(smooth = cv_settings$grid, stringsAsFactors = FALSE)
  out$cv_reps <- reps
  out$cv_reps_ranked <- vapply(out$smooth, function(sm) sum(is.finite(long$rank_in_rep[long$smooth == sm])), integer(1))
  out$cv_score <- vapply(out$smooth, function(sm) {
    z <- long$rank_in_rep[long$smooth == sm]
    if (any(is.finite(z))) mean(z, na.rm = TRUE) else NA_real_
  }, numeric(1))
  out$cv_reps_accepted <- vapply(out$smooth, function(sm) sum(long$accepted & long$smooth == sm), integer(1))
  out$cv_chisq <- vapply(out$smooth, function(sm) {
    z <- long$score[long$accepted & long$smooth == sm]
    if (length(z) >= need) stats::median(z) else NA_real_
  }, numeric(1))
  out$cv_chisq_treepl <- NA_real_
  out$cv_folds <- NA_integer_
  out$cv_folds_used <- NA_integer_
  list(scores = out, replicates = long)
}

## Lowest score wins (cv_score for replicated random CV, else cv_chisq);
## exact ties favour the larger smoothing.
.treepl_best_smooth <- function(sc) {
  v <- if (!is.null(sc$cv_score)) sc$cv_score else sc$cv_chisq
  ok <- is.finite(v)
  if (!any(ok)) return(NA_real_)
  z <- sc[ok, , drop = FALSE]
  z <- z[order(v[ok], -z$smooth), , drop = FALSE]
  z$smooth[1L]
}

## Step 2 with edge extension. Run CV on `smoothing`. If the best value is the
## largest (or smallest) value, run CV again on that value plus two more values
## further out in the same direction (same ratio); repeat while the outermost
## value keeps winning, up to `max_extend` extra runs. The anchor value is
## rescored in each extension run, so every choice is made inside one CV run,
## and the search never turns back. With cv_type "random" and reps > 1 each
## run is replicated (treepl_cv_replicated). cfg_path, log_path and cvoutfile
## are functions of (run, rep); rep is NULL for an unreplicated run.
## Returns all score rows (with cv_run), the replicate rows, and the selection.
treepl_cv_search <- function(treepl_bin, tree_path, node_bounds, smoothing, cv_type,
                             max_extend, cfg_path, log_path, cvoutfile, ...,
                             reps = 1L, seed = 1L, limits = c(1e-8, 1e8)) {
  settings <- treepl_cv_settings(smoothing)
  step <- 1 / settings$cvmultstep
  out <- list()
  replicates <- list()
  run <- 1L
  direction <- 0
  at_edge <- FALSE
  replicated <- identical(cv_type, "random") && reps > 1L
  repeat {
    if (replicated) {
      res <- treepl_cv_replicated(treepl_bin, tree_path, node_bounds, settings, reps, seed,
                                  cfg_path = function(r) cfg_path(run, r),
                                  log_path = function(r) log_path(run, r),
                                  cvoutfile = function(r) cvoutfile(run, r), ...)
      sc <- res$scores
      res$replicates$cv_run <- run
      replicates[[run]] <- res$replicates
    } else {
      sc <- treepl_cross_validate(treepl_bin, tree_path, node_bounds, settings, cv_type,
                                  cfg_path = cfg_path(run, NULL), log_path = log_path(run, NULL),
                                  cvoutfile = cvoutfile(run, NULL), seed = seed, ...)
    }
    sc$cv_run <- run
    out[[run]] <- sc
    best <- .treepl_best_smooth(sc)
    if (!is.finite(best)) break
    hi <- isTRUE(all.equal(best, max(sc$smooth)))
    lo <- isTRUE(all.equal(best, min(sc$smooth)))
    if (direction == 0) direction <- if (hi) 1 else if (lo) -1 else 0
    at_edge <- (direction > 0 && hi) || (direction < 0 && lo)
    if (!at_edge || run > max_extend) break
    nxt <- if (direction > 0) best * step^(0:2) else best / step^(0:2)
    if (any(nxt < limits[1L] | nxt > limits[2L])) break
    settings <- treepl_cv_settings(nxt)
    run <- run + 1L
  }
  list(scores = do.call(rbind, lapply(out, function(x) x[, union(names(out[[1L]]), names(x)), drop = FALSE])),
       replicates = if (length(replicates)) do.call(rbind, replicates) else NULL,
       selected = .treepl_best_smooth(out[[length(out)]]), at_edge = isTRUE(at_edge))
}

## Final penalized-likelihood objective of a dating run (last
## "after opt calc:" line of the treePL log). Lower is better.
treepl_parse_final_objective <- function(log_path) {
  if (!file.exists(log_path)) return(NA_real_)
  txt <- readLines(log_path, warn = FALSE)
  hit <- grep("^after opt calc:", txt, value = TRUE)
  if (!length(hit)) return(NA_real_)
  suppressWarnings(as.numeric(sub("^after opt calc:\\s*", "", hit[length(hit)])))
}

## Step 3. Date the full tree once per optimizer setting at one smoothing.
## `variants` is a named list of extra config lines (optimizer lines); `seeds`
## is a vector of the same length. Every dated tree is checked against the
## calibration bounds (tolerance 0.001); a tree outside them is not usable.
## The usable run with the lowest objective is selected, and the ages of every
## usable run are compared with it node by node. `jobs` > 1 runs the fits in
## parallel (forked workers, Unix only).
treepl_date_variants <- function(treepl_bin, tree_path, phy, node_bounds, smooth, numsites,
                                 variants, seeds, out_tree, cfg_path, log_path,
                                 thorough = TRUE, opt = NULL, plsimaniter = NULL,
                                 pliter = NULL, omp_threads = 1L, jobs = 1L) {
  ref_sig <- .treepl_signature_index(phy)
  cal_idx <- match(as.character(node_bounds$node), names(ref_sig))
  fit_one <- function(v) {
    name <- names(variants)[v]
    unlink(c(out_tree(name), paste0(out_tree(name), ".r8s")), force = TRUE)
    treepl_write_cfg(cfg_path(name), tree_path, out_tree(name), smooth, numsites, node_bounds,
                     thorough = thorough, prime = FALSE, opt = opt,
                     plsimaniter = plsimaniter, pliter = pliter,
                     extra_lines = variants[[v]], seed = seeds[v])
    started <- Sys.time()
    exit <- treepl_run_one(treepl_bin, cfg_path(name), log_path(name), omp_threads = omp_threads)
    seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    tr <- .treepl_read_valid_tree(out_tree(name))
    ok <- !is.null(tr) && ape::Ntip(tr) == ape::Ntip(phy)
    a <- if (ok) .treepl_ages_by_signature(tr, ref_sig) else NULL
    ok <- ok && all(is.finite(a))
    in_bounds <- ok && all(a[cal_idx] >= node_bounds$age_min - 1e-3 &
                             a[cal_idx] <= node_bounds$age_max + 1e-3)
    list(ages = if (in_bounds) a else NULL, row = data.frame(
      variant = name, seed = seeds[v],
      opt_lines = paste(variants[[v]], collapse = "; "),
      exit_code = as.integer(exit),
      within_bounds = in_bounds,
      objective = if (in_bounds) treepl_parse_final_objective(log_path(name)) else NA_real_,
      root_age = if (ok) max(a) else NA_real_,
      seconds = round(seconds, 2),
      tree_file = if (in_bounds) out_tree(name) else "",
      stringsAsFactors = FALSE
    ))
  }
  idx <- seq_along(variants)
  res <- if (jobs > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(idx, fit_one, mc.cores = jobs, mc.preschedule = FALSE)
  } else {
    lapply(idx, fit_one)
  }
  failed <- vapply(res, function(r) !is.list(r) || is.null(r$row), logical(1))
  if (any(failed)) stop("treePL optimizer fits crashed for ", sum(failed), " setting(s).")
  tab <- do.call(rbind, lapply(res, `[[`, "row"))
  ages <- lapply(res, `[[`, "ages")
  usable <- which(is.finite(tab$objective))
  if (!length(usable)) stop("treePL wrote no valid dated tree within the bounds under any optimizer setting; see ",
                            log_path(names(variants)[1L]))
  best <- usable[which.min(tab$objective[usable])]
  tab$selected <- seq_len(nrow(tab)) == best
  tab$max_abs_age_diff <- NA_real_
  tab$median_rel_age_diff <- NA_real_
  for (v in usable) {
    d <- abs(ages[[v]] - ages[[best]])
    tab$max_abs_age_diff[v] <- max(d)
    tab$median_rel_age_diff[v] <- stats::median(d / pmax(ages[[best]], 1e-12))
  }
  names(ages) <- tab$variant
  list(table = tab, selected = tab$variant[best], lines = variants[[best]],
       seed = seeds[best], tree_file = tab$tree_file[best], ages = ages)
}

## The full optimizer grid of A. Pyron's VGP benchmark: opt 0-5 x optad 0-5 x
## moredetail x moredetailad (144 settings), as config lines named o-a-m-d.
treepl_optimizer_grid <- function() {
  g <- expand.grid(opt = 0:5, optad = 0:5, md = 0:1, mdad = 0:1)
  out <- lapply(seq_len(nrow(g)), function(i) c(
    paste0("opt = ", g$opt[i]), paste0("optad = ", g$optad[i]),
    if (g$md[i] == 1L) "moredetail" else NULL,
    if (g$mdad[i] == 1L) "moredetailad" else NULL
  ))
  names(out) <- sprintf("o%d-a%d-m%d-d%d", g$opt, g$optad, g$md, g$mdad)
  out
}

## Rank optimizer settings from a screen table (one row per setting x seed;
## `setting` names the setting). A setting is usable only if every seed gave
## a tree within the bounds. Usable settings are ranked by median, then worst,
## then best objective, then total run time.
treepl_rank_settings <- function(tab, n_seeds) {
  sp <- split(tab, tab$setting)
  rk <- do.call(rbind, lapply(names(sp), function(k) {
    z <- sp[[k]]
    v <- z$objective[is.finite(z$objective)]
    data.frame(setting = k, n_fits = nrow(z), n_valid = length(v),
               usable = length(v) == n_seeds,
               median = if (length(v)) stats::median(v) else NA_real_,
               worst = if (length(v)) max(v) else NA_real_,
               best = if (length(v)) min(v) else NA_real_,
               seconds = sum(z$seconds), stringsAsFactors = FALSE)
  }))
  rk <- rk[order(!rk$usable, rk$median, rk$worst, rk$best, rk$seconds), , drop = FALSE]
  rk$rank <- ifelse(rk$usable, cumsum(rk$usable), NA_integer_)
  rk
}

.treepl_read_valid_tree <- function(path) {
  if (!file.exists(path)) return(NULL)
  tr <- try(ape::read.tree(path), silent = TRUE)
  if (inherits(tr, "try-error") || !inherits(tr, "phylo")) return(NULL)
  if (is.null(tr$edge.length) || any(!is.finite(tr$edge.length)) || any(tr$edge.length < 0)) {
    return(NULL)
  }
  tr
}

.treepl_poisson_perturb_phy <- function(phy, n_sites, min_edge = 1e-12) {
  if (!is.finite(n_sites) || n_sites <= 0) {
    stop("n_sites must be a positive finite number for treePL bootstrap.")
  }
  out <- phy
  lambda <- pmax(as.numeric(phy$edge.length) * n_sites, 0)
  out$edge.length <- stats::rpois(length(lambda), lambda) / n_sites
  out$edge.length <- pmax(out$edge.length, min_edge)
  out
}

## Internal-node signatures (sorted descendant tip labels), named by node number.
.treepl_signature_index <- function(phy) {
  parts <- ape::prop.part(phy)
  labs <- attr(parts, "labels")
  sig <- vapply(parts, function(p) paste(sort(labs[p]), collapse = "\r"), character(1))
  names(sig) <- as.character(seq.int(ape::Ntip(phy) + 1L, length.out = length(parts)))
  sig
}

.treepl_node_age_vector <- function(dated_tree) {
  depths <- ape::node.depth.edgelength(dated_tree)
  root_to_tip <- max(depths[seq_len(ape::Ntip(dated_tree))], na.rm = TRUE)
  ages <- root_to_tip - depths
  storage.mode(ages) <- "double"
  ages
}

.treepl_ages_by_signature <- function(dated_tree, ref_sig) {
  tree_sig <- .treepl_signature_index(dated_tree)
  node_age <- .treepl_node_age_vector(dated_tree)
  age_by_sig <- stats::setNames(as.numeric(node_age[as.integer(names(tree_sig))]), tree_sig)
  out <- unname(age_by_sig[ref_sig])
  storage.mode(out) <- "double"
  out
}

.treepl_bootstrap_one <- function(i, phy, ref_sig, treepl_bin, smooth, numsites,
                                  node_bounds, prime_lines, thorough, opt,
                                  plsimaniter, pliter, omp_threads, min_edge,
                                  n_sites, workdir, prefix,
                                  return_tree = FALSE, seed = NULL) {
  if (is.finite(seed)) set.seed(as.integer(seed))
  phy_boot <- try(.treepl_poisson_perturb_phy(phy, n_sites, min_edge = min_edge), silent = TRUE)
  if (inherits(phy_boot, "try-error")) {
    return(list(ok = FALSE, ages = NULL, tree = NULL))
  }
  in_tree <- file.path(workdir, sprintf("%s_boot_%03d_in.tre", prefix, i))
  cfg_path <- file.path(workdir, sprintf("%s_boot_%03d.cfg", prefix, i))
  out_tree <- file.path(workdir, sprintf("%s_boot_%03d.tre", prefix, i))
  log_path <- file.path(workdir, sprintf("%s_boot_%03d.log", prefix, i))
  ape::write.tree(phy_boot, in_tree)
  treepl_write_cfg(
    cfg_path, in_tree, out_tree, smooth, numsites, node_bounds,
    thorough = thorough, prime = FALSE, opt = opt,
    plsimaniter = plsimaniter, pliter = pliter, extra_lines = prime_lines,
    seed = if (is.finite(seed)) seed %% 2147483647L else NULL
  )
  treepl_run_one(treepl_bin, cfg_path, log_path, omp_threads = omp_threads)
  tr_i <- .treepl_read_valid_tree(out_tree)
  if (is.null(tr_i)) {
    return(list(ok = FALSE, ages = NULL, tree = NULL))
  }
  ages_i <- .treepl_ages_by_signature(tr_i, ref_sig)
  if (!all(is.finite(ages_i))) {
    return(list(ok = FALSE, ages = NULL, tree = NULL))
  }
  list(ok = TRUE, ages = ages_i, tree = if (isTRUE(return_tree)) tr_i else NULL)
}

## node_bounds: the shared merged bounds (node, taxonA, taxonB, age_min,
## age_max) mapped on `phy`. smooth and prime_lines come from the main run,
## so every replicate uses the same settings as the point estimate.
treepl_bootstrap_ci <- function(treepl_bin, phy, node_bounds, smooth, dated_tree,
                                prime_lines = character(0), B = 100L,
                                n_sites = 1000L, numsites = as.integer(n_sites),
                                thorough = TRUE,
                                opt = NULL, plsimaniter = NULL, pliter = NULL,
                                quiet = FALSE, omp_threads = 1L,
                                jobs = 1L,
                                min_edge = 1e-12, trees = FALSE,
                                workdir = tempdir(), prefix = "treepl_bootstrap") {
  B <- as.integer(B)
  if (!is.finite(B) || B < 1L) stop("B must be a positive integer.")
  if (!is.finite(n_sites) || n_sites <= 0) {
    stop("n_sites must be a positive finite number for treePL bootstrap.")
  }
  if (!is.finite(numsites) || numsites < 1L) {
    stop("numsites must be a positive integer for treePL.")
  }
  jobs <- as.integer(jobs)
  if (!is.finite(jobs) || jobs < 1L) stop("jobs must be a positive integer.")
  if (!nzchar(treepl_bin) || !file.exists(treepl_bin)) {
    stop("A valid treePL binary is required for treePL bootstrap.")
  }
  if (!inherits(dated_tree, "phylo")) stop("treePL bootstrap needs the dated point-estimate tree.")

  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  ref_sig <- .treepl_signature_index(phy)
  int_nodes <- as.integer(names(ref_sig))
  base_age <- .treepl_ages_by_signature(dated_tree, ref_sig)
  if (!all(is.finite(base_age))) {
    stop("treePL bootstrap could not map all base node ages onto the reference topology.")
  }

  boot_age <- matrix(NA_real_, nrow = length(int_nodes), ncol = B)
  boot_trees <- if (isTRUE(trees)) vector("list", B) else NULL
  ## Some empirical datasets have a low bootstrap convergence rate in treePL
  ## even when the selected full-tree solution is stable. Allow a much larger
  ## retry budget so CI generation can complete without changing the target B.
  max_attempts <- max(as.integer(B * 20L), as.integer(B + 100L))
  rep_seeds <- sample.int(.Machine$integer.max, max_attempts)
  n_ok <- 0L
  attempt_i <- 0L

  record_bootstrap_success <- function(res_i) {
    if (!is.list(res_i) || !isTRUE(res_i$ok) || n_ok >= B) return(FALSE)
    n_ok <<- n_ok + 1L
    boot_age[, n_ok] <<- res_i$ages
    if (isTRUE(trees)) boot_trees[[n_ok]] <<- res_i$tree
    TRUE
  }

  run_one <- function(i) {
    .treepl_bootstrap_one(
      i = i, phy = phy, ref_sig = ref_sig, treepl_bin = treepl_bin,
      smooth = smooth, numsites = numsites, node_bounds = node_bounds,
      prime_lines = prime_lines, thorough = thorough, opt = opt,
      plsimaniter = plsimaniter, pliter = pliter, omp_threads = omp_threads,
      min_edge = min_edge, n_sites = n_sites, workdir = workdir, prefix = prefix,
      return_tree = trees, seed = rep_seeds[i]
    )
  }

  if (jobs > 1L && .Platform$OS.type == "unix" && !quiet) {
    cat("Running treePL bootstrap in parallel with ", jobs, " workers\n", sep = "")
  }
  while (n_ok < B && attempt_i < max_attempts) {
    if (jobs <= 1L || .Platform$OS.type != "unix") {
      attempt_i <- attempt_i + 1L
      res_list <- list(run_one(attempt_i))
    } else {
      batch_ids <- seq.int(attempt_i + 1L, length.out = min(jobs, max_attempts - attempt_i))
      attempt_i <- max(batch_ids)
      res_list <- parallel::mclapply(batch_ids, run_one, mc.cores = jobs,
                                     mc.preschedule = FALSE, mc.set.seed = FALSE)
    }
    for (res_i in res_list) {
      record_bootstrap_success(res_i)
      if (n_ok >= B) break
    }
    if (!quiet) {
      cat("\rRunning treePL bootstrap: ", min(n_ok + 1L, B), "/", B,
          " (attempt ", attempt_i, ")", sep = "")
    }
  }
  if (!quiet) cat("\n")
  if (n_ok < B) {
    stop(
      "treePL bootstrap produced only ", n_ok,
      " successful replicates out of ", attempt_i,
      " attempts."
    )
  }

  qfun <- function(probs) {
    apply(
      boot_age,
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
    age = base_age,
    ci_lo = qfun(0.025),
    ci_hi = qfun(0.975),
    q50_lo = qfun(0.25),
    q50_hi = qfun(0.75),
    se = NA_real_,
    stringsAsFactors = FALSE
  )

  out <- list(
    tree = dated_tree,
    ci = ci_df,
    bounds = node_bounds,
    n_bootstrap_ok = n_ok,
    n_bootstrap_attempts = attempt_i
  )
  if (isTRUE(trees)) out$trees <- boot_trees[seq_len(n_ok)]
  out
}

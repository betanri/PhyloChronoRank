##############################################################
## treepl_helpers.R
##
## Repo-local treePL support (Smith & O'Meara 2012). The dating
## engine runs treePL in the standard order:
##
##   1. prime  — treePL tests optimizers and prints the best
##               opt/optad/optcvad/moredetail lines, then exits
##   2. cv     — cross-validation over a geometric smoothing series
##               (cvstart -> cvstop by cvmultstep); the smoothing
##               with the lowest CV chi-square is selected
##   3. final  — one thorough run with the selected smoothing and
##               the primed optimizer lines
##
## treePL's penalized-likelihood objective is not comparable across
## smoothing values (the penalty grows with smoothing), so it is
## never used to pick smoothing.
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
                             extra_lines = NULL) {
  lines <- c(
    paste0("treefile = ", tree_path),
    paste0("outfile = ", out_tree_path),
    paste0("numsites = ", as.integer(numsites)),
    paste0("smooth = ", .treepl_fmt_num(smooth))
  )
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
                         omp_threads = 1L) {
  junk_out <- tempfile(pattern = "treepl_prime_", fileext = ".tre")
  on.exit(unlink(c(junk_out, paste0(junk_out, ".r8s")), force = TRUE), add = TRUE)
  treepl_write_cfg(
    cfg_path, tree_path, junk_out, smooth, numsites, node_bounds,
    thorough = thorough, prime = TRUE,
    opt = opt, plsimaniter = plsimaniter, pliter = pliter
  )
  treepl_run_one(treepl_bin, cfg_path, log_path, omp_threads = omp_threads)
  lines <- treepl_parse_prime_lines(log_path)
  if (!length(lines)) {
    stop("treePL prime step did not report optimizer settings; see ", log_path)
  }
  lines
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

## Step 2. Cross-validate the smoothing grid. Returns one row per grid value.
treepl_cross_validate <- function(treepl_bin, tree_path, node_bounds, cv_settings,
                                  cv_type, numsites, prime_lines,
                                  cfg_path, log_path, cvoutfile, thorough = TRUE,
                                  opt = NULL, plsimaniter = NULL, pliter = NULL,
                                  omp_threads = 1L) {
  junk_out <- tempfile(pattern = "treepl_cv_", fileext = ".tre")
  on.exit(unlink(c(junk_out, paste0(junk_out, ".r8s")), force = TRUE), add = TRUE)
  if (file.exists(cvoutfile)) unlink(cvoutfile)
  treepl_write_cfg(
    cfg_path, tree_path, junk_out, cv_settings$cvstart, numsites, node_bounds,
    thorough = thorough, prime = FALSE,
    opt = opt, plsimaniter = plsimaniter, pliter = pliter,
    extra_lines = c(prime_lines, treepl_cv_lines(cv_settings, cv_type, cvoutfile))
  )
  treepl_run_one(treepl_bin, cfg_path, log_path, omp_threads = omp_threads)
  scores <- treepl_parse_cv_scores(cvoutfile, cv_settings$grid)
  out <- data.frame(smooth = cv_settings$grid, stringsAsFactors = FALSE)
  out$cv_chisq <- scores$cv_chisq[match(out$smooth, scores$smooth)]
  out
}

treepl_select_smoothing <- function(cv_scores) {
  ok <- is.finite(cv_scores$cv_chisq)
  if (!any(ok)) return(NA_real_)
  z <- cv_scores[ok, , drop = FALSE]
  z$smooth[which.min(z$cv_chisq)]
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
    plsimaniter = plsimaniter, pliter = pliter, extra_lines = prime_lines
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

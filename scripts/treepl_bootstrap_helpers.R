##############################################################
## treepl_bootstrap_helpers.R
##
## Repo-local treePL bootstrap CI support.
## Branch lengths are perturbed with the same Poisson resampling
## scheme used in the shared PCR uncertainty layer, treePL is rerun
## on each perturbed phylogram, and node ages are summarized with
## empirical quantiles after matching internal nodes by descendant
## tip-set signatures.
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

.treepl_children_lookup <- function(phy) {
  total <- ape::Ntip(phy) + phy$Nnode
  ch <- vector("list", total)
  for (k in seq_len(nrow(phy$edge))) {
    ch[[phy$edge[k, 1L]]] <- c(ch[[phy$edge[k, 1L]]], phy$edge[k, 2L])
  }
  ch
}

.treepl_descendant_tip_labels <- function(phy, node, children = NULL) {
  if (is.null(children)) children <- .treepl_children_lookup(phy)
  n_tip <- ape::Ntip(phy)
  stack <- node
  tips <- integer(0)
  while (length(stack)) {
    cur <- stack[[1L]]
    stack <- stack[-1L]
    if (cur <= n_tip) {
      tips <- c(tips, cur)
    } else {
      stack <- c(children[[cur]], stack)
    }
  }
  sort(phy$tip.label[unique(tips)])
}

.treepl_select_node_pair <- function(phy, node, children = NULL) {
  if (is.null(children)) children <- .treepl_children_lookup(phy)
  kids <- children[[node]]
  if (length(kids) < 2L) {
    stop("Could not identify two child clades for treePL node ", node, ".")
  }
  child_min_tips <- vapply(
    kids,
    function(nd) .treepl_descendant_tip_labels(phy, nd, children)[1L],
    FUN.VALUE = character(1)
  )
  kids <- kids[order(child_min_tips)]
  c(
    .treepl_descendant_tip_labels(phy, kids[1L], children)[1L],
    .treepl_descendant_tip_labels(phy, kids[2L], children)[1L]
  )
}

.treepl_build_node_bounds <- function(phy, calibration_df) {
  if (!exists("reltime_merge_calibration_bounds", mode = "function")) {
    stop("treePL bootstrap helper requires reltime_merge_calibration_bounds() from reltime_helpers.R.")
  }
  bounds <- reltime_merge_calibration_bounds(phy, calibration_df)
  children <- .treepl_children_lookup(phy)
  pairs <- t(vapply(
    bounds$node,
    function(nd) .treepl_select_node_pair(phy, as.integer(nd), children),
    FUN.VALUE = c("", "")
  ))
  data.frame(
    node = bounds$node,
    taxonA = pairs[, 1L],
    taxonB = pairs[, 2L],
    age_min = bounds$age_min,
    age_max = bounds$age_max,
    stringsAsFactors = FALSE
  )
}

treepl_write_cfg <- function(cfg_path, tree_path, out_tree_path, smooth, numsites,
                             node_bounds, thorough = TRUE, prime = TRUE,
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

.treepl_signature_index <- function(phy) {
  n_tip <- ape::Ntip(phy)
  int_nodes <- seq.int(n_tip + 1L, n_tip + phy$Nnode)
  children <- .treepl_children_lookup(phy)
  sig <- vapply(
    int_nodes,
    function(nd) paste(.treepl_descendant_tip_labels(phy, nd, children), collapse = "\r"),
    FUN.VALUE = character(1)
  )
  names(sig) <- as.character(int_nodes)
  sig
}

.treepl_node_age_vector <- function(dated_tree) {
  depths <- ape::node.depth.edgelength(dated_tree)
  if (length(depths) != ape::Ntip(dated_tree) + dated_tree$Nnode) {
    stop("treePL bootstrap could not compute node depths for the dated tree.")
  }
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
                                  n_sites, workdir, prefix, quiet,
                                  return_tree = FALSE, seed = NULL) {
  if (is.finite(seed)) set.seed(as.integer(seed))
  if (!quiet) cat("\rRunning treePL bootstrap:", i)
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

treepl_bootstrap_ci <- function(treepl_bin, phy, calibration_df, smooth,
                                dated_tree = NULL, B = 100L,
                                n_sites = 1000L, numsites = as.integer(n_sites),
                                thorough = TRUE, prime = TRUE,
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

  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  node_bounds <- .treepl_build_node_bounds(phy, calibration_df)
  ref_sig <- .treepl_signature_index(phy)

  base_tree_path <- file.path(workdir, paste0(prefix, "_base_input.tre"))
  ape::write.tree(phy, base_tree_path)

  prime_lines <- character(0)
  if (isTRUE(prime)) {
    prime_cfg <- file.path(workdir, paste0(prefix, "_prime.cfg"))
    prime_out <- file.path(workdir, paste0(prefix, "_prime_out.tre"))
    prime_log <- file.path(workdir, paste0(prefix, "_prime.log"))
    treepl_write_cfg(
      prime_cfg, base_tree_path, prime_out, smooth, numsites, node_bounds,
      thorough = thorough, prime = TRUE, opt = opt,
      plsimaniter = plsimaniter, pliter = pliter
    )
    treepl_run_one(treepl_bin, prime_cfg, prime_log, omp_threads = omp_threads)
    prime_lines <- treepl_parse_prime_lines(prime_log)
    if (file.exists(prime_out)) unlink(prime_out)
    if (file.exists(paste0(prime_out, ".r8s"))) unlink(paste0(prime_out, ".r8s"))
  }

  if (is.null(dated_tree)) {
    base_cfg <- file.path(workdir, paste0(prefix, "_base.cfg"))
    base_out <- file.path(workdir, paste0(prefix, "_base.tre"))
    base_log <- file.path(workdir, paste0(prefix, "_base.log"))
    treepl_write_cfg(
      base_cfg, base_tree_path, base_out, smooth, numsites, node_bounds,
      thorough = thorough, prime = FALSE, opt = opt,
      plsimaniter = plsimaniter, pliter = pliter, extra_lines = prime_lines
    )
    treepl_run_one(treepl_bin, base_cfg, base_log, omp_threads = omp_threads)
    dated_tree <- .treepl_read_valid_tree(base_out)
    if (is.null(dated_tree)) stop("treePL bootstrap could not generate a valid base dated tree.")
  } else if (is.character(dated_tree)) {
    dated_tree <- .treepl_read_valid_tree(dated_tree)
    if (is.null(dated_tree)) stop("treePL bootstrap could not read the supplied dated tree.")
  }

  int_nodes <- as.integer(names(ref_sig))
  base_age <- .treepl_ages_by_signature(dated_tree, ref_sig)
  if (!all(is.finite(base_age))) {
    stop("treePL bootstrap could not map all base node ages onto the reference topology.")
  }

  boot_age <- matrix(NA_real_, nrow = length(int_nodes), ncol = B)
  boot_trees <- if (isTRUE(trees)) vector("list", B) else NULL
  ok <- logical(B)
  rep_seeds <- sample.int(.Machine$integer.max, B)

  if (jobs <= 1L || .Platform$OS.type != "unix") {
    for (i in seq_len(B)) {
      if (!quiet) cat("\rRunning treePL bootstrap:", i, "/", B)
      res_i <- .treepl_bootstrap_one(
        i = i,
        phy = phy,
        ref_sig = ref_sig,
        treepl_bin = treepl_bin,
        smooth = smooth,
        numsites = numsites,
        node_bounds = node_bounds,
        prime_lines = prime_lines,
        thorough = thorough,
        opt = opt,
        plsimaniter = plsimaniter,
        pliter = pliter,
        omp_threads = omp_threads,
        min_edge = min_edge,
        n_sites = n_sites,
        workdir = workdir,
        prefix = prefix,
        quiet = TRUE,
        return_tree = trees,
        seed = rep_seeds[i]
      )
      if (!isTRUE(res_i$ok)) next
      boot_age[, i] <- res_i$ages
      ok[i] <- TRUE
      if (isTRUE(trees)) boot_trees[[i]] <- res_i$tree
    }
  } else {
    if (!quiet) {
      cat("Running treePL bootstrap in parallel with ", jobs, " workers\n", sep = "")
    }
    res_list <- parallel::mclapply(
      seq_len(B),
      function(i) {
        .treepl_bootstrap_one(
          i = i,
          phy = phy,
          ref_sig = ref_sig,
          treepl_bin = treepl_bin,
          smooth = smooth,
          numsites = numsites,
          node_bounds = node_bounds,
          prime_lines = prime_lines,
          thorough = thorough,
          opt = opt,
          plsimaniter = plsimaniter,
          pliter = pliter,
          omp_threads = omp_threads,
          min_edge = min_edge,
          n_sites = n_sites,
          workdir = workdir,
          prefix = prefix,
          quiet = TRUE,
          return_tree = trees,
          seed = rep_seeds[i]
        )
      },
      mc.cores = jobs,
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
    for (i in seq_along(res_list)) {
      res_i <- res_list[[i]]
      if (!is.list(res_i) || !isTRUE(res_i$ok)) next
      boot_age[, i] <- res_i$ages
      ok[i] <- TRUE
      if (isTRUE(trees)) boot_trees[[i]] <- res_i$tree
    }
  }
  if (!quiet) cat("\n")
  if (!any(ok)) stop("treePL bootstrap produced no successful replicates")

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
    bootstrap_ok = ok,
    n_bootstrap_ok = sum(ok)
  )
  if (isTRUE(trees)) out$trees <- boot_trees[ok]
  out
}

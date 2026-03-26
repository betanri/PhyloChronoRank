#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
kv <- list()
for (a in args) {
  if (!grepl("^--", a)) next
  parts <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  key <- parts[1]
  val <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "TRUE"
  kv[[key]] <- val
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/run_dating_grid.R",
      "--phylogram=TARGET.tre",
      "--outdir=OUT",
      "[--calibrations-csv=CALS.csv | --reference-time-tree=REF_TIME.tre]",
      "[--methods=chronos,treepl,reltime]",
      "[--chronos-lambdas=0.01,0.1,1,10,100]",
      "[--treepl-smoothing=0.01,0.1,1,10,100]",
      "[--chronos-discrete-k=2,3,5]",
      "[--chronos-attempt-timeout=90]",
      "[--treepl-bin=/path/to/treePL]",
      "[--megacc-bin=/path/to/megacc]",
      "[--root-age=123.4]",
      "[--ci-sites=1000]",
      "[--chronos-ci-type=parametric]",
      "[--chronos-ci-reps=100]",
      "[--treepl-bootstrap-reps=0]",
      "[--treepl-bootstrap-jobs=1]",
      "[--reltime-bootstrap-reps=100]",
      "[--chronos-pml-rds=FIT.rds]",
      "[--calibration-tag=NAME]",
      sep = " "
    ),
    "\n\n",
    "Required arguments:\n",
    "  --phylogram             Rooted target phylogram.\n",
    "  --outdir                Output directory.\n",
    "  Exactly one of:\n",
    "    --calibrations-csv    Pairwise calibration CSV with taxonA,taxonB,age_min[,age_max].\n",
    "    --reference-time-tree Reference ultrametric tree for congruification.\n\n",
    "Optional arguments:\n",
    "  --methods               Comma list from chronos,treepl,reltime,reltime_mega. Default: chronos,treepl,reltime.\n",
    "  --megacc-bin            megacc executable for MEGA RelTime. Default: MEGACC_BIN env, PATH megacc.\n",
    "  --phylogram-tree-name   Select named tree from a named-Newick or multi-tree file.\n",
    "  --phylogram-tree-index  1-based fallback tree index. Default: 1.\n",
    "  --reference-tree-name   Optional named tree selector for the reference time tree.\n",
    "  --reference-tree-index  1-based fallback index for the reference tree. Default: 1.\n",
    "  --calibration-tag       If the calibration CSV has a candidate column, keep rows tagged with this value plus blank/all rows.\n",
    "  --root-age              Optional exact root age appended to the shared calibration set for all methods.\n",
    "  --out-prefix            Prefix for output files. Default: target tree id.\n",
    "  --chronos-lambdas       Lambda grid. Default: 0.01,0.1,1,10,100.\n",
    "  --chronos-models        Clock models. Default: clock,correlated,relaxed,discrete.\n",
    "  --chronos-discrete-k    nb.rate.cat grid for the discrete model. Default: 2,3,5.\n",
    "  --chronos-retries       Retries per chronos setting. Default: 2.\n",
    "  --chronos-attempt-timeout Maximum elapsed seconds allowed per chronos attempt. Default: 0 (no timeout).\n",
    "  --treepl-smoothing      Smoothing grid. Default: 0.01,0.1,1,10,100.\n",
    "  --treepl-bin            treePL executable. Default: TREEPL_BIN env, PATH treePL, then ../treePL.\n",
    "  --treepl-numsites       numsites value written to treePL configs. Default: 1000.\n",
    "  --treepl-threads        OMP_NUM_THREADS for treePL. Default: 1.\n",
    "  --treepl-thorough       Use thorough optimization in treePL (strongly recommended). Default: TRUE.\n",
    "  --treepl-prime          Prime treePL optimizer before full run. Default: TRUE.\n",
    "  --treepl-opt            treePL optimization algorithm (1=LM, 2=stochastic, 4=auto). Default: not set.\n",
    "  --treepl-plsimaniter    Simulated-annealing iterations for treePL. Default: not set (treePL default).\n",
    "  --treepl-pliter         Penalty-likelihood iterations for treePL. Default: not set (treePL default).\n",
    "  --chronos-iter-max      Maximum iterations for chronos optimizer. Default: 10000 (ape default).\n",
    "  --chronos-tol           Convergence tolerance for chronos. Default: 1e-8 (ape default).\n",
    "  --ci-sites              Alignment length used for optional uncertainty summaries. Default: 1000.\n",
    "  --chronos-ci-type       chronos bootstrap type: nonparametric, semiparametric, or parametric. Default: parametric.\n",
    "  --chronos-ci-reps       chronos bootstrap replicates. Default: 100.\n",
    "  --treepl-bootstrap-reps treePL bootstrap replicates for the shared uncertainty layer. Default: 0 (skip treePL CI).\n",
    "  --treepl-bootstrap-jobs Parallel workers for treePL bootstrap replicates. Default: 1.\n",
    "  --reltime-bootstrap-reps RelTime bootstrap replicates for the shared uncertainty layer. Default: 100.\n",
    "  --chronos-pml-rds       Optional RDS containing a phangorn pml.output for nonparametric/semiparametric chronos bootstrap.\n",
    "  --reltime-sites         Backward-compatible alias for --ci-sites.\n",
    "  --help                  Print this message.\n",
    sep = ""
  )
}

if ("help" %in% names(kv) || !length(args)) {
  usage()
  quit(save = "no", status = 0)
}

need <- c("phylogram", "outdir")
miss <- need[!need %in% names(kv)]
if (length(miss)) {
  usage()
  stop("Missing required args: ", paste(miss, collapse = ", "))
}

has_csv <- "calibrations-csv" %in% names(kv)
has_ref <- "reference-time-tree" %in% names(kv)
if (identical(has_csv, has_ref)) {
  usage()
  stop("Provide exactly one of --calibrations-csv or --reference-time-tree.")
}

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package ape is required.")
}
suppressPackageStartupMessages(library(ape))

## quadprog is no longer required -- calibration bounds are now enforced via
## MEGA-style local rescaling (see .reltime_project_node_ages_local).

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file_arg)) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
repo_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "reltime_helpers.R"))
source(file.path(script_dir, "chronos_ci_helpers.R"))
source(file.path(script_dir, "treepl_bootstrap_helpers.R"))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

parse_num_grid <- function(x, arg_name) {
  vals <- as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1L]]))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) stop(arg_name, " must contain at least one numeric value.")
  unique(vals)
}

parse_int_grid <- function(x, arg_name) {
  vals <- suppressWarnings(as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1L]])))
  vals <- vals[is.finite(vals) & vals >= 1L]
  if (!length(vals)) stop(arg_name, " must contain at least one positive integer.")
  unique(vals)
}

parse_chr_grid <- function(x, arg_name) {
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  vals <- vals[nzchar(vals)]
  if (!length(vals)) stop(arg_name, " must contain at least one value.")
  unique(vals)
}

sanitize_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) x <- "dating_run"
  x
}

fmt_num <- function(x) {
  format(as.numeric(x), scientific = FALSE, trim = TRUE, digits = 15L)
}

fmt_token <- function(x) {
  out <- fmt_num(x)
  out <- gsub("-", "m", out, fixed = TRUE)
  gsub("\\.", "p", out)
}

bind_rows_fill <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (!length(lst)) return(data.frame())
  cols <- unique(unlist(lapply(lst, names), use.names = FALSE))
  padded <- lapply(lst, function(df) {
    miss <- setdiff(cols, names(df))
    for (nm in miss) df[[nm]] <- NA
    df[, cols, drop = FALSE]
  })
  do.call(rbind, padded)
}

build_ci_inputs <- function(phy, node_bounds) {
  root_node <- Ntip(phy) + 1L
  calib_df <- data.frame(
    node = node_bounds$node,
    age.min = node_bounds$age_min,
    age.max = node_bounds$age_max,
    stringsAsFactors = FALSE
  )
  cal_min <- as.list(stats::setNames(calib_df$age.min, calib_df$node))
  cal_max <- as.list(stats::setNames(calib_df$age.max, calib_df$node))

  root_var <- 0
  root_row <- calib_df[calib_df$node == root_node, , drop = FALSE]
  if (nrow(root_row) &&
      is.finite(root_row$age.min[1L]) &&
      is.finite(root_row$age.max[1L]) &&
      abs(root_row$age.max[1L] - root_row$age.min[1L]) > 1e-10) {
    root_var <- ((root_row$age.max[1L] - root_row$age.min[1L]) / (2 * 1.96))^2
  }

  list(
    root_var = root_var,
    cal_min = cal_min,
    cal_max = cal_max,
    calib_bundle = list(chronos_calib = calib_df)
  )
}

make_unscored_uncertainty_row <- function(candidate, source = "not_scored") {
  ci_width_summary(
    candidate,
    data.frame(ci_lo = numeric(), ci_hi = numeric(), stringsAsFactors = FALSE),
    source = source
  )
}

load_chronos_pml_output <- function(path) {
  if (!nzchar(path)) return(NULL)
  obj <- readRDS(normalizePath(path, winslash = "/", mustWork = TRUE))
  if (is.list(obj) && !is.null(obj$pml_output)) obj <- obj$pml_output
  if (!is.list(obj) || is.null(obj$tree)) {
    stop("--chronos-pml-rds must contain a phangorn pml.output object or a list with $pml_output.")
  }
  obj
}

is_nexus_file <- function(path, max_lines = 50L) {
  con <- file(path, "r")
  on.exit(close(con), add = TRUE)
  for (i in seq_len(max_lines)) {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (!length(ln)) break
    s <- trimws(ln)
    if (!nzchar(s)) next
    su <- toupper(s)
    return(grepl("^#?NEXUS\\b", su) || grepl("\\bBEGIN\\s+TREES\\b", su))
  }
  FALSE
}

read_named_newick <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- x[grepl("=", x, fixed = TRUE)]
  if (!length(x)) stop("No named-Newick lines found in: ", path)
  nm <- trimws(sub("\\s*=.*$", "", x))
  nw <- trimws(sub("^.*=\\s*", "", x))
  tr <- lapply(nw, function(s) ape::read.tree(text = s))
  names(tr) <- nm
  class(tr) <- "multiPhylo"
  tr
}

read_any_tree <- function(path) {
  if (!file.exists(path)) stop("Tree file not found: ", path)
  tr <- NULL
  if (is_nexus_file(path)) {
    tr <- try(ape::read.nexus(path), silent = TRUE)
    if (inherits(tr, "try-error")) tr <- try(ape::read.tree(path), silent = TRUE)
  } else {
    tr <- try(ape::read.tree(path), silent = TRUE)
    if (inherits(tr, "try-error")) tr <- try(ape::read.nexus(path), silent = TRUE)
  }
  if (inherits(tr, "try-error")) stop("Could not parse tree file: ", path)
  tr
}

load_tree_from_file <- function(path, tree_name = "", tree_index = 1L) {
  x <- try(read_named_newick(path), silent = TRUE)
  if (!inherits(x, "try-error")) {
    if (nzchar(tree_name)) {
      if (!(tree_name %in% names(x))) stop("Named tree not found: ", tree_name)
      return(list(tree = x[[tree_name]], tree_id = tree_name))
    }
    if (tree_index < 1L || tree_index > length(x)) stop("Tree index out of range.")
    return(list(tree = x[[tree_index]], tree_id = names(x)[tree_index]))
  }

  y <- read_any_tree(path)
  if (inherits(y, "multiPhylo")) {
    if (nzchar(tree_name) && !is.null(names(y)) && tree_name %in% names(y)) {
      return(list(tree = y[[tree_name]], tree_id = tree_name))
    }
    if (tree_index < 1L || tree_index > length(y)) stop("Tree index out of range for multiPhylo.")
    tree_id <- if (!is.null(names(y)) && nzchar(names(y)[tree_index])) names(y)[tree_index] else paste0("tree_", tree_index)
    return(list(tree = y[[tree_index]], tree_id = tree_id))
  }
  list(tree = y, tree_id = tools::file_path_sans_ext(basename(path)))
}

normalize_calibration_df <- function(df, calibration_tag = NULL) {
  names(df) <- tolower(names(df))
  if (!("taxona" %in% names(df)) && "taxon_a" %in% names(df)) df$taxona <- df$taxon_a
  if (!("taxonb" %in% names(df)) && "taxon_b" %in% names(df)) df$taxonb <- df$taxon_b
  if (!("age_min" %in% names(df)) && "minage" %in% names(df)) df$age_min <- df$minage
  if (!("age_max" %in% names(df)) && "maxage" %in% names(df)) df$age_max <- df$maxage
  if (!("age_min" %in% names(df)) && "age" %in% names(df)) df$age_min <- df$age
  if (!("age_max" %in% names(df)) && "age" %in% names(df)) df$age_max <- df$age
  if (!("age_max" %in% names(df))) df$age_max <- Inf

  need_cal <- c("taxona", "taxonb", "age_min")
  miss_cal <- setdiff(need_cal, names(df))
  if (length(miss_cal)) {
    stop("Calibration CSV missing required columns: ", paste(miss_cal, collapse = ", "))
  }

  if ("candidate" %in% names(df)) {
    tag <- trimws(as.character(df$candidate))
    specific <- sort(unique(tag[nzchar(tag) & !(tolower(tag) %in% c("all", "na"))]))
    if (!is.null(calibration_tag) && nzchar(calibration_tag)) {
      keep <- !nzchar(tag) | tolower(tag) == "all" | tag == calibration_tag
      df <- df[keep, , drop = FALSE]
    } else if (length(specific) > 1L) {
      stop(
        "Calibration CSV contains candidate-specific rows. Pass --calibration-tag=<name> ",
        "or provide a method-agnostic calibration CSV."
      )
    }
  }

  out <- data.frame(
    taxonA = as.character(df$taxona),
    taxonB = as.character(df$taxonb),
    age_min = as.numeric(df$age_min),
    age_max = as.numeric(df$age_max),
    calibration_source = "csv",
    stringsAsFactors = FALSE
  )
  out$age_max[is.na(out$age_max)] <- Inf
  out <- out[is.finite(out$age_min) & !is.na(out$age_max), , drop = FALSE]
  if (!nrow(out)) stop("No usable calibration rows after normalization.")
  out
}

build_calib_pairs_congruif <- function(target_tree, ref_tree) {
  if (!requireNamespace("geiger", quietly = TRUE)) {
    stop("Package geiger is required for congruification runs.")
  }
  if (!ape::is.ultrametric(ref_tree)) {
    if (!requireNamespace("phytools", quietly = TRUE)) {
      stop("Package phytools is required to force the reference tree ultrametric for congruification.")
    }
    ref_tree <- phytools::force.ultrametric(ref_tree, method = "extend")
  }
  cg <- geiger::congruify.phylo(ref_tree, target_tree)
  cc <- cg$calibrations
  if (is.null(cc) || !nrow(cc)) stop("Congruification returned zero calibrations.")
  out <- data.frame(
    taxonA = as.character(cc$taxonA),
    taxonB = as.character(cc$taxonB),
    age_min = as.numeric(cc$MinAge),
    age_max = as.numeric(cc$MaxAge),
    calibration_source = "congruify",
    stringsAsFactors = FALSE
  )
  out$age_max[is.na(out$age_max)] <- Inf
  out <- out[is.finite(out$age_min) & !is.na(out$age_max), , drop = FALSE]
  if (!nrow(out)) stop("Congruification produced no usable calibration rows.")
  out
}

build_children_lookup <- function(phy) {
  total <- ape::Ntip(phy) + phy$Nnode
  ch <- vector("list", total)
  for (k in seq_len(nrow(phy$edge))) {
    ch[[phy$edge[k, 1L]]] <- c(ch[[phy$edge[k, 1L]]], phy$edge[k, 2L])
  }
  ch
}

descendant_tip_labels <- function(phy, node, children = NULL) {
  if (is.null(children)) children <- build_children_lookup(phy)
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

select_root_pair <- function(phy) {
  if (!ape::is.rooted(phy)) stop("The target phylogram must be rooted to inject a root-age calibration.")
  root_node <- ape::Ntip(phy) + 1L
  children <- build_children_lookup(phy)
  kids <- children[[root_node]]
  if (length(kids) < 2L) stop("Could not identify two root children for the root-age calibration.")
  child_min_tips <- vapply(kids, function(nd) descendant_tip_labels(phy, nd, children)[1L], FUN.VALUE = character(1))
  kids <- kids[order(child_min_tips)]
  c(
    descendant_tip_labels(phy, kids[1L], children)[1L],
    descendant_tip_labels(phy, kids[2L], children)[1L]
  )
}

append_root_calibration <- function(pair_df, phy, root_age) {
  pair <- select_root_pair(phy)
  rbind(
    pair_df,
    data.frame(
      taxonA = pair[1L],
      taxonB = pair[2L],
      age_min = as.numeric(root_age),
      age_max = as.numeric(root_age),
      calibration_source = "root_age_arg",
      stringsAsFactors = FALSE
    )
  )
}

map_pairwise_calibrations_to_nodes <- function(phy, pair_df) {
  mapped_rows <- vector("list", nrow(pair_df))
  missing_taxa <- 0L
  bad_mrca <- 0L
  kept <- 0L

  for (i in seq_len(nrow(pair_df))) {
    a <- pair_df$taxonA[i]
    b <- pair_df$taxonB[i]
    if (!(a %in% phy$tip.label) || !(b %in% phy$tip.label)) {
      missing_taxa <- missing_taxa + 1L
      next
    }
    nd <- ape::getMRCA(phy, c(a, b))
    if (is.null(nd) || !is.finite(nd)) {
      bad_mrca <- bad_mrca + 1L
      next
    }
    kept <- kept + 1L
    mapped_rows[[kept]] <- data.frame(
      pair_index = i,
      node = as.integer(nd),
      taxonA = a,
      taxonB = b,
      age_min = as.numeric(pair_df$age_min[i]),
      age_max = as.numeric(pair_df$age_max[i]),
      calibration_source = as.character(pair_df$calibration_source[i] %||% "unknown"),
      stringsAsFactors = FALSE
    )
  }

  mapped_rows <- mapped_rows[seq_len(kept)]
  if (!length(mapped_rows)) {
    stop("No calibration rows mapped onto the target phylogram.")
  }

  mapped_df <- do.call(rbind, mapped_rows)
  merged_rows <- list()
  dropped_rows <- list()
  keep_n <- 0L
  drop_n <- 0L

  for (nd in sort(unique(mapped_df$node))) {
    z <- mapped_df[mapped_df$node == nd, , drop = FALSE]
    mn <- max(z$age_min, na.rm = TRUE)
    mx <- min(z$age_max, na.rm = TRUE)
    if (is.finite(mn) && !is.na(mx) && mn <= mx) {
      keep_n <- keep_n + 1L
      merged_rows[[keep_n]] <- data.frame(
        node = as.integer(nd),
        taxonA = z$taxonA[1L],
        taxonB = z$taxonB[1L],
        age_min = mn,
        age_max = mx,
        n_merged = nrow(z),
        pair_indices = paste(z$pair_index, collapse = ";"),
        calibration_sources = paste(unique(z$calibration_source), collapse = ";"),
        stringsAsFactors = FALSE
      )
    } else {
      drop_n <- drop_n + 1L
      dropped_rows[[drop_n]] <- data.frame(
        node = as.integer(nd),
        age_min_max = mn,
        age_max_min = mx,
        n_conflicting = nrow(z),
        pair_indices = paste(z$pair_index, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }

  merged_df <- if (keep_n) do.call(rbind, merged_rows) else data.frame()
  dropped_df <- if (drop_n) do.call(rbind, dropped_rows) else data.frame()
  if (!nrow(merged_df)) {
    stop("All mapped calibration nodes became inconsistent after duplicate-node intersection.")
  }

  list(
    mapped = mapped_df,
    merged = merged_df,
    dropped_inconsistent = dropped_df,
    dropped_missing_taxa = missing_taxa,
    dropped_bad_mrca = bad_mrca
  )
}

build_chronos_calib_from_node_bounds <- function(phy, node_bounds) {
  ## Cap Inf age_max to a safe finite value so older ape versions don't choke.
  ## 10x the largest finite age bound is large enough to be effectively unconstrained.
  age_max <- node_bounds$age_max
  cap <- 10 * max(c(node_bounds$age_min, age_max[is.finite(age_max)]), na.rm = TRUE)
  age_max[!is.finite(age_max)] <- cap
  ape::makeChronosCalib(
    phy,
    node = node_bounds$node,
    age.min = node_bounds$age_min,
    age.max = age_max
  )
}

fit_score_chronos <- function(tr) {
  ph <- attr(tr, "PHIIC")
  if (is.list(ph) && is.finite(ph$PHIIC)) return(ph$PHIIC)
  pl <- attr(tr, "ploglik")
  if (is.finite(pl)) return(-pl)
  Inf
}

write_chronos_helper_script <- function(path) {
  lines <- c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "params <- readRDS(args[1])",
    "out_file <- args[2]",
    "suppressPackageStartupMessages(library(ape))",
    "`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a",
    "if (!is.null(params$seed) && is.finite(params$seed)) set.seed(as.integer(params$seed))",
    "ctrl <- ape::chronos.control()",
    "ctrl$iter.max <- as.integer(params$iter_max)",
    "ctrl$eval.max <- as.integer(params$iter_max)",
    "ctrl$tol <- as.numeric(params$tol)",
    "if (identical(params$model, 'discrete') && is.finite(params$nb_rate_cat)) {",
    "  ctrl$nb.rate.cat <- as.integer(params$nb_rate_cat)",
    "}",
    "tr <- try(",
    "  ape::chronos(",
    "    params$phy,",
    "    lambda = params$lambda,",
    "    model = params$model,",
    "    calibration = params$calib,",
    "    quiet = TRUE,",
    "    control = ctrl",
    "  ),",
    "  silent = TRUE",
    ")",
    "if (inherits(tr, 'try-error')) {",
    "  res <- list(ok = FALSE, error = conditionMessage(attr(tr, 'condition') %||% simpleError(as.character(tr))))",
    "} else if (!inherits(tr, 'phylo') || is.null(tr$edge.length) || any(!is.finite(tr$edge.length)) || any(tr$edge.length < 0)) {",
    "  res <- list(ok = FALSE, error = 'chronos returned an invalid tree.')",
    "} else {",
    "  ph <- attr(tr, 'PHIIC')",
    "  phiic <- if (is.list(ph) && is.finite(ph$PHIIC)) ph$PHIIC else NA_real_",
    "  ploglik <- attr(tr, 'ploglik') %||% NA_real_",
    "  fit_score <- if (is.finite(phiic)) phiic else Inf",
    "  res <- list(ok = TRUE, tree = tr, phiic = phiic, ploglik = ploglik, fit_score = fit_score)",
    "}",
    "saveRDS(res, out_file)"
  )
  writeLines(lines, path)
}

run_chronos_subprocess <- function(phy, calib, model, lambda, nb_rate_cat, iter_max, tol, timeout, attempt) {
  helper_file <- tempfile(pattern = "chronos_attempt_", fileext = ".R")
  input_file <- tempfile(pattern = "chronos_attempt_", fileext = ".rds")
  output_file <- tempfile(pattern = "chronos_attempt_", fileext = ".rds")
  on.exit(unlink(c(helper_file, input_file, output_file), force = TRUE), add = TRUE)

  write_chronos_helper_script(helper_file)
  saveRDS(list(
    phy = phy,
    calib = calib,
    model = model,
    lambda = lambda,
    nb_rate_cat = nb_rate_cat,
    iter_max = iter_max,
    tol = tol,
    seed = as.integer((as.numeric(Sys.time()) * 1e6 + Sys.getpid() + attempt) %% .Machine$integer.max)
  ), input_file)

  timeout_bin <- Sys.which("timeout")
  if (!nzchar(timeout_bin)) timeout_bin <- Sys.which("gtimeout")
  if (!nzchar(timeout_bin)) {
    stop("A timeout-enabled chronos attempt was requested, but no timeout binary was found.")
  }

  args <- c(as.character(ceiling(timeout)), "Rscript", helper_file, input_file, output_file)
  suppressWarnings(system2(timeout_bin, args = args, stdout = TRUE, stderr = TRUE))

  if (file.exists(output_file)) {
    return(readRDS(output_file))
  }

  list(ok = FALSE, error = paste0("chronos attempt timed out after ", ceiling(timeout), " seconds."))
}

run_chronos_one <- function(phy, calib, model, lambda, nb_rate_cat = NA_integer_,
                            retries = 2L, iter_max = 10000L, tol = 1e-8,
                            attempt_timeout = 0) {
  best <- NULL
  last_error <- NA_character_
  retries <- max(1L, as.integer(retries))
  timeout <- suppressWarnings(as.numeric(attempt_timeout))
  if (!is.finite(timeout) || timeout <= 0) timeout <- 0

  for (attempt in seq_len(retries)) {
    ctrl <- ape::chronos.control()
    ctrl$iter.max <- as.integer(iter_max)
    ctrl$eval.max <- as.integer(iter_max)
    ctrl$tol <- as.numeric(tol)
    if (identical(model, "discrete") && is.finite(nb_rate_cat)) {
      ctrl$nb.rate.cat <- as.integer(nb_rate_cat)
    }
    if (timeout > 0) {
      subres <- run_chronos_subprocess(
        phy = phy,
        calib = calib,
        model = model,
        lambda = lambda,
        nb_rate_cat = nb_rate_cat,
        iter_max = iter_max,
        tol = tol,
        timeout = timeout,
        attempt = attempt
      )
      if (!isTRUE(subres$ok)) {
        last_error <- subres$error %||% paste0("chronos attempt ", attempt, " failed.")
        next
      }
      res <- list(
        tree = subres$tree,
        phiic = subres$phiic,
        ploglik = subres$ploglik,
        fit_score = subres$fit_score,
        attempt = attempt
      )
    } else {
      tr <- try(
        ape::chronos(
          phy,
          lambda = lambda,
          model = model,
          calibration = calib,
          quiet = TRUE,
          control = ctrl
        ),
        silent = TRUE
      )
      if (inherits(tr, "try-error")) {
        last_error <- conditionMessage(attr(tr, "condition") %||% simpleError(as.character(tr)))
        next
      }
      if (!inherits(tr, "phylo") || is.null(tr$edge.length) ||
          any(!is.finite(tr$edge.length)) || any(tr$edge.length < 0)) {
        last_error <- "chronos returned an invalid tree."
        next
      }

      res <- list(
        tree = tr,
        phiic = if (is.list(attr(tr, "PHIIC"))) attr(tr, "PHIIC")$PHIIC else NA_real_,
        ploglik = attr(tr, "ploglik") %||% NA_real_,
        fit_score = fit_score_chronos(tr),
        attempt = attempt
      )
    }
    if (is.null(best) || res$fit_score < best$fit_score) best <- res
  }

  if (is.null(best)) {
    return(list(
      status = "FAILED",
      error = last_error %||% "chronos failed without a reported error.",
      tree = NULL,
      phiic = NA_real_,
      ploglik = NA_real_,
      fit_score = NA_real_,
      attempt = NA_integer_
    ))
  }

  c(best, list(status = "OK", error = NA_character_))
}

tail_text <- function(path, n = 12L) {
  if (!file.exists(path)) return(NA_character_)
  txt <- readLines(path, warn = FALSE)
  if (!length(txt)) return(NA_character_)
  paste(tail(txt, n), collapse = "\n")
}

parse_treepl_final_objective <- function(log_path) {
  if (!file.exists(log_path)) return(NA_real_)
  txt <- readLines(log_path, warn = FALSE)
  for (pat in c("after opt calc2:", "after opt calc1:", "exit siman:")) {
    hits <- grep(pat, txt, fixed = TRUE, value = TRUE)
    if (!length(hits)) next
    val <- suppressWarnings(as.numeric(trimws(sub("^.*:\\s*", "", tail(hits, 1L)))))
    if (is.finite(val)) return(val)
  }
  NA_real_
}

distance_from_one <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x) & x > 0, abs(log10(x)), Inf)
}

build_representative_candidates <- function(all_df) {
  reps <- list()
  add_i <- 0L

  chronos_df <- all_df[
    all_df$method == "chronos" &
      all_df$status == "OK" &
      nzchar(all_df$tree_file),
    ,
    drop = FALSE
  ]
  if (nrow(chronos_df)) {
    for (mdl in unique(as.character(chronos_df$model))) {
      z <- chronos_df[chronos_df$model == mdl, , drop = FALSE]
      if (!nrow(z)) next
      ord <- order(
        ifelse(is.finite(z$fit_score), z$fit_score, Inf),
        distance_from_one(z$lambda),
        ifelse(is.finite(z$lambda), z$lambda, Inf),
        ifelse(is.finite(z$nb_rate_cat), z$nb_rate_cat, Inf),
        z$candidate
      )
      best <- z[ord[1L], , drop = FALSE]
      add_i <- add_i + 1L
      ## Build descriptive candidate name including selected lambda and k
      cand_name <- paste0("chronos_", mdl)
      lam_str <- fmt_token(best$lambda[1L])
      cand_name <- paste0(cand_name, "_lambda", lam_str)
      if (identical(mdl, "discrete") && is.finite(best$nb_rate_cat[1L])) {
        cand_name <- paste0(cand_name, "_k", best$nb_rate_cat[1L])
      }
      reps[[add_i]] <- data.frame(
        candidate = cand_name,
        tree_file = best$tree_file[1L],
        method = "chronos",
        model = mdl,
        lambda = best$lambda[1L],
        nb_rate_cat = best$nb_rate_cat[1L],
        smooth = NA_real_,
        selected_from_candidate = best$candidate[1L],
        selection_score = best$fit_score[1L],
        selection_rule = "min_fit_score_then_lambda_closest_to_1",
        stringsAsFactors = FALSE
      )
    }
  }

  treepl_df <- all_df[
    all_df$method == "treePL" &
      grepl("^OK", all_df$status) &
      nzchar(all_df$tree_file),
    ,
    drop = FALSE
  ]
  if (nrow(treepl_df)) {
    if (!("treepl_objective" %in% names(treepl_df))) treepl_df$treepl_objective <- NA_real_
    missing_obj <- !is.finite(treepl_df$treepl_objective)
    if (any(missing_obj) && "log_file" %in% names(treepl_df)) {
      treepl_df$treepl_objective[missing_obj] <- vapply(
        treepl_df$log_file[missing_obj],
        parse_treepl_final_objective,
        FUN.VALUE = numeric(1)
      )
    }
    ord <- order(
      ifelse(is.finite(treepl_df$treepl_objective), treepl_df$treepl_objective, Inf),
      distance_from_one(treepl_df$smooth),
      ifelse(is.finite(treepl_df$smooth), treepl_df$smooth, Inf),
      treepl_df$candidate
    )
    best <- treepl_df[ord[1L], , drop = FALSE]
    add_i <- add_i + 1L
    ## Build descriptive candidate name including selected smoothing
    smooth_str <- fmt_token(best$smooth[1L])
    treepl_cand_name <- paste0("treePL_smooth", smooth_str)
    reps[[add_i]] <- data.frame(
      candidate = treepl_cand_name,
      tree_file = best$tree_file[1L],
      method = "treePL",
      model = "treePL",
      lambda = NA_real_,
      nb_rate_cat = NA_integer_,
      smooth = best$smooth[1L],
      selected_from_candidate = best$candidate[1L],
      selection_score = best$treepl_objective[1L],
      selection_rule = "min_treepl_objective_then_smoothing_closest_to_1",
      stringsAsFactors = FALSE
    )
  }

  rel_df <- all_df[
    all_df$method == "RelTime" &
      all_df$status == "OK" &
      nzchar(all_df$tree_file),
    ,
    drop = FALSE
  ]
  if (nrow(rel_df)) {
    best <- rel_df[1L, , drop = FALSE]
    add_i <- add_i + 1L
    reps[[add_i]] <- data.frame(
      candidate = "RelTime",
      tree_file = best$tree_file[1L],
      method = "RelTime",
      model = "RelTime",
      lambda = NA_real_,
      nb_rate_cat = NA_integer_,
      smooth = NA_real_,
      selected_from_candidate = best$candidate[1L],
      selection_score = NA_real_,
      selection_rule = "single_reltime_run",
      stringsAsFactors = FALSE
    )
  }

  mega_df <- all_df[
    all_df$method == "RelTime_MEGA" &
      all_df$status == "OK" &
      nzchar(all_df$tree_file),
    ,
    drop = FALSE
  ]
  if (nrow(mega_df)) {
    best <- mega_df[1L, , drop = FALSE]
    add_i <- add_i + 1L
    reps[[add_i]] <- data.frame(
      candidate = "RelTime_MEGA",
      tree_file = best$tree_file[1L],
      method = "RelTime_MEGA",
      model = "RelTime_MEGA",
      lambda = NA_real_,
      nb_rate_cat = NA_integer_,
      smooth = NA_real_,
      selected_from_candidate = best$candidate[1L],
      selection_score = NA_real_,
      selection_rule = "mega_reltime_run",
      stringsAsFactors = FALSE
    )
  }

  bind_rows_fill(reps)
}

copy_if_present <- function(src, dst) {
  if (!is.character(src) || length(src) != 1L || !nzchar(src) || !file.exists(src)) return(FALSE)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  file.copy(src, dst, overwrite = TRUE)
}

## In-memory version: takes phylo object + CI data.frame (with matching node numbers)
## and writes annotated NEXUS directly. Avoids write/read round-trip node reordering.
annotate_tree_with_ci_from_objects <- function(tr, ci_df, tree_dst, centered = TRUE) {
  if (!inherits(tr, "phylo") || !is.data.frame(ci_df)) return(FALSE)
  if (!all(c("node", "ci_lo", "ci_hi") %in% names(ci_df))) return(FALSE)
  n_tip <- ape::Ntip(tr)
  int_nodes <- seq.int(n_tip + 1L, n_tip + tr$Nnode)
  d <- ape::node.depth.edgelength(tr)
  root_age <- max(d, na.rm = TRUE)
  node_age <- root_age - d
  idx <- match(int_nodes, ci_df$node)
  lo <- ci_df$ci_lo[idx]
  hi <- ci_df$ci_hi[idx]
  age <- node_age[int_nodes]
  bad <- !is.finite(lo) | !is.finite(hi) | hi < lo
  lo[bad] <- age[bad]
  hi[bad] <- age[bad]

  ## Recenter CIs on point estimates for display (Paradis et al. 2023).
  ## Raw bootstrap percentile CIs can be off-center — the point estimate
  ## may sit near one edge of the interval. Recentering preserves the CI
  ## width but shifts the bar so the node age sits at the center.
  ## Raw (uncentered) values are kept in the CSV.
  if (isTRUE(centered)) {
    w95 <- pmax(0, hi - lo)
    lo <- age - (w95 / 2)
    hi <- age + (w95 / 2)
    lo <- pmax(0, lo)  ## ages can't go negative
  }

  ## 50% CI (IQR) — include if available
  has_q50 <- all(c("q50_lo", "q50_hi") %in% names(ci_df))
  if (has_q50) {
    q50_lo <- ci_df$q50_lo[idx]
    q50_hi <- ci_df$q50_hi[idx]
    bad50 <- !is.finite(q50_lo) | !is.finite(q50_hi) | q50_hi < q50_lo
    q50_lo[bad50] <- age[bad50]
    q50_hi[bad50] <- age[bad50]
    if (isTRUE(centered)) {
      w50 <- pmax(0, q50_hi - q50_lo)
      q50_lo <- age - (w50 / 2)
      q50_hi <- age + (w50 / 2)
      q50_lo <- pmax(0, q50_lo)
    }
  }

  tokens <- sprintf("CODXNODE%04d", seq_along(int_nodes))
  tr$node.label <- tokens
  txt <- ape::write.tree(tr)
  fmt <- function(x) format(as.numeric(x), scientific = FALSE, trim = TRUE, digits = 15L)
  for (i in seq_along(tokens)) {
    ann <- paste0("height_95%_HPD={", fmt(lo[i]), ",", fmt(hi[i]), "}")
    if (has_q50) {
      ann <- paste0(ann, ",height_50%_HPD={", fmt(q50_lo[i]), ",", fmt(q50_hi[i]), "}")
    }
    comment_i <- paste0("[&", ann, "]")
    txt <- sub(tokens[i], comment_i, txt, fixed = TRUE)
  }
  .write_figtree_nexus(txt, tree_dst)
  TRUE
}

## File-based version: reads tree + CI CSV, matches by descendant tip sets to
## avoid node-numbering mismatches from write.tree → read.tree round-trips.
annotate_tree_with_ci_comments <- function(tree_src, ci_src, tree_dst, centered = TRUE) {
  if (!is.character(tree_src) || length(tree_src) != 1L || !nzchar(tree_src) || !file.exists(tree_src)) return(FALSE)
  if (!is.character(ci_src) || length(ci_src) != 1L || !nzchar(ci_src) || !file.exists(ci_src)) return(FALSE)
  tr <- try(ape::read.tree(tree_src), silent = TRUE)
  ci <- try(read.csv(ci_src, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(tr, "try-error") || !inherits(tr, "phylo")) return(FALSE)
  if (inherits(ci, "try-error") || !all(c("node", "ci_lo", "ci_hi") %in% names(ci))) return(FALSE)
  n_tip <- ape::Ntip(tr)
  int_nodes <- seq.int(n_tip + 1L, n_tip + tr$Nnode)
  d <- ape::node.depth.edgelength(tr)
  root_age <- max(d, na.rm = TRUE)
  node_age <- root_age - d

  lo <- rep(NA_real_, length(int_nodes))
  hi <- rep(NA_real_, length(int_nodes))
  age <- node_age[int_nodes]
  ci_nodes <- suppressWarnings(as.integer(ci$node))
  exact_match_ok <- length(ci_nodes) == nrow(ci) &&
    any(is.finite(ci_nodes)) &&
    !anyDuplicated(ci_nodes[is.finite(ci_nodes)]) &&
    sum(ci_nodes %in% int_nodes, na.rm = TRUE) >= max(1L, floor(length(int_nodes) * 0.5))

  if (isTRUE(exact_match_ok)) {
    idx <- match(int_nodes, ci_nodes)
    lo <- ci$ci_lo[idx]
    hi <- ci$ci_hi[idx]
  } else if ("age" %in% names(ci)) {
    ## Fall back to closest-age greedy matching only when the CI table does
    ## not carry trustworthy node IDs for this read-back tree.
    ci_ages <- ci$age
    ## Sort tree nodes by age to get consistent greedy ordering
    ord <- order(age, decreasing = TRUE)
    ci_ord <- order(ci_ages, decreasing = TRUE)
    used <- logical(length(ci_ages))
    for (jj in seq_along(ord)) {
      j <- ord[jj]
      best_k <- NA_integer_; best_d <- Inf
      for (kk in seq_along(ci_ord)) {
        k <- ci_ord[kk]
        if (used[k]) next
        dd <- abs(age[j] - ci_ages[k])
        if (dd < best_d) { best_d <- dd; best_k <- k }
        if (dd < 1e-8) break  # exact match
      }
      if (!is.na(best_k) && best_d < 1e-4) {
        lo[j] <- ci$ci_lo[best_k]
        hi[j] <- ci$ci_hi[best_k]
        used[best_k] <- TRUE
      }
    }
  } else {
    return(FALSE)
  }
  bad <- !is.finite(lo) | !is.finite(hi) | hi < lo
  lo[bad] <- age[bad]
  hi[bad] <- age[bad]

  ## Recenter CIs on point estimates for display (same logic as object-based version)
  if (isTRUE(centered)) {
    w95 <- pmax(0, hi - lo)
    lo <- age - (w95 / 2)
    hi <- age + (w95 / 2)
    lo <- pmax(0, lo)
  }

  ## 50% CI (IQR) — include if available
  has_q50 <- all(c("q50_lo", "q50_hi") %in% names(ci))
  if (has_q50) {
    q50_lo <- rep(NA_real_, length(int_nodes))
    q50_hi <- rep(NA_real_, length(int_nodes))
    if (isTRUE(exact_match_ok)) {
      q50_lo <- ci$q50_lo[idx]
      q50_hi <- ci$q50_hi[idx]
    }
    bad50 <- !is.finite(q50_lo) | !is.finite(q50_hi) | q50_hi < q50_lo
    q50_lo[bad50] <- age[bad50]
    q50_hi[bad50] <- age[bad50]
    if (isTRUE(centered)) {
      w50 <- pmax(0, q50_hi - q50_lo)
      q50_lo <- age - (w50 / 2)
      q50_hi <- age + (w50 / 2)
      q50_lo <- pmax(0, q50_lo)
    }
  }

  tokens <- sprintf("CODXNODE%04d", seq_along(int_nodes))
  tr$node.label <- tokens
  txt <- ape::write.tree(tr)
  fmt <- function(x) format(as.numeric(x), scientific = FALSE, trim = TRUE, digits = 15L)
  for (i in seq_along(tokens)) {
    ann <- paste0("height_95%_HPD={", fmt(lo[i]), ",", fmt(hi[i]), "}")
    if (has_q50) {
      ann <- paste0(ann, ",height_50%_HPD={", fmt(q50_lo[i]), ",", fmt(q50_hi[i]), "}")
    }
    comment_i <- paste0("[&", ann, "]")
    txt <- sub(tokens[i], comment_i, txt, fixed = TRUE)
  }
  .write_figtree_nexus(txt, tree_dst)
  TRUE
}

.write_figtree_nexus <- function(newick_txt, dst) {
  nexus_txt <- c(
    "#NEXUS",
    "",
    "begin trees;",
    paste0("\ttree tree_1 = [&R] ", newick_txt),
    "end;",
    "",
    "begin figtree;",
    "\tset branchLabels.displayAttribute=\"Branch times\";",
    "\tset branchLabels.isShown=false;",
    "\tset layout.layoutType=\"RECTILINEAR\";",
    "\tset nodeBars.barWidth=4.0;",
    "\tset nodeBars.displayAttribute=\"height_95%_HPD\";",
    "\tset nodeBars.isShown=true;",
    "\tset nodeLabels.displayAttribute=\"Node ages\";",
    "\tset nodeLabels.isShown=false;",
    "\tset scaleBar.isShown=true;",
    "\tset tipLabels.displayAttribute=\"Names\";",
    "\tset tipLabels.isShown=true;",
    "\tset trees.rooting=true;",
    "\tset trees.rootingType=\"User Selection\";",
    "\tset trees.transform=false;",
    "end;"
  )
  writeLines(nexus_txt, dst)
}

materialize_key_outputs <- function(outdir, candidates_df, fit_files = character(0)) {
  main_tree_dir <- file.path(outdir, "MAIN_OUTPUT_TREES")
  model_fit_dir <- file.path(outdir, "model_fits")
  dir.create(main_tree_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(model_fit_dir, recursive = TRUE, showWarnings = FALSE)

  if (nrow(candidates_df)) {
    for (i in seq_len(nrow(candidates_df))) {
      src <- candidates_df$tree_file[i]
      ext <- tools::file_ext(src)
      ext <- if (nzchar(ext)) paste0(".", ext) else ".tre"
      dst <- file.path(main_tree_dir, paste0(candidates_df$candidate[i], ext))
      copy_if_present(src, dst)
      candidates_df$tree_file[i] <- dst
    }
  }

  if (length(fit_files)) {
    for (nm in names(fit_files)) {
      src <- fit_files[[nm]]
      if (!is.character(src) || length(src) != 1L || !nzchar(src) || !file.exists(src)) next
      ext <- tools::file_ext(src)
      ext <- if (nzchar(ext)) paste0(".", ext) else ""
      dst <- file.path(model_fit_dir, paste0(nm, ext))
      copy_if_present(src, dst)
    }
  }

  list(
    main_tree_dir = main_tree_dir,
    model_fit_dir = model_fit_dir,
    candidates_df = candidates_df
  )
}

build_representative_uncertainty <- function(candidates_df, uncertainty_df) {
  if (!nrow(candidates_df)) return(data.frame())
  if (is.null(uncertainty_df) || !nrow(uncertainty_df) || !("candidate" %in% names(uncertainty_df))) {
    return(do.call(rbind, lapply(
      candidates_df$candidate,
      make_unscored_uncertainty_row
    )))
  }
  out <- vector("list", nrow(candidates_df))
  for (i in seq_len(nrow(candidates_df))) {
    src <- candidates_df$selected_from_candidate[i] %||% candidates_df$candidate[i]
    hit <- uncertainty_df[uncertainty_df$candidate == src, , drop = FALSE]
    if (!nrow(hit)) {
      out[[i]] <- make_unscored_uncertainty_row(candidates_df$candidate[i])
      next
    }
    row_i <- hit[1L, , drop = FALSE]
    row_i$candidate <- candidates_df$candidate[i]
    out[[i]] <- row_i
  }
  bind_rows_fill(out)
}

copy_representative_ci_files <- function(outdir, candidates_df, all_df) {
  if (!nrow(candidates_df) || is.null(all_df) || !nrow(all_df)) {
    return(list(ci_dir = "", main_tree_dir = ""))
  }
  main_dir <- file.path(outdir, "MAIN_OUTPUT_TREES")
  ci_dir   <- file.path(main_dir, "ci")
  dir.create(ci_dir, recursive = TRUE, showWarnings = FALSE)
  copied_ci <- FALSE

  for (i in seq_len(nrow(candidates_df))) {
    src_id <- candidates_df$selected_from_candidate[i] %||% candidates_df$candidate[i]
    hit <- all_df[all_df$candidate == src_id, , drop = FALSE]
    if (!nrow(hit)) next
    tree_src <- candidates_df$tree_file[i]
    cand_name <- candidates_df$candidate[i]
    ci_method <- candidates_df$method[i]

    ## -- Bootstrap / main CI file (skip RelTime — handled below) --
    if (!(ci_method %in% c("RelTime", "RelTime_MEGA")) &&
        "ci_file" %in% names(hit) && nzchar(hit$ci_file[1L]) && file.exists(hit$ci_file[1L])) {
      ## CSV into ci/ subfolder
      csv_dst <- file.path(ci_dir, paste0(cand_name, "_ci.csv"))
      copy_if_present(hit$ci_file[1L], csv_dst)
      copied_ci <- TRUE

      ## Overwrite the plain tree in MAIN_OUTPUT_TREES with CI-annotated version
      ci_tree_src <- file.path(outdir, "selected_ci",
                               paste0(prefix, "_", cand_name, "_with_CI.tre"))
      if (file.exists(ci_tree_src)) {
        main_dst <- file.path(main_dir, paste0(cand_name, "_with_CI.tre"))
        copy_if_present(ci_tree_src, main_dst)
        ## Remove the plain (non-CI) tree to avoid confusion
        plain_dst <- file.path(main_dir, paste0(cand_name, ".tre"))
        if (file.exists(plain_dst)) file.remove(plain_dst)
      }
    }

    ## -- Tao analytical CI file (RelTime R and MEGA) --
    if ("tao_ci_file" %in% names(hit) && nzchar(hit$tao_ci_file[1L]) && file.exists(hit$tao_ci_file[1L])) {
      csv_dst <- file.path(ci_dir, paste0(cand_name, "_tao_ci.csv"))
      copy_if_present(hit$tao_ci_file[1L], csv_dst)
      copied_ci <- TRUE
    }

    ## -- RelTime R: copy both CI-annotated trees from reltime/ dir --
    if (identical(ci_method, "RelTime")) {
      rt_dir <- file.path(outdir, "reltime")
      boot_ci_tre <- file.path(rt_dir, paste0(prefix, "_RelTime_with_bootstrap_CI.tre"))
      tao_ci_tre  <- file.path(rt_dir, paste0(prefix, "_RelTime_with_tao_CI.tre"))
      if (file.exists(boot_ci_tre)) {
        copy_if_present(boot_ci_tre, file.path(main_dir, "RelTime_R_with_bootstrap_CI.tre"))
      }
      if (file.exists(tao_ci_tre)) {
        copy_if_present(tao_ci_tre, file.path(main_dir, "RelTime_R_with_tao_CI.tre"))
      }
      ## Remove the plain RelTime tree
      plain_dst <- file.path(main_dir, paste0(cand_name, ".tre"))
      if (file.exists(plain_dst)) file.remove(plain_dst)
    }

    ## -- RelTime MEGA: copy Tao CI-annotated tree --
    if (identical(ci_method, "RelTime_MEGA")) {
      mega_dir <- file.path(outdir, "reltime_mega")
      mega_ci_tre <- file.path(mega_dir, paste0(prefix, "_RelTime_MEGA_with_tao_CI.tre"))
      if (file.exists(mega_ci_tre)) {
        copy_if_present(mega_ci_tre, file.path(main_dir, "RelTime_MEGA_with_tao_CI.tre"))
      }
      plain_dst <- file.path(main_dir, paste0(cand_name, ".tre"))
      if (file.exists(plain_dst)) file.remove(plain_dst)
    }
  }

  list(
    ci_dir = if (copied_ci) ci_dir else "",
    main_tree_dir = main_dir
  )
}

resolve_treepl_bin <- function(user_bin, repo_dir) {
  candidates <- c(
    user_bin %||% NA_character_,
    Sys.getenv("TREEPL_BIN", unset = ""),
    Sys.which("treePL"),
    file.path(dirname(repo_dir), "treePL")
  )
  candidates <- unique(trimws(candidates))
  candidates <- candidates[nzchar(candidates)]
  for (cand in candidates) {
    if (file.exists(cand)) return(normalizePath(cand, winslash = "/", mustWork = TRUE))
  }
  ""
}

## ---------------------------------------------------------------------------
## MEGA-CC RelTime helpers
## ---------------------------------------------------------------------------

resolve_megacc_bin <- function(user_bin) {
  candidates <- c(
    user_bin %||% NA_character_,
    Sys.getenv("MEGACC_BIN", unset = ""),
    Sys.which("megacc")
  )
  candidates <- unique(trimws(candidates))
  candidates <- candidates[nzchar(candidates)]
  for (cand in candidates) {
    if (file.exists(cand)) return(normalizePath(cand, winslash = "/", mustWork = TRUE))
  }
  ""
}

write_mega_mao <- function(mao_path) {
  lines <- c(
    "; MEGA RelTime from branch lengths settings",
    "[ MEGAinfo ]",
    "ver=0",
    "[ DataSettings ]",
    "datatype=Newick Tree",
    "[ ProcessTypes ]",
    "ppRelTimeBLens=true",
    "[ AnalysisSettings ]"
  )
  writeLines(lines, mao_path)
}

write_mega_outgroup <- function(outgroup_path, phy) {
  ## Identify outgroup: tips descended from the smaller root child
  root_node <- Ntip(phy) + 1L
  kids <- phy$edge[phy$edge[, 1] == root_node, 2]
  count_tips <- function(nd) {
    if (nd <= Ntip(phy)) return(1L)
    length(extract.clade(phy, nd)$tip.label)
  }
  n_tips <- vapply(kids, count_tips, integer(1))
  ## Outgroup = side with fewer taxa
  og_kid <- kids[which.min(n_tips)]
  if (og_kid <= Ntip(phy)) {
    og_tips <- phy$tip.label[og_kid]
  } else {
    og_tips <- extract.clade(phy, og_kid)$tip.label
  }
  writeLines(paste0(og_tips, "=outgroup"), outgroup_path)
}

write_mega_calibrations <- function(cal_path, node_bounds, root_node = NULL) {
  ## MEGA RelTime cannot calibrate the root node — only internal nodes below
  ## the root.  Exclude the root calibration if present.
  keep <- rep(TRUE, nrow(node_bounds))
  if (!is.null(root_node)) {
    keep <- node_bounds$node != root_node
    n_dropped <- sum(!keep)
    if (n_dropped > 0L)
      message("Note: excluded ", n_dropped,
              " root-node calibration(s) — MEGA RelTime cannot constrain the root.")
  }
  nb <- node_bounds[keep, , drop = FALSE]
  lines <- character(nrow(nb))
  for (i in seq_len(nrow(nb))) {
    tag <- paste0("cal_", nb$node[i])
    lines[i] <- sprintf(
      "!MRCA='%s' TaxonA='%s' TaxonB='%s' MinTime=%.11f MaxTime=%.11f calibrationName='%s';",
      tag, nb$taxonA[i], nb$taxonB[i],
      nb$age_min[i], nb$age_max[i], tag
    )
  }
  writeLines(lines, cal_path)
}

write_mega_input_tree <- function(tree_path, phy) {
  ## MEGA-CC needs underscores in tip labels replaced with spaces, wrapped in quotes
  ## Actually: MEGA replaces underscores with spaces internally, so we write as-is
  ## and it reads them fine.  Write plain Newick.
  ape::write.tree(phy, file = tree_path)
}

run_megacc_reltime <- function(megacc_bin, phy, node_bounds, work_dir, prefix) {
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  mao_path      <- file.path(work_dir, "reltime_blens.mao")
  tree_path     <- file.path(work_dir, paste0(prefix, "_phylogram_for_mega.nwk"))
  outgroup_path <- file.path(work_dir, "outgroup.txt")
  cal_path      <- file.path(work_dir, "mega_calibrations.txt")
  out_prefix    <- file.path(work_dir, paste0(prefix, "_mega_reltime"))
  log_path      <- file.path(work_dir, paste0(prefix, "_mega_reltime.log"))

  ## Graft a mock outgroup tip sister to the entire ingroup so MEGA can
  ## treat it as the single outgroup.  This gives the root node a proper
  ## outgroup branch, avoiding the near-polytomy collapse at the root.
  mock_tip <- "MOCK_OUTGROUP_FOR_MEGA"
  d <- node.depth.edgelength(phy)
  root_depth <- max(d)
  nwk <- sub(";$", "", ape::write.tree(phy))
  new_nwk <- paste0("(", nwk, ":0.001,", mock_tip, ":",
                    format(root_depth + 0.001, scientific = FALSE), ");")
  phy_mega <- ape::read.tree(text = new_nwk)

  write_mega_mao(mao_path)
  write_mega_input_tree(tree_path, phy_mega)
  ## Outgroup = the mock tip only
  writeLines(paste0(mock_tip, "=outgroup"), outgroup_path)
  ## Write calibrations using original phy node numbering; the MRCA pairs
  ## reference tip names that exist in phy_mega too, so MEGA can find them.
  write_mega_calibrations(cal_path, node_bounds)

  ## Run megacc
  exit_code <- suppressWarnings(system2(
    megacc_bin,
    args = c("-a", mao_path, "-t", tree_path,
             "-c", cal_path, "-g", outgroup_path,
             "-o", out_prefix),
    stdout = log_path,
    stderr = log_path
  ))

  ## Parse output tree and drop the mock outgroup
  exact_tree_file <- paste0(out_prefix, "_exactTimes.nwk")
  if (file.exists(exact_tree_file)) {
    mega_tree <- ape::read.tree(exact_tree_file)
    ## Fix MEGA's label mangling: strip quotes and restore underscores
    mega_tree$tip.label <- gsub("^'|'$", "", mega_tree$tip.label)
    mega_tree$tip.label <- gsub(" ", "_", mega_tree$tip.label)
    ## Drop mock outgroup
    mock_idx <- match(mock_tip, mega_tree$tip.label)
    if (!is.na(mock_idx)) {
      mega_tree <- ape::drop.tip(mega_tree, mock_tip)
    }
    ## Write cleaned tree
    out_tree <- file.path(work_dir, paste0(prefix, "_RelTime_MEGA.tre"))
    ape::write.tree(mega_tree, file = out_tree)
    list(tree = mega_tree, tree_file = out_tree, status = "OK", error = NA_character_)
  } else {
    log_txt <- if (file.exists(log_path)) paste(readLines(log_path, warn = FALSE), collapse = "\n") else ""
    list(tree = NULL, tree_file = "", status = "FAILED",
         error = paste0("megacc exit code ", exit_code, "; ", log_txt))
  }
}

megacc_tao_ci <- function(megacc_bin, phy, node_bounds, work_dir, prefix,
                          n_sites = 1000L) {
  ## Run base MEGA RelTime
  base <- run_megacc_reltime(megacc_bin, phy, node_bounds, work_dir, prefix)
  if (base$status != "OK") return(base)

  mega_tree <- base$tree
  n_tip  <- Ntip(mega_tree)
  n_node <- Nnode(mega_tree)
  int_nodes <- seq.int(n_tip + 1L, n_tip + n_node)

  ## ---- Parse native Tao 95% CIs from MEGA's NEXUS output ----
  nexus_files <- list.files(work_dir, pattern = "_nexus\\.tre$", full.names = TRUE)
  if (!length(nexus_files)) {
    warning("No MEGA NEXUS output found — skipping Tao CIs")
    base$ci <- NULL
    return(base)
  }
  nexus_txt <- paste(readLines(nexus_files[1L], warn = FALSE), collapse = "\n")
  tree_line <- grep("^tree TREE1", strsplit(nexus_txt, "\n")[[1L]], value = TRUE)
  if (!length(tree_line)) {
    warning("No 'tree TREE1' line in MEGA NEXUS — skipping Tao CIs")
    base$ci <- NULL
    return(base)
  }

  ## Extract divtime and divtime_95%_CI for each internal node
  ## MEGA annotates internal nodes as: [&rate=...,divtime=T,divtime_95%_CI={hi,lo}]
  ## The NEXUS tree has numbered tips (1..N) and annotations on internal nodes.
  ## We parse the annotations in the order they appear (postorder traversal),
  ## then match to the cleaned tree via a reconstructed mapping.

  ci_raw <- regmatches(tree_line, gregexpr("divtime_95%_CI=[{][^}]+[}]", tree_line))[[1]]
  dt_raw <- regmatches(tree_line, gregexpr("divtime=([0-9eE.+-]+)", tree_line))[[1]]
  dt_vals <- as.numeric(gsub("divtime=", "", dt_raw))

  if (length(ci_raw) != length(dt_raw) || length(ci_raw) == 0L) {
    warning("Could not parse MEGA Tao CIs (", length(ci_raw), " CIs, ",
            length(dt_raw), " divtimes)")
    base$ci <- NULL
    return(base)
  }

  ## Parse CI bounds
  ci_parsed <- gsub("divtime_95%_CI=[{]|[}]", "", ci_raw)
  ci_split  <- strsplit(ci_parsed, ",")
  ci_lo <- vapply(ci_split, function(x) min(as.numeric(x)), numeric(1))
  ci_hi <- vapply(ci_split, function(x) max(as.numeric(x)), numeric(1))

  ## MEGA node ages from the cleaned tree
  dist_from_root <- node.depth.edgelength(mega_tree)
  root_h <- max(dist_from_root)
  mega_ages <- root_h - dist_from_root

  ## Match parsed annotations to tree nodes by divtime value (closest match)
  ## There are n_node-1 annotations (root excluded) or n_node.
  ## Build lookup: for each annotation, find the internal node with closest age.
  matched_ci_lo <- rep(NA_real_, n_node)
  matched_ci_hi <- rep(NA_real_, n_node)
  int_ages <- mega_ages[int_nodes]
  used <- logical(length(dt_vals))

  for (j in seq_along(int_nodes)) {
    best_k <- NA_integer_
    best_d <- Inf
    for (k in seq_along(dt_vals)) {
      if (used[k]) next
      d <- abs(int_ages[j] - dt_vals[k])
      if (d < best_d) { best_d <- d; best_k <- k }
    }
    if (!is.na(best_k) && best_d < 1e-6) {
      matched_ci_lo[j] <- ci_lo[best_k]
      matched_ci_hi[j] <- ci_hi[best_k]
      used[best_k] <- TRUE
    }
  }

  ## For the root (no CI annotation typically), set CI = point estimate
  root_idx <- which(int_nodes == n_tip + 1L)
  if (is.na(matched_ci_lo[root_idx])) {
    matched_ci_lo[root_idx] <- int_ages[root_idx]
    matched_ci_hi[root_idx] <- int_ages[root_idx]
  }

  ## Clamp CIs for near-zero-age nodes (MEGA's collapsed backbone branches).
  ## These nodes have age ≈ 0 but Tao CIs spanning the full tree depth,
  ## producing absurd FigTree bars.  Set their CIs to the point estimate.
  near_zero <- int_ages < 1e-4
  matched_ci_lo[near_zero] <- int_ages[near_zero]
  matched_ci_hi[near_zero] <- int_ages[near_zero]

  ci_df <- data.frame(
    node  = int_nodes,
    age   = int_ages,
    ci_lo = matched_ci_lo,
    ci_hi = matched_ci_hi,
    se    = (matched_ci_hi - matched_ci_lo) / (2 * 1.96),
    stringsAsFactors = FALSE
  )

  base$ci <- ci_df
  base
}

write_treepl_cfg <- function(cfg_path, tree_path, out_tree_path, smooth, numsites,
                             node_bounds, thorough = TRUE, prime = TRUE,
                             opt = NULL, plsimaniter = NULL, pliter = NULL,
                             extra_lines = NULL) {
  lines <- c(
    paste0("treefile = ", tree_path),
    paste0("outfile = ", out_tree_path),
    paste0("numsites = ", as.integer(numsites)),
    paste0("smooth = ", fmt_num(smooth))
  )
  ## Optimization settings (Smith & O'Meara 2012 strongly recommend thorough + prime)
  if (isTRUE(thorough)) lines <- c(lines, "thorough")
  if (isTRUE(prime))    lines <- c(lines, "prime")
  if (!is.null(opt) && nzchar(opt))
    lines <- c(lines, paste0("opt = ", as.integer(opt)))
  if (!is.null(plsimaniter) && nzchar(plsimaniter))
    lines <- c(lines, paste0("plsimaniter = ", as.integer(plsimaniter)))
  if (!is.null(pliter) && nzchar(pliter))
    lines <- c(lines, paste0("pliter = ", as.integer(pliter)))
  if (length(extra_lines)) {
    extra_lines <- trimws(as.character(extra_lines))
    extra_lines <- extra_lines[nzchar(extra_lines)]
    lines <- c(lines, extra_lines)
  }
  ## Calibration constraints
  for (i in seq_len(nrow(node_bounds))) {
    tag <- paste0("N", node_bounds$node[i])
    lines <- c(
      lines,
      paste0("mrca = ", tag, " ", node_bounds$taxonA[i], " ", node_bounds$taxonB[i]),
      paste0("min = ", tag, " ", fmt_num(node_bounds$age_min[i]))
    )
    if (is.finite(node_bounds$age_max[i])) {
      lines <- c(lines, paste0("max = ", tag, " ", fmt_num(node_bounds$age_max[i])))
    }
  }
  writeLines(lines, cfg_path)
}

parse_treepl_prime_lines <- function(log_path) {
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

run_treepl_one <- function(treepl_bin, cfg_path, log_path, omp_threads = 1L) {
  if (file.exists(log_path)) unlink(log_path)
  suppressWarnings(system2(
    treepl_bin,
    args = cfg_path,
    stdout = log_path,
    stderr = log_path,
    env = paste0("OMP_NUM_THREADS=", as.integer(omp_threads))
  ))
}

methods <- tolower(parse_chr_grid(kv[["methods"]] %||% "chronos,treepl,reltime", "--methods"))
valid_methods <- c("chronos", "treepl", "reltime", "reltime_mega")
bad_methods <- setdiff(methods, valid_methods)
if (length(bad_methods)) {
  stop("Unknown methods: ", paste(bad_methods, collapse = ", "))
}

chronos_lambdas <- parse_num_grid(kv[["chronos-lambdas"]] %||% "0.01,0.1,1,10,100", "--chronos-lambdas")
chronos_models <- parse_chr_grid(kv[["chronos-models"]] %||% "clock,correlated,relaxed,discrete", "--chronos-models")
chronos_models <- tolower(chronos_models)
bad_models <- setdiff(chronos_models, c("clock", "correlated", "relaxed", "discrete"))
if (length(bad_models)) stop("Unsupported chronos model(s): ", paste(bad_models, collapse = ", "))
chronos_discrete_k <- parse_int_grid(kv[["chronos-discrete-k"]] %||% "2,3,5", "--chronos-discrete-k")
chronos_retries <- as.integer(kv[["chronos-retries"]] %||% "2")
if (!is.finite(chronos_retries) || chronos_retries < 1L) stop("--chronos-retries must be >= 1.")
chronos_attempt_timeout <- as.numeric(kv[["chronos-attempt-timeout"]] %||% "0")
if (!is.finite(chronos_attempt_timeout) || chronos_attempt_timeout < 0) {
  stop("--chronos-attempt-timeout must be >= 0.")
}
chronos_iter_max <- as.integer(kv[["chronos-iter-max"]] %||% "10000")
chronos_tol <- as.numeric(kv[["chronos-tol"]] %||% "1e-8")

treepl_smoothing <- parse_num_grid(kv[["treepl-smoothing"]] %||% "0.01,0.1,1,10,100", "--treepl-smoothing")
treepl_numsites <- as.integer(kv[["treepl-numsites"]] %||% "1000")
treepl_threads <- as.integer(kv[["treepl-threads"]] %||% "1")
if (!is.finite(treepl_numsites) || treepl_numsites < 1L) stop("--treepl-numsites must be >= 1.")
if (!is.finite(treepl_threads) || treepl_threads < 1L) stop("--treepl-threads must be >= 1.")
treepl_thorough <- !identical(tolower(kv[["treepl-thorough"]] %||% "TRUE"), "false")
treepl_prime <- !identical(tolower(kv[["treepl-prime"]] %||% "TRUE"), "false")
treepl_opt <- kv[["treepl-opt"]] %||% NULL
treepl_plsimaniter <- kv[["treepl-plsimaniter"]] %||% NULL
treepl_pliter <- kv[["treepl-pliter"]] %||% NULL

ci_sites <- as.integer(kv[["ci-sites"]] %||% kv[["reltime-sites"]] %||% "1000")
if (!is.finite(ci_sites) || ci_sites < 1L) stop("--ci-sites must be >= 1.")
chronos_ci_type <- kv[["chronos-ci-type"]] %||% "parametric"
chronos_ci_reps <- as.integer(kv[["chronos-ci-reps"]] %||% "100")
if (!is.finite(chronos_ci_reps) || chronos_ci_reps < 1L) stop("--chronos-ci-reps must be >= 1.")
reltime_bootstrap_reps <- as.integer(kv[["reltime-bootstrap-reps"]] %||% "100")
if (!is.finite(reltime_bootstrap_reps) || reltime_bootstrap_reps < 1L) {
  stop("--reltime-bootstrap-reps must be >= 1.")
}
treepl_bootstrap_reps <- as.integer(kv[["treepl-bootstrap-reps"]] %||% "0")
if (!is.finite(treepl_bootstrap_reps) || treepl_bootstrap_reps < 0L) {
  stop("--treepl-bootstrap-reps must be >= 0.")
}
treepl_bootstrap_jobs <- as.integer(kv[["treepl-bootstrap-jobs"]] %||% "1")
if (!is.finite(treepl_bootstrap_jobs) || treepl_bootstrap_jobs < 1L) {
  stop("--treepl-bootstrap-jobs must be >= 1.")
}
chronos_pml_output <- load_chronos_pml_output(kv[["chronos-pml-rds"]] %||% "")
if (!nzchar(chronos_ci_type)) stop("--chronos-ci-type must be non-empty.")
if (!(chronos_ci_type %in% c("nonparametric", "semiparametric", "parametric"))) {
  stop("--chronos-ci-type must be one of nonparametric, semiparametric, parametric.")
}
if (chronos_ci_type %in% c("nonparametric", "semiparametric") && is.null(chronos_pml_output)) {
  stop("--chronos-pml-rds is required for --chronos-ci-type=", chronos_ci_type, ".")
}

phy_file <- normalizePath(kv[["phylogram"]], winslash = "/", mustWork = TRUE)
phy_tree_name <- kv[["phylogram-tree-name"]] %||% ""
phy_tree_index <- as.integer(kv[["phylogram-tree-index"]] %||% "1")
if (!is.finite(phy_tree_index) || phy_tree_index < 1L) stop("--phylogram-tree-index must be >= 1.")
phy_loaded <- load_tree_from_file(phy_file, tree_name = phy_tree_name, tree_index = phy_tree_index)
phy <- ape::ladderize(phy_loaded$tree)
if (is.null(phy$edge.length) || any(!is.finite(phy$edge.length)) || any(phy$edge.length <= 0)) {
  stop("Target phylogram must have strictly positive branch lengths.")
}

outdir <- kv[["outdir"]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir, winslash = "/", mustWork = TRUE)
prefix <- sanitize_id(kv[["out-prefix"]] %||% phy_loaded$tree_id)

log_file <- file.path(outdir, paste0(prefix, "_run_dating_grid.log"))
msg <- function(...) {
  s <- paste0(...)
  cat(s, "\n")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s), file = log_file, append = TRUE)
}

dir.create(file.path(outdir, "chronos", "trees"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "chronos", "ci"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "treepl", "trees"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "treepl", "ci"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "treepl", "configs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "treepl", "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "reltime"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "reltime_mega"), recursive = TRUE, showWarnings = FALSE)

target_tree_copy <- file.path(outdir, paste0(prefix, "_phylogram_used.tre"))
ape::write.tree(phy, file = target_tree_copy)

calibration_tag <- kv[["calibration-tag"]] %||% NULL
pair_df <- NULL
ref_tree_copy <- NA_character_
if (has_csv) {
  cal_csv <- normalizePath(kv[["calibrations-csv"]], winslash = "/", mustWork = TRUE)
  pair_df <- normalize_calibration_df(read.csv(cal_csv, stringsAsFactors = FALSE), calibration_tag = calibration_tag)
} else {
  ref_file <- normalizePath(kv[["reference-time-tree"]], winslash = "/", mustWork = TRUE)
  ref_tree_name <- kv[["reference-tree-name"]] %||% ""
  ref_tree_index <- as.integer(kv[["reference-tree-index"]] %||% "1")
  if (!is.finite(ref_tree_index) || ref_tree_index < 1L) stop("--reference-tree-index must be >= 1.")
  ref_loaded <- load_tree_from_file(ref_file, tree_name = ref_tree_name, tree_index = ref_tree_index)
  ref_tree <- ape::ladderize(ref_loaded$tree)
  ref_tree_copy <- file.path(outdir, paste0(prefix, "_reference_time_tree_used.tre"))
  ape::write.tree(ref_tree, file = ref_tree_copy)
  pair_df <- build_calib_pairs_congruif(phy, ref_tree)
}

if ("root-age" %in% names(kv)) {
  root_age <- as.numeric(kv[["root-age"]])
  if (!is.finite(root_age) || root_age <= 0) stop("--root-age must be a positive number.")
  pair_df <- append_root_calibration(pair_df, phy, root_age)
} else {
  root_age <- NA_real_
}

pair_df$pair_id <- seq_len(nrow(pair_df))
pair_map <- map_pairwise_calibrations_to_nodes(phy, pair_df)
node_bounds <- pair_map$merged

pair_csv <- file.path(outdir, paste0(prefix, "_calibration_pairs_used.csv"))
mapped_csv <- file.path(outdir, paste0(prefix, "_calibration_pairs_mapped.csv"))
bounds_csv <- file.path(outdir, paste0(prefix, "_calibration_bounds_used.csv"))
drop_csv <- file.path(outdir, paste0(prefix, "_calibration_bounds_dropped.csv"))
write.csv(pair_df, pair_csv, row.names = FALSE)
write.csv(pair_map$mapped, mapped_csv, row.names = FALSE)
write.csv(node_bounds, bounds_csv, row.names = FALSE)
if (nrow(pair_map$dropped_inconsistent)) {
  write.csv(pair_map$dropped_inconsistent, drop_csv, row.names = FALSE)
}

msg("Target tree: ", phy_file, " | tree_id=", phy_loaded$tree_id)
msg("Output prefix: ", prefix)
msg("Calibration rows input: ", nrow(pair_df))
msg("Calibration rows mapped: ", nrow(pair_map$mapped))
msg("Merged node bounds retained: ", nrow(node_bounds))
msg("Dropped calibrations (missing taxa): ", pair_map$dropped_missing_taxa)
msg("Dropped calibrations (bad MRCA): ", pair_map$dropped_bad_mrca)
if (nrow(pair_map$dropped_inconsistent)) {
  msg("Dropped calibration nodes (inconsistent after merge): ", nrow(pair_map$dropped_inconsistent))
}

candidate_rows <- list()
all_rows <- list()
uncertainty_rows <- list()
ci_inputs <- build_ci_inputs(phy, node_bounds)
chronos_runs_csv <- ""
treepl_runs_csv <- ""
reltime_run_csv <- ""

if ("chronos" %in% methods) {
  msg("Running chronos grid...")
  chronos_calib <- build_chronos_calib_from_node_bounds(phy, node_bounds)
  chronos_rows <- list()
  chronos_trees_cache <- list()   ## cache in-memory chronos trees for post-selection CI
  row_i <- 0L

  for (mdl in chronos_models) {
    k_grid <- if (identical(mdl, "discrete")) chronos_discrete_k else NA_integer_
    if (all(!is.finite(k_grid))) k_grid <- NA_integer_
    for (lam in chronos_lambdas) {
      for (kcat in k_grid) {
        candidate <- paste0("chronos_", mdl, "_lambda", fmt_token(lam))
        if (identical(mdl, "discrete") && is.finite(kcat)) {
          candidate <- paste0(candidate, "_k", as.integer(kcat))
        }
        out_tree <- file.path(outdir, "chronos", "trees", paste0(prefix, "_", candidate, ".tre"))
        ## Checkpoint: skip if tree already exists from a previous interrupted run
        if (file.exists(out_tree)) {
          cached_tree <- tryCatch(ape::read.tree(out_tree), error = function(e) NULL)
          if (!is.null(cached_tree) && inherits(cached_tree, "phylo")) {
            run <- list(status = "OK", tree = cached_tree, error = NA_character_)
            chronos_trees_cache[[candidate]] <- cached_tree
            msg("  checkpoint: reusing ", candidate)
          } else {
            run <- list(status = "FAIL", tree = NULL, error = "cached tree unreadable")
          }
        } else {
        run <- run_chronos_one(
          phy = phy,
          calib = chronos_calib,
          model = mdl,
          lambda = lam,
          nb_rate_cat = kcat,
          retries = chronos_retries,
          iter_max = chronos_iter_max,
          tol = chronos_tol,
          attempt_timeout = chronos_attempt_timeout
        )
        }
        ## If discrete model fails, try progressively dropping calibrations
        ## (chronos discrete init can fail with many fixed-point calibrations)
        discrete_dropped_cals <- character(0)
        if (!identical(run$status, "OK") && identical(mdl, "discrete")) {
          cal_subset <- chronos_calib
          max_drops <- min(nrow(chronos_calib) - 2L, ceiling(nrow(chronos_calib) * 0.6))
          for (.drop_round in seq_len(max_drops)) {
            found <- FALSE
            ## Try dropping each remaining cal; keep first that works
            for (.j in seq_len(nrow(cal_subset))) {
              test_cal <- cal_subset[-.j, , drop = FALSE]
              test_run <- run_chronos_one(phy = phy, calib = test_cal, model = mdl,
                lambda = lam, nb_rate_cat = kcat, retries = 3L,
                iter_max = chronos_iter_max, tol = chronos_tol,
                attempt_timeout = chronos_attempt_timeout)
              if (identical(test_run$status, "OK")) {
                discrete_dropped_cals <- c(discrete_dropped_cals,
                  sprintf("node=%d age=%.2f", cal_subset$node[.j], cal_subset$age.min[.j]))
                cal_subset <- test_cal
                run <- test_run
                found <- TRUE
                break
              }
            }
            if (found) break
            ## No single drop worked — remove shallowest cal and retry
            shallow_idx <- which.min(cal_subset$age.min)
            discrete_dropped_cals <- c(discrete_dropped_cals,
              sprintf("node=%d age=%.2f (shallowest)", cal_subset$node[shallow_idx], cal_subset$age.min[shallow_idx]))
            cal_subset <- cal_subset[-shallow_idx, , drop = FALSE]
          }
          if (identical(run$status, "OK") && length(discrete_dropped_cals)) {
            cat(sprintf("  discrete model: succeeded after dropping %d calibration(s):\n",
              length(discrete_dropped_cals)))
            for (.dc in discrete_dropped_cals) cat("    -", .dc, "\n")
            cat(sprintf("  Using %d / %d calibrations.\n", nrow(cal_subset), nrow(chronos_calib)))
          }
        }
        if (identical(run$status, "OK")) {
          ape::write.tree(run$tree, file = out_tree)
          chronos_trees_cache[[candidate]] <- run$tree  ## keep in-memory for post-selection CI
          out_ci <- ""
          candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
            candidate = candidate,
            tree_file = out_tree,
            stringsAsFactors = FALSE
          )
        } else {
          out_tree <- ""
          out_ci <- ""
        }
        row_i <- row_i + 1L
        chronos_rows[[row_i]] <- data.frame(
          candidate = candidate,
          method = "chronos",
          model = mdl,
          lambda = lam,
          nb_rate_cat = if (is.finite(kcat)) as.integer(kcat) else NA_integer_,
          status = run$status,
          phiic = run$phiic,
          ploglik = run$ploglik,
          fit_score = run$fit_score,
          attempts_used = run$attempt,
          tree_file = out_tree,
          ci_file = out_ci,
          n_sites = if (nzchar(out_ci)) ci_sites else NA_integer_,
          error = run$error %||% NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  chronos_df <- do.call(rbind, chronos_rows)
  chronos_runs_csv <- file.path(outdir, "chronos", paste0(prefix, "_chronos_runs.csv"))
  write.csv(chronos_df, chronos_runs_csv, row.names = FALSE)
  all_rows[[length(all_rows) + 1L]] <- chronos_df
  msg("chronos finished: ", sum(chronos_df$status == "OK"), "/", nrow(chronos_df), " successful.")
}

if ("treepl" %in% methods) {
  treepl_bin <- treepl_resolve_bin(kv[["treepl-bin"]] %||% NULL, repo_dir)
  if (!nzchar(treepl_bin)) {
    stop("treePL requested, but no executable was found. Pass --treepl-bin or set TREEPL_BIN.")
  }
  msg("Running treePL grid with binary: ", treepl_bin)
  treepl_rows <- list()
  treepl_input_tree <- file.path(outdir, "treepl", paste0(prefix, "_phylogram_for_treepl.tre"))
  ape::write.tree(phy, file = treepl_input_tree)
  for (i in seq_along(treepl_smoothing)) {
    smooth <- treepl_smoothing[i]
    smooth_tag <- fmt_token(smooth)
    candidate <- paste0("treePL_smooth", smooth_tag)
    cfg_path <- file.path(outdir, "treepl", "configs", paste0(prefix, "_", candidate, ".cfg"))
    out_tree <- file.path(outdir, "treepl", "trees", paste0(prefix, "_", candidate, ".tre"))
    log_path <- file.path(outdir, "treepl", "logs", paste0(prefix, "_", candidate, ".log"))
    prime_cfg_path <- file.path(outdir, "treepl", "configs", paste0(prefix, "_", candidate, ".prime.cfg"))
    prime_log_path <- file.path(outdir, "treepl", "logs", paste0(prefix, "_", candidate, ".prime.log"))
    ## Checkpoint: skip if tree already exists from a previous interrupted run
    if (file.exists(out_tree)) {
      cached_tree <- tryCatch(ape::read.tree(out_tree), error = function(e) NULL)
      if (!is.null(cached_tree) && inherits(cached_tree, "phylo")) {
        msg("  checkpoint: reusing ", candidate)
        treepl_rows[[length(treepl_rows) + 1L]] <- data.frame(
          candidate = candidate, method = "treePL", model = "treePL",
          lambda = NA, nb_rate_cat = NA, smooth = smooth,
          treepl_objective = NA_real_, status = "OK", exit_code = 0L,
          tree_file = normalizePath(out_tree), ci_file = "", n_sites = NA,
          cfg_file = cfg_path, log_file = log_path, error = NA_character_,
          stringsAsFactors = FALSE)
        candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
          candidate = candidate, tree_file = normalizePath(out_tree),
          method = "treePL", model = "treePL", lambda = NA, nb_rate_cat = NA,
          smooth = smooth, stringsAsFactors = FALSE)
        next
      }
    }
    if (file.exists(out_tree)) unlink(out_tree)
    if (file.exists(paste0(out_tree, ".r8s"))) unlink(paste0(out_tree, ".r8s"))

    prime_lines <- character(0)
    if (isTRUE(treepl_prime)) {
      treepl_write_cfg(
        prime_cfg_path, treepl_input_tree, out_tree, smooth, treepl_numsites, node_bounds,
        thorough = treepl_thorough, prime = TRUE,
        opt = treepl_opt, plsimaniter = treepl_plsimaniter, pliter = treepl_pliter
      )
      run_treepl_one(treepl_bin, prime_cfg_path, prime_log_path, omp_threads = treepl_threads)
      prime_lines <- treepl_parse_prime_lines(prime_log_path)
      if (file.exists(out_tree)) unlink(out_tree)
      if (file.exists(paste0(out_tree, ".r8s"))) unlink(paste0(out_tree, ".r8s"))
    }

    treepl_write_cfg(
      cfg_path, treepl_input_tree, out_tree, smooth, treepl_numsites, node_bounds,
      thorough = treepl_thorough, prime = FALSE,
      opt = treepl_opt, plsimaniter = treepl_plsimaniter, pliter = treepl_pliter,
      extra_lines = prime_lines
    )
    exit_code <- run_treepl_one(treepl_bin, cfg_path, log_path, omp_threads = treepl_threads)

    status <- "FAILED"
    err <- tail_text(log_path)
    out_ci <- ""
    ## treePL can exit non-zero (e.g. convergence warnings) yet still write
    ## a valid output tree.  Check the tree file regardless of exit code.
    if (file.exists(out_tree)) {
      tr <- try(ape::read.tree(out_tree), silent = TRUE)
      if (!inherits(tr, "try-error") && inherits(tr, "phylo") &&
          !is.null(tr$edge.length) && all(is.finite(tr$edge.length)) && all(tr$edge.length >= 0)) {
        status <- if (identical(exit_code, 0L)) "OK" else "OK_NONZERO_EXIT"
        err <- NA_character_
        candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
          candidate = candidate,
          tree_file = out_tree,
          stringsAsFactors = FALSE
        )
      }
    }

    treepl_rows[[i]] <- data.frame(
      candidate = candidate,
      method = "treePL",
      model = "treePL",
      lambda = NA_real_,
      nb_rate_cat = NA_integer_,
      smooth = smooth,
      treepl_objective = if (grepl("^OK", status)) parse_treepl_final_objective(log_path) else NA_real_,
      status = status,
      exit_code = exit_code,
      tree_file = if (grepl("^OK", status)) out_tree else "",
      ci_file = out_ci,
      n_sites = if (nzchar(out_ci)) ci_sites else NA_integer_,
      cfg_file = cfg_path,
      log_file = log_path,
      error = err %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }

  treepl_df <- do.call(rbind, treepl_rows)
  treepl_runs_csv <- file.path(outdir, "treepl", paste0(prefix, "_treepl_runs.csv"))
  write.csv(treepl_df, treepl_runs_csv, row.names = FALSE)
  all_rows[[length(all_rows) + 1L]] <- treepl_df
  msg("treePL finished: ", sum(grepl("^OK", treepl_df$status)), "/", nrow(treepl_df), " successful.")
}

if ("reltime" %in% methods) {
  msg("Running RelTime...")
  rel_candidate <- "RelTime"
  rel_tree_file <- file.path(outdir, "reltime", paste0(prefix, "_RelTime.tre"))
  rel_boot_ci_file <- file.path(outdir, "reltime", paste0(prefix, "_RelTime_bootstrap_ci.csv"))
  rel_tao_ci_file <- file.path(outdir, "reltime", paste0(prefix, "_RelTime_ci.csv"))
  rel_bounds_file <- file.path(outdir, "reltime", paste0(prefix, "_RelTime_bounds_used.csv"))
  rel_summary_file <- file.path(outdir, "reltime", paste0(prefix, "_RelTime_run.csv"))
  reltime_run_csv <- rel_summary_file

  rel_run <- try(
    reltime_bootstrap_ci(
      phy = phy,
      calibration_df = pair_df[, c("taxonA", "taxonB", "age_min", "age_max"), drop = FALSE],
      root_age = if (is.finite(root_age)) root_age else NULL,
      B = reltime_bootstrap_reps,
      n_sites = ci_sites,
      quiet = TRUE
    ),
    silent = TRUE
  )

  if (!inherits(rel_run, "try-error")) {
    ape::write.tree(rel_run$tree, file = rel_tree_file)
    write.csv(rel_run$ci, rel_boot_ci_file, row.names = FALSE)
    write.csv(rel_run$bounds, rel_bounds_file, row.names = FALSE)
    rel_tao_ci <- reltime_tao_ci_from_bounds_run(
      phy = phy,
      rel_run = rel_run,
      n_sites = ci_sites
    )
    write.csv(rel_tao_ci, rel_tao_ci_file, row.names = FALSE)
    ## Write FigTree-annotated NEXUS trees with HPD bars (use in-memory tree
    ## to avoid node-number reordering from write.tree → read.tree round-trip)
    rel_boot_annotated <- file.path(outdir, "reltime",
                                     paste0(prefix, "_RelTime_with_bootstrap_CI.tre"))
    annotate_tree_with_ci_from_objects(rel_run$tree, rel_run$ci, rel_boot_annotated)
    rel_tao_annotated <- file.path(outdir, "reltime",
                                    paste0(prefix, "_RelTime_with_tao_CI.tre"))
    annotate_tree_with_ci_from_objects(rel_run$tree, rel_tao_ci, rel_tao_annotated)
    rel_df <- data.frame(
      candidate = rel_candidate,
      method = "RelTime",
      model = "RelTime",
      lambda = NA_real_,
      nb_rate_cat = NA_integer_,
      smooth = NA_real_,
      status = "OK",
      tree_file = rel_tree_file,
      ci_file = rel_boot_ci_file,
      tao_ci_file = rel_tao_ci_file,
      bounds_file = rel_bounds_file,
      n_sites = ci_sites,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
    uncertainty_rows[[length(uncertainty_rows) + 1L]] <- ci_width_summary(
      rel_candidate,
      rel_run$ci,
      source = "RelTime-bootstrap"
    )
    candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
      candidate = rel_candidate,
      tree_file = rel_tree_file,
      stringsAsFactors = FALSE
    )
    msg("RelTime finished successfully.")
  } else {
    rel_df <- data.frame(
      candidate = rel_candidate,
      method = "RelTime",
      model = "RelTime",
      lambda = NA_real_,
      nb_rate_cat = NA_integer_,
      smooth = NA_real_,
      status = "FAILED",
      tree_file = "",
      ci_file = "",
      tao_ci_file = "",
      bounds_file = "",
      n_sites = ci_sites,
      error = conditionMessage(attr(rel_run, "condition")),
      stringsAsFactors = FALSE
    )
    msg("RelTime failed: ", rel_df$error[1L])
  }

  write.csv(rel_df, rel_summary_file, row.names = FALSE)
  all_rows[[length(all_rows) + 1L]] <- rel_df
}

## ---------------------------------------------------------------------------
## MEGA-CC RelTime (optional — runs only if megacc binary is available)
## ---------------------------------------------------------------------------
megacc_run_csv <- ""
if ("reltime_mega" %in% methods) {
  megacc_bin <- resolve_megacc_bin(kv[["megacc-bin"]] %||% NULL)
  if (!nzchar(megacc_bin)) {
    stop("reltime_mega requested, but no megacc executable was found. ",
         "Pass --megacc-bin or set MEGACC_BIN or place megacc on your PATH.")
  }
  msg("Running MEGA-CC RelTime with binary: ", megacc_bin)

  mega_work_dir <- file.path(outdir, "reltime_mega")
  mega_candidate <- "RelTime_MEGA"
  mega_summary_file <- file.path(mega_work_dir, paste0(prefix, "_RelTime_MEGA_run.csv"))
  mega_tao_ci_file  <- file.path(mega_work_dir, paste0(prefix, "_RelTime_MEGA_tao_ci.csv"))
  megacc_run_csv <- mega_summary_file

  ## Checkpoint: skip MEGA if output tree already exists
  mega_tree_file <- file.path(mega_work_dir, paste0(prefix, "_RelTime_MEGA.tre"))
  if (file.exists(mega_tree_file) && file.exists(mega_tao_ci_file)) {
    msg("  checkpoint: reusing MEGA RelTime from previous run")
    mega_result <- list(
      status = "OK",
      tree_file = mega_tree_file,
      ci = tryCatch(read.csv(mega_tao_ci_file), error = function(e) NULL)
    )
  } else {
  msg("Running MEGA-CC RelTime with Tao analytical CIs...")
  mega_result <- try(
    megacc_tao_ci(
      megacc_bin = megacc_bin,
      phy = phy,
      node_bounds = node_bounds,
      work_dir = mega_work_dir,
      prefix = prefix,
      n_sites = ci_sites
    ),
    silent = TRUE
  )
  }  ## end of MEGA checkpoint else

  if (!inherits(mega_result, "try-error") && mega_result$status == "OK") {
    has_ci <- !is.null(mega_result$ci) && is.data.frame(mega_result$ci) && nrow(mega_result$ci) > 0
    if (has_ci) {
      write.csv(mega_result$ci, mega_tao_ci_file, row.names = FALSE)
      msg("MEGA-CC RelTime Tao analytical CI written: ", mega_tao_ci_file)
      ## Write FigTree-annotated NEXUS tree with HPD bars
      mega_annotated_tree <- file.path(mega_work_dir,
                                        paste0(prefix, "_RelTime_MEGA_with_tao_CI.tre"))
      annotate_tree_with_ci_comments(mega_result$tree_file, mega_tao_ci_file,
                                      mega_annotated_tree)
      msg("MEGA-CC RelTime annotated tree: ", mega_annotated_tree)
    }
    mega_df <- data.frame(
      candidate = mega_candidate,
      method = "RelTime_MEGA",
      model = "RelTime_MEGA",
      lambda = NA_real_,
      nb_rate_cat = NA_integer_,
      smooth = NA_real_,
      status = "OK",
      tree_file = if (has_ci && file.exists(mega_annotated_tree)) mega_annotated_tree else mega_result$tree_file,
      ci_file = if (has_ci) mega_tao_ci_file else "",
      n_sites = if (has_ci) ci_sites else NA_integer_,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
    if (has_ci) {
      uncertainty_rows[[length(uncertainty_rows) + 1L]] <- ci_width_summary(
        mega_candidate,
        mega_result$ci,
        source = "RelTime_MEGA-Tao"
      )
    }
    candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
      candidate = mega_candidate,
      tree_file = if (has_ci && file.exists(mega_annotated_tree)) mega_annotated_tree else mega_result$tree_file,
      stringsAsFactors = FALSE
    )
    msg("MEGA-CC RelTime finished successfully.")
  } else {
    err_msg <- if (inherits(mega_result, "try-error")) {
      conditionMessage(attr(mega_result, "condition"))
    } else {
      mega_result$error
    }
    mega_df <- data.frame(
      candidate = mega_candidate,
      method = "RelTime_MEGA",
      model = "RelTime_MEGA",
      lambda = NA_real_,
      nb_rate_cat = NA_integer_,
      smooth = NA_real_,
      status = "FAILED",
      tree_file = "",
      ci_file = "",
      n_sites = NA_integer_,
      error = err_msg,
      stringsAsFactors = FALSE
    )
    msg("MEGA-CC RelTime failed: ", mega_df$error[1L])
  }

  write.csv(mega_df, mega_summary_file, row.names = FALSE)
  all_rows[[length(all_rows) + 1L]] <- mega_df
}

all_df <- bind_rows_fill(all_rows)
all_runs_csv <- file.path(outdir, paste0(prefix, "_all_runs_summary.csv"))
write.csv(all_df, all_runs_csv, row.names = FALSE)

uncertainty_df <- bind_rows_fill(uncertainty_rows)
full_candidates_df <- if (length(candidate_rows)) {
  do.call(rbind, candidate_rows)
} else {
  data.frame(candidate = character(0), tree_file = character(0), stringsAsFactors = FALSE)
}
if (nrow(full_candidates_df)) {
  full_candidates_df <- full_candidates_df[!duplicated(full_candidates_df$candidate), , drop = FALSE]
}
full_grid_dir <- file.path(outdir, "full_grid")
dir.create(full_grid_dir, recursive = TRUE, showWarnings = FALSE)
full_candidates_csv <- file.path(full_grid_dir, "candidates.csv")
write.csv(full_candidates_df, full_candidates_csv, row.names = FALSE, quote = TRUE)
full_uncertainty_csv <- file.path(full_grid_dir, "uncertainty_summary_long.csv")
if (nrow(uncertainty_df)) {
  write.csv(uncertainty_df, full_uncertainty_csv, row.names = FALSE)
}

candidates_df <- build_representative_candidates(all_df)
if (!nrow(candidates_df)) {
  candidates_df <- full_candidates_df
}

## ---- Post-selection CI computation (only for the 7 final trees) ----
msg("Computing CIs for selected candidates...")
ci_output_dir <- file.path(outdir, "selected_ci")
dir.create(ci_output_dir, recursive = TRUE, showWarnings = FALSE)

for (ci_i in seq_len(nrow(candidates_df))) {
  ci_cand   <- candidates_df$candidate[ci_i]
  ci_method <- candidates_df$method[ci_i]
  ci_tree_f <- candidates_df$tree_file[ci_i]
  if (!nzchar(ci_tree_f) || !file.exists(ci_tree_f)) next

  ## RelTime already has CIs — skip (they were computed inline, not in a grid)
  if (ci_method %in% c("RelTime", "RelTime_MEGA")) next

  ci_csv_out  <- file.path(ci_output_dir, paste0(prefix, "_", ci_cand, "_ci.csv"))
  ci_tree_out <- file.path(ci_output_dir, paste0(prefix, "_", ci_cand, "_with_CI.tre"))

  ## Checkpoint: skip if CI CSV already exists from a previous interrupted run
  if (file.exists(ci_csv_out)) {
    msg("  checkpoint: reusing CI for ", ci_cand)
    if (length(src_idx <- which(all_df$candidate == ci_cand)))
      all_df$ci_file[src_idx[1L]] <- ci_csv_out
    next
  }

  ci_tree     <- try(ape::read.tree(ci_tree_f), silent = TRUE)
  if (inherits(ci_tree, "try-error")) next

  if (identical(ci_method, "chronos")) {
    ## Look up model/lambda/nb_rate_cat from the source candidate row
    src_cand <- candidates_df$selected_from_candidate[ci_i] %||% ci_cand
    src_row  <- all_df[all_df$candidate == src_cand & all_df$status == "OK", , drop = FALSE]
    ci_model <- if (nrow(src_row)) src_row$model[1L] else "relaxed"
    ci_lambda <- if (nrow(src_row)) src_row$lambda[1L] else 1
    ci_nrc <- if (nrow(src_row) && "nb_rate_cat" %in% names(src_row)) src_row$nb_rate_cat[1L] else NA_integer_

    ## Use cached in-memory chronos tree but strip the call attribute so
    ## chronos_ci() rebuilds it cleanly — the original call can cause the
    ## bootstrap to return the starting ages unchanged (zero-width CIs).
    cached_tree <- chronos_trees_cache[[src_cand]]
    if (is.null(cached_tree)) {
      msg("  chronos CI skipped for ", ci_cand, " — no cached in-memory tree")
      next
    }
    attr(cached_tree, "call") <- NULL

    msg("  chronos CI for ", ci_cand, " (", chronos_ci_reps, " reps)...")
    ci_ctrl <- ape::chronos.control()
    ci_ctrl$iter.max <- as.integer(chronos_iter_max)
    ci_ctrl$eval.max <- as.integer(chronos_iter_max)
    ci_ctrl$tol      <- as.numeric(chronos_tol)
    if (identical(ci_model, "discrete") && is.finite(ci_nrc)) {
      ci_ctrl$nb.rate.cat <- as.integer(ci_nrc)
    }
    ci_res <- try(
      chronos_ci(
        phy, cached_tree,
        pml_output = chronos_pml_output,
        calibration = chronos_calib,
        B = chronos_ci_reps,
        type = chronos_ci_type,
        n_sites = ci_sites,
        model = ci_model,
        lambda = ci_lambda,
        nb_rate_cat = ci_nrc,
        control = ci_ctrl,
        quiet = TRUE
      ),
      silent = TRUE
    )
    if (!inherits(ci_res, "try-error")) {
      write.csv(ci_res, ci_csv_out, row.names = FALSE)
      annotate_tree_with_ci_from_objects(cached_tree, ci_res, ci_tree_out)
      uncertainty_rows[[length(uncertainty_rows) + 1L]] <- ci_width_summary(
        ci_cand, ci_res, source = "chronos-bootstrap"
      )
      src_idx <- which(all_df$candidate == src_cand & all_df$status == "OK")
      if (length(src_idx)) all_df$ci_file[src_idx[1L]] <- ci_csv_out
      msg("  -> ", ci_csv_out)
    } else {
      msg("  chronos CI failed for ", ci_cand, ": ", conditionMessage(attr(ci_res, "condition")))
    }

  } else if (identical(ci_method, "treePL")) {
    if (treepl_bootstrap_reps > 0L) {
      ci_smooth <- candidates_df$smooth[ci_i]
      if (!is.finite(ci_smooth)) ci_smooth <- 1
      msg("  treePL CI for ", ci_cand, " (smooth=", ci_smooth, ", ", treepl_bootstrap_reps, " reps)...")
      ci_res <- try(
        treepl_bootstrap_ci(
          treepl_bin = treepl_bin,
          phy = phy,
          calibration_df = pair_df[, c("taxonA", "taxonB", "age_min", "age_max"), drop = FALSE],
          smooth = ci_smooth,
          dated_tree = ci_tree,
          B = treepl_bootstrap_reps,
          n_sites = ci_sites,
          numsites = treepl_numsites,
          thorough = treepl_thorough,
          prime = treepl_prime,
          opt = treepl_opt,
          plsimaniter = treepl_plsimaniter,
          pliter = treepl_pliter,
          quiet = TRUE,
          omp_threads = treepl_threads,
          jobs = treepl_bootstrap_jobs,
          workdir = tempfile(pattern = paste0(prefix, "_", ci_cand, "_treepl_boot_")),
          prefix = paste0(prefix, "_", ci_cand)
        ),
        silent = TRUE
      )
      if (!inherits(ci_res, "try-error")) {
        write.csv(ci_res$ci, ci_csv_out, row.names = FALSE)
        annotate_tree_with_ci_from_objects(ci_tree, ci_res$ci, ci_tree_out)
        uncertainty_rows[[length(uncertainty_rows) + 1L]] <- ci_width_summary(
          ci_cand, ci_res$ci, source = "treePL-bootstrap"
        )
        src_cand <- candidates_df$selected_from_candidate[ci_i] %||% ci_cand
        src_idx <- which(all_df$candidate == src_cand & grepl("^OK", all_df$status))
        if (length(src_idx)) all_df$ci_file[src_idx[1L]] <- ci_csv_out
        msg("  -> ", ci_csv_out)
      } else {
        msg("  treePL CI failed for ", ci_cand, ": ", conditionMessage(attr(ci_res, "condition")))
      }
    }
  }
}

## Rebuild uncertainty_df with the new post-selection entries
uncertainty_df <- bind_rows_fill(uncertainty_rows)
msg("Post-selection CI computation done.")

key_dirs <- materialize_key_outputs(
  outdir,
  candidates_df,
  fit_files = c(
    all_runs_summary = all_runs_csv,
    chronos_runs = chronos_runs_csv,
    treepl_runs = treepl_runs_csv,
    reltime_run = reltime_run_csv,
    reltime_mega_run = megacc_run_csv
  )
)
candidates_df <- key_dirs$candidates_df
selected_uncertainty_df <- build_representative_uncertainty(candidates_df, uncertainty_df)
uncertainty_csv <- file.path(outdir, "uncertainty_summary_long.csv")
if (nrow(selected_uncertainty_df)) {
  write.csv(selected_uncertainty_df, uncertainty_csv, row.names = FALSE)
}
selected_ci_paths <- copy_representative_ci_files(outdir, candidates_df, all_df)
candidates_csv <- file.path(outdir, "candidates.csv")
write.csv(candidates_df, candidates_csv, row.names = FALSE, quote = TRUE)

meta_lines <- c(
  paste0("Target phylogram: ", phy_file),
  paste0("Target tree id: ", phy_loaded$tree_id),
  paste0("Output prefix: ", prefix),
  paste0("Methods requested: ", paste(methods, collapse = ",")),
  paste0("Calibration source: ", if (has_csv) "csv" else "congruify"),
  if (has_csv) paste0("Calibration CSV: ", normalizePath(kv[["calibrations-csv"]], winslash = "/", mustWork = TRUE)) else paste0("Reference time tree: ", normalizePath(kv[["reference-time-tree"]], winslash = "/", mustWork = TRUE)),
  if (is.finite(root_age)) paste0("Root age appended: ", fmt_num(root_age)) else "Root age appended: none",
  paste0("chronos bootstrap type: ", chronos_ci_type),
  paste0("chronos bootstrap replicates: ", chronos_ci_reps),
  paste0("chronos bootstrap pml source: ", if (!is.null(chronos_pml_output)) "external_pml_rds" else "parametric_proxy_from_phylogram"),
  paste0("treePL bootstrap replicates: ", treepl_bootstrap_reps),
  paste0("treePL bootstrap jobs: ", treepl_bootstrap_jobs),
  paste0("RelTime bootstrap replicates: ", reltime_bootstrap_reps),
  "RelTime Tao analytical CI: written as supplemental file only; not used in uncertainty_summary_long.csv",
  if ("reltime_mega" %in% methods) paste0("MEGA-CC binary: ", if (exists("megacc_bin") && nzchar(megacc_bin)) megacc_bin else "not found") else NULL,
  if ("reltime_mega" %in% methods) "MEGA-CC RelTime CI: Tao analytical (delta-method variance)" else NULL,
  paste0("Calibration pairs input: ", nrow(pair_df)),
  paste0("Calibration pairs mapped: ", nrow(pair_map$mapped)),
  paste0("Calibration bounds retained: ", nrow(node_bounds)),
  paste0("Calibration bounds dropped as inconsistent: ", nrow(pair_map$dropped_inconsistent)),
  paste0("Successful full-grid candidates written: ", nrow(full_candidates_df)),
  paste0("Representative candidates written: ", nrow(candidates_df)),
  paste0("Combined runs summary: ", all_runs_csv),
  paste0("Full-grid candidates CSV: ", full_candidates_csv),
  paste0("Full-grid uncertainty summary CSV: ", if (nrow(uncertainty_df)) full_uncertainty_csv else "none"),
  paste0("Candidates CSV: ", candidates_csv),
  paste0("MAIN_OUTPUT_TREES folder: ", key_dirs$main_tree_dir),
  paste0("CI CSVs folder: ", if (nzchar(selected_ci_paths$ci_dir)) selected_ci_paths$ci_dir else "none"),
  paste0("Model fits folder: ", key_dirs$model_fit_dir),
  paste0("Representative uncertainty summary CSV: ", if (nrow(selected_uncertainty_df)) uncertainty_csv else "none")
)
writeLines(meta_lines, file.path(outdir, paste0(prefix, "_run_metadata.txt")))

msg("Saved combined runs summary: ", all_runs_csv)
msg("Saved full-grid candidates CSV: ", full_candidates_csv)
msg("Saved candidates CSV: ", candidates_csv)
msg("Saved representative tree copies: ", key_dirs$main_tree_dir)
msg("Saved model-fit summaries: ", key_dirs$model_fit_dir)
msg("Done.")

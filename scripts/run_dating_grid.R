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
      "[--methods=chronos,treepl,reltime_mega]",
      "[--chronos-lambdas=0.01,0.1,1,10,100]",
      "[--treepl-smoothing=0.01,0.1,1,10,100]",
      "[--treepl-cv=loo]",
      "[--chronos-discrete-k=2,3,5]",
      "[--chronos-ci-models=clock,correlated,discrete,relaxed]",
      "[--chronos-attempt-timeout=90]",
      "[--chronos-age-start-tree=SEED_DATED_TREE.tre]",
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
      "[--seed=1]",
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
    "  --methods               Comma list from chronos,treepl,reltime_mega. Default: chronos,treepl,reltime_mega.\n",
    "                          RelTime runs only through MEGA-CC; 'reltime' is accepted as an alias for reltime_mega.\n",
    "  --megacc-bin            megacc executable for MEGA RelTime. Default: MEGACC_BIN env, PATH megacc.\n",
    "  --phylogram-tree-name   Select named tree from a named-Newick or multi-tree file.\n",
    "  --phylogram-tree-index  1-based fallback tree index. Default: 1.\n",
    "  --reference-tree-name   Optional named tree selector for the reference time tree.\n",
    "  --reference-tree-index  1-based fallback index for the reference tree. Default: 1.\n",
    "  --calibration-tag       If the calibration CSV has a candidate column, keep rows tagged with this value plus blank/all rows.\n",
    "  --additional-calibrations-csv  Optional CSV with additional calibrations (e.g., fossil min-ages) to merge\n",
    "                          with the primary calibration set (from --calibrations-csv or congruification).\n",
    "                          Format: taxonA,taxonB,age_min[,age_max]. Min-only constraints leave age_max blank.\n",
    "  --root-age              Optional exact root age appended to the shared calibration set for all methods.\n",
    "  --out-prefix            Prefix for output files. Default: target tree id.\n",
    "  --chronos-lambdas       Lambda grid. Default: 0.01,0.1,1,10,100.\n",
    "  --chronos-models        Clock models. Default: clock,correlated,relaxed,discrete.\n",
    "  --chronos-discrete-k    nb.rate.cat grid for the discrete model. Default: 2,3,5.\n",
    "  --chronos-ci-models     Which selected chronos models get bootstrap CIs. Default: clock,correlated,discrete,relaxed.\n",
    "  --chronos-retries       Retries per chronos setting. Default: 2.\n",
    "  --chronos-attempt-timeout Maximum elapsed seconds allowed per chronos attempt. Default: 0 (no timeout).\n",
    "  --chronos-age-start-tree Optional dated tree used only to seed chronos starting ages.\n",
    "  --treepl-smoothing      Smoothing grid for treePL cross-validation; must be a geometric series\n",
    "                          (constant ratio). Default: 0.01,0.1,1,10,100.\n",
    "  --treepl-cv             treePL cross-validation type: loo (leave-one-out, treePL default) or random\n",
    "                          (random subsample and replicate CV; faster on large trees). Default: loo.\n",
    "  --treepl-bin            treePL executable. Default: TREEPL_BIN env, PATH treePL, then ../treePL.\n",
    "  --treepl-numsites       numsites value written to treePL configs; set to the alignment length. Default: 1000.\n",
    "  --treepl-threads        OMP_NUM_THREADS for treePL. Default: 1.\n",
    "  --treepl-thorough       Use thorough optimization in treePL (strongly recommended). Default: TRUE.\n",
    "  --treepl-prime          Prime treePL optimizer settings before CV and dating. Default: TRUE.\n",
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
    "  --reltime-bootstrap-reps MEGA RelTime bootstrap replicates (megacc rerun per replicate). Default: 100.\n",
    "  --chronos-pml-rds       Optional RDS containing a phangorn pml.output for nonparametric/semiparametric chronos bootstrap.\n",
    "  --subset-tips           Trees with more tips than this tune chronos/treePL on a ~200-tip subset. Default: 500.\n",
    "  --seed                  Random seed (tuning subset, bootstraps). Default: 1.\n",
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

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file_arg)) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
repo_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "chronos_ci_helpers.R"))
source(file.path(script_dir, "treepl_helpers.R"))
source(file.path(script_dir, "reltime_mega_helpers.R"))

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
    for (nm in miss) df[[nm]] <- rep(NA, nrow(df))
    df[, cols, drop = FALSE]
  })
  do.call(rbind, padded)
}

read_csv_if_exists <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path) || !file.exists(path)) {
    return(data.frame())
  }
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
}

dedupe_rows_last <- function(df, key_cols = "candidate") {
  if (!is.data.frame(df) || !nrow(df)) return(df)
  key_cols <- intersect(key_cols, names(df))
  if (!length(key_cols)) return(df)
  key <- do.call(
    paste,
    c(unname(df[key_cols]), sep = "\r")
  )
  df[!duplicated(key, fromLast = TRUE), , drop = FALSE]
}

scalar_num_or_na <- function(x) {
  if (is.null(x) || length(x) != 1L) return(NA_real_)
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return(NA_real_)
  x
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
  ## NEXUS files may contain lines like "Dimensions NTax=..." that fool the
  ## named-Newick reader into fabricating bogus one-tip trees.  Only try the
  ## named-Newick path for non-NEXUS files.
  if (!is_nexus_file(path)) {
    x <- try(read_named_newick(path), silent = TRUE)
    if (!inherits(x, "try-error")) {
      if (nzchar(tree_name)) {
        if (!(tree_name %in% names(x))) stop("Named tree not found: ", tree_name)
        return(list(tree = x[[tree_name]], tree_id = tree_name))
      }
      if (tree_index < 1L || tree_index > length(x)) stop("Tree index out of range.")
      return(list(tree = x[[tree_index]], tree_id = names(x)[tree_index]))
    }
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

## ---------------------------------------------------------------------------
## build_calibration_preserving_subset
## Prunes a large phylogram to ~target_n tips while preserving:
##   1. All calibration-defining taxa (tip pairs for each MRCA calibration)
##   2. Phylogenetic diversity (tips spread across tree depth)
##   3. Heterotachous extremes (tips on longest and shortest terminal branches)
## ---------------------------------------------------------------------------
build_calibration_preserving_subset <- function(phy, node_bounds, target_n = 200L) {
  n_tip <- ape::Ntip(phy)
  if (n_tip <= target_n) return(list(subset_phy = phy, subset_bounds = node_bounds,
    n_original = n_tip, n_subset = n_tip))

  ## 1. Calibration-defining tips: two tips per calibrated node that preserve MRCA
  cal_tips <- character(0)
  for (i in seq_len(nrow(node_bounds))) {
    nd <- node_bounds$node[i]
    kids <- phy$edge[phy$edge[, 1] == nd, 2]
    for (k in kids) {
      if (k <= n_tip) {
        cal_tips <- union(cal_tips, phy$tip.label[k])
      } else {
        clade_tips <- sort(ape::extract.clade(phy, k)$tip.label)
        cal_tips <- union(cal_tips, clade_tips[1])
      }
    }
  }

  ## 2. Heterotachous extremes: tips with longest and shortest terminal branches
  tip_edges <- which(phy$edge[, 2] <= n_tip)
  tip_bl <- phy$edge.length[tip_edges]
  tip_ids <- phy$tip.label[phy$edge[tip_edges, 2]]
  names(tip_bl) <- tip_ids
  n_extreme <- max(10L, target_n %/% 10L)
  extreme_long  <- names(sort(tip_bl, decreasing = TRUE))[seq_len(min(n_extreme, length(tip_bl)))]
  extreme_short <- names(sort(tip_bl, decreasing = FALSE))[seq_len(min(n_extreme, length(tip_bl)))]
  must_keep <- unique(c(cal_tips, extreme_long, extreme_short))

  ## 3. Fill remaining slots with phylo-diversity-spread tips
  remaining <- setdiff(phy$tip.label, must_keep)
  n_add <- max(0L, target_n - length(must_keep))
  if (n_add > 0 && length(remaining) > 0) {
    ## Use node depth to spread tips across the tree
    depths <- ape::node.depth.edgelength(phy)[seq_len(n_tip)]
    names(depths) <- phy$tip.label
    rem_depths <- depths[remaining]
    ## Stratified: divide depth range into n_add bins, pick one per bin
    breaks <- seq(min(rem_depths), max(rem_depths) + 1e-12, length.out = n_add + 1L)
    selected <- character(0)
    for (b in seq_len(n_add)) {
      in_bin <- names(rem_depths[rem_depths >= breaks[b] & rem_depths < breaks[b + 1]])
      in_bin <- setdiff(in_bin, selected)
      if (length(in_bin)) selected <- c(selected, sample(in_bin, 1))
    }
    must_keep <- unique(c(must_keep, selected))
  }

  ## 4. Prune
  tips_to_keep <- intersect(must_keep, phy$tip.label)
  subset_phy <- ape::drop.tip(phy, setdiff(phy$tip.label, tips_to_keep))

  ## 5. Remap calibrations
  keep_cal <- logical(nrow(node_bounds))
  new_nodes <- integer(nrow(node_bounds))
  new_taxonA <- rep(NA_character_, nrow(node_bounds))
  new_taxonB <- rep(NA_character_, nrow(node_bounds))
  for (i in seq_len(nrow(node_bounds))) {
    clade_tips <- ape::extract.clade(phy, node_bounds$node[i])$tip.label
    sub_tips <- intersect(clade_tips, subset_phy$tip.label)
    if (length(sub_tips) >= 2) {
      new_nd <- ape::getMRCA(subset_phy, sub_tips)
      if (!is.null(new_nd)) {
        kids <- subset_phy$edge[subset_phy$edge[, 1] == new_nd, 2]
        if (length(kids) >= 2) {
          pick_tip <- function(node_id) {
            if (node_id <= ape::Ntip(subset_phy)) {
              subset_phy$tip.label[node_id]
            } else {
              sort(ape::extract.clade(subset_phy, node_id)$tip.label)[1]
            }
          }
          new_nodes[i] <- new_nd
          new_taxonA[i] <- pick_tip(kids[1])
          new_taxonB[i] <- pick_tip(kids[2])
          keep_cal[i] <- nzchar(new_taxonA[i]) && nzchar(new_taxonB[i]) &&
            !identical(new_taxonA[i], new_taxonB[i])
        }
      }
    }
  }
  subset_bounds <- node_bounds[keep_cal, , drop = FALSE]
  subset_bounds$node <- new_nodes[keep_cal]
  if ("taxonA" %in% names(subset_bounds)) subset_bounds$taxonA <- new_taxonA[keep_cal]
  if ("taxonB" %in% names(subset_bounds)) subset_bounds$taxonB <- new_taxonB[keep_cal]

  list(subset_phy = subset_phy, subset_bounds = subset_bounds,
       n_original = n_tip, n_subset = ape::Ntip(subset_phy))
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

tree_node_ages <- function(tr) {
  d <- ape::node.depth.edgelength(tr)
  ages <- max(d) - d
  ages[seq_len(ape::Ntip(tr))] <- 0
  ages
}

map_tree_node_ages_by_clade <- function(source_tree, target_tree) {
  if (!inherits(source_tree, "phylo") || !inherits(target_tree, "phylo")) return(NULL)
  if (ape::Ntip(source_tree) != ape::Ntip(target_tree)) return(NULL)
  if (!setequal(source_tree$tip.label, target_tree$tip.label)) return(NULL)

  src_int <- seq.int(ape::Ntip(source_tree) + 1L, ape::Ntip(source_tree) + source_tree$Nnode)
  tgt_int <- seq.int(ape::Ntip(target_tree) + 1L, ape::Ntip(target_tree) + target_tree$Nnode)
  src_age <- tree_node_ages(source_tree)

  src_key <- vapply(src_int, function(nd) {
    paste(descendant_tip_labels(source_tree, nd), collapse = "\r")
  }, character(1))
  tgt_key <- vapply(tgt_int, function(nd) {
    paste(descendant_tip_labels(target_tree, nd), collapse = "\r")
  }, character(1))

  idx <- match(tgt_key, src_key)
  if (anyNA(idx)) return(NULL)

  out <- rep(NA_real_, ape::Ntip(target_tree) + target_tree$Nnode)
  out[seq_len(ape::Ntip(target_tree))] <- 0
  out[tgt_int] <- src_age[src_int[idx]]
  out
}

build_seeded_chronos_calib_from_node_bounds <- function(phy, node_bounds, age_start_tree = NULL) {
  calib <- build_chronos_calib_from_node_bounds(phy, node_bounds)
  if (is.null(age_start_tree) || !inherits(age_start_tree, "phylo")) return(calib)

  seed_age <- tryCatch(
    map_tree_node_ages_by_clade(age_start_tree, phy),
    error = function(e) NULL
  )
  if (is.null(seed_age)) return(calib)

  total <- ape::Ntip(phy) + phy$Nnode
  lower <- rep(0, total)
  upper <- rep(Inf, total)
  lower[node_bounds$node] <- node_bounds$age_min
  upper[node_bounds$node] <- node_bounds$age_max
  seed_age <- project_seed_ages_to_bounds(phy, seed_age, lower, upper)
  if (is.null(seed_age)) return(calib)
  calib$age.start <- seed_age[node_bounds$node]
  calib
}

## Clamp seed ages into the node bounds, then raise parents (postorder) so
## every parent is older than its children. Returns NULL if the seed cannot
## be made consistent with the bounds; chronos then picks its own start ages.
project_seed_ages_to_bounds <- function(phy, node_age, lower, upper, eps = 1e-6) {
  n_tip <- ape::Ntip(phy)
  node_age <- pmin(pmax(node_age, lower), upper)
  node_age[seq_len(n_tip)] <- 0
  po <- ape::reorder.phylo(phy, "postorder")$edge
  for (k in seq_len(nrow(po))) {
    p <- po[k, 1L]
    ch <- po[k, 2L]
    if (node_age[p] < node_age[ch] + eps) node_age[p] <- min(node_age[ch] + eps, upper[p])
  }
  ordered <- all(node_age[phy$edge[, 1L]] > node_age[phy$edge[, 2L]])
  in_bounds <- all(node_age >= lower - 1e-9 & node_age <= upper + 1e-9)
  if (!ordered || !in_bounds) return(NULL)
  node_age
}

fit_score_chronos <- function(tr) {
  ph <- attr(tr, "PHIIC")
  phiic <- if (is.list(ph)) scalar_num_or_na(ph$PHIIC) else scalar_num_or_na(ph)
  if (is.finite(phiic)) return(phiic)
  pl <- scalar_num_or_na(attr(tr, "ploglik"))
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
    "scalar_num_or_na <- function(x) {",
    "  if (is.null(x) || length(x) != 1L) return(NA_real_)",
    "  x <- suppressWarnings(as.numeric(x))",
    "  if (!is.finite(x)) return(NA_real_)",
    "  x",
    "}",
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
    "  phiic <- if (is.list(ph)) scalar_num_or_na(ph$PHIIC) else scalar_num_or_na(ph)",
    "  ploglik <- scalar_num_or_na(attr(tr, 'ploglik'))",
    "  fit_score <- if (is.finite(phiic)) phiic else if (is.finite(ploglik)) -ploglik else Inf",
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
        phiic = {
          ph <- attr(tr, "PHIIC")
          if (is.list(ph)) scalar_num_or_na(ph$PHIIC) else scalar_num_or_na(ph)
        },
        ploglik = scalar_num_or_na(attr(tr, "ploglik")),
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
    ## Only the CV-selected smoothing is dated, so normally one row is left.
    if (!("cv_chisq" %in% names(treepl_df))) treepl_df$cv_chisq <- NA_real_
    ord <- order(
      ifelse(is.finite(treepl_df$cv_chisq), treepl_df$cv_chisq, Inf),
      treepl_df$candidate
    )
    best <- treepl_df[ord[1L], , drop = FALSE]
    add_i <- add_i + 1L
    reps[[add_i]] <- data.frame(
      candidate = best$candidate[1L],
      tree_file = best$tree_file[1L],
      method = "treePL",
      model = "treePL",
      lambda = NA_real_,
      nb_rate_cat = NA_integer_,
      smooth = best$smooth[1L],
      selected_from_candidate = best$candidate[1L],
      selection_score = best$cv_chisq[1L],
      selection_rule = "min_treepl_cv_chisq",
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

tree_tip_count <- function(tree_path) {
  if (!is.character(tree_path) || length(tree_path) != 1L || !nzchar(tree_path) || !file.exists(tree_path)) {
    return(NA_integer_)
  }
  tr <- try(read_any_tree(tree_path), silent = TRUE)
  if (inherits(tr, "try-error")) return(NA_integer_)
  if (inherits(tr, "multiPhylo")) tr <- tr[[1L]]
  if (!inherits(tr, "phylo")) return(NA_integer_)
  as.integer(ape::Ntip(tr))
}

tree_has_expected_tips <- function(tree_path, expected_tips) {
  is.finite(expected_tips) && identical(tree_tip_count(tree_path), as.integer(expected_tips))
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

  existing_main_trees <- list.files(main_tree_dir, pattern = "\\.tre$", full.names = TRUE)
  if (length(existing_main_trees)) unlink(existing_main_trees, force = TRUE)

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

## Fill MAIN_OUTPUT_TREES with the best available version of each selected
## candidate: the CI-annotated tree when its CI exists, else the plain tree.
## Later CI stages run in the same output folder and replace plain trees.
copy_representative_ci_files <- function(outdir, candidates_df, all_df) {
  if (!nrow(candidates_df)) {
    return(list(ci_dir = "", main_tree_dir = "", candidates_df = candidates_df))
  }
  main_dir <- file.path(outdir, "MAIN_OUTPUT_TREES")
  ci_dir   <- file.path(main_dir, "ci")
  dir.create(ci_dir, recursive = TRUE, showWarnings = FALSE)
  existing_ci_csv <- list.files(ci_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(existing_ci_csv)) unlink(existing_ci_csv, force = TRUE)
  copied_ci <- FALSE

  ci_expected <- function(method, model) {
    if (identical(method, "treePL")) return(treepl_bootstrap_reps > 0L)
    if (identical(method, "chronos")) {
      return(chronos_ci_reps > 0L && tolower(model) %in% chronos_ci_models)
    }
    FALSE
  }

  for (i in seq_len(nrow(candidates_df))) {
    cand_name <- candidates_df$candidate[i]
    cand_method <- candidates_df$method[i]

    if (identical(cand_method, "RelTime_MEGA")) {
      mega_dir <- file.path(outdir, "reltime_mega")
      boot_csv <- file.path(mega_dir, paste0(prefix, "_RelTime_MEGA_bootstrap_ci.csv"))
      boot_tre <- file.path(mega_dir, paste0(prefix, "_RelTime_MEGA_with_bootstrap_CI.tre"))
      tao_csv  <- file.path(mega_dir, paste0(prefix, "_RelTime_MEGA_tao_ci.csv"))
      tao_tre  <- file.path(mega_dir, paste0(prefix, "_RelTime_MEGA_with_tao_CI.tre"))
      plain    <- file.path(mega_dir, paste0(prefix, "_RelTime_MEGA.tre"))
      dst <- ""
      if (copy_if_present(tao_tre, file.path(main_dir, "RelTime_MEGA_with_tao_CI.tre"))) {
        copied_ci <- copy_if_present(tao_csv, file.path(ci_dir, "RelTime_MEGA_tao_ci.csv")) || copied_ci
        dst <- file.path(main_dir, "RelTime_MEGA_with_tao_CI.tre")
      }
      if (copy_if_present(boot_tre, file.path(main_dir, "RelTime_MEGA_with_bootstrap_CI.tre"))) {
        copied_ci <- copy_if_present(boot_csv, file.path(ci_dir, "RelTime_MEGA_ci.csv")) || copied_ci
        dst <- file.path(main_dir, "RelTime_MEGA_with_bootstrap_CI.tre")
      }
      if (!nzchar(dst) && copy_if_present(plain, file.path(main_dir, "RelTime_MEGA.tre"))) {
        dst <- file.path(main_dir, "RelTime_MEGA.tre")
      }
      if (nzchar(dst)) candidates_df$tree_file[i] <- dst
      next
    }

    src_id <- candidates_df$selected_from_candidate[i] %||% cand_name
    hit <- all_df[all_df$candidate == src_id, , drop = FALSE]
    hit_ci_file <- if (nrow(hit) && "ci_file" %in% names(hit)) as.character(hit$ci_file[1L]) else ""
    if (is.na(hit_ci_file)) hit_ci_file <- ""
    ci_csv_src <- if (nzchar(hit_ci_file) && file.exists(hit_ci_file)) {
      hit_ci_file
    } else {
      file.path(outdir, "selected_ci", paste0(prefix, "_", cand_name, "_ci.csv"))
    }
    ci_tree_src <- file.path(outdir, "selected_ci", paste0(prefix, "_", cand_name, "_with_CI.tre"))

    if (file.exists(ci_csv_src) && file.exists(ci_tree_src)) {
      dst <- file.path(main_dir, paste0(cand_name, "_with_CI.tre"))
      copy_if_present(ci_tree_src, dst)
      copy_if_present(ci_csv_src, file.path(ci_dir, paste0(cand_name, "_ci.csv")))
      copied_ci <- TRUE
    } else {
      if (ci_expected(cand_method, candidates_df$model[i])) {
        msg("  no CI available for ", cand_name, "; MAIN_OUTPUT_TREES keeps the plain tree")
      }
      dst <- file.path(main_dir, paste0(cand_name, ".tre"))
      copy_if_present(candidates_df$tree_file[i], dst)
    }
    candidates_df$tree_file[i] <- dst
  }

  list(
    ci_dir = if (copied_ci) ci_dir else "",
    main_tree_dir = main_dir,
    candidates_df = candidates_df
  )
}

methods <- tolower(parse_chr_grid(kv[["methods"]] %||% "chronos,treepl,reltime_mega", "--methods"))
methods <- unique(sub("^reltime$", "reltime_mega", methods))
valid_methods <- c("chronos", "treepl", "reltime_mega")
bad_methods <- setdiff(methods, valid_methods)
if (length(bad_methods)) {
  stop("Unknown methods: ", paste(bad_methods, collapse = ", "))
}

chronos_lambdas <- parse_num_grid(kv[["chronos-lambdas"]] %||% "0.01,0.1,1,10,100", "--chronos-lambdas")
chronos_models <- parse_chr_grid(kv[["chronos-models"]] %||% "clock,correlated,relaxed,discrete", "--chronos-models")
chronos_models <- tolower(chronos_models)
bad_models <- setdiff(chronos_models, c("clock", "correlated", "relaxed", "discrete"))
if (length(bad_models)) stop("Unsupported chronos model(s): ", paste(bad_models, collapse = ", "))
chronos_ci_models <- parse_chr_grid(
  kv[["chronos-ci-models"]] %||% "clock,correlated,discrete,relaxed",
  "--chronos-ci-models"
)
chronos_ci_models <- tolower(chronos_ci_models)
bad_ci_models <- setdiff(chronos_ci_models, c("clock", "correlated", "relaxed", "discrete"))
if (length(bad_ci_models)) stop("Unsupported chronos CI model(s): ", paste(bad_ci_models, collapse = ", "))
chronos_discrete_k <- parse_int_grid(kv[["chronos-discrete-k"]] %||% "2,3,5", "--chronos-discrete-k")
chronos_retries <- as.integer(kv[["chronos-retries"]] %||% "2")
if (!is.finite(chronos_retries) || chronos_retries < 1L) stop("--chronos-retries must be >= 1.")
chronos_attempt_timeout <- as.numeric(kv[["chronos-attempt-timeout"]] %||% "0")
if (!is.finite(chronos_attempt_timeout) || chronos_attempt_timeout < 0) {
  stop("--chronos-attempt-timeout must be >= 0.")
}
chronos_iter_max <- as.integer(kv[["chronos-iter-max"]] %||% "10000")
chronos_tol <- as.numeric(kv[["chronos-tol"]] %||% "1e-8")
chronos_allow_cal_drop <- !identical(
  tolower(kv[["chronos-allow-cal-drop"]] %||% "TRUE"),
  "false"
)

treepl_smoothing <- parse_num_grid(kv[["treepl-smoothing"]] %||% "0.01,0.1,1,10,100", "--treepl-smoothing")
treepl_cv_type <- tolower(kv[["treepl-cv"]] %||% "loo")
if (!(treepl_cv_type %in% c("loo", "random"))) stop("--treepl-cv must be loo or random.")
if ("treepl" %in% methods) treepl_cv <- treepl_cv_settings(treepl_smoothing)
treepl_numsites <- as.integer(kv[["treepl-numsites"]] %||% "1000")
treepl_threads <- as.integer(kv[["treepl-threads"]] %||% "1")
if (!is.finite(treepl_numsites) || treepl_numsites < 1L) stop("--treepl-numsites must be >= 1.")
if (!is.finite(treepl_threads) || treepl_threads < 1L) stop("--treepl-threads must be >= 1.")
treepl_thorough <- !identical(tolower(kv[["treepl-thorough"]] %||% "TRUE"), "false")
treepl_prime <- !identical(tolower(kv[["treepl-prime"]] %||% "TRUE"), "false")
treepl_opt <- kv[["treepl-opt"]] %||% NULL
treepl_plsimaniter <- kv[["treepl-plsimaniter"]] %||% NULL
treepl_pliter <- kv[["treepl-pliter"]] %||% NULL

ci_sites <- as.integer(kv[["ci-sites"]] %||% "1000")
if (!is.finite(ci_sites) || ci_sites < 1L) stop("--ci-sites must be >= 1.")
chronos_ci_type <- kv[["chronos-ci-type"]] %||% "parametric"
chronos_ci_reps <- as.integer(kv[["chronos-ci-reps"]] %||% "100")
if (!is.finite(chronos_ci_reps) || chronos_ci_reps < 0L) stop("--chronos-ci-reps must be >= 0.")
chronos_age_start_tree_path <- kv[["chronos-age-start-tree"]] %||% ""
reltime_bootstrap_reps <- as.integer(kv[["reltime-bootstrap-reps"]] %||% "100")
if (!is.finite(reltime_bootstrap_reps) || reltime_bootstrap_reps < 0L) {
  stop("--reltime-bootstrap-reps must be >= 0.")
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
chronos_age_start_tree <- NULL
if (nzchar(chronos_age_start_tree_path)) {
  chronos_seed_loaded <- load_tree_from_file(
    normalizePath(chronos_age_start_tree_path, winslash = "/", mustWork = TRUE),
    tree_index = 1L
  )
  chronos_age_start_tree <- ape::ladderize(chronos_seed_loaded$tree)
}

outdir <- kv[["outdir"]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir, winslash = "/", mustWork = TRUE)
prefix <- sanitize_id(kv[["out-prefix"]] %||% phy_loaded$tree_id)

existing_all_runs_csv <- file.path(outdir, paste0(prefix, "_all_runs_summary.csv"))
existing_full_candidates_csv <- file.path(outdir, "full_grid", "candidates.csv")
existing_full_uncertainty_csv <- file.path(outdir, "full_grid", "uncertainty_summary_long.csv")
existing_selected_uncertainty_csv <- file.path(outdir, "uncertainty_summary_long.csv")
existing_chronos_runs_csv <- file.path(outdir, "chronos", paste0(prefix, "_chronos_runs.csv"))
existing_treepl_runs_csv <- file.path(outdir, "treepl", paste0(prefix, "_treepl_runs.csv"))

log_file <- file.path(outdir, paste0(prefix, "_run_dating_grid.log"))
msg <- function(...) {
  s <- paste0(...)
  cat(s, "\n")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s), file = log_file, append = TRUE)
}

if ("chronos" %in% methods) {
  dir.create(file.path(outdir, "chronos", "trees"), recursive = TRUE, showWarnings = FALSE)
}
if ("treepl" %in% methods) {
  for (sub_dir in c("trees", "configs", "logs")) {
    dir.create(file.path(outdir, "treepl", sub_dir), recursive = TRUE, showWarnings = FALSE)
  }
}
if ("reltime_mega" %in% methods) {
  dir.create(file.path(outdir, "reltime_mega"), recursive = TRUE, showWarnings = FALSE)
}

existing_all_df <- dedupe_rows_last(read_csv_if_exists(existing_all_runs_csv), c("candidate", "method"))
existing_full_candidates_df <- dedupe_rows_last(read_csv_if_exists(existing_full_candidates_csv), "candidate")
existing_uncertainty_df <- bind_rows_fill(list(
  read_csv_if_exists(existing_full_uncertainty_csv),
  read_csv_if_exists(existing_selected_uncertainty_csv)
))
existing_uncertainty_df <- dedupe_rows_last(existing_uncertainty_df, "candidate")
existing_chronos_df <- dedupe_rows_last(read_csv_if_exists(existing_chronos_runs_csv), "candidate")
existing_treepl_df <- dedupe_rows_last(read_csv_if_exists(existing_treepl_runs_csv), "candidate")

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

## Merge additional manual calibrations (e.g., fossil min-age constraints)
## with congruification-derived or CSV-supplied calibrations.
if ("additional-calibrations-csv" %in% names(kv)) {
  add_csv <- normalizePath(kv[["additional-calibrations-csv"]], winslash = "/", mustWork = TRUE)
  add_df <- normalize_calibration_df(read.csv(add_csv, stringsAsFactors = FALSE), calibration_tag = calibration_tag)
  if (!is.null(add_df) && nrow(add_df) > 0L) {
    ## Ensure consistent columns; fill missing age_max with NA (min-only constraint)
    for (col in setdiff(names(pair_df), names(add_df))) add_df[[col]] <- NA
    for (col in setdiff(names(add_df), names(pair_df))) pair_df[[col]] <- NA
    pair_df <- rbind(pair_df, add_df[, names(pair_df), drop = FALSE])
    message(sprintf("  Appended %d additional calibration(s) from %s", nrow(add_df), basename(add_csv)))
  }
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

candidate_rows <- if (nrow(existing_full_candidates_df)) list(existing_full_candidates_df) else list()
all_rows <- if (nrow(existing_all_df)) list(existing_all_df) else list()
uncertainty_rows <- if (nrow(existing_uncertainty_df)) list(existing_uncertainty_df) else list()
chronos_runs_csv <- ""
chronos_trees_cache <- list()   ## in-memory chronos trees for post-selection CI
treepl_bin <- ""
treepl_runs_csv <- ""
treepl_prime_lines_full <- character(0)
megacc_run_csv <- ""

## Fixed seed so every stage of a staged run builds the same tuning subset.
set.seed(as.integer(kv[["seed"]] %||% "1"))

## ---------------------------------------------------------------------------
## Large-tree subset tuning (auto for >500 tips, override with --subset-tips)
## When active: tune lambda/k/smooth grids on a ~200-tip subset, then rerun
## only the winning settings on the full tree.
## ---------------------------------------------------------------------------
subset_threshold <- suppressWarnings(as.integer(kv[["subset-tips"]] %||% "500"))
if (!is.finite(subset_threshold) || subset_threshold < 50L) subset_threshold <- 500L
use_subset_tuning <- Ntip(phy) > subset_threshold

if (use_subset_tuning) {
  subset_target_n <- min(200L, Ntip(phy) %/% 2L)
  msg("Large tree detected (", Ntip(phy), " tips > ", subset_threshold, "). ",
      "Tuning grids on ~", subset_target_n, "-tip subset, then rerunning winners on full tree.")
  sub <- build_calibration_preserving_subset(phy, node_bounds, target_n = subset_target_n)
  msg("  Subset: ", sub$n_subset, " tips (", nrow(sub$subset_bounds), " calibrations retained)")
  phy_tune <- sub$subset_phy
  node_bounds_tune <- sub$subset_bounds
} else {
  phy_tune <- phy
  node_bounds_tune <- node_bounds
}

if ("chronos" %in% methods) {
  msg("Running chronos grid", if (use_subset_tuning) " (on subset)" else "", "...")
  ## Use subset tree for grid tuning if large-tree mode is active
  chronos_calib <- build_seeded_chronos_calib_from_node_bounds(
    phy_tune,
    node_bounds_tune,
    age_start_tree = chronos_age_start_tree
  )
  chronos_phy_grid <- phy_tune  ## tree used for grid search
  chronos_rows <- list()
  row_i <- 0L
  chronos_runs_csv <- file.path(outdir, "chronos", paste0(prefix, "_chronos_runs.csv"))
  update_incremental_chronos_outputs <- function() {
    chronos_df_partial <- bind_rows_fill(chronos_rows)
    if (!nrow(chronos_df_partial)) return(invisible(NULL))

    write.csv(chronos_df_partial, chronos_runs_csv, row.names = FALSE)

    incremental_all_df <- dedupe_rows_last(
      bind_rows_fill(list(existing_all_df, chronos_df_partial)),
      c("candidate", "method")
    )
    incremental_candidates_df <- build_representative_candidates(incremental_all_df)
    if (!nrow(incremental_candidates_df)) return(invisible(NULL))
    ## Subset-tuning trees must never reach MAIN_OUTPUT_TREES.
    if (use_subset_tuning) return(invisible(NULL))

    main_output_dir <- file.path(outdir, "MAIN_OUTPUT_TREES")
    dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)
    chronos_reps <- incremental_candidates_df[
      incremental_candidates_df$method == "chronos",
      ,
      drop = FALSE
    ]
    if (nrow(chronos_reps)) {
      for (.rr in seq_len(nrow(chronos_reps))) {
        src <- as.character(chronos_reps$tree_file[.rr])
        if (length(src) == 1L && nzchar(src) && file.exists(src)) {
          candidate_name <- as.character(chronos_reps$candidate[.rr])
          dst <- file.path(main_output_dir, paste0(candidate_name, ".tre"))
          ok <- file.copy(src, dst, overwrite = TRUE)
          if (isTRUE(ok)) {
            incremental_candidates_df$tree_file[
              incremental_candidates_df$candidate == candidate_name
            ] <- dst
            msg("  promoted current chronos winner to MAIN_OUTPUT_TREES: ", basename(dst))
          }
        }
      }
    }

    write.csv(
      incremental_candidates_df,
      file.path(outdir, "candidates.csv"),
      row.names = FALSE,
      quote = TRUE
    )
    invisible(NULL)
  }

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
        run <- NULL
        if (file.exists(out_tree)) {
          cached_tree <- tryCatch(ape::read.tree(out_tree), error = function(e) NULL)
          if (!is.null(cached_tree) && inherits(cached_tree, "phylo")) {
            cached_row <- existing_chronos_df[
              existing_chronos_df$candidate == candidate,
              ,
              drop = FALSE
            ]
            cached_phiic <- if (nrow(cached_row) && "phiic" %in% names(cached_row)) {
              scalar_num_or_na(cached_row$phiic[1L])
            } else {
              scalar_num_or_na(tryCatch(attr(cached_tree, "PHIIC"), error = function(e) NA_real_))
            }
            cached_ploglik <- if (nrow(cached_row) && "ploglik" %in% names(cached_row)) {
              scalar_num_or_na(cached_row$ploglik[1L])
            } else {
              scalar_num_or_na(tryCatch(attr(cached_tree, "ploglik"), error = function(e) NA_real_))
            }
            cached_fit_score <- if (nrow(cached_row) && "fit_score" %in% names(cached_row)) {
              scalar_num_or_na(cached_row$fit_score[1L])
            } else {
              cached_phiic
            }
            cached_attempts <- if (nrow(cached_row) && "attempts_used" %in% names(cached_row)) {
              suppressWarnings(as.integer(cached_row$attempts_used[1L]))
            } else {
              0L
            }
            if (!is.finite(cached_attempts)) cached_attempts <- 0L
            if (is.finite(cached_fit_score)) {
              run <- list(status = "OK", tree = cached_tree, error = NA_character_,
                          phiic = cached_phiic, ploglik = cached_ploglik,
                          fit_score = cached_fit_score, attempt = cached_attempts)
              chronos_trees_cache[[candidate]] <- cached_tree
              msg("  checkpoint: reusing ", candidate)
            } else {
              msg("  checkpoint tree lacks saved chronos score; rerunning ", candidate)
            }
          } else {
            msg("  checkpoint tree unreadable; rerunning ", candidate)
          }
        }
        if (is.null(run)) {
          run <- run_chronos_one(
            phy = chronos_phy_grid,
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
        if (!identical(run$status, "OK") && identical(mdl, "discrete") && isTRUE(chronos_allow_cal_drop)) {
          cal_subset <- chronos_calib
          max_drops <- min(nrow(chronos_calib) - 2L, ceiling(nrow(chronos_calib) * 0.6))
          for (.drop_round in seq_len(max_drops)) {
            found <- FALSE
            ## Try dropping each remaining cal; keep first that works
            for (.j in seq_len(nrow(cal_subset))) {
              test_cal <- cal_subset[-.j, , drop = FALSE]
              test_run <- run_chronos_one(phy = chronos_phy_grid, calib = test_cal, model = mdl,
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
        update_incremental_chronos_outputs()
      }
    }
  }

  chronos_df <- bind_rows_fill(chronos_rows)
  write.csv(chronos_df, chronos_runs_csv, row.names = FALSE)
  all_rows[[length(all_rows) + 1L]] <- chronos_df
  msg("chronos finished: ", sum(chronos_df$status == "OK"), "/", nrow(chronos_df), " successful.")
}

## ---------------------------------------------------------------------------
## treePL: prime -> cross-validate smoothing -> one thorough dating run.
## Priming and CV use the tuning tree (a subset for large trees); the final
## run always uses the full tree, primed on the full tree.
## ---------------------------------------------------------------------------
if ("treepl" %in% methods) {
  treepl_bin <- treepl_resolve_bin(kv[["treepl-bin"]] %||% NULL, repo_dir)
  if (!nzchar(treepl_bin)) {
    stop("treePL requested, but no executable was found. Pass --treepl-bin or set TREEPL_BIN.")
  }
  msg("Running treePL with binary: ", treepl_bin)
  tp_dir <- file.path(outdir, "treepl")
  tp_cfg <- function(name) file.path(tp_dir, "configs", paste0(prefix, "_treePL_", name, ".cfg"))
  tp_log <- function(name) file.path(tp_dir, "logs", paste0(prefix, "_treePL_", name, ".log"))
  tp_prime <- function(tree_path, bounds, name) {
    treepl_run_prime(
      treepl_bin, tree_path, bounds,
      smooth = exp(mean(log(treepl_cv$grid))), numsites = treepl_numsites,
      cfg_path = tp_cfg(name), log_path = tp_log(name), thorough = treepl_thorough,
      opt = treepl_opt, plsimaniter = treepl_plsimaniter, pliter = treepl_pliter,
      omp_threads = treepl_threads
    )
  }
  full_tree_path <- file.path(tp_dir, paste0(prefix, "_phylogram_for_treepl.tre"))
  ape::write.tree(phy, file = full_tree_path)
  cv_tree_path <- full_tree_path
  if (use_subset_tuning) {
    cv_tree_path <- file.path(tp_dir, paste0(prefix, "_phylogram_for_treepl_cv_subset.tre"))
    ape::write.tree(phy_tune, file = cv_tree_path)
  }
  cv_scores_csv <- file.path(tp_dir, paste0(prefix, "_treePL_cv_scores.csv"))
  cv_outfile <- file.path(tp_dir, paste0(prefix, "_treePL_cv.out"))
  final_status <- "FAILED"
  final_error <- NA_character_
  final_exit <- NA_integer_
  best_smooth <- NA_real_
  final_candidate <- ""
  final_tree <- ""

  treepl_res <- try({
    ## Step 1-2: prime + cross-validation (checkpointed by the CV score table)
    cv_scores <- read_csv_if_exists(cv_scores_csv)
    cv_reusable <- nrow(cv_scores) &&
      all(c("smooth", "cv_chisq", "cv_type") %in% names(cv_scores)) &&
      isTRUE(all.equal(sort(cv_scores$smooth), sort(treepl_cv$grid))) &&
      identical(unique(as.character(cv_scores$cv_type)), treepl_cv_type) &&
      any(is.finite(cv_scores$cv_chisq))
    if (cv_reusable) {
      msg("  checkpoint: reusing treePL cross-validation scores")
    } else {
      cv_prime_lines <- character(0)
      if (isTRUE(treepl_prime)) {
        msg("  priming treePL optimizer", if (use_subset_tuning) " on the tuning subset" else "", "...")
        cv_prime_lines <- tp_prime(cv_tree_path, node_bounds_tune,
                                   if (use_subset_tuning) "prime_cv_subset" else "prime")
      }
      msg("  treePL ", treepl_cv_type, " cross-validation, smoothing ",
          paste(fmt_num(treepl_cv$grid), collapse = ","), "...")
      cv_scores <- treepl_cross_validate(
        treepl_bin, cv_tree_path, node_bounds_tune, treepl_cv, treepl_cv_type,
        numsites = treepl_numsites, prime_lines = cv_prime_lines,
        cfg_path = tp_cfg("cv"), log_path = tp_log("cv"), cvoutfile = cv_outfile,
        thorough = treepl_thorough, opt = treepl_opt,
        plsimaniter = treepl_plsimaniter, pliter = treepl_pliter,
        omp_threads = treepl_threads
      )
      cv_scores$cv_type <- treepl_cv_type
      write.csv(cv_scores, cv_scores_csv, row.names = FALSE)
    }
    for (j in seq_len(nrow(cv_scores))) {
      msg("    smooth=", fmt_num(cv_scores$smooth[j]), "  CV chisq=", fmt_num(cv_scores$cv_chisq[j]))
    }
    best_smooth <- treepl_select_smoothing(cv_scores)
    if (!is.finite(best_smooth)) {
      stop("treePL cross-validation produced no usable scores; see ", tp_log("cv"))
    }
    msg("  CV-selected smoothing: ", fmt_num(best_smooth))
    if (best_smooth %in% range(treepl_cv$grid)) {
      msg("  WARNING: the CV-selected smoothing is at the edge of --treepl-smoothing; ",
          "the CV optimum may lie outside the grid. Consider extending the grid.")
    }

    ## Step 3: final dating run on the full tree
    if (isTRUE(treepl_prime)) {
      treepl_prime_lines_full <- treepl_parse_prime_lines(tp_log("prime"))
      if (!length(treepl_prime_lines_full)) {
        msg("  priming treePL optimizer on the full tree...")
        treepl_prime_lines_full <- tp_prime(full_tree_path, node_bounds, "prime")
      }
    }
    final_candidate <- paste0("treePL_smooth", fmt_token(best_smooth))
    final_tree <- file.path(tp_dir, "trees", paste0(prefix, "_", final_candidate, ".tre"))
    if (tree_has_expected_tips(final_tree, Ntip(phy)) && !is.null(.treepl_read_valid_tree(final_tree))) {
      msg("  checkpoint: reusing ", final_candidate)
      final_exit <- 0L
    } else {
      unlink(c(final_tree, paste0(final_tree, ".r8s")), force = TRUE)
      treepl_write_cfg(
        tp_cfg(final_candidate), full_tree_path, final_tree, best_smooth, treepl_numsites, node_bounds,
        thorough = treepl_thorough, prime = FALSE,
        opt = treepl_opt, plsimaniter = treepl_plsimaniter, pliter = treepl_pliter,
        extra_lines = treepl_prime_lines_full
      )
      final_exit <- treepl_run_one(treepl_bin, tp_cfg(final_candidate), tp_log(final_candidate),
                                   omp_threads = treepl_threads)
    }
    ## treePL can exit non-zero (e.g. convergence warnings) yet still write
    ## a valid output tree, so the tree file decides success.
    if (is.null(.treepl_read_valid_tree(final_tree))) {
      stop("treePL final run wrote no valid tree; ", tail_text(tp_log(final_candidate)))
    }
    final_status <- if (identical(final_exit, 0L)) "OK" else "OK_NONZERO_EXIT"
    cv_scores
  }, silent = TRUE)

  if (inherits(treepl_res, "try-error")) {
    final_error <- conditionMessage(attr(treepl_res, "condition"))
    msg("treePL failed: ", final_error)
    cv_scores <- read_csv_if_exists(cv_scores_csv)
    if (!nrow(cv_scores)) cv_scores <- data.frame(smooth = treepl_cv$grid, cv_chisq = NA_real_)
  } else {
    cv_scores <- treepl_res
  }

  treepl_rows <- lapply(seq_len(nrow(cv_scores)), function(j) {
    s <- cv_scores$smooth[j]
    is_sel <- is.finite(best_smooth) && isTRUE(all.equal(s, best_smooth))
    cand_j <- paste0("treePL_smooth", fmt_token(s))
    data.frame(
      candidate = cand_j,
      method = "treePL",
      model = "treePL",
      lambda = NA_real_,
      nb_rate_cat = NA_integer_,
      smooth = s,
      cv_type = treepl_cv_type,
      cv_chisq = cv_scores$cv_chisq[j],
      cv_selected = is_sel,
      status = if (is_sel) final_status else if (is.na(final_error)) "CV_ONLY" else "FAILED",
      exit_code = if (is_sel) final_exit else NA_integer_,
      tree_file = if (is_sel && grepl("^OK", final_status)) final_tree else "",
      ci_file = "",
      n_sites = NA_integer_,
      cfg_file = if (is_sel) tp_cfg(cand_j) else tp_cfg("cv"),
      log_file = if (is_sel) tp_log(cand_j) else tp_log("cv"),
      error = if (is_sel || !is.na(final_error)) final_error else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  treepl_df <- bind_rows_fill(treepl_rows)
  if (grepl("^OK", final_status)) {
    candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
      candidate = final_candidate,
      tree_file = final_tree,
      stringsAsFactors = FALSE
    )
  }
  treepl_runs_csv <- file.path(tp_dir, paste0(prefix, "_treepl_runs.csv"))
  write.csv(treepl_df, treepl_runs_csv, row.names = FALSE)
  all_rows[[length(all_rows) + 1L]] <- treepl_df
  msg("treePL finished: ", if (grepl("^OK", final_status)) paste0(final_candidate, " dated") else "FAILED")
}

## ---------------------------------------------------------------------------
## RelTime through MEGA-CC: point estimate + Tao analytical CI, then a Poisson
## branch-length bootstrap that reruns megacc on each perturbed phylogram.
## Always runs on the full tree.
## ---------------------------------------------------------------------------
if ("reltime_mega" %in% methods) {
  megacc_bin <- resolve_megacc_bin(kv[["megacc-bin"]] %||% NULL)
  if (!nzchar(megacc_bin)) {
    stop("reltime_mega requested, but no megacc executable was found. ",
         "Pass --megacc-bin or set MEGACC_BIN or place megacc on your PATH.")
  }
  msg("Running MEGA-CC RelTime with binary: ", megacc_bin)

  mega_work_dir <- file.path(outdir, "reltime_mega")
  mega_candidate <- "RelTime_MEGA"
  mega_file <- function(suffix) file.path(mega_work_dir, paste0(prefix, "_RelTime_MEGA", suffix))
  mega_tree_file <- mega_file(".tre")
  mega_tao_ci_file <- mega_file("_tao_ci.csv")
  mega_boot_ci_file <- mega_file("_bootstrap_ci.csv")
  megacc_run_csv <- mega_file("_run.csv")
  mega_status <- "OK"
  mega_error <- NA_character_

  if (file.exists(mega_tree_file) && file.exists(mega_tao_ci_file)) {
    msg("  checkpoint: reusing MEGA RelTime point estimate")
  } else {
    mega_result <- try(
      megacc_tao_ci(megacc_bin, phy, node_bounds, mega_work_dir, prefix),
      silent = TRUE
    )
    if (inherits(mega_result, "try-error")) {
      mega_status <- "FAILED"
      mega_error <- conditionMessage(attr(mega_result, "condition"))
    } else if (!identical(mega_result$status, "OK")) {
      mega_status <- "FAILED"
      mega_error <- mega_result$error
    } else if (is.data.frame(mega_result$ci) && nrow(mega_result$ci)) {
      write.csv(mega_result$ci, mega_tao_ci_file, row.names = FALSE)
      annotate_tree_with_ci_from_objects(mega_result$tree, mega_result$ci, mega_file("_with_tao_CI.tre"))
      msg("  MEGA Tao analytical CI written: ", mega_tao_ci_file)
    }
  }

  if (identical(mega_status, "OK") && reltime_bootstrap_reps > 0L) {
    if (file.exists(mega_boot_ci_file) && file.exists(mega_file("_with_bootstrap_CI.tre"))) {
      msg("  checkpoint: reusing MEGA RelTime bootstrap CI")
    } else {
      msg("  MEGA RelTime bootstrap CI (", reltime_bootstrap_reps, " megacc reruns)...")
      mega_tree <- ape::read.tree(mega_tree_file)
      boot_ci <- try(
        megacc_bootstrap_ci(megacc_bin, phy, node_bounds, mega_tree,
                            B = reltime_bootstrap_reps, n_sites = ci_sites, prefix = prefix),
        silent = TRUE
      )
      if (inherits(boot_ci, "try-error")) {
        msg("  MEGA RelTime bootstrap CI failed: ", conditionMessage(attr(boot_ci, "condition")))
      } else {
        write.csv(boot_ci, mega_boot_ci_file, row.names = FALSE)
        annotate_tree_with_ci_from_objects(mega_tree, boot_ci, mega_file("_with_bootstrap_CI.tre"))
        msg("  -> ", mega_boot_ci_file)
      }
    }
  }

  if (identical(mega_status, "OK")) {
    ## The bootstrap CI shares the Poisson scheme used for chronos and treePL,
    ## so it feeds uncertainty_summary_long.csv; Tao is the fallback.
    if (file.exists(mega_boot_ci_file)) {
      uncertainty_rows[[length(uncertainty_rows) + 1L]] <- ci_width_summary(
        mega_candidate, read.csv(mega_boot_ci_file), source = "RelTime_MEGA-bootstrap"
      )
    } else if (file.exists(mega_tao_ci_file)) {
      uncertainty_rows[[length(uncertainty_rows) + 1L]] <- ci_width_summary(
        mega_candidate, read.csv(mega_tao_ci_file), source = "RelTime_MEGA-Tao"
      )
    }
    candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
      candidate = mega_candidate,
      tree_file = mega_tree_file,
      stringsAsFactors = FALSE
    )
    msg("MEGA-CC RelTime finished successfully.")
  } else {
    msg("MEGA-CC RelTime failed: ", mega_error)
  }

  mega_df <- data.frame(
    candidate = mega_candidate,
    method = "RelTime_MEGA",
    model = "RelTime_MEGA",
    lambda = NA_real_,
    nb_rate_cat = NA_integer_,
    smooth = NA_real_,
    status = mega_status,
    tree_file = if (identical(mega_status, "OK")) mega_tree_file else "",
    ci_file = if (file.exists(mega_boot_ci_file)) mega_boot_ci_file else "",
    tao_ci_file = if (file.exists(mega_tao_ci_file)) mega_tao_ci_file else "",
    n_sites = ci_sites,
    error = mega_error,
    stringsAsFactors = FALSE
  )
  write.csv(mega_df, megacc_run_csv, row.names = FALSE)
  all_rows[[length(all_rows) + 1L]] <- mega_df
}

all_df <- dedupe_rows_last(bind_rows_fill(all_rows), c("candidate", "method"))
all_runs_csv <- file.path(outdir, paste0(prefix, "_all_runs_summary.csv"))
write.csv(all_df, all_runs_csv, row.names = FALSE)

uncertainty_df <- dedupe_rows_last(bind_rows_fill(uncertainty_rows), "candidate")
full_candidates_df <- if (length(candidate_rows)) {
  dedupe_rows_last(bind_rows_fill(candidate_rows), "candidate")
} else {
  data.frame(candidate = character(0), tree_file = character(0), stringsAsFactors = FALSE)
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

## ---- Subset tuning: rerun chronos winners on the full tree ----
## treePL already dates the full tree after subset CV, and RelTime always
## runs on the full tree, so only chronos needs a rerun here.
if (use_subset_tuning) {
  msg("Checking full-tree versions of the selected chronos candidates (", Ntip(phy), " tips)...")
  full_chronos_calib <- build_seeded_chronos_calib_from_node_bounds(
    phy,
    node_bounds,
    age_start_tree = chronos_age_start_tree
  )
  expected_full_tips <- as.integer(Ntip(phy))
  for (ri in seq_len(nrow(candidates_df))) {
    cand <- candidates_df$candidate[ri]
    if (!identical(candidates_df$method[ri], "chronos")) next
    src <- candidates_df$selected_from_candidate[ri] %||% cand
    src_row <- all_df[all_df$candidate == src & all_df$status == "OK", , drop = FALSE]
    if (!nrow(src_row)) next
    mdl <- src_row$model[1]; lam <- src_row$lambda[1]; kcat <- src_row$nb_rate_cat[1]
    full_tree_path <- file.path(outdir, "chronos", "trees",
      paste0(prefix, "_", cand, "_fulltree.tre"))
    if (file.exists(full_tree_path) && !tree_has_expected_tips(full_tree_path, expected_full_tips)) {
      unlink(full_tree_path, force = TRUE)
    }
    if (file.exists(full_tree_path)) {
      candidates_df$tree_file[ri] <- normalizePath(full_tree_path)
      msg("  checkpoint: reusing full-tree ", cand)
      next
    }
    if (!("chronos" %in% methods)) next
    full_run <- run_chronos_one(phy = phy, calib = full_chronos_calib,
      model = mdl, lambda = lam, nb_rate_cat = if (is.finite(kcat)) kcat else NA_integer_,
      retries = chronos_retries, iter_max = chronos_iter_max, tol = chronos_tol,
      attempt_timeout = chronos_attempt_timeout)
    if (identical(full_run$status, "OK")) {
      ape::write.tree(full_run$tree, full_tree_path)
      if (tree_has_expected_tips(full_tree_path, expected_full_tips)) {
        candidates_df$tree_file[ri] <- normalizePath(full_tree_path)
        chronos_trees_cache[[cand]] <- full_run$tree
        msg("  ", cand, ": rerun OK on full tree")
      } else {
        unlink(full_tree_path, force = TRUE)
        msg("  ", cand, ": rerun produced a non-full tree — dropping from final outputs")
      }
    } else {
      msg("  ", cand, ": rerun FAILED on full tree — dropping from final outputs")
    }
  }

  final_tip_counts <- vapply(candidates_df$tree_file, tree_tip_count, integer(1))
  keep_full <- is.finite(final_tip_counts) & final_tip_counts == expected_full_tips
  if (any(!keep_full)) {
    dropped <- candidates_df$candidate[!keep_full]
    msg("Dropping non-full candidates from final outputs: ", paste(dropped, collapse = ", "))
    candidates_df <- candidates_df[keep_full, , drop = FALSE]
  }
}

## ---- Post-selection bootstrap CIs for the selected chronos and treePL trees ----
## RelTime_MEGA CIs are computed in the MEGA block above.
msg("Computing CIs for selected candidates...")
ci_output_dir <- file.path(outdir, "selected_ci")
dir.create(ci_output_dir, recursive = TRUE, showWarnings = FALSE)

for (ci_i in seq_len(nrow(candidates_df))) {
  ci_cand   <- candidates_df$candidate[ci_i]
  ci_method <- candidates_df$method[ci_i]
  ci_tree_f <- candidates_df$tree_file[ci_i]
  if (!nzchar(ci_tree_f) || !file.exists(ci_tree_f)) next
  if (!(ci_method %in% c("chronos", "treePL"))) next

  ci_csv_out  <- file.path(ci_output_dir, paste0(prefix, "_", ci_cand, "_ci.csv"))
  ci_tree_out <- file.path(ci_output_dir, paste0(prefix, "_", ci_cand, "_with_CI.tre"))

  ## Checkpoint: skip if CI CSV already exists from a previous interrupted run
  if (file.exists(ci_csv_out)) {
    msg("  checkpoint: reusing CI for ", ci_cand)
    if (length(src_idx <- which(all_df$candidate == ci_cand)))
      all_df$ci_file[src_idx[1L]] <- ci_csv_out
    next
  }

  run_chronos_ci <- identical(ci_method, "chronos") && "chronos" %in% methods && chronos_ci_reps > 0L
  run_treepl_ci <- identical(ci_method, "treePL") && "treepl" %in% methods && treepl_bootstrap_reps > 0L
  if (!run_chronos_ci && !run_treepl_ci) next

  ci_tree     <- try(ape::read.tree(ci_tree_f), silent = TRUE)
  if (inherits(ci_tree, "try-error")) next

  if (run_chronos_ci) {
    ## Look up model/lambda/nb_rate_cat from the source candidate row
    src_cand <- candidates_df$selected_from_candidate[ci_i] %||% ci_cand
    src_row  <- all_df[all_df$candidate == src_cand & all_df$status == "OK", , drop = FALSE]
    ci_model <- if (nrow(src_row)) src_row$model[1L] else "relaxed"
    ci_lambda <- if (nrow(src_row)) src_row$lambda[1L] else 1
    ci_nrc <- if (nrow(src_row) && "nb_rate_cat" %in% names(src_row)) src_row$nb_rate_cat[1L] else NA_integer_

    if (!(tolower(ci_model) %in% chronos_ci_models)) {
      msg("  chronos CI skipped for ", ci_cand, " (model ", ci_model, " not requested in --chronos-ci-models)")
      next
    }

    ## Prefer the file-backed full-tree chronogram for CI generation.
    ## In-memory cached chronos objects can carry state that makes downstream
    ## chronosCI bootstraps brittle on some large empirical trees. The
    ## file-backed full-tree chronogram is the safer source for stage-5 CI.
    ci_source_tree <- if (inherits(ci_tree, "phylo")) ci_tree else chronos_trees_cache[[src_cand]]
    if (is.null(ci_source_tree) || !inherits(ci_source_tree, "phylo")) {
      msg("  chronos CI skipped for ", ci_cand, " — no usable full-tree chronogram")
      next
    }
    attr(ci_source_tree, "call") <- NULL

    msg("  chronos CI for ", ci_cand, " (", chronos_ci_reps, " reps)...")
    ci_ctrl <- ape::chronos.control()
    ci_ctrl$iter.max <- as.integer(chronos_iter_max)
    ci_ctrl$eval.max <- as.integer(chronos_iter_max)
    ci_ctrl$tol      <- as.numeric(chronos_tol)
    if (identical(ci_model, "discrete") && is.finite(ci_nrc)) {
      ci_ctrl$nb.rate.cat <- as.integer(ci_nrc)
    }
    ci_calibration <- if (use_subset_tuning) full_chronos_calib else chronos_calib
    ci_res <- try(
      chronos_ci(
        phy, ci_source_tree,
        pml_output = chronos_pml_output,
        calibration = ci_calibration,
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
      annotate_tree_with_ci_from_objects(ci_source_tree, ci_res, ci_tree_out)
      uncertainty_rows[[length(uncertainty_rows) + 1L]] <- ci_width_summary(
        ci_cand, ci_res, source = "chronos-bootstrap"
      )
      src_idx <- which(all_df$candidate == src_cand & all_df$status == "OK")
      if (length(src_idx)) all_df$ci_file[src_idx[1L]] <- ci_csv_out
      msg("  -> ", ci_csv_out)
    } else {
      msg("  chronos CI failed for ", ci_cand, ": ", conditionMessage(attr(ci_res, "condition")))
    }

  } else if (run_treepl_ci) {
    ci_smooth <- candidates_df$smooth[ci_i]
    msg("  treePL CI for ", ci_cand, " (smooth=", fmt_num(ci_smooth), ", ", treepl_bootstrap_reps, " reps)...")
    ci_res <- try(
      treepl_bootstrap_ci(
        treepl_bin = treepl_bin,
        phy = phy,
        node_bounds = node_bounds,
        smooth = ci_smooth,
        dated_tree = ci_tree,
        prime_lines = treepl_prime_lines_full,
        B = treepl_bootstrap_reps,
        n_sites = ci_sites,
        numsites = treepl_numsites,
        thorough = treepl_thorough,
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

## Rebuild uncertainty_df with the new post-selection entries
uncertainty_df <- dedupe_rows_last(bind_rows_fill(uncertainty_rows), "candidate")
msg("Post-selection CI computation done.")

key_dirs <- materialize_key_outputs(
  outdir,
  candidates_df,
  fit_files = c(
    all_runs_summary = all_runs_csv,
    chronos_runs = chronos_runs_csv,
    treepl_runs = treepl_runs_csv,
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
candidates_df <- selected_ci_paths$candidates_df
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
  paste0("Subset tuning: ", if (use_subset_tuning) paste0("yes (", Ntip(phy_tune), " of ", Ntip(phy), " tips)") else "no"),
  paste0("chronos bootstrap type: ", chronos_ci_type),
  paste0("chronos bootstrap replicates: ", chronos_ci_reps),
  paste0("chronos bootstrap CI models: ", paste(chronos_ci_models, collapse = ",")),
  paste0("chronos bootstrap pml source: ", if (!is.null(chronos_pml_output)) "external_pml_rds" else "parametric_proxy_from_phylogram"),
  if ("treepl" %in% methods) paste0("treePL smoothing selection: ", treepl_cv_type, " cross-validation over ",
                                    paste(fmt_num(treepl_cv$grid), collapse = ","), " (lowest CV chi-square wins)") else NULL,
  paste0("treePL numsites: ", treepl_numsites),
  paste0("treePL bootstrap replicates: ", treepl_bootstrap_reps),
  paste0("treePL bootstrap jobs: ", treepl_bootstrap_jobs),
  if ("reltime_mega" %in% methods) paste0("MEGA-CC binary: ", megacc_bin) else NULL,
  paste0("RelTime (MEGA) bootstrap replicates: ", reltime_bootstrap_reps, " (megacc rerun per replicate)"),
  "RelTime (MEGA) Tao analytical CI: supplemental; uncertainty_summary_long.csv uses the bootstrap CI when present",
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

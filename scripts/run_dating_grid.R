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
      "[--chronos-discrete-k=5]",
      "[--treepl-bin=/path/to/treePL]",
      "[--root-age=123.4]",
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
    "  --methods               Comma list from chronos,treepl,reltime. Default: all three.\n",
    "  --phylogram-tree-name   Select named tree from a named-Newick or multi-tree file.\n",
    "  --phylogram-tree-index  1-based fallback tree index. Default: 1.\n",
    "  --reference-tree-name   Optional named tree selector for the reference time tree.\n",
    "  --reference-tree-index  1-based fallback index for the reference tree. Default: 1.\n",
    "  --calibration-tag       If the calibration CSV has a candidate column, keep rows tagged with this value plus blank/all rows.\n",
    "  --root-age              Optional exact root age appended to the shared calibration set for all methods.\n",
    "  --out-prefix            Prefix for output files. Default: target tree id.\n",
    "  --chronos-lambdas       Lambda grid. Default: 0.01,0.1,1,10,100.\n",
    "  --chronos-models        Clock models. Default: clock,correlated,relaxed,discrete.\n",
    "  --chronos-discrete-k    nb.rate.cat grid for the discrete model. Default: 5.\n",
    "  --chronos-retries       Retries per chronos setting. Default: 2.\n",
    "  --treepl-smoothing      Smoothing grid. Default: 0.01,0.1,1,10,100.\n",
    "  --treepl-bin            treePL executable. Default: TREEPL_BIN env, PATH treePL, then ../treePL.\n",
    "  --treepl-numsites       numsites value written to treePL configs. Default: 1000.\n",
    "  --treepl-threads        OMP_NUM_THREADS for treePL. Default: 1.\n",
    "  --reltime-sites         Alignment length used for RelTime CI widths. Default: 1000.\n",
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
source(file.path(script_dir, "reltime_helpers.R"))

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
  ape::makeChronosCalib(
    phy,
    node = node_bounds$node,
    age.min = node_bounds$age_min,
    age.max = node_bounds$age_max
  )
}

fit_score_chronos <- function(tr) {
  ph <- attr(tr, "PHIIC")
  if (is.list(ph) && is.finite(ph$PHIIC)) return(ph$PHIIC)
  pl <- attr(tr, "ploglik")
  if (is.finite(pl)) return(-pl)
  Inf
}

run_chronos_one <- function(phy, calib, model, lambda, nb_rate_cat = NA_integer_, retries = 2L) {
  best <- NULL
  last_error <- NA_character_
  retries <- max(1L, as.integer(retries))

  for (attempt in seq_len(retries)) {
    ctrl <- ape::chronos.control()
    if (identical(model, "discrete") && is.finite(nb_rate_cat)) {
      ctrl$nb.rate.cat <- as.integer(nb_rate_cat)
    }
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
        any(!is.finite(tr$edge.length)) || any(tr$edge.length <= 0)) {
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

write_treepl_cfg <- function(cfg_path, tree_path, out_tree_path, smooth, numsites, node_bounds) {
  lines <- c(
    paste0("treefile = ", tree_path),
    paste0("outfile = ", out_tree_path),
    paste0("numsites = ", as.integer(numsites)),
    paste0("smooth = ", fmt_num(smooth))
  )
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

run_treepl_one <- function(treepl_bin, cfg_path, log_path, omp_threads = 1L) {
  if (file.exists(log_path)) unlink(log_path)
  cmd <- sprintf(
    "OMP_NUM_THREADS=%d %s %s > %s 2>&1",
    as.integer(omp_threads),
    shQuote(treepl_bin),
    shQuote(cfg_path),
    shQuote(log_path)
  )
  suppressWarnings(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE))
}

methods <- tolower(parse_chr_grid(kv[["methods"]] %||% "chronos,treepl,reltime", "--methods"))
valid_methods <- c("chronos", "treepl", "reltime")
bad_methods <- setdiff(methods, valid_methods)
if (length(bad_methods)) {
  stop("Unknown methods: ", paste(bad_methods, collapse = ", "))
}

chronos_lambdas <- parse_num_grid(kv[["chronos-lambdas"]] %||% "0.01,0.1,1,10,100", "--chronos-lambdas")
chronos_models <- parse_chr_grid(kv[["chronos-models"]] %||% "clock,correlated,relaxed,discrete", "--chronos-models")
chronos_models <- tolower(chronos_models)
bad_models <- setdiff(chronos_models, c("clock", "correlated", "relaxed", "discrete"))
if (length(bad_models)) stop("Unsupported chronos model(s): ", paste(bad_models, collapse = ", "))
chronos_discrete_k <- parse_int_grid(kv[["chronos-discrete-k"]] %||% "5", "--chronos-discrete-k")
chronos_retries <- as.integer(kv[["chronos-retries"]] %||% "2")
if (!is.finite(chronos_retries) || chronos_retries < 1L) stop("--chronos-retries must be >= 1.")

treepl_smoothing <- parse_num_grid(kv[["treepl-smoothing"]] %||% "0.01,0.1,1,10,100", "--treepl-smoothing")
treepl_numsites <- as.integer(kv[["treepl-numsites"]] %||% "1000")
treepl_threads <- as.integer(kv[["treepl-threads"]] %||% "1")
if (!is.finite(treepl_numsites) || treepl_numsites < 1L) stop("--treepl-numsites must be >= 1.")
if (!is.finite(treepl_threads) || treepl_threads < 1L) stop("--treepl-threads must be >= 1.")

reltime_sites <- as.integer(kv[["reltime-sites"]] %||% "1000")
if (!is.finite(reltime_sites) || reltime_sites < 1L) stop("--reltime-sites must be >= 1.")

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
dir.create(file.path(outdir, "treepl", "trees"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "treepl", "configs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "treepl", "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "reltime"), recursive = TRUE, showWarnings = FALSE)

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

if ("chronos" %in% methods) {
  msg("Running chronos grid...")
  chronos_calib <- build_chronos_calib_from_node_bounds(phy, node_bounds)
  chronos_rows <- list()
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
        run <- run_chronos_one(
          phy = phy,
          calib = chronos_calib,
          model = mdl,
          lambda = lam,
          nb_rate_cat = kcat,
          retries = chronos_retries
        )
        if (identical(run$status, "OK")) {
          ape::write.tree(run$tree, file = out_tree)
          candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
            candidate = candidate,
            tree_file = out_tree,
            stringsAsFactors = FALSE
          )
        } else {
          out_tree <- ""
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
          error = run$error %||% NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  chronos_df <- do.call(rbind, chronos_rows)
  write.csv(chronos_df, file.path(outdir, "chronos", paste0(prefix, "_chronos_runs.csv")), row.names = FALSE)
  all_rows[[length(all_rows) + 1L]] <- chronos_df
  msg("chronos finished: ", sum(chronos_df$status == "OK"), "/", nrow(chronos_df), " successful.")
}

if ("treepl" %in% methods) {
  treepl_bin <- resolve_treepl_bin(kv[["treepl-bin"]] %||% NULL, repo_dir)
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
    if (file.exists(out_tree)) unlink(out_tree)
    write_treepl_cfg(cfg_path, treepl_input_tree, out_tree, smooth, treepl_numsites, node_bounds)
    exit_code <- run_treepl_one(treepl_bin, cfg_path, log_path, omp_threads = treepl_threads)

    status <- "FAILED"
    err <- tail_text(log_path)
    if (identical(exit_code, 0L) && file.exists(out_tree)) {
      tr <- try(ape::read.tree(out_tree), silent = TRUE)
      if (!inherits(tr, "try-error") && inherits(tr, "phylo") &&
          !is.null(tr$edge.length) && all(is.finite(tr$edge.length)) && all(tr$edge.length > 0)) {
        status <- "OK"
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
      status = status,
      exit_code = exit_code,
      tree_file = if (identical(status, "OK")) out_tree else "",
      cfg_file = cfg_path,
      log_file = log_path,
      error = err %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }

  treepl_df <- do.call(rbind, treepl_rows)
  write.csv(treepl_df, file.path(outdir, "treepl", paste0(prefix, "_treepl_runs.csv")), row.names = FALSE)
  all_rows[[length(all_rows) + 1L]] <- treepl_df
  msg("treePL finished: ", sum(treepl_df$status == "OK"), "/", nrow(treepl_df), " successful.")
}

if ("reltime" %in% methods) {
  msg("Running RelTime...")
  rel_candidate <- "RelTime"
  rel_tree_file <- file.path(outdir, "reltime", paste0(prefix, "_RelTime.tre"))
  rel_ci_file <- file.path(outdir, "reltime", paste0(prefix, "_RelTime_ci.csv"))
  rel_bounds_file <- file.path(outdir, "reltime", paste0(prefix, "_RelTime_bounds_used.csv"))
  rel_summary_file <- file.path(outdir, "reltime", paste0(prefix, "_RelTime_run.csv"))

  rel_run <- try(
    run_reltime_with_bounds_ci(
      phy = phy,
      calibration_df = pair_df[, c("taxonA", "taxonB", "age_min", "age_max"), drop = FALSE],
      root_age = if (is.finite(root_age)) root_age else NULL,
      n_sites = reltime_sites
    ),
    silent = TRUE
  )

  if (!inherits(rel_run, "try-error")) {
    ape::write.tree(rel_run$tree, file = rel_tree_file)
    write.csv(rel_run$ci, rel_ci_file, row.names = FALSE)
    write.csv(rel_run$bounds, rel_bounds_file, row.names = FALSE)
    rel_df <- data.frame(
      candidate = rel_candidate,
      method = "RelTime",
      model = "RelTime",
      lambda = NA_real_,
      nb_rate_cat = NA_integer_,
      smooth = NA_real_,
      status = "OK",
      tree_file = rel_tree_file,
      ci_file = rel_ci_file,
      bounds_file = rel_bounds_file,
      n_sites = reltime_sites,
      error = NA_character_,
      stringsAsFactors = FALSE
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
      bounds_file = "",
      n_sites = reltime_sites,
      error = conditionMessage(attr(rel_run, "condition")),
      stringsAsFactors = FALSE
    )
    msg("RelTime failed: ", rel_df$error[1L])
  }

  write.csv(rel_df, rel_summary_file, row.names = FALSE)
  all_rows[[length(all_rows) + 1L]] <- rel_df
}

all_df <- bind_rows_fill(all_rows)
all_runs_csv <- file.path(outdir, paste0(prefix, "_all_runs_summary.csv"))
write.csv(all_df, all_runs_csv, row.names = FALSE)

candidates_df <- if (length(candidate_rows)) do.call(rbind, candidate_rows) else data.frame(candidate = character(0), tree_file = character(0), stringsAsFactors = FALSE)
if (nrow(candidates_df)) {
  candidates_df <- candidates_df[!duplicated(candidates_df$candidate), , drop = FALSE]
}
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
  paste0("Calibration pairs input: ", nrow(pair_df)),
  paste0("Calibration pairs mapped: ", nrow(pair_map$mapped)),
  paste0("Calibration bounds retained: ", nrow(node_bounds)),
  paste0("Calibration bounds dropped as inconsistent: ", nrow(pair_map$dropped_inconsistent)),
  paste0("Successful candidates written: ", nrow(candidates_df)),
  paste0("Combined runs summary: ", all_runs_csv),
  paste0("Candidates CSV: ", candidates_csv)
)
writeLines(meta_lines, file.path(outdir, paste0(prefix, "_run_metadata.txt")))

msg("Saved combined runs summary: ", all_runs_csv)
msg("Saved candidates CSV: ", candidates_csv)
msg("Done.")

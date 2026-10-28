##############################################################
## reltime_mega_helpers.R
##
## RelTime support through MEGA-CC (`megacc`). MEGA is the only
## RelTime engine in this repo.
##
## - run_megacc_reltime(): one RelTime run from branch lengths
## - megacc_tao_ci():      run + parse MEGA's analytical (Tao et al.
##                         2020) 95% CIs from the NEXUS output
## - megacc_bootstrap_ci(): Poisson branch-length bootstrap; each
##                         replicate reruns megacc on a perturbed
##                         phylogram (same scheme as chronos/treePL)
##############################################################

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

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

write_mega_outgroup <- function(outgroup_path, outgroup_tip = "MOCK_OUTGROUP_FOR_MEGA") {
  writeLines(paste0(outgroup_tip, "=outgroup"), outgroup_path)
}

## All node bounds are written, including the ingroup root. The mock outgroup
## added by add_mock_mega_outgroup() makes the ingroup root an ordinary
## internal node, so MEGA can apply the root calibration.
write_mega_calibrations <- function(cal_path, node_bounds) {
  lines <- character(nrow(node_bounds))
  for (i in seq_len(nrow(node_bounds))) {
    tag <- paste0("cal_", node_bounds$node[i])
    if (is.finite(node_bounds$age_max[i])) {
      lines[i] <- sprintf(
        "!MRCA='%s' TaxonA='%s' TaxonB='%s' MinTime=%.11f MaxTime=%.11f calibrationName='%s';",
        tag, node_bounds$taxonA[i], node_bounds$taxonB[i],
        node_bounds$age_min[i], node_bounds$age_max[i], tag
      )
    } else {
      lines[i] <- sprintf(
        "!MRCA='%s' TaxonA='%s' TaxonB='%s' MinTime=%.11f calibrationName='%s';",
        tag, node_bounds$taxonA[i], node_bounds$taxonB[i],
        node_bounds$age_min[i], tag
      )
    }
  }
  writeLines(lines, cal_path)
}

add_mock_mega_outgroup <- function(phy, mock_tip = "MOCK_OUTGROUP_FOR_MEGA") {
  root_depth <- max(ape::node.depth.edgelength(phy))
  nwk <- sub(";$", "", ape::write.tree(phy))
  new_nwk <- paste0("(", nwk, ":0.001,", mock_tip, ":",
                    format(root_depth + 0.001, scientific = FALSE), ");")
  ape::read.tree(text = new_nwk)
}

run_megacc_reltime <- function(megacc_bin, phy, node_bounds, work_dir, prefix) {
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  ## MEGA-CC is fragile with long/deep file paths on some empirical trees.
  ## Run from a short temp directory, then copy artifacts back into work_dir.
  mega_work <- tempfile(pattern = paste0("megacc_", prefix, "_"), tmpdir = tempdir())
  dir.create(mega_work, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(mega_work, recursive = TRUE, force = TRUE), add = TRUE)

  mao_path      <- file.path(mega_work, "reltime_blens.mao")
  tree_path     <- file.path(mega_work, paste0(prefix, "_phylogram_for_mega.nwk"))
  outgroup_path <- file.path(mega_work, "outgroup.txt")
  cal_path      <- file.path(mega_work, "mega_calibrations.txt")
  out_prefix    <- file.path(mega_work, paste0(prefix, "_mega_reltime"))
  log_path      <- file.path(mega_work, paste0(prefix, "_mega_reltime.log"))
  mock_tip      <- "MOCK_OUTGROUP_FOR_MEGA"

  write_mega_mao(mao_path)
  ape::write.tree(add_mock_mega_outgroup(phy, mock_tip = mock_tip), file = tree_path)
  write_mega_outgroup(outgroup_path, outgroup_tip = mock_tip)
  write_mega_calibrations(cal_path, node_bounds)

  exit_code <- suppressWarnings(system2(
    megacc_bin,
    args = c("-a", mao_path, "-t", tree_path,
             "-c", cal_path, "-g", outgroup_path,
             "-o", out_prefix),
    stdout = log_path,
    stderr = log_path
  ))

  staged_files <- list.files(mega_work, all.files = FALSE, full.names = TRUE, no.. = TRUE)
  if (length(staged_files)) {
    invisible(file.copy(staged_files, work_dir, overwrite = TRUE, recursive = TRUE))
  }

  exact_tree_file <- paste0(out_prefix, "_exactTimes.nwk")
  if (!file.exists(exact_tree_file)) {
    log_txt <- if (file.exists(log_path)) paste(readLines(log_path, warn = FALSE), collapse = "\n") else ""
    return(list(tree = NULL, tree_file = "", nexus_file = "", status = "FAILED",
                error = paste0("megacc exit code ", exit_code, "; ", log_txt)))
  }

  mega_tree <- ape::read.tree(exact_tree_file)
  ## MEGA wraps labels in quotes and turns underscores into spaces.
  mega_tree$tip.label <- gsub("^'|'$", "", mega_tree$tip.label)
  mega_tree$tip.label <- gsub(" ", "_", mega_tree$tip.label)
  if (mock_tip %in% mega_tree$tip.label) {
    mega_tree <- ape::drop.tip(mega_tree, mock_tip)
  }
  out_tree <- file.path(work_dir, paste0(prefix, "_RelTime_MEGA.tre"))
  ape::write.tree(mega_tree, file = out_tree)
  list(
    tree = mega_tree,
    tree_file = out_tree,
    nexus_file = file.path(work_dir, paste0(basename(out_prefix), "_nexus.tre")),
    status = "OK",
    error = NA_character_
  )
}

megacc_tao_ci <- function(megacc_bin, phy, node_bounds, work_dir, prefix) {
  base <- run_megacc_reltime(megacc_bin, phy, node_bounds, work_dir, prefix)
  if (base$status != "OK") return(base)

  mega_tree <- base$tree
  n_tip  <- ape::Ntip(mega_tree)
  n_node <- ape::Nnode(mega_tree)
  int_nodes <- seq.int(n_tip + 1L, n_tip + n_node)

  ## ---- Parse native Tao 95% CIs from MEGA's NEXUS output ----
  if (!file.exists(base$nexus_file)) {
    warning("No MEGA NEXUS output found — skipping Tao CIs")
    base$ci <- NULL
    return(base)
  }
  nexus_lines <- readLines(base$nexus_file, warn = FALSE)
  tree_line <- grep("^tree TREE1", nexus_lines, value = TRUE)
  if (!length(tree_line)) {
    warning("No 'tree TREE1' line in MEGA NEXUS — skipping Tao CIs")
    base$ci <- NULL
    return(base)
  }

  ## MEGA annotates internal nodes as [&rate=...,divtime=T,divtime_95%_CI={hi,lo}].
  ## Annotations are matched to the cleaned tree by divtime value.
  ci_raw <- regmatches(tree_line, gregexpr("divtime_95%_CI=[{][^}]+[}]", tree_line))[[1]]
  dt_raw <- regmatches(tree_line, gregexpr("divtime=([0-9eE.+-]+)", tree_line))[[1]]
  dt_vals <- as.numeric(gsub("divtime=", "", dt_raw))

  if (length(ci_raw) != length(dt_raw) || length(ci_raw) == 0L) {
    warning("Could not parse MEGA Tao CIs (", length(ci_raw), " CIs, ",
            length(dt_raw), " divtimes)")
    base$ci <- NULL
    return(base)
  }

  ci_split <- strsplit(gsub("divtime_95%_CI=[{]|[}]", "", ci_raw), ",")
  ci_lo <- vapply(ci_split, function(x) min(as.numeric(x)), numeric(1))
  ci_hi <- vapply(ci_split, function(x) max(as.numeric(x)), numeric(1))

  int_ages <- .mega_node_ages(mega_tree)[int_nodes]
  matched_ci_lo <- rep(NA_real_, n_node)
  matched_ci_hi <- rep(NA_real_, n_node)
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

  ## The root usually carries no CI annotation: set CI = point estimate.
  root_idx <- which(int_nodes == n_tip + 1L)
  if (is.na(matched_ci_lo[root_idx])) {
    matched_ci_lo[root_idx] <- int_ages[root_idx]
    matched_ci_hi[root_idx] <- int_ages[root_idx]
  }

  ## Near-zero-age nodes (MEGA's collapsed backbone branches) can get Tao CIs
  ## that span the full tree depth. Set their CIs to the point estimate.
  near_zero <- int_ages < 1e-4
  matched_ci_lo[near_zero] <- int_ages[near_zero]
  matched_ci_hi[near_zero] <- int_ages[near_zero]

  base$ci <- data.frame(
    node  = int_nodes,
    age   = int_ages,
    ci_lo = matched_ci_lo,
    ci_hi = matched_ci_hi,
    se    = (matched_ci_hi - matched_ci_lo) / (2 * 1.96),
    stringsAsFactors = FALSE
  )
  base
}

.mega_node_ages <- function(tr) {
  d <- ape::node.depth.edgelength(tr)
  ages <- max(d[seq_len(ape::Ntip(tr))]) - d
  ages[seq_len(ape::Ntip(tr))] <- 0
  ages
}

## Internal-node signatures (sorted descendant tip labels), named by node number.
.mega_signature_index <- function(tr) {
  parts <- ape::prop.part(tr)
  labs <- attr(parts, "labels")
  sig <- vapply(parts, function(p) paste(sort(labs[p]), collapse = "\r"), character(1))
  names(sig) <- as.character(seq.int(ape::Ntip(tr) + 1L, length.out = length(parts)))
  sig
}

.mega_ages_by_signature <- function(tr, ref_sig) {
  sig <- .mega_signature_index(tr)
  ages <- .mega_node_ages(tr)[as.integer(names(sig))]
  unname(stats::setNames(ages, sig)[ref_sig])
}

.mega_poisson_perturb_phy <- function(phy, n_sites, min_edge = 1e-12) {
  out <- phy
  lambda <- pmax(as.numeric(phy$edge.length) * n_sites, 0)
  out$edge.length <- pmax(stats::rpois(length(lambda), lambda) / n_sites, min_edge)
  out
}

## Poisson branch-length bootstrap (Paradis et al. 2023) for MEGA RelTime.
## Each replicate perturbs the phylogram, reruns megacc with the same node
## bounds, and collects node ages by descendant tip set. CI nodes are numbered
## on `mega_tree` (the point-estimate tree as read from file).
megacc_bootstrap_ci <- function(megacc_bin, phy, node_bounds, mega_tree,
                                B = 100L, n_sites = 1000L, min_edge = 1e-12,
                                prefix = "reltime_boot", quiet = TRUE) {
  B <- as.integer(B)
  if (!is.finite(B) || B < 1L) stop("B must be a positive integer.")
  if (!is.finite(n_sites) || n_sites <= 0) stop("n_sites must be positive for the RelTime bootstrap.")

  ref_sig <- .mega_signature_index(mega_tree)
  int_nodes <- as.integer(names(ref_sig))
  base_age <- .mega_node_ages(mega_tree)[int_nodes]

  boot_age <- matrix(NA_real_, nrow = length(int_nodes), ncol = B)
  max_attempts <- max(5L * B, B + 20L)
  n_ok <- 0L
  attempt_i <- 0L
  while (n_ok < B && attempt_i < max_attempts) {
    attempt_i <- attempt_i + 1L
    if (!quiet) cat("\rRunning MEGA RelTime bootstrap: ", n_ok + 1L, "/", B, sep = "")
    phy_boot <- .mega_poisson_perturb_phy(phy, n_sites, min_edge = min_edge)
    rep_dir <- tempfile(pattern = sprintf("%s_boot_%03d_", prefix, attempt_i))
    res <- try(run_megacc_reltime(megacc_bin, phy_boot, node_bounds, rep_dir, prefix), silent = TRUE)
    unlink(rep_dir, recursive = TRUE, force = TRUE)
    if (inherits(res, "try-error") || !identical(res$status, "OK")) next
    ages <- .mega_ages_by_signature(res$tree, ref_sig)
    if (!all(is.finite(ages))) next
    n_ok <- n_ok + 1L
    boot_age[, n_ok] <- ages
  }
  if (!quiet) cat("\n")
  if (n_ok < B) {
    stop("MEGA RelTime bootstrap produced only ", n_ok, " successful replicates out of ",
         attempt_i, " attempts.")
  }

  qfun <- function(probs) {
    apply(boot_age, 1L, stats::quantile, probs = probs, names = FALSE, na.rm = TRUE, type = 8)
  }
  data.frame(
    node = int_nodes,
    age = base_age,
    ci_lo = qfun(0.025),
    ci_hi = qfun(0.975),
    q50_lo = qfun(0.25),
    q50_hi = qfun(0.75),
    se = NA_real_,
    stringsAsFactors = FALSE
  )
}

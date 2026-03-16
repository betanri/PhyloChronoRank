#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
kv <- list()
for (a in args) {
  if (!grepl("^--", a)) next
  parts <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  key <- parts[1]
  val <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "TRUE"
  kv[[key]] <- val
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else getwd()
if (!nzchar(script_dir) || !dir.exists(script_dir)) script_dir <- getwd()
base_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

source(file.path(script_dir, "reltime_helpers.R"))
source(file.path(script_dir, "chronos_ci_helpers.R"))
source(file.path(script_dir, "treepl_bootstrap_helpers.R"))

resolve_vgp_dir <- function() {
  raw <- kv[["vgp-dir"]] %||% Sys.getenv("VGP_DIR", unset = "")
  if (!nzchar(raw)) {
    sibling <- file.path(base_dir, "..", "chronos", "VGP")
    if (dir.exists(sibling)) raw <- sibling
  }
  if (!nzchar(raw)) {
    stop("Pass --vgp-dir=/PATH/TO/VGP_FOLDER or set VGP_DIR=/PATH/TO/VGP_FOLDER.")
  }
  normalizePath(raw, winslash = "/", mustWork = TRUE)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
seed <- suppressWarnings(as.integer(kv[["seed"]] %||% "20260315"))
if (is.finite(seed)) set.seed(seed)
ci_sites <- as.integer(kv[["ci-sites"]] %||% "1000")
if (!is.finite(ci_sites) || ci_sites < 1L) stop("--ci-sites must be >= 1.")
chronos_ci_reps <- as.integer(kv[["chronos-ci-reps"]] %||% "100")
if (!is.finite(chronos_ci_reps) || chronos_ci_reps < 1L) stop("--chronos-ci-reps must be >= 1.")
reltime_bootstrap_reps <- as.integer(kv[["reltime-bootstrap-reps"]] %||% "100")
if (!is.finite(reltime_bootstrap_reps) || reltime_bootstrap_reps < 1L) {
  stop("--reltime-bootstrap-reps must be >= 1.")
}
treepl_bootstrap_reps <- as.integer(kv[["treepl-bootstrap-reps"]] %||% "100")
if (!is.finite(treepl_bootstrap_reps) || treepl_bootstrap_reps < 1L) {
  stop("--treepl-bootstrap-reps must be >= 1.")
}
treepl_bootstrap_jobs <- as.integer(kv[["treepl-bootstrap-jobs"]] %||% "4")
if (!is.finite(treepl_bootstrap_jobs) || treepl_bootstrap_jobs < 1L) {
  stop("--treepl-bootstrap-jobs must be >= 1.")
}
treepl_numsites <- as.integer(kv[["treepl-numsites"]] %||% "1000")
if (!is.finite(treepl_numsites) || treepl_numsites < 1L) stop("--treepl-numsites must be >= 1.")
treepl_threads <- as.integer(kv[["treepl-threads"]] %||% "1")
if (!is.finite(treepl_threads) || treepl_threads < 1L) stop("--treepl-threads must be >= 1.")
treepl_thorough <- !identical(tolower(kv[["treepl-thorough"]] %||% "TRUE"), "false")
treepl_prime <- !identical(tolower(kv[["treepl-prime"]] %||% "TRUE"), "false")
treepl_opt <- kv[["treepl-opt"]] %||% NULL
treepl_plsimaniter <- kv[["treepl-plsimaniter"]] %||% NULL
treepl_pliter <- kv[["treepl-pliter"]] %||% NULL
treepl_bin <- treepl_resolve_bin(kv[["treepl-bin"]] %||% NULL, base_dir)
if (!nzchar(treepl_bin)) stop("treePL binary not found. Pass --treepl-bin or set TREEPL_BIN.")

build_ci_inputs_local <- function(phy, calibration_df) {
  bounds <- reltime_merge_calibration_bounds(phy, calibration_df)
  root_node <- Ntip(phy) + 1L
  calib_df <- data.frame(
    node = bounds$node,
    age.min = bounds$age_min,
    age.max = bounds$age_max,
    stringsAsFactors = FALSE
  )
  root_var <- 0
  root_row <- calib_df[calib_df$node == root_node, , drop = FALSE]
  if (nrow(root_row) &&
      is.finite(root_row$age.min[1L]) &&
      is.finite(root_row$age.max[1L]) &&
      abs(root_row$age.max[1L] - root_row$age.min[1L]) > 1e-10) {
    root_var <- ((root_row$age.max[1L] - root_row$age.min[1L]) / (2 * 1.96))^2
  }
  list(
    bounds = bounds,
    root_var = root_var,
    cal_min = as.list(stats::setNames(calib_df$age.min, calib_df$node)),
    cal_max = as.list(stats::setNames(calib_df$age.max, calib_df$node)),
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

append_candidate_row <- function(cand, candidate_name, tree_file) {
  if (!all(c("candidate", "tree_file") %in% names(cand))) {
    stop("Candidate table must contain candidate and tree_file columns.")
  }
  out <- cand[cand$candidate != candidate_name, , drop = FALSE]
  new_row <- out[0, , drop = FALSE]
  new_row[1, ] <- NA
  new_row$candidate <- candidate_name
  new_row$tree_file <- tree_file
  rbind(out, new_row)
}

sanitize_tree_file_column <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if ("tree_file" %in% names(df)) {
    df$tree_file <- file.path("/PATH/TO/FOLDER", basename(df$tree_file))
  }
  write.csv(df, path, row.names = FALSE)
}

write_reltime_outputs <- function(ref_tree, calibration_df, tree_file,
                                  bootstrap_ci_file, tao_ci_file, bounds_file) {
  rel_boot <- reltime_bootstrap_ci(
    phy = ref_tree,
    calibration_df = calibration_df,
    B = reltime_bootstrap_reps,
    n_sites = ci_sites,
    quiet = TRUE
  )
  rel_tao <- reltime_tao_ci_from_bounds_run(
    phy = ref_tree,
    rel_run = rel_boot,
    n_sites = ci_sites
  )
  write.tree(rel_boot$tree, tree_file)
  write.csv(rel_boot$ci, bootstrap_ci_file, row.names = FALSE)
  write.csv(rel_tao, tao_ci_file, row.names = FALSE)
  write.csv(rel_boot$bounds, bounds_file, row.names = FALSE)
  list(tree = rel_boot$tree, bootstrap_ci = rel_boot$ci, tao_ci = rel_tao, bounds = rel_boot$bounds)
}

write_treepl_outputs <- function(ref_tree, calibration_df, tree_file,
                                 bootstrap_ci_file, smooth, prefix) {
  treepl_boot <- treepl_bootstrap_ci(
    treepl_bin = treepl_bin,
    phy = ref_tree,
    calibration_df = calibration_df,
    smooth = smooth,
    dated_tree = tree_file,
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
    workdir = tempfile(pattern = paste0(prefix, "_treepl_boot_")),
    prefix = prefix
  )
  write.csv(treepl_boot$ci, bootstrap_ci_file, row.names = FALSE)
  treepl_boot
}

load_terap_chronos_settings <- function() {
  fit_path <- file.path(base_dir, "examples", "terapontoid", "summary_terap_empirical_model_fits.csv")
  fit_df <- read.csv(fit_path, stringsAsFactors = FALSE)
  fit_df$tree_key <- basename(fit_df$dated_tree_file)
  split(
    data.frame(
      model = fit_df$model,
      lambda = as.numeric(fit_df$best_phiic_lambda),
      nb_rate_cat = as.numeric(fit_df$best_phiic_nb_rate_cat),
      stringsAsFactors = FALSE
    ),
    fit_df$tree_key
  )
}

load_vgp_chronos_settings <- function(vgp_dir) {
  fit_path <- file.path(vgp_dir, "chronos_empirical_out_full_per_model", "summary_fulltree_per_model_selected.csv")
  fit_df <- read.csv(fit_path, stringsAsFactors = FALSE)
  fit_df$tree_key <- basename(fit_df$out_tree)
  split(
    data.frame(
      model = fit_df$model,
      lambda = as.numeric(fit_df$lambda),
      nb_rate_cat = as.numeric(fit_df$nb_rate_cat),
      stringsAsFactors = FALSE
    ),
    fit_df$tree_key
  )
}

get_chronos_settings <- function(settings_map, tree_path) {
  tree_key <- basename(tree_path)
  hit <- settings_map[[tree_key]]
  if (is.null(hit) || !nrow(hit)) {
    stop("No chronos settings found for tree: ", tree_key)
  }
  hit[1, , drop = FALSE]
}

write_terap_outputs <- function() {
  ex_dir <- file.path(base_dir, "examples", "terapontoid")
  ref_tree <- read.tree(file.path(ex_dir, "Terapontoid_ML_MAIN_phylogram_used.tree"))
  cal_df <- read.csv(file.path(ex_dir, "Terapontoid_ML_MAIN_calibrations_used.csv"), stringsAsFactors = FALSE)
  ci_inputs <- build_ci_inputs_local(ref_tree, cal_df)
  rel_tree_file <- file.path(ex_dir, "Terapontoid_ML_MAIN_RelTime_full_bounds.tre")
  rel_boot_ci_file <- file.path(ex_dir, "Terapontoid_ML_MAIN_RelTime_full_bounds_bootstrap_ci.csv")
  rel_tao_ci_file <- file.path(ex_dir, "Terapontoid_ML_MAIN_RelTime_full_bounds_ci.csv")
  rel_bounds_file <- file.path(ex_dir, "Terapontoid_ML_MAIN_RelTime_bounds_used.csv")
  treepl_tree_file <- file.path(ex_dir, "Terapontoid_ML_MAIN_treePL_congruify.tre")
  treepl_boot_ci_file <- file.path(ex_dir, "Terapontoid_ML_MAIN_treePL_congruify_bootstrap_ci.csv")
  rel_outputs <- write_reltime_outputs(
    ref_tree = ref_tree,
    calibration_df = cal_df,
    tree_file = rel_tree_file,
    bootstrap_ci_file = rel_boot_ci_file,
    tao_ci_file = rel_tao_ci_file,
    bounds_file = rel_bounds_file
  )
  treepl_outputs <- write_treepl_outputs(
    ref_tree = ref_tree,
    calibration_df = cal_df,
    tree_file = treepl_tree_file,
    bootstrap_ci_file = treepl_boot_ci_file,
    smooth = 0.01,
    prefix = "terap_treepl"
  )
  cand <- read.csv(file.path(ex_dir, "candidates.csv"), stringsAsFactors = FALSE)
  chronos_settings <- load_terap_chronos_settings()

  uncertainty_rows <- list()
  for (i in seq_len(nrow(cand))) {
    candidate <- cand$candidate[i]
    tree_path <- file.path(ex_dir, cand$tree_file[i])
    if (candidate == "RelTime") {
      ci_df <- rel_outputs$bootstrap_ci
      source_tag <- "RelTime-bootstrap"
    } else if (identical(candidate, "treePL")) {
      ci_df <- treepl_outputs$ci
      source_tag <- "treePL-bootstrap"
    } else if (startsWith(candidate, "chronos_")) {
      tr <- read.tree(tree_path)
      meta <- get_chronos_settings(chronos_settings, tree_path)
      ci_path <- sub("\\.tre$", "_ci.csv", tree_path)
      ci_try <- try(
        chronos_ci(
          ref_tree, tr,
          calibration = ci_inputs$calib_bundle$chronos_calib,
          B = chronos_ci_reps,
          type = kv[["chronos-ci-type"]] %||% "parametric",
          n_sites = ci_sites,
          model = meta$model[1],
          lambda = meta$lambda[1],
          nb_rate_cat = meta$nb_rate_cat[1],
          quiet = TRUE
        ),
        silent = TRUE
      )
      source_tag <- "chronos-bootstrap"
      if (inherits(ci_try, "try-error")) {
        ci_df <- data.frame(ci_lo = numeric(), ci_hi = numeric())
      } else {
        ci_df <- ci_try
        write.csv(ci_df, ci_path, row.names = FALSE)
      }
    } else {
      ci_df <- data.frame(ci_lo = numeric(), ci_hi = numeric())
      source_tag <- "not_scored"
    }
    if (!nrow(ci_df)) {
      uncertainty_rows[[length(uncertainty_rows) + 1L]] <- make_unscored_uncertainty_row(candidate, source = source_tag)
    } else {
      uncertainty_rows[[length(uncertainty_rows) + 1L]] <- ci_width_summary(candidate, ci_df, source = source_tag)
    }
  }

  uncertainty_file <- file.path(ex_dir, "uncertainty_summary_long.csv")
  write.csv(do.call(rbind, uncertainty_rows), uncertainty_file, row.names = FALSE)

  outdir <- tempfile(pattern = "pcr_terap_uncertainty_")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    "Rscript",
    c(
      file.path(base_dir, "scripts", "run_pcr.R"),
      paste0("--ref-tree=", file.path(ex_dir, "Terapontoid_ML_MAIN_phylogram_used.tree")),
      paste0("--candidates-csv=", file.path(ex_dir, "candidates.csv")),
      paste0("--uncertainty-csv=", uncertainty_file),
      paste0("--outdir=", outdir)
    )
  )
  if (!identical(status, 0L)) stop("Terapontoidei PCR rerun failed.")

  file.copy(
    file.path(outdir, "summary_pcr_metrics.csv"),
    file.path(ex_dir, "summary_terap_empirical_postfit_metrics.csv"),
    overwrite = TRUE
  )
}

write_vgp_outputs <- function() {
  vgp_dir <- resolve_vgp_dir()
  ref_tree_file <- file.path(vgp_dir, "roadies_v1.1.16b.numbers.nwk")
  cal_file <- file.path(vgp_dir, "roadies_manual_calibrations.csv")
  rel_tree_file <- file.path(vgp_dir, "pcr_rerun_20260310", "roadies_v1.1.16b.numbers_RelTime_full_bounds.tre")
  rel_boot_ci_file <- file.path(vgp_dir, "pcr_rerun_20260310", "roadies_v1.1.16b.numbers_RelTime_full_bounds_bootstrap_ci.csv")
  rel_tao_ci_file <- file.path(vgp_dir, "pcr_rerun_20260310", "roadies_v1.1.16b.numbers_RelTime_full_bounds_ci.csv")
  rel_bounds_file <- file.path(vgp_dir, "pcr_rerun_20260310", "roadies_v1.1.16b.numbers_RelTime_bounds_used.csv")

  ref_tree <- read.tree(ref_tree_file)
  cal_df <- read.csv(cal_file, stringsAsFactors = FALSE)
  rel_outputs <- write_reltime_outputs(
    ref_tree = ref_tree,
    calibration_df = cal_df,
    tree_file = rel_tree_file,
    bootstrap_ci_file = rel_boot_ci_file,
    tao_ci_file = rel_tao_ci_file,
    bounds_file = rel_bounds_file
  )
  ci_inputs <- build_ci_inputs_local(ref_tree, cal_df)
  cand_file_raw <- kv[["vertebrate-candidates-csv"]] %||% ""
  if (nzchar(cand_file_raw)) {
    cand <- read.csv(normalizePath(cand_file_raw, winslash = "/", mustWork = TRUE), stringsAsFactors = FALSE)
  } else {
    cand_file <- file.path(vgp_dir, "pcr_rerun_20260310", "candidates_core5.csv")
    cand <- read.csv(cand_file, stringsAsFactors = FALSE)
    cand <- append_candidate_row(cand, "RelTime", rel_tree_file)
  }
  treepl_row <- cand[grepl("^treepl_best-smooth-", cand$candidate, ignore.case = TRUE), , drop = FALSE]
  if (!nrow(treepl_row)) {
    stop("Could not identify the selected treePL candidate in the vertebrate manifest.")
  }
  treepl_smooth <- as.numeric(sub("^treepl_best-smooth-", "", treepl_row$candidate[1]))
  if (!is.finite(treepl_smooth)) stop("Could not parse vertebrate treePL smoothing value.")
  treepl_boot_ci_file <- file.path(vgp_dir, "pcr_rerun_20260310", "roadies_v1.1.16b.numbers_treepl_dated_best_bootstrap_ci.csv")
  treepl_outputs <- write_treepl_outputs(
    ref_tree = ref_tree,
    calibration_df = cal_df,
    tree_file = treepl_row$tree_file[1],
    bootstrap_ci_file = treepl_boot_ci_file,
    smooth = treepl_smooth,
    prefix = "vgp_treepl"
  )
  chronos_settings <- load_vgp_chronos_settings(vgp_dir)

  uncertainty_rows <- list()
  for (i in seq_len(nrow(cand))) {
    candidate <- cand$candidate[i]
    tree_path <- cand$tree_file[i]
    if (candidate == "RelTime") {
      ci_df <- rel_outputs$bootstrap_ci
      source_tag <- "RelTime-bootstrap"
    } else if (grepl("^treepl_best-smooth-", candidate, ignore.case = TRUE)) {
      ci_df <- treepl_outputs$ci
      source_tag <- "treePL-bootstrap"
    } else if (startsWith(candidate, "chronos_")) {
      tr <- read.tree(tree_path)
      meta <- get_chronos_settings(chronos_settings, tree_path)
      ci_try <- try(
        chronos_ci(
          ref_tree, tr,
          calibration = ci_inputs$calib_bundle$chronos_calib,
          B = chronos_ci_reps,
          type = kv[["chronos-ci-type"]] %||% "parametric",
          n_sites = ci_sites,
          model = meta$model[1],
          lambda = meta$lambda[1],
          nb_rate_cat = meta$nb_rate_cat[1],
          quiet = TRUE
        ),
        silent = TRUE
      )
      if (inherits(ci_try, "try-error")) {
        ci_df <- data.frame(ci_lo = numeric(), ci_hi = numeric())
        source_tag <- "not_scored"
      } else {
        ci_df <- ci_try
        source_tag <- "chronos-bootstrap"
      }
    } else {
      ci_df <- data.frame(ci_lo = numeric(), ci_hi = numeric())
      source_tag <- "not_scored"
    }
    if (!nrow(ci_df)) {
      uncertainty_rows[[length(uncertainty_rows) + 1L]] <- make_unscored_uncertainty_row(candidate, source = source_tag)
    } else {
      uncertainty_rows[[length(uncertainty_rows) + 1L]] <- ci_width_summary(candidate, ci_df, source = source_tag)
    }
  }

  out_uncertainty_file <- file.path(
    base_dir, "examples", "unpublished_vertebrate", "postfit_metrics", "uncertainty_summary_long.csv"
  )
  write.csv(do.call(rbind, uncertainty_rows), out_uncertainty_file, row.names = FALSE)

  cand_tmp <- tempfile(pattern = "vgp_candidates_", fileext = ".csv")
  write.csv(cand, cand_tmp, row.names = FALSE, quote = FALSE)

  outdir <- tempfile(pattern = "pcr_vgp_uncertainty_")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    "Rscript",
    c(
      file.path(base_dir, "scripts", "run_pcr.R"),
      paste0("--ref-tree=", ref_tree_file),
      paste0("--candidates-csv=", cand_tmp),
      paste0("--calibrations-csv=", cal_file),
      paste0("--uncertainty-csv=", out_uncertainty_file),
      paste0("--outdir=", outdir)
    )
  )
  if (!identical(status, 0L)) stop("Unpublished vertebrate PCR rerun failed.")

  file.copy(
    file.path(outdir, "summary_pcr_metrics.csv"),
    file.path(base_dir, "examples", "unpublished_vertebrate", "postfit_metrics", "summary_unpublished_vertebrate_postfit_metrics.csv"),
    overwrite = TRUE
  )
  sanitize_tree_file_column(
    file.path(base_dir, "examples", "unpublished_vertebrate", "postfit_metrics", "summary_unpublished_vertebrate_postfit_metrics.csv")
  )
}

write_terap_outputs()
write_vgp_outputs()

system2("Rscript", file.path(base_dir, "scripts", "make_terapontoid_postfit_figures.R"))
system2("Rscript", file.path(base_dir, "scripts", "make_unpublished_vertebrate_postfit_figure.R"))

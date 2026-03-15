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

write_terap_outputs <- function() {
  ex_dir <- file.path(base_dir, "examples", "terapontoid")
  ref_tree <- read.tree(file.path(ex_dir, "Terapontoid_ML_MAIN_phylogram_used.tree"))
  cal_df <- read.csv(file.path(ex_dir, "Terapontoid_ML_MAIN_calibrations_used.csv"), stringsAsFactors = FALSE)
  ci_inputs <- build_ci_inputs_local(ref_tree, cal_df)
  cand <- read.csv(file.path(ex_dir, "candidates.csv"), stringsAsFactors = FALSE)

  uncertainty_rows <- list()
  for (i in seq_len(nrow(cand))) {
    candidate <- cand$candidate[i]
    tree_path <- file.path(ex_dir, cand$tree_file[i])
    if (candidate == "RelTime") {
      ci_path <- file.path(ex_dir, "Terapontoid_ML_MAIN_RelTime_full_bounds_ci.csv")
      ci_df <- read.csv(ci_path, stringsAsFactors = FALSE)
      source_tag <- "RelTime-CI"
    } else {
      tr <- read.tree(tree_path)
      ci_path <- sub("\\.tre$", "_ci.csv", tree_path)
      if (candidate == "treePL") {
        ci_df <- treepl_ci(ref_tree, tr, ci_inputs$calib_bundle, n_sites = 1000L)
        source_tag <- "TreePL-CI"
      } else {
        ci_df <- chronos_ci(
          ref_tree, tr,
          n_sites = 1000L,
          root_var = ci_inputs$root_var,
          cal_min = ci_inputs$cal_min,
          cal_max = ci_inputs$cal_max
        )
        source_tag <- "ChronosCI"
      }
      write.csv(ci_df, ci_path, row.names = FALSE)
    }
    uncertainty_rows[[length(uncertainty_rows) + 1L]] <- ci_width_summary(candidate, ci_df, source = source_tag)
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
  rel_ci_file <- file.path(vgp_dir, "pcr_rerun_20260310", "roadies_v1.1.16b.numbers_RelTime_full_bounds_ci.csv")

  ref_tree <- read.tree(ref_tree_file)
  cal_df <- read.csv(cal_file, stringsAsFactors = FALSE)
  ci_inputs <- build_ci_inputs_local(ref_tree, cal_df)
  cand_file_raw <- kv[["vertebrate-candidates-csv"]] %||% ""
  if (nzchar(cand_file_raw)) {
    cand <- read.csv(normalizePath(cand_file_raw, winslash = "/", mustWork = TRUE), stringsAsFactors = FALSE)
  } else {
    cand_file <- file.path(vgp_dir, "pcr_rerun_20260310", "candidates_core5.csv")
    cand <- read.csv(cand_file, stringsAsFactors = FALSE)
    cand <- append_candidate_row(cand, "RelTime", rel_tree_file)
  }

  uncertainty_rows <- list()
  for (i in seq_len(nrow(cand))) {
    candidate <- cand$candidate[i]
    tree_path <- cand$tree_file[i]
    if (candidate == "RelTime") {
      rel_ci_path <- kv[["vertebrate-reltime-ci"]] %||% sub("\\.tre$", "_ci.csv", tree_path)
      if (!file.exists(rel_ci_path)) rel_ci_path <- rel_ci_file
      ci_df <- read.csv(rel_ci_path, stringsAsFactors = FALSE)
      source_tag <- "RelTime-CI"
    } else {
      tr <- read.tree(tree_path)
      if (grepl("^treepl", candidate, ignore.case = TRUE)) {
        ci_df <- treepl_ci(ref_tree, tr, ci_inputs$calib_bundle, n_sites = 1000L)
        source_tag <- "TreePL-CI"
      } else {
        ci_df <- chronos_ci(
          ref_tree, tr,
          n_sites = 1000L,
          root_var = ci_inputs$root_var,
          cal_min = ci_inputs$cal_min,
          cal_max = ci_inputs$cal_max
        )
        source_tag <- "ChronosCI"
      }
    }
    uncertainty_rows[[length(uncertainty_rows) + 1L]] <- ci_width_summary(candidate, ci_df, source = source_tag)
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

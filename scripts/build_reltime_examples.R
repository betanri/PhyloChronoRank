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

script_dir <- dirname(normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  winslash = "/",
  mustWork = FALSE
))
if (!nzchar(script_dir) || !dir.exists(script_dir)) script_dir <- getwd()
base_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "reltime_helpers.R"))

vgp_dir <- if ("vgp-dir" %in% names(kv)) {
  normalizePath(kv[["vgp-dir"]], winslash = "/", mustWork = TRUE)
} else {
  normalizePath("/Users/ricardobetancur/Desktop/Proxy_Misplaced/chronos/VGP", winslash = "/", mustWork = TRUE)
}

write_terap_outputs <- function() {
  ex_dir <- file.path(base_dir, "examples", "terapontoid")
  ref_tree <- read.tree(file.path(ex_dir, "Terapontoid_ML_MAIN_phylogram_used.tree"))
  cal_df <- read.csv(file.path(ex_dir, "Terapontoid_ML_MAIN_calibrations_used.csv"), stringsAsFactors = FALSE)
  rel <- run_reltime_with_bounds_ci(ref_tree, cal_df)

  tree_file <- file.path(ex_dir, "Terapontoid_ML_MAIN_RelTime_full_bounds.tre")
  ci_file <- file.path(ex_dir, "Terapontoid_ML_MAIN_RelTime_full_bounds_ci.csv")
  bounds_file <- file.path(ex_dir, "Terapontoid_ML_MAIN_RelTime_bounds_used.csv")
  write.tree(rel$tree, tree_file)
  write.csv(rel$ci, ci_file, row.names = FALSE)
  write.csv(rel$bounds, bounds_file, row.names = FALSE)

  cand_file <- file.path(ex_dir, "candidates.csv")
  cand <- read.csv(cand_file, stringsAsFactors = FALSE)
  rel_row <- data.frame(
    candidate = "RelTime",
    tree_file = basename(tree_file),
    stringsAsFactors = FALSE
  )
  cand <- cand[cand$candidate != "RelTime", , drop = FALSE]
  cand <- rbind(cand, rel_row)
  write.csv(cand, cand_file, row.names = FALSE, quote = FALSE)

  outdir <- tempfile(pattern = "pcr_terap_reltime_")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    "Rscript",
    c(
      file.path(base_dir, "scripts", "run_pcr.R"),
      paste0("--ref-tree=", file.path(ex_dir, "Terapontoid_ML_MAIN_phylogram_used.tree")),
      paste0("--candidates-csv=", cand_file),
      paste0("--outdir=", outdir)
    )
  )
  if (!identical(status, 0L)) stop("Terapontoid PCR rerun failed.")

  summary_file <- file.path(ex_dir, "summary_terap_empirical_postfit_metrics.csv")
  file.copy(file.path(outdir, "summary_pcr_metrics.csv"), summary_file, overwrite = TRUE)
}

write_vgp_outputs <- function() {
  ref_tree_file <- file.path(vgp_dir, "roadies_v1.1.16b.numbers.nwk")
  cal_file <- file.path(vgp_dir, "roadies_manual_calibrations.csv")
  cand_file <- file.path(vgp_dir, "pcr_rerun_20260310", "candidates_core5.csv")

  ref_tree <- read.tree(ref_tree_file)
  cal_df <- read.csv(cal_file, stringsAsFactors = FALSE)
  rel <- run_reltime_with_bounds_ci(ref_tree, cal_df)

  rel_tree_tmp <- file.path(vgp_dir, "pcr_rerun_20260310", "roadies_v1.1.16b.numbers_RelTime_full_bounds.tre")
  rel_ci_tmp <- file.path(vgp_dir, "pcr_rerun_20260310", "roadies_v1.1.16b.numbers_RelTime_full_bounds_ci.csv")
  rel_bounds_tmp <- file.path(vgp_dir, "pcr_rerun_20260310", "roadies_v1.1.16b.numbers_RelTime_bounds_used.csv")
  write.tree(rel$tree, rel_tree_tmp)
  write.csv(rel$ci, rel_ci_tmp, row.names = FALSE)
  write.csv(rel$bounds, rel_bounds_tmp, row.names = FALSE)

  cand <- read.csv(cand_file, stringsAsFactors = FALSE)
  cand <- cand[cand$candidate != "RelTime", , drop = FALSE]
  cand <- rbind(cand, data.frame(candidate = "RelTime", tree_file = rel_tree_tmp, stringsAsFactors = FALSE))
  cand_tmp <- tempfile(pattern = "vgp_candidates_", fileext = ".csv")
  write.csv(cand, cand_tmp, row.names = FALSE, quote = FALSE)

  outdir <- tempfile(pattern = "pcr_vgp_reltime_")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    "Rscript",
    c(
      file.path(base_dir, "scripts", "run_pcr.R"),
      paste0("--ref-tree=", ref_tree_file),
      paste0("--candidates-csv=", cand_tmp),
      paste0("--calibrations-csv=", cal_file),
      paste0("--outdir=", outdir)
    )
  )
  if (!identical(status, 0L)) stop("Unpublished vertebrate PCR rerun failed.")

  summary_file <- file.path(
    base_dir, "examples", "unpublished_vertebrate", "postfit_metrics",
    "summary_unpublished_vertebrate_postfit_metrics.csv"
  )
  file.copy(file.path(outdir, "summary_pcr_metrics.csv"), summary_file, overwrite = TRUE)
}

write_terap_outputs()
write_vgp_outputs()

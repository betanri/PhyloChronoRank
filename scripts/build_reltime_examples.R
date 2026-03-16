#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else getwd()
if (!nzchar(script_dir) || !dir.exists(script_dir)) script_dir <- getwd()

message(
  "build_reltime_examples.R delegates to build_uncertainty_examples.R so ",
  "the bundled RelTime tree, bootstrap summary, supplemental Tao file, and ",
  "public example tables stay in sync."
)

status <- system2(
  "Rscript",
  c(file.path(script_dir, "build_uncertainty_examples.R"), args)
)
if (!identical(status, 0L)) {
  stop("Delegated uncertainty-example rebuild failed.")
}

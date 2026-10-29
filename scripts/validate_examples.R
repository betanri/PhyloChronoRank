## Check that the PCR tables in 2_PCR_POSTFIT_METRICS/README.md match their
## summary_pcr_metrics.csv files. Thin wrapper around readme_tables.R --check.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
status <- system2("Rscript", c(file.path(script_dir, "readme_tables.R"), "--check"))
quit(status = status)

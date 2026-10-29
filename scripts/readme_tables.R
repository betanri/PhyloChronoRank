## Build or check the PCR result tables in 2_PCR_POSTFIT_METRICS/README.md.
##
## Each table is generated from one summary_pcr_metrics.csv. The README marks
## it with a comment on the line above the table:
##   <!-- pcr-table: path/to/summary_pcr_metrics.csv -->
## (path relative to the repo root, or absolute).
##
##   Rscript scripts/readme_tables.R --check   compare every table cell with its CSV
##   Rscript scripts/readme_tables.R --write   regenerate every table from its CSV
##
## Tables whose CSV is not in the repo (unpublished datasets) are skipped.
args <- commandArgs(trailingOnly = TRUE)
mode <- if ("--write" %in% args) "write" else "check"
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
repo_dir <- normalizePath(file.path(script_dir, ".."))
readme <- file.path(repo_dir, "2_PCR_POSTFIT_METRICS", "README.md")
lines <- readLines(readme, warn = FALSE)

## README column -> CSV column, digits shown
cols <- list(
  c("burst loss", "burst_loss", 3),
  c("int. concordance", "internode_concordance", 3),
  c("compression", "compression_score", 3),
  c("depth R²", "depth_r2", 3),
  c("rate irregularity", "rate_irregularity", 3),
  c("mean relative gap", "mean_relative_gap", 3),
  c("uncertainty width", "uncertainty_mean_width_ma", 2),
  c("Fam 1", "family_radiation_score", 2),
  c("Fam 2", "family_global_score", 2),
  c("Fam 3", "family_ratecal_score", 2),
  c("core rank", "rank_mean_core", 2)
)
fmt <- function(x, d) ifelse(is.finite(x), formatC(round(x, d), format = "f", digits = d), "NA")

make_table <- function(df, unit_label) {
  keep <- vapply(cols, function(cc) cc[2] %in% names(df) && any(is.finite(df[[cc[2]]])), logical(1))
  use <- cols[keep]
  df <- df[order(df$rank_mean_core, df$candidate), , drop = FALSE]
  best <- min(df$rank_mean_core, na.rm = TRUE)
  head <- c("candidate", vapply(use, function(cc) if (cc[1] == "uncertainty width") paste0(cc[1], " (", unit_label, ")") else cc[1], ""))
  out <- c(paste0("| ", paste(head, collapse = " | "), " |"),
           paste0("| --- | ", paste(rep("---:", length(use)), collapse = " | "), " |"))
  for (i in seq_len(nrow(df))) {
    cells <- vapply(use, function(cc) fmt(df[[cc[2]]][i], as.integer(cc[3])), "")
    if (isTRUE(abs(df$rank_mean_core[i] - best) < 1e-9)) cells[length(cells)] <- paste0("**", cells[length(cells)], "**")
    out <- c(out, paste0("| `", df$candidate[i], "` | ", paste(cells, collapse = " | "), " |"))
  }
  out
}

markers <- grep("^<!-- pcr-table: .* -->$", lines)
if (!length(markers)) stop("No <!-- pcr-table: ... --> markers in ", readme)
problems <- 0L
checked <- 0L
new_lines <- lines
offset <- 0L
for (m in markers) {
  spec <- sub("^<!-- pcr-table: (.*) -->$", "\\1", lines[m])
  parts <- strsplit(spec, "\\s*\\|\\s*")[[1]]
  path <- parts[1]
  unit <- if (length(parts) > 1) sub("^unit=", "", parts[2]) else "Ma"
  csv <- if (grepl("^/", path)) path else file.path(repo_dir, path)
  start <- m + 1L
  end <- start
  while (end <= length(lines) && startsWith(trimws(lines[end]), "|")) end <- end + 1L
  end <- end - 1L
  if (!file.exists(csv)) {
    cat("skip  (CSV not in repo):", path, "\n")
    next
  }
  df <- read.csv(csv, stringsAsFactors = FALSE)
  tab <- make_table(df, unit)
  if (mode == "write") {
    a <- start + offset
    b <- end + offset
    new_lines <- c(new_lines[seq_len(a - 1L)], tab, if (b < length(new_lines)) new_lines[(b + 1L):length(new_lines)])
    offset <- offset + length(tab) - (end - start + 1L)
    cat("wrote", path, "\n")
  } else {
    checked <- checked + 1L
    current <- lines[start:end]
    if (!identical(trimws(current), trimws(tab))) {
      problems <- problems + 1L
      cat("FAIL ", path, "(README line", start, ")\n")
      d <- which(trimws(current) != trimws(tab)[seq_along(current)])
      for (k in head(d, 3L)) cat("   README:", current[k], "\n   CSV   :", tab[k], "\n")
    } else {
      cat("ok   ", path, "\n")
    }
  }
}
if (mode == "write") {
  writeLines(new_lines, readme)
} else {
  cat(checked, "table(s) checked,", problems, "mismatch(es)\n")
  if (problems) quit(status = 1)
}

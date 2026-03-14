args <- commandArgs(trailingOnly = TRUE)
quiet <- "--quiet" %in% args

msg <- function(...) {
  if (!quiet) cat(..., "\n", sep = "")
}

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

trim <- function(x) trimws(x, which = "both")

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else NA_character_
if (is.na(script_dir) || !nzchar(script_dir) || !dir.exists(script_dir)) script_dir <- getwd()
repo_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

readme_path <- file.path(repo_dir, "2_PCR_POSTFIT_METRICS", "README.md")
if (!file.exists(readme_path)) fail("2_PCR_POSTFIT_METRICS/README.md not found.")
readme_lines <- readLines(readme_path, warn = FALSE)

normalize_candidate <- function(x) {
  x <- gsub("`", "", x, fixed = TRUE)
  x <- trim(x)
  treepl_match <- regexec("^treePL \\(smooth *= *([^)]+)\\)$", x)
  treepl_parts <- regmatches(x, treepl_match)
  has_treepl <- lengths(treepl_parts) == 2L
  if (any(has_treepl)) {
    smooth_vals <- vapply(treepl_parts[has_treepl], `[`, character(1), 2L)
    x[has_treepl] <- paste0("treepl_best-smooth-", smooth_vals)
  }
  x
}

parse_md_table <- function(lines, start_idx) {
  table_lines <- character(0)
  for (i in start_idx:length(lines)) {
    ln <- lines[i]
    if (!startsWith(trim(ln), "|")) break
    table_lines <- c(table_lines, ln)
  }
  if (length(table_lines) < 3) fail("Markdown table parsing failed near line ", start_idx, ".")

  split_row <- function(line) {
    parts <- strsplit(line, "|", fixed = TRUE)[[1]]
    parts <- trim(parts)
    parts[parts != ""]
  }

  header <- split_row(table_lines[1])
  body <- table_lines[-c(1, 2)]
  rows <- lapply(body, split_row)
  max_cols <- length(header)
  rows <- lapply(rows, function(r) {
    if (length(r) != max_cols) fail("Unexpected markdown table width in README.")
    r
  })
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- header
  out[] <- lapply(out, function(col) trim(gsub("`", "", col, fixed = TRUE)))
  out
}

extract_table_after_heading <- function(example_heading, table_heading = "### Ranked post-fit results (lower is better)") {
  ex_idx <- grep(paste0("## ", example_heading), readme_lines, fixed = TRUE)
  if (!length(ex_idx)) fail("Could not find example heading: ", example_heading)
  tbl_head_idx <- grep(table_heading, readme_lines, fixed = TRUE)
  tbl_head_idx <- tbl_head_idx[tbl_head_idx > ex_idx[1]]
  if (!length(tbl_head_idx)) fail("Could not find table heading after ", example_heading)
  start_idx <- grep("^\\| candidate \\|", readme_lines)
  start_idx <- start_idx[start_idx > tbl_head_idx[1]]
  if (!length(start_idx)) fail("Could not find markdown table after ", example_heading)
  parse_md_table(readme_lines, start_idx[1])
}

fmt_num <- function(x, digits) sprintf(paste0("%.", digits, "f"), as.numeric(x))

assert_table_matches <- function(name, table_df, csv_df, mapping, digits, candidate_transform = identity) {
  table_df$candidate <- normalize_candidate(table_df$candidate)
  csv_df$candidate <- normalize_candidate(candidate_transform(csv_df$candidate))
  table_df <- table_df[order(table_df$candidate), , drop = FALSE]
  csv_df <- csv_df[order(csv_df$candidate), , drop = FALSE]

  if (!identical(table_df$candidate, csv_df$candidate)) {
    fail(name, ": README candidates do not match CSV candidates.\nREADME: ",
         paste(table_df$candidate, collapse = ", "),
         "\nCSV: ", paste(csv_df$candidate, collapse = ", "))
  }

  for (col_label in names(mapping)) {
    csv_col <- mapping[[col_label]]
    digs <- digits[[col_label]]
    left <- table_df[[col_label]]
    right <- fmt_num(csv_df[[csv_col]], digs)
    if (!identical(left, right)) {
      fail(name, ": mismatch in column '", col_label, "'.\nREADME: ",
           paste(left, collapse = ", "),
           "\nCSV: ", paste(right, collapse = ", "))
    }
  }
  msg("[ok] ", name, ": README table matches shipped CSV.")
}

compare_frames <- function(name, got, expected, cols, digits = 8) {
  got$candidate <- normalize_candidate(got$candidate)
  expected$candidate <- normalize_candidate(expected$candidate)
  got <- got[order(got$candidate), , drop = FALSE]
  expected <- expected[order(expected$candidate), , drop = FALSE]
  if (!identical(got$candidate, expected$candidate)) {
    fail(name, ": candidate mismatch.\nGot: ",
         paste(got$candidate, collapse = ", "),
         "\nExpected: ", paste(expected$candidate, collapse = ", "))
  }
  for (col in cols) {
    if (!all(col %in% names(got))) fail(name, ": missing rerun column ", col)
    if (!all(col %in% names(expected))) fail(name, ": missing expected column ", col)
    g <- round(as.numeric(got[[col]]), digits)
    e <- round(as.numeric(expected[[col]]), digits)
    if (!identical(g, e)) {
      fail(name, ": mismatch in rerun column '", col, "'.\nGot: ",
           paste(g, collapse = ", "),
           "\nExpected: ", paste(e, collapse = ", "))
    }
  }
  msg("[ok] ", name, ": rerun matches shipped metrics for shared columns.")
}

run_pcr <- function(args_vec, outdir) {
  unlink(outdir, recursive = TRUE, force = TRUE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  status <- system2("Rscript", c(file.path(repo_dir, "scripts", "run_pcr.R"), args_vec, paste0("--outdir=", outdir)))
  if (!identical(status, 0L)) fail("PCR rerun failed for outdir ", outdir)
  read.csv(file.path(outdir, "summary_pcr_metrics.csv"), stringsAsFactors = FALSE, check.names = FALSE)
}

# README vs shipped CSVs
table1 <- extract_table_after_heading("Example 1: Empirical dataset with two competing chronograms (Syngnatharia)")
csv1 <- read.csv(file.path(repo_dir, "examples", "syngnatharia", "postfit_metrics", "syngnatharia_postfit_metrics.csv"), stringsAsFactors = FALSE, check.names = FALSE)
assert_table_matches(
  "Example 1",
  table1,
  csv1,
  mapping = c(
    "burst loss" = "burst_loss",
    "pulse preservation (burst)" = "pulse_burst_selector_error",
    "pulse preservation (overall)" = "pulse_default_selector_error",
    "mean relative gap" = "mean_relative_gap",
    "rate irregularity" = "rate_irregularity",
    "uncertainty width (mean HPD width, Ma)" = "uncertainty_mean_width_ma",
    "core overall mean rank (pulse = 1/3)" = "rank_mean_3families"
  ),
  digits = c(
    "burst loss" = 4,
    "pulse preservation (burst)" = 4,
    "pulse preservation (overall)" = 4,
    "mean relative gap" = 4,
    "rate irregularity" = 4,
    "uncertainty width (mean HPD width, Ma)" = 2,
    "core overall mean rank (pulse = 1/3)" = 2
  )
)

table2 <- extract_table_after_heading("Example 2: Empirical dataset with six competing chronograms (Terapontoidei)")
csv2 <- read.csv(file.path(repo_dir, "examples", "terapontoid", "summary_terap_empirical_postfit_metrics.csv"), stringsAsFactors = FALSE, check.names = FALSE)
assert_table_matches(
  "Example 2",
  table2,
  csv2,
  mapping = c(
    "burst loss" = "burst_loss",
    "pulse preservation (burst)" = "pulse_burst_selector_error",
    "pulse preservation (overall)" = "pulse_default_selector_error",
    "rate irregularity" = "rate_irregularity",
    "overall mean rank (pulse = 1/2)" = "rank_mean_core"
  ),
  digits = c(
    "burst loss" = 4,
    "pulse preservation (burst)" = 4,
    "pulse preservation (overall)" = 4,
    "rate irregularity" = 4,
    "overall mean rank (pulse = 1/2)" = 2
  )
)

table3 <- extract_table_after_heading("Example 3: Unpublished vertebrate dataset (derived outputs only)")
csv3 <- read.csv(file.path(repo_dir, "examples", "unpublished_vertebrate", "postfit_metrics", "summary_unpublished_vertebrate_postfit_metrics.csv"), stringsAsFactors = FALSE, check.names = FALSE)
assert_table_matches(
  "Example 3",
  table3,
  csv3,
  mapping = c(
    "burst loss" = "burst_loss",
    "pulse preservation (burst)" = "pulse_burst_selector_error",
    "pulse preservation (overall)" = "pulse_default_selector_error",
    "mean relative gap" = "mean_relative_gap",
    "rate irregularity" = "rate_irregularity",
    "core overall mean rank (pulse = 1/3)" = "rank_mean_core"
  ),
  digits = c(
    "burst loss" = 4,
    "pulse preservation (burst)" = 4,
    "pulse preservation (overall)" = 4,
    "mean relative gap" = 4,
    "rate irregularity" = 4,
    "core overall mean rank (pulse = 1/3)" = 2
  ),
  candidate_transform = function(x) sub("^treepl_best-smooth-[^,]+$", "treePL", x)
)

# Rerun public examples
rerun1 <- run_pcr(
  c(
    paste0("--ref-tree=", file.path(repo_dir, "examples", "syngnatharia", "backbone_Raxml_besttree_matrix75.tre")),
    paste0("--candidates-csv=", file.path(repo_dir, "examples", "syngnatharia", "candidates.csv")),
    paste0("--calibrations-csv=", file.path(repo_dir, "examples", "syngnatharia", "calibrations_by_candidate.csv")),
    paste0("--uncertainty-csv=", file.path(repo_dir, "examples", "syngnatharia", "uncertainty_summary_long.csv"))
  ),
  outdir = file.path(tempdir(), "pcr_validate_syng")
)
compare_frames(
  "Example 1",
  rerun1,
  csv1,
  cols = c(
    "pulse_default_selector_error",
    "burst_loss",
    "pulse_burst_selector_error",
    "rate_irregularity",
    "mean_relative_gap",
    "uncertainty_mean_width_ma"
  )
)

rerun2 <- run_pcr(
  c(
    paste0("--ref-tree=", file.path(repo_dir, "examples", "terapontoid", "Terapontoid_ML_MAIN_phylogram_used.tree")),
    paste0("--candidates-csv=", file.path(repo_dir, "examples", "terapontoid", "candidates.csv"))
  ),
  outdir = file.path(tempdir(), "pcr_validate_terap")
)
compare_frames(
  "Example 2",
  rerun2,
  csv2,
  cols = c(
    "pulse_default_selector_error",
    "burst_loss",
    "pulse_burst_selector_error",
    "rate_irregularity",
    "rank_mean_core"
  )
)

msg("[ok] Example 3: README table matches derived CSV (rerun intentionally skipped; raw inputs are withheld).")
msg("")
msg("All bundled example validations passed.")

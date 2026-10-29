## Draw the Part 1 example figure (figures/terapontoid_pipeline_output.png):
## the input phylogram and every candidate chronogram of one staged run,
## with 95% node-age intervals where the tree carries them.
##
##   Rscript scripts/make_pipeline_example_figure.R [run_dir] [phylogram] [out_png]
## Defaults: examples/terapontoid/run_dense19, examples/terapontoid/phylogram.tre,
## figures/terapontoid_pipeline_output.png
suppressPackageStartupMessages({ library(ape); library(treeio) })
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
repo_dir <- normalizePath(file.path(script_dir, ".."))
a <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(a) >= 1) a[1] else file.path(repo_dir, "examples/terapontoid/run_dense19")
phy_file <- if (length(a) >= 2) a[2] else file.path(repo_dir, "examples/terapontoid/phylogram.tre")
out_png <- if (length(a) >= 3) a[3] else file.path(repo_dir, "figures/terapontoid_pipeline_output.png")
main_dir <- file.path(run_dir, "MAIN_OUTPUT_TREES")

panels <- list(
  c("^chronos_clock", "chronos clock"),
  c("^chronos_correlated", "chronos correlated"),
  c("^chronos_discrete", "chronos discrete"),
  c("^chronos_relaxed", "chronos relaxed"),
  c("^treePL", "treePL"),
  c("^RelTime_MEGA_with_tao_CI", "RelTime MEGA (CIs: Tao et al.)"),
  c("^RelTime_MEGA_with_bootstrap_CI", "RelTime MEGA (CIs: Poisson bootstrap)"),
  c("^MDCat", "MD-Cat (no CIs)")
)
files <- list.files(main_dir, pattern = "\\.tre$", full.names = TRUE)

## Read a candidate with its 95% intervals (NEXUS with height_95%_HPD) or as
## plain Newick. Returns the tree and a node x interval matrix (ages).
read_candidate <- function(f) {
  first <- readLines(f, n = 1L, warn = FALSE)
  if (grepl("^#NEXUS", first)) {
    b <- treeio::read.beast(f)
    tr <- treeio::as.phylo(b)
    d <- treeio::get.data(b)
    hpd <- d[[intersect(c("height_95%_HPD", "height_0.95_HPD"), names(d))[1]]]
    ci <- matrix(NA_real_, ape::Ntip(tr) + tr$Nnode, 2L)
    if (!is.null(hpd)) {
      for (k in seq_len(nrow(d))) {
        v <- suppressWarnings(as.numeric(hpd[[k]]))
        if (length(v) == 2L && all(is.finite(v))) ci[as.integer(d$node[k]), ] <- sort(v)
      }
    }
  } else {
    tr <- ape::read.tree(f)
    ci <- matrix(NA_real_, ape::Ntip(tr) + tr$Nnode, 2L)
  }
  list(tree = tr, ci = ci)
}

phy <- ape::ladderize(ape::read.tree(phy_file), right = FALSE)
used <- lapply(panels, function(p) { hit <- grep(p[1], basename(files)); if (length(hit)) files[hit[1]] else NA })
keep <- !is.na(unlist(used))
panels <- panels[keep]
used <- unlist(used)[keep]
ages <- vapply(used, function(f) { t <- read_candidate(f)$tree; max(ape::node.depth.edgelength(t)) }, numeric(1))
xmax <- max(ages) * 1.02

n_pan <- length(panels) + 1L
ncol <- if (n_pan <= 8L) 4L else 3L
nrow <- ceiling(n_pan / ncol)
png(out_png, width = 800 * ncol, height = 800 * nrow, res = 150)
par(mfrow = c(nrow, ncol), mar = c(4, 1, 2.5, 1))
plot(phy, show.tip.label = FALSE, main = "Input phylogram", cex.main = 1.3)
add.scale.bar(cex = 0.8)
for (i in seq_along(panels)) {
  cand <- read_candidate(used[i])
  tr <- ape::ladderize(cand$tree, right = FALSE)
  ## ladderize keeps node numbers, so the interval rows still match
  root_age <- max(ape::node.depth.edgelength(tr))
  plot(tr, show.tip.label = FALSE, main = panels[[i]][2], cex.main = 1.3,
       x.lim = c(root_age - xmax, root_age))
  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  nodes <- (ape::Ntip(tr) + 1L):(ape::Ntip(tr) + tr$Nnode)
  for (nd in nodes) {
    if (all(is.finite(cand$ci[nd, ]))) {
      rect(root_age - cand$ci[nd, 2], lp$yy[nd] - 0.35, root_age - cand$ci[nd, 1], lp$yy[nd] + 0.35,
           col = grDevices::adjustcolor("#3f2fb5", 0.55), border = NA)
    }
  }
  ticks <- pretty(c(0, xmax))
  ticks <- ticks[ticks <= xmax]
  axis(1, at = root_age - ticks, labels = ticks)
  mtext("Ma", side = 1, line = 2.5, cex = 0.8, font = 2)
}
dev.off()
message("Wrote ", out_png)

suppressPackageStartupMessages({ library(ape) })

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
if (length(file_arg)) {
  script_dir <- dirname(normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE))
} else {
  script_dir <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
}
base_dir <- normalizePath(file.path(script_dir, '..'), winslash = '/', mustWork = TRUE)
out_fig <- file.path(base_dir, "figures")
out_csv <- normalizePath(file.path(base_dir, "..", "2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE", "EXAMPLE_FILES", "OUTPUT_DEMO"),
                         winslash = "/", mustWork = TRUE)

tempo <- read.csv(file.path(out_csv, "summary_terap_empirical_model_fits.csv"), stringsAsFactors = FALSE)
tempo <- tempo[order(tempo$tempo_composite, tempo$tempo_mae_all), ]
best2 <- tempo$model[1:2]
worst2 <- tempo$model[(nrow(tempo)-1):nrow(tempo)]

read_tree <- function(model) {
  p <- file.path(out_csv, sprintf("Terapontoid_ML_MAIN_chronos_dated_model%s.tre", model))
  read.tree(p)
}

phy <- read.tree(file.path(out_csv, "Terapontoid_ML_MAIN_phylogram_used.tree"))
mods <- list(
  discrete = read_tree("discrete"),
  clock = read_tree("clock"),
  correlated = read_tree("correlated"),
  relaxed = read_tree("relaxed")
)

normalize_tree <- function(tr) {
  d <- node.depth.edgelength(tr)
  m <- max(d[1:Ntip(tr)])
  tr$edge.length <- tr$edge.length / m
  tr
}

get_norm_height <- function(tr) {
  d <- node.depth.edgelength(tr)
  m <- max(d[1:Ntip(tr)])
  (m - d) / m
}

# 1) Clean tree panel + arrows on top row only
p_phy <- normalize_tree(phy)
p_mods <- lapply(mods, normalize_tree)

order_models <- c(best2, worst2)
order_models <- unique(order_models)
if (!all(c("discrete", "clock", "correlated", "relaxed") %in% order_models)) {
  order_models <- c("discrete", "clock", "correlated", "relaxed")
}

clade_key <- function(tips) paste(sort(tips), collapse = "|")
internal_node_keys <- function(tr) {
  n <- Ntip(tr)
  nodes <- (n + 1):(n + tr$Nnode)
  keys <- sapply(nodes, function(nd) clade_key(extract.clade(tr, nd)$tip.label))
  names(keys) <- as.character(nodes)
  keys
}
find_nodes_by_keys <- function(tr, keys) {
  kmap <- internal_node_keys(tr)
  idx <- match(keys, unname(kmap))
  as.integer(names(kmap)[idx])
}
pick_focal_nodes <- function(tr) {
  n <- Ntip(tr)
  h <- get_norm_height(tr)
  nodes <- (n + 1):(n + tr$Nnode)
  early <- nodes[h[nodes] >= quantile(h[nodes], 0.75, na.rm = TRUE)]
  tips_list <- lapply(early, function(nd) sort(extract.clade(tr, nd)$tip.label))
  sizes <- sapply(tips_list, length)
  ord <- order(sizes, decreasing = TRUE)
  selected <- integer(0)
  selected_tips <- list()
  for (k in ord) {
    cand_nd <- early[k]
    cand_tips <- tips_list[[k]]
    nested <- FALSE
    if (length(selected_tips) > 0) {
      for (st in selected_tips) {
        if (all(cand_tips %in% st) || all(st %in% cand_tips)) {
          nested <- TRUE
          break
        }
      }
    }
    if (!nested) {
      selected <- c(selected, cand_nd)
      selected_tips[[length(selected_tips) + 1]] <- cand_tips
    }
    if (length(selected) == 2) break
  }
  if (length(selected) < 2) selected <- early[ord[1:2]]
  selected
}

focal_phy_nodes <- pick_focal_nodes(p_phy)
focal_keys <- sapply(focal_phy_nodes, function(nd) clade_key(extract.clade(p_phy, nd)$tip.label))

draw_arrows <- function(tr, keyset) {
  n <- Ntip(tr)
  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  nodes <- find_nodes_by_keys(tr, keyset)
  nodes <- nodes[is.finite(nodes)]
  if (!length(nodes)) return(invisible(NULL))
  dy <- 0.08 * diff(lp$y.lim)
  for (i in seq_along(nodes)) {
    nd <- nodes[i]
    xt <- lp$xx[nd]
    yt <- lp$yy[nd - n]
    x0 <- max(lp$x.lim[1] + 0.05, xt - 0.11)
    y0 <- yt + ifelse(i == 1, -dy, dy)
    arrows(x0, y0, xt - 0.004, yt, length = 0.09, col = "#e31a1c", lwd = 2.7)
  }
}

draw_tree_panel <- function() {
  layout(matrix(c(1,2,3,4,5,6), nrow = 2, byrow = TRUE))
  par(mar = c(0.9, 0.9, 3.2, 0.9), oma = c(0.3, 0.3, 2.4, 0.3), xpd = NA)

  plot(p_phy, show.tip.label = FALSE, no.margin = TRUE, main = "")
  title(main = "Input phylogram", cex.main = 2.1, line = 0.95, font.main = 2)
  draw_arrows(p_phy, focal_keys)

  for (m in order_models) {
    rec <- tempo[tempo$model == m, ][1, ]
    ttl <- sprintf("%s | all %.3f | early %.3f", tools::toTitleCase(m), rec$tempo_mae_all, rec$tempo_mae_early_q75)
    plot(p_mods[[m]], show.tip.label = FALSE, no.margin = TRUE, main = "")
    title(main = ttl, cex.main = 2.1, line = 0.95, font.main = 2)
    if (m %in% c("discrete", "clock")) draw_arrows(p_mods[[m]], focal_keys)
  }

  plot.new()
  legend("center",
         legend = c("Best two by composite: discrete, clock",
                    "Worst two by composite: relaxed, correlated",
                    "Manual annotation ready (no auto arrows/dots)."),
         bty = "n", cex = 1.75, text.font = 2)

  mtext("Figure A: Tree-shape comparison (best vs worst)",
        side = 3, outer = TRUE, line = 0.25, cex = 2.8, font = 2)
}

png(file.path(out_fig, "branching_tempo_tree_panel_clean_v3.png"), width = 3000, height = 2100, res = 220)
draw_tree_panel()
dev.off()

pdf(file.path(out_fig, "branching_tempo_tree_panel_clean_v3.pdf"), width = 15, height = 10.5, onefile = FALSE)
draw_tree_panel()
dev.off()

# 2) Node-height ECDF vs phylogram
get_internal_h <- function(tr) {
  h <- get_norm_height(tr)
  n <- Ntip(tr)
  h[(n + 1):(n + tr$Nnode)]
}

h_phy <- get_internal_h(phy)
hs <- lapply(mods, get_internal_h)

png(file.path(out_fig, "branching_tempo_node_height_ecdf.png"), width = 1400, height = 950, res = 180)
par(mar = c(5,5,3,1))
plot(ecdf(h_phy), do.points = FALSE, verticals = FALSE, lwd = 3, col = "black",
     xlab = "Normalized node height (0 tip-side, 1 root-side)", ylab = "ECDF",
     main = "Branching-tempo shape vs input phylogram")
cols <- c(discrete = "#1b9e77", clock = "#2c7fb8", correlated = "#d95f0e", relaxed = "#7570b3")
for (m in names(hs)) {
  lines(ecdf(hs[[m]]), do.points = FALSE, verticals = FALSE, lwd = 2.5, col = cols[m], lty = ifelse(m %in% best2, 1, 2))
}
legend("bottomright",
       legend = c("phylogram", "discrete", "clock", "correlated", "relaxed"),
       col = c("black", cols), lwd = c(3,2.5,2.5,2.5,2.5), lty = c(1,1,1,2,2), bty = "n")
dev.off()

# 3) Tempo metrics bar chart
show <- tempo[match(c("discrete","clock","correlated","relaxed"), tempo$model), ]

png(file.path(out_fig, "branching_tempo_metric_bars.png"), width = 1400, height = 950, res = 180)
par(mar = c(8,5,4,1))
mat <- rbind(show$tempo_mae_all, show$tempo_mae_early_q75)
colnames(mat) <- show$model
ylim_top <- max(mat, na.rm = TRUE) * 1.22
barplot(mat, beside = TRUE, col = c("#6baed6", "#fb6a4a"), las = 2,
        ylab = "Tempo MAE (lower is better)", main = "Branching-tempo metric by model",
        ylim = c(0, ylim_top))
legend("top", inset = c(0, -0.02), horiz = TRUE,
       legend = c("overall MAE", "early_q75 MAE"),
       fill = c("#6baed6", "#fb6a4a"), bty = "n", cex = 1.1)
dev.off()

write.csv(tempo, file.path(out_csv, "summary_terap_empirical_model_fits_ranked.csv"), row.names = FALSE)

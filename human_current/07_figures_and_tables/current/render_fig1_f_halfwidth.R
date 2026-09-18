if (!nzchar(Sys.getenv("NMF_WORK_ROOT"))) stop("Set NMF_WORK_ROOT; original source directories must remain read-only.")
#!/usr/bin/env Rscript
options(future.globals.maxSize = 50 * 1024^3)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggraph)
  library(igraph)
  library(tidygraph)
  library(dplyr)
  library(scales)
  library(grid)
  library(svglite)
})

# Native half-width redraw of the confirmed Fig1 f network.  The source tables
# are read-only archived outputs; no correlation matrix or network statistic is
# recomputed upstream.  Edge/node rules match the GITHUB v15 panel-j producer.
PROJ <- paste0(Sys.getenv("NMF_WORK_ROOT"), "/inputs/cortex_nmf_program/archived")
RES  <- file.path(PROJ, "results/crossregion_v1")
INT  <- file.path(PROJ, "figures/fig1/_intermediate")
OUT  <- paste0(Sys.getenv("NMF_WORK_ROOT"), "/figures/human_revision/panels/Fig1_f_halfwidth.pdf")

CLASS_COL <- c(
  excitatory = "#2166AC",
  inhibitory = "#B2182B",
  `non-neuronal` = "#1B7837"
)

# ---- canonical retained-54 functional names from Table S3
pn <- read.delim(
  paste0(Sys.getenv("NMF_WORK_ROOT"), "/tables/TableS3_program_annotation.tsv"),
  quote = "", comment.char = "", stringsAsFactors = FALSE,
  check.names = FALSE
)
pn$program <- as.integer(pn$cnmf_component)
pn$lab_short <- as.character(pn$functional_name)
short_map <- setNames(pn$lab_short, pn$program)
.lsh <- function(p) {
  key <- as.character(p)
  ifelse(key %in% names(short_map), short_map[key], paste0("Component ", p))
}

spec <- read.csv(file.path(INT, "program_specificity.csv"), stringsAsFactors = FALSE)
spec$program <- as.integer(spec$program)

# ---- original panel-j graph rules (do not change)
DROP <- c(4L, 18L, 35L, 52L)
cm <- read.csv(
  file.path(INT, "program_program_corr.csv"),
  check.names = FALSE, row.names = 1
)
cm <- as.matrix(cm)
diag(cm) <- 0
progs_all <- as.integer(colnames(cm))
keep <- !(progs_all %in% DROP)
cm <- cm[keep, keep, drop = FALSE]
progs_j <- as.integer(colnames(cm))

nn <- nrow(cm)
KTOP <- 3L
MINW <- 0.10
ekeep <- matrix(FALSE, nn, nn)
for (i in seq_len(nn)) {
  ord <- order(cm[i, ], decreasing = TRUE)
  ord <- ord[ord != i][seq_len(KTOP)]
  ord <- ord[cm[i, ord] > MINW]
  ekeep[i, ord] <- TRUE
}
ekeep <- ekeep | t(ekeep)
idx <- which(upper.tri(ekeep) & ekeep, arr.ind = TRUE)
edges <- data.frame(
  from = as.character(progs_j[idx[, 1]]),
  to = as.character(progs_j[idx[, 2]]),
  weight = cm[idx],
  stringsAsFactors = FALSE
)
nodes <- data.frame(program = progs_j) %>%
  left_join(spec[, c("program", "dominant_class", "gini")], by = "program")
nodes$name <- as.character(nodes$program)
g <- tbl_graph(
  nodes = nodes, edges = edges, directed = FALSE, node_key = "name"
) %>%
  activate(nodes) %>%
  mutate(deg = centrality_degree(weights = NULL))
nd <- g %>% activate(nodes) %>% as_tibble()
N_HUB <- 12L
deg_cut <- sort(nd$deg, decreasing = TRUE)[min(N_HUB, length(nd$deg))]
is_hub <- nd$deg >= deg_cut & nd$deg > 0
sl <- .lsh(nd$program)
wrap_label <- function(x) {
  if (!nzchar(x)) return(x)
  paste(strwrap(x, width = 18), collapse = "\n")
}
wrapped_sl <- vapply(sl, wrap_label, character(1))
g <- g %>% activate(nodes) %>%
  mutate(
    hub_lab = ifelse(is_hub & nzchar(sl), wrapped_sl, NA_character_),
    is_hub = is_hub
  )

# ---- deterministic original FR layout, rescaled for the new near-square half
set.seed(7)
lay <- create_layout(g, layout = "fr", niter = 2000)
lay$x <- scales::rescale(lay$x, to = c(0, 0.96))
lay$y <- scales::rescale(lay$y, to = c(0, 1.0))

p <- ggraph(lay) +
  geom_edge_link(aes(width = weight, alpha = weight), colour = "grey60") +
  geom_node_point(
    aes(fill = dominant_class, size = deg + 1),
    shape = 21, colour = "white", stroke = 0.3
  ) +
  ggrepel::geom_text_repel(
    data = subset(lay, !is.na(hub_lab)),
    aes(x = x, y = y, label = hub_lab),
    size = 1.85, fontface = "bold", colour = "grey10",
    bg.color = "white", bg.r = 0.10,
    box.padding = 0.22, point.padding = 0.10,
    min.segment.length = 0, segment.size = 0.18,
    segment.colour = "grey55", force = 0.70,
    force_pull = 1.1, max.overlaps = Inf, max.iter = 4000,
    seed = 7, na.rm = TRUE
  ) +
  scale_edge_width(range = c(0.12, 0.9), guide = "none") +
  scale_edge_alpha(range = c(0.12, 0.6), guide = "none") +
  scale_fill_manual(values = CLASS_COL, name = "dominant class") +
  scale_size_continuous(range = c(1.2, 5), name = "degree") +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  scale_y_continuous(expand = expansion(mult = 0.08)) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(title = "Co-activity network") +
  theme_void(base_size = 7, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.text = element_text(size = 4.3),
    legend.title = element_text(size = 4.8),
    legend.key.size = unit(1.8, "mm"),
    legend.key.width = unit(4.5, "mm"),
    legend.spacing.x = unit(0.6, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(2.5, 2.5, 1.0, 2.5)
  )

ggsave(
  OUT, p,
  width = 175 / 72, height = 196 / 72, units = "in",
  device = cairo_pdf, bg = "white", limitsize = FALSE
)

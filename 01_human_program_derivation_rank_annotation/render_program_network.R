PROJECT_ROOT <- Sys.getenv("CORTEX_PROGRAM_ROOT", unset = "/DATA/cortex_nmf_program")
suppressPackageStartupMessages({
    library(ggplot2)
    library(ComplexHeatmap)
    library(circlize)
    library(ggrepel)
    library(ggridges)
    library(ggraph)
    library(igraph)
    library(tidygraph)
    library(dplyr)
    library(tidyr)
    library(grid)
    library(scales)
    library(RColorBrewer)
    library(svglite)
    library(ggrastr)
    library(png)
})
set.seed(42)
PROJ <- PROJECT_ROOT
RES <- file.path(PROJ, "results/crossregion_v1")
INT <- file.path(PROJ, "figures/fig1/_intermediate")
OUT <- file.path(PROJ, "figures/fig1")
SVGD <- file.path(OUT, "svg_panels")
BOX <- list(j = c(w = 99.077, h = 61.923))
mm2in <- function(x) x/25.4
CLASS_COL <- c(excitatory = "#2166AC", inhibitory = "#B2182B", `non-neuronal` = "#1B7837")
pn <- read.delim(file.path(RES, "program_names.tsv"), quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
pn$program <- as.integer(pn$cnmf_component)
retained_raw <- pn$program[pn$new_P != "EXCLUDED"]
pn$fdr <- suppressWarnings(as.numeric(pn$fdr))
pn$lab_short <- ifelse(!is.na(pn$fdr) & pn$fdr >= 0.25, paste0(pn$name_short, "*"), pn$name_short)
short_map <- setNames(pn$lab_short, pn$program)
.lsh <- function(p) ifelse(as.character(p) %in% names(short_map), short_map[as.character(p)], paste0("P", p))
scc <- read.csv(file.path(INT, "subclass_class.csv"), stringsAsFactors = FALSE)
class_map <- setNames(scc$class, scc$subclass)
spec <- read.csv(file.path(INT, "program_specificity.csv"))
spec$program <- as.integer(spec$program)
save_gg_svg <- function(p, panel, scale_in = 1) {
    b <- BOX[[panel]]
    w <- mm2in(b["w"]) * scale_in
    h <- mm2in(b["h"]) * scale_in
    f <- file.path(SVGD, sprintf("fig1_%s.svg", panel))
    svglite(f, width = as.numeric(w), height = as.numeric(h), bg = "white")
    plot(p)
    dev.off()
}
theme_nat <- function(base = 6.6) {
    theme_classic(base_size = base, base_family = "Helvetica") + theme(panel.grid = element_blank(), plot.background = element_rect(fill = "white", colour = NA))
}
{
    cm <- read.csv(file.path(INT, "program_program_corr.csv"), check.names = FALSE, row.names = 1)
    cm <- as.matrix(cm)
    diag(cm) <- 0
    progs_all <- as.integer(colnames(cm))
    cm <- cm[as.character(retained_raw), as.character(retained_raw), drop = FALSE]
    progs_j <- as.integer(colnames(cm))
    nn <- nrow(cm)
    KTOP <- 3L
    MINW <- 0.1
    ekeep <- matrix(FALSE, nn, nn)
    for (i in seq_len(nn)) {
        ord <- order(cm[i, ], decreasing = TRUE)
        ord <- ord[ord != i][seq_len(KTOP)]
        ord <- ord[cm[i, ord] > MINW]
        ekeep[i, ord] <- TRUE
    }
    ekeep <- ekeep | t(ekeep)
    idx <- which(upper.tri(ekeep) & ekeep, arr.ind = TRUE)
    edges <- data.frame(from = as.character(progs_j[idx[, 1]]), to = as.character(progs_j[idx[, 2]]), weight = cm[idx], stringsAsFactors = FALSE)
    nodes <- data.frame(program = progs_j) %>% left_join(spec[, c("program", "dominant_class", "gini")], by = "program")
    nodes$name <- as.character(nodes$program)
    g <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE, node_key = "name") %>% activate(nodes) %>% mutate(deg = centrality_degree(weights = NULL))
    nd <- g %>% activate(nodes) %>% as_tibble()
    N_HUB <- 12L
    deg_cut <- sort(nd$deg, decreasing = TRUE)[min(N_HUB, length(nd$deg))]
    is_hub <- nd$deg >= deg_cut & nd$deg > 0
    sl <- .lsh(nd$program)
    g <- g %>% activate(nodes) %>% mutate(hub_lab = ifelse(is_hub & nzchar(sl), sl, NA_character_), is_hub = is_hub)
    set.seed(7)
    FF_AR <- as.numeric(BOX[["j"]]["w"]/BOX[["j"]]["h"])
    lay <- create_layout(g, layout = "fr", niter = 2000)
    lay$x <- scales::rescale(lay$x, to = c(0, FF_AR))
    lay$y <- scales::rescale(lay$y, to = c(0, 1))
    p_j <- ggraph(lay) + geom_edge_link(aes(width = weight, alpha = weight), colour = "grey60") + geom_node_point(aes(fill = dominant_class, size = deg + 1), shape = 21, colour = "white", stroke = 0.3) + ggrepel::geom_text_repel(data = subset(lay, !is.na(hub_lab)), aes(x = x, y = y, label = hub_lab), size = 1.95, fontface = "bold", colour = "grey10", bg.color = "white", bg.r = 0.12, box.padding = 0.18, point.padding = 0.1, min.segment.length = 0, segment.size = 0.2, segment.colour = "grey55", force = 0.4, 
        force_pull = 1.2, max.overlaps = 20, max.iter = 4000, seed = 7, na.rm = TRUE) + scale_edge_width(range = c(0.12, 0.9), guide = "none") + scale_edge_alpha(range = c(0.12, 0.6), guide = "none") + scale_fill_manual(values = CLASS_COL, name = "dominant class") + scale_size_continuous(range = c(1.2, 5), name = "degree") + scale_x_continuous(expand = expansion(mult = 0.06)) + scale_y_continuous(expand = expansion(mult = 0.08)) + coord_fixed(ratio = 1, clip = "off") + labs(title = "Co-activity network") + 
        theme_void(base_size = 7, base_family = "Helvetica") + theme(plot.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 5.5), legend.title = element_text(size = 6), legend.key.size = unit(3.2, "mm"), legend.position = "right", plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(0, 0, 0, 0))
    save_gg_svg(p_j, "j")
}

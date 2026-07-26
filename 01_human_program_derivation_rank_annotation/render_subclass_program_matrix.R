.annotation_group <- function(value) ifelse(endsWith(as.character(value), "sig"), "class_a", "class_b")
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
dir.create(SVGD, recursive = TRUE, showWarnings = FALSE)
BOX <- list(a = c(w = 173, h = 40), b = c(w = 105.861, h = 102.553), c = c(w = 62, h = 62), d = c(w = 220, h = 52), e = c(w = 64, h = 92), f = c(w = 90, h = 90), g = c(w = 78, h = 46), h = c(w = 55.139, h = 42.414), i = c(w = 61.923, h = 61.923), j = c(w = 99.077, h = 61.923))
mm2in <- function(x) x/25.4
CLASS_COL <- c(excitatory = "#3B6FB6", inhibitory = "#C0392B", `non-neuronal` = "#27865E")
CLASS_LV <- c("excitatory", "inhibitory", "non-neuronal")
theme_nat <- function(base = 6.6) {
    theme_classic(base_size = base, base_family = "Helvetica") + theme(axis.line = element_line(linewidth = 0.35, colour = "black"), axis.ticks = element_line(linewidth = 0.35, colour = "black"), axis.title = element_text(size = base), axis.text = element_text(size = base - 0.6, colour = "black"), legend.title = element_text(size = base - 0.4), legend.text = element_text(size = base - 1), legend.key.size = unit(3.2, "mm"), strip.text = element_text(size = base - 0.2, face = "bold"), strip.background = element_blank(), 
        plot.title = element_text(size = base + 1, face = "bold"), plot.subtitle = element_text(size = base - 0.6, colour = "grey30"), panel.grid = element_blank(), plot.background = element_rect(fill = "white", colour = NA))
}
theme_set(theme_nat())
pn <- read.delim(file.path(RES, "program_names.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "", check.names = FALSE)
pn$program <- as.integer(pn$cnmf_component)
retained_raw <- pn$program[pn$new_P != "EXCLUDED"]
pn$fdr <- suppressWarnings(as.numeric(pn$fdr))
pn$lab_short <- ifelse(.annotation_group(pn$confidence) == "class_b", paste0(pn$name_short, "*"), pn$name_short)
lab_map <- setNames(pn$name_full, pn$program)
short_map <- setNames(pn$lab_short, pn$program)
.lab <- function(p) ifelse(as.character(p) %in% names(lab_map), lab_map[as.character(p)], paste0("P", p))
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
{
    m <- read.csv(file.path(INT, "program_x_subclass_mean.csv"), check.names = FALSE)
    progs_all <- as.integer(m$program)
    mat <- as.matrix(m[, -1])
    rownames(mat) <- progs_all
    mat <- mat[progs_all %in% retained_raw, , drop = FALSE]
    mat <- mat[match(retained_raw, rownames(mat)), , drop = FALSE]
    n_prog_f <- nrow(mat)
    frac <- mat/rowSums(mat)
    subs <- colnames(frac)
    S <- matrix(0, length(subs), length(subs), dimnames = list(subs, subs))
    for (i in seq_len(nrow(frac))) {
        w <- frac[i, ]
        S <- S + outer(w, w)
    }
    diag(S) <- 0
    sub_cls <- class_map[subs]
    ord <- order(match(sub_cls, CLASS_LV), subs)
    subs_o <- subs[ord]
    cls_o <- sub_cls[ord]
    S <- S[subs_o, subs_o]
    offdiag <- S[upper.tri(S)]
    thr <- as.numeric(quantile(offdiag, 0.8))
    Sthr <- S
    Sthr[Sthr < thr] <- 0
    grid.col <- setNames(unname(CLASS_COL[cls_o]), subs_o)
    col_mat <- matrix(NA_character_, nrow(Sthr), ncol(Sthr), dimnames = dimnames(Sthr))
    for (i in seq_len(nrow(Sthr))) {
        cc <- grDevices::adjustcolor(unname(CLASS_COL[cls_o[i]]), alpha.f = 0.45)
        col_mat[i, ] <- cc
    }
    col_mat[Sthr == 0] <- "#00000000"
    bf <- BOX[["f"]]
    fsvg <- file.path(SVGD, "fig1_f.svg")
    svglite(fsvg, width = mm2in(bf["w"]), height = mm2in(bf["h"]), bg = "white")
    par(mar = c(0, 0, 0, 0), xpd = NA)
    circlize::circos.clear()
    cls_run <- rle(as.character(cls_o))$lengths
    gvec <- unlist(lapply(seq_along(cls_run), function(k) c(rep(1.6, cls_run[k] - 1), 7)))
    names(gvec) <- subs_o
    circlize::circos.par(gap.after = gvec, start.degree = 90, points.overflow.warning = FALSE, cell.padding = c(0, 0, 0, 0), canvas.xlim = c(-1.45, 1.45), canvas.ylim = c(-1.45, 1.45))
    circlize::chordDiagram(Sthr, order = subs_o, grid.col = grid.col, col = col_mat, transparency = 0.5, directional = 0, symmetric = TRUE, reduce = 0, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.18), link.lwd = 0.3, link.border = NA, scale = FALSE)
    circlize::circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
        s <- circlize::get.cell.meta.data("sector.index")
        xlm <- circlize::get.cell.meta.data("xlim")
        circlize::circos.text(mean(xlm), circlize::get.cell.meta.data("ylim")[2] + 0.4, labels = s, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.46, col = "black")
    })
    circlize::circos.clear()
    text(0, 1.4, "Program sharing", cex = 0.8, font = 2, adj = c(0.5, 1))
    ly <- -1.36
    lx <- -0.66
    for (k in seq_along(CLASS_LV)) {
        points(lx + (k - 1) * 0.66, ly, pch = 22, cex = 0.9, bg = unname(CLASS_COL[CLASS_LV[k]]), col = "grey30", lwd = 0.4)
        text(lx + (k - 1) * 0.66 + 0.05, ly, CLASS_LV[k], adj = c(0, 0.5), cex = 0.46)
    }
    dev.off()
    n_links <- sum(Sthr[upper.tri(Sthr)] > 0)
}

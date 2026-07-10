#!/usr/bin/env Rscript
# =============================================================================
# Figure B (LARGE / INTEGRATED) — Program x Program spatial co-organization
# ONE CNS-grade full-page figure = abstract 8-panel logic + spatial tissue
# exemplars ported to NATIVE ggplot patchwork (rule 8: same object, multi-modality;
# unified theme/fonts; vector cairo_pdf). NO PNG montage.
#
# Framing: HONEST. co-org = real biology; gene-similarity = honest annotation
# (NOT dismissal). bin50 = 25um, NOT dense mixing -> NO co-density-floor argument.
#
# Order: PAN-REGIONAL commonality first
#   a  60x60 co-org heatmap (r=25um)
#   b  co-org vs gene-cosine scatter
#   SP1 spatial grounding grid: 3 chips x [Myelination P45 | Oligodendrocyte dev. P13 | Neurofilament cytoskeleton P8]
#       (2 oligodendroglial programs co-distribute; neurofilament-cytoskeleton P8 avoids)
#   c  low-gene-overlap niche network (microglia-synapse hub P36/P45/P54)
#   SP2 the SAME hub IN TISSUE: P36 micro map | P26 glutamate-receptor map | co-distribution r~0.79
#   d  strongest avoidances lollipop
#   e  rigor (reproducibility + donor-LOO)
# then REGIONAL difference
#   f  region small-multiples | g region clustering | h cross-region span
# =============================================================================
suppressPackageStartupMessages({
  library(reticulate); library(ggplot2); library(patchwork)
  library(ComplexHeatmap); library(circlize); library(grid)
  library(igraph); library(ggraph); library(ggridges)
  library(ggrepel); library(dendextend); library(data.table)
  library(viridisLite); library(scales)
  library(svglite); library(ggrastr)
})
# RELAYOUT: dense per-bin spatial point layers (sp1 9 maps + sp2 maps + sp2
# scatter) -> 100k-1M square markers blow past rsvg's XML element cap during
# ink-crop. Rasterise JUST the point layers (400dpi); axes/text/legends stay
# vector. Scoped helper for the spatial geoms only (abstract panels stay vector).
geom_point_r <- function(...) ggrastr::rasterise(ggplot2::geom_point(...), dpi=400, dev="ragg")
np <- import("numpy")

ROOT  <- "CORTEX_PROGRAM_ROOT"
MC    <- file.path(ROOT, "results/crossregion_v1/markcorr_v2")
FIN   <- file.path(MC, "final"); SIM <- file.path(MC, "similarity"); RIG <- file.path(MC, "rigor")
SPR   <- file.path(MC, "spatial_for_R")
# BETWEENCHIP NULL UPDATE (2026-06-25):
# Load between-chip Stouffer q for progprog (1208 headline pairs).
BTWNDIR <- file.path(ROOT, "results/crossregion_v1/markcorr_betweenchip_v1")
sq_pp <- read.delim(file.path(BTWNDIR, "betweenchip_progprog_stouffer_q.tsv"),
                    stringsAsFactors=FALSE)
sq_pp_hl <- sq_pp[sq_pp$is_headline, ]  # 1208 headline pairs
sq_pp_hl$pa <- as.integer(sub("program_", "", sq_pp_hl$A_name))
sq_pp_hl$pb <- as.integer(sub("program_", "", sq_pp_hl$B_name))
bt_keys_from_sq <- paste(sq_pp_hl$pa, sq_pp_hl$pb)
OUT   <- file.path(ROOT, "figures/markcorr_B")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

BASE <- 6
theme_b <- theme_minimal(base_size = BASE) + theme(
  text          = element_text(size = BASE, colour = "grey15"),
  plot.title    = element_text(size = BASE + 1, face = "bold", hjust = 0, margin = margin(b = 2)),
  plot.subtitle = element_text(size = BASE - 0.5, colour = "grey40", hjust = 0, margin = margin(b = 3)),
  axis.title    = element_text(size = BASE),
  axis.text     = element_text(size = BASE - 1, colour = "grey25"),
  legend.title  = element_text(size = BASE - 0.5),
  legend.text   = element_text(size = BASE - 1),
  legend.key.size = unit(3.2, "mm"),
  legend.margin = margin(0,0,0,0),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(linewidth = 0.18, colour = "grey90"),
  plot.margin   = margin(4,4,4,4)
)
theme_set(theme_b)

CLASS_COL <- c(exc="#D7263D", inh="#1B6CA8", glia="#2E933C",
               nonneuron="#8E44AD", vascular="#E08A00")
DIV <- colorRampPalette(c("#2166AC","#4393C3","#D1E5F0","#F7F7F7",
                          "#FDDBC7","#D6604D","#B2182B"))(101)

# =============================================================================
# LOAD (abstract)
# =============================================================================
lab <- read.delim(file.path(ROOT,"results/crossregion_v1/program_labels.tsv"),
                   check.names=FALSE, stringsAsFactors=FALSE)
ann <- read.delim(file.path(ROOT,"results/crossregion_v1/program_annotation.tsv"),
                  check.names=FALSE, stringsAsFactors=FALSE)

# AUTHORITATIVE FDR table = results/crossregion_v1/program_names.tsv (brain_term_NES FDR;
# same table Figure A's editor used). program_annotation_gobp.tsv is STALE — do NOT use.
# Star = confidence of the DISPLAYED curated functional name -> bind to that program's FDR.
nmfdr <- read.delim(file.path(ROOT,"results/crossregion_v1/program_names.tsv"),
                    check.names=FALSE, stringsAsFactors=FALSE)
{
  .prog_col <- if ("new_P" %in% names(nmfdr)) nmfdr$new_P else nmfdr$program
  nmfdr$program <- .prog_col  # ensure $program alias exists
  nmfdr$pid <- as.integer(sub("^P", "", .prog_col))  # strip "P" prefix if present
}
# STAR AUTHORITY 20260615: weak-annotation star = program_names.tsv confidence column
# (the 24 "brain-weak" programs), NOT the raw fdr>0.05 numeric. The fdr>0.05 set
# diverged from the curated brain-weak whitelist (over/under-starred); confidence is
# the single source of truth for which displayed functional names are weak.
FDR_WEAK <- setNames(nmfdr$confidence == "brain-weak", nmfdr$pid)   # TRUE if brain-weak (weak annotation)
star_pid <- function(p) {
  # p = old pid; FDR_WEAK keyed by new pid
  np <- op2np(p)
  if (is.na(np)) return("")
  ifelse(isTRUE(FDR_WEAK[[as.character(np)]]), "*", "")
}  # trailing * if weak (RENUMBER 20260620)
# label a "Pn" token (e.g. used in pair strings "P8 x P13") with its FDR star
star_lab <- function(p) {
  np <- op2np(p)
  if (is.na(np)) return(paste0("P", p, "(excl)"))
  paste0("P", np, star_pid(p))
}  # display new pid (RENUMBER 20260620)
ann$pid <- as.integer(sub("^P","",ann$program))
lab$pid <- as.integer(lab$program)
meta <- merge(lab[,c("pid","label_short","confidence")], ann[,c("pid","class")], by="pid")
meta <- meta[order(meta$pid),]

# RENUMBER 20260620: load old->new pid map for display label conversion.
# All matrix indices (HUBS, pid, g25) remain old pid throughout.
# Only display strings (P# labels in figure text) switch to new pid.
rmap <- read.delim(file.path(ROOT,"results/crossregion_v1/program_renumber_map.tsv"),
                   check.names=FALSE, stringsAsFactors=FALSE)
rmap$old_pid_int <- as.integer(rmap$old_P)
rmap$new_pid_int <- suppressWarnings(as.integer(rmap$new_P))  # NA for EXCLUDED
old2new <- setNames(rmap$new_pid_int, as.character(rmap$old_pid_int))
# helper: given old pid integer, return new pid integer (or NA if excluded)
op2np <- function(op) { r <- old2new[as.character(as.integer(op))]; if(is.na(r)) NA_integer_ else as.integer(r) }
# NAMEFIX 20260614: authoritative display name = program_names.tsv name_short
# (NOT stale program_labels.tsv label_short). Strip any trailing '*' so the
# weak-FDR star is added ONLY by star_pid() dynamically (avoid double star).
meta$short <- sub("\\*+$", "", nmfdr$name_short[match(vapply(meta$pid, op2np, integer(1)), nmfdr$pid)])  # RENUMBER: match via new pid
meta$plab  <- paste0("P", vapply(meta$pid, op2np, integer(1)), " ", meta$short)
CLASS_ORD <- c("exc","inh","glia","nonneuron","vascular")
meta$class <- factor(meta$class, levels=CLASS_ORD)
N <- nrow(meta)

# COHORT-TECHNICAL EXCLUSION (consistency P0): six programs are cohort/batch-driven
# technical artefacts, NOT biological members — they must NOT be shown as members in
# any Fig.8 panel (network nodes, headline co-orgs, exemplars, avoidances, span).
# Filter them out of every data-driven selection below.
COHORT_TECH <- c(9, 18, 19, 35, 52, 57)
keep_pid <- function(p) !(as.integer(p) %in% COHORT_TECH)

zg  <- np$load(file.path(FIN,"progprog_median_iqr.npz"), allow_pickle=TRUE)
L2  <- zg$f[["log2_median_g"]]; FSS <- zg$f[["frac_same_sign"]]
RING<- as.numeric(zg$f[["ring_edges_um"]])
g25 <- L2[,,1]; fss25 <- FSS[,,1]
rownames(g25) <- meta$pid; colnames(g25) <- meta$pid

gcos <- as.matrix(read.csv(file.path(SIM,"gene_cosine.csv"), row.names=1, check.names=FALSE))
simp <- read.delim(file.path(SIM,"program_similarity_pairs.tsv"), stringsAsFactors=FALSE)
simp$pa <- as.integer(sub("program_","",simp$A)); simp$pb <- as.integer(sub("program_","",simp$B))

loo <- read.delim(file.path(RIG,"progprog_donor_loo.tsv"), stringsAsFactors=FALSE)
jac <- read.delim(file.path(RIG,"progprog_jaccard.tsv"), stringsAsFactors=FALSE)

zb  <- np$load(file.path(FIN,"progprog_byarea_median_iqr.npz"), allow_pickle=TRUE)
L2b <- zb$f[["log2_median_g"]]
AREA<- as.character(zb$f[["area_names"]])
UNST<- as.integer(zb$f[["unstable"]])
STABLE <- AREA[UNST==0]

pid2lab <- setNames(meta$plab, meta$pid)
pid2short<- setNames(meta$short, meta$pid)
pid2cls <- setNames(as.character(meta$class), meta$pid)

cluster_within_class_order <- function(mat, class_vec, class_levels){
  unlist(lapply(class_levels, function(cl){
    idx <- which(as.character(class_vec) == cl)
    if (length(idx) <= 2) return(idx)
    subm <- mat[idx, idx, drop=FALSE]
    idx[hclust(dist(subm))$order]
  }))
}

# =============================================================================
# PANEL a — 60x60 heatmap
# =============================================================================
ord  <- cluster_within_class_order(g25, meta$class, CLASS_ORD)
M    <- g25[ord, ord]
cls  <- as.character(meta$class)[ord]; pids <- meta$pid[ord]
cls_split <- factor(cls, levels=CLASS_ORD)
rng  <- max(abs(quantile(M, c(.01,.99), na.rm=TRUE)))
colf <- colorRamp2(seq(-rng, rng, length.out=101), DIV)
cbar <- HeatmapAnnotation(class=cls, col=list(class=CLASS_COL),
                          simple_anno_size=unit(1.4,"mm"), show_legend=FALSE,
                          annotation_name_gp=gpar(fontsize=BASE-1),
                          show_annotation_name=FALSE)
rbar <- rowAnnotation(class=cls, col=list(class=CLASS_COL),
                      simple_anno_size=unit(1.4,"mm"),
                      annotation_legend_param=list(title="Program class",
                        title_gp=gpar(fontsize=BASE-0.5), labels_gp=gpar(fontsize=BASE-1),
                        grid_height=unit(3,"mm"), grid_width=unit(3,"mm")),
                      show_annotation_name=FALSE)
ht <- Heatmap(M, name="log2 g\n(25 µm)", col=colf, rect_gp=gpar(col=NA),
              cluster_rows=FALSE, cluster_columns=FALSE,
              row_split=cls_split, column_split=cls_split,
              row_gap=unit(1.0, "mm"), column_gap=unit(1.0, "mm"),
              top_annotation=cbar, left_annotation=rbar,
              row_title=NULL,
              show_row_names=FALSE, show_column_names=FALSE,
              width=unit(46,"mm"), height=unit(46,"mm"),
              column_title="Program x program spatial co-organization (r = 25 µm)",
              column_title_gp=gpar(fontsize=BASE+1, fontface="bold"),
              heatmap_legend_param=list(title_gp=gpar(fontsize=BASE-0.5),
                labels_gp=gpar(fontsize=BASE-1), legend_height=unit(14,"mm"),
                grid_width=unit(3,"mm")))
panel_a <- grid.grabExpr(draw(ht, heatmap_legend_side="right",
                              annotation_legend_side="right", merge_legend=TRUE,
                              padding=unit(c(1,1,1,1),"mm")), wrap.grobs=TRUE)
panel_a <- wrap_elements(full=panel_a)

# =============================================================================
# PANEL b — co-org vs gene-similarity scatter
# =============================================================================
sb <- simp[simp$pa < simp$pb, ]
sb <- sb[keep_pid(sb$pa) & keep_pid(sb$pb), ]   # drop cohort-technical programs (P0)
ij <- cbind(match(sb$pa, meta$pid), match(sb$pb, meta$pid))
sb$g  <- g25[ij]; sb$absg <- abs(sb$g)
rr <- cor(sb$gene_cosine, sb$g, use="complete.obs")
sb$indep <- sb$gene_cosine < 0.25 & sb$g > 0.32
nb_indep <- sum(sb$indep)
set.seed(7)
yhi <- max(sb$g, na.rm=TRUE); ylo <- min(sb$g, na.rm=TRUE)
panel_b <- ggplot(sb, aes(gene_cosine, g)) +
  geom_hline(yintercept=c(-0.32,0.32), linetype=3, linewidth=0.2, colour="grey60") +
  geom_vline(xintercept=0.25, linetype=3, linewidth=0.2, colour="grey60") +
  geom_point(data=subset(sb,!indep), colour="grey72", size=0.35, alpha=0.45, stroke=0) +
  geom_point(data=subset(sb, indep), aes(gene_cosine), colour="#B2182B",
             size=0.55, alpha=0.55, stroke=0,
             position=position_jitter(width=0.012, height=0, seed=7)) +
  geom_smooth(method="lm", se=FALSE, colour="#1B6CA8", linewidth=0.4) +
  annotate("text", x=0.40, y=ylo + (yhi-ylo)*0.06,
           label=sprintf("Pearson r = %.2f", rr), hjust=0, size=2.0, colour="#1B6CA8") +
  annotate("text", x=0.40, y=ylo + (yhi-ylo)*0.20,
           label=sprintf("low-gene-overlap\nstrong co-orgs: %d", nb_indep),
           hjust=0, vjust=0, size=1.9, colour="#B2182B", lineheight=0.9) +
  scale_x_continuous(expand=expansion(mult=c(0.01,0.03))) +
  scale_y_continuous(expand=expansion(mult=c(0.02,0.02))) +
  labs(x="Gene-loading cosine similarity", y="Spatial co-org  log2 g (25 µm)",
       title="Co-organization rises with gene similarity",
       subtitle="strong co-orgs are mostly gene-similar co-regulated programs;\na distinct subset is low-gene-overlap") +
  coord_cartesian(clip="off")

# =============================================================================
# SPATIAL helper — single tissue map panel (ggplot geom_point, coord_fixed)
# masked-invalid greyed; per-panel vmin/vmax; viridis/magma; y inverted.
# =============================================================================
map_panel <- function(d, title, cmap="viridis", legend_lab="SCT",
                      title_size=BASE-0.5, pt=0.18, show_legend=TRUE){
  vmin <- d$vmin[1]; vmax <- d$vmax[1]
  pal  <- if (cmap=="magma") viridisLite::magma(256) else viridisLite::viridis(256)
  inv  <- d[!d$valid, ]; val <- d[d$valid & is.finite(d$val), ]
  ggplot() +
    {if (nrow(inv)>0) geom_point_r(data=inv, aes(x, y), colour="#dcdcdc",
                                 size=pt, shape=15, stroke=0)} +
    geom_point_r(data=val, aes(x, y, colour=val), size=pt, shape=15, stroke=0) +
    scale_colour_gradientn(colours=pal, limits=c(vmin, vmax), oob=scales::squish,
                           name=legend_lab,
                           guide=if (show_legend) guide_colourbar(barwidth=unit(1.4,"mm"),
                                                 barheight=unit(9,"mm")) else "none") +
    scale_y_reverse() +
    coord_fixed(clip="off") +
    labs(title=title) +
    theme_void(base_size=BASE) +
    theme(plot.title=element_text(size=title_size, hjust=0.5, margin=margin(b=1)),
          legend.title=element_text(size=BASE-1.5), legend.text=element_text(size=BASE-1.5),
          legend.key.height=unit(9,"mm"), legend.key.width=unit(1.4,"mm"),
          legend.margin=margin(0,0,0,0), legend.box.spacing=unit(1,"mm"),
          plot.margin=margin(1,1,1,1))
}

# shared compact reference bar for the SP1 grid (per-map relative scaling).
# RELAYOUT FIX: HORIZONTAL twin colorbars, placed in a thin band UNDER the 3x3
# grid (was a floating tall column on the right -> dead whitespace). Same content
# (per-map relative SCT; viridis = co-distribute programs, magma = avoid program).
shared_bar_grob <- function(){
  bdf <- data.frame(x=seq(0,1,length.out=256))
  mk <- function(pal, lab) ggplot(bdf, aes(x=x, y=1, fill=x)) +
    geom_raster() +
    scale_fill_gradientn(colours=pal, guide="none") +
    scale_x_continuous(breaks=c(0.02,0.98), labels=c("low","high"),
                       expand=expansion(0,0)) +
    scale_y_continuous(expand=expansion(0,0)) +
    labs(subtitle=lab) +
    theme_void(base_size=BASE) +
    theme(axis.text.x=element_text(size=BASE-1.5, vjust=1,
            margin=margin(t=0.4, unit="mm")),
          axis.ticks.length=unit(0,"mm"),
          plot.subtitle=element_text(size=BASE-1.5, hjust=0.5, margin=margin(b=0.8)),
          plot.margin=margin(0.5,2,0.5,2))
  vb <- mk(viridisLite::viridis(256), paste0("P", op2np(45), " / P", op2np(13), star_pid(13), " co-distribute"))
  mb <- mk(viridisLite::magma(256),  "P8 avoids")
  patchwork::wrap_plots(plot_spacer(), vb, plot_spacer(), mb, plot_spacer(),
                        nrow=1, widths=c(0.16,0.31,0.06,0.31,0.16)) +
    plot_annotation(title="SCT score (per-map relative scaling)",
      theme=theme(plot.title=element_text(size=BASE-1, face="bold", hjust=0.5,
                                          margin=margin(b=1)),
                  plot.margin=margin(1,2,0,2)))
}

# =============================================================================
# SP1 — spatial grounding grid (3 chips x 3 programs)
# =============================================================================
sp1 <- fread(file.path(SPR,"sp1_long.csv"))
REG_LAB <- c("DLPFC (prefrontal)", "V1 (occipital)", "M1 (motor)")  # row order 0,1,2
sp1_panels <- vector("list", 9)
for (rI in 0:2) for (cI in 0:2) {
  d <- sp1[sp1$row==rI & sp1$col==cI, ]
  cmap <- d$cmap[1]
  base_t <- d$ptitle[1]
  # NAMEFIX 20260614: ptitle now carries authoritative program_names name_short
  # (P8 Neurofilament cytoskeleton (pan-neuronal); P13 Oligodendrocyte development;
  # P45 Myelination). Weak-FDR star added dynamically via star_pid() (no hard-coded *).
  # RENUMBER 20260620: replace old P-numbers in ptitle with new P-numbers
  for(op_str in names(old2new)) {
    np_val <- old2new[[op_str]]
    if (!is.na(np_val)) {
      base_t <- gsub(paste0("(P", op_str, ")"), paste0("(P", np_val, ")"), base_t, fixed=TRUE)
    }
  }
  base_t <- sub(paste0("(P", op2np(13), ")"),
                paste0("(P", op2np(13), star_pid(13), ")"), base_t, fixed=TRUE)  # P13->new: star if brain-weak
  # add region prefix on first column only
  ttl <- if (cI==0) paste0(REG_LAB[rI+1], "\n", base_t) else base_t
  sp1_panels[[rI*3 + cI + 1]] <- map_panel(d, ttl, cmap=cmap, legend_lab="SCT",
                                           title_size=BASE-1, pt=0.22,
                                           show_legend=FALSE) +
    theme(plot.margin=margin(1,6,1,6))  # horizontal gutter so map titles do not collide (widened micro-touch)
}
sp1_maps <- wrap_plots(sp1_panels, nrow=3, ncol=3)
# RELAYOUT FIX: STACK the 3x3 map grid OVER a thin horizontal colorbar band
# (no more floating right column / dead whitespace). Maps reclaim the full width.
sp1_grid <- (sp1_maps / shared_bar_grob()) + plot_layout(heights=c(1, 0.07)) +
  plot_annotation(title=sprintf("Two oligodendroglial programs (P%d, P%d%s) co-distribute in tissue; neurofilament-cytoskeleton program (P8) avoids",
                                op2np(45), op2np(13), star_pid(13)),
                  theme=theme(plot.title=element_text(size=BASE+0.5, face="bold", hjust=0,
                                                      margin=margin(b=1)),
                              plot.margin=margin(2,2,2,2)))
sp1_block <- wrap_elements(full=sp1_grid)

# =============================================================================
# PANEL c — low-gene-overlap niche network (microglia-synapse hub)
# =============================================================================
ed <- sb[sb$gene_cosine < 0.25 & abs(sb$g) > 0.32, c("pa","pb","g","gene_cosine")]
ed <- ed[ed$g > 0, ]
deg <- table(c(ed$pa, ed$pb))
usei <- as.integer(names(deg))
vdf <- data.frame(pid=usei, deg=as.integer(deg),
                  class=pid2cls[as.character(usei)], short=pid2short[as.character(usei)])
vdf$absg_strength <- sapply(vdf$pid, function(p) sum(abs(ed$g[ed$pa==p | ed$pb==p])))
gph <- graph_from_data_frame(d=data.frame(from=ed$pa, to=ed$pb, w=ed$g, aw=abs(ed$g)),
                             vertices=vdf, directed=FALSE)
HUBS <- c(40,49,60)
hub_partner <- function(h, n=6){
  pr <- ed[ed$pa==h | ed$pb==h, ]; pr$partner <- ifelse(pr$pa==h, pr$pb, pr$pa)
  pr <- pr[order(-pr$g), ]; head(unique(pr$partner), n)
}
PART <- unique(unlist(lapply(HUBS, hub_partner, n=4)))
PART <- setdiff(PART, HUBS)
set.seed(11)
lay <- create_layout(gph, layout="stress")
xr <- range(lay$x); yr <- range(lay$y)
sx <- function(f) xr[1] + diff(xr)*f
sy <- function(f) yr[1] + diff(yr)*f
hub_xy <- list("40"=c(sx(0.18), sy(0.80)), "49"=c(sx(0.82), sy(0.78)),
               "60"=c(sx(0.50), sy(0.18)))
for(h in names(hub_xy)){
  idx <- which(as.character(lay$name) == h)
  if(length(idx)==1){ lay$x[idx] <- hub_xy[[h]][1]; lay$y[idx] <- hub_xy[[h]][2] }
}
place_partners <- function(h, angle_start, angle_end, radius){
  partners <- hub_partner(as.integer(h), n=6)
  partners <- partners[partners %in% as.integer(lay$name)]
  if (!length(partners)) return(invisible(NULL))
  ang <- seq(angle_start, angle_end, length.out=length(partners))
  hx <- hub_xy[[as.character(h)]][1]
  hy <- hub_xy[[as.character(h)]][2]
  for (i in seq_along(partners)){
    idx <- which(as.integer(lay$name) == partners[i])
    if (!length(idx)) next
    lay$x[idx] <<- hx + cos(ang[i]) * radius
    lay$y[idx] <<- hy + sin(ang[i]) * radius
  }
}
place_partners(40, pi * 0.70, pi * 1.35, diff(xr) * 0.12)
place_partners(49, -pi * 0.35, pi * 0.30, diff(xr) * 0.12)
place_partners(60, pi * 1.55, pi * 1.95, diff(xr) * 0.12)
lay$is_hub  <- as.integer(lay$name) %in% HUBS
lay$is_part <- as.integer(lay$name) %in% PART
lay$lab <- NA_character_
lay$lab[lay$is_hub]  <- sapply(lay$name[lay$is_hub],  star_lab)  # P# + FDR star if weak
lay$lab[lay$is_part] <- sapply(lay$name[lay$is_part], star_lab)
lay$lab_bold <- lay$is_hub
ename <- as_edgelist(gph, names=TRUE)
exy <- data.frame(
  x   = lay$x[match(ename[,1], lay$name)], y   = lay$y[match(ename[,1], lay$name)],
  xend= lay$x[match(ename[,2], lay$name)], yend= lay$y[match(ename[,2], lay$name)],
  aw  = abs(E(gph)$w),
  hub = (as.integer(ename[,1]) %in% HUBS) | (as.integer(ename[,2]) %in% HUBS)
)
exy <- exy[order(exy$hub), ]
panel_c <- ggplot() +
  geom_segment(data=subset(exy,!hub),
               aes(x=x,y=y,xend=xend,yend=yend, linewidth=aw, alpha=aw),
               colour="#C9A227", lineend="round") +
  geom_segment(data=subset(exy, hub),
               aes(x=x,y=y,xend=xend,yend=yend, linewidth=aw, alpha=aw),
               colour="#8E44AD", lineend="round") +
  scale_linewidth(range=c(0.10,0.85), guide="none") +
  scale_alpha(range=c(0.12,0.62), guide="none") +
  geom_point(data=subset(lay, is_hub), aes(x,y), size=10, colour="#8E44AD", alpha=0.14) +
  geom_point(data=subset(lay, is_hub), aes(x,y), size=7,  colour="#8E44AD", alpha=0.20) +
  geom_point(data=lay, aes(x,y, size=deg, fill=class), shape=21, colour="white", stroke=0.2) +
  scale_size(range=c(0.9,6.2), name="co-org degree",
             guide=guide_legend(override.aes=list(fill="grey50"))) +
  scale_fill_manual(values=CLASS_COL, name="class",
                    guide=guide_legend(override.aes=list(size=2.6, colour="white", stroke=0.2))) +
  geom_text_repel(data=subset(lay, lab_bold), aes(x,y,label=lab),
                  size=2.6, fontface="bold", colour="#4A235A",
                  max.overlaps=Inf, min.segment.length=0,
                  segment.size=0.2, segment.colour="grey55",
                  box.padding=0.5, point.padding=0.3, seed=3) +
  geom_text_repel(data=subset(lay, is_part), aes(x,y,label=lab),
                  size=1.9, colour="grey15",
                  max.overlaps=Inf, min.segment.length=0,
                  segment.size=0.12, segment.colour="grey75",
                  box.padding=0.28, seed=3) +
  coord_equal(clip="off") +
  labs(title="Low-gene-overlap co-organization niches",
       subtitle="edges: gene cosine < 0.25 & |log2 g| > 0.32 – the microglia–synapse hub",
       caption=paste0("Microglial programs P", op2np(40), " / P", op2np(49), " / P", op2np(60),
                       " (purple, pinned) form low-gene-overlap hubs co-locating with synaptic programs (",
                       star_lab(17), ", ", star_lab(29), ", ", star_lab(47), ");\n",
                       "the adjacent inflammatory program P", op2np(54),
                       " (complement/MHC) forms no hub (degree 1).  *=weak functional annotation (brain-weak).")) +
  theme_void(base_size=BASE) +
  theme(plot.title=element_text(size=BASE+1,face="bold"),
        plot.subtitle=element_text(size=BASE-0.5,colour="grey40"),
        plot.caption=element_text(size=BASE-1.5, colour="grey35", hjust=0,
                                  lineheight=1.0, margin=margin(t=3)),
        legend.title=element_text(size=BASE-0.5), legend.text=element_text(size=BASE-1),
        legend.key.size=unit(3,"mm"), legend.position="right",
        legend.box="vertical",
        plot.margin=margin(4,6,4,6))

# =============================================================================
# SP2 — the SAME hub IN TISSUE: P40 micro | P29 synapse | co-distribution scatter
# =============================================================================
sp2m <- fread(file.path(SPR,"sp2_maps.csv"))
sp2s <- fread(file.path(SPR,"sp2_scatter.csv"))
m2 <- readLines(file.path(SPR,"sp2_meta.txt"))
getm <- function(k) sub(paste0("^",k,"\t"), "", grep(paste0("^",k,"\t"), m2, value=TRUE))
sp2_r   <- as.numeric(getm("r"))
sp2_reg <- getm("region")
progs2 <- unique(sp2m$prog)
# NAMEFIX 20260614: authoritative SP2 display names from program_names.tsv name_short
# (P40 = Microglial immune activation; P29 = Glutamate-receptor (postsyn.));
# star added dynamically by star_pid().
nm_P40 <- paste0("P", op2np(40), " ", pid2short[["40"]], star_pid(40))   # "P36 Microglial immune activation (RENUMBER)"
nm_P29 <- paste0("P", op2np(29), " ", pid2short[["29"]], star_pid(29))   # "P26 Glutamate-receptor (postsyn.) (RENUMBER)"
mapA <- map_panel(sp2m[sp2m$prog==progs2[1], ],
                  paste0(sp2_reg, "\n", nm_P40),
                  cmap="viridis", title_size=BASE-1, show_legend=FALSE)
mapB <- map_panel(sp2m[sp2m$prog==progs2[2], ],
                  nm_P29, cmap="viridis", title_size=BASE-1,
                  show_legend=TRUE)
scat2 <- ggplot(sp2s, aes(micro, syn)) +
  geom_point_r(colour="#444444", size=0.18, alpha=0.13, stroke=0) +
  geom_smooth(method="lm", se=FALSE, colour="#B2182B", linewidth=0.4) +
  annotate("text", x=-Inf, y=Inf, hjust=-0.12, vjust=1.4, size=2.1, colour="#B2182B",
           label=sprintf("r = %.2f", sp2_r)) +
  labs(x=nm_P40, y=nm_P29,
       title=sprintf("Representative hub edge: P%d x P%d (per-bin r = %.2f)", op2np(40), op2np(29), sp2_r)) +
  theme(plot.title=element_text(size=BASE-1, face="plain", hjust=0.5),
        axis.text=element_text(size=BASE-1.5), axis.title=element_text(size=BASE-1),
        legend.position="none",
        panel.grid.major=element_line(linewidth=0.12, colour="grey92"),
        plot.margin=margin(2,3,2,2))
sp2_grid <- (mapA | mapB | scat2) + plot_layout(widths=c(1,1,1)) +
  plot_annotation(title=sprintf("The microglia-synapse hub in tissue (V1, low-gene-overlap); scatter = one representative edge P%d x P%d, r = %.2f", op2np(40), op2np(29), sp2_r),
                  theme=theme(plot.title=element_text(size=BASE+0.5, face="bold", hjust=0,
                                                      margin=margin(b=1)),
                              plot.margin=margin(2,2,2,2)))
sp2_block <- wrap_elements(full=sp2_grid)

# =============================================================================
# PANEL d — strongest avoidances lollipop
# =============================================================================
neur <- c(8,24); olig <- c(36,26,37,13)   # P57 removed (cohort-technical, P0)
av <- sb[ (sb$pa %in% neur & sb$pb %in% olig) | (sb$pa %in% olig & sb$pb %in% neur), ]
av$fss <- fss25[cbind(match(av$pa,meta$pid), match(av$pb,meta$pid))]
av$pair <- paste0(sapply(av$pa, star_lab)," x ",sapply(av$pb, star_lab))  # FDR star per pid
av <- av[order(av$g), ]; av$pair <- factor(av$pair, levels=av$pair)
panel_d <- ggplot(av, aes(g, pair)) +
  geom_vline(xintercept=0, linewidth=0.2, colour="grey70") +
  geom_segment(aes(x=0, xend=g, yend=pair), colour="#2166AC", linewidth=0.4) +
  geom_point(aes(fill=fss), shape=21, size=2.0, colour="white", stroke=0.2) +
  scale_fill_gradient(low="#9ecae1", high="#08306b", name="frac\nsame-sign",
                      limits=c(0.85,1), breaks=c(0.9,1)) +
  labs(x="Spatial co-org  log2 g (25 µm)", y=NULL,
       title="Strongest spatial avoidances",
       subtitle=sprintf("Neurofilament cytoskeleton (%s,%s) vs\noligodendrocyte/myelin (%s,%s,%s,%s)",
                       star_lab(8), star_lab(24), star_lab(13), star_lab(26), star_lab(36), star_lab(37)),
       caption="Reflects a gray-/white-matter spatial partition (not program-level repulsion of neurons by oligodendrocytes):\noligodendrocyte-diff. P38 instead CO-occurs (log2 g +0.2) as reference control.  *=weak functional annotation (brain-weak).") +
  theme(axis.text.y=element_text(size=BASE-1.5),
        plot.caption=element_text(size=BASE-1.5, colour="grey35", hjust=0,
                                  lineheight=1.0, margin=margin(t=3)))

# =============================================================================
# PANEL e — rigor
# =============================================================================
ut <- which(upper.tri(g25), arr.ind=TRUE)
# drop any co-org pair involving a cohort-technical program (P0) before all
# data-driven top-N selections (panel e headline hist, g region-clustering, h span)
ut <- ut[keep_pid(meta$pid[ut[,1]]) & keep_pid(meta$pid[ut[,2]]), ]
# BETWEENCHIP NULL UPDATE (2026-06-25): use sq_pp is_headline (1208 pairs) instead of |g|>0.32
hl_all <- data.frame(g=g25[ut], fss=fss25[ut],
                     pa=meta$pid[ut[,1]], pb=meta$pid[ut[,2]])
hl_all$pair_key <- paste(hl_all$pa, hl_all$pb)
hl <- hl_all[hl_all$pair_key %in% bt_keys_from_sq, ]  # 1208 between-chip headline pairs
e1 <- ggplot(hl, aes(fss)) +
  geom_histogram(bins=28, fill="#1B6CA8", colour="white", linewidth=0.1) +
  geom_vline(xintercept=median(hl$fss), linetype=2, linewidth=0.3, colour="#B2182B") +
  annotate("text", x=median(hl$fss), y=Inf, vjust=1.4, hjust=1.1, size=1.9, colour="#B2182B",
           label=sprintf("median %.2f", median(hl$fss))) +
  labs(x="frac same-sign across 44 chips", y="headline co-orgs",
       title="Reproducibility of 1,208 headline co-orgs (between-chip null)",
       subtitle="direction agreement across chips + donor-LOO stability;\nsanity check: GM/WM compartment recovery (ground-truth)") +
  coord_cartesian(clip="off")
loo$absfull <- abs(loo$full_log2)
lh <- loo[loo$absfull > 0.32, ]
lh <- lh[keep_pid(lh$A_idx + 1) & keep_pid(lh$B_idx + 1), ]   # drop cohort-technical (P0)
lh <- lh[order(-lh$absfull), ]; lh <- head(lh, 22)
lh$pair <- paste0(sapply(lh$A_idx+1, star_lab)," x ",sapply(lh$B_idx+1, star_lab))  # FDR star per pid
lh$pair <- factor(lh$pair, levels=rev(lh$pair))
e2 <- ggplot(lh, aes(y=pair)) +
  geom_vline(xintercept=c(-0.32,0.32), linetype=3, linewidth=0.2, colour="grey70") +
  geom_segment(aes(x=log2(loo_min_median_g), xend=log2(loo_max_median_g), yend=pair),
               linewidth=0.5, colour="grey55") +
  geom_point(aes(x=full_log2, colour=stable_90), size=0.9) +
  scale_colour_manual(values=c("True"="#2E933C","False"="grey55"), name="LOO stable (90%)",
                      labels=c("True"="yes","False"="no")) +
  labs(x="log2 g  (point=full; bar=donor-LOO min-max)", y=NULL, title=NULL, subtitle=NULL) +
  theme(axis.text.y=element_text(size=BASE-2))
panel_e <- wrap_elements(full = (e1 + e2 + plot_layout(widths=c(1,1.15))))

# =============================================================================
# PANEL f — region small-multiples
# =============================================================================
ex <- list(c(37,38), c(8,36), c(40,29), c(17,47))   # P17xP52 (cohort-technical P52) -> P17xP47 valid synaptic co-org (P0)
exlab <- c(paste0(star_lab(37),"x",star_lab(38)," oligo/myelin"),
           paste0(star_lab(8),"x",star_lab(36)," neuron-oligo avoid"),
           paste0(star_lab(40),"x",star_lab(29)," microglia-glutR"),
           paste0(star_lab(17),"x",star_lab(47)," synaptic assembly"))  # weak-annotation star per pid via star_lab (brain-weak whitelist)
ai <- match(STABLE, AREA)
frows <- do.call(rbind, lapply(seq_along(ex), function(k){
  pa<-ex[[k]][1]; pb<-ex[[k]][2]
  data.frame(area=STABLE, ex=exlab[k], g=sapply(ai, function(a) L2b[a, pa, pb, 1]))
}))
frows$ex <- factor(frows$ex, levels=exlab)
panel_f <- ggplot(frows, aes(reorder(area,g), g)) +
  geom_hline(yintercept=0, linewidth=0.2, colour="grey70") +
  geom_col(aes(fill=g>0), width=0.7, show.legend=FALSE) +
  scale_fill_manual(values=c("TRUE"="#B2182B","FALSE"="#2166AC")) +
  scale_y_continuous(expand=expansion(mult=c(0.04,0.06))) +
  facet_wrap(~ex, nrow=1) +
  labs(x=NULL, y="log2 g (25 µm)",
       title="Exemplar co-orgs across 9 stable regions",
       subtitle="per-area median g; unstable regions (n<3 chips) excluded") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=BASE-1),
        axis.title.y=element_text(size=BASE),
        panel.spacing=unit(2,"mm"), strip.text=element_text(size=BASE-1),
        plot.margin=margin(4,4,4,8))

# =============================================================================
# PANEL g — region clustering
# =============================================================================
hp <- data.frame(a=meta$pid[ut[,1]], b=meta$pid[ut[,2]], g=g25[ut])
hp <- hp[order(-abs(hp$g)), ]; hp <- head(hp, 40)
Mg <- sapply(seq_len(nrow(hp)), function(k){
  sapply(ai, function(a) L2b[a, hp$a[k], hp$b[k], 1])
})
rownames(Mg) <- STABLE
colnames(Mg) <- paste0(sapply(hp$a, star_lab),"x",sapply(hp$b, star_lab))  # FDR star per pid (shown in panel h y-axis)
SENS <- c("V1","S1","M1")
gcol <- ifelse(STABLE %in% SENS, "#E08A00", "#1B6CA8")
rng2 <- max(abs(Mg), na.rm=TRUE)
htg <- Heatmap(Mg, name="log2 g", col=colorRamp2(seq(-rng2,rng2,length.out=101), DIV),
               cluster_columns=TRUE, cluster_rows=TRUE, show_column_names=FALSE,
               row_names_gp=gpar(fontsize=BASE-0.5, col=gcol), row_names_side="left",
               row_dend_width=unit(6,"mm"), column_dend_height=unit(5,"mm"),
               width=unit(52,"mm"), height=unit(26,"mm"),
               column_title="Region clustering over 40 headline co-orgs",
               column_title_gp=gpar(fontsize=BASE+1, fontface="bold"),
               heatmap_legend_param=list(title_gp=gpar(fontsize=BASE-0.5),
                 labels_gp=gpar(fontsize=BASE-1), legend_height=unit(12,"mm"),
                 grid_width=unit(3,"mm")))
panel_g <- grid.grabExpr(draw(htg, padding=unit(c(1,1,1,1),"mm")), wrap.grobs=TRUE)
panel_g <- wrap_elements(full=panel_g)

# =============================================================================
# PANEL h — cross-region span
# =============================================================================
span <- data.frame(pair=colnames(Mg), mn=apply(Mg,2,min,na.rm=TRUE), mx=apply(Mg,2,max,na.rm=TRUE))
span$span <- span$mx - span$mn
span <- span[order(-span$span), ]; span <- head(span, 22)
span$pair <- factor(span$pair, levels=rev(span$pair))
panel_h <- ggplot(span, aes(y=pair)) +
  geom_vline(xintercept=0, linewidth=0.2, colour="grey70") +
  geom_segment(aes(x=mn, xend=mx, yend=pair), colour="grey55", linewidth=0.4) +
  geom_point(aes(x=mn), colour="#2166AC", size=0.9) +
  geom_point(aes(x=mx), colour="#B2182B", size=0.9) +
  labs(x="log2 g range across stable regions (blue=min, red=max)", y=NULL,
       title="Most region-modulated co-orgs",
       subtitle="per co-org min->max across 9 stable regions, sorted by span") +
  theme(axis.text.y=element_text(size=BASE-2))

# =============================================================================
# HUB BLOCK (tag d) — microglia-synapse hub: network (top, wide) OVER the SAME
# hub IN TISSUE (P36 micro | P26 synapse | representative-edge scatter). One
# tight composite -> one tag. Dual-modality (graph + tissue + scatter) intact.
# =============================================================================
sp2_row <- (mapA | mapB | scat2) + plot_layout(widths=c(1,1,1.05))
hub_block <- (panel_c / sp2_row) + plot_layout(heights=c(1.45, 1)) +
  plot_annotation(theme=theme(plot.margin=margin(0,0,0,0)))
hub_block <- wrap_elements(full=hub_block)

# =============================================================================
# SVGUTILS RELAYOUT export tail — write 9 standalone ZERO-MARGIN panel SVGs
# (content FROZEN; only layout/packing changes). Each panel -> tight svglite
# SVG at ~natural size; ink-crop trims to true bbox; make_template_nested packs
# them gutter-0 tangram; compose_svgutils_skill stitches to vector PDF.
#   a=heatmap b=scatter c=sp1 spatial grid d=hub(network+tissue) e=avoid lollipop
#   f=rigor g=region exemplars h=region-clustering heatmap i=span
# =============================================================================
# Worker override: keep live source/data root, but let this run emit panels into
# a temp svg_panels directory so layout trials do not touch canonical outputs.
SVGD <- Sys.getenv("FIGB_SVGD", file.path(ROOT, "scripts/figmarkcorr_B/svg_panels"))
dir.create(SVGD, recursive=TRUE, showWarnings=FALSE)
mm2in <- function(x) as.numeric(x)/25.4
svgf  <- function(id) file.path(SVGD, sprintf("figB_%s.svg", id))
zero_m <- theme(plot.margin = margin(0,0,0,0))
emit <- function(id, obj, w_mm, h_mm){
  svglite(svgf(id), width=mm2in(w_mm), height=mm2in(h_mm), bg="white")
  print(obj); invisible(dev.off())
  cat(sprintf("panel %s SVG (%.0fx%.0f mm)\n", id, w_mm, h_mm))
}

# t : figure title band (FROZEN content: honest title + FDR-star legend). Full-
# width thin band, top of figure (was plot_annotation title/caption pre-relayout).
ttl_b <- "Low-gene-overlap spatial co-organization of cortical gene programs"
cap_b <- "*=weak functional annotation (brain-weak).  Co-org = spatial mark-correlation g(r=25 µm), 44 chips; bin50 = 25 µm (multi-cell, not single cells)."
ptitle_b <- ggplot() + theme_void() +
  annotate("text", x=0.012, y=0.78, hjust=0, vjust=1, size=9/.pt, fontface="bold", label=ttl_b) +
  annotate("text", x=0.012, y=0.30, hjust=0, vjust=1, size=5.5/.pt, colour="grey30",
           lineheight=1.0, label=cap_b) +
  scale_x_continuous(limits=c(0,1), expand=expansion(0,0)) +
  scale_y_continuous(limits=c(0,1), expand=expansion(0,0)) +
  theme(plot.margin=margin(0,1,0,1))
# CRITICAL: svglite's default fix_text_size=TRUE writes textLength +
# lengthAdjust="spacingAndGlyphs" on every <text>. That renders fine standalone,
# but once svgutils wraps the band in a scale() transform, rsvg mis-lays the
# glyphs and collapses the whole title to a single stray "(" in the composite.
# fix_text_size=FALSE emits plain <text> (no textLength) -> composes correctly.
svglite(svgf("t"), width=mm2in(172), height=mm2in(13), bg="white", fix_text_size=FALSE)
print(ptitle_b); invisible(dev.off())
cat("panel t SVG (172x13 mm, fix_text_size=FALSE)\n")

# a : 60x60 heatmap (grabbed grob, tight padding)
pa_tight <- grid.grabExpr(draw(ht, heatmap_legend_side="right",
        annotation_legend_side="right", merge_legend=TRUE,
        padding=unit(c(1,1,1,1),"mm")), wrap.grobs=TRUE)
emit("a", wrap_elements(full=pa_tight), 78, 70)

# b : co-org vs gene-cosine scatter
emit("b", panel_b + zero_m, 60, 58)

# c : SP1 spatial program grid (3x3 maps over horizontal colorbar band)
emit("c", sp1_block & zero_m, 120, 112)

# d : microglia-synapse hub (network + tissue maps + scatter)
emit("d", hub_block & zero_m, 120, 120)

# e/f/g row is the embed pinch-point. Make e/f slightly taller and narrow g a
# touch so the row gains vertical height without changing panel order/content.
# e : strongest avoidances lollipop
emit("e", panel_d + zero_m, 70, 64)

# f : rigor (reproducibility hist + donor-LOO dumbbell)
emit("f", panel_e & zero_m, 96, 64)

# g : region exemplars (4-facet)
emit("g", panel_f + zero_m, 110, 54)

# h : region-clustering heatmap (grob)
emit("h", panel_g, 78, 56)

# i : cross-region span
emit("i", panel_h + zero_m, 86, 56)

cat("\nALL SVG PANELS WRITTEN (incl. title band t) ->", SVGD, "\n")
cat("DONE. Pearson r:", round(rr,3),
    "| sp2 micro-syn r:", round(sp2_r,3),
    "| stable regions:", length(STABLE), "\n")

#!/usr/bin/env Rscript
# Cross-region program conservation vs rewiring — 7+ panels, zero-margin vector SVG each.
# Panels: a lollipop conservation rank, b 14x14 region-pair heatmap, c expr-vs-neigh scatter
# (rewiring quadrant), d conservation-distance curve, e hit partner networks (ggraph),
# f hit spatial fields (bin50), g multi-step filter funnel.
# Display labels map source component IDs to contiguous retained-program IDs.
suppressMessages({
  library(ggplot2); library(data.table); library(ggrepel)
  library(ggraph); library(igraph); library(scales); library(svglite)
})
ROOT <- Sys.getenv("CORTEX_NMF_ROOT", "CORTEX_PROGRAM_ROOT")
B   <- file.path(ROOT, "results/crossregion_v1")
OUT <- Sys.getenv("XREGION_OUTPUT_DIR", file.path(B, "xregion_auroc"))
FIG <- Sys.getenv("FIGS6_OUTPUT_DIR", file.path(ROOT, "figures/xregion_auroc"))
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
PANELS <- file.path(FIG, "panels"); dir.create(PANELS, showWarnings=FALSE)

BASefont <- 7
th <- theme_minimal(base_size=BASefont) +
  theme(plot.margin=margin(0,0,0,0), panel.grid.minor=element_blank(),
        panel.grid.major=element_line(linewidth=.2, colour="grey90"),
        axis.text=element_text(size=BASefont-1, colour="grey20"),
        axis.title=element_text(size=BASefont, colour="grey10"),
        plot.title=element_text(size=BASefont+1, face="bold"),
        legend.title=element_text(size=BASefont-1), legend.text=element_text(size=BASefont-2),
        legend.key.size=unit(3,"mm"))
sv <- function(p, name, w, h) {
  svglite(file.path(PANELS, paste0(name,".svg")), width=w, height=h, bg="white")
  print(p); invisible(dev.off())
  cat("panel", name, "->", w, "x", h, "in\n")
}
# Map source component IDs to retained display labels.
nm <- fread(Sys.getenv("RETAIN_MAP", file.path(ROOT,"tables/TableS3_program_annotation.tsv")))
setnames(nm, "functional_name", "name_short")
# nm$new_P is "P1".."P54" or "EXCLUDED"; nm$cnmf_component is the old integer 1..60
# build lab from new_P + name_short (confidence stars preserved)
nm[, lab := ifelse(new_P == "EXCLUDED", NA_character_,
             ifelse(confidence=="brain-weak",
                    paste0(new_P," ",name_short," *"),
                    paste0(new_P," ",name_short)))]
# namemap: key = old cnmf_component integer (as character), value = new-P label
namemap <- setNames(nm$lab, sub("^P", "", nm$new_P))
# old2new: key = old integer (as character), value = new_P string e.g. "P21"
old2new <- setNames(nm$new_P, as.character(nm$cnmf_component))
stopifnot(nrow(nm) == 54, identical(sort(as.integer(names(namemap))), 1:54))
cat("namemap sample (retained P21):", namemap["21"], "\n")
# -----------------------------------------------------------------------

## ---------- a: lollipop conservation ranking (expr vs neigh, top+bottom) ----------
m1 <- fread(file.path(OUT,"m1_expr_conservation_per_program.tsv"))
m2b<- fread(file.path(OUT,"m2b_neighborhood_conservation_per_program.tsv"))
M  <- merge(m1, m2b, by="program")
M[, lab := namemap[as.character(program)]]
M[, rewgap := expr_cons_auroc - neigh_cons_auroc]
# show top12 + bottom12 by neighborhood conservation
ord <- M[order(neigh_cons_auroc)]
sel <- rbind(head(ord,12), tail(ord,12))
sel[, lab := factor(lab, levels=lab[order(neigh_cons_auroc)])]
pa <- ggplot(sel) +
  geom_segment(aes(y=lab, yend=lab, x=neigh_cons_auroc, xend=expr_cons_auroc), colour="grey70", linewidth=.5) +
  geom_point(aes(y=lab, x=expr_cons_auroc, colour="Expression"), size=1.6) +
  geom_point(aes(y=lab, x=neigh_cons_auroc, colour="Co-activation"), size=1.6) +
  scale_colour_manual(values=c(Expression="#2C6FB0", `Co-activation`="#D1495B"), name=NULL) +
  scale_x_continuous("AUROC", limits=c(0.5,1), breaks=seq(.5,1,.1)) +
  labs(y=NULL, title="Conservation") +
  th + theme(legend.position=c(.18,.5), legend.background=element_rect(fill=alpha("white",.7),colour=NA))
sv(pa, "a_lollipop", 2.41, 2.54)

## ---------- b: 14x14 region-pair expression conservation heatmap (square cells) ----------
pmdf <- as.data.frame(fread(file.path(OUT,"m1_expr_pairwise_auroc_matrix.tsv")))
rownames(pmdf) <- pmdf[[1]]; pmdf[[1]] <- NULL
pm <- as.matrix(pmdf)
storage.mode(pm) <- "numeric"
# order regions by AP using region_pc_coords PC1
pc <- fread(file.path(B,"region_pc_coords.tsv"))
roword <- pc[order(PC1)]$region
pm <- pm[roword, roword]
dd <- as.data.table(reshape2::melt(pm)); setnames(dd, c("rA","rB","auroc"))
dd[, auroc := as.numeric(auroc)]
dd[, rA:=factor(rA,levels=roword)][, rB:=factor(rB,levels=rev(roword))]
pb <- ggplot(dd, aes(rA, rB, fill=auroc)) +
  geom_tile(colour="white", linewidth=.3) +
  scale_fill_gradientn("AUC", colours=c("#3B4CC0","#EAEAEA","#B40426"),
                       values=rescale(c(.85,.93,1)), limits=c(.85,1)) +
  coord_fixed() + labs(x=NULL,y=NULL,title="Expression conservation\nby region pair") +
  th + theme(plot.title=element_text(size=BASefont, face="bold", lineheight=.95),
             axis.text.x=element_text(angle=90, vjust=.5, hjust=1, size=BASefont-2),
             axis.text.y=element_text(size=BASefont-2),
             plot.margin=margin(1,2,1,1))  # editfig_r1: 2-line title + small right pad
sv(pb, "b_pairheat", 2.79, 2.54)

## ---------- c: expression vs neighborhood scatter (rewiring quadrant) ----------
hits <- fread(file.path(OUT,"m4_rewiring_hits.tsv"))
M[, is_hit := program %in% hits$program]
xm <- mean(M$expr_cons_auroc); ym <- mean(M$neigh_cons_auroc)
pc3 <- ggplot(M, aes(expr_cons_auroc, neigh_cons_auroc)) +
  annotate("rect", xmin=xm, xmax=1.0, ymin=0.55, ymax=ym, fill="#F4D03F", alpha=.18) +
  annotate("text", x=0.995, y=0.60, label="Rewiring\nquadrant", hjust=1, size=2.0, colour="#9A7D0A", lineheight=.9) +
  geom_hline(yintercept=ym, linetype=2, colour="grey60", linewidth=.3) +
  geom_vline(xintercept=xm, linetype=2, colour="grey60", linewidth=.3) +
  geom_point(aes(colour=is_hit, size=is_hit)) +
  ggrepel::geom_text_repel(data=M[is_hit==TRUE], aes(label=lab), size=2.0, box.padding=.4,
                           min.segment.length=0, max.overlaps=20, colour="#7B241C") +
  scale_colour_manual(values=c(`FALSE`="grey65",`TRUE`="#C0392B"), guide="none") +
  scale_size_manual(values=c(`FALSE`=1.1,`TRUE`=2.2), guide="none") +
  coord_fixed(ratio=1, xlim=c(.74,1), ylim=c(.55,.95)) +
  labs(x="Expression AUROC", y="Co-activation AUROC",
       title="Expression vs\nneighborhood conservation") +
  th + theme(plot.title=element_text(size=BASefont, face="bold", lineheight=.95))
sv(pc3, "c_scatter", 1.89, 2.54)

## ---------- d: conservation vs region distance curve ----------
pts <- fread(file.path(OUT,"m3_distance_curve_points.tsv"))
st  <- fread(file.path(OUT,"m3_distance_curve_stats.tsv"))
long <- rbind(
  data.table(dist=pts$trans_dist, auroc=pts$expr_auroc, metric="Expression"),
  data.table(dist=pts$trans_dist, auroc=pts$neigh_auroc, metric="Co-activation"))
rE <- st[conservation_metric=="expr"&distance=="transcriptomic"]
rN <- st[conservation_metric=="neigh"&distance=="transcriptomic"]
subt <- sprintf("Expr rho=%.2f, Co-act rho=%.2f (perm p<1e-3)", rE$spearman_rho, rN$spearman_rho)
pd4 <- ggplot(long, aes(dist, auroc, colour=metric)) +
  geom_point(size=1, alpha=.6) +
  geom_smooth(method="lm", se=TRUE, linewidth=.6, alpha=.15) +
  scale_colour_manual(values=c(Expression="#2C6FB0",`Co-activation`="#D1495B"), name=NULL) +
  scale_x_continuous(expand=expansion(mult=c(.02,.08))) +
  labs(x="Transcriptomic region distance", y="Pairwise conservation AUROC",
       title="Conservation decays with region distance", subtitle=subt) +
  th + theme(plot.subtitle=element_text(size=BASefont-2, colour="grey30"),
             legend.position=c(.78,.88), aspect.ratio=.68,
             legend.background=element_rect(fill=alpha("white",.7),colour=NA))
sv(pd4, "d_distcurve", 3.50, 2.40)

## ---------- e-g: strict-hit network pairs and spatial fields ----------
ed <- fread(file.path(OUT,"fig_partner_network_edges.tsv"))
hits_tab <- fread(file.path(OUT,"m4_rewiring_hits.tsv"))
panel_plan <- fread(file.path(OUT,"fig_hit_panel_plan.tsv"))
stopifnot(nrow(panel_plan) == 3, nrow(hits_tab) > 0)
jacc_of <- setNames(hits_tab$mean_partner_jaccard, hits_tab$program)
TOPK <- 8
EDGE_BLUE <- "#7FB3D5"

# build a two-subnetwork data frame for one hit (radial layout per region, side-by-side)
mkpair_df <- function(h, ref, tgt) {
  get_top <- function(reg) {
    s <- ed[hit==h & region==reg & corr>0.05][order(-corr)]
    s[seq_len(min(TOPK,.N))]
  }
  sref <- get_top(ref); stgt <- get_top(tgt)
  pref <- sref$partner; ptgt <- stgt$partner
  shared <- intersect(pref, ptgt)
  # radial coords: hub at center, partners on a ring (deterministic order by corr)
  radial <- function(sub, cx) {
    n <- nrow(sub)
    ang <- seq(90, 90-360+360/n, length.out=n) * pi/180  # start top, clockwise
    data.table(partner=sub$partner, corr=sub$corr,
               px=cx + 1.0*cos(ang), py=1.0*sin(ang))
  }
  rr <- radial(sref, 0); rt <- radial(stgt, 3.2)   # two clusters, x-offset apart
  rr[, reg := ref]; rt[, reg := tgt]
  hubs <- data.table(reg=c(ref,tgt), hx=c(0,3.2), hy=c(0,0))
  list(
    nodes = rbind(
      cbind(rr, role=ifelse(rr$partner %in% shared,"shared","ref_only")),
      cbind(rt, role=ifelse(rt$partner %in% shared,"shared","tgt_only"))),
    hubs = hubs, shared = shared,
    edges = rbind(
      data.table(reg=ref, x0=0,   y0=0, x1=rr$px, y1=rr$py, w=rr$corr),
      data.table(reg=tgt, x0=3.2, y0=0, x1=rt$px, y1=rt$py, w=rt$corr)),
    nshare=length(shared), nref=length(pref), ntgt=length(ptgt))
}

mknet_pair <- function(h, ref, tgt) {
  hubname <- namemap[as.character(h)]
  hub_newP <- paste0("P", h)
  D <- mkpair_df(h, ref, tgt)
  nd <- D$nodes; eg <- D$edges; hb <- D$hubs
  turn <- 1 - D$nshare/((D$nref+D$ntgt)/2)   # mean turnover over the two regions
  rolecol <- c(shared="#7D3C98", ref_only="#2C6FB0", tgt_only="#D1495B")
  nd[, partner_label := paste0("P", partner)]
  ggplot() +
    geom_segment(data=eg, aes(x=x0,y=y0,xend=x1,yend=y1, linewidth=w, alpha=w),
                 colour=EDGE_BLUE, lineend="round") +
    geom_point(data=hb, aes(hx,hy), size=5.2, colour="#922B21") +
    geom_text(data=hb, aes(hx,hy,label=hub_newP), size=2.0, colour="white", fontface="bold") +
    geom_point(data=nd, aes(px,py, fill=role), shape=21, colour="white", stroke=.3, size=3.4) +
    geom_text(data=nd, aes(px,py,label=partner_label), size=1.6, colour="grey15") +
    geom_text(data=hb, aes(hx, -1.45, label=reg), size=2.4, fontface="bold", colour="grey10") +
    scale_fill_manual(values=rolecol, breaks=c("shared","ref_only","tgt_only"),
                      labels=c("shared","ref-specific","target-specific"), name=NULL) +
    scale_linewidth(range=c(.25,1.5), guide="none") + scale_alpha(range=c(.35,.95), guide="none") +
    coord_fixed(clip="off", xlim=c(-1.4,4.8), ylim=c(-1.7,1.4)) +
    labs(title=sprintf("%s\n(%s vs %s)", hubname, ref, tgt),
         subtitle=sprintf("shared %d/%d  |  Jaccard %.2f  |  turnover %.0f%%",
                          D$nshare, max(D$nref,D$ntgt), jacc_of[as.character(h)], 100*turn)) +
    theme_void(base_size=BASefont) +
    theme(plot.title=element_text(size=BASefont-1, face="bold", hjust=.5, lineheight=.9),
          plot.subtitle=element_text(size=BASefont-2, colour="grey35", hjust=.5),
          legend.position="bottom", legend.text=element_text(size=BASefont-2),
          legend.key.size=unit(2.6,"mm"), legend.box.spacing=unit(0,"mm"),
          plot.margin=margin(1,1,1,1))
}
for (i in seq_len(nrow(panel_plan))) {
  p <- panel_plan[i]
  sv(mknet_pair(p$program, p$reference, p$target), sprintf("e_netpair_%d", p$slot), 2.35, 1.68)
}

sf <- fread(file.path(OUT,"fig_hit_spatial_fields.tsv"))
mkfield <- function(reg, prog) {
  s <- sf[region==reg]
  col <- paste0("program_",prog)
  s[, val := get(col)]
  # robust scale 2-98 pct for color
  lim <- quantile(s$val, c(.02,.98), na.rm=TRUE)
  ggplot(s, aes(x, y, colour=pmin(pmax(val,lim[1]),lim[2]))) +
    geom_point(size=.25, shape=16) +
    scale_colour_viridis_c(option="magma", name=NULL) +
    coord_fixed() +
    # region label placed BELOW the field (caption) to avoid colliding with the
    # adjacent panel-d x-axis title in the abutting (gutter-0) row above.
    labs(caption=sprintf("%s", reg)) +
    theme_void(base_size=BASefont) +
    theme(plot.caption=element_text(size=BASefont, hjust=.5, face="bold",
                                    colour="grey10", margin=margin(t=1)),
          legend.position="none", plot.margin=margin(0,0,0,0))
}
for (i in seq_len(nrow(panel_plan))) {
  p <- panel_plan[i]
  sv(mkfield(p$field_region, p$program), sprintf("f_field_%d", p$slot), 1.55, 1.62)
}

## ---------- g: multi-step filter funnel ----------
fn <- fread(file.path(OUT,"m4_filter_funnel.tsv"))
stage_labels <- c(
  all_programs="All retained programs",
  gateA_expr_conserved="Gate A: expression conserved",
  `+gateB_neigh_low`="+ Gate B: neighborhood low",
  `+gateC_coact_spatial`="+ Gate C: spatial agreement",
  `+gateD_anchor_HITS`="+ Gate D: anchored hit"
)
fn[, stage := stage_labels[stage]]
fn[, stage := factor(stage, levels=rev(stage))]
fn[, lab := sprintf("%s (n=%d)", stage, n)]
pg <- ggplot(fn, aes(x=n, y=stage)) +
  geom_col(aes(fill=n), width=.65) +
  geom_text(aes(label=n), hjust=-.3, size=2.2) +
  scale_fill_gradient(low="#AED6F1", high="#1A5276", guide="none") +
  scale_x_continuous("Programs passing", expand=expansion(mult=c(0,.18))) +
  labs(y=NULL, title="Multi-step rewiring filter") +
  th + theme(axis.text.y=element_text(size=BASefont-1))
sv(pg, "g_funnel", 2.50, 1.62)

## ---------- spacer: thin gutter row between middle (d/network) and bottom (fields) ----------
## gutter-0 packing abuts d's x-axis title onto the bottom row; a thin near-blank
## spacer with a single faint hairline survives ink-crop and provides breathing room.
psp <- ggplot(data.frame(x=c(0,1), y=c(0,1))) +
  # two near-white hairlines top+bottom span full width AND set the ink height so
  # ink-crop keeps a band of real height (true blank would crop to zero). Invisible
  # on white paper; purpose = a ~3-4mm vertical gutter between middle and bottom rows.
  geom_segment(aes(x=0, xend=1, y=0,  yend=0 ), colour="#FEFEFE", linewidth=.1) +
  geom_segment(aes(x=0, xend=1, y=1,  yend=1 ), colour="#FEFEFE", linewidth=.1) +
  coord_cartesian(xlim=c(0,1), ylim=c(-.05,1.05), clip="off") +
  theme_void() + theme(plot.margin=margin(0,0,0,0))
# wide+short render; the y-extent between the two hairlines fixes the crop AR
# (~ width/height) so the spacer becomes a thin band (target ~3-4mm tall).
svglite(file.path(PANELS,"spacer.svg"), width=7.0, height=0.18, bg="white")
print(psp); invisible(dev.off())
cat("panel spacer -> 7.0 x 0.16 in\n")

cat("ALL PANELS DONE\n")

#!/usr/bin/env Rscript
# =============================================================================
# Figure B — Program x Program spatial co-organization (EDITOR PASS v2)
# Native single-script patchwork, shared theme, cairo_pdf vector.
# Panels a-e = COMMON (global), f-h = REGIONAL (per-area, 9 stable regions).
# Framing: HONEST. gene-similar co-orgs = real co-regulated modules (NOT dismissed);
# similarity used as annotation layer. NO "dense-tissue co-density floor" reasoning.
# =============================================================================
suppressPackageStartupMessages({
  library(reticulate); library(ggplot2); library(patchwork)
  library(ComplexHeatmap); library(circlize); library(grid)
  library(igraph); library(ggraph); library(ggridges)
  library(ggrepel); library(dendextend)
})
np <- import("numpy")

ROOT  <- "CORTEX_PROGRAM_ROOT"
MC    <- file.path(ROOT, "results/crossregion_v1/markcorr_v2")
FIN   <- file.path(MC, "final"); SIM <- file.path(MC, "similarity"); RIG <- file.path(MC, "rigor")
OUT   <- file.path(ROOT, "figures/markcorr_B")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- global theme: fonts >=5pt uniform ----
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
TAG <- theme(plot.tag = element_text(size = BASE + 3, face = "bold"),
             plot.tag.position = c(0.003, 0.992))

CLASS_COL <- c(exc="#D7263D", inh="#1B6CA8", glia="#2E933C",
               nonneuron="#8E44AD", vascular="#E08A00")
DIV <- colorRampPalette(c("#2166AC","#4393C3","#D1E5F0","#F7F7F7",
                          "#FDDBC7","#D6604D","#B2182B"))(101)

# =============================================================================
# LOAD
# =============================================================================
# NAME AUTHORITY (single source): program_names.tsv name_short.
# Star is re-judged from confidence (brain-weak => '*'); stale label_short/inherited stars NOT used.
lab <- read.delim(file.path(ROOT,"results/crossregion_v1/program_names.tsv"),
                   check.names=FALSE, stringsAsFactors=FALSE)        # program,name_full,name_short,confidence,...
ann <- read.delim(file.path(ROOT,"results/crossregion_v1/program_annotation.tsv"),
                  check.names=FALSE, stringsAsFactors=FALSE)         # P01.., class
ann$pid <- as.integer(sub("^P","",ann$program))
lab$pid <- as.integer(lab$program)
lab$name_short <- trimws(sub("\\*+$", "", lab$name_short))             # strip any inherited stale star
meta <- merge(lab[,c("pid","name_short","confidence")], ann[,c("pid","class")], by="pid")
meta <- meta[order(meta$pid),]
meta$short <- ifelse(meta$confidence=="brain-weak",
                     paste0(meta$name_short, "*"), meta$name_short)  # star = brain-weak per authority
meta$plab  <- paste0("P", meta$pid, " ", meta$short)
CLASS_ORD <- c("exc","inh","glia","nonneuron","vascular")
meta$class <- factor(meta$class, levels=CLASS_ORD)
N <- nrow(meta)                                                      # 60

# global npz (r=25um = ring index 0)
zg  <- np$load(file.path(FIN,"progprog_median_iqr.npz"), allow_pickle=TRUE)
L2  <- zg$f[["log2_median_g"]]                                       # 60x60x10
FSS <- zg$f[["frac_same_sign"]]
RING<- as.numeric(zg$f[["ring_edges_um"]])                          # 0,25,...,500
g25 <- L2[,,1]                                                       # r=25um slice (outer edge 25)
fss25 <- FSS[,,1]
rownames(g25) <- meta$pid; colnames(g25) <- meta$pid

# similarity (60x60, header 'Pn name')
gcos <- as.matrix(read.csv(file.path(SIM,"gene_cosine.csv"), row.names=1, check.names=FALSE))
apear<- as.matrix(read.csv(file.path(SIM,"activity_pearson.csv"), row.names=1, check.names=FALSE))
simp <- read.delim(file.path(SIM,"program_similarity_pairs.tsv"), stringsAsFactors=FALSE)
simp$pa <- as.integer(sub("program_","",simp$A)); simp$pb <- as.integer(sub("program_","",simp$B))

# rigor
loo <- read.delim(file.path(RIG,"progprog_donor_loo.tsv"), stringsAsFactors=FALSE)
jac <- read.delim(file.path(RIG,"progprog_jaccard.tsv"), stringsAsFactors=FALSE)

# byarea
zb  <- np$load(file.path(FIN,"progprog_byarea_median_iqr.npz"), allow_pickle=TRUE)
L2b <- zb$f[["log2_median_g"]]                                       # 14x60x60x10
AREA<- as.character(zb$f[["area_names"]])
UNST<- as.integer(zb$f[["unstable"]])
STABLE <- AREA[UNST==0]                                              # 9 stable

pid2lab <- setNames(meta$plab, meta$pid)
pid2short<- setNames(meta$short, meta$pid)
pid2cls <- setNames(as.character(meta$class), meta$pid)

# =============================================================================
# PANEL a — global 60x60 heatmap (r=25um), class-ordered, class color bars
# =============================================================================
ord  <- order(meta$class, meta$pid)
M    <- g25[ord, ord]
cls  <- as.character(meta$class)[ord]; pids <- meta$pid[ord]
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
ht <- Heatmap(M, name="log2 g\n(25 um)", col=colf, rect_gp=gpar(col=NA),
              cluster_rows=FALSE, cluster_columns=FALSE,
              top_annotation=cbar, left_annotation=rbar,
              show_row_names=FALSE, show_column_names=FALSE,
              width=unit(46,"mm"), height=unit(46,"mm"),
              column_title="Program x program spatial co-organization (r = 25 um)",
              column_title_gp=gpar(fontsize=BASE+1, fontface="bold"),
              heatmap_legend_param=list(title_gp=gpar(fontsize=BASE-0.5),
                labels_gp=gpar(fontsize=BASE-1), legend_height=unit(14,"mm"),
                grid_width=unit(3,"mm")))
panel_a <- grid.grabExpr(draw(ht, heatmap_legend_side="right",
                              annotation_legend_side="right", merge_legend=TRUE,
                              padding=unit(c(1,1,1,1),"mm")), wrap.grobs=TRUE)
panel_a <- wrap_elements(full=panel_a)

# =============================================================================
# PANEL b — co-org vs gene-similarity scatter (honest annotation)
# FIXES (P1-1,P1-3,P2-2): wrap subtitle to 2 lines; x-jitter+alpha spread red blob;
#   crop dead whitespace; recompute Pearson r (verified 0.78).
# =============================================================================
sb <- simp[simp$pa < simp$pb, ]
ij <- cbind(match(sb$pa, meta$pid), match(sb$pb, meta$pid))
sb$g  <- g25[ij]
sb$absg <- abs(sb$g)
rr <- cor(sb$gene_cosine, sb$g, use="complete.obs")                  # P2-2: recomputed, == 0.78
# gene-INDEPENDENT strong co-orgs: low cosine, high positive g
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
           label=sprintf("gene-independent\nstrong co-orgs: %d", nb_indep),
           hjust=0, vjust=0, size=1.9, colour="#B2182B", lineheight=0.9) +
  scale_x_continuous(expand=expansion(mult=c(0.01,0.03))) +
  scale_y_continuous(expand=expansion(mult=c(0.02,0.02))) +
  labs(x="Gene-loading cosine similarity", y="Spatial co-org  log2 g (25 um)",
       title="Co-organization rises with gene similarity",
       subtitle="strong co-orgs are mostly gene-similar co-regulated modules;\na distinct subset is gene-independent") +
  coord_cartesian(clip="off")

# =============================================================================
# PANEL c — gene-INDEPENDENT niche network (the microglia-synapse hub)  [P0 REBUILD]
# FIXES: stress layout + manual pin of hubs P40/P49/P60 to spread coords;
#   halo ring underlay on hubs; label only hubs (bold/large) + their top synaptic
#   partners with leader lines; edge alpha/width mapped to |log2 g|;
#   annotation note moved to margin caption (kept "P54 forms no hub" point, deg=1);
#   size legend kept HERE only (deduped from panel d). subtitle en-dash.
# =============================================================================
ed <- sb[sb$gene_cosine < 0.25 & abs(sb$g) > 0.32, c("pa","pb","g","gene_cosine")]
ed <- ed[ed$g > 0, ]                                                # positive co-orgs (co-occurrence niches)
deg <- table(c(ed$pa, ed$pb))
usei <- as.integer(names(deg))
vdf <- data.frame(pid=usei, deg=as.integer(deg),
                  class=pid2cls[as.character(usei)], short=pid2short[as.character(usei)])
vdf$absg_strength <- sapply(vdf$pid, function(p) sum(abs(ed$g[ed$pa==p | ed$pb==p])))
gph <- graph_from_data_frame(d=data.frame(from=ed$pa, to=ed$pb, w=ed$g, aw=abs(ed$g)),
                             vertices=vdf, directed=FALSE)
HUBS <- c(40,49,60)                                                 # microglial programs (nonneuron)
# top synaptic partners of the hubs to label (highest g, synaptic/neuronal)
hub_partner <- function(h, n=4){
  pr <- ed[ed$pa==h | ed$pb==h, ]; pr$partner <- ifelse(pr$pa==h, pr$pb, pr$pa)
  pr <- pr[order(-pr$g), ]; head(unique(pr$partner), n)
}
PART <- unique(unlist(lapply(HUBS, hub_partner, n=4)))
PART <- setdiff(PART, HUBS)

# ---- stress layout, then PIN hubs to well-separated coords ----
set.seed(11)
lay <- create_layout(gph, layout="stress")
xr <- range(lay$x); yr <- range(lay$y)
sx <- function(f) xr[1] + diff(xr)*f
sy <- function(f) yr[1] + diff(yr)*f
hub_xy <- list("40"=c(sx(0.18), sy(0.80)),
               "49"=c(sx(0.82), sy(0.78)),
               "60"=c(sx(0.50), sy(0.18)))
for(h in names(hub_xy)){
  idx <- which(as.character(lay$name) == h)
  if(length(idx)==1){ lay$x[idx] <- hub_xy[[h]][1]; lay$y[idx] <- hub_xy[[h]][2] }
}
# node attrs onto layout
lay$is_hub  <- as.integer(lay$name) %in% HUBS
lay$is_part <- as.integer(lay$name) %in% PART
lay$lab <- NA_character_
lay$lab[lay$is_hub]  <- paste0("P", lay$name[lay$is_hub])
lay$lab[lay$is_part] <- paste0("P", lay$name[lay$is_part])
lay$lab_bold <- lay$is_hub

# edge endpoints for manual segment draw (so we follow pinned coords)
ename <- as_edgelist(gph, names=TRUE)
exy <- data.frame(
  x   = lay$x[match(ename[,1], lay$name)], y   = lay$y[match(ename[,1], lay$name)],
  xend= lay$x[match(ename[,2], lay$name)], yend= lay$y[match(ename[,2], lay$name)],
  aw  = abs(E(gph)$w),
  hub = (as.integer(ename[,1]) %in% HUBS) | (as.integer(ename[,2]) %in% HUBS)
)
exy <- exy[order(exy$hub), ]                                         # hub edges drawn on top

panel_c <- ggplot() +
  # edges: width+alpha mapped to |log2 g|; hub edges pop
  geom_segment(data=subset(exy,!hub),
               aes(x=x,y=y,xend=xend,yend=yend, linewidth=aw, alpha=aw),
               colour="#C9A227", lineend="round") +
  geom_segment(data=subset(exy, hub),
               aes(x=x,y=y,xend=xend,yend=yend, linewidth=aw, alpha=aw),
               colour="#8E44AD", lineend="round") +
  scale_linewidth(range=c(0.10,0.85), guide="none") +
  scale_alpha(range=c(0.12,0.62), guide="none") +
  # HALO under hubs
  geom_point(data=subset(lay, is_hub), aes(x,y), size=10, colour="#8E44AD", alpha=0.14) +
  geom_point(data=subset(lay, is_hub), aes(x,y), size=7,  colour="#8E44AD", alpha=0.20) +
  # nodes
  geom_point(data=lay, aes(x,y, size=deg, fill=class), shape=21, colour="white", stroke=0.2) +
  scale_size(range=c(0.9,6.2), name="co-org degree",
             guide=guide_legend(override.aes=list(fill="grey50"))) +
  scale_fill_manual(values=CLASS_COL, name="class", guide="none") +
  # labels: hubs bold/larger, partners regular, leader lines, no silent drops
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
  labs(title="Gene-independent co-organization niches",
       subtitle="edges: gene cosine < 0.25 & |log2 g| > 0.32 – the microglia–synapse hub",
       caption="Microglial programs P40 / P49 / P60 (purple, pinned) form gene-independent hubs co-locating with synaptic programs (P17, P29, P47, P52);\nthe adjacent microglial complement/MHC program P54 forms no hub (degree 1).") +
  theme_void(base_size=BASE) +
  theme(plot.title=element_text(size=BASE+1,face="bold"),
        plot.subtitle=element_text(size=BASE-0.5,colour="grey40"),
        plot.caption=element_text(size=BASE-1.5, colour="grey35", hjust=0,
                                  lineheight=1.0, margin=margin(t=3)),
        legend.title=element_text(size=BASE-0.5), legend.text=element_text(size=BASE-1),
        legend.key.size=unit(3,"mm"), legend.position="right",
        plot.margin=margin(4,6,4,6))

# =============================================================================
# PANEL d — strongest AVOIDANCES (neuronal proton/Ca2+ vs oligo/myelin) lollipop
# FIXES (P1-1, P2-1): wrap subtitle to 2 lines; DROP duplicate size legend (n/a here);
#   panel d has no size legend so just subtitle wrap + keep fss legend.
# =============================================================================
neur <- c(8,24,57); olig <- c(36,26,37,13)
av <- sb[ (sb$pa %in% neur & sb$pb %in% olig) | (sb$pa %in% olig & sb$pb %in% neur), ]
av$fss <- fss25[cbind(match(av$pa,meta$pid), match(av$pb,meta$pid))]
av$pair <- paste0("P",av$pa," x P",av$pb)
av <- av[order(av$g), ]
av$pair <- factor(av$pair, levels=av$pair)
panel_d <- ggplot(av, aes(g, pair)) +
  geom_vline(xintercept=0, linewidth=0.2, colour="grey70") +
  geom_segment(aes(x=0, xend=g, yend=pair), colour="#2166AC", linewidth=0.4) +
  geom_point(aes(fill=fss), shape=21, size=2.0, colour="white", stroke=0.2) +
  scale_fill_gradient(low="#9ecae1", high="#08306b", name="frac\nsame-sign",
                      limits=c(0.85,1), breaks=c(0.9,1)) +
  labs(x="Spatial co-org  log2 g (25 um)", y=NULL,
       title="Strongest spatial avoidances",
       subtitle="neuronal proton/Ca2+ (P8,P24,P57) vs\noligodendrocyte/myelin (P13,P26,P36,P37)") +
  theme(axis.text.y=element_text(size=BASE-1.5))

# =============================================================================
# PANEL e — rigor: donor-LOO + jaccard + frac_same_sign  (GOOD; unchanged)
# =============================================================================
ut <- which(upper.tri(g25), arr.ind=TRUE)
hl <- data.frame(g=g25[ut], fss=fss25[ut])
hl <- hl[abs(hl$g) > 0.32, ]
e1 <- ggplot(hl, aes(fss)) +
  geom_histogram(bins=28, fill="#1B6CA8", colour="white", linewidth=0.1) +
  geom_vline(xintercept=median(hl$fss), linetype=2, linewidth=0.3, colour="#B2182B") +
  annotate("text", x=median(hl$fss), y=Inf, vjust=1.4, hjust=1.1, size=1.9, colour="#B2182B",
           label=sprintf("median %.2f", median(hl$fss))) +
  labs(x="frac same-sign across 44 chips", y="headline co-orgs",
       title="Reproducibility of headline co-orgs",
       subtitle="direction agreement across chips + donor-LOO stability") +
  coord_cartesian(clip="off")
loo$absfull <- abs(loo$full_log2)
lh <- loo[loo$absfull > 0.32, ]
lh <- lh[order(-lh$absfull), ]
lh <- head(lh, 22)
lh$pair <- paste0("P",lh$A_idx+1," x P",lh$B_idx+1)
lh$pair <- factor(lh$pair, levels=rev(lh$pair))
e2 <- ggplot(lh, aes(y=pair)) +
  geom_vline(xintercept=c(-0.32,0.32), linetype=3, linewidth=0.2, colour="grey70") +
  geom_segment(aes(x=log2(loo_min_median_g), xend=log2(loo_max_median_g), yend=pair),
               linewidth=0.5, colour="grey55") +
  geom_point(aes(x=full_log2, colour=stable_90), size=0.9) +
  scale_colour_manual(values=c("True"="#2E933C","False"="grey55"), name="LOO stable (90%)",
                      labels=c("True"="yes","False"="no")) +
  labs(x="log2 g  (point=full; bar=donor-LOO min-max)", y=NULL,
       title=NULL, subtitle=NULL) +
  theme(axis.text.y=element_text(size=BASE-2))
panel_e <- e1 + e2 + plot_layout(widths=c(1,1.15))

# =============================================================================
# PANEL f — region small-multiples  [P1-2: more left room + >=5pt y labels]
# =============================================================================
ex <- list(c(37,38), c(8,36), c(40,29), c(17,52))
exlab <- c("P37xP38 oligo/myelin", "P8xP36 neuron-oligo avoid", "P40xP29 microglia-glutR", "P17xP52 GABA-synapse")
ai <- match(STABLE, AREA)
frows <- do.call(rbind, lapply(seq_along(ex), function(k){
  pa<-ex[[k]][1]; pb<-ex[[k]][2]
  data.frame(area=STABLE, ex=exlab[k],
             g=sapply(ai, function(a) L2b[a, pa, pb, 1]))
}))
frows$ex <- factor(frows$ex, levels=exlab)
panel_f <- ggplot(frows, aes(reorder(area,g), g)) +
  geom_hline(yintercept=0, linewidth=0.2, colour="grey70") +
  geom_col(aes(fill=g>0), width=0.7, show.legend=FALSE) +
  scale_fill_manual(values=c("TRUE"="#B2182B","FALSE"="#2166AC")) +
  scale_y_continuous(expand=expansion(mult=c(0.04,0.06))) +
  facet_wrap(~ex, nrow=1) +
  labs(x=NULL, y="log2 g (25 um)",
       title="Exemplar co-orgs across 9 stable regions",
       subtitle="per-area median g; unstable regions (n<3 chips) excluded") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=BASE-1),
        axis.title.y=element_text(size=BASE),
        panel.spacing=unit(2,"mm"), strip.text=element_text(size=BASE-1),
        plot.margin=margin(4,4,4,8))

# =============================================================================
# PANEL g — region clustering: 9 stable regions x headline co-orgs (GOOD; unchanged)
# =============================================================================
hp <- data.frame(a=meta$pid[ut[,1]], b=meta$pid[ut[,2]], g=g25[ut])
hp <- hp[order(-abs(hp$g)), ]; hp <- head(hp, 40)
Mg <- sapply(seq_len(nrow(hp)), function(k){
  sapply(ai, function(a) L2b[a, hp$a[k], hp$b[k], 1])
})
rownames(Mg) <- STABLE
colnames(Mg) <- paste0("P",hp$a,"xP",hp$b)
SENS <- c("V1","S1","M1")
gcol <- ifelse(STABLE %in% SENS, "#E08A00", "#1B6CA8")
rng2 <- max(abs(Mg), na.rm=TRUE)
htg <- Heatmap(Mg, name="log2 g", col=colorRamp2(seq(-rng2,rng2,length.out=101), DIV),
               cluster_columns=TRUE, cluster_rows=TRUE,
               show_column_names=FALSE,
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
# PANEL h — cross-region span (GOOD; unchanged)
# =============================================================================
span <- data.frame(pair=colnames(Mg),
                   mn=apply(Mg,2,min,na.rm=TRUE), mx=apply(Mg,2,max,na.rm=TRUE))
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
# ASSEMBLE  (180mm wide x ~235mm portrait)
# =============================================================================
design <- "
AAB
AAB
CCD
EEE
FFF
GGH
GGH
"
fig <- (panel_a) + (panel_b) + (panel_c) + (panel_d) +
       (panel_e) + (panel_f) + (panel_g) + (panel_h) +
  plot_layout(design=design, heights=c(1,1,1.12,0.9,1,1,1)) +
  plot_annotation(tag_levels=list(c("a","b","c","d","e","f","g","h"))) &
  theme(plot.tag=element_text(size=BASE+3, face="bold"),
        plot.tag.position=c(0.004,0.985))

W_IN <- 180/25.4; H_IN <- 235/25.4
ggsave(file.path(OUT,"figB_full.pdf"), fig, width=W_IN, height=H_IN,
       device=cairo_pdf, bg="white")
ggsave(file.path(OUT,"figB_full.png"), fig, width=W_IN, height=H_IN,
       dpi=300, bg="white")
cat("DONE. Pearson r:", round(rr,3), "| indep co-orgs (panel b):", nb_indep,
    "| network nodes:", vcount(gph), "edges:", ecount(gph),
    "| hub deg P40/49/60:", deg["40"], deg["49"], deg["60"],
    "| stable regions:", length(STABLE), "\n")

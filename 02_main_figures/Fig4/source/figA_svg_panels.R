#!/usr/bin/env Rscript
# =====================================================================
# Figure A LARGE : cell-type x program spatial co-localization
# Integrates the 8 abstract panels (a-h) + spatial tissue exemplars
# (HERO 3chip x 3field grid + co-loc/avoidance ZOOM) into ONE native
# patchwork. cairo_pdf vector. ~180mm wide portrait.
#
# Spatial panels ported to ggplot geom_point (rule 8: native patchwork,
# unified theme/fonts, NOT png montage). Data pre-exported per-chip by
# scripts/figmarkcorr_spatial/export_spatial_for_R.py.
# Chip swap: DLPFC B01012B2 (oligo-dead) -> SPL B02221E6 (oligo-rich).
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ComplexHeatmap)
  library(circlize); library(grid); library(ggrepel); library(dendextend)
  library(data.table); library(scales); library(viridisLite)
})

ROOT <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
FIN  <- file.path(ROOT, "markcorr_v2/final")
RIG  <- file.path(ROOT, "markcorr_v2/rigor")
SPDIR<- "CORTEX_PROGRAM_ROOT/scripts/figmarkcorr_spatial/export_R"
OUT  <- "CORTEX_PROGRAM_ROOT/figures/markcorr_A"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

BASE_FONT <- "sans"
FS <- 6
TAG_SZ <- 9
HEADLINE_THR <- 0.32

shorten_lab <- function(x){
  x <- gsub(" x ", " · ", x, fixed=TRUE)
  x <- gsub("Oligodendrocyte", "Oligo", x)
  x <- gsub("OLIGO", "Oligo", x)
  x <- gsub("CHANDELIER", "Chand.", x)
  x <- gsub("LINC00507", "", x)
  x <- gsub("Glutamate-receptor signaling", "Glu-receptor sig.", x, fixed=TRUE)
  x <- gsub("Glutamate receptor signaling", "Glu-receptor sig.", x, fixed=TRUE)
  x <- gsub("Cation channel (interneuron)", "Cation channel (IN)", x, fixed=TRUE)
  x <- gsub("\\s+", " ", x); trimws(x)
}

pairs <- fread(file.path(FIN, "cellprog_pairs_median_iqr.tsv"))
barea <- fread(file.path(FIN, "cellprog_byarea_median_iqr.tsv"))
logo  <- fread(file.path(RIG, "cellprog_logo.tsv"))
jacc  <- fread(file.path(RIG, "cellprog_jaccard.tsv"))
nm    <- fread(file.path(ROOT, "program_names.tsv"))

nm    <- nm[new_P != "EXCLUDED"]
nm[, pid := cnmf_component]
nm[, new_pid := as.integer(sub("^P", "", new_P))]
nm[, weak := (confidence == "brain-weak")]
nm[, plab := sprintf("P%d %s%s", new_pid, name_short, ifelse(weak, "*", ""))]
plab_map <- setNames(nm$plab, paste0("program_", nm$pid))
pshort   <- setNames(nm$name_short, paste0("program_", nm$pid))
pstar    <- setNames(ifelse(nm$weak, "*", ""), paste0("program_", nm$pid))
pshort_s <- function(B) paste0(pshort[B], pstar[B])  # short name + trailing * if weak

ct_class <- c(
  AST="Astro", OLIGO="Oligo", OPC="Oligo", MICRO="Immune", ENDO="Vascular",
  VLMC="Vascular", `L2-L3 IT LINC00507`="Exc-IT", `L3-L4 IT RORB`="Exc-IT",
  `L4-L5 IT RORB`="Exc-IT", `L6 IT`="Exc-IT", `L6 IT CAR3`="Exc-IT",
  ET="Exc-deep", NP="Exc-deep", `L6 CT`="Exc-deep", L6B="Exc-deep",
  `L6 CAR3`="Exc-deep", PVALB="Inh-MGE", SST="Inh-MGE", CHANDELIER="Inh-MGE",
  LAMP5="Inh-CGE", VIP="Inh-CGE", NDNF="Inh-CGE", PAX6="Inh-CGE")

r25 <- pairs[ring_um == 25]
r25[, plab := plab_map[B]]
r25[, ct := A]

# Between-chip null-model input:
# Use between-chip Stouffer q (318 headline pairs, q<0.05 and |log2g|>0.32)
# to define headline_progs while displaying the current 54-program namespace.
BTWNDIR <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1/markcorr_betweenchip_v1"
sq <- fread(file.path(BTWNDIR, "betweenchip_cellprog_stouffer_q.tsv"))
sq_hl <- sq[is_headline == TRUE]  # 318 headline pairs
surv_keys_all <- paste(logo[survives==TRUE]$A, logo[survives==TRUE]$B)
r25[, sk := paste(A,B)]
.posB <- sq_hl[median_log2g>0][order(-median_log2g)][1:min(12,.N)]$B_name
.negB <- sq_hl[median_log2g<0][order(median_log2g)][1:min(5,.N)]$B_name
.lgB  <- {l<-copy(logo)[!is.na(orig_log2g)&is.finite(logo_log2g)][abs(orig_log2g)>=HEADLINE_THR][order(-abs(orig_log2g))][1:min(14,.N)]; l$B}
.spB  <- {s<-barea[ring_um==25 & area %in% c("AG","DLPFC","FPPFC","M1","S1","SMG","SPL","V1","VLPFC")][,.(lo=min(log2_median_g,na.rm=TRUE),hi=max(log2_median_g,na.rm=TRUE),n=.N),by=.(A,B)][n>=7]; s[,span:=hi-lo]; s[order(-span)][1:10]$B}
headline_progs <- unique(c(.posB,.negB,.lgB,.spB))
headline_progs <- headline_progs[!is.na(headline_progs)]
sq_hl_keys <- paste(sq_hl$A_name, sq_hl$B_name)  # "CELLTYPE program_N"

th <- theme_minimal(base_size = FS, base_family = BASE_FONT) +
  theme(
    plot.tag = element_text(face="bold", size=TAG_SZ),
    plot.title = element_text(size=FS+1, face="bold", hjust=0),
    plot.subtitle = element_text(size=FS-0.5, color="grey35"),
    axis.title = element_text(size=FS),
    axis.text  = element_text(size=FS-0.5, color="grey20"),
    legend.title = element_text(size=FS-0.5),
    legend.text  = element_text(size=FS-1),
    legend.key.size = unit(3.2, "mm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth=0.2, color="grey90"),
    plot.margin = margin(3,3,3,3)
  )
theme_set(th)

pal_div <- colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F7F7F7",
                              "#F4A582","#D6604D","#B2182B"))(101)
col_pos <- "#B2182B"; col_neg <- "#2166AC"; col_grey <- "grey75"

# =====================================================================
# Panel a : master heatmap
# =====================================================================
M <- dcast(r25, A ~ B, value.var="log2_median_g")
rn <- M$A; M$A <- NULL
Mm <- as.matrix(M); rownames(Mm) <- rn
COHORT_TECH_A <- paste0("program_", c(9,18,19,35,52,57))
ord_cols <- paste0("program_", 1:60)
ord_cols <- ord_cols[ord_cols %in% colnames(Mm) & !(ord_cols %in% COHORT_TECH_A)]
Mm <- Mm[, ord_cols]
col_pid <- as.integer(sub("program_", "", ord_cols))
col_wk  <- ifelse(nm$weak[match(col_pid, nm$pid)], "*", "")
col_new_pid <- nm$new_pid[match(col_pid, nm$pid)]
col_lab <- ifelse(ord_cols %in% headline_progs,
                  sprintf("P%d%s", col_new_pid, col_wk), sprintf("%d%s", col_new_pid, col_wk))
col_fontsz <- ifelse(ord_cols %in% headline_progs, 5.2, 4.2)
colnames(Mm) <- ord_cols
cls <- ifelse(rownames(Mm) %in% names(ct_class), ct_class[rownames(Mm)], "Other")
cls_col <- c(Astro="#1B9E77", Oligo="#7570B3", Immune="#E7298A",
             Vascular="#666666", `Exc-IT`="#D95F02", `Exc-deep`="#A6761D",
             `Inh-MGE`="#1F78B4", `Inh-CGE`="#66A61E", Other="grey80")
# force factor so every present class is enumerated in the legend (bar<->legend match)
cls <- factor(cls, levels=names(cls_col)[names(cls_col) %in% cls])
mx <- max(abs(Mm), na.rm=TRUE); mx <- min(mx, 1.6)
ha <- rowAnnotation(Class=cls, col=list(Class=cls_col),
        annotation_name_gp=gpar(fontsize=FS-0.5),
        annotation_legend_param=list(Class=list(title_gp=gpar(fontsize=FS-0.5),
          labels_gp=gpar(fontsize=FS-1), grid_height=unit(3,"mm"), grid_width=unit(3,"mm"))))
ht <- Heatmap(Mm, name="log2g",
        col=colorRamp2(seq(-mx,mx,length.out=101), pal_div),
        cluster_rows=TRUE, cluster_columns=TRUE,
        show_row_dend=TRUE, show_column_dend=TRUE,
        row_dend_width=unit(5,"mm"), column_dend_height=unit(4,"mm"),
        row_names_gp=gpar(fontsize=5.2),
        column_labels=col_lab,
        column_names_gp=gpar(fontsize=col_fontsz),
        column_names_rot=90, column_names_side="bottom",
        column_dend_side="top",
        left_annotation=ha,
        column_names_max_height=unit(12,"mm"),
        row_names_max_width=unit(28,"mm"),
        rect_gp=gpar(col="white", lwd=0.15),
        width=unit(54*1.85,"mm"), height=unit(42,"mm"),  # 54 biological programs x 1.85mm; 22 rows = ~1.9mm/row
        heatmap_legend_param=list(title="log2 g", title_gp=gpar(fontsize=FS-0.5),
          labels_gp=gpar(fontsize=FS-1), legend_height=unit(14,"mm"),
          grid_width=unit(3,"mm")))
pa <- grid.grabExpr(draw(ht, heatmap_legend_side="right",
        annotation_legend_side="right", merge_legend=TRUE,
        padding=unit(c(6,2,5,1),"mm")))  # b,l,t,r : bottom pad keeps col labels inside grab; top pad clears figure subtitle
pa <- wrap_elements(full = pa)

# =====================================================================
# Panel b : headline niche forest
# =====================================================================
r25[, surv := { key <- paste(A,B); s <- setNames(logo$survives, paste(logo$A,logo$B));
                ifelse(key %in% names(s), s[key], NA) }]
hl <- r25[abs(log2_median_g) >= HEADLINE_THR]
pos <- hl[log2_median_g > 0][order(-log2_median_g)]
surv_keys <- paste(logo[survives==TRUE]$A, logo[survives==TRUE]$B)
pos[, sk := paste(A,B)]
pos_s <- pos[sk %in% surv_keys][1:min(8,.N)]
micro_niche_force <- c("VIP program_33", "L6 IT program_54")
pos_add <- pos[sk %in% micro_niche_force & !(sk %in% pos_s$sk)]
pos_s <- rbind(pos_s, pos_add, fill=TRUE)
neg <- hl[log2_median_g < 0][order(log2_median_g)][1:min(3,.N)]
fb <- rbind(pos_s, neg, fill=TRUE)
fb[, lab := shorten_lab(sprintf("%s x %s", ct, plab))]
fb[, dir := ifelse(log2_median_g>0,"co-localize","avoid")]
fb[, log2q1 := log2(q1_g)]; fb[, log2q3 := log2(q3_g)]
fb <- fb[order(log2_median_g)]
fb[, lab := make.unique(lab)]
fb[, lab := factor(lab, levels=lab)]
pb <- ggplot(fb, aes(log2_median_g, lab, color=dir)) +
  geom_vline(xintercept=0, linewidth=0.3, color="grey60") +
  geom_errorbarh(aes(xmin=log2q1, xmax=log2q3), height=0, linewidth=0.5, alpha=0.6) +
  geom_point(size=1.4) +
  scale_color_manual(values=c(`co-localize`=col_pos, avoid=col_neg), name=NULL) +
  labs(x=expression(log[2]~median~g~"(r=25"*mu*"m)"), y=NULL,
       title="Headline niches (LOGO-surviving)") +
  theme(axis.text.y=element_text(size=5.2, hjust=1), legend.position="bottom",
        legend.margin=margin(0,0,0,0), plot.margin=margin(3,3,3,16))

# =====================================================================
# Panel c : g(r) curves
# =====================================================================
ex <- list(c("OLIGO","program_45"), c("OLIGO","program_37"),
           c("PVALB","program_7"), c("CHANDELIER","program_50"),
           c("SST","program_20"), c("L2-L3 IT LINC00507","program_45"))
exdt <- rbindlist(lapply(ex, function(e){
  d <- pairs[A==e[1] & B==e[2]]; if(nrow(d)==0) return(NULL)
  d[, lab := shorten_lab(sprintf("%s x %s", A, pshort_s(B)))]; d }))
broad_ct <- c("OLIGO"); exdt[, arch := ifelse(A %in% broad_ct, "domain-scale (broad)", "micro-niche (single-bin)")]
exdt[, log2g := log2_median_g]
exdt[, lo := log2(q1_g)]; exdt[, hi := log2(q3_g)]
c_series <- unique(exdt$lab)
c_lty <- setNames(rep(c("solid","longdash","dotted"), length.out=length(c_series)), c_series)
c_shp <- setNames(rep(c(16,17,15,18), length.out=length(c_series)), c_series)
pc <- ggplot(exdt, aes(ring_um, log2g, color=lab, fill=lab)) +
  geom_hline(yintercept=0, linewidth=0.3, color="grey70") +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.10, color=NA) +
  geom_line(aes(linetype=lab), linewidth=0.6) +
  geom_point(aes(shape=lab), size=1.0) +
  facet_wrap(~arch, ncol=1) +
  scale_color_brewer(palette="Dark2", name=NULL) +
  scale_fill_brewer(palette="Dark2", guide="none") +
  scale_linetype_manual(values=c_lty, name=NULL) +
  scale_shape_manual(values=c_shp, name=NULL) +
  labs(x=expression(distance~r~"("*mu*"m)"), y=expression(log[2]~median~g(r)),
       title="Distance decay: two archetypes") +
  theme(legend.position="bottom", legend.text=element_text(size=5),
        strip.text=element_text(size=FS-0.5, face="bold"),
        legend.key.size=unit(3.4,"mm"), legend.box.margin=margin(2,0,0,0),
        legend.spacing.y=unit(0.5,"mm")) +
  guides(color=guide_legend(ncol=2, byrow=TRUE),
         linetype=guide_legend(ncol=2, byrow=TRUE),
         shape=guide_legend(ncol=2, byrow=TRUE))

# =====================================================================
# Panel d : volcano
# =====================================================================
vd <- copy(r25)
vd[, tier := ifelse(abs(log2_median_g)>=HEADLINE_THR,
              ifelse(log2_median_g>0,"headline +","headline -"),"sub-threshold")]
vd[, lab := shorten_lab(sprintf("%s x %s", ct, pshort_s(B)))]
set.seed(7)
vd[, yj := frac_same_sign + runif(.N, -0.012, 0.012)]
toplab <- vd[order(-abs(log2_median_g))][1:5]
dxr <- max(abs(vd$log2_median_g), na.rm=TRUE)*1.05
pd <- ggplot(vd, aes(log2_median_g, yj, color=tier)) +
  geom_vline(xintercept=c(-HEADLINE_THR,HEADLINE_THR), linetype=2, linewidth=0.25, color="grey60") +
  geom_point(size=0.7, alpha=0.7) +
  geom_text_repel(data=toplab, aes(label=lab), size=1.8, max.overlaps=Inf,
                  segment.size=0.2, min.segment.length=0, color="grey20",
                  box.padding=0.6, point.padding=0.3, force=6, force_pull=0.5,
                  seed=7) +
  scale_color_manual(values=c(`headline +`=col_pos, `headline -`=col_neg,
                              `sub-threshold`="grey80"), name=NULL) +
  scale_x_continuous(limits=c(-dxr, dxr)) +
  labs(x=expression(log[2]~median~g~"(r=25"*mu*"m)"), y="fraction same-sign chips",
       title="Effect x consistency") +
  theme(legend.position="bottom", legend.key.size=unit(2.6,"mm"))

# =====================================================================
# Panel e : LOGO dumbbell
# =====================================================================
lg <- copy(logo)
lg[, plab := plab_map[B]]
lg[, lab := shorten_lab(sprintf("%s x %s", A, pshort_s(B)))]
lg <- lg[!is.na(orig_log2g) & is.finite(logo_log2g)]
lge <- lg[order(-abs(orig_log2g))][abs(orig_log2g)>=HEADLINE_THR][1:min(10,.N)]
lge[, lab := make.unique(lab)]
lge <- lge[order(orig_log2g)]
lge[, lab := factor(lab, levels=lab)]
lge[, surv := ifelse(survives,"survives","weakened")]
pe <- ggplot(lge) +
  geom_segment(aes(x=orig_log2g, xend=logo_log2g, y=lab, yend=lab),
               linewidth=0.5, color="grey70") +
  geom_point(aes(orig_log2g, lab, shape="full"), size=1.5, color="grey25") +
  geom_point(aes(logo_log2g, lab, color=surv, shape="LOGO"), size=1.5) +
  scale_color_manual(values=c(survives=col_pos, weakened="#E08214"), name=NULL) +
  scale_shape_manual(values=c(full=16, LOGO=17), name=NULL) +
  geom_text(aes(x=pmax(orig_log2g,logo_log2g)+0.05, y=lab,
                label=sprintf("J=%.2f", jaccard)), size=1.6, hjust=0, color="grey40") +
  labs(x=expression(log[2]~median~g), y=NULL, title="LOGO gene-removal retention") +
  theme(axis.text.y=element_text(size=5.2, hjust=1), legend.position="bottom",
        legend.key.size=unit(2.6,"mm"), plot.margin=margin(3,3,3,16)) +
  expand_limits(x=max(lge$orig_log2g)+0.35)

# =====================================================================
# Panel f : region small-multiples
# =====================================================================
stable <- c("AG","DLPFC","FPPFC","M1","S1","SMG","SPL","V1","VLPFC")
unstable <- c("ACC","ITG","PoCG","S1E","STG")
exf <- list(c("OLIGO","program_37"), c("PVALB","program_7"),
            c("AST","program_15"), c("SST","program_20"))
fdt <- rbindlist(lapply(exf, function(e){
  d <- barea[A==e[1] & B==e[2] & ring_um==25]; if(nrow(d)==0) return(NULL)
  d[, lab := shorten_lab(sprintf("%s x %s", A, pshort_s(B)))]; d}))
fdt[, region := factor(area, levels=unique(c(stable,unstable)))]
fdt[, stab := ifelse(area %in% stable, "stable","unstable")]
pf <- ggplot(fdt, aes(region, log2_median_g, fill=stab)) +
  geom_hline(yintercept=0, linewidth=0.3, color="grey70") +
  geom_col(width=0.7) +
  facet_wrap(~lab, ncol=2, scales="free_y",
             labeller=label_wrap_gen(width=20)) +
  scale_fill_manual(values=c(stable="#3182BD", unstable="grey80"), name=NULL) +
  labs(x=NULL, y=expression(log[2]~median~g), title="Region variation (exemplars)") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=FS-2),
        strip.text=element_text(size=FS-1.5, face="bold", lineheight=0.9,
                                margin=margin(1.5,1,1.5,1)),
        legend.position="bottom", legend.key.size=unit(2.6,"mm"),
        plot.margin=margin(3,3,3,6))

# =====================================================================
# Panel g : region clustering heatmap
# =====================================================================
hn <- pos_s[1:min(10,nrow(pos_s))]
hn_keys <- paste(hn$A, hn$B)
gdt <- barea[ring_um==25 & area %in% stable & paste(A,B) %in% hn_keys]
gdt[, niche := shorten_lab(sprintf("%s x %s", A, plab_map[B]))]
Gm <- dcast(gdt, niche ~ area, value.var="log2_median_g")
gr <- Gm$niche; Gm$niche <- NULL; Gm <- as.matrix(Gm); rownames(Gm) <- gr
gmx <- max(abs(Gm), na.rm=TRUE)
htg <- Heatmap(Gm, name="log2g_g",
        col=colorRamp2(seq(-gmx,gmx,length.out=101), pal_div),
        cluster_rows=TRUE, cluster_columns=TRUE,
        show_row_dend=FALSE, column_dend_height=unit(4,"mm"),
        row_names_gp=gpar(fontsize=5.2), column_names_gp=gpar(fontsize=5.2),
        column_names_rot=45,
        rect_gp=gpar(col="white", lwd=0.3),
        width=unit(ncol(Gm)*5,"mm"), height=unit(nrow(Gm)*5,"mm"),
        column_title="region clustering (stable)", column_title_gp=gpar(fontsize=FS),
        heatmap_legend_param=list(title="log2 g", title_gp=gpar(fontsize=FS-0.5),
          labels_gp=gpar(fontsize=FS-1), legend_height=unit(10,"mm"), grid_width=unit(3,"mm")))
pg <- grid.grabExpr(draw(htg, heatmap_legend_side="right", padding=unit(c(6,4,2,2),"mm")))
pg <- wrap_elements(full=pg)

# =====================================================================
# Panel h : cross-region span dumbbell
# =====================================================================
sp <- barea[ring_um==25 & area %in% stable]
sp <- sp[, .(lo=min(log2_median_g,na.rm=TRUE), hi=max(log2_median_g,na.rm=TRUE),
             md=median(log2_median_g,na.rm=TRUE), n=.N), by=.(A,B)]
sp <- sp[n>=7]
sp[, span := hi-lo]
sp[, lab := shorten_lab(sprintf("%s x %s", A, pshort_s(B)))]
sph <- sp[order(-span)][1:10]
sph[, lab := make.unique(lab)]
sph <- sph[order(span)]
sph[, lab := factor(lab, levels=lab)]
ph <- ggplot(sph) +
  geom_segment(aes(x=lo, xend=hi, y=lab, yend=lab, color=span), linewidth=0.9) +
  geom_point(aes(lo, lab), size=1.1, color="grey40") +
  geom_point(aes(hi, lab), size=1.1, color="grey20") +
  geom_point(aes(md, lab), size=0.9, shape=18, color="white") +
  scale_color_viridis_c(option="magma", direction=-1, name="span") +
  labs(x=expression(log[2]~median~g~"(min"%->%"max across stable regions)"),
       y=NULL, title="Most region-modulated niches") +
  theme(axis.text.y=element_text(size=5.2), legend.position="right",
        legend.key.size=unit(3,"mm"), plot.margin=margin(3,3,3,5))

# =====================================================================
# SPATIAL panels (ggplot port) -- HERO grid + ZOOM
# =====================================================================
# RASTERIZE dense per-bin point layers: 140k-410k square markers blow past
# rsvg's 1e6 XML-element cap during ink-crop. Wrap geom_point -> ggrastr
# rasterise (400dpi) for ALL subsequent (spatial) layers ONLY. Abstract
# panels a-l above already built with vector geom_point (unaffected).
suppressPackageStartupMessages(library(ggrastr))
geom_point <- function(...) ggrastr::rasterise(ggplot2::geom_point(...), dpi=400, dev="ragg")
nmv <- fread(file.path(SPDIR,"names.tsv"))
nmm <- setNames(nmv$value, nmv$key)
myelin_full <- nmm[["myelin_full"]]
zoom_chip_lab <- nmm[["zoom_chip"]]

hero <- fread(file.path(SPDIR,"hero_bins.csv.gz"))
chip_order <- c("B02221E6","D00865B3","A01186A4")  # SPL / V1 / M1
reg_levels <- c("SPL (parietal)","V1 (occipital)","M1 (motor)")
hero[, region_label := factor(region_label, levels=reg_levels)]

PT_SZ <- 0.18   # tissue dot size (square markers, dense)

# base spatial theme (no axes, equal aspect)
sp_theme <- theme_void(base_size=FS, base_family=BASE_FONT) +
  theme(plot.title=element_text(size=FS-0.5, hjust=0.5, margin=margin(b=1), lineheight=0.9),
        strip.text=element_text(size=FS-0.5, face="bold"),
        legend.position="bottom",
        legend.title=element_text(size=FS-1),
        legend.text=element_text(size=FS-1),
        legend.key.height=unit(2.4,"mm"), legend.key.width=unit(7,"mm"),
        legend.box.spacing=unit(1,"mm"),
        panel.spacing=unit(0.3,"mm"),
        plot.margin=margin(1,0.5,1,0.5))

# helper: one spatial column ggplot, faceted by chip (3 rows), shared scale via _s
spatial_col <- function(dat, sval, low_col, high_col, title, lab, viridis=FALSE){
  dv <- dat[valid==TRUE]
  di <- dat[valid==FALSE]
  g <- ggplot() +
    geom_point(data=di, aes(x,y), color="#dcdcdc", size=PT_SZ, shape=15) +
    geom_point(data=dv, aes(x,y, color=.data[[sval]]), size=PT_SZ, shape=15)
  if(viridis){
    g <- g + scale_color_viridis_c(option="viridis", limits=c(0,1), name=lab,
              breaks=c(0,0.5,1), labels=c("lo","","hi"))
  } else {
    g <- g + scale_color_gradient(low="#fff5f0", high=high_col, limits=c(0,1), name=lab,
              breaks=c(0,0.5,1), labels=c("lo","","hi"))
  }
  g + facet_grid(region_label ~ ., switch="y") +
    coord_fixed() + scale_y_reverse() +
    guides(color=guide_colorbar(title.position="left", title.vjust=1)) +
    labs(title=title) + sp_theme +
    theme(strip.placement="outside",
          strip.clip="off",
          panel.spacing.y=unit(0.6,"mm"),
          strip.text.y.left=element_text(angle=90, size=FS-1, face="bold",
                                         margin=margin(0,1.2,0,0.5)),
          plot.margin=margin(1,0.5,1,0.5))
}

# column titles (functional names)
sc_oligo  <- spatial_col(hero, "oligo_s",  NULL, "#cb181d", "OLIGO weight",        "rel.\nweight", viridis=FALSE)
sc_myelin <- spatial_col(hero, "myelin_s", NULL, NULL,      "Myelination (P41)",   "rel.\nSCT",    viridis=TRUE)
sc_it     <- spatial_col(hero, "it_s",     NULL, "#cb181d", "L2-L3 IT weight",     "rel.\nweight", viridis=FALSE)
# only leftmost keeps the region strip labels
sc_myelin <- sc_myelin + theme(strip.text.y.left=element_blank())
sc_it     <- sc_it     + theme(strip.text.y.left=element_blank())

# leftmost col carries the region strip (extra width); give it a touch more so the
# three MAP areas come out equal-sized & square. guides collected to one bottom row.
hero_block <- (sc_oligo | sc_myelin | sc_it) +
  plot_layout(guides="collect", widths=c(1.06,1,1)) +
  plot_annotation(theme=theme(plot.margin=margin(0,0,0,0))) &
  theme(legend.position="bottom")
hero_block <- wrap_elements(full = hero_block)

# ---- ZOOM : co-localization vs avoidance (V1) ----
zoom <- fread(file.path(SPDIR,"zoom_bins.csv.gz"))
zhi  <- zoom[myelin_hi==TRUE]
zoom_panel <- function(dat, sval, high_col, title){
  dv <- dat[valid==TRUE]; di <- dat[valid==FALSE]
  ggplot() +
    geom_point(data=di, aes(x,y), color="#f0f0f0", size=PT_SZ, shape=15) +
    geom_point(data=zhi, aes(x,y), color="#6baed6", size=PT_SZ*0.85, shape=15, alpha=0.40) +
    geom_point(data=dv[get(sval)>0.45], aes(x,y, color=.data[[sval]]), size=PT_SZ, shape=15) +
    scale_color_gradient(low="#fcae91", high=high_col, limits=c(0.45,1), name="rel. weight",
              breaks=c(0.5,1), labels=c("mid","hi"), oob=scales::squish) +
    coord_fixed() + scale_y_reverse() + labs(title=title) + sp_theme +
    theme(plot.title=element_text(size=FS-0.5, hjust=0.5, lineheight=0.9, face="bold"))
}
z_oligo <- zoom_panel(zoom, "oligo_s", "#cb181d",
                      "OLIGO ⊕ Myelin\n(co-localize)")
z_it    <- zoom_panel(zoom, "it_s", "#cb181d",
                      "L2-L3 IT ⊘ Myelin\n(avoid)")
# explanatory caption + two maps; pad with spacers so maps spread, titles never collide
# short single-line caption; detail lives in the two panel titles + shared legend
# short caption (V1 only, no double-paren); narrow column so the two maps fill width
zcap <- ggplot()+theme_void()+
  annotate("text",x=0.5,y=0.5,hjust=0.5,vjust=0.5,size=1.75,fontface="bold",
           lineheight=1.2,
           label="Co-localization vs\nmutual exclusion (V1)\nblue = Myelin high-zone")+
  xlim(0,1)+ylim(0,1)
zoom_block <- (zcap | z_oligo | z_it) +
  plot_layout(widths=c(0.55,1,1), guides="collect") &
  theme(legend.position="bottom",
        legend.title=element_text(size=FS-1), legend.text=element_text(size=FS-1),
        legend.key.height=unit(2.4,"mm"), legend.key.width=unit(7,"mm"))
zoom_block <- wrap_elements(full = zoom_block)

# =====================================================================
# NEW AVOIDANCE 1 : P47 Reactive-astro/vascular  (-)  upper-layer IT
#   region-contrast, COMPACT to 2 rows (M1 strong | V1 absent), 2 cols each:
#   [P47 SCT field]  [upper-IT weight field + P47 high-zone overlay (blue)]
#   In M1 the P47 high-zone sits where IT is empty (exclusion); in V1 they overlap.
# =====================================================================
avn  <- fread(file.path(SPDIR,"avoid_names.tsv"))
avnm <- setNames(avn$value, avn$key)
p51_short <- avnm[["p51_short"]]      # "Reactive astrocyte/vascular"
p8_short  <- avnm[["p8_short"]]       # "Mid-deep IT neuropil/cytoskeleton" (relabeled, from program_names.tsv)
av2_r     <- avnm[["avoid2_r"]]

av1 <- fread(file.path(SPDIR,"avoid1_bins.csv.gz"))
# shorten region strip labels (long labels garble the rotated y-strip)
av1[, region_label := ifelse(grepl("^M1", region_label), "M1 (strong)", "V1 (absent)")]
a1_levels <- c("M1 (strong)","V1 (absent)")
av1[, region_label := factor(region_label, levels=a1_levels)]

# one faceted (2-region) column for avoid-1
av1_col <- function(dat, sval, low_col, high_col, title, lab, viridis=FALSE, overlay=FALSE){
  dv <- dat[valid==TRUE]; di <- dat[valid==FALSE]; dhi <- dat[p51_hi==TRUE]
  g <- ggplot() +
    geom_point(data=di, aes(x,y), color="#dcdcdc", size=PT_SZ, shape=15) +
    geom_point(data=dv, aes(x,y, color=.data[[sval]]), size=PT_SZ, shape=15)
  if(overlay) g <- g + geom_point(data=dhi, aes(x,y), color="#08519c",
                                  size=PT_SZ*0.75, shape=15, alpha=0.42)
  if(viridis){
    g <- g + scale_color_viridis_c(option="viridis", limits=c(0,1), name=lab,
              breaks=c(0,0.5,1), labels=c("lo","","hi"))
  } else {
    g <- g + scale_color_gradient(low="#fff5f0", high=high_col, limits=c(0,1), name=lab,
              breaks=c(0,0.5,1), labels=c("lo","","hi"))
  }
  g + facet_grid(region_label ~ ., switch="y") +
    coord_fixed() + scale_y_reverse() +
    guides(color=guide_colorbar(title.position="left", title.vjust=1)) +
    labs(title=title) + sp_theme +
    theme(strip.placement="outside", strip.clip="off",
          panel.spacing.y=unit(0.6,"mm"),
          strip.text.y.left=element_text(angle=90, size=FS-1, face="bold",
                                         margin=margin(0,1.2,0,0.5)))
}
a1_p51 <- av1_col(av1, "p51_s", NULL, NULL, paste0("Reactive astro/vasc (P47)"),
                  "rel.\nSCT", viridis=TRUE)
a1_it  <- av1_col(av1, "upit_s", NULL, "#cb181d", "Upper-IT weight\n(+P47 hi-zone, blue)",
                  "rel.\nweight") + theme(strip.text.y.left=element_blank())
avcap1 <- ggplot()+theme_void()+
  annotate("text",x=0.5,y=0.5,hjust=0.5,vjust=0.5,size=1.75,fontface="bold",lineheight=1.2,
    label="Region-specific\nexclusion (M1 vs V1)\nblue = P47 high-zone")+xlim(0,1)+ylim(0,1)
avoid1_block <- (avcap1 | a1_p51 | a1_it) +
  plot_layout(widths=c(0.42,1,1), guides="collect") &
  theme(legend.position="bottom",
        legend.title=element_text(size=FS-1), legend.text=element_text(size=FS-1),
        legend.key.height=unit(2.2,"mm"), legend.key.width=unit(5,"mm"))
avoid1_block <- wrap_elements(full = avoid1_block)

# =====================================================================
# NEW AVOIDANCE 2 : OPC  (-)  dense neurofilament/microtubule cytoskeleton (P8) — single M1 chip
#   [OPC weight]  [P8 cytoskeleton SCT]  [OPC weight + P8 high-zone overlay (blue)]
#   (scatter dropped to save space; r reported in caption)
#   P8 = pan-neuronal neurofilament/microtubule cytoskeleton (renamed from "IT neuropil").
# =====================================================================
av2 <- fread(file.path(SPDIR,"avoid2_bins.csv.gz"))
av2_panel <- function(dat, sval, high_col, title, lab, viridis=FALSE, overlay=FALSE){
  dv <- dat[valid==TRUE]; di <- dat[valid==FALSE]; dhi <- dat[p8_hi==TRUE]
  g <- ggplot() +
    geom_point(data=di, aes(x,y), color="#dcdcdc", size=PT_SZ, shape=15) +
    geom_point(data=dv, aes(x,y, color=.data[[sval]]), size=PT_SZ, shape=15)
  if(overlay) g <- g + geom_point(data=dhi, aes(x,y), color="#08519c",
                                  size=PT_SZ*0.75, shape=15, alpha=0.42)
  if(viridis){
    g <- g + scale_color_viridis_c(option="viridis", limits=c(0,1), name=lab,
              breaks=c(0,0.5,1), labels=c("lo","","hi"))
  } else {
    g <- g + scale_color_gradient(low="#fff5f0", high=high_col, limits=c(0,1), name=lab,
              breaks=c(0,0.5,1), labels=c("lo","","hi"))
  }
  g + coord_fixed() + scale_y_reverse() +
    guides(color=guide_colorbar(title.position="left", title.vjust=1)) +
    labs(title=title) + sp_theme +
    theme(plot.title=element_text(size=FS-0.5, hjust=0.5, lineheight=0.9, face="bold"))
}
a2_opc  <- av2_panel(av2, "opc_s", "#cb181d", "M1\nOPC weight", "rel.\nweight")
a2_p8   <- av2_panel(av2, "p8_s",  NULL,      "Neurofilament cytoskel. (P8)\nSCT score", "rel.\nSCT", viridis=TRUE)
a2_ov   <- av2_panel(av2, "opc_s", "#cb181d", "OPC weight\n(+P8 hi-zone, blue)", "rel.\nweight",
                     overlay=TRUE)
avcap2 <- ggplot()+theme_void()+
  annotate("text",x=0.5,y=0.5,hjust=0.5,vjust=0.5,size=1.75,fontface="bold",lineheight=1.2,
    label=sprintf("OPC tiles territory\ncomplementary to dense\nneuropil (r=%s)", av2_r))+
  xlim(0,1)+ylim(0,1)
avoid2_block <- (avcap2 | a2_opc | a2_p8 | a2_ov) +
  plot_layout(widths=c(0.42,1,1,1), guides="collect") &
  theme(legend.position="bottom",
        legend.title=element_text(size=FS-1), legend.text=element_text(size=FS-1),
        legend.key.height=unit(2.2,"mm"), legend.key.width=unit(5,"mm"))
avoid2_block <- wrap_elements(full = avoid2_block)

# =====================================================================
# UNUSED LEGACY CO-LOCALIZATION HEROES (positive niches) retained for traceability.
# The public Fig. 4 composition uses panels f/g/h/i.
#   m : L6 IT (+) P49 Microglial complement/MHC   (M1 A01186A4, +0.41, fss=1.0, LOGO J=0.01)
#   n : VIP   (+) P30 Cholinergic synapse         (AG C00841F3, +0.75, fss=1.0, LOGO survives)
#   each : [L6 IT weight] [program SCT] [L6 IT weight + program hi-zone overlay (blue)]
#   recipe = avoid2 sibling but co-localization (overlay sits ON the cell territory).
# =====================================================================
cln  <- fread(file.path(SPDIR,"coloc_names.tsv"))
clnm <- setNames(cln$value, cln$key)
p54_short <- clnm[["p54_short"]]   # Microglial complement/MHC
p33_short <- clnm[["p33_short"]]   # Cholinergic synapse
coloc1_r  <- clnm[["coloc1_r"]]
coloc2_r  <- clnm[["coloc2_r"]]

coloc_panel <- function(dat, sval, hi_col, high_col, title, lab, viridis=FALSE, overlay=FALSE){
  dv <- dat[valid==TRUE]; di <- dat[valid==FALSE]
  dhi <- if(overlay) dat[get(hi_col)==TRUE] else dat[0]
  g <- ggplot() +
    geom_point(data=di, aes(x,y), color="#dcdcdc", size=PT_SZ, shape=15) +
    geom_point(data=dv, aes(x,y, color=.data[[sval]]), size=PT_SZ, shape=15)
  if(overlay) g <- g + geom_point(data=dhi, aes(x,y), color="#08519c",
                                  size=PT_SZ*0.75, shape=15, alpha=0.42)
  if(viridis){
    g <- g + scale_color_viridis_c(option="viridis", limits=c(0,1), name=lab,
              breaks=c(0,0.5,1), labels=c("lo","","hi"))
  } else {
    g <- g + scale_color_gradient(low="#fff5f0", high=high_col, limits=c(0,1), name=lab,
              breaks=c(0,0.5,1), labels=c("lo","","hi"))
  }
  g + coord_fixed() + scale_y_reverse() +
    guides(color=guide_colorbar(title.position="left", title.vjust=1)) +
    labs(title=title) + sp_theme +
    theme(plot.title=element_text(size=FS-0.5, hjust=0.5, lineheight=0.9, face="bold"))
}

# ---- m : L6 IT (+) P49 (M1; raw component 54) ----
co1 <- fread(file.path(SPDIR,"coloc1_bins.csv.gz"))
m_l6  <- coloc_panel(co1, "l6it_s", NULL,     "#cb181d", "M1\nL6 IT weight", "rel.\nweight")
m_p54 <- coloc_panel(co1, "p54_s",  NULL,     NULL,      sprintf("Microglial compl./MHC (P49)\nSCT score"),
                     "rel.\nSCT", viridis=TRUE)
m_ov  <- coloc_panel(co1, "l6it_s", "p54_hi", "#cb181d", "L6 IT weight\n(+P49 hi-zone, blue)", "rel.\nweight",
                     overlay=TRUE)
mcap <- ggplot()+theme_void()+
  annotate("text",x=0.5,y=0.5,hjust=0.5,vjust=0.5,size=1.75,fontface="bold",lineheight=1.2,
    label=sprintf("L6 IT co-localizes with microglial\ncomplement/MHC (P49)\nlog2 g=+0.41, all 44 chips +"))+
  xlim(0,1)+ylim(0,1)
coloc1_block <- (mcap | m_l6 | m_p54 | m_ov) +
  plot_layout(widths=c(0.42,1,1,1), guides="collect") &
  theme(legend.position="bottom",
        legend.title=element_text(size=FS-1), legend.text=element_text(size=FS-1),
        legend.key.height=unit(2.2,"mm"), legend.key.width=unit(5,"mm"))
coloc1_block <- wrap_elements(full = coloc1_block)

# ---- n : VIP (+) P30 Cholinergic synapse (AG; raw component 33) ----
# REPLACES prior L6 IT x P17 (trivial IT x IT after P17 rename). VIP (CCK+) interneurons
# co-localise with the cholinergic-synapse program (dom PAX6) -> upper-layer disinhibitory
# microcircuit assembling in space. log2 g=+0.75, frac_same_sign=1.0 (44/44 chips), LOGO survives.
co2 <- fread(file.path(SPDIR,"coloc2_bins.csv.gz"))
n_l6  <- coloc_panel(co2, "vip_s", NULL,     "#cb181d", "AG\nVIP weight", "rel.\nweight")
n_p17 <- coloc_panel(co2, "p33_s",  NULL,     NULL,      "Cholinergic synapse (P30)\nSCT score",
                     "rel.\nSCT", viridis=TRUE)
n_ov  <- coloc_panel(co2, "vip_s", "p33_hi", "#cb181d", "VIP weight\n(+P30 hi-zone, blue)", "rel.\nweight",
                     overlay=TRUE)
ncap <- ggplot()+theme_void()+
  annotate("text",x=0.5,y=0.5,hjust=0.5,vjust=0.5,size=1.75,fontface="bold",lineheight=1.2,
    label=sprintf("VIP interneurons co-organize\nwith cholinergic synapse (P30)\nlog2 g=+0.75, all 44 chips +"))+
  xlim(0,1)+ylim(0,1)
coloc2_block <- (ncap | n_l6 | n_p17 | n_ov) +
  plot_layout(widths=c(0.42,1,1,1), guides="collect") &
  theme(legend.position="bottom",
        legend.title=element_text(size=FS-1), legend.text=element_text(size=FS-1),
        legend.key.height=unit(2.2,"mm"), legend.key.width=unit(5,"mm"))
coloc2_block <- wrap_elements(full = coloc2_block)

# =====================================================================
# SVG export tail: write standalone zero-margin SVGs from the panel objects
# defined above. Each panel is rendered at its natural size, ink-cropped, and
# assembled into the vector composite by the published compose step.
#   emit-ids -> final printed tag (compose --tag-sequential, reading order, skip t):
#   a=pa(54-program heatmap)->a  b=pb(forest)->b  c=pc(g(r))->c  d=pd(volcano)->d
#   e=pe(LOGO)->e  f=hero_block->f  g=zoom_block->g  h=avoid1->h
#   i=avoid2->i  j=pf(region-var)->j  k=pg(region-heatmap)->k  l=ph(span)->l  t=title(no tag)
# =====================================================================
suppressPackageStartupMessages(library(svglite))
SVGD <- "CORTEX_PROGRAM_ROOT/scripts/figmarkcorr_A/svg_panels"
dir.create(SVGD, recursive=TRUE, showWarnings=FALSE)
mm2in <- function(x) as.numeric(x)/25.4
svgf  <- function(id) file.path(SVGD, sprintf("figA_%s.svg", id))
published_svgf <- function(id) file.path(SVGD, sprintf("%s.svg", id))

# strip outer whitespace on ggplot panels at export time (ink-crop finishes job)
zero_m  <- theme(plot.margin = margin(0,0,0,0))

emit <- function(id, obj, w_mm, h_mm){
  svglite(svgf(id), width=mm2in(w_mm), height=mm2in(h_mm), bg="white")
  print(obj); invisible(dev.off())
  target <- published_svgf(id)
  if (!file.exists(target) || !nzchar(Sys.readlink(target))) {
    file.copy(svgf(id), target, overwrite=TRUE)
  }
  cat(sprintf("panel %s SVG (%.0fx%.0f mm)\n", id, w_mm, h_mm))
}

# ---- a (master heatmap) + d (volcano): current 54-program overview anchors ----
# Panel a is generated by the overview ComplexHeatmap source so it retains the
# annotation bars, separators, LOGO dots, and legends.
overview_a <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1/markcorr_v2/figures/markcorr_overview/figA_cellprog_heatmap.svg"
if (!file.exists(overview_a)) {
  stop("Missing overview panel a SVG: run results/crossregion_v1/markcorr_v2/scripts/figmarkcorr_overview/fig_overview.R first")
}
file.copy(overview_a, svgf("a"), overwrite=TRUE)
file.copy(overview_a, published_svgf("a"), overwrite=TRUE)
cat("panel a SVG copied from the current 54-program overview\n")
emit("d", pd + zero_m,                70,  43)  # effect x cross-chip reproducibility

# ---- b,c,e, j, k, l : abstract ggplot / heatmap panels ----
emit("b", pb + zero_m,                52, 62)   # headline forest
emit("c", pc + zero_m,                52, 70)   # distance curves (2 facets stacked)
emit("e", pe + zero_m,                58, 62)   # LOGO dumbbell
emit("j", pf + zero_m,                70, 64)   # region-variation exemplars (4 facet)
# pg / ph are heatmap-grob (k) and ggplot dumbbell (l)
emit("k", pg,                         62, 56)   # region-clustering heatmap (grob)
emit("l", ph + zero_m,                96, 40)   # span dumbbell (wide)

# ---- f,g,h,i : current spatial blocks (patchwork comps) -> ONE tight SVG each ----
# tighten: kill block outer margins; sizes chosen ~ to natural map AR so
# ink-crop leaves near-zero inter-map gap (the main dead-space fix).
emit("f", hero_block   & zero_m,      120, 96)  # HERO 3chip x 3field square grid
emit("g", zoom_block   & zero_m,      120, 42)  # OLIGO/myelin co-localization and IT/myelin avoidance
emit("h", avoid1_block & zero_m,      96,  52)  # P47 x upper-IT region contrast
emit("i", avoid2_block & zero_m,      120, 42)  # OPC x P8 IT-neuropil avoidance (cap|3 maps)

# ---- t : title + subtitle band (FROZEN content: honest title, LOGO note,
#          *=FDR legend, GM/WM framing). Full-width thin band, top of figure. ----
ttl <- "Spatially reproducible cell-type-program co-localization across human cortex"
sub <- paste0(
  "a, 22 cell types × 54 programs, log2 median g(r=25µm), 44 chips; labels use the current P1–P54 program namespace. ",
  "b, headline niches (LOGO-surviving). c, distance decay. ",
  "d, effect × cross-chip reproducibility. ",
  "e, leave-one-gene-out gene removal. ",
  "f, SPL/V1/M1 tissue exemplars: OLIGO weight, P41 myelination and L2-L3 IT weight. ",
  "g, zoomed OLIGO-myelination co-localization and L2-L3 IT-myelination avoidance. ",
  "h, region-specific reactive-astro/vascular P47 excludes upper-layer IT in M1 but not V1. ",
  "OPC–neurofilament/microtubule cytoskeleton P8 (i). ",
  "j–l, regional modulation. Headline co-localizations verified by leave-one-gene-out (e), ",
  "excluding shared-gene artifacts. *=brain-weak functional annotation.")
ptitle <- ggplot() + theme_void() +
  annotate("text", x=0, y=1, hjust=0, vjust=1, size=8/.pt, fontface="bold", label=ttl) +
  annotate("text", x=0, y=0.66, hjust=0, vjust=1, size=5.0/.pt, color="grey30",
           lineheight=1.0,
           label=paste(strwrap(sub, width=155), collapse="\n")) +
  xlim(0,1) + ylim(0,1) +
  theme(plot.margin=margin(0,0,0,0))
emit("t", ptitle, 180, 20)

cat("\nSVG PANELS WRITTEN (a,b,c,d,e,f,g,h,i,j,k,l,t) ->", SVGD,
    "\n(current 54-program svgutils relayout; legacy m/n not composed)\n")

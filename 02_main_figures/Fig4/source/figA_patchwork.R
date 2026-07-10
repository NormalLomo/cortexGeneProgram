#!/usr/bin/env Rscript
# =====================================================================
# Figure A: cell-type x program spatial co-localization.
# Standalone patchwork assembly. 8 panels (a-h). cairo_pdf vector.
# Data: cortex_nmf_program markcorr_v2 median+IQR products (ALL REAL).
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ComplexHeatmap)
  library(circlize); library(grid); library(ggrepel); library(dendextend)
  library(data.table); library(scales)
})

ROOT <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
FIN  <- file.path(ROOT, "markcorr_v2/final")
RIG  <- file.path(ROOT, "markcorr_v2/rigor")
OUT  <- "CORTEX_PROGRAM_ROOT/figures/markcorr_A"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

BASE_FONT <- "sans"      # maps to Liberation Sans (Arial-metric) on this box
FS <- 6                  # base font pt (>=5 enforced)
TAG_SZ <- 9
HEADLINE_THR <- 0.32     # |log2 median g| at r=25um

# ---- label shortener (rule 6: drop " x ", Oligodendrocyte->Oligo, keep P{n}) ----
shorten_lab <- function(x){
  x <- gsub(" x ", " · ", x, fixed=TRUE)   # " x " -> middot
  x <- gsub("Oligodendrocyte", "Oligo", x)
  x <- gsub("OLIGO", "Oligo", x)
  x <- gsub("CHANDELIER", "Chand.", x)
  x <- gsub("LINC00507", "", x)
  x <- gsub("\\s+", " ", x); trimws(x)
}

# ---- load ----
pairs <- fread(file.path(FIN, "cellprog_pairs_median_iqr.tsv"))
barea <- fread(file.path(FIN, "cellprog_byarea_median_iqr.tsv"))
logo  <- fread(file.path(RIG, "cellprog_logo.tsv"))
jacc  <- fread(file.path(RIG, "cellprog_jaccard.tsv"))
nm    <- fread(file.path(ROOT, "program_names.tsv"))

# program short label "P{n} name", weak (fdr>0.05 brain-weak) marked with *
nm[, pid := as.integer(program)]
nm[, weak := (confidence == "brain-weak") | (fdr > 0.05)]
nm[, plab := sprintf("P%d %s%s", pid, name_short, ifelse(weak, "*", ""))]
plab_map <- setNames(nm$plab, paste0("program_", nm$pid))
pshort   <- setNames(nm$name_short, paste0("program_", nm$pid))

# program -> NMF class proxy from prev_top_term not available; use brain_term as class color later
# cell-type class grouping (broad) for the master-heatmap class bar
ct_class <- c(
  AST="Astro", OLIGO="Oligo", OPC="Oligo", MICRO="Immune", ENDO="Vascular",
  VLMC="Vascular", `L2-L3 IT LINC00507`="Exc-IT", `L3-L4 IT RORB`="Exc-IT",
  `L4-L5 IT RORB`="Exc-IT", `L6 IT`="Exc-IT", `L6 IT CAR3`="Exc-IT",
  ET="Exc-deep", NP="Exc-deep", `L6 CT`="Exc-deep", L6B="Exc-deep",
  `L6 CAR3`="Exc-deep", PVALB="Inh-MGE", SST="Inh-MGE", CHANDELIER="Inh-MGE",
  LAMP5="Inh-CGE", VIP="Inh-CGE", NDNF="Inh-CGE", PAX6="Inh-CGE")
# (some keys may not match exactly; fallback "Other")

r25 <- pairs[ring_um == 25]
r25[, plab := plab_map[B]]
r25[, ct := A]

# ---- headline-program set: ONLY programs that surface as a named row in b/e/g/h ----
# (replicate those panels' top-N selections up front so panel a labels match exactly)
surv_keys_all <- paste(logo[survives==TRUE]$A, logo[survives==TRUE]$B)
r25[, sk := paste(A,B)]
.hl  <- r25[abs(log2_median_g) >= HEADLINE_THR]
.posB <- .hl[log2_median_g>0][sk %in% surv_keys_all][order(-log2_median_g)][1:min(12,.N)]$B  # panel b/g pos
.negB <- .hl[log2_median_g<0][order(log2_median_g)][1:min(5,.N)]$B                            # panel b neg
.lgB  <- {l<-copy(logo)[!is.na(orig_log2g)&is.finite(logo_log2g)][abs(orig_log2g)>=HEADLINE_THR][order(-abs(orig_log2g))][1:min(14,.N)]; l$B}  # panel e
.spB  <- {s<-barea[ring_um==25 & area %in% c("AG","DLPFC","FPPFC","M1","S1","SMG","SPL","V1","VLPFC")][,.(lo=min(log2_median_g,na.rm=TRUE),hi=max(log2_median_g,na.rm=TRUE),n=.N),by=.(A,B)][n>=7]; s[,span:=hi-lo]; s[order(-span)][1:10]$B}  # panel h
headline_progs <- unique(c(.posB,.negB,.lgB,.spB))
headline_progs <- headline_progs[!is.na(headline_progs)]

# ---------- theme ----------
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
# Panel a : master heatmap 22 ct x 60 prog @ r=25, ComplexHeatmap grabExpr
# =====================================================================
M <- dcast(r25, A ~ B, value.var="log2_median_g")
rn <- M$A; M$A <- NULL
Mm <- as.matrix(M); rownames(Mm) <- rn
# order columns by program number
ord_cols <- paste0("program_", 1:60); ord_cols <- ord_cols[ord_cols %in% colnames(Mm)]
Mm <- Mm[, ord_cols]
# selective column labels: headline programs get "P{n} name", rest bare number
col_pid <- as.integer(sub("program_", "", ord_cols))
# panel a is a LOCATOR: named cols show just "P{n}" bold-ish; full names live in b/e.
# (avoids grabExpr clipping of long rotated text; numbers stay legible)
col_lab <- ifelse(ord_cols %in% headline_progs,
                  sprintf("P%d", col_pid),
                  sprintf("%d", col_pid))
col_fontsz <- ifelse(ord_cols %in% headline_progs, 5.2, 4.2)
colnames(Mm) <- ord_cols   # keep ids as colnames referencely
cls <- ifelse(rownames(Mm) %in% names(ct_class), ct_class[rownames(Mm)], "Other")
cls_col <- c(Astro="#1B9E77", Oligo="#7570B3", Immune="#E7298A",
             Vascular="#666666", `Exc-IT`="#D95F02", `Exc-deep`="#A6761D",
             `Inh-MGE`="#1F78B4", `Inh-CGE`="#66A61E", Other="grey80")
mx <- max(abs(Mm), na.rm=TRUE); mx <- min(mx, 1.6)
ha <- rowAnnotation(Class=cls, col=list(Class=cls_col),
        annotation_name_gp=gpar(fontsize=FS-0.5),
        annotation_legend_param=list(Class=list(title_gp=gpar(fontsize=FS-0.5),
          labels_gp=gpar(fontsize=FS-1), grid_height=unit(3,"mm"), grid_width=unit(3,"mm"))))
ht <- Heatmap(Mm, name="log2g",
        col=colorRamp2(seq(-mx,mx,length.out=101), pal_div),
        cluster_rows=TRUE, cluster_columns=TRUE,
        show_row_dend=TRUE, show_column_dend=TRUE,
        row_dend_width=unit(5,"mm"), column_dend_height=unit(5,"mm"),
        row_names_gp=gpar(fontsize=5.2),
        column_labels=col_lab,
        column_names_gp=gpar(fontsize=col_fontsz),
        column_names_rot=90, column_names_side="top",  # labels on top -> never clipped by panel b
        column_dend_side="bottom",                      # dend below body, labels above
        left_annotation=ha,
        column_names_max_height=unit(10,"mm"),
        row_names_max_width=unit(28,"mm"),
        rect_gp=gpar(col="white", lwd=0.15),
        width=unit(60*1.85,"mm"), height=unit(22*1.35,"mm"),
        column_title="22 cell types x 60 programs  |  log2 median g(r=25um)  (P-prefixed cols = headline programs, named in b/e/g/h)",
        column_title_gp=gpar(fontsize=FS),
        heatmap_legend_param=list(title="log2 g", title_gp=gpar(fontsize=FS-0.5),
          labels_gp=gpar(fontsize=FS-1), legend_height=unit(14,"mm"),
          grid_width=unit(3,"mm")))
pa <- grid.grabExpr(draw(ht, heatmap_legend_side="right",
        annotation_legend_side="right", merge_legend=TRUE,
        padding=unit(c(3,1,1,1),"mm")))   # short P{n}/num labels -> small pad, no clip/float
pa <- wrap_elements(full = pa)

# =====================================================================
# Panel b : headline niche forest (LOGO-surviving) - point=log2g, whisker=IQR
# =====================================================================
r25[, surv := { key <- paste(A,B); s <- setNames(logo$survives, paste(logo$A,logo$B));
                ifelse(key %in% names(s), s[key], NA) }]
hl <- r25[abs(log2_median_g) >= HEADLINE_THR]
pos <- hl[log2_median_g > 0][order(-log2_median_g)]
# keep LOGO-surviving positives
surv_keys <- paste(logo[survives==TRUE]$A, logo[survives==TRUE]$B)
pos[, sk := paste(A,B)]
pos_s <- pos[sk %in% surv_keys][1:min(12,.N)]
neg <- hl[log2_median_g < 0][order(log2_median_g)][1:min(5,.N)]
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
  theme(axis.text.y=element_text(size=5.2), legend.position="bottom",
        legend.margin=margin(0,0,0,0),
        plot.margin=margin(3,3,3,5))

# =====================================================================
# Panel c : g(r) distance-decay curves, IQR ribbon, 2 archetypes
# =====================================================================
# exemplars: oligo broad plateau vs interneuron micro-niche
ex <- list(c("OLIGO","program_45"), c("OLIGO","program_37"),
           c("PVALB","program_7"), c("CHANDELIER","program_50"),
           c("SST","program_20"), c("L2-L3 IT LINC00507","program_45"))
exdt <- rbindlist(lapply(ex, function(e){
  d <- pairs[A==e[1] & B==e[2]]; if(nrow(d)==0) return(NULL)
  d[, lab := shorten_lab(sprintf("%s x %s", A, pshort[B]))]; d }))
# archetype tag
broad_ct <- c("OLIGO"); exdt[, arch := ifelse(A %in% broad_ct, "domain-scale (broad)", "micro-niche (single-bin)")]
exdt[, log2g := log2_median_g]
exdt[, lo := log2(q1_g)]; exdt[, hi := log2(q3_g)]
# per-series linetype + shape so near-zero lines stay distinguishable beyond color
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
        legend.key.size=unit(3.4,"mm"),
        legend.box.margin=margin(2,0,0,0),
        legend.spacing.y=unit(0.5,"mm")) +
  guides(color=guide_legend(ncol=2, byrow=TRUE),
         linetype=guide_legend(ncol=2, byrow=TRUE),
         shape=guide_legend(ncol=2, byrow=TRUE))

# =====================================================================
# Panel d : effect x consistency volcano
# =====================================================================
vd <- copy(r25)
vd[, tier := ifelse(abs(log2_median_g)>=HEADLINE_THR,
              ifelse(log2_median_g>0,"headline +","headline -"),"sub-threshold")]
vd[, lab := shorten_lab(sprintf("%s x %s", ct, pshort[B]))]
set.seed(7)
vd[, yj := frac_same_sign + runif(.N, -0.012, 0.012)]   # light jitter de-stripe
# callouts: only a FEW extreme pairs (top by |effect|)
toplab <- vd[order(-abs(log2_median_g))][1:8]
dxr <- max(abs(vd$log2_median_g), na.rm=TRUE)*1.05         # symmetric x range -> fill cell
pd <- ggplot(vd, aes(log2_median_g, yj, color=tier)) +
  geom_vline(xintercept=c(-HEADLINE_THR,HEADLINE_THR), linetype=2, linewidth=0.25, color="grey60") +
  geom_point(size=0.7, alpha=0.7) +
  geom_text_repel(data=toplab, aes(label=lab), size=1.8, max.overlaps=8,
                  segment.size=0.2, min.segment.length=0, color="grey20",
                  box.padding=0.3, force=2) +
  scale_color_manual(values=c(`headline +`=col_pos, `headline -`=col_neg,
                              `sub-threshold`="grey80"), name=NULL) +
  scale_x_continuous(limits=c(-dxr, dxr)) +
  labs(x=expression(log[2]~median~g~"(r=25"*mu*"m)"), y="fraction same-sign chips",
       title="Effect x consistency") +
  theme(legend.position="bottom", legend.key.size=unit(2.6,"mm"))

# =====================================================================
# Panel e : LOGO retention dumbbell
# =====================================================================
lg <- copy(logo)
lg[, plab := plab_map[B]]
lg[, lab := shorten_lab(sprintf("%s x %s", A, pshort[B]))]
lg <- lg[!is.na(orig_log2g) & is.finite(logo_log2g)]
lge <- lg[order(-abs(orig_log2g))][abs(orig_log2g)>=HEADLINE_THR][1:min(14,.N)]
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
  labs(x=expression(log[2]~median~g), y=NULL,
       title="LOGO gene-removal retention") +
  theme(axis.text.y=element_text(size=5.2), legend.position="bottom",
        legend.key.size=unit(2.6,"mm"),
        plot.margin=margin(3,3,3,5)) +
  expand_limits(x=max(lge$orig_log2g)+0.35)

# =====================================================================
# Panel f : region small-multiples (exemplar niches across 9 stable regions)
# =====================================================================
stable <- c("AG","DLPFC","FPPFC","M1","S1","SMG","SPL","V1","VLPFC")
unstable <- c("ACC","ITG","PoCG","S1E","STG")
exf <- list(c("OLIGO","program_37"), c("PVALB","program_7"),
            c("AST","program_15"), c("SST","program_20"))
fdt <- rbindlist(lapply(exf, function(e){
  d <- barea[A==e[1] & B==e[2] & ring_um==25]; if(nrow(d)==0) return(NULL)
  d[, lab := shorten_lab(sprintf("%s x %s", A, pshort[B]))]; d}))
fdt[, region := factor(area, levels=unique(c(stable,unstable)))]
fdt[, stab := ifelse(area %in% stable, "stable","unstable")]
pf <- ggplot(fdt, aes(region, log2_median_g, fill=stab)) +
  geom_hline(yintercept=0, linewidth=0.3, color="grey70") +
  geom_col(width=0.7) +
  facet_wrap(~lab, ncol=2, scales="free_y") +
  scale_fill_manual(values=c(stable="#3182BD", unstable="grey80"), name=NULL) +
  labs(x=NULL, y=expression(log[2]~median~g), title="Region variation (exemplars)") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=FS-2),
        strip.text=element_text(size=FS-1, face="bold"),
        legend.position="bottom", legend.key.size=unit(2.6,"mm"))

# =====================================================================
# Panel g : region clustering heatmap (9 stable x headline niches)
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
        width=unit(ncol(Gm)*5,"mm"), height=unit(nrow(Gm)*5,"mm"),  # square 5mm cells
        column_title="region clustering (stable)", column_title_gp=gpar(fontsize=FS),
        heatmap_legend_param=list(title="log2 g", title_gp=gpar(fontsize=FS-0.5),
          labels_gp=gpar(fontsize=FS-1), legend_height=unit(10,"mm"), grid_width=unit(3,"mm")))
pg <- grid.grabExpr(draw(htg, heatmap_legend_side="right",
        padding=unit(c(6,4,2,2),"mm")))   # bottom+left pad so rotated x-labels not clipped
pg <- wrap_elements(full=pg)

# =====================================================================
# Panel h : cross-region span dumbbell (min->max across stable regions)
# =====================================================================
sp <- barea[ring_um==25 & area %in% stable]
sp <- sp[, .(lo=min(log2_median_g,na.rm=TRUE), hi=max(log2_median_g,na.rm=TRUE),
             md=median(log2_median_g,na.rm=TRUE), n=.N), by=.(A,B)]
sp <- sp[n>=7]               # present in most stable regions
sp[, span := hi-lo]
sp[, lab := shorten_lab(sprintf("%s x %s", A, pshort[B]))]
sph <- sp[order(-span)][1:10]          # top 10 only -> readable gaps
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
        legend.key.size=unit(3,"mm"),
        plot.margin=margin(3,3,3,5))

# =====================================================================
# assemble : portrait ~180 x 230 mm
# layout: a spans top full width; then b,c ; d,e ; f,g ; h
# =====================================================================
design <- "
AAAA
AAAA
AAAA
BBCC
BBCC
DDEE
DDEE
FFGG
FFGG
FFGG
HHHH
HHHH
"
fig <- pa + pb + pc + pd + pe + pf + pg + ph +
  plot_layout(design=design,
    heights=c(0.85,0.85,0.85, 1,1, 1,1, 1,1,1, 1,1)) +
  plot_annotation(tag_levels="a",
    title="Figure A  |  Cell-type x gene-program spatial co-localization across human cortex",
    subtitle="Mark-correlation g(r) median+IQR over 44 chips; LOGO-rigorous headline niches; pan-region commonality (a-e) and regional modulation (f-h)",
    theme=theme(plot.title=element_text(size=FS+2, face="bold"),
                plot.subtitle=element_text(size=FS, color="grey35")))

pdf_path <- file.path(OUT, "figA_full.pdf")
png_path <- file.path(OUT, "figA_full.png")
cairo_pdf(pdf_path, width=180/25.4, height=230/25.4, family=BASE_FONT)
print(fig); dev.off()

ragg_ok <- requireNamespace("ragg", quietly=TRUE)
if (ragg_ok) {
  ragg::agg_png(png_path, width=180/25.4, height=230/25.4, units="in", res=200)
  print(fig); dev.off()
} else {
  png(png_path, width=180/25.4, height=230/25.4, units="in", res=200, type="cairo")
  print(fig); dev.off()
}
cat("DONE\n", pdf_path, "\n", png_path, "\n")

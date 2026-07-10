suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scatterpie); library(ggridges); library(ComplexHeatmap); library(circlize)
  library(grid)
})
work <- "CORTEX_PROGRAM_ROOT/scripts/fig2/"
outd <- "CORTEX_PROGRAM_ROOT/figures/fig2/"

LAYER_ORDER <- c("ARACHNOID","L1","L2","L3","L4","L5","L6","WM")
LAYER_COLORS <- c(ARACHNOID="#9E9E9E",L1="#3B4CC0",L2="#5A78D6",L3="#7DA0E0",
                  L4="#36A66B",L5="#E8A93B",L6="#D65A3B",WM="#7B3294")
CELLTYPES <- c("AST","CHANDELIER","ENDO","ET","L2-L3 IT LINC00507","L3-L4 IT RORB",
 "L4-L5 IT RORB","L6 CAR3","L6 CT","L6 IT","L6B","LAMP5","MICRO","NDNF","NP",
 "OLIGO","OPC","PAX6","PVALB","SST","VIP","VLMC")
# 22-class palette (qualitative, color-blind aware blend)
CT_COLORS <- c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df8a",
 "#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b","#c49c94","#e377c2",
 "#f7b6d2","#7f7f7f","#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5",
 "#393b79","#637939")
names(CT_COLORS) <- CELLTYPES

theme_pub <- function(base=7) theme_classic(base_size=base) +
  theme(axis.text=element_text(color="black"),
        axis.line=element_line(linewidth=0.3),
        axis.ticks=element_line(linewidth=0.3),
        plot.title=element_text(size=base+1, face="plain", hjust=0),
        legend.key.size=unit(3,"mm"), legend.text=element_text(size=base-1))

ch <- readLines(paste0(work,"_choices.txt"))
REP <- sub("REP=","",grep("^REP=",ch,value=TRUE)[1])
ex_lines <- grep("^EX\t",ch,value=TRUE)
EX <- setNames(sub(".*\t.*\t","",ex_lines), sub("^EX\t","",sub("\t[^\t]*$","",ex_lines)))

## ---- program -> functional GO:BP label lookup (from program_names.tsv) ----
## label = "P{n} {name_short}"  (trailing " P{n}" stripped; "*" suffix if brain-weak)
NM <- fread("CORTEX_PROGRAM_ROOT/results/crossregion_v1/program_names.tsv")
.mk_label <- function(n) {
  r <- NM[program==as.integer(n)]
  if (nrow(r)==0) return(paste0("P",n))
  ns <- trimws(sub("\\s+P[0-9]+$","",as.character(r$name_short[1])))   # strip trailing " P{n}"
  star <- if (as.character(r$confidence[1])=="brain-weak") "*" else ""
  paste0("P",n," ",ns,star)
}
PLABEL <- setNames(vapply(NM$program,.mk_label,character(1)), as.character(NM$program))
plabel <- function(prog_key){   # accepts "program_37" or "37"
  n <- sub("program_","",as.character(prog_key))
  out <- PLABEL[n]; out[is.na(out)] <- paste0("P",n[is.na(out)]); unname(out)
}
## anatomical/cell-type context tag per exemplar (kept as spatial anchor)
EX_CONTEXT <- c(program_5="L2/3",program_35="L4",program_53="L6",
                program_37="OLIGO / WM",program_56="ENDO / vascular",
                program_29="Inhibitory")
ex_title <- function(prog_key){
  ctx <- EX_CONTEXT[as.character(prog_key)]
  lab <- plabel(prog_key)
  ifelse(is.na(ctx), lab, paste0(ctx," · ",lab))
}

## ================= PANEL C: scatterpie of RCTD composition on coarse grid =====
rctd <- fread(paste0(work,"repchip_rctd.tsv"))
meta <- fread(paste0(work,"repchip_meta.tsv"))
m <- merge(meta[,.(bin,x,y)], rctd[rctd_pass_mask==TRUE], by="bin")
# coarse grid: keep pies <=800; bin to ~700um cells (x in um, 50um native)
gx <- 700
m[, gx2 := floor(x/gx)*gx]
m[, gy2 := floor(y/gx)*gx]
agg <- m[, c(lapply(.SD, mean), .(n=.N)), by=.(gx2,gy2), .SDcols=CELLTYPES]
agg <- agg[n>=4]  # require >=4 bins per pie for stability
# normalize comp to sum=1
comp <- as.matrix(agg[,..CELLTYPES]); comp <- comp/rowSums(comp)
agg[,(CELLTYPES):=as.data.table(comp)]
# radius scaled by sqrt(n) modest
agg[, r := gx*0.42]
agg[, yy := max(gy2)-gy2]   # flip y to match panel a
cat("panel c pies:", nrow(agg), "\n")
pc <- ggplot() +
  geom_scatterpie(aes(x=gx2, y=yy, r=r), data=agg, cols=CELLTYPES,
                  color=NA, linewidth=0) +
  scale_fill_manual(values=CT_COLORS, name="cell type") +
  coord_equal() +
  labs(title=sprintf("RCTD composition (chip %s, %d-µm grid)", REP, gx)) +
  theme_void(base_size=7) +
  theme(plot.title=element_text(size=7.2,hjust=0,margin=margin(b=4)),
        plot.margin=margin(t=12,r=4,b=2,l=4),
        legend.key.size=unit(2.6,"mm"), legend.text=element_text(size=5.6),
        legend.title=element_text(size=6.5)) +
  guides(fill=guide_legend(ncol=1))
ggsave(paste0(outd,"fig2_c.pdf"), pc, width=5.2, height=3.8, device=cairo_pdf)
cat("panel c done\n")

## ================= PANEL D: ridgeline of selected programs across depth ========
pcl <- fread(paste0(work,"prog_x_layer_per_chip.tsv"))
selprog <- names(EX)
pd <- pcl[program %in% selprog & n>=50]
pd[, majorDomain := factor(majorDomain, levels=rev(LAYER_ORDER))] # ridge stacks bottom->top
pd[, prog_label := ex_title(program)]
pd[, prog_label := factor(prog_label, levels=ex_title(selprog))]
# ridgeline: distribution of per-chip mean_z across chips, per layer, per program
pdr <- ggplot(pd, aes(x=mean_z, y=majorDomain, fill=after_stat(x))) +
  ggridges::geom_density_ridges_gradient(scale=1.6, rel_min_height=0.01,
      linewidth=0.2, color="grey30") +
  scale_fill_gradient2(low="#2166AC", mid="#F7F7F7", high="#B2182B", midpoint=0,
      name="mean z") +
  facet_wrap(~prog_label, nrow=2) +
  labs(title="Program activity along cortical depth (pia→WM, across 44 chips)",
       x="mean program z-activity", y="cortical depth") +
  theme_pub(7) +
  theme(strip.background=element_blank(), strip.text=element_text(size=6.5),
        panel.spacing=unit(1.5,"mm"))
ggsave(paste0(outd,"fig2_d.pdf"), pdr, width=6.4, height=3.8, device=cairo_pdf)
cat("panel d done\n")

## ================= PANEL E: program x layer ComplexHeatmap (clustered) =========
G <- fread(paste0(work,"prog_x_layer_global.tsv"))
G <- G[match(LAYER_ORDER, majorDomain)]
progcols <- paste0("program_",1:60)
M <- t(as.matrix(G[,..progcols]))         # 60 programs x 8 layers
rownames(M) <- sub("program_","P",progcols)
colnames(M) <- G$majorDomain
col_fun <- colorRamp2(c(-1.2,-0.4,0,0.4,1.2),
   c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#B2182B"))
# column (layer) annotation = depth color
col_ann <- HeatmapAnnotation(layer=colnames(M),
   col=list(layer=LAYER_COLORS), show_legend=FALSE,
   annotation_name_gp=gpar(fontsize=6),
   simple_anno_size=unit(2.5,"mm"))
# label peak cell-type assoc per program from corr file
cc <- fread(paste0(work,"program_celltype_corr.tsv"))
setnames(cc, 1, "program")
peakct <- apply(as.matrix(cc[,..CELLTYPES]),1,function(r) CELLTYPES[which.max(r)])
names(peakct) <- cc$program
# LEGIBILITY: ~60 programs => full GO names overflow. Spotlight programs (those
# featured in panels b/d/f/h) get full functional name; the rest stay "P{n} · {peakct}".
SPOTLIGHT <- c(5,29,35,37,53,56,26)
pnum <- as.integer(sub("program_","",progcols))
rowlab_all <- ifelse(pnum %in% SPOTLIGHT,
                 paste0(plabel(progcols), " · ", peakct[progcols]),
                 paste0(rownames(M), " · ", peakct[progcols]))
names(rowlab_all) <- rownames(M)
# bold + slightly larger for spotlight rows; >=5pt everywhere
face_all <- ifelse(pnum %in% SPOTLIGHT, "bold", "plain"); names(face_all) <- rownames(M)
size_all <- ifelse(pnum %in% SPOTLIGHT, 5.6, 5.2);        names(size_all) <- rownames(M)

# LAYOUT FIX: 60 rows in one column => too tall+narrow. Cluster once, then split the
# clustered order into two equal halves drawn SIDE BY SIDE (30 rows each, half height).
hc <- hclust(dist(M, method="euclidean"), method="ward.D2")
ord <- hc$order                       # row indices in cluster (dendrogram) order
half <- ceiling(length(ord)/2)        # 30
idx1 <- ord[1:half]; idx2 <- ord[(half+1):length(ord)]

mk_half <- function(idx, with_legend, ttl, names_side){
  Mi <- M[idx,,drop=FALSE]
  rn <- rownames(Mi)
  Heatmap(Mi, name="mean z", col=col_fun,
     cluster_columns=FALSE, column_order=LAYER_ORDER,
     cluster_rows=FALSE, row_order=seq_len(nrow(Mi)),   # keep global cluster order
     row_labels=rowlab_all[rn],
     row_names_side=names_side,                         # labels point OUTWARD (no collision)
     row_names_gp=gpar(fontsize=size_all[rn], fontface=face_all[rn]),
     column_names_gp=gpar(fontsize=6), column_names_rot=45,   # rotate to avoid overlap
     top_annotation=HeatmapAnnotation(layer=colnames(Mi), col=list(layer=LAYER_COLORS),
        show_legend=FALSE, annotation_name_gp=gpar(fontsize=5.5),
        simple_anno_size=unit(2.2,"mm")),
     show_heatmap_legend=with_legend,
     heatmap_legend_param=list(title_gp=gpar(fontsize=6.5),labels_gp=gpar(fontsize=5.5),
         legend_height=unit(15,"mm"),grid_width=unit(2.5,"mm")),
     width=unit(22,"mm"), height=unit(76,"mm"),
     column_title=ttl, column_title_gp=gpar(fontsize=7))
}
# left half: labels on LEFT; right half: labels on RIGHT -> both face outward, never clip mid-gap
ht1 <- mk_half(idx1, FALSE, "Program × cortical layer", "left")
ht2 <- mk_half(idx2, TRUE,  "(continued)",              "right")
# concatenate horizontally: 2 half-height columns => wide ~1.3:1 block
ht_list <- ht1 + ht2
pdf(paste0(outd,"fig2_e.pdf"), width=6.6, height=4.2)
draw(ht_list, heatmap_legend_side="right", ht_gap=unit(6,"mm"),
     padding=unit(c(2,2,2,2),"mm"))
dev.off()
cat("panel e done (2-col split)\n")

## ================= PANEL G: cross-chip reproducibility =========================
gsum <- fread(paste0(work,"panelg_summary.tsv"))
grep_ <- fread(paste0(work,"panelg_reproducibility.tsv"))
ord <- gsum[order(median)]$program
grep_[, program := factor(program, levels=ord)]
grep_[, pid := as.integer(sub("program_","",as.character(program)))]
# boxplot per program, colored by median repro
medmap <- setNames(gsum$median, gsum$program)
grep_[, medv := medmap[as.character(program)]]
pg <- ggplot(grep_, aes(x=program, y=corr_to_mean)) +
  geom_hline(yintercept=c(0,0.5,1), color="grey85", linewidth=0.25) +
  geom_boxplot(aes(fill=medv), outlier.size=0.2, outlier.colour="grey60",
               linewidth=0.15, width=0.72) +
  scale_fill_gradientn(colours=c("#D73027","#FEE08B","#1A9850"),
      limits=c(min(grep_$medv),1), name="median r") +
  scale_x_discrete(labels=function(x) sub("program_","P",x)) +
  labs(title="Cross-chip reproducibility of program layer profiles (n=44 chips)",
       x="program (sorted by median)", y="correlation to mean layer profile") +
  coord_cartesian(ylim=c(min(0,quantile(grep_$corr_to_mean,0.01)),1)) +
  theme_pub(7) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=3.6))
ggsave(paste0(outd,"fig2_g.pdf"), pg, width=6.6, height=2.8, device=cairo_pdf)
cat("panel g done\n")

## ================= PANEL H: WM-distance radial profile =========================
rings <- fread(paste0(work,"panelh_rings.tsv"))
rings[, ring := as.numeric(as.character(ring))]
wmp <- c("program_37","program_26","program_5","program_56")
# functional label + retained anatomical role hint (spatial anchor for this panel)
.role <- c(program_37="OLIGO/WM", program_26="OLIGO/WM",
           program_5="L2/3 ctrl", program_56="vascular")
plab <- setNames(paste0(plabel(wmp), " (", .role[wmp], ")"), wmp)
long <- melt(rings, id.vars=c("ring","chip"), measure.vars=wmp,
             variable.name="program", value.name="mean_z")
long[, prog_label := plab[as.character(program)]]
long[, prog_label := factor(prog_label, levels=plab[wmp])]
prog_pal <- setNames(c("#7B3294","#C2A5CF","#5A78D6","#D6604D"), plab[wmp])
ph <- ggplot(long, aes(x=ring, y=mean_z, color=prog_label, group=interaction(prog_label,chip))) +
  geom_hline(yintercept=0, color="grey80", linewidth=0.3) +
  geom_line(aes(linetype=chip), linewidth=0.5) +
  geom_point(size=0.6) +
  scale_color_manual(values=prog_pal, name="program") +
  scale_linetype_manual(values=c("solid","31"), name="chip") +
  labs(title="Program activity vs distance from white-matter seeds (OLIGO≥0.5)",
       x="distance from WM (µm)", y="mean program z-activity") +
  theme_pub(7) +
  theme(legend.position="right", legend.box="vertical")
ggsave(paste0(outd,"fig2_h.pdf"), ph, width=5.0, height=2.8, device=cairo_pdf)
cat("panel h done\n")
cat("R PANELS DONE\n")

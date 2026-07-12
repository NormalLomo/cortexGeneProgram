#!/usr/bin/env Rscript
# Fig. 7 panel a: retained-program cross-species laminar heatmap.
# Labels come from program_names_curated.tsv; rows are clustered on
# the plotted matrix and the center-of-mass statistic remains an annotation.
suppressMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# Compact heatmap and legend spacing.
#   - TITLE_PADDING: gap above/below column titles (default 5.5pt) -> shrink so main title hugs body.
#   - legend_gap / legend_grid_*: tighten the bottom legend band so it hugs the body and itself.
#   - HEATMAP_LEGEND_PADDING / ANNOTATION_LEGEND_PADDING: gap between body and the legend band.
ht_opt(
  TITLE_PADDING = unit(0.8, "mm"),
  HEATMAP_LEGEND_PADDING = unit(0.5, "mm"),
  ANNOTATION_LEGEND_PADDING = unit(0.5, "mm"),
  legend_gap = unit(2.5, "mm")
)

args <- commandArgs(trailingOnly=TRUE)
root_idx <- match("--canonical-root", args)
canonical_root <- if (!is.na(root_idx) && length(args) > root_idx) args[[root_idx + 1]] else Sys.getenv("CORTEX_PROGRAM_CANONICAL_ROOT")
if (!nzchar(canonical_root)) stop("Set --canonical-root or CORTEX_PROGRAM_CANONICAL_ROOT.")

AGG <- file.path(canonical_root, "results/xspecies_humanmap_v1/spatial_xspecies/_aggregate")
OUT_SVG <- file.path(dirname(AGG), "figures/Fig_spatial_univ/panels/panel_a_complex_heatmap.svg")
# Cross-region curated display names keyed by program.
CURATED_TSV <- file.path(dirname(AGG), "figures/Fig_spatial_univ/program_names_curated.tsv")
dir.create(dirname(OUT_SVG), recursive = TRUE, showWarnings = FALSE)

LAYERS <- c("L1","L2","L4","L5","L6")   # dropL3 20260617: 5 common layers (L3 absent in mouse mapping), matches legend/Methods/LOCKED
SPECIES <- c("human","macaque","mouse")
SP_FILE <- c(human="human", macaque="monkey", mouse="mouse")

read_layer <- function(sp, suffix){
  f <- file.path(AGG, paste0(SP_FILE[sp], "_program_x_layer", suffix, ".tsv"))
  m <- as.matrix(read.table(f, sep="\t", header=TRUE, row.names=1,
                            check.names=FALSE, colClasses="character", na.strings=""))
  m <- m[, LAYERS, drop=FALSE]
  apply(m, 2, as.numeric) -> mm
  rownames(mm) <- rownames(m)
  mm
}
# heatmap BODY (mat) = baseline-removed profile (visual decision, unchanged).
read_blr <- function(sp) read_layer(sp, "_baselinerm")
H  <- read_blr("human")
MK <- read_blr("macaque")
MO <- read_blr("mouse")
# RIGHT consistency strip = LOCKED per-program within-species z layer-profile
# The consistency strip uses per-program within-species z profiles (median
# Spearman 0.70 Hu-Mq / 0.90 Hu-Mo). Separate from the body's baseline-removed mat.
read_z <- function(sp) read_layer(sp, "_perprogz")
Hz  <- read_z("human")
MKz <- read_z("macaque")
MOz <- read_z("mouse")

lk <- read.table(file.path(AGG,"laminar_consistency_LOCKED.tsv"), sep="\t",
                 header=TRUE, check.names=FALSE, quote="", stringsAsFactors=FALSE)
rownames(lk) <- lk$program

progs <- rownames(H)
# Drop the six excluded source components so the
# laminar-universality heatmap is computed over the 54 biologically-interpreted programs. Applied
# immediately after rownames(H) so EVERYTHING downstream (matrix, functional-class anno, A/B/C
# laminar-class anno, COM, consistency strips, ward.D2 row clustering + dendrogram, curated row
# names) is rebuilt on 54 automatically. The A/B/C counts therefore become 38/7/9 by data, with
# no hard-coded number to edit. Nothing else in the science/encoding changes.
DROP6 <- paste0("program_", c(9, 18, 19, 35, 52, 57))
progs <- setdiff(progs, DROP6)   # 60 -> 54; all downstream indexing is by-name on progs/progs_o
stopifnot(all(progs %in% rownames(lk)), length(progs) == 54)

func_map <- c(
  program_1="Gene regulation", program_29="Gene regulation", program_51="Gene regulation",
  program_53="Gene regulation", program_57="Gene regulation",
  program_7="Neuronal signaling", program_9="Neuronal signaling", program_17="Neuronal signaling",
  program_19="Neuronal signaling", program_25="Neuronal signaling", program_33="Neuronal signaling",
  program_45="Glia / myelin",
  program_30="Vascular", program_46="Vascular", program_48="Vascular", program_56="Vascular",
  program_27="Immune", program_35="Immune", program_37="Immune", program_40="Immune",
  program_47="Immune", program_49="Immune", program_54="Immune", program_58="Immune", program_60="Immune",
  program_10="ECM / adhesion", program_13="ECM / adhesion", program_22="ECM / adhesion",
  program_31="ECM / adhesion", program_39="ECM / adhesion", program_50="ECM / adhesion",
  program_4="Transport / metabolism", program_5="Transport / metabolism", program_12="Transport / metabolism",
  program_24="Transport / metabolism", program_28="Transport / metabolism", program_44="Transport / metabolism",
  program_52="Transport / metabolism", program_59="Transport / metabolism",
  program_2="Development", program_3="Development", program_8="Development", program_15="Development",
  program_20="Development", program_26="Development", program_34="Development", program_36="Development",
  program_38="Development", program_43="Development", program_55="Development",
  program_18="Sensory perception", program_21="Sensory perception", program_41="Sensory perception",
  program_6="Secretion / signaling", program_14="Secretion / signaling", program_16="Secretion / signaling",
  program_11="Contractile", program_32="Contractile", program_42="Contractile"
)
func_class <- func_map[progs]
func_class[is.na(func_class)] <- "Other"

pair_cor <- function(a, b){
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  if (sd(a[ok])==0 || sd(b[ok])==0) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method="spearman"))
}
# Consistency-strip Spearman values use per-program within-species z profiles.
# (Hz/MKz/MOz), NOT the baseline-removed body matrices (H/MK/MO).
consist_hmq <- sapply(progs, function(p) pair_cor(Hz[p,], MKz[p,]))
consist_hmo <- sapply(progs, function(p) pair_cor(Hz[p,], MOz[p,]))
names(consist_hmq) <- progs; names(consist_hmo) <- progs

com <- lk[progs, "human_pref"]
names(com) <- progs
# Rows are ordered by the dendrogram; COM remains a left annotation.
progs_o <- progs

col_ids <- character(0); col_layer <- character(0); col_sp <- character(0)
mat_list <- list()
for (L in LAYERS){
  for (sp in SPECIES){
    M <- switch(sp, human=H, macaque=MK, mouse=MO)
    col_ids  <- c(col_ids,  paste0(L,".",sp))
    col_layer<- c(col_layer, L)
    col_sp   <- c(col_sp,   sp)
    mat_list[[length(mat_list)+1]] <- M[progs_o, L]
  }
}
mat <- do.call(cbind, mat_list)
rownames(mat) <- progs_o
colnames(mat) <- col_ids

LIM <- 0.35
mat_c <- mat
mat_c[is.finite(mat_c) & mat_c >  LIM] <-  LIM
mat_c[is.finite(mat_c) & mat_c < -LIM] <- -LIM

# Hierarchical clustering of retained programs on the plotted matrix.
#   distance = euclidean, linkage = ward.D2 (reveals block structure most clearly).
#   NA cells (layer not mapped / undefined) -> 0 (neutral baseline-removed score) for the
#   DISTANCE computation ONLY; the displayed matrix keeps NA (grey), values unchanged.
mat_clust <- mat_c
mat_clust[!is.finite(mat_clust)] <- 0
row_hc <- hclust(dist(mat_clust, method = "euclidean"), method = "ward.D2")
col_fun <- colorRamp2(
  c(-LIM, -0.18, -0.07, -0.025, -0.010, 0, 0.010, 0.025, 0.07, 0.18, LIM),
  c("#155E63","#2E9AA0","#74C2C2","#A9D6D2","#DCEDEA","#FBFBF7",
    "#F7E2C4","#F2C386","#EB9B4E","#D9772A","#9C4509")
)

sp_cols  <- c(human="#3B5BA5", macaque="#C0504D", mouse="#4E9E5A")
func_cats <- c("Neuronal signaling","Glia / myelin","Vascular","Secretion / signaling",
               "Transport / metabolism","Gene regulation","ECM / adhesion","Immune",
               "Development","Sensory perception","Contractile","Other")
func_cols <- setNames(
  c("#C44E52","#8172B3","#4C72B0","#55A868","#CCB974","#937860",
    "#DA8BC3","#DD8452","#64B5CD","#B07AA1","#9C755F","#BBBBBB"),
  func_cats)
abc_cols <- c(A="#2E7D32", B="#F9A825", C="#C62828")

func_class_o <- func_class[progs_o]
abc_o <- lk[progs_o, "class_ABC"]
com_o <- com[progs_o]
consist_hmq_o <- consist_hmq[progs_o]
consist_hmo_o <- consist_hmo[progs_o]

# top annotation: species color strip (per column) ; legend HORIZONTAL
top_anno <- HeatmapAnnotation(
  Species = col_sp,
  col = list(Species = sp_cols),
  annotation_label = "Species",
  annotation_name_gp = gpar(fontsize = 8),
  annotation_name_side = "left",
  simple_anno_size = unit(3.6, "mm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    title="Species", labels=c("Human","Macaque","Mouse"),
    at=c("human","macaque","mouse"),
    title_gp=gpar(fontsize=7,fontface="bold"),
    labels_gp=gpar(fontsize=6.5),
    direction="horizontal", nrow=1
  )
)

# left row annotation: functional class + A/B/C + COM ; legends HORIZONTAL
com_col_fun <- colorRamp2(c(min(com_o), median(com_o), max(com_o)),
                          c("#F2E9C8","#C9C9E8","#3A0CA3"))
left_anno <- rowAnnotation(
  `Functional class` = func_class_o,
  `Laminar class` = abc_o,
  `COM` = com_o,
  col = list(`Functional class`=func_cols, `Laminar class`=abc_cols, `COM`=com_col_fun),
  na_col = "grey85",
  annotation_name_gp = gpar(fontsize = 8),
  annotation_name_side = "top",
  annotation_name_rot = 45,
  annotation_name_offset = unit(1.0, "mm"),
  simple_anno_size = unit(3.6, "mm"),
  annotation_legend_param = list(
    `Functional class`=list(title="Functional class",
                            title_gp=gpar(fontsize=7,fontface="bold"),
                            labels_gp=gpar(fontsize=6.5),
                            at=func_cats, direction="horizontal", nrow=2),
    `Laminar class`=list(title="Laminar class",
                       title_gp=gpar(fontsize=7,fontface="bold"),
                       labels_gp=gpar(fontsize=6.5),
                       at=c("A","B","C"),
                       labels=c("A universal","B partial","C divergent"),
                       direction="horizontal", nrow=1),
    `COM`=list(title="Human COM (upper -> deep)",
               title_gp=gpar(fontsize=7,fontface="bold"),
               labels_gp=gpar(fontsize=6.5),
               legend_width=unit(34,"mm"),
               direction="horizontal")
  )
)

# RIGHT row annotation = MAIN narrative: cross-species laminar consistency, WIDENED bars.
cons_col_fun <- colorRamp2(c(-1, -0.5, 0, 0.5, 1),
                           c("#762A83","#C2A5CF","#F2F2F2","#7FBF7B","#1B7837"))
right_anno <- rowAnnotation(
  `Hu-Mq` = consist_hmq_o,
  `Hu-Mo` = consist_hmo_o,
  col = list(`Hu-Mq`=cons_col_fun, `Hu-Mo`=cons_col_fun),
  na_col = "grey85",
  annotation_name_gp = gpar(fontsize = 9, fontface="bold"),
  annotation_name_side = "top",
  annotation_name_rot = 0,
  annotation_name_offset = unit(1.2, "mm"),
  simple_anno_size = unit(6.5, "mm"),     # WIDENED -> main narrative
  gap = unit(1.2, "mm"),
  annotation_legend_param = list(
    `Hu-Mq`=list(
      title="Cross-species laminar consistency (Spearman r): Hu-Mq / Hu-Mo",
      title_gp=gpar(fontsize=8,fontface="bold"),
      labels_gp=gpar(fontsize=6.5),
      at=c(-1,-0.5,0,0.5,1),
      direction="horizontal",
      legend_width=unit(60,"mm"),
      grid_height=unit(5,"mm"),
      border="grey40"),
    `Hu-Mo`=list(show=FALSE)
  )
)

# Row labels use the cross-region curated display-name table.
#   Keyed by `program` (program_1..program_60), reindexed to the EXISTING clustered/display row
#   vector progs_o so the row->label correspondence is preserved 1:1; matrix/order untouched.
cur <- read.table(CURATED_TSV, sep="\t", header=TRUE, check.names=FALSE,
                  quote="", stringsAsFactors=FALSE)
rownames(cur) <- cur$program
stopifnot(all(progs_o %in% rownames(cur)))   # every plotted program has a curated name
row_lab <- cur[progs_o, "curated_name"]

ht <- Heatmap(
  mat_c,
  name = "Baseline-removed laminar score",
  col = col_fun,
  na_col = "grey85",
  cluster_rows = row_hc, cluster_columns = FALSE,
  show_row_dend = TRUE,
  row_dend_side = "left",
  row_dend_width = unit(11, "mm"),                  # moderate dend width, doesn't dominate
  row_dend_gp = gpar(lwd = 0.5),
  column_order = seq_len(ncol(mat_c)),
  column_split = factor(col_layer, levels = LAYERS),
  column_gap = unit(1.8, "mm"),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  show_column_names = FALSE,
  row_names_side = "right",
  row_labels = row_lab,
  row_names_gp = gpar(fontsize = 6),
  top_annotation = top_anno,
  left_annotation = left_anno,
  right_annotation = right_anno,
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.15),
  width  = unit(136, "mm"),
  height = unit(236, "mm"),
  heatmap_legend_param = list(
    title = "Baseline-removed laminar score (shared scale, all 3 species)",
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5),
    at = c(-0.35,-0.18,-0.07,0,0.07,0.18,0.35),
    direction = "horizontal",
    legend_width = unit(66, "mm"),
    border = "grey40"
  )
)

# Device size follows the dendrogram and legend footprint.
svg(OUT_SVG, width = 6.20 + 136/25.4, height = 12.30 + (236-210)/25.4 - 1.28)
draw(ht,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legend = TRUE,
     column_title = "Laminar programs are universal across human, macaque and mouse cortex",
     column_title_gp = gpar(fontsize = 13, fontface = "bold"),
     padding = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))
# subtitle (two lines, left-anchored): universality read from strips
# Keep the subtitle close to the main title.
grid.text(
  "Universality is read from the cross-species consistency strips (per-program profile r); the heatmap shows",
  x = unit(4, "mm"), y = unit(1, "npc") - unit(4.0, "mm"), just = c("left","top"),
  gp = gpar(fontsize = 8, fontface = "italic", col = "#2E7D32"))
grid.text(
  "baseline-removed score on a shared scale - human amplitude is genuinely smaller.   Upper-vs-deep preference: human-mouse 94% / human-macaque 83% sign-concordant.",
  x = unit(4, "mm"), y = unit(1, "npc") - unit(7.6, "mm"), just = c("left","top"),
  gp = gpar(fontsize = 8, fontface = "italic", col = "#444444"))
# footnote
grid.text(
  paste0("Each column block = one cortical layer (L1, L2, L4, L5, L6; L3 absent in mouse); within a block, columns are Human / Macaque / Mouse. ",
         "Cell colour = baseline-removed laminar score on one shared diverging scale (no per-column normalization; ",
         "values clipped at +/-0.35). Rows are hierarchically clustered (euclidean distance, Ward.D2 linkage) on the ",
         "cross-species laminar profile (dendrogram at left); the COM strip shows each program's human centre-of-mass (upper to deep). ",
         "Right strips (main evidence): cross-species laminar consistency = Spearman r of the 5-layer ",
         "baseline-removed profile (Hu-Mq, Hu-Mo). Macaque uses cell-bin mapping (caveat). Grey = layer not mapped / r undefined."),
  x = unit(4, "mm"), y = unit(2.5, "mm"), just = c("left","bottom"),
  gp = gpar(fontsize = 5.6, col = "grey25"))
dev.off()
cat("WROTE", OUT_SVG, "\n")
cat("rows", nrow(mat_c), "cols", ncol(mat_c), "\n")
cat("classA/B/C:", paste(table(abc_o), collapse="/"), "\n")
cat("Hu-Mq consist: median=", round(median(consist_hmq, na.rm=TRUE),3),
    " nNA=", sum(is.na(consist_hmq)), "\n", sep="")
cat("Hu-Mo consist: median=", round(median(consist_hmo, na.rm=TRUE),3),
    " nNA=", sum(is.na(consist_hmo)), "\n", sep="")
# Cluster diagnostic: report COM ranges for several dendrogram cuts.
cat("--- ROW CLUSTERING (euclidean + ward.D2) ---\n")
ord_idx <- row_hc$order
for (K in c(3,4,5)) {
  cl <- cutree(row_hc, k = K)
  cl_disp <- cl[ord_idx]                 # cluster ids in displayed (top->bottom) order
  rle_blocks <- rle(cl_disp)
  cat("k=", K, " block sizes (top->bottom): ", paste(rle_blocks$lengths, collapse="/"), "\n", sep="")
  for (b in unique(cl_disp)) {
    members <- progs_o[ord_idx][cl_disp == b]
    com_b <- com_o[match(members, progs_o)]
    cat("   cluster", b, " n=", length(members),
        " COM[min/med/max]=", round(min(com_b,na.rm=TRUE),2), "/",
        round(median(com_b,na.rm=TRUE),2), "/", round(max(com_b,na.rm=TRUE),2), "\n", sep="")
  }
}

#!/usr/bin/env Rscript
# Panel A v9 = v8 with ROW LABELS ONLY swapped to crossregion-curated brain names (羅老師 v17 brief).
#   THE ONLY CHANGE vs v8 is the DISPLAYED ROW-LABEL STRINGS: they are now read from
#   program_names_curated_forfig_v18.tsv (curated_name; brain-relevant GO:BP terms, P{n} prefix,
#   trailing * on weak-evidence programs) INSTEAD OF the raw _aggregate/program_v9_names.tsv
#   v9_name column (which carried non-brain generic GO terms: Substantia Nigra / Eye / Kidney /
#   Cardiac Muscle / Skeletal / Osteoclast / Endoderm / Taste ...). NOTHING ELSE CHANGES:
#   the 60x15 numeric matrix, cmap, shared scale, NA-grey, black cell border, 5-layer column
#   split, within-block 3 species, species colors, functional-class anno, A-B-C anno, COM anno,
#   Hu-Mq/Hu-Mo consistency anno, ward.D2 row clustering + dendrogram, row ORDER, and all
#   spacing/canvas are BYTE-IDENTICAL to v8. The label is a pure display string; swapping the
#   text does not touch any data/order/colour/clustering. Output -> panelA_complexhm_v10.svg.
# --- v8 header below ---
# Panel A v8 = v7 with reference+OUTER WHITESPACE KILLED AT SOURCE (羅老師 v15 brief).
#   ONLY whitespace/spacing changes vs v7 (SCIENCE / DATA / CMAP / SHARED SCALE / black cell
#   border / 5-layer split / within-block 3 species / species colors / functional-class anno /
#   A-B-C anno / COM anno / Hu-Mq+Hu-Mo consistency anno / 60 program 6pt names / ward.D2 row
#   clustering + dendrogram are BYTE-IDENTICAL to v7). Diagnosed dead bands on the v7 render:
#     (w1) device width carried +0.7in headroom for the dendrogram -> removed (dend is only 11mm).
#     (w2) device height carried ~16mm bottom headroom above the footnote -> shrunk to hug content.
#     (w3) top band: main title -> subtitle gap of ~9mm -> subtitles pulled up tight under title,
#          ht_opt TITLE_PADDING shrunk, draw padding top trimmed.
#     (w4) heatmap body -> bottom legend band ~4mm gap + legend reference gaps -> tightened via
#          ht_opt(legend_gap) and the draw legend padding.
#     (w5) footnote pulled up tight under the legend band (was floating with ~16mm white above).
#   Net: panel a SVG viewBox hugs the ink on all four sides with no reference dead bands.
# --- v7 header below ---
# Panel A v7 = v6 with ROW ORDER changed from COM-sort to HIERARCHICAL CLUSTERING (羅老師 v14 brief).
#   ONLY the row ordering changes vs v6: rows are now ordered by hclust(dist(plot matrix),
#   euclidean + ward.D2) with the row dendrogram SHOWN on the left, so programs with similar
#   cross-species laminar profiles cluster into blocks (universal-upper / universal-deep etc.).
#   COM is NOT dropped: it stays as a left annotation column (already in v6's left_anno), now
#   following the new clustered row order. EVERYTHING ELSE (data / cmap / shared scale / black
#   cell border / 5-layer column split / within-block 3 species / species colors / functional-class
#   anno / A-B-C anno / Hu-Mq+Hu-Mo consistency anno / 60 program 6pt names) byte-identical to v6.
# --- v6 header below ---
# Panel A v6 = v5 + BLACK cell border (羅老師 v12 brief). ONLY rect_gp changed vs v5
#   (col white->black, lwd 0.4->0.15). SCIENCE/DATA/CMAP/ORDER/ANNO/NAMES byte-identical to v5.
# --- v5 header below ---
# Panel A v5 = v4 with TIGHTENED whitespace (羅老師 v11 brief).
# ONLY layout/whitespace changed vs v4 (SCIENCE / DATA / NARRATIVE / ENCODING UNCHANGED):
#   (w1) draw padding c(6,4,12,6)mm -> c(1,1,2,1)mm  (kills the top 12mm + outer dead margins)
#   (w2) heatmap body width 120->136 mm, height 210->236 mm (eats the freed whitespace; cells bigger; rows still 6pt)
#   (w3) svg() canvas shrunk to hug the bigger body + bottom legend band (was 14.6x12.8 inch)
# Everything else is BYTE-IDENTICAL to v4.
# --- original v4 header below ---
# Panel A v4 (CNS-grade ComplexHeatmap) = v3 ENLARGED for v9 headline.
# CHANGES vs v3 (SCIENCE / DATA / NARRATIVE UNCHANGED):
#   (i)   legends moved to a horizontal BAND BELOW the heatmap (heatmap_legend_side="bottom",
#         annotation_legend_side="bottom", every legend direction="horizontal") -> frees the
#         RIGHT width that v3 spent on a tall legend column, giving it to the heatmap body + row names.
#   (ii)  heatmap body enlarged: width 86 -> 120 mm, height 165 -> 210 mm so 60 rows breathe.
#   (iii) row-name fontsize 5 -> 6 (clearer, still >=5pt). Functional-class anno name fontsize bumped slightly.
#   (iv)  canvas (svg width/height) enlarged to host the bigger body + bottom legend band.
# KEPT: baseline-removed shared scale, consistency strips Hu-Mq/Hu-Mo on the RIGHT (main narrative),
#       5-layer split, within-block 3 species, species colors, functional class, A/B/C, COM order,
#       NA grey, subtitles, footnote. Human K60 only. All-English, 0 CJK.
suppressMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# v8: global spacing options to kill the diagnosed dead bands (spacing-only; no science change).
#   - TITLE_PADDING: gap above/below column titles (default 5.5pt) -> shrink so main title hugs body.
#   - legend_gap / legend_grid_*: tighten the bottom legend band so it hugs the body and itself.
#   - HEATMAP_LEGEND_PADDING / ANNOTATION_LEGEND_PADDING: gap between body and the legend band.
ht_opt(
  TITLE_PADDING = unit(0.8, "mm"),
  HEATMAP_LEGEND_PADDING = unit(0.5, "mm"),
  ANNOTATION_LEGEND_PADDING = unit(0.5, "mm"),
  legend_gap = unit(2.5, "mm")
)

AGG <- "__PRIVATE_CANONICAL_ROOT__/results/xspecies_humanmap_v1/spatial_xspecies/_aggregate"
OUT_SVG <- file.path(dirname(AGG), "figures/Fig_spatial_univ/panels/panelA_complexhm_v10.svg")
# v9: crossregion-curated brain-relevant display names (the ONLY change vs v8). Keyed by program.
CURATED_TSV <- file.path(dirname(AGG), "figures/Fig_spatial_univ/program_names_curated_forfig_v18.tsv")
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
# (stripfix 20260617): matches body L119 / legend L461 / Methods L259 (median
# Spearman 0.70 Hu-Mq / 0.90 Hu-Mo). Separate from the body's baseline-removed mat.
read_z <- function(sp) read_layer(sp, "_perprogz")
Hz  <- read_z("human")
MKz <- read_z("macaque")
MOz <- read_z("mouse")

lk <- read.table(file.path(AGG,"laminar_consistency_LOCKED.tsv"), sep="\t",
                 header=TRUE, check.names=FALSE, quote="", stringsAsFactors=FALSE)
rownames(lk) <- lk$program
v9 <- read.table(file.path(AGG,"program_v9_names.tsv"), sep="\t",
                 header=TRUE, check.names=FALSE, quote="", stringsAsFactors=FALSE)
rownames(v9) <- v9$program

progs <- rownames(H)
# fig54甲 (2026-06-20): drop the six cohort-technical programs (P9, P18, P19, P35, P52, P57) so the
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
# stripfix 20260617: strip Spearman from LOCKED per-program within-species z
# (Hz/MKz/MOz), NOT the baseline-removed body matrices (H/MK/MO).
consist_hmq <- sapply(progs, function(p) pair_cor(Hz[p,], MKz[p,]))
consist_hmo <- sapply(progs, function(p) pair_cor(Hz[p,], MOz[p,]))
names(consist_hmq) <- progs; names(consist_hmo) <- progs

com <- lk[progs, "human_pref"]
names(com) <- progs
# v7: rows no longer manually ordered by COM. Keep original program order in the matrix;
# ComplexHeatmap will reorder rows via the row dendrogram (hclust below). COM is kept as a
# left-annotation column and follows the clustered order automatically.
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

# v7: hierarchical clustering of rows (programs) on the 60x15 PLOTTED matrix.
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
               legend_width=unit(34,"mm"),   # v8: widened to push first legend row toward the right edge
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
      legend_width=unit(60,"mm"),   # v8: widened 46->60mm to push the bottom band toward the right edge (kills right dead band)
      grid_height=unit(5,"mm"),
      border="grey40"),
    `Hu-Mo`=list(show=FALSE)
  )
)

# row labels = crossregion-curated brain names (v17): display-string swap ONLY (no data/order change).
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
  cluster_rows = row_hc, cluster_columns = FALSE,   # v7: rows clustered (euclidean + ward.D2); columns stay fixed
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
  rect_gp = gpar(col = "black", lwd = 0.15),   # v6: black cell border (was white 0.4), thin so 60x15 cells stay crisp
  width  = unit(136, "mm"),    # v5: 120 -> 136 (eat freed whitespace; cells bigger)
  height = unit(236, "mm"),    # v5: 210 -> 236 (eat freed whitespace; rows breathe, still 6pt)
  heatmap_legend_param = list(
    title = "Baseline-removed laminar score (shared scale, all 3 species)",
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5),
    at = c(-0.35,-0.18,-0.07,0,0.07,0.18,0.35),
    direction = "horizontal",
    legend_width = unit(66, "mm"),   # v8: widened 52->66mm to fill the bottom band toward the right edge
    border = "grey40"
  )
)

# v8 canvas: drop the +0.7in width headroom (dend is only 11mm) and ~0.5in bottom-height headroom
#   (diagnosed ~16mm dead white above the footnote). Device sized to hug content; ink-crop still
#   trims any residual outer margin to a tight bbox.
svg(OUT_SVG, width = 6.20 + 136/25.4, height = 12.30 + (236-210)/25.4 - 1.28)
draw(ht,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legend = TRUE,
     column_title = "Laminar programs are universal across human, macaque and mouse cortex",
     column_title_gp = gpar(fontsize = 13, fontface = "bold"),
     padding = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))   # v8: c(1,1,2,1) -> c(0.5x4); kill outer dead margin
# subtitle (two lines, left-anchored): universality read from strips
# v8: pulled up tight under the main title (was 5.6 / 10.0 mm) to kill the ~9mm top dead band.
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
# v7 cluster diagnostic: cut the dendrogram into k=4 blocks, report COM range per block (top->bottom display order)
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

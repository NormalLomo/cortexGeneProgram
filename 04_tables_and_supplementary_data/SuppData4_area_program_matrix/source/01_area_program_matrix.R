#!/usr/bin/env Rscript
# =====================================================================
# Extended Data Fig. 5  —  FULL region x ALL-54-program z-heatmap
# Transparency companion to Fig.3 (which shows only the curated
# variable subset). Here: the COMPLETE 14-region x 54-program matrix.
# RENUMBERED 2026-06-20: excluded P9/18/19/35/52/57 removed; remaining
# 54 programs relabelled with new P1-P54 numbers from program_renumber_map.tsv.
# Native R / ComplexHeatmap, single script, vector PDF + PNG.
# =====================================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize)
  library(ggplot2); library(patchwork); library(grid)
  library(RColorBrewer); library(dplyr)
})

set.seed(42)

## ---- paths ----------------------------------------------------------
RES <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
OUT <- "CORTEX_PROGRAM_ROOT/figures/extended"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- shared anatomy / palettes (consistent with Fig.3) --------------
REGION_ORDER <- c("V1","S1","S1E","PoCG","M1","STG","SPL","SMG",
                  "AG","ITG","VLPFC","DLPFC","FPPFC","ACC")
LOBE <- c(V1="Occipital",
          S1="Parietal", S1E="Parietal", PoCG="Parietal",
          SPL="Parietal", SMG="Parietal", AG="Parietal",
          STG="Temporal", ITG="Temporal",
          M1="Frontal/PFC", VLPFC="Frontal/PFC",
          DLPFC="Frontal/PFC", FPPFC="Frontal/PFC",
          ACC="Limbic")
LOBE_ORDER <- c("Occipital","Parietal","Temporal","Frontal/PFC","Limbic")
lobe_pal <- c(Occipital   = "#4C6EB1",
              Parietal    = "#33A089",
              Temporal    = "#E2A22C",
              `Frontal/PFC`= "#C44E52",
              Limbic      = "#8064A2")
class_pal <- c(variable = "#C44E52", stable = "#9AA0A6")
div_pal   <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(256)

## ---- Nature theme (matches Fig.3 theme_nat) -------------------------
theme_nat <- function(base = 6.6) {
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      axis.line   = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.35, colour = "black"),
      axis.title  = element_text(size = base),
      axis.text   = element_text(size = base - 0.6, colour = "black"),
      legend.title= element_text(size = base - 0.4),
      legend.text = element_text(size = base - 1.0),
      legend.key.size = unit(3.2, "mm"),
      plot.title  = element_text(size = base + 1, face = "bold"),
      plot.subtitle = element_text(size = base - 0.6, colour = "grey30"),
      panel.grid  = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}
theme_set(theme_nat())

## ---- load -----------------------------------------------------------
z_df   <- read.table(file.path(RES, "program_region_zscore.tsv"),
                     sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
var_df <- read.table(file.path(RES, "program_variability.tsv"),
                     sep = "\t", header = TRUE, check.names = FALSE)
nm_df  <- read.table(file.path(RES, "program_names.tsv"),
                     sep = "\t", header = TRUE, quote = "", comment.char = "",
                     stringsAsFactors = FALSE, check.names = FALSE)
## ---- RENUMBER 2026-06-20: load mapping, exclude, relabel -----------
rmap <- read.table(file.path(RES, "program_renumber_map.tsv"),
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# old_P (int) -> new_P ("EXCLUDED" or "1".."54")
EXCLUDED_OLD <- as.character(rmap$old_P[rmap$new_P == "EXCLUDED"])  # "9","18","19","35","52","57"
kept_old     <- as.character(rmap$old_P[rmap$new_P != "EXCLUDED"])   # 54 old indices
old2new      <- setNames(rmap$new_P[rmap$new_P != "EXCLUDED"],
                         as.character(rmap$old_P[rmap$new_P != "EXCLUDED"]))

## keep only 54 kept columns, rename to new index (string)
M_z_full <- as.matrix(z_df[, as.character(1:60)])   # full 60-col matrix
M_z      <- M_z_full[, kept_old]                    # 14 x 54 (old col names)
colnames(M_z) <- old2new[kept_old]                  # rename to new P numbers (string)
progs    <- colnames(M_z)                            # "1".."54" (new indices as string)
stopifnot(nrow(M_z) == 14L, ncol(M_z) == 54L)
stopifnot(all(rownames(M_z) %in% REGION_ORDER))

## var_df: filter out excluded, relabel program col to new P number
var_df$program <- as.character(var_df$program)
var_df$class   <- factor(var_df$class, levels = c("variable","stable"))
var_df <- var_df[!var_df$program %in% EXCLUDED_OLD, ]
var_df$program <- old2new[var_df$program]  # relabel to new P number (string)

## nm_df: use cnmf_component (= old P) to build old->new label; keep only kept
nm_df$old_prog <- as.character(nm_df$cnmf_component)
nm_df <- nm_df[!nm_df$old_prog %in% EXCLUDED_OLD, ]
nm_df$program  <- old2new[nm_df$old_prog]   # new P number (string)

## ---- label maker (P{n} short-name; weak-evidence -> trailing *) -----
## 'p' here is the new P number (string "1".."54")
.mk_label <- function(p, ns, conf) {
  ns  <- trimws(ns)
  suf <- paste0(" P", p)
  if (grepl(paste0(suf, "$"), ns)) ns <- trimws(sub(paste0(suf, "$"), "", ns))
  star <- ifelse(conf == "brain-weak", "*", "")
  paste0("P", p, " ", ns, star)
}
prog_lab <- setNames(
  mapply(.mk_label, nm_df$program, nm_df$name_short, nm_df$confidence),
  nm_df$program)
prog_pn  <- setNames(paste0("P", nm_df$program), nm_df$program)
.lab <- function(p) ifelse(p %in% names(prog_lab), prog_lab[p], paste0("P", p))

## =====================================================================
## COLUMN ORDER = programs sorted by variability (eta2_region desc).
## This places the 14 variable programs (high eta2) at the LEFT and the
## 46 stable ones to the RIGHT, giving a transparent eta2 gradient across
## the full matrix without hierarchical column reshuffling.
## =====================================================================
ord_prog  <- var_df %>% arrange(desc(eta2_region)) %>% pull(program)
M         <- M_z[, ord_prog]                       # region x program (eta2-ordered cols)

prog_class <- setNames(as.character(var_df$class[match(ord_prog, var_df$program)]), ord_prog)
prog_eta   <- setNames(var_df$eta2_region[match(ord_prog, var_df$program)], ord_prog)
n_var      <- sum(prog_class == "variable")        # dynamic (was 14 in 60-prog)
n_stab     <- sum(prog_class == "stable")          # dynamic (was 46 in 60-prog)
cat(sprintf("variable / stable count: %d / %d\n", n_var, n_stab))

## column labels: fully label the 14 variable programs (figure focus),
## bare "P{n}" for the 46 stable ones (keeps the dense axis legible while
## every program stays traceable by number). >=4pt rotated 90deg.
col_lab <- ifelse(prog_class[ord_prog] == "variable", .lab(ord_prog), prog_pn[ord_prog])

## column split: variable | stable (slit gap marks the eta2 boundary)
col_split <- factor(ifelse(prog_class[ord_prog] == "variable",
                           sprintf("variable (%d)", n_var), sprintf("stable (%d)", n_stab)),
                    levels = c(sprintf("variable (%d)", n_var), sprintf("stable (%d)", n_stab)))

## =====================================================================
## PANEL a — FULL 14 region x 54 program z-heatmap (square-ish cells)
##   rows = 14 regions, clustered (ward.D2 on 54-program z-vector)
##   cols = 54 programs, eta2-ordered, split variable|stable
##   top annotation: class strip + eta2 barplot
##   left annotation: lobe strip
## Square-ish cells: explicit per-cell width/height so the 14x54 grid
## reads as a clean matrix (cells ~3.0 mm; 54 cols wide banner, but each
## CELL is near-square, satisfying the square-cell rule for heatmaps).
## =====================================================================
CELL <- unit(3.0, "mm")                            # near-square cell edge
top_anno <- HeatmapAnnotation(
  class = prog_class[ord_prog],
  `eta2 (region)` = anno_barplot(
      prog_eta[ord_prog], gp = gpar(fill = "#6E7B8B", col = NA),
      height = unit(9, "mm"),
      axis_param = list(gp = gpar(fontsize = 5))),
  col = list(class = class_pal),
  annotation_name_gp = gpar(fontsize = 6),
  annotation_legend_param = list(class = list(
      title = "class", title_gp = gpar(fontsize = 6),
      labels_gp = gpar(fontsize = 5.5))),
  simple_anno_size = unit(2.6, "mm"),
  show_legend = c(class = TRUE)
)
left_anno <- rowAnnotation(
  lobe = LOBE[rownames(M)],
  col = list(lobe = lobe_pal),
  annotation_name_gp = gpar(fontsize = 6),
  annotation_legend_param = list(lobe = list(
      title = "lobe", title_gp = gpar(fontsize = 6),
      labels_gp = gpar(fontsize = 5.5))),
  simple_anno_size = unit(2.6, "mm")
)
ht <- Heatmap(
  M, name = "z-score",
  col = colorRamp2(seq(-2.5, 2.5, length.out = 256), div_pal),
  width  = CELL * ncol(M),                          # 54 near-square cells
  height = CELL * nrow(M),                          # 14 near-square cells
  top_annotation = top_anno, left_annotation = left_anno,
  cluster_rows = TRUE, cluster_columns = FALSE,     # cols fixed = eta2 order
  column_split = col_split, column_gap = unit(1.6, "mm"),
  cluster_column_slices = FALSE,
  show_row_dend = TRUE, row_dend_width = unit(7, "mm"),
  row_names_gp = gpar(fontsize = 6),
  row_names_side = "left",
  column_labels = col_lab, column_names_gp = gpar(fontsize = 4.2),
  column_names_rot = 90, column_names_side = "bottom",
  column_title_gp = gpar(fontsize = 6.5, fontface = "bold"),
  row_title = "cortical region (clustered)", row_title_gp = gpar(fontsize = 6.5),
  border = TRUE, border_gp = gpar(col = "grey60", lwd = 0.4),
  rect_gp = gpar(col = NA),
  heatmap_legend_param = list(
      title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5.5),
      legend_height = unit(16, "mm"))
)
cat("panel a (full 14x54 heatmap) defined\n")

## =====================================================================
## PANEL b — companion: FULL eta2(region) variability profile, all 54
## programs, on the SAME column axis as the heatmap (eta2-ordered, split
## variable|stable). Vertical bars; variable programs get full
## functional labels (rotated 90deg), stable ones keep "P{n}".
## Placed directly under the heatmap so columns align conceptually.
## =====================================================================
b_df <- var_df %>%
  mutate(program = factor(program, levels = ord_prog)) %>%   # heatmap col order
  arrange(program) %>%
  mutate(xlab = ifelse(class == "variable", .lab(as.character(program)),
                       prog_pn[as.character(program)]),
         xlab = factor(xlab, levels = xlab))
b_split_x <- n_var + 0.5                                     # variable|stable divider
p_b <- ggplot(b_df, aes(x = xlab, y = eta2_region, fill = class)) +
  geom_col(width = 0.78, colour = NA) +
  geom_vline(xintercept = b_split_x, linetype = "22",
             linewidth = 0.3, colour = "grey45") +
  annotate("text", x = b_split_x - 0.4, y = max(b_df$eta2_region)*0.98,
           label = sprintf("variable (%d)", n_var), hjust = 1, vjust = 1, size = 1.9,
           colour = class_pal["variable"], fontface = "bold") +
  annotate("text", x = b_split_x + 0.4, y = max(b_df$eta2_region)*0.98,
           label = sprintf("stable (%d)", n_stab), hjust = 0, vjust = 1, size = 1.9,
           colour = class_pal["stable"], fontface = "bold") +
  scale_fill_manual(values = class_pal, name = "class") +
  scale_y_continuous(name = expression(eta^2~"(region)"),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(expand = expansion(add = c(0.6, 0.6))) +
  labs(x = NULL, title = NULL) +
  theme_nat() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 4.2, colour = "black"),
        axis.ticks.x = element_line(linewidth = 0.2),
        legend.position = c(0.93, 0.78),
        legend.key.size = unit(2.6, "mm"),
        plot.margin = margin(t = 2, r = 6, b = 2, l = 6, unit = "pt"))
gb_b <- ggplotGrob(p_b)
cat("panel b (full eta2 ranking, heatmap col order) done\n")

## =====================================================================
## ASSEMBLE — NATIVE grid two-row layout (no patchwork outer stack, which
## collapses a fixed-mm ComplexHeatmap grob). Top row = the heatmap drawn
## by ComplexHeatmap into its own viewport; bottom row = the eta2 ggplot
## grob. Page width = heatmap's true content width so there is zero dead
## horizontal whitespace; row heights sized to each panel's content.
## =====================================================================
# true content sizes (mm) -> inches
HM_BODY_W_MM <- 54*3.0 + 1.6                 # 54 cells @3mm + split gap
HM_BODY_H_MM <- 14*3.0                       # 14 cells @3mm
# heatmap full footprint: body + left(rowname ~14 + dend 7) + top(class 2.6 +
# eta2 bar 9 + titles ~10) + bottom(90deg labels ~42) + right legends ~26
HM_W_IN <- (HM_BODY_W_MM + 14 + 7 + 26) / 25.4
HM_H_IN <- (HM_BODY_H_MM + 2.6 + 9 + 12 + 44) / 25.4
# eta2 bar panel: same body width; height for bars + rotated labels
BAR_H_IN <- 2.7
PAD_IN   <- 0.30                              # outer margin
TITLE_IN <- 0.62                              # top title/subtitle band

FIG_W <- HM_W_IN + 2*PAD_IN
FIG_H <- TITLE_IN + HM_H_IN + BAR_H_IN + 2*PAD_IN + 0.25

png_path <- file.path(OUT, "ed_fig5_full_region_heatmap.png")
pdf_path <- file.path(OUT, "ed_fig5_full_region_heatmap.pdf")

draw_fig <- function() {
  grid.newpage()
  # outer master viewport with padding
  pushViewport(viewport(x = unit(PAD_IN, "in"), y = unit(PAD_IN, "in"),
                        width  = unit(FIG_W - 2*PAD_IN, "in"),
                        height = unit(FIG_H - 2*PAD_IN, "in"),
                        just = c("left","bottom")))
  # vertical layout: [title band] / [heatmap] / [gap] / [eta2 bar]
  lay <- grid.layout(nrow = 4, ncol = 1,
                     heights = unit.c(unit(TITLE_IN, "in"),
                                      unit(HM_H_IN, "in"),
                                      unit(0.18, "in"),
                                      unit(BAR_H_IN, "in")))
  pushViewport(viewport(layout = lay))

  ## --- title band (row 1) ---
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text("Extended Data Fig. 5  Full region × program activity matrix (all 54 cNMF programs)",
            x = unit(0, "npc"), y = unit(0.72, "npc"), just = c("left","center"),
            gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Helvetica"))
  grid.text(paste0("Complete 14-region × 54-program z-score matrix and full ",
                   "variability ranking; Fig.3 shows only the curated variable subset."),
            x = unit(0, "npc"), y = unit(0.30, "npc"), just = c("left","center"),
            gp = gpar(fontsize = 7.2, col = "grey30", fontfamily = "Helvetica"))
  popViewport()

  ## --- panel a: heatmap (row 2) ---
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.text("a", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left","top"),
            gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "Helvetica"))
  # nested viewport offset right so the panel tag 'a' is not overdrawn
  pushViewport(viewport(x = unit(2.2, "mm"), width = unit(1, "npc") - unit(2.2, "mm"),
                        just = "left"))
  draw(ht, newpage = FALSE,
       heatmap_legend_side = "right", annotation_legend_side = "right",
       merge_legend = TRUE,
       column_title = "All 54 cortical gene programs (cNMF K=60, 6 cohort-technical excluded), ordered by cross-region variability (η²)",
       column_title_gp = gpar(fontsize = 7, fontface = "bold"))
  popViewport(2)

  ## --- panel b: eta2 bar (row 4) ---
  pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
  grid.text("b", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left","top"),
            gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "Helvetica"))
  grid.text("Cross-region variability of all 54 programs (η² ranking)",
            x = unit(3, "mm"), y = unit(1, "npc") - unit(0.5, "mm"),
            just = c("left","top"),
            gp = gpar(fontsize = 7.6, fontface = "bold", fontfamily = "Helvetica"))
  pushViewport(viewport(y = 0, height = unit(1, "npc") - unit(5, "mm"), just = "bottom"))
  grid.draw(gb_b)
  popViewport(2)

  popViewport(2)  # layout + master
}

# cairo devices: render the Unicode η (U+03B7) glyph correctly AND embed
# fonts as subsetted Type1C in the PDF (base pdf() fails on mbcs->sbcs for η).
cairo_pdf(pdf_path, width = FIG_W, height = FIG_H, family = "Helvetica")
draw_fig(); dev.off()
png(png_path, width = FIG_W, height = FIG_H, units = "in", res = 220,
    bg = "white", type = "cairo")
draw_fig(); dev.off()
cat(sprintf("assembly done  [page %.2f x %.2f in]\n", FIG_W, FIG_H))

## ---- report numbers -------------------------------------------------
cat("\n=== ED FIG 5 SUMMARY ===\n")
cat(sprintf("matrix dims        : %d regions x %d programs\n", nrow(M), ncol(M)))
cat(sprintf("column order       : eta2(region) descending\n"))
cat(sprintf("row order          : ward.D2 clustering of regions on 54-program z\n"))
cat(sprintf("variable / stable  : %d / %d (column split shown)\n", n_var, n_stab))
cat(sprintf("page (in)          : %.1f x %.1f\n", FIG_W, FIG_H))
hc <- hclust(dist(M, method = "euclidean"), method = "ward.D2")
cat("clustered region order :", paste(rownames(M)[hc$order], collapse = " "), "\n")
cat("top-5 variable (eta2)  :", paste(head(ord_prog, 5), collapse = " "), "\n")
cat("\nALL DONE\n")

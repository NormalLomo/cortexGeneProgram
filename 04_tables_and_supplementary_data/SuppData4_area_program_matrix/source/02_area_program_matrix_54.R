#!/usr/bin/env Rscript
# =====================================================================
# Extended Data Fig. 5 (54-program version) — FULL region x 54-program z-heatmap
# Removes 6 excluded cohort-technical programs from the original 60.
# Excluded old-component ids: 9, 18, 19, 35, 52, 57
# All other settings identical to original ed_fig5_full_region_heatmap.R
# Outputs: figures/extended/ed_fig5_full_region_heatmap_54.{pdf,png}
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

## ---- excluded component ids (old numbering) -------------------------
EXCLUDED_OLD <- c(9L, 18L, 19L, 35L, 52L, 57L)

## ---- load renumber map to get new_P labels --------------------------
remap <- read.table(file.path(RES, "program_renumber_map.tsv"),
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# kept rows only; new_P is e.g. "1","2"... (numeric strings) for kept
kept_map <- remap[remap$status == "kept", ]
stopifnot(nrow(kept_map) == 54L)
# old_P (integer) -> new_P (integer string)
old_to_new <- setNames(as.character(kept_map$new_P), as.character(kept_map$old_P))

## ---- shared anatomy / palettes (consistent with Fig.2) --------------
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

## ---- Nature theme ---------------------------------------------------
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

# all 60 programs (old numbering as integers)
all_progs_old <- as.character(1:60)
kept_progs_old <- all_progs_old[!(all_progs_old %in% as.character(EXCLUDED_OLD))]
stopifnot(length(kept_progs_old) == 54L)

# filter matrix to 54 kept columns
M_z_full <- as.matrix(z_df[, all_progs_old])
stopifnot(nrow(M_z_full) == 14L, ncol(M_z_full) == 60L)
M_z <- M_z_full[, kept_progs_old]   # 14 x 54
stopifnot(ncol(M_z) == 54L)
cat("Kept programs (old ids):", paste(kept_progs_old, collapse=" "), "\n")

# filter variability table to 54
var_df$program <- as.character(var_df$program)
var_df_kept <- var_df[var_df$program %in% kept_progs_old, ]
var_df_kept$class <- factor(var_df_kept$class, levels = c("variable","stable"))
stopifnot(nrow(var_df_kept) == 54L)

# program names: nm_df has new_P and cnmf_component columns
nm_df$cnmf_component <- as.character(nm_df$cnmf_component)
nm_df_kept <- nm_df[nm_df$new_P != "EXCLUDED", ]
stopifnot(nrow(nm_df_kept) == 54L)
# map: old integer id -> new_P label
nm_df_kept$old_id <- nm_df_kept$cnmf_component

## ---- label maker using new_P numbers --------------------------------
.mk_label <- function(old_id, ns, conf) {
  # old_id: character integer like "1", "10", etc.
  new_p <- old_to_new[old_id]   # new P number string e.g. "1", "54"
  ns    <- trimws(ns)
  # remove trailing " P{n}" if present
  suf   <- paste0(" P", new_p)
  if (grepl(paste0(suf, "$"), ns)) ns <- trimws(sub(paste0(suf, "$"), "", ns))
  star  <- ifelse(conf == "brain-weak", "*", "")
  paste0("P", new_p, " ", ns, star)
}
.mk_pn <- function(old_id) {
  new_p <- old_to_new[old_id]
  paste0("P", new_p)
}
prog_lab <- setNames(
  mapply(.mk_label, nm_df_kept$old_id, nm_df_kept$name_short, nm_df_kept$confidence),
  nm_df_kept$old_id)
prog_pn <- setNames(
  sapply(nm_df_kept$old_id, .mk_pn),
  nm_df_kept$old_id)
.lab <- function(p) ifelse(p %in% names(prog_lab), prog_lab[p], paste0("P", p))

## ---- column order = variability descending --------------------------
ord_prog <- var_df_kept %>% arrange(desc(eta2_region)) %>% pull(program)
M        <- M_z[, ord_prog]   # 14 x 54, eta2-ordered columns
stopifnot(ncol(M) == 54L, nrow(M) == 14L, all(rownames(M) %in% REGION_ORDER))

prog_class <- setNames(as.character(var_df_kept$class[match(ord_prog, var_df_kept$program)]), ord_prog)
prog_eta   <- setNames(var_df_kept$eta2_region[match(ord_prog, var_df_kept$program)], ord_prog)
n_var      <- sum(prog_class == "variable")
n_stab     <- sum(prog_class == "stable")
cat(sprintf("variable=%d, stable=%d (total=%d)\n", n_var, n_stab, n_var + n_stab))

# column labels: fully label variable programs, bare P{new} for stable
col_lab <- ifelse(prog_class[ord_prog] == "variable", .lab(ord_prog), prog_pn[ord_prog])
col_split <- factor(ifelse(prog_class[ord_prog] == "variable",
                           sprintf("variable (%d)", n_var),
                           sprintf("stable (%d)", n_stab)),
                    levels = c(sprintf("variable (%d)", n_var),
                               sprintf("stable (%d)", n_stab)))

## ---- panel a: heatmap -----------------------------------------------
CELL <- unit(3.0, "mm")
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
  width  = CELL * ncol(M),   # 54 cells
  height = CELL * nrow(M),   # 14 cells
  top_annotation = top_anno, left_annotation = left_anno,
  cluster_rows = TRUE, cluster_columns = FALSE,
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

## ---- panel b: eta2 bar ----------------------------------------------
b_df <- var_df_kept %>%
  mutate(program = factor(program, levels = ord_prog)) %>%
  arrange(program) %>%
  mutate(new_label = ifelse(class == "variable", .lab(as.character(program)),
                            prog_pn[as.character(program)]),
         new_label = factor(new_label, levels = new_label))
b_split_x <- n_var + 0.5
p_b <- ggplot(b_df, aes(x = new_label, y = eta2_region, fill = class)) +
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

## ---- assemble -------------------------------------------------------
HM_BODY_W_MM <- 54*3.0 + 1.6   # 54 cells @3mm + split gap
HM_BODY_H_MM <- 14*3.0
HM_W_IN <- (HM_BODY_W_MM + 14 + 7 + 26) / 25.4
HM_H_IN <- (HM_BODY_H_MM + 2.6 + 9 + 12 + 44) / 25.4
BAR_H_IN <- 2.7
PAD_IN   <- 0.30
TITLE_IN <- 0.62

FIG_W <- HM_W_IN + 2*PAD_IN
FIG_H <- TITLE_IN + HM_H_IN + BAR_H_IN + 2*PAD_IN + 0.25

png_path <- file.path(OUT, "ed_fig5_full_region_heatmap_54.png")
pdf_path <- file.path(OUT, "ed_fig5_full_region_heatmap_54.pdf")

draw_fig <- function() {
  grid.newpage()
  pushViewport(viewport(x = unit(PAD_IN, "in"), y = unit(PAD_IN, "in"),
                        width  = unit(FIG_W - 2*PAD_IN, "in"),
                        height = unit(FIG_H - 2*PAD_IN, "in"),
                        just = c("left","bottom")))
  lay <- grid.layout(nrow = 4, ncol = 1,
                     heights = unit.c(unit(TITLE_IN, "in"),
                                      unit(HM_H_IN, "in"),
                                      unit(0.18, "in"),
                                      unit(BAR_H_IN, "in")))
  pushViewport(viewport(layout = lay))

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text("Fig. S6 |  Full region × program activity matrix (54 cNMF programs)",
            x = unit(0, "npc"), y = unit(0.72, "npc"), just = c("left","center"),
            gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Helvetica"))
  grid.text(paste0("Complete 14-region x 54-program z-score matrix and variability ranking. ",
                   "Six cohort-technical programs excluded (Methods). ",
                   "Fig.2 shows only the curated variable subset."),
            x = unit(0, "npc"), y = unit(0.30, "npc"), just = c("left","center"),
            gp = gpar(fontsize = 7.2, col = "grey30", fontfamily = "Helvetica"))
  popViewport()

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.text("a", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left","top"),
            gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "Helvetica"))
  pushViewport(viewport(x = unit(2.2, "mm"), width = unit(1, "npc") - unit(2.2, "mm"),
                        just = "left"))
  draw(ht, newpage = FALSE,
       heatmap_legend_side = "right", annotation_legend_side = "right",
       merge_legend = TRUE,
       column_title = sprintf("54 cortical gene programs (cNMF K=60, 6 cohort-technical excluded), ordered by cross-region variability (%s²)",
                              "η"),
       column_title_gp = gpar(fontsize = 7, fontface = "bold"))
  popViewport(2)

  pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
  grid.text("b", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left","top"),
            gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "Helvetica"))
  grid.text(sprintf("Cross-region variability of 54 programs (%s² ranking)", "η"),
            x = unit(3, "mm"), y = unit(1, "npc") - unit(0.5, "mm"),
            just = c("left","top"),
            gp = gpar(fontsize = 7.6, fontface = "bold", fontfamily = "Helvetica"))
  pushViewport(viewport(y = 0, height = unit(1, "npc") - unit(5, "mm"), just = "bottom"))
  grid.draw(gb_b)
  popViewport(2)

  popViewport(2)
}

cairo_pdf(pdf_path, width = FIG_W, height = FIG_H, family = "Helvetica")
draw_fig(); dev.off()
png(png_path, width = FIG_W, height = FIG_H, units = "in", res = 220,
    bg = "white", type = "cairo")
draw_fig(); dev.off()
cat(sprintf("assembly done  [page %.2f x %.2f in]\n", FIG_W, FIG_H))

cat("\n=== SUMMARY ===\n")
cat(sprintf("matrix dims: %d regions x %d programs\n", nrow(M), ncol(M)))
cat(sprintf("excluded old-component ids: %s\n", paste(sort(EXCLUDED_OLD), collapse=",")))
cat(sprintf("variable / stable: %d / %d\n", n_var, n_stab))
cat(sprintf("page (in): %.1f x %.1f\n", FIG_W, FIG_H))
hc <- hclust(dist(M, method = "euclidean"), method = "ward.D2")
cat("clustered region order:", paste(rownames(M)[hc$order], collapse=" "), "\n")
cat("top-5 variable (eta2):", paste(head(ord_prog, 5), collapse=" "), "\n")
cat("\nALL DONE\n")

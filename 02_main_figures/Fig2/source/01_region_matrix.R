#!/usr/bin/env Rscript
# =====================================================================
# Fig.3  Cross-region gene-program variation  (human cortex, cNMF K=60)
# 9 panels (a-i). Submission-grade, Nature-family style.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ComplexHeatmap)
  library(circlize); library(ggalluvial); library(ggrepel)
  library(dplyr); library(tidyr); library(grid)
  library(scales); library(RColorBrewer)
})

set.seed(42)

## ---- paths ----------------------------------------------------------
RES <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
OUT <- "CORTEX_PROGRAM_ROOT/figures/fig3"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- fixed anatomical region order & lobes --------------------------
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

## ---- palettes (unified family across panels) ------------------------
lobe_pal <- c(Occipital   = "#4C6EB1",   # blue
              Parietal    = "#33A089",   # teal-green
              Temporal    = "#E2A22C",   # amber
              `Frontal/PFC`= "#C44E52",  # red
              Limbic      = "#8064A2")   # purple
class_pal <- c(variable = "#C44E52", stable = "#9AA0A6")
# diverging z palette
div_pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(256)

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
      strip.text  = element_text(size = base - 0.2, face = "bold"),
      strip.background = element_blank(),
      plot.title  = element_text(size = base + 1, face = "bold"),
      plot.subtitle = element_text(size = base - 0.6, colour = "grey30"),
      panel.grid  = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}
theme_set(theme_nat())

## ---- load data ------------------------------------------------------
mean_df <- read.table(file.path(RES, "region_program_mean.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
z_df    <- read.table(file.path(RES, "program_region_zscore.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
var_df  <- read.table(file.path(RES, "program_variability.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE)
grad_pre<- read.table(file.path(RES, "program_gradient.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE)
gsub_df <- read.table(file.path(RES, "panel_g_subsample.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE)

progs <- as.character(1:60)
# matrices: regions (rows) x programs (cols). transpose to programs x regions
M_mean <- as.matrix(mean_df[, progs])            # region x program
M_z    <- as.matrix(z_df[,    progs])            # region x program
stopifnot(all(rownames(M_z) %in% REGION_ORDER))

var_df$program <- as.character(var_df$program)
var_df$class   <- factor(var_df$class, levels = c("variable","stable"))
variable_progs <- var_df$program[var_df$class == "variable"]

# region-level mean activity & CV per program (across 14 regions)
prog_mean_act <- colMeans(M_mean)                       # mean activity across regions
prog_cv       <- apply(M_mean, 2, function(x) sd(x)/mean(x))

## ---- program -> functional label lookup -----------------------------
## display label = "P{n} {name_short}"; strip redundant trailing " P{n}"
## suffix (disambiguation artifact); brain-weak programs get trailing "*".
nm_df <- read.table(file.path(RES, "program_names.tsv"),
                    sep = "\t", header = TRUE, quote = "", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE)
# program_names.tsv has columns (new_P, cnmf_component, name_short,
# confidence, ...). All upstream data tables (M_z / M_mean / var_df) key on the LEGACY
# cnmf_component string (1..60). So our lookups must be keyed by cnmf_component but DISPLAY
# the new_P tag ("P{new_int}").
nm_df$cnmf_component <- as.character(nm_df$cnmf_component)
nm_df$new_P <- as.character(nm_df$new_P)
# drop EXCLUDED entries; keep only the 54 biologically-interpreted programs in the lookup
nm_df_keep <- nm_df[toupper(nm_df$new_P) != "EXCLUDED", ]
.mk_label_current <- function(newP, ns, conf) {
  ns  <- trimws(ns)
  suf <- paste0(" ", newP)
  if (grepl(paste0(suf, "$"), ns)) ns <- trimws(sub(paste0(suf, "$"), "", ns))
  star <- ifelse(conf == "brain-weak", "*", "")
  paste0(newP, " ", ns, star)
}
# keyed by LEGACY cnmf_component string ("1", "14", ...) -> "P{new_int} name_short[*]"
prog_lab <- setNames(
  mapply(.mk_label_current, nm_df_keep$new_P, nm_df_keep$name_short, nm_df_keep$confidence),
  nm_df_keep$cnmf_component)
# bare new_P tag (used for non-highlighted items in dense panels + alluvial strata)
prog_pn  <- setNames(nm_df_keep$new_P, nm_df_keep$cnmf_component)
# safety: any cnmf id missing from name map falls back to legacy "P{cnmf}"
.lab <- function(p) ifelse(p %in% names(prog_lab), prog_lab[p], paste0("P", p))

## =====================================================================
## helper to compute gradient (r & slope) of region mean vs anat rank
## =====================================================================
anat_rank <- setNames(seq_along(REGION_ORDER), REGION_ORDER)
grad_tab <- bind_rows(lapply(progs, function(p) {
  v <- M_mean[REGION_ORDER, p]
  ct <- suppressWarnings(cor.test(anat_rank[REGION_ORDER], v, method = "pearson"))
  sl <- coef(lm(v ~ anat_rank[REGION_ORDER]))[2]
  data.frame(program = p, grad_r = unname(ct$estimate),
             grad_slope = unname(sl), grad_p = ct$p.value)
}))

## =====================================================================
## PANEL a — ComplexHeatmap: region x program activity (z), clustered
## =====================================================================
# matrix: regions (rows) x programs (cols); color = z-score
Ha <- t(M_z[REGION_ORDER, ])                  # program x region NOT used; keep region rows
Ha <- M_z[REGION_ORDER, progs]                # region x program
# annotations
prog_class <- setNames(as.character(var_df$class[match(progs, var_df$program)]), progs)
prog_eta   <- setNames(var_df$eta2_region[match(progs, var_df$program)], progs)
# column labels: 60 programs is too dense for full names on every column,
# so fully label only the *variable* programs (the figure's focus) and keep
# stable ones as bare "P{n}". Rotated 90deg, small font (legibility rule).
col_lab_a <- ifelse(prog_class[progs] == "variable", .lab(progs), prog_pn[progs])

top_anno <- HeatmapAnnotation(
  class = prog_class,
  eta2  = anno_barplot(prog_eta, gp = gpar(fill = "#6E7B8B", col = NA),
                       height = unit(7, "mm"), axis_param = list(gp = gpar(fontsize = 5))),
  col = list(class = class_pal),
  annotation_name_gp = gpar(fontsize = 6),
  annotation_legend_param = list(class = list(title = "class",
                                  title_gp = gpar(fontsize = 6),
                                  labels_gp = gpar(fontsize = 5.5))),
  simple_anno_size = unit(2.5, "mm"),
  show_legend = c(class = TRUE)
)
left_anno <- rowAnnotation(
  lobe = LOBE[REGION_ORDER],
  col = list(lobe = lobe_pal),
  annotation_name_gp = gpar(fontsize = 6),
  annotation_legend_param = list(lobe = list(title = "lobe",
                                  title_gp = gpar(fontsize = 6),
                                  labels_gp = gpar(fontsize = 5.5))),
  simple_anno_size = unit(2.5, "mm")
)
ht_a <- Heatmap(
  Ha, name = "z-score",
  col = colorRamp2(seq(-2.5, 2.5, length.out = 256), div_pal),
  top_annotation = top_anno, left_annotation = left_anno,
  cluster_rows = TRUE, cluster_columns = TRUE,
  show_row_dend = TRUE, show_column_dend = TRUE,
  row_dend_width = unit(6, "mm"), column_dend_height = unit(6, "mm"),
  row_names_gp = gpar(fontsize = 6),
  column_labels = col_lab_a, column_names_gp = gpar(fontsize = 4),
  column_names_rot = 90,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                              labels_gp = gpar(fontsize = 5.5),
                              legend_height = unit(15, "mm")),
  row_title = "cortical region", column_title = "program (cNMF K=60)",
  row_title_gp = gpar(fontsize = 6.5), column_title_gp = gpar(fontsize = 6.5)
)
# Panel a was too flat (~2.3:1 banner). Export taller (7.0 x 4.3 in, ~1.6:1)
# so the 14 region rows get taller cells. ComplexHeatmap annotation tracks
# (col dendrogram 6mm, class strip 2.5mm, eta2 barplot 7mm, lobe strip 2.5mm,
# row dendrogram 6mm) are FIXED in mm, so the extra ~1.3 in of canvas height
# is absorbed entirely by the heatmap BODY -> the 14 z-score rows grow taller
# (not the text/tracks, which stay at their fixed mm sizes & fontsizes).
PANEL_A_W <- 7.0; PANEL_A_H <- 4.3
pdf(file.path(OUT, "panel_a.pdf"), width = PANEL_A_W, height = PANEL_A_H, family = "Helvetica")
draw(ht_a, heatmap_legend_side = "right", annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()
# capture as grob for assembly
gb_a <- grid.grabExpr(draw(ht_a, heatmap_legend_side = "right",
                           annotation_legend_side = "right", merge_legend = TRUE),
                      width = PANEL_A_W, height = PANEL_A_H)
p_a <- patchwork::wrap_elements(full = gb_a)
cat("panel a done\n")

## =====================================================================
## PANEL b — variability ranking lollipop (eta2_region desc)
## =====================================================================
b_df <- var_df %>% arrange(desc(eta2_region)) %>%
  mutate(rank = row_number(),
         program = factor(program, levels = rev(program)),
         ypos = as.integer(program))            # 1..60; rank-1 -> ypos 60 (top)
# Label only the top-8 variable programs. The crowded y-axis tick labels (and
# earlier auto-repel attempts) collided because the top ranks bunch near the
# top of the axis. Instead of letting an optimiser fight that geometry, we
# DETERMINISTICALLY assign each of the 8 callouts a fixed, evenly-spaced y in
# the right margin and draw an explicit leader line from each point -> its
# label. Even spacing across the panel height => zero overlap, guaranteed.
b_xmax <- max(b_df$eta2_region)
nlab   <- 8L
b_lab  <- b_df %>% filter(rank <= nlab) %>% arrange(rank) %>%
  mutate(lab   = .lab(as.character(program)),
         lab_x = b_xmax * 1.05,                 # common x for all label text
         # spread label-y evenly from high (top) to low across the panel so
         # adjacent labels never touch; order follows rank (top label = rank 1)
         lab_y = seq(58, 12, length.out = nlab))
p_b <- ggplot(b_df, aes(x = eta2_region, y = ypos, colour = class)) +
  geom_segment(aes(x = 0, xend = eta2_region, yend = ypos), linewidth = 0.4) +
  geom_point(size = 1.4) +
  # leader lines: from each labelled point to its fixed label anchor
  geom_segment(data = b_lab, inherit.aes = FALSE,
               aes(x = eta2_region, y = ypos, xend = lab_x, yend = lab_y),
               colour = "grey55", linewidth = 0.18) +
  geom_text(data = b_lab, inherit.aes = FALSE,
            aes(x = lab_x, y = lab_y, label = lab),
            hjust = 0, size = 1.8, colour = "black") +
  scale_colour_manual(values = class_pal, name = "class") +
  # wide right expansion so the full label text (longest ~30 chars) fits in
  # the margin without clipping; labels start at lab_x = b_xmax * 1.05.
  scale_x_continuous(expand = expansion(mult = c(0, 1.9)),
                     name = expression(eta^2~"(region)")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(y = NULL, title = "Program variability ranking") +
  theme_nat() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.82, 0.30),
        # 1-col cell in assembled layout is narrow & tall -> force a tall
        # aspect so the lollipop content stays legible (not squashed wide).
        aspect.ratio = 1.6)
ggsave(file.path(OUT, "panel_b.pdf"), p_b, width = 3.0, height = 3.6, device = cairo_pdf)
cat("panel b done\n")

## =====================================================================
## PANEL c — mean activity vs CV scatter + marginal densities
## =====================================================================
c_df <- data.frame(program = progs,
                   mean_act = prog_mean_act[progs],
                   cv = prog_cv[progs],
                   class = var_df$class[match(progs, var_df$program)],
                   eta2 = prog_eta[progs])
top_lab <- var_df %>% arrange(desc(eta2_region)) %>% slice_head(n = 5) %>% pull(program)
# label only the top-5 most variable programs; stronger repulsion keeps the
# few callouts from overlapping each other or the dense point cloud.
c_df$lab <- ifelse(c_df$program %in% top_lab, .lab(c_df$program), "")

p_c_main <- ggplot(c_df, aes(mean_act, cv, colour = class)) +
  geom_point(aes(size = eta2), alpha = 0.85) +
  geom_text_repel(aes(label = lab), size = 1.9, colour = "black",
                  max.overlaps = Inf, segment.size = 0.2, min.segment.length = 0,
                  box.padding = 0.55, point.padding = 0.3, force = 4,
                  seed = 42) +
  scale_colour_manual(values = class_pal, name = "class", guide = "none") +
  scale_size_continuous(range = c(0.4, 2.6), name = expression(eta^2)) +
  scale_x_continuous(name = "mean activity (across regions)") +
  scale_y_continuous(name = "CV across regions") +
  theme_nat() + theme(legend.position = c(0.85, 0.82))
p_c_top <- ggplot(c_df, aes(mean_act, fill = class, colour = class)) +
  geom_density(alpha = 0.35, linewidth = 0.3) +
  scale_fill_manual(values = class_pal, guide = "none") +
  scale_colour_manual(values = class_pal, guide = "none") +
  theme_void()
p_c_right <- ggplot(c_df, aes(cv, fill = class, colour = class)) +
  geom_density(alpha = 0.35, linewidth = 0.3) +
  scale_fill_manual(values = class_pal, guide = "none") +
  scale_colour_manual(values = class_pal, guide = "none") +
  coord_flip() + theme_void()
p_c_grp <- (p_c_top + plot_spacer() + p_c_main + p_c_right) +
  plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4))
ggsave(file.path(OUT, "panel_c.pdf"), p_c_grp, width = 3.0, height = 3.0, device = cairo_pdf)
# wrap as single element so patchwork tags it once in the assembled figure
p_c <- wrap_elements(full = p_c_grp)
cat("panel c done\n")

## =====================================================================
## PANEL d — REPLACEMENT: gradient along PROGRAM-CLUSTERED region axis
## (region order = ward.D2 / Euclidean clustering of regions on their
##  60-program z-vector; replaces invalid PCA-PC1 cohort-confounded axis)
## =====================================================================
# cluster the 14 regions on their 60-program z-vector
d_reg_dist <- dist(M_z[rownames(M_z), progs], method = "euclidean")
hc_reg     <- hclust(d_reg_dist, method = "ward.D2")
leaf_order <- hc_reg$labels[hc_reg$order]
axis_rank  <- setNames(seq_along(leaf_order), leaf_order)

# persist region cluster order + clustered-axis gradient (results dir)
reg_order_tab <- data.frame(region = leaf_order,
                            leaf_order = seq_along(leaf_order),
                            axis_rank = as.integer(axis_rank[leaf_order]))
write.table(reg_order_tab, file.path(RES, "region_cluster_order.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
grad_clust <- bind_rows(lapply(progs, function(p) {
  v  <- M_mean[leaf_order, p]; x <- axis_rank[leaf_order]
  ct <- suppressWarnings(cor.test(x, v, method = "pearson"))
  sl <- coef(lm(v ~ x))[2]
  data.frame(program = p, axis_slope = unname(sl),
             axis_r = unname(ct$estimate), axis_p = ct$p.value)
}))
write.table(grad_clust, file.path(RES, "program_gradient_clustered.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# pick variable programs with strongest gradient on the clustered axis
var_grad <- grad_clust %>% filter(program %in% variable_progs) %>%
  arrange(desc(abs(axis_r)))
nshow <- min(8L, nrow(var_grad))
show_progs <- var_grad$program[seq_len(nshow)]

d_long <- as.data.frame(M_z[leaf_order, show_progs]) %>%
  mutate(region = factor(leaf_order, levels = leaf_order)) %>%
  pivot_longer(-region, names_to = "program", values_to = "z")
d_grad <- grad_clust %>% filter(program %in% show_progs)
d_lab  <- d_long %>% group_by(program) %>% slice_max(as.integer(region), n = 1) %>%
  left_join(d_grad, by = "program") %>%
  # full functional label + gradient r on a second line (keeps horizontal
  # extent small so the 8 line-end labels fit in the right margin)
  mutate(txt = sprintf("%s\n(r=%.2f)", .lab(program), axis_r))
ord_progs  <- d_grad$program[order(-d_grad$axis_r)]
prog_cols8 <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(show_progs)),
                       ord_progs)
# lobe color strip under axis (anatomy reference along the clustered order)
lobe_strip <- data.frame(region = factor(leaf_order, levels = leaf_order),
                         lobe = factor(LOBE[leaf_order], levels = LOBE_ORDER),
                         x = seq_along(leaf_order))
# push the lobe strip further below the lowest trace (extra gap above strip)
ystrip <- min(c(d_long$z, 0)) - 0.85; ystep <- 0.32

p_d <- ggplot(d_long, aes(region, z, group = program, colour = program)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
  geom_line(linewidth = 0.55) +
  geom_point(size = 0.7) +
  geom_tile(data = lobe_strip, inherit.aes = FALSE,
            aes(x = x, y = ystrip, fill = lobe),
            width = 1, height = ystep, colour = "white", linewidth = 0.2) +
  geom_text_repel(data = d_lab, aes(label = txt), hjust = 0, direction = "y",
                  nudge_x = 0.4, size = 1.55, segment.size = 0.2, lineheight = 0.85,
                  xlim = c(14.3, NA), max.overlaps = Inf, force = 3,
                  box.padding = 0.4, point.padding = 0.2, seed = 42) +
  scale_colour_manual(values = prog_cols8, guide = "none") +
  scale_fill_manual(values = lobe_pal, name = "lobe",
                    guide = guide_legend(override.aes = list(colour = NA))) +
  scale_x_discrete(expand = expansion(mult = c(0.04, 0.62))) +
  scale_y_continuous(expand = expansion(mult = c(0.16, 0.06))) +
  labs(x = "region (program-clustered order)", y = "region z-score",
       title = "Program gradients along clustered region axis",
       subtitle = "axis = ward.D2 clustering of regions on 60-program profile") +
  theme_nat() +
  # lobe legend moved OUT of the plotting area (to the right margin) so it
  # no longer collides with the line traces in the lower-left corner.
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
        legend.position = "right",
        legend.key.size = unit(2.4, "mm"))
ggsave(file.path(OUT, "panel_d.pdf"), p_d, width = 4.0, height = 3.2, device = cairo_pdf)
file.copy(file.path(OUT, "panel_d.pdf"), file.path(OUT, "fig3_d.pdf"), overwrite = TRUE)
cat("panel d (clustered) done\n")

## =====================================================================
## PANEL e — region similarity biplot (PCA of region x program z)
## =====================================================================
pca <- prcomp(M_z[REGION_ORDER, progs], center = TRUE, scale. = FALSE)
ve  <- (pca$sdev^2) / sum(pca$sdev^2) * 100
scores <- as.data.frame(pca$x[, 1:2]); scores$region <- rownames(scores)
scores$lobe <- factor(LOBE[scores$region], levels = LOBE_ORDER)
# loadings for top variable programs (arrows)
load_progs <- var_df %>% arrange(desc(eta2_region)) %>% slice_head(n = 6) %>% pull(program)
ld <- as.data.frame(pca$rotation[load_progs, 1:2]); ld$program <- rownames(ld)
sc_fac <- 0.55 * max(abs(scores$PC1), abs(scores$PC2)) / max(abs(ld$PC1), abs(ld$PC2))
ld$PC1s <- ld$PC1 * sc_fac; ld$PC2s <- ld$PC2 * sc_fac
ld$lab  <- .lab(ld$program)            # full functional label for loading arrows

p_e <- ggplot(scores, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
  geom_segment(data = ld, aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
               arrow = arrow(length = unit(1.4, "mm")), colour = "grey45",
               linewidth = 0.3, inherit.aes = FALSE) +
  geom_text_repel(data = ld, aes(PC1s, PC2s, label = lab), colour = "grey25",
            size = 1.55, fontface = "italic", segment.size = 0.15,
            segment.colour = "grey60", min.segment.length = 0,
            max.overlaps = 30, box.padding = 0.2, inherit.aes = FALSE) +
  geom_point(aes(fill = lobe), size = 2.1, shape = 21, colour = "white", stroke = 0.3) +
  geom_text_repel(aes(label = region, colour = lobe), size = 2.0,
                  max.overlaps = 20, segment.size = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = lobe_pal, name = "lobe") +
  scale_colour_manual(values = lobe_pal, guide = "none") +
  labs(x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]),
       title = "Region similarity (program profile)") +
  coord_cartesian(clip = "off") +
  theme_nat() + theme(legend.position = c(0.86, 0.22),
                      plot.margin = margin(t = 1, r = 2, b = 6, l = 2, unit = "pt"),
                      # 1-col cell -> tall aspect keeps the PCA biplot legible.
                      aspect.ratio = 1.6)
ggsave(file.path(OUT, "panel_e.pdf"), p_e, width = 3.0, height = 3.0, device = cairo_pdf)
cat("panel e done\n")

## =====================================================================
## PANEL f — full-width row of 6 horizontal-lollipop facets:
## each facet = one top-variable program's z-score across the 14 regions
## (regions on shared y-axis; region tick labels only on leftmost facet).
## Replaces the prior coord_polar petals (illegible at final size).
## =====================================================================
top6 <- var_df %>% arrange(desc(eta2_region)) %>% slice_head(n = 6) %>% pull(program)
# facet title = full functional label; split "P{n}" onto its own line so the
# long name fits within a narrow facet strip without truncation/overlap.
f_lab6 <- setNames(sub("^(P[0-9]+) ", "\\1\n", .lab(top6)), top6)
# regions ordered top->bottom by anatomical REGION_ORDER (V1 at top)
f_reg_lvls <- rev(REGION_ORDER)
f_df <- as.data.frame(M_z[REGION_ORDER, top6]) %>%
  mutate(region = factor(REGION_ORDER, levels = f_reg_lvls)) %>%
  pivot_longer(-region, names_to = "program", values_to = "z") %>%
  mutate(program = factor(f_lab6[program], levels = f_lab6[top6]),
         lobe = factor(LOBE[as.character(region)], levels = LOBE_ORDER))
p_f <- ggplot(f_df, aes(x = z, y = region, fill = lobe)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_segment(aes(x = 0, xend = z, yend = region, colour = lobe), linewidth = 0.45) +
  geom_point(shape = 21, colour = "white", stroke = 0.25, size = 1.5) +
  facet_grid(~ program) +   # cols only -> single shared y-axis on left facet
  scale_fill_manual(values = lobe_pal, name = "lobe") +
  scale_colour_manual(values = lobe_pal, guide = "none") +
  scale_x_continuous(name = "region z-score",
                     expand = expansion(mult = c(0.06, 0.08))) +
  labs(y = NULL, title = "Regional activity of top-variable programs") +
  theme_nat() +
  theme(axis.text.y = element_text(size = 5.0, colour = "black"),
        axis.text.x = element_text(size = 5.0),
        strip.text = element_text(size = 5.4, face = "bold", lineheight = 0.9),
        panel.spacing.x = unit(2.2, "mm"),
        legend.position = "bottom", legend.key.size = unit(2.6, "mm"))
ggsave(file.path(OUT, "panel_f.pdf"), p_f, width = 7.0, height = 2.4, device = cairo_pdf)
cat("panel f done\n")

## =====================================================================
## PANEL g — spotlight violins (programs 14,6,1,19) from parquet subsample
## =====================================================================
# facet labels = full functional labels for the 4 spotlight programs
g_facet <- c(p14 = .lab("14"), p6 = .lab("6"), p1 = .lab("1"), p19 = .lab("19"))
g_long <- gsub_df %>%
  pivot_longer(c(p1, p6, p14, p19), names_to = "program", values_to = "act") %>%
  mutate(program = recode(program, !!!g_facet),
         program = factor(program, levels = unname(g_facet)),
         region = factor(region, levels = REGION_ORDER),
         lobe = factor(LOBE[as.character(region)], levels = LOBE_ORDER))
g_med <- g_long %>% group_by(program, region, lobe) %>%
  summarise(med = median(act), .groups = "drop")
# cap y for visibility (99th pct per program)
g_cap <- g_long %>% group_by(program) %>%
  summarise(ymax = quantile(act, 0.99), .groups = "drop")
g_long <- g_long %>% left_join(g_cap, by = "program")

p_g <- ggplot(g_long, aes(region, act, fill = lobe)) +
  geom_violin(scale = "width", width = 0.9, linewidth = 0.15, colour = "grey30") +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.2,
               fill = "white", colour = "grey20") +
  geom_point(data = g_med, aes(region, med), inherit.aes = FALSE,
             size = 0.4, colour = "black") +
  facet_wrap(~ program, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = lobe_pal, name = "lobe") +
  labs(x = NULL, y = "per-cell program activity",
       title = "Spotlight programs across regions",
       subtitle = "P14 ↔ L3-L4 IT RORB (within-subclass η²=0.76)") +
  theme_nat() +
  # lobe legend dropped here: panel f (immediately above g in the assembled
  # layout) carries the identical shared lobe legend, so g need not repeat it.
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.0),
        strip.text = element_text(size = 5.2, face = "bold"),
        legend.position = "none")
ggsave(file.path(OUT, "panel_g.pdf"), p_g, width = 4.2, height = 3.2, device = cairo_pdf)
cat("panel g done\n")

## =====================================================================
## PANEL h — effect-size bubble matrix (14 variable programs x region)
## =====================================================================
h_progs <- var_df %>% filter(class == "variable") %>%
  arrange(desc(eta2_region)) %>% pull(program)
h_df <- as.data.frame(M_z[REGION_ORDER, h_progs]) %>%
  mutate(region = factor(REGION_ORDER, levels = REGION_ORDER)) %>%
  pivot_longer(-region, names_to = "program", values_to = "z") %>%
  mutate(program = factor(program, levels = rev(h_progs)))
# y-axis = full functional label (14 rows on the short axis: fits cleanly)
h_ylab <- setNames(.lab(rev(h_progs)), rev(h_progs))
p_h <- ggplot(h_df, aes(region, program)) +
  geom_point(aes(size = abs(z), fill = z), shape = 21, colour = "grey40", stroke = 0.2) +
  scale_fill_gradientn(colours = div_pal, limits = c(-2.6, 2.6),
                       oob = scales::squish, name = "z-score") +
  scale_size_continuous(range = c(0.3, 3.4), name = "|z|") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(labels = h_ylab) +
  labs(x = NULL, y = "variable program (η² desc)",
       title = "Region-specific effect sizes") +
  theme_nat() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 5.5),
        axis.text.y = element_text(size = 5.0, colour = "black"),
        legend.position = "right",
        panel.grid.major = element_line(linewidth = 0.15, colour = "grey92"),
        # dotplot was slightly too wide for its cell -> make it a touch taller.
        aspect.ratio = 0.7)
ggsave(file.path(OUT, "panel_h.pdf"), p_h, width = 3.5, height = 3.2, device = cairo_pdf)
cat("panel h done\n")

## =====================================================================
## PANEL i — rank-shift alluvial across lobes (top ~10 variable programs)
## =====================================================================
# Draw the top 10 from the 54 retained programs.
# (drop the 6 cohort-technical EXCLUDED programs: cnmf 9, 18, 19, 35, 52, 57).
EXCLUDED_CNMF <- as.character(c(9, 18, 19, 35, 52, 57))
top10 <- var_df %>% filter(!program %in% EXCLUDED_CNMF) %>%
  arrange(desc(eta2_region)) %>% slice_head(n = 10) %>% pull(program)
# mean activity per lobe (mean over regions in lobe), then rank programs within lobe
lobe_of_region <- LOBE[REGION_ORDER]
i_long <- as.data.frame(M_mean[REGION_ORDER, top10]) %>%
  mutate(region = REGION_ORDER, lobe = lobe_of_region[REGION_ORDER]) %>%
  pivot_longer(c(-region, -lobe), names_to = "program", values_to = "act") %>%
  group_by(lobe, program) %>% summarise(act = mean(act), .groups = "drop") %>%
  group_by(lobe) %>% mutate(rank = rank(-act, ties.method = "first")) %>% ungroup() %>%
  mutate(lobe = factor(lobe, levels = LOBE_ORDER),
         program = factor(program, levels = top10))
i_pal <- setNames(colorRampPalette(brewer.pal(10, "Spectral"))(length(top10)), top10)
# rank-slot alluvial: stratum = rank position (1..10), alluvium = program.
# Each program flows between its rank slot in each lobe -> crossing ribbons.
nP <- length(top10)
# stratum factor: level order controls vertical stacking. ggalluvial stacks the
# FIRST level at the BOTTOM, so put rank nP first (bottom) and rank 1 last (top).
i_long <- i_long %>% mutate(rankf = factor(rank, levels = nP:1))  # rank1 -> top
# slot centers: with nP unit-height strata stacked, the j-th level from the
# bottom is centered at y = j - 0.5.  level nP:1 means bottom->top = rank nP..1,
# so rank r sits at y = (nP - r) + 0.5.
rank_breaks <- (nP - (1:nP)) + 0.5            # y center of each rank slot
rank_labels <- paste0("rank ", 1:nP)           # rank 1 (top) .. rank nP (bottom)
# Convert source component IDs to retained display IDs for strata text;
# remap to new_P display ("P13", "P9"...) using prog_pn lookup. The `program` factor levels stay as
# cnmf strings (so alluvium fill/data join is intact); add a new `lab_pn` column for the strata text.
i_long <- i_long %>% mutate(lab_pn = unname(prog_pn[as.character(program)]))
p_i <- ggplot(i_long, aes(x = lobe, y = 1, stratum = rankf, alluvium = program,
                          fill = program, label = lab_pn)) +
  geom_alluvium(width = 0.26, alpha = 0.62, colour = "white", linewidth = 0.12,
                curve_type = "sigmoid") +
  geom_stratum(width = 0.26, fill = NA, colour = "grey80", linewidth = 0.15) +
  geom_text(stat = "stratum", size = 1.75, colour = "grey10", fontface = "bold") +
  # strata text = new_P tag ("P{new_int}"); legend = full functional labels
  scale_fill_manual(values = i_pal, name = "program", labels = .lab(top10)) +
  scale_y_continuous(name = "activity rank within lobe (1 = highest)",
                     breaks = rank_breaks, labels = rank_labels,
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, title = "Program rank shifts across lobes") +
  theme_nat() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5),
        axis.text.y = element_text(size = 5.0, colour = "grey30"),
        axis.ticks.y = element_line(linewidth = 0.25),
        legend.position = "right", legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 4.4))
ggsave(file.path(OUT, "panel_i.pdf"), p_i, width = 4.4, height = 3.2, device = cairo_pdf)
cat("panel i done\n")

## =====================================================================
## ASSEMBLE
## =====================================================================
design <- "
AAAAAA
BCCDDE
FFFFFF
GGGHHH
GGGIII
"
fig <- p_a + p_b + p_c + p_d + p_e + p_f + p_g + p_h + p_i +
  plot_layout(design = design, heights = c(1.85, 1.35, 0.78, 1.0, 1.0)) +
  plot_annotation(tag_levels = "a",
                  title = "Fig. 3  Cross-region variation of cortical gene programs",
                  theme = theme(plot.title = element_text(size = 10, face = "bold",
                                                          family = "Helvetica"))) &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(face = "bold", size = 9, family = "Helvetica",
                                hjust = 0, vjust = 1))

ggsave(file.path(OUT, "fig3_crossregion.pdf"), fig,
       width = 13.0, height = 13.0, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT, "fig3_crossregion.png"), fig,
       width = 13.0, height = 13.0, dpi = 200, bg = "white", limitsize = FALSE)
file.copy(file.path(OUT, "fig3_crossregion.png"), "/tmp/fig3_aspect.png", overwrite = TRUE)
# task-requested alias name for the reassembled full figure
file.copy(file.path(OUT, "fig3_crossregion.pdf"),
          file.path(OUT, "fig3_full.pdf"), overwrite = TRUE)
cat("assembly done\n")

## ---- emit gradient numbers for report (clustered axis) -------------
rep_grad <- var_grad %>% slice_head(n = nshow) %>%
  left_join(var_df[, c("program","eta2_region")], by = "program") %>%
  select(program, eta2_region, axis_slope, axis_r, axis_p)
cat("\n=== CLUSTERED REGION ORDER ===\n")
cat(paste(seq_along(leaf_order), leaf_order, sep = ":"), sep = "  "); cat("\n")
cat("=== GRADIENT (shown variable programs, clustered axis) ===\n")
print(rep_grad, row.names = FALSE, digits = 3)
write.table(rep_grad, file.path(OUT, "panel_d_gradient_numbers.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nALL PANELS COMPLETE\n")

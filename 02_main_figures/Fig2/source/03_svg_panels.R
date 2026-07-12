#!/usr/bin/env Rscript
# =====================================================================
# Fig.3  Cross-region gene-program variation  (human cortex, cNMF K=60)
# 9 panels (a-i). Submission-grade, Nature-family style.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ComplexHeatmap)
  library(circlize); library(ggalluvial); library(ggrepel)
  library(dplyr); library(tidyr); library(grid)
  library(scales); library(RColorBrewer); library(svglite)
})

set.seed(42)

## ---- paths ----------------------------------------------------------
RES <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
OUT <- "CORTEX_PROGRAM_ROOT/figures/fig3"
SVGD <- file.path(OUT, "svg_panels")            # individual SVG panels for svgutils compose
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(SVGD, recursive = TRUE, showWarnings = FALSE)
mm2in <- function(x) as.numeric(x) / 25.4
svgf  <- function(id) file.path(SVGD, sprintf("fig3_%s.svg", id))

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
# cohort-robustness tier (WITHIN the 14 ANOVA-variable programs): which survive
# cohort + donor ANCOVA control vs which drop out. Robust-7 get a heavier,
# saturated red (the figure's highlight); cohort-sensitive get a muted warm tone
# (de-emphasised); the stable bulk stays grey. Source (real files):
#   results/crossregion_v1/program_variability.tsv and the encoded robust/sensitive sets
ROBUST7    <- as.character(c(1, 3, 4, 6, 8, 10, 14))
SENSITIVE7 <- as.character(c(9, 18, 19, 35, 37, 52, 57))
tier_of <- function(p) factor(
  ifelse(p %in% ROBUST7,    "cohort-robust",
  ifelse(p %in% SENSITIVE7, "cohort-sensitive", "stable")),
  levels = c("cohort-robust", "cohort-sensitive", "stable"))
tier_pal <- c(`cohort-robust`    = "#B2182B",   # strong dark red  (highlight robust-7)
              `cohort-sensitive` = "#F4A582",   # muted warm       (de-emphasised)
              stable             = "#9AA0A6")    # grey bulk
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

# Fig.54 (甲案): physically exclude the six cohort-technical programs from EVERY
# panel of this figure. The whole figure is restricted to the 54 biologically-
# interpreted programs; the six (P9/P18/P19/P35/P52/P57) drop out of the heatmap
# columns, the variability ranking, the scatter, the gradient traces and the
# region-variable counts. After exclusion the region-variable set is the seven
# cohort-robust programs plus P37 = 8 (P37 is retained-biological but cohort-
# sensitive); the six removed programs are exactly the cohort-technical ones.
EXCLUDE6 <- as.character(c(9, 18, 19, 35, 52, 57))
progs <- setdiff(as.character(1:60), EXCLUDE6)    # 54 biologically-interpreted
# matrices: regions (rows) x programs (cols). transpose to programs x regions
M_mean <- as.matrix(mean_df[, progs])            # region x program
M_z    <- as.matrix(z_df[,    progs])            # region x program
stopifnot(all(rownames(M_z) %in% REGION_ORDER))

var_df$program <- as.character(var_df$program)
var_df <- var_df[!(var_df$program %in% EXCLUDE6), ]   # 54-program multiverse
var_df$class   <- factor(var_df$class, levels = c("variable","stable"))
variable_progs <- var_df$program[var_df$class == "variable"]  # = ROBUST7 + P37 (8)

# region-level mean activity & CV per program (across 14 regions)
prog_mean_act <- colMeans(M_mean)                       # mean activity across regions
prog_cv       <- apply(M_mean, 2, function(x) sd(x)/mean(x))

## ---- program -> functional label lookup -----------------------------
## display label = "P{n} {name_short}"; strip redundant trailing " P{n}"
## suffix (disambiguation artifact); brain-weak programs get trailing "*".
nm_df <- read.table(file.path(RES, "program_names.tsv"),
                    sep = "\t", header = TRUE, quote = "", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE)
# program_names.tsv uses cnmf_component as the original component key.
# (old_P integer 1-60), and new_P column holds "P1"-"P54" or "EXCLUDED".
# We keep old_P as the lookup key (matching data matrix column names), but
# generate labels that display the new_P number so every figure panel shows
# the post-exclusion P-numbers (P1-P54).
nm_df$program   <- as.character(nm_df$cnmf_component)   # old_P key (1-60)
nm_df$new_P_num <- suppressWarnings(as.integer(sub("^P", "", nm_df$new_P)))  # 1-54 or NA for EXCLUDED
.mk_label <- function(p_old, new_p_num, ns, conf) {
  if (is.na(new_p_num)) return(paste0("EXCLUDED_P", p_old))
  ns  <- trimws(ns)
  suf <- paste0(" P", new_p_num)
  if (grepl(paste0(suf, "$"), ns)) ns <- trimws(sub(paste0(suf, "$"), "", ns))
  star <- ifelse(conf == "brain-weak", "*", "")
  paste0("P", new_p_num, " ", ns, star)
}
prog_lab <- setNames(
  mapply(.mk_label, nm_df$program, nm_df$new_P_num, nm_df$name_short, nm_df$confidence),
  nm_df$program)
# bare "P{n}" fallback using new_P number (used for non-highlighted items in dense panels)
prog_pn  <- setNames(
  ifelse(is.na(nm_df$new_P_num), paste0("P", nm_df$program), paste0("P", nm_df$new_P_num)),
  nm_df$program)
# safety: any program id missing from name map falls back to "P{n}" (old_P)
.lab <- function(p) ifelse(p %in% names(prog_lab), prog_lab[p], paste0("P", p))
# compact line-end label: "P{n} {short}" with the short name capped at a word
# boundary (no mid-word break) + ellipsis if truncated. Used ONLY for panel-d
# line-end labels to keep the right-margin block narrow (smaller panel-d
# footprint). Does NOT change underlying names or any other panel's labels.
.lab_short <- function(p, maxchar = 17L) {
  full <- .lab(p)                                   # "P{n} short name*"
  m    <- regmatches(full, regexec("^(P[0-9]+)\\s+(.*)$", full))[[1]]
  if (length(m) < 3) return(full)
  pn <- m[2]; nm <- m[3]
  star <- ifelse(grepl("\\*$", nm), "*", ""); nm <- sub("\\*$", "", nm)
  if (nchar(nm) > maxchar) {
    cut <- substr(nm, 1, maxchar)
    sp  <- regexpr("\\s[^\\s]*$", cut)              # last space within the cap
    if (sp > 4) cut <- substr(cut, 1, sp - 1)       # trim back to a word boundary
    nm <- paste0(sub("[\\s[:punct:]]+$", "", cut), "…")
  }
  paste0(pn, " ", nm, star)
}

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
prog_tier  <- setNames(as.character(tier_of(progs)), progs)
prog_eta   <- setNames(var_df$eta2_region[match(progs, var_df$program)], progs)
# Keep the panel biologically readable, but trim long labels to the shortest
# unambiguous "P{n} short-name" form so the heatmap can carry clearer group
# structure without a label wall.
col_lab_a <- vapply(progs, .lab_short, character(1), maxchar = 14L)
prog_tier_fac <- factor(prog_tier[progs],
                        levels = c("cohort-robust", "cohort-sensitive", "stable"))

top_anno <- HeatmapAnnotation(
  tier = prog_tier,
  eta2  = anno_barplot(prog_eta, gp = gpar(fill = "#6E7B8B", col = NA),
                       height = unit(7, "mm"), axis_param = list(gp = gpar(fontsize = 5))),
  col = list(tier = tier_pal),
  annotation_name_gp = gpar(fontsize = 6),
  annotation_legend_param = list(tier = list(title = "cohort robustness",
                                  title_gp = gpar(fontsize = 6),
                                  labels_gp = gpar(fontsize = 5.5))),
  simple_anno_size = unit(2.5, "mm"),
  show_legend = c(tier = TRUE)
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
CELL_A <- 3.4                          # per-cell size in mm (1:1 square cells); widened 2.6->3.4 for fully-expanded col labels (Fig 2A redraw 2026-06-23)
ht_a <- Heatmap(
  Ha, name = "z-score",
  col = colorRamp2(seq(-2.5, 2.5, length.out = 256), div_pal),
  top_annotation = top_anno, left_annotation = left_anno,
  cluster_rows = TRUE, cluster_columns = TRUE,
  column_split = prog_tier_fac,
  column_gap = unit(c(1.2, 2.0), "mm"),
  show_row_dend = TRUE, show_column_dend = TRUE,
  row_dend_width = unit(6, "mm"), column_dend_height = unit(6, "mm"),
  row_names_gp = gpar(fontsize = 6),
  column_labels = col_lab_a, column_names_gp = gpar(fontsize = 4.8),
  column_names_rot = 90,
  column_names_max_height = unit(38, "mm"),
  column_title = "programs (n=54, cNMF K=60 backbone)",
  column_title_gp = gpar(fontsize = 6.5, fontface = "bold"),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                              labels_gp = gpar(fontsize = 5.5),
                              legend_height = unit(20, "mm")),
  row_title = "cortical region",
  row_title_gp = gpar(fontsize = 6.5),
  # SQUARE cells (user FB): each grid cell is exactly s x s mm, so the body is
  # ncol*s wide x nrow*s tall. 54 programs x 14 regions -> body AR = 54/14 ~ 3.9
  # (a WIDE SHORT banner). width/height pin the BODY; tracks/labels stay fixed mm.
  width  = ncol(Ha) * unit(CELL_A, "mm"),
  height = nrow(Ha) * unit(CELL_A, "mm")
)
# SQUARE-cell banner: body = 54*2.6 = 140.4mm wide x 14*2.6 = 36.4mm tall.
# Canvas must be WIDE-SHORT and large enough to also hold the fixed-mm tracks
# (col dend 6mm, class strip 2.5mm, eta2 barplot 7mm, lobe strip 2.5mm, row
# dend 6mm), the rotated column labels, and the right-side legend. 9.6x3.1in
# comfortably contains it; ink-crop trims the slack -> resulting AR ~3.5-4.5.
PANEL_A_W <- 13.5; PANEL_A_H <- 6.0  # widened 9.6->13.5 + raised 4.2->6.0 for fully-expanded 54-prog col labels (Fig 2A redraw 2026-06-23, 羅老師 Direction A)
pdf(file.path(OUT, "panel_a.pdf"), width = PANEL_A_W, height = PANEL_A_H, family = "Helvetica")
draw(ht_a, heatmap_legend_side = "right", annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()
# SVG twin (for svgutils tangram compose)
svglite(svgf("a"), width = PANEL_A_W, height = PANEL_A_H, bg = "white")
draw(ht_a, heatmap_legend_side = "right", annotation_legend_side = "right",
     merge_legend = TRUE)
invisible(dev.off())
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
         tier = tier_of(program),               # cohort robustness (within the 14)
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
p_b <- ggplot(b_df, aes(x = eta2_region, y = ypos, colour = tier)) +
  geom_segment(aes(x = 0, xend = eta2_region, yend = ypos), linewidth = 0.4) +
  geom_point(aes(size = tier == "cohort-robust")) +
  scale_size_manual(values = c(`TRUE` = 1.9, `FALSE` = 1.1), guide = "none") +
  # leader lines: from each labelled point to its fixed label anchor
  geom_segment(data = b_lab, inherit.aes = FALSE,
               aes(x = eta2_region, y = ypos, xend = lab_x, yend = lab_y),
               colour = "grey55", linewidth = 0.18) +
  geom_text(data = b_lab, inherit.aes = FALSE,
            aes(x = lab_x, y = lab_y, label = lab),
            hjust = 0, size = 1.8, colour = "black") +
  scale_colour_manual(values = tier_pal, name = "cohort robustness") +
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
svglite(svgf("b"), width = 3.0, height = 3.6, bg = "white"); print(p_b); invisible(dev.off())
cat("panel b done\n")

## =====================================================================
## PANEL c — mean activity vs CV scatter + marginal densities
## =====================================================================
c_df <- data.frame(program = progs,
                   mean_act = prog_mean_act[progs],
                   cv = prog_cv[progs],
                   class = var_df$class[match(progs, var_df$program)],
                   tier = tier_of(progs),
                   eta2 = prog_eta[progs])
top_lab <- var_df %>% arrange(desc(eta2_region)) %>% slice_head(n = 5) %>% pull(program)
# label only the top-5 most variable programs; stronger repulsion keeps the
# few callouts from overlapping each other or the dense point cloud.
c_df$lab <- ifelse(c_df$program %in% top_lab, .lab(c_df$program), "")

p_c_main <- ggplot(c_df, aes(mean_act, cv, colour = tier)) +
  geom_point(aes(size = eta2), alpha = 0.85) +
  geom_text_repel(aes(label = lab), size = 2.0, colour = "black",
                  max.overlaps = Inf, segment.size = 0.2, min.segment.length = 0,
                  box.padding = 0.8, point.padding = 0.4, force = 8, force_pull = 0.5,
                  seed = 42) +
  scale_colour_manual(values = tier_pal, name = "cohort robustness", guide = "none") +
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
svglite(svgf("c"), width = 3.0, height = 3.0, bg = "white"); print(p_c_grp); invisible(dev.off())
# wrap as single element so patchwork tags it once in the assembled figure
p_c <- wrap_elements(full = p_c_grp)
cat("panel c done\n")

## =====================================================================
## PANEL d — REPLACEMENT: gradient along PROGRAM-CLUSTERED region axis
## (region order = ward.D2 / Euclidean clustering of regions on their
##  60-program z-vector; replaces invalid PCA-PC1 cohort-confounded axis)
## =====================================================================
# cluster the 14 regions on their 54-program z-vector
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
  # line-end label: short functional name + gradient r on a second line. Cap the
  # name length (no mid-word break) so the right-margin label block stays narrow
  # -> shrinks panel-d footprint without changing any data/name/highlight.
  mutate(txt = sprintf("%s\n(r=%.2f)", .lab_short(program), axis_r))
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
                  nudge_x = 0.3, size = 1.78, segment.size = 0.18, lineheight = 0.85,
                  xlim = c(14.25, NA), max.overlaps = Inf, force = 3,
                  box.padding = 0.32, point.padding = 0.2, seed = 42) +
  scale_colour_manual(values = prog_cols8, guide = "none") +
  scale_fill_manual(values = lobe_pal, name = "lobe",
                    guide = guide_legend(override.aes = list(colour = NA),
                                         nrow = 1)) +
  # shrink right expansion: short capped labels need far less right margin than
  # the old full names did (0.62 -> 0.30) -> panel-d gets much less wide.
  scale_x_discrete(expand = expansion(mult = c(0.04, 0.30))) +
  scale_y_continuous(expand = expansion(mult = c(0.16, 0.06))) +
  labs(x = "region (program-clustered order)", y = "region z-score",
       title = "Program gradients along clustered region axis",
       subtitle = "axis = ward.D2 clustering of regions on 54-program profile") +
  theme_nat() +
  # lobe legend moved to a single-row strip along the BOTTOM (was right margin):
  # frees the right side so the panel reads near-square (lower natural aspect),
  # while staying >=5pt and not colliding with the traces.
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 5.6),
        legend.text = element_text(size = 5.2),
        legend.key.size = unit(2.4, "mm"),
        legend.margin = margin(t = 0, b = 0))
# squarer canvas (was 4.0x3.2, AR 1.25): narrower width + slightly taller (bottom
# legend) drives the cropped natural aspect toward ~1.0-1.15.
ggsave(file.path(OUT, "panel_d.pdf"), p_d, width = 3.3, height = 3.1, device = cairo_pdf)
svglite(svgf("d"), width = 3.3, height = 3.1, bg = "white"); print(p_d); invisible(dev.off())
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
                  max.overlaps = Inf, segment.size = 0.2, min.segment.length = 0,
                  box.padding = 0.5, point.padding = 0.3, force = 6, seed = 42,
                  show.legend = FALSE) +
  scale_fill_manual(values = lobe_pal, name = "lobe") +
  scale_colour_manual(values = lobe_pal, guide = "none") +
  labs(x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]),
       title = "Region similarity (program profile)") +
  coord_cartesian(clip = "off") +
  # NO forced aspect.ratio: let the biplot panel fill the canvas so the
  # plot is naturally near-square (was aspect.ratio=1.6 -> tall-narrow strip
  # with ~70pt dead margin L+R, which read as the d-e gap). Compact inset
  # legend bottom-right, tight margins.
  theme_nat() + theme(legend.position = c(0.90, 0.18),
                      legend.background = element_rect(fill = "white", colour = NA),
                      legend.key.size = unit(2.4, "mm"),
                      legend.title = element_text(size = 5.6),
                      legend.text  = element_text(size = 5.2),
                      plot.margin = margin(t = 1, r = 3, b = 2, l = 2, unit = "pt"))
ggsave(file.path(OUT, "panel_e.pdf"), p_e, width = 3.4, height = 3.0, device = cairo_pdf)
svglite(svgf("e"), width = 3.4, height = 3.0, bg = "white"); print(p_e); invisible(dev.off())
cat("panel e done\n")

## =====================================================================
## PANEL f — full-width row of 6 horizontal-lollipop facets:
## each facet = one top-variable program's z-score across the 14 regions
## (regions on shared y-axis; region tick labels only on leftmost facet).
## Replaces the prior coord_polar petals (illegible at final size).
## =====================================================================
# Restrict panel-f examples to the region-variable programs.
# to COHORT-ROBUST programs only (ROBUST7 = {old_P 1,3,4,6,8,10,14} = new_P {1,3,4,6,8,9,13});
# the prior unfiltered top-6-by-eta2 leaked cohort-technical P19/P18/P57 (excluded
# from biological interpretation). Pick the 6 most region-variable among
# the cohort-robust set -> new_P P6,P1,P13,P4,P8,P9. Names via .lab()=name_short.
top6 <- var_df %>% filter(as.character(program) %in% ROBUST7) %>%
  arrange(desc(eta2_region)) %>% slice_head(n = 6) %>%
  pull(program) %>% as.character()
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
svglite(svgf("f"), width = 7.0, height = 2.4, bg = "white"); print(p_f); invisible(dev.off())
cat("panel f done\n")

## =====================================================================
## PANEL g — spotlight violins (programs 14,6,1) from parquet subsample
## =====================================================================
# Fig.54 (甲案): the former fourth spotlight P19 is a cohort-technical program
# (excluded from biological interpretation), so it is dropped here; the spotlight
# now shows the three cohort-robust star programs new_P13/P6/P1 (old_P14/6/1). (P19's per-cell
# subsample column p19 in panel_g_subsample.tsv is simply not pivoted.)
g_ids   <- c("14", "6", "1")
g_facet <- setNames(
  sprintf("%s · cohort-%s", .lab(g_ids),
          ifelse(g_ids %in% ROBUST7, "robust", "sensitive")),
  paste0("p", g_ids))
g_long <- gsub_df %>%
  pivot_longer(c(p1, p6, p14), names_to = "program", values_to = "act") %>%
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
  # boxplots (was violins: ugly with the heavy zero mass). lobe fill kept;
  # thin lines + tiny faint outliers so the 40-69% zero data reads cleanly.
  geom_boxplot(width = 0.7, linewidth = 0.3, colour = "grey25",
               outlier.size = 0.3, outlier.colour = "grey55",
               outlier.alpha = 0.4, outlier.stroke = 0,
               fatten = 1.4) +
  geom_point(data = g_med, aes(region, med), inherit.aes = FALSE,
             size = 0.4, colour = "black") +
  facet_wrap(~ program, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = lobe_pal, name = "lobe") +
  labs(x = NULL, y = "per-cell program activity",
       title = "Spotlight programs across regions",
       subtitle = "P13 ↔ L3-L4 IT RORB (per-cell region η²=0.220); facet tag = cohort robustness") +
  theme_nat() +
  # lobe legend dropped here: panel f (immediately above g in the assembled
  # layout) carries the identical shared lobe legend, so g need not repeat it.
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.0),
        strip.text = element_text(size = 5.2, face = "bold"),
        legend.position = "none")
ggsave(file.path(OUT, "panel_g.pdf"), p_g, width = 4.2, height = 3.2, device = cairo_pdf)
svglite(svgf("g"), width = 4.2, height = 3.2, bg = "white"); print(p_g); invisible(dev.off())
cat("panel g done\n")

## =====================================================================
## PANEL h — REGION-SIGNATURE heatmap (region-anchored; replaces the prior
## program-centric effect-size dotplot, which overlapped panel f in intent
## and did not let a reader read off "region X = high A / low B").
## rows = 14 cortical regions (anatomical order V1..ACC); columns = the UNION
## of each region's top-K distinctive programs (K up + K down per region),
## ordered so columns sweep V1-programs (left) -> ACC-programs (right) =>
## the signature forms a readable diagonal band. z-coloured, SQUARE cells
## (skill rule 1: ncol*s x nrow*s -> wide-short banner, no stretch). Each
## region's #1-up program is marked with a black ring so the eye anchors the
## row's defining program. A reader scans a region row = its named signature
## (reds = high programs, blues = low).
## =====================================================================
SIG_K <- 2L                                     # top-2 up + top-2 down per region
# per-region top-K up / down on the region x program z matrix
.sig_topk <- function(reg, k, up = TRUE) {
  v <- M_z[reg, progs]
  names(sort(v, decreasing = up))[seq_len(k)]
}
# union, preserving first-seen order across regions (V1 first)
sig_union <- character(0)
sig_top1  <- setNames(character(length(REGION_ORDER)), REGION_ORDER)  # each region #1-up
for (reg in REGION_ORDER) {
  ups <- .sig_topk(reg, SIG_K, TRUE); dns <- .sig_topk(reg, SIG_K, FALSE)
  sig_top1[reg] <- ups[1]
  for (p in c(ups, dns)) if (!(p %in% sig_union)) sig_union <- c(sig_union, p)
}
# column order: by the region (anatomical rank) where each program peaks (argmax z),
# tie-break by that peak z desc -> columns ordered V1-ish .. ACC-ish (diagonal band)
.colkey <- function(p) {
  br <- names(which.max(M_z[REGION_ORDER, p])); c(anat_rank[[br]], -M_z[br, p])
}
sig_cols <- sig_union[order(sapply(sig_union, function(p) .colkey(p)[1]),
                            sapply(sig_union, function(p) .colkey(p)[2]))]
# matrix: region (rows, V1..ACC) x program (cols, ordered)
Hh <- M_z[REGION_ORDER, sig_cols, drop = FALSE]
colnames(Hh) <- sig_cols
h_collab <- .lab(sig_cols)                       # "P{n} functional name"
# mark each region's #1-up program: which (row,col) cells get a ring
sig_mark <- matrix(FALSE, nrow = nrow(Hh), ncol = ncol(Hh),
                   dimnames = dimnames(Hh))
for (reg in REGION_ORDER) sig_mark[reg, sig_top1[reg]] <- TRUE
CELL_H <- 2.6                                    # square cell size (mm), matches panel a
left_anno_h <- rowAnnotation(
  lobe = LOBE[REGION_ORDER],
  col = list(lobe = lobe_pal),
  annotation_name_gp = gpar(fontsize = 6),
  show_legend = FALSE,                           # lobe legend already on panels a/f
  simple_anno_size = unit(2.5, "mm"))
ht_h <- Heatmap(
  Hh, name = "z-score",
  col = colorRamp2(seq(-2.5, 2.5, length.out = 256), div_pal),
  left_annotation = left_anno_h,
  cluster_rows = FALSE, cluster_columns = FALSE,   # fixed anatomical / diagonal order
  row_names_side = "left", row_names_gp = gpar(fontsize = 5.6),
  column_names_side = "top", column_names_rot = 90,
  column_labels = h_collab, column_names_gp = gpar(fontsize = 5.0),
  rect_gp = gpar(col = "white", lwd = 0.4),
  # ring each region's #1-up program cell (anchors the row's defining program)
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (sig_mark[i, j])
      grid.rect(x, y, w, h, gp = gpar(col = "black", lwd = 0.9, fill = NA))
  },
  row_title = "cortical region", row_title_gp = gpar(fontsize = 6.5),
  column_title = "region-distinctive program (top-2 up/down per region; ring = region's #1 up)",
  column_title_side = "bottom", column_title_gp = gpar(fontsize = 5.6),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                              labels_gp = gpar(fontsize = 5.5),
                              legend_height = unit(15, "mm")),
  width  = ncol(Hh) * unit(CELL_H, "mm"),
  height = nrow(Hh) * unit(CELL_H, "mm"))
# body = 40*2.6 = 104mm wide x 14*2.6 = 36.4mm tall (wide-short banner). Canvas
# wide-short + tall enough for the rotated 90deg column labels (~34 char) on top,
# the bottom caption, left lobe strip + region names, and the right legend.
PANEL_H_W <- 6.6; PANEL_H_H <- 3.6
pdf(file.path(OUT, "panel_h.pdf"), width = PANEL_H_W, height = PANEL_H_H, family = "Helvetica")
draw(ht_h, heatmap_legend_side = "right", padding = unit(c(0, 0, 0, 0), "mm"))
dev.off()
svglite(svgf("h"), width = PANEL_H_W, height = PANEL_H_H, bg = "white")
draw(ht_h, heatmap_legend_side = "right", padding = unit(c(0, 0, 0, 0), "mm"))
invisible(dev.off())
# grob wrapper so the local patchwork preview (ASSEMBLE block) still builds
# (panel h is now a ComplexHeatmap grob, not a ggplot object).
gb_h <- grid.grabExpr(draw(ht_h, heatmap_legend_side = "right",
                           padding = unit(c(0, 0, 0, 0), "mm")),
                      width = PANEL_H_W, height = PANEL_H_H)
p_h <- patchwork::wrap_elements(full = gb_h)
cat(sprintf("panel h done (region-signature heatmap: %d regions x %d programs)\n",
            nrow(Hh), ncol(Hh)))

## ---- SUPP TABLE: 14 regions x top-5 up + top-5 down programs ---------
sig_supp <- lapply(REGION_ORDER, function(reg) {
  ups <- .sig_topk(reg, 5L, TRUE); dns <- .sig_topk(reg, 5L, FALSE)
  rbind(
    data.frame(region = reg, lobe = unname(LOBE[reg]),
               direction = "up",   rank = 1:5, program = ups,
               program_name = .lab(ups), z = round(M_z[reg, ups], 3)),
    data.frame(region = reg, lobe = unname(LOBE[reg]),
               direction = "down", rank = 1:5, program = dns,
               program_name = .lab(dns), z = round(M_z[reg, dns], 3)))
}) %>% bind_rows()
write.table(sig_supp, file.path(RES, "supp_table_region_signatures.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("wrote supp_table_region_signatures.tsv (14 regions x top-5 up/down)\n")

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
# Convert source component IDs to retained display IDs for strata text.
i_long <- i_long %>% mutate(lab_pn = unname(prog_pn[as.character(program)]))
p_i <- ggplot(i_long, aes(x = lobe, y = 1, stratum = rankf, alluvium = program,
                          fill = program, label = lab_pn)) +
  geom_alluvium(width = 0.26, alpha = 0.62, colour = "white", linewidth = 0.12,
                curve_type = "sigmoid") +
  geom_stratum(width = 0.26, fill = NA, colour = "grey80", linewidth = 0.15) +
  geom_text(stat = "stratum", size = 2.7, colour = "grey10", fontface = "bold") +
  # strata text = rank slot (kept numeric); legend = full functional labels
  scale_fill_manual(values = i_pal, name = NULL, labels = .lab(top10)) +
  scale_y_continuous(name = "activity rank within lobe (1 = highest)",
                     breaks = rank_breaks, labels = rank_labels,
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, title = "Program rank shifts across lobes") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                             keywidth = unit(3.0, "mm"),
                             keyheight = unit(3.0, "mm"))) +
  theme_nat() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5),
        axis.text.y = element_text(size = 5.0, colour = "grey30"),
        axis.ticks.y = element_line(linewidth = 0.25),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal", legend.key.size = unit(3.0, "mm"),
        legend.text = element_text(size = 5.0),
        legend.margin = margin(t = 1, b = 0))
ggsave(file.path(OUT, "panel_i.pdf"), p_i, width = 4.4, height = 3.7, device = cairo_pdf)
svglite(svgf("i"), width = 4.4, height = 3.7, bg = "white"); print(p_i); invisible(dev.off())
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

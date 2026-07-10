#!/usr/bin/env Rscript
# =====================================================================
# Fig.3 panel d — REPLACEMENT
# Region axis derived by hierarchical clustering of regions on their
# 60-program z-score profile (ward.D2 / Euclidean) -> dendrogram leaf
# order = "program-defined region axis". Gradient recomputed vs this
# clustered axis rank (replaces the invalid PCA-PC1 axis that tracked
# sequencing cohort / sample size).
# Style/palette matched to existing fig3_region.R (theme_nat, Dark2).
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel)
  library(dplyr); library(tidyr); library(grid)
  library(scales); library(RColorBrewer)
})

set.seed(42)

## ---- paths ----------------------------------------------------------
RES <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
OUT <- "CORTEX_PROGRAM_ROOT/figures/fig3"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- lobes (for color annotation only; NOT used to define the axis) -
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

## ---- Nature theme (identical to fig3_region.R) ----------------------
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

progs <- as.character(1:60)
M_mean <- as.matrix(mean_df[, progs])            # region x program (raw mean activity)
M_z    <- as.matrix(z_df[,    progs])            # region x program (z across regions)

var_df$program <- as.character(var_df$program)
var_df$class   <- factor(var_df$class, levels = c("variable","stable"))
variable_progs <- var_df$program[var_df$class == "variable"]

## =====================================================================
## STEP 1 — cluster the 14 REGIONS on their 60-program z-vector
## ward.D2 on Euclidean distance -> dendrogram leaf order = region axis
## =====================================================================
regions <- rownames(M_z)                          # 14 regions
d_reg   <- dist(M_z[regions, progs], method = "euclidean")
hc      <- hclust(d_reg, method = "ward.D2")
leaf_order <- hc$labels[hc$order]                 # dendrogram leaf order

cat("\n=== CLUSTERED REGION ORDER (ward.D2 / Euclidean on 60-program z) ===\n")
cat(paste(seq_along(leaf_order), leaf_order, sep = ":"), sep = "  ")
cat("\nlobe sequence: ",
    paste(LOBE[leaf_order], collapse = " -> "), "\n", sep = "")

## axis rank = position along the clustered leaf order (1..14)
axis_rank <- setNames(seq_along(leaf_order), leaf_order)

reg_order_tab <- data.frame(
  region    = leaf_order,
  leaf_order= seq_along(leaf_order),
  axis_rank = as.integer(axis_rank[leaf_order]),
  lobe      = unname(LOBE[leaf_order])
)
write.table(reg_order_tab[, c("region","leaf_order","axis_rank")],
            file.path(RES, "region_cluster_order.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("wrote region_cluster_order.tsv\n")

## =====================================================================
## STEP 2 — recompute gradient (Pearson r + slope) of each program's
## 14 region-mean values vs the clustered axis rank
## =====================================================================
grad_tab <- bind_rows(lapply(progs, function(p) {
  v  <- M_mean[leaf_order, p]                     # region-mean activity in clustered order
  x  <- axis_rank[leaf_order]
  ct <- suppressWarnings(cor.test(x, v, method = "pearson"))
  sl <- coef(lm(v ~ x))[2]
  data.frame(program   = p,
             axis_slope = unname(sl),
             axis_r     = unname(ct$estimate),
             axis_p     = ct$p.value)
}))

write.table(grad_tab, file.path(RES, "program_gradient_clustered.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("wrote program_gradient_clustered.tsv\n")

## =====================================================================
## STEP 3 — redraw PANEL d along the clustered region order
## show the top variable programs (by |axis_r| on the new axis), among
## the class=="variable" set, with per-program gradient annotation.
## =====================================================================
# rank variable programs by strength of the gradient on the clustered axis
var_grad <- grad_tab %>%
  filter(program %in% variable_progs) %>%
  left_join(var_df[, c("program","eta2_region")], by = "program") %>%
  arrange(desc(abs(axis_r)))
nshow <- min(8L, nrow(var_grad))
show_progs <- var_grad$program[seq_len(nshow)]

d_long <- as.data.frame(M_z[leaf_order, show_progs]) %>%
  mutate(region = factor(leaf_order, levels = leaf_order)) %>%
  pivot_longer(-region, names_to = "program", values_to = "z")

d_grad <- grad_tab %>% filter(program %in% show_progs)
d_lab  <- d_long %>% group_by(program) %>% slice_max(as.integer(region), n = 1) %>%
  left_join(d_grad, by = "program") %>%
  mutate(txt = sprintf("p%s  (r=%.2f)", program, axis_r))

# color programs by sign/strength of gradient (Dark2, matching panel-d family)
ord_progs <- d_grad$program[order(-d_grad$axis_r)]
prog_cols <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(show_progs)),
                      ord_progs)

# lobe color strip under the x-axis to expose anatomy along the clustered order
lobe_strip <- data.frame(region = factor(leaf_order, levels = leaf_order),
                         lobe   = factor(LOBE[leaf_order], levels = LOBE_ORDER),
                         x      = seq_along(leaf_order))
ystrip <- min(c(d_long$z, 0)) - 0.55
ystep  <- 0.32

p_d <- ggplot(d_long, aes(region, z, group = program, colour = program)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
  geom_line(linewidth = 0.55) +
  geom_point(size = 0.7) +
  # lobe color strip beneath the axis (anatomy reference along clustered order)
  geom_tile(data = lobe_strip, inherit.aes = FALSE,
            aes(x = x, y = ystrip, fill = lobe),
            width = 1, height = ystep, colour = "white", linewidth = 0.2) +
  geom_text_repel(data = d_lab, aes(label = txt), hjust = 0, direction = "y",
                  nudge_x = 0.4, size = 1.85, segment.size = 0.2,
                  xlim = c(14.3, NA), max.overlaps = 20) +
  scale_colour_manual(values = prog_cols, guide = "none") +
  scale_fill_manual(values = lobe_pal, name = "lobe",
                    guide = guide_legend(override.aes = list(colour = NA))) +
  scale_x_discrete(expand = expansion(mult = c(0.04, 0.30))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.06))) +
  labs(x = "region (program-clustered order)", y = "region z-score",
       title = "Program gradients along clustered region axis",
       subtitle = "axis = ward.D2 clustering of regions on 60-program profile") +
  theme_nat() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
        legend.position = c(0.18, 0.16),
        legend.key.size = unit(2.4, "mm"),
        legend.background = element_rect(fill = scales::alpha("white", 0.6), colour = NA))

ggsave(file.path(OUT, "panel_d.pdf"), p_d, width = 3.4, height = 3.0, device = cairo_pdf)
# task-requested alias name
file.copy(file.path(OUT, "panel_d.pdf"), file.path(OUT, "fig3_d.pdf"), overwrite = TRUE)
cat("panel d (clustered) done -> panel_d.pdf + fig3_d.pdf\n")

## ---- emit numbers for report ---------------------------------------
rep_grad <- var_grad %>% slice_head(n = nshow) %>%
  mutate(shown = TRUE) %>%
  select(program, eta2_region, axis_slope, axis_r, axis_p)
write.table(rep_grad, file.path(OUT, "panel_d_gradient_numbers.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== SHOWN VARIABLE PROGRAMS (top by |axis_r|, n=", nshow, ") ===\n", sep = "")
print(rep_grad, row.names = FALSE, digits = 3)

cat("\n=== ALL VARIABLE PROGRAMS gradient on clustered axis ===\n")
print(var_grad %>% select(program, eta2_region, axis_slope, axis_r, axis_p),
      row.names = FALSE, digits = 3)

cat("\nPANEL D REPLACEMENT COMPLETE\n")

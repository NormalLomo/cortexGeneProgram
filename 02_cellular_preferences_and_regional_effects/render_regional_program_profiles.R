.annotation_group <- function(value) ifelse(endsWith(as.character(value), "sig"), "class_a", "class_b")
PROJECT_ROOT <- Sys.getenv("CORTEX_PROGRAM_ROOT", unset = "/DATA/cortex_nmf_program")
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(ComplexHeatmap)
    library(circlize)
    library(ggalluvial)
    library(ggrepel)
    library(dplyr)
    library(tidyr)
    library(grid)
    library(scales)
    library(RColorBrewer)
    library(svglite)
})
set.seed(42)
RES <- file.path(PROJECT_ROOT, "results/crossregion_v1")
OUT <- file.path(PROJECT_ROOT, "figures/fig2")
SVGD <- file.path(OUT, "svg_panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(SVGD, recursive = TRUE, showWarnings = FALSE)
mm2in <- function(x) as.numeric(x)/25.4
svgf <- function(id) file.path(SVGD, sprintf("fig2_%s.svg", id))
REGION_ORDER <- c("V1", "S1", "S1E", "PoCG", "M1", "STG", "SPL", "SMG", "AG", "ITG", "VLPFC", "DLPFC", "FPPFC", "ACC")
LOBE <- c(V1 = "Occipital", S1 = "Parietal", S1E = "Parietal", PoCG = "Parietal", SPL = "Parietal", SMG = "Parietal", AG = "Parietal", STG = "Temporal", ITG = "Temporal", M1 = "Frontal/PFC", VLPFC = "Frontal/PFC", DLPFC = "Frontal/PFC", FPPFC = "Frontal/PFC", ACC = "Limbic")
LOBE_ORDER <- c("Occipital", "Parietal", "Temporal", "Frontal/PFC", "Limbic")
lobe_pal <- c(Occipital = "#4C6EB1", Parietal = "#33A089", Temporal = "#E2A22C", `Frontal/PFC` = "#C44E52", Limbic = "#8064A2")
class_pal <- c(variable = "#C44E52", stable = "#9AA0A6")
ROBUST7 <- as.character(c(1, 3, 4, 6, 8, 10, 14))
SENSITIVE7 <- as.character(c(9, 18, 19, 35, 37, 52, 57))
tier_of <- function(p) factor(ifelse(p %in% ROBUST7, "cohort-robust", ifelse(p %in% SENSITIVE7, "cohort-sensitive", "stable")), levels = c("cohort-robust", "cohort-sensitive", "stable"))
tier_pal <- c(`cohort-robust` = "#B2182B", `cohort-sensitive` = "#F4A582", stable = "#9AA0A6")
div_pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(256)
theme_nat <- function(base = 6.6) {
    theme_classic(base_size = base, base_family = "Helvetica") + theme(axis.line = element_line(linewidth = 0.35, colour = "black"), axis.ticks = element_line(linewidth = 0.35, colour = "black"), axis.title = element_text(size = base), axis.text = element_text(size = base - 0.6, colour = "black"), legend.title = element_text(size = base - 0.4), legend.text = element_text(size = base - 1), legend.key.size = unit(3.2, "mm"), strip.text = element_text(size = base - 0.2, face = "bold"), strip.background = element_blank(), 
        plot.title = element_text(size = base + 1, face = "bold"), plot.subtitle = element_text(size = base - 0.6, colour = "grey30"), panel.grid = element_blank(), plot.background = element_rect(fill = "white", colour = NA))
}
theme_set(theme_nat())
mean_df <- read.table(file.path(RES, "region_program_mean.tsv"), sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
z_df <- read.table(file.path(RES, "program_region_zscore.tsv"), sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
var_df <- read.table(file.path(RES, "program_variability.tsv"), sep = "\t", header = TRUE, check.names = FALSE)
gsub_df <- read.table(file.path(RES, "panel_g_subsample.tsv"), sep = "\t", header = TRUE, check.names = FALSE)
remap <- read.table(file.path(RES, "program_renumber_map.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
remap$new_P_num <- suppressWarnings(as.integer(sub("^P", "", remap$new_P)))
kept_map <- remap[remap$status == "kept", ]
kept_map <- kept_map[order(kept_map$new_P_num), ]
progs <- as.character(kept_map$old_P)
M_mean <- as.matrix(mean_df[, progs])
M_z <- as.matrix(z_df[, progs])
var_df$program <- as.character(var_df$program)
var_df <- var_df[var_df$program %in% progs, ]
var_df$class <- factor(var_df$class, levels = c("variable", "stable"))
variable_progs <- var_df$program[var_df$class == "variable"]
prog_mean_act <- colMeans(M_mean)
prog_cv <- apply(M_mean, 2, function(x) sd(x)/mean(x))
nm_df <- read.table(file.path(RES, "program_names.tsv"), sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
nm_df$program <- as.character(nm_df$cnmf_component)
nm_df$new_P_num <- suppressWarnings(as.integer(sub("^P", "", nm_df$new_P)))
.mk_label <- function(p_old, new_p_num, ns, conf) {
    if (is.na(new_p_num)) 
        return(paste0("EXCLUDED_P", p_old))
    ns <- trimws(ns)
    suf <- paste0(" P", new_p_num)
    if (grepl(paste0(suf, "$"), ns)) 
        ns <- trimws(sub(paste0(suf, "$"), "", ns))
    star <- ifelse(.annotation_group(conf) == "class_b", "*", "")
    paste0("P", new_p_num, " ", ns, star)
}
prog_lab <- setNames(mapply(.mk_label, nm_df$program, nm_df$new_P_num, nm_df$name_short, nm_df$confidence), nm_df$program)
prog_pn <- setNames(ifelse(is.na(nm_df$new_P_num), paste0("P", nm_df$program), paste0("P", nm_df$new_P_num)), nm_df$program)
.lab <- function(p) ifelse(p %in% names(prog_lab), prog_lab[p], paste0("P", p))
.lab_short <- function(p, maxchar = 17L) {
    full <- .lab(p)
    m <- regmatches(full, regexec("^(P[0-9]+)\\s+(.*)$", full))[[1]]
    if (length(m) < 3) 
        return(full)
    pn <- m[2]
    nm <- m[3]
    star <- ifelse(grepl("\\*$", nm), "*", "")
    nm <- sub("\\*$", "", nm)
    if (nchar(nm) > maxchar) {
        cut <- substr(nm, 1, maxchar)
        sp <- regexpr("\\s[^\\s]*$", cut)
        if (sp > 4) 
            cut <- substr(cut, 1, sp - 1)
        nm <- paste0(sub("[\\s[:punct:]]+$", "", cut), "…")
    }
    paste0(pn, " ", nm, star)
}
anat_rank <- setNames(seq_along(REGION_ORDER), REGION_ORDER)
Ha <- t(M_z[REGION_ORDER, ])
Ha <- M_z[REGION_ORDER, progs]
prog_class <- setNames(as.character(var_df$class[match(progs, var_df$program)]), progs)
prog_tier <- setNames(as.character(tier_of(progs)), progs)
prog_eta <- setNames(var_df$eta2_region[match(progs, var_df$program)], progs)
col_lab_a <- vapply(progs, .lab_short, character(1), maxchar = 14L)
prog_tier_fac <- factor(prog_tier[progs], levels = c("cohort-robust", "cohort-sensitive", "stable"))
top_anno <- HeatmapAnnotation(tier = prog_tier, eta2 = anno_barplot(prog_eta, gp = gpar(fill = "#6E7B8B", col = NA), height = unit(7, "mm"), axis_param = list(gp = gpar(fontsize = 5))), col = list(tier = tier_pal), annotation_name_gp = gpar(fontsize = 6), annotation_legend_param = list(tier = list(title = "cohort robustness", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5.5))), simple_anno_size = unit(2.5, "mm"), show_legend = c(tier = TRUE))
left_anno <- rowAnnotation(lobe = LOBE[REGION_ORDER], col = list(lobe = lobe_pal), annotation_name_gp = gpar(fontsize = 6), annotation_legend_param = list(lobe = list(title = "lobe", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5.5))), simple_anno_size = unit(2.5, "mm"))
CELL_A <- 3.4
ht_a <- Heatmap(Ha, name = "z-score", col = colorRamp2(seq(-2.5, 2.5, length.out = 256), div_pal), top_annotation = top_anno, left_annotation = left_anno, cluster_rows = TRUE, cluster_columns = TRUE, column_split = prog_tier_fac, column_gap = unit(c(1.2, 2), "mm"), show_row_dend = TRUE, show_column_dend = TRUE, row_dend_width = unit(6, "mm"), column_dend_height = unit(6, "mm"), row_names_gp = gpar(fontsize = 6), column_labels = col_lab_a, column_names_gp = gpar(fontsize = 4.8), column_names_rot = 90, 
    column_names_max_height = unit(38, "mm"), column_title = "retained programs (n=54)", column_title_gp = gpar(fontsize = 6.5, fontface = "bold"), heatmap_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5.5), legend_height = unit(20, "mm")), row_title = "cortical region", row_title_gp = gpar(fontsize = 6.5), width = ncol(Ha) * unit(CELL_A, "mm"), height = nrow(Ha) * unit(CELL_A, "mm"))
PANEL_A_W <- 13.5
PANEL_A_H <- 6
pdf(file.path(OUT, "panel_a.pdf"), width = PANEL_A_W, height = PANEL_A_H, family = "Helvetica")
draw(ht_a, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
dev.off()
svglite(svgf("a"), width = PANEL_A_W, height = PANEL_A_H, bg = "white")
draw(ht_a, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
invisible(dev.off())
gb_a <- grid.grabExpr(draw(ht_a, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE), width = PANEL_A_W, height = PANEL_A_H)
p_a <- patchwork::wrap_elements(full = gb_a)
b_df <- var_df %>% arrange(desc(eta2_region)) %>% mutate(rank = row_number(), tier = tier_of(program), program = factor(program, levels = rev(program)), ypos = as.integer(program))
b_xmax <- max(b_df$eta2_region)
nlab <- 8L
b_lab <- b_df %>% filter(rank <= nlab) %>% arrange(rank) %>% mutate(lab = .lab(as.character(program)), lab_x = b_xmax * 1.05, lab_y = seq(58, 12, length.out = nlab))
p_b <- ggplot(b_df, aes(x = eta2_region, y = ypos, colour = tier)) + geom_segment(aes(x = 0, xend = eta2_region, yend = ypos), linewidth = 0.4) + geom_point(aes(size = tier == "cohort-robust")) + scale_size_manual(values = c(`TRUE` = 1.9, `FALSE` = 1.1), guide = "none") + geom_segment(data = b_lab, inherit.aes = FALSE, aes(x = eta2_region, y = ypos, xend = lab_x, yend = lab_y), colour = "grey55", linewidth = 0.18) + geom_text(data = b_lab, inherit.aes = FALSE, aes(x = lab_x, y = lab_y, label = lab), 
    hjust = 0, size = 1.8, colour = "black") + scale_colour_manual(values = tier_pal, name = "cohort robustness") + scale_x_continuous(expand = expansion(mult = c(0, 1.9)), name = expression(eta^2 ~ "(region)")) + scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) + labs(y = NULL, title = "Program variability ranking") + theme_nat() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = c(0.82, 0.3), aspect.ratio = 1.6)
ggsave(file.path(OUT, "panel_b.pdf"), p_b, width = 3, height = 3.6, device = cairo_pdf)
svglite(svgf("b"), width = 3, height = 3.6, bg = "white")
plot(p_b)
invisible(dev.off())
c_df <- data.frame(program = progs, mean_act = prog_mean_act[progs], cv = prog_cv[progs], class = var_df$class[match(progs, var_df$program)], tier = tier_of(progs), eta2 = prog_eta[progs])
top_lab <- var_df %>% arrange(desc(eta2_region)) %>% slice_head(n = 5) %>% pull(program)
c_df$lab <- ifelse(c_df$program %in% top_lab, .lab(c_df$program), "")
p_c_main <- ggplot(c_df, aes(mean_act, cv, colour = tier)) + geom_point(aes(size = eta2), alpha = 0.85) + geom_text_repel(aes(label = lab), size = 2, colour = "black", max.overlaps = Inf, segment.size = 0.2, min.segment.length = 0, box.padding = 0.8, point.padding = 0.4, force = 8, force_pull = 0.5, seed = 42) + scale_colour_manual(values = tier_pal, name = "cohort robustness", guide = "none") + scale_size_continuous(range = c(0.4, 2.6), name = expression(eta^2)) + scale_x_continuous(name = "mean activity (across regions)") + 
    scale_y_continuous(name = "CV across regions") + theme_nat() + theme(legend.position = c(0.85, 0.82))
p_c_top <- ggplot(c_df, aes(mean_act, fill = class, colour = class)) + geom_density(alpha = 0.35, linewidth = 0.3) + scale_fill_manual(values = class_pal, guide = "none") + scale_colour_manual(values = class_pal, guide = "none") + theme_void()
p_c_right <- ggplot(c_df, aes(cv, fill = class, colour = class)) + geom_density(alpha = 0.35, linewidth = 0.3) + scale_fill_manual(values = class_pal, guide = "none") + scale_colour_manual(values = class_pal, guide = "none") + coord_flip() + theme_void()
p_c_grp <- (p_c_top + plot_spacer() + p_c_main + p_c_right) + plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4))
ggsave(file.path(OUT, "panel_c.pdf"), p_c_grp, width = 3, height = 3, device = cairo_pdf)
svglite(svgf("c"), width = 3, height = 3, bg = "white")
plot(p_c_grp)
invisible(dev.off())
p_c <- wrap_elements(full = p_c_grp)
d_reg_dist <- dist(M_z[rownames(M_z), progs], method = "euclidean")
hc_reg <- hclust(d_reg_dist, method = "ward.D2")
leaf_order <- hc_reg$labels[hc_reg$order]
axis_rank <- setNames(seq_along(leaf_order), leaf_order)
reg_order_tab <- data.frame(region = leaf_order, leaf_order = seq_along(leaf_order), axis_rank = as.integer(axis_rank[leaf_order]))
write.table(reg_order_tab, file.path(RES, "region_cluster_order.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
grad_clust <- bind_rows(lapply(progs, function(p) {
    v <- M_mean[leaf_order, p]
    x <- axis_rank[leaf_order]
    ct <- suppressWarnings(cor.test(x, v, method = "pearson"))
    sl <- coef(lm(v ~ x))[2]
    data.frame(program = p, axis_slope = unname(sl), axis_r = unname(ct$estimate), axis_p = ct$p.value)
}))
write.table(grad_clust, file.path(RES, "program_gradient_clustered.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
var_grad <- grad_clust %>% filter(program %in% variable_progs) %>% arrange(desc(abs(axis_r)))
nshow <- min(8L, nrow(var_grad))
show_progs <- var_grad$program[seq_len(nshow)]
d_long <- as.data.frame(M_z[leaf_order, show_progs]) %>% mutate(region = factor(leaf_order, levels = leaf_order)) %>% pivot_longer(-region, names_to = "program", values_to = "z")
d_grad <- grad_clust %>% filter(program %in% show_progs)
d_lab <- d_long %>% group_by(program) %>% slice_max(as.integer(region), n = 1) %>% left_join(d_grad, by = "program") %>% mutate(txt = sprintf("%s\n(r=%.2f)", .lab_short(program), axis_r))
ord_progs <- d_grad$program[order(-d_grad$axis_r)]
prog_cols8 <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(show_progs)), ord_progs)
lobe_strip <- data.frame(region = factor(leaf_order, levels = leaf_order), lobe = factor(LOBE[leaf_order], levels = LOBE_ORDER), x = seq_along(leaf_order))
ystrip <- min(c(d_long$z, 0)) - 0.85
ystep <- 0.32
p_d <- ggplot(d_long, aes(region, z, group = program, colour = program)) + geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") + geom_line(linewidth = 0.55) + geom_point(size = 0.7) + geom_tile(data = lobe_strip, inherit.aes = FALSE, aes(x = x, y = ystrip, fill = lobe), width = 1, height = ystep, colour = "white", linewidth = 0.2) + geom_text_repel(data = d_lab, aes(label = txt), hjust = 0, direction = "y", nudge_x = 0.3, size = 1.78, segment.size = 0.18, lineheight = 0.85, xlim = c(14.25, 
    NA), max.overlaps = Inf, force = 3, box.padding = 0.32, point.padding = 0.2, seed = 42) + scale_colour_manual(values = prog_cols8, guide = "none") + scale_fill_manual(values = lobe_pal, name = "lobe", guide = guide_legend(override.aes = list(colour = NA), nrow = 1)) + scale_x_discrete(expand = expansion(mult = c(0.04, 0.3))) + scale_y_continuous(expand = expansion(mult = c(0.16, 0.06))) + labs(x = "region (program-clustered order)", y = "region z-score", title = "Program gradients along clustered region axis", 
    subtitle = "axis = ward.D2 clustering of regions on 54-program profile") + theme_nat() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5), legend.position = "bottom", legend.box = "horizontal", legend.title = element_text(size = 5.6), legend.text = element_text(size = 5.2), legend.key.size = unit(2.4, "mm"), legend.margin = margin(t = 0, b = 0))
ggsave(file.path(OUT, "panel_d.pdf"), p_d, width = 3.3, height = 3.1, device = cairo_pdf)
svglite(svgf("d"), width = 3.3, height = 3.1, bg = "white")
plot(p_d)
invisible(dev.off())
file.copy(file.path(OUT, "panel_d.pdf"), file.path(OUT, "fig2_d.pdf"), overwrite = TRUE)
pca <- prcomp(M_z[REGION_ORDER, progs], center = TRUE, scale. = FALSE)
ve <- (pca$sdev^2)/sum(pca$sdev^2) * 100
scores <- as.data.frame(pca$x[, 1:2])
scores$region <- rownames(scores)
scores$lobe <- factor(LOBE[scores$region], levels = LOBE_ORDER)
load_progs <- var_df %>% arrange(desc(eta2_region)) %>% slice_head(n = 6) %>% pull(program)
ld <- as.data.frame(pca$rotation[load_progs, 1:2])
ld$program <- rownames(ld)
sc_fac <- 0.55 * max(abs(scores$PC1), abs(scores$PC2))/max(abs(ld$PC1), abs(ld$PC2))
ld$PC1s <- ld$PC1 * sc_fac
ld$PC2s <- ld$PC2 * sc_fac
ld$lab <- .lab(ld$program)
p_e <- ggplot(scores, aes(PC1, PC2)) + geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") + geom_segment(data = ld, aes(x = 0, y = 0, xend = PC1s, yend = PC2s), arrow = arrow(length = unit(1.4, "mm")), colour = "grey45", linewidth = 0.3, inherit.aes = FALSE) + geom_text_repel(data = ld, aes(PC1s, PC2s, label = lab), colour = "grey25", size = 1.55, fontface = "italic", segment.size = 0.15, segment.colour = "grey60", min.segment.length = 0, 
    max.overlaps = 30, box.padding = 0.2, inherit.aes = FALSE) + geom_point(aes(fill = lobe), size = 2.1, shape = 21, colour = "white", stroke = 0.3) + geom_text_repel(aes(label = region, colour = lobe), size = 2, max.overlaps = Inf, segment.size = 0.2, min.segment.length = 0, box.padding = 0.5, point.padding = 0.3, force = 6, seed = 42, show.legend = FALSE) + scale_fill_manual(values = lobe_pal, name = "lobe") + scale_colour_manual(values = lobe_pal, guide = "none") + labs(x = sprintf("PC1 (%.0f%%)", 
    ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]), title = "Region similarity (program profile)") + coord_cartesian(clip = "off") + theme_nat() + theme(legend.position = c(0.9, 0.18), legend.background = element_rect(fill = "white", colour = NA), legend.key.size = unit(2.4, "mm"), legend.title = element_text(size = 5.6), legend.text = element_text(size = 5.2), plot.margin = margin(t = 1, r = 3, b = 2, l = 2, unit = "pt"))
ggsave(file.path(OUT, "panel_e.pdf"), p_e, width = 3.4, height = 3, device = cairo_pdf)
svglite(svgf("e"), width = 3.4, height = 3, bg = "white")
plot(p_e)
invisible(dev.off())
top6 <- var_df %>% filter(as.character(program) %in% ROBUST7) %>% arrange(desc(eta2_region)) %>% slice_head(n = 6) %>% pull(program) %>% as.character()
f_lab6 <- setNames(sub("^(P[0-9]+) ", "\\1\n", .lab(top6)), top6)
f_reg_lvls <- rev(REGION_ORDER)
f_df <- as.data.frame(M_z[REGION_ORDER, top6]) %>% mutate(region = factor(REGION_ORDER, levels = f_reg_lvls)) %>% pivot_longer(-region, names_to = "program", values_to = "z") %>% mutate(program = factor(f_lab6[program], levels = f_lab6[top6]), lobe = factor(LOBE[as.character(region)], levels = LOBE_ORDER))
p_f <- ggplot(f_df, aes(x = z, y = region, fill = lobe)) + geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") + geom_segment(aes(x = 0, xend = z, yend = region, colour = lobe), linewidth = 0.45) + geom_point(shape = 21, colour = "white", stroke = 0.25, size = 1.5) + facet_grid(~program) + scale_fill_manual(values = lobe_pal, name = "lobe") + scale_colour_manual(values = lobe_pal, guide = "none") + scale_x_continuous(name = "region z-score", expand = expansion(mult = c(0.06, 0.08))) + 
    labs(y = NULL, title = "Regional activity of top-variable programs") + theme_nat() + theme(axis.text.y = element_text(size = 5, colour = "black"), axis.text.x = element_text(size = 5), strip.text = element_text(size = 5.4, face = "bold", lineheight = 0.9), panel.spacing.x = unit(2.2, "mm"), legend.position = "bottom", legend.key.size = unit(2.6, "mm"))
ggsave(file.path(OUT, "panel_f.pdf"), p_f, width = 7, height = 2.4, device = cairo_pdf)
svglite(svgf("f"), width = 7, height = 2.4, bg = "white")
plot(p_f)
invisible(dev.off())
g_ids <- c("14", "6", "1")
g_subclass <- c(`14` = "L3-L4 IT RORB", `6` = "L2-L3 IT LINC00507", `1` = "L6 IT")
g_facet <- setNames(sprintf("%s\n%s · cohort-%s", .lab(g_ids), g_subclass[g_ids], ifelse(g_ids %in% ROBUST7, "robust", "sensitive")), g_ids)
g_long <- gsub_df %>% mutate(program = as.character(program)) %>% filter(program %in% g_ids, subclass == g_subclass[program]) %>% mutate(program = recode(program, !!!g_facet), program = factor(program, levels = unname(g_facet)), region = factor(region, levels = REGION_ORDER), lobe = factor(LOBE[as.character(region)], levels = LOBE_ORDER))
g_med <- g_long %>% group_by(program, region, lobe) %>% summarise(med = median(act), .groups = "drop")
g_cap <- g_long %>% group_by(program) %>% summarise(ymax = quantile(act, 0.99), .groups = "drop")
g_long <- g_long %>% left_join(g_cap, by = "program")
p_g <- ggplot(g_long, aes(region, act, fill = lobe)) + geom_boxplot(width = 0.7, linewidth = 0.3, colour = "grey25", outlier.size = 0.3, outlier.colour = "grey55", outlier.alpha = 0.4, outlier.stroke = 0, fatten = 1.4) + geom_point(data = g_med, aes(region, med), inherit.aes = FALSE, size = 0.4, colour = "black") + facet_wrap(~program, ncol = 2, scales = "free_y") + scale_fill_manual(values = lobe_pal, name = "lobe") + labs(x = NULL, y = "per-cell program activity", title = "Spotlight programs within selected subclasses") + 
    theme_nat() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), strip.text = element_text(size = 5.2, face = "bold"), legend.position = "none")
ggsave(file.path(OUT, "panel_g.pdf"), p_g, width = 4.2, height = 3.2, device = cairo_pdf)
svglite(svgf("g"), width = 4.2, height = 3.2, bg = "white")
plot(p_g)
invisible(dev.off())
SIG_K <- 2L
.sig_topk <- function(reg, k, up = TRUE) {
    v <- M_z[reg, progs]
    names(sort(v, decreasing = up))[seq_len(k)]
}
sig_union <- character(0)
sig_top1 <- setNames(character(length(REGION_ORDER)), REGION_ORDER)
for (reg in REGION_ORDER) {
    ups <- .sig_topk(reg, SIG_K, TRUE)
    dns <- .sig_topk(reg, SIG_K, FALSE)
    sig_top1[reg] <- ups[1]
    for (p in c(ups, dns)) if (!(p %in% sig_union)) 
        sig_union <- c(sig_union, p)
}
.colkey <- function(p) {
    br <- names(which.max(M_z[REGION_ORDER, p]))
    c(anat_rank[[br]], -M_z[br, p])
}
sig_cols <- sig_union[order(sapply(sig_union, function(p) .colkey(p)[1]), sapply(sig_union, function(p) .colkey(p)[2]))]
Hh <- M_z[REGION_ORDER, sig_cols, drop = FALSE]
colnames(Hh) <- sig_cols
h_collab <- .lab(sig_cols)
sig_mark <- matrix(FALSE, nrow = nrow(Hh), ncol = ncol(Hh), dimnames = dimnames(Hh))
for (reg in REGION_ORDER) sig_mark[reg, sig_top1[reg]] <- TRUE
CELL_H <- 2.6
left_anno_h <- rowAnnotation(lobe = LOBE[REGION_ORDER], col = list(lobe = lobe_pal), annotation_name_gp = gpar(fontsize = 6), show_legend = FALSE, simple_anno_size = unit(2.5, "mm"))
ht_h <- Heatmap(Hh, name = "z-score", col = colorRamp2(seq(-2.5, 2.5, length.out = 256), div_pal), left_annotation = left_anno_h, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 5.6), column_names_side = "top", column_names_rot = 90, column_labels = h_collab, column_names_gp = gpar(fontsize = 5), rect_gp = gpar(col = "white", lwd = 0.4), cell_fun = function(j, i, x, y, w, h, fill) {
    if (sig_mark[i, j]) 
        grid.rect(x, y, w, h, gp = gpar(col = "black", lwd = 0.9, fill = NA))
}, row_title = "cortical region", row_title_gp = gpar(fontsize = 6.5), column_title = "region-distinctive program (top-2 up/down per region; ring = region's #1 up)", column_title_side = "bottom", column_title_gp = gpar(fontsize = 5.6), heatmap_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5.5), legend_height = unit(15, "mm")), width = ncol(Hh) * unit(CELL_H, "mm"), height = nrow(Hh) * unit(CELL_H, "mm"))
PANEL_H_W <- 6.6
PANEL_H_H <- 3.6
pdf(file.path(OUT, "panel_h.pdf"), width = PANEL_H_W, height = PANEL_H_H, family = "Helvetica")
draw(ht_h, heatmap_legend_side = "right", padding = unit(c(0, 0, 0, 0), "mm"))
dev.off()
svglite(svgf("h"), width = PANEL_H_W, height = PANEL_H_H, bg = "white")
draw(ht_h, heatmap_legend_side = "right", padding = unit(c(0, 0, 0, 0), "mm"))
invisible(dev.off())
gb_h <- grid.grabExpr(draw(ht_h, heatmap_legend_side = "right", padding = unit(c(0, 0, 0, 0), "mm")), width = PANEL_H_W, height = PANEL_H_H)
p_h <- patchwork::wrap_elements(full = gb_h)
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
display_map <- read.delim(file.path(dirname(dirname(normalizePath(script_file))), "08_external_annotations", "inputs", "program_display_map.tsv"), stringsAsFactors = FALSE)
sig_supp <- lapply(REGION_ORDER, function(reg) {
    raw_z <- z_df[reg, as.character(1:60)]
    raw_ranked <- function(up) {
        data.frame(raw_program = names(sort(raw_z, decreasing = up))[seq_len(5L)], rank = seq_len(5L), stringsAsFactors = FALSE)
    }
    selected <- bind_rows(transform(raw_ranked(TRUE), direction = "up"), transform(raw_ranked(FALSE), direction = "down")) %>% filter(raw_program %in% progs)
    selected <- selected %>% left_join(display_map, by = c(raw_program = "program")) %>% mutate(region = reg, lobe = unname(LOBE[reg]), program = as.integer(display_program), program_name = paste0(new_P, " ", name_short, ifelse(as.integer(display_marker) == 1L, "*", "")), z = round(raw_z[raw_program], 3))
    selected[, c("region", "lobe", "direction", "rank", "program", "program_name", "z")]
}) %>% bind_rows()
write.table(sig_supp, file.path(RES, "supp_table_region_signatures.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
EXCLUDED_CNMF <- as.character(remap$old_P[remap$status != "kept"])
top10 <- var_df %>% filter(!program %in% EXCLUDED_CNMF) %>% arrange(desc(eta2_region)) %>% slice_head(n = 10) %>% pull(program)
lobe_of_region <- LOBE[REGION_ORDER]
i_long <- as.data.frame(M_mean[REGION_ORDER, top10]) %>% mutate(region = REGION_ORDER, lobe = lobe_of_region[REGION_ORDER]) %>% pivot_longer(c(-region, -lobe), names_to = "program", values_to = "act") %>% group_by(lobe, program) %>% summarise(act = mean(act), .groups = "drop") %>% group_by(lobe) %>% mutate(rank = rank(-act, ties.method = "first")) %>% ungroup() %>% mutate(lobe = factor(lobe, levels = LOBE_ORDER), program = factor(program, levels = top10))
i_pal <- setNames(colorRampPalette(brewer.pal(10, "Spectral"))(length(top10)), top10)
nP <- length(top10)
i_long <- i_long %>% mutate(rankf = factor(rank, levels = nP:1))
rank_breaks <- (nP - (1:nP)) + 0.5
rank_labels <- paste0("rank ", 1:nP)
i_long <- i_long %>% mutate(lab_pn = unname(prog_pn[as.character(program)]))
p_i <- ggplot(i_long, aes(x = lobe, y = 1, stratum = rankf, alluvium = program, fill = program, label = lab_pn)) + geom_alluvium(width = 0.26, alpha = 0.62, colour = "white", linewidth = 0.12, curve_type = "sigmoid") + geom_stratum(width = 0.26, fill = NA, colour = "grey80", linewidth = 0.15) + geom_text(stat = "stratum", size = 2.7, colour = "grey10", fontface = "bold") + scale_fill_manual(values = i_pal, name = NULL, labels = .lab(top10)) + scale_y_continuous(name = "activity rank within lobe (1 = highest)", 
    breaks = rank_breaks, labels = rank_labels, expand = expansion(mult = c(0.02, 0.02))) + labs(x = NULL, title = "Program rank shifts across lobes") + guides(fill = guide_legend(nrow = 2, byrow = TRUE, keywidth = unit(3, "mm"), keyheight = unit(3, "mm"))) + theme_nat() + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5), axis.text.y = element_text(size = 5, colour = "grey30"), axis.ticks.y = element_line(linewidth = 0.25), legend.position = "bottom", legend.direction = "horizontal", 
    legend.box = "horizontal", legend.key.size = unit(3, "mm"), legend.text = element_text(size = 5), legend.margin = margin(t = 1, b = 0))
ggsave(file.path(OUT, "panel_i.pdf"), p_i, width = 4.4, height = 3.7, device = cairo_pdf)
svglite(svgf("i"), width = 4.4, height = 3.7, bg = "white")
plot(p_i)
invisible(dev.off())
design <- "\nAAAAAA\nBCCDDE\nFFFFFF\nGGGHHH\nGGGIII\n"
fig <- p_a + p_b + p_c + p_d + p_e + p_f + p_g + p_h + p_i + plot_layout(design = design, heights = c(1.85, 1.35, 0.78, 1, 1)) + plot_annotation(tag_levels = "a", title = "Fig. 3  Cross-region variation of cortical gene programs", theme = theme(plot.title = element_text(size = 10, face = "bold", family = "Helvetica"))) & theme(plot.tag.position = c(0, 1), plot.tag = element_text(face = "bold", size = 9, family = "Helvetica", hjust = 0, vjust = 1))
ggsave(file.path(OUT, "fig2_crossregion.pdf"), fig, width = 13, height = 13, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT, "fig2_crossregion.png"), fig, width = 13, height = 13, dpi = 200, bg = "white", limitsize = FALSE)
file.copy(file.path(OUT, "fig2_crossregion.pdf"), file.path(OUT, "fig2_full.pdf"), overwrite = TRUE)
rep_grad <- var_grad %>% slice_head(n = nshow) %>% left_join(var_df[, c("program", "eta2_region")], by = "program") %>% select(program, eta2_region, axis_slope, axis_r, axis_p)
write.table(rep_grad, file.path(OUT, "panel_d_gradient_numbers.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

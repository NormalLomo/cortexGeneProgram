PROJECT_ROOT <- Sys.getenv("CORTEX_PROGRAM_ROOT", unset = "/DATA/cortex_nmf_program")
suppressPackageStartupMessages({
    library(ggplot2)
    library(svglite)
})
base <- file.path(PROJECT_ROOT, "results/xspecies_humanmap_v1/figures/cross_species_nuclear_projection")
data_d <- file.path(base, "data")
out_d <- file.path(base, "svg_panels")
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)
pairs <- read.csv(file.path(data_d, "panelA_pairs_long.csv"), stringsAsFactors = FALSE)
meds <- read.csv(file.path(data_d, "panelA_medians.csv"), stringsAsFactors = FALSE)
ord <- c("human-mouse", "human-macaque", "mouse-macaque")
lab <- c(`human-mouse` = "human–mouse", `human-macaque` = "human–macaque", `mouse-macaque` = "mouse–macaque")
pairs$pair <- factor(pairs$pair, levels = ord)
meds$pair <- factor(meds$pair, levels = ord)
meds <- meds[order(meds$pair), ]
pair_cols <- c(`human-mouse` = "#3E8E7E", `human-macaque` = "#C8743C", `mouse-macaque` = "#9A7B4F")
ff <- "Liberation Sans"
meds$ylab <- meds$median
lab_txt <- sprintf("%.2f", meds$median)
set.seed(11)
p <- ggplot(pairs, aes(x = pair, y = cosine, fill = pair, colour = pair)) + geom_violin(width = 0.85, alpha = 0.85, linewidth = 0.4, trim = FALSE, scale = "width", show.legend = FALSE) + geom_jitter(width = 0.07, height = 0, size = 0.55, alpha = 0.5, stroke = 0, show.legend = FALSE) + geom_boxplot(width = 0.12, fill = "white", colour = "#3A3A3A", linewidth = 0.4, outlier.shape = NA, show.legend = FALSE) + geom_text(data = meds, aes(x = pair, y = 1, label = sprintf("%.2f", median)), inherit.aes = FALSE, 
    family = ff, size = 5.5/.pt, colour = "#3A3A3A", vjust = 0) + scale_fill_manual(values = pair_cols) + scale_colour_manual(values = pair_cols) + scale_x_discrete(labels = lab, expand = expansion(add = c(0.55, 0.55))) + scale_y_continuous(limits = c(0, 1.06), breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(mult = c(0.01, 0))) + labs(y = "Cosine similarity (conservation)") + theme_minimal(base_family = ff) + theme(text = element_text(family = ff, colour = "#3A3A3A"), axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 6.5, colour = "#3A3A3A", margin = margin(r = 1)), axis.text.x = element_text(size = 5.5, colour = "#3A3A3A"), axis.text.y = element_text(size = 5.5, colour = "#3A3A3A"), axis.line = element_line(linewidth = 0.4, colour = "#3A3A3A"), axis.ticks = element_line(linewidth = 0.4, colour = "#3A3A3A"), axis.ticks.length = unit(1.2, "pt"), panel.grid.major.y = element_line(linewidth = 0.25, colour = "#E8E8E8"), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), 
    legend.position = "none", plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent", colour = NA), plot.margin = margin(0, 0, 0, 0))
W <- 2.35
H <- 1.85
svglite(file.path(out_d, "a.svg"), width = W, height = H, bg = "transparent", system_fonts = list(sans = ff))
invisible(dev.off())

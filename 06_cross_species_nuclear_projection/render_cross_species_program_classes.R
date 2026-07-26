PROJECT_ROOT <- Sys.getenv("CORTEX_PROGRAM_ROOT", unset = "/DATA/cortex_nmf_program")
suppressMessages({
    library(ggplot2)
    library(data.table)
    library(scales)
    library(svglite)
})
CJK <- "Liberation Sans"
BASE <- file.path(PROJECT_ROOT, "results/xspecies_humanmap_v1")
D <- file.path(BASE, "figures", "cross_species_nuclear_projection", "data")
OUTD <- file.path(BASE, "figures", "cross_species_nuclear_projection", "svg_panels")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
PT <- 1/2.845276
gtsize <- function(pt) pt * 0.3528
sp_cols <- c(mouse = "#3E8E7E", macaque = "#C8743C")
med_cols <- c(mouse = "#2F6B5E", macaque = "#A85C2A")
GREY_AX <- "#3A3A3A"
GREY_NOTE <- "#5A5A5A"
GREY_REF <- "#E0E0E0"
GREY_DUMB <- "#B8B8B8"
zt <- theme_minimal(base_size = 6, base_family = CJK) + theme(text = element_text(family = CJK, color = GREY_AX), plot.margin = margin(0, 0, 0, 0), panel.grid.minor = element_blank(), panel.grid.major.x = element_line(linewidth = 0.18, color = "#EDEDED"), panel.grid.major.y = element_blank(), axis.line.x = element_line(linewidth = 0.4, color = GREY_AX), axis.ticks.x = element_line(linewidth = 0.3, color = GREY_AX), axis.ticks.length = unit(1.2, "pt"), axis.text.x = element_text(size = 5.5, color = "grey20"), 
    axis.text.y = element_text(size = 5.5, color = "grey15"), axis.title.x = element_text(size = 6.5, color = GREY_AX), plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = c(0.985, 0.14), legend.justification = c(1, 0), legend.title = element_text(size = 5.5, face = "bold"), legend.text = element_text(size = 5, color = "grey20"), legend.key.size = unit(4, "pt"), legend.spacing.y = unit(0.5, "pt"), legend.background = element_rect(fill = "white", color = NA), legend.margin = margin(1, 
        2, 1, 2))
cc <- fread(file.path(D, "panelC_v2_strat.csv"))
cmed <- fread(file.path(D, "panelC_v2_inner_med.csv"))
cst <- fread(file.path(D, "panelC_v2_stat.csv"))
inner_lv <- cmed$inner
cc$inner <- factor(cc$inner, levels = inner_lv)
cmed$inner <- factor(cmed$inner, levels = inner_lv)
cc$species <- ifelse(cc$species == "鼠", "mouse", "macaque")
cc$species <- factor(cc$species, levels = c("mouse", "macaque"))
ybase <- setNames(seq_along(inner_lv), inner_lv)
off <- 0.17
cc$ynum <- ybase[as.character(cc$inner)] + ifelse(cc$species == "mouse", -off, off)
cmed$ynum_m <- ybase[as.character(cmed$inner)] - off
cmed$ynum_q <- ybase[as.character(cmed$inner)] + off
cc$is_unres <- cc$inner == "unresolved"
fmt_p <- function(p) if (p < 0.001) sprintf("%.1e", p) else sprintf("%.3f", p)
pC <- ggplot() + geom_hline(yintercept = ybase, linewidth = 0.22, color = GREY_REF) + geom_segment(data = cmed, aes(x = mou_med, xend = mac_med, y = ynum_m, yend = ynum_q), linewidth = 0.5, color = GREY_DUMB, lineend = "round") + geom_point(data = cc[low_conf == FALSE & is_unres == FALSE], aes(cosine, ynum, color = species), size = 0.85, alpha = 0.8) + geom_point(data = cc[is_unres == TRUE], aes(cosine, ynum, color = species), size = 0.85, alpha = 0.4, shape = 17) + geom_point(data = cc[low_conf == 
    TRUE & is_unres == FALSE], aes(cosine, ynum, color = species), size = 0.85, alpha = 0.42, shape = 17) + geom_point(data = cmed, aes(mou_med, ynum_m, fill = "mouse"), shape = 23, size = 1.9, color = med_cols["mouse"], stroke = 0.4) + geom_point(data = cmed, aes(mac_med, ynum_q, fill = "macaque"), shape = 23, size = 1.9, color = med_cols["macaque"], stroke = 0.4) + scale_color_manual(values = sp_cols, name = "Species", guide = guide_legend(override.aes = list(size = 1.5, alpha = 1))) + scale_fill_manual(values = sp_cols, 
    guide = "none") + scale_y_continuous(breaks = ybase, labels = inner_lv, limits = c(0.45, length(inner_lv) + 0.62), expand = c(0, 0)) + scale_x_continuous(limits = c(0, 0.83), breaks = c(0, 0.2, 0.4, 0.6, 0.8), expand = expansion(mult = c(0, 0.004))) + annotate("text", x = 0.82, y = length(inner_lv) + 0.42, label = "more conserved →", size = gtsize(5), color = "grey45", hjust = 1, fontface = "italic") + annotate("text", x = 0.005, y = 3.4, hjust = 0, vjust = 0.5, size = gtsize(5), color = GREY_NOTE, 
    lineheight = 1.05, label = sprintf("neuron vs non-neuron\nMann–Whitney\nhuman–mouse p=%s ***\nhuman–macaque p=%s (trend)", fmt_p(cst$mou_p[1]), fmt_p(cst$mac_p[1]))) + labs(x = "Cross-species alignment cosine", y = NULL) + coord_cartesian(clip = "off") + zt
ggsave(file.path(OUTD, "c.svg"), pC, width = 2.55, height = 2.32, units = "in", device = svglite::svglite)

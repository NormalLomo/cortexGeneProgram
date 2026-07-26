.annotation_group <- function(value) ifelse(endsWith(as.character(value), "sig"), "class_a", "class_b")
PROJECT_ROOT <- Sys.getenv("CORTEX_PROGRAM_ROOT", unset = "/DATA/cortex_nmf_program")
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(readr)
    library(svglite)
    library(stringr)
})
base_dir <- file.path(PROJECT_ROOT, "results/xspecies_humanmap_v1/figures/cross_species_nuclear_projection")
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(base_dir, "svg_panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
names_authority <- file.path(PROJECT_ROOT, "results/crossregion_v1/program_names.tsv")
nm <- read_tsv(names_authority, show_col_types = FALSE)
new2name <- nm %>% transmute(new_id = new_P, name_short = name_short)
auth_map <- new2name %>% transmute(program = new_id, new_id, name_short)
lol <- read_csv(file.path(data_dir, "panelB_lollipop.csv"), show_col_types = FALSE)
curc <- read_csv(file.path(data_dir, "func_class_curated_current.csv"), show_col_types = FALSE)
df <- lol %>% select(program, h_mac_cosine, order) %>% left_join(curc %>% select(program, function_class), by = "program") %>% left_join(auth_map, by = "program")
df <- df %>% filter(!is.na(new_id))
agg_map <- c(neuron_synapse = "Synaptic", neuron_ion_channel = "Synaptic", neuron_axon_guidance = "Synaptic", neuron_neuropeptide = "Neuropeptide/IEG", neuron_activity_IEG = "Neuropeptide/IEG", glia_astrocyte = "Glia (oligo/astro)", glia_oligo_myelin = "Glia (oligo/astro)", glia_microglia_immune = "Microglia/immune", vascular = "Vascular", cytoskeleton = "Cytoskeleton", metabolic_housekeeping = "Cytoskeleton", other_unresolved = "Unresolved")
df$big_class <- agg_map[df$function_class]
df$big_class[is.na(df$big_class)] <- "Unresolved"
class_cols <- c(Synaptic = "#3B6FB6", `Neuropeptide/IEG` = "#6B4E9E", `Glia (oligo/astro)` = "#3E8E7E", Vascular = "#C8743C", `Microglia/immune` = "#B0506A", Cytoskeleton = "#8A7A52", Unresolved = "#9AA0A6")
class_order_top <- c("Microglia/immune", "Cytoskeleton", "Glia (oligo/astro)", "Synaptic", "Vascular", "Neuropeptide/IEG", "Unresolved")
class_n <- df %>% group_by(big_class) %>% summarise(n = n(), .groups = "drop")
class_levels <- class_order_top[class_order_top %in% df$big_class]
df$big_class <- factor(df$big_class, levels = rev(class_levels))
n_by <- setNames(class_n$n, class_n$big_class)
class_lab <- setNames(sprintf("%s  (n=%d)", class_levels, n_by[class_levels]), class_levels)
set.seed(11)
df <- df %>% group_by(big_class) %>% mutate(yj = as.numeric(big_class) + (runif(n()) - 0.5) * 0.52) %>% ungroup()
ord_b <- df$program[order(df$h_mac_cosine)]
named_bot <- as.character(head(ord_b, 6))
named_top <- as.character(tail(ord_b, 6))
named_set <- c(named_top, named_bot)
df$is_named <- df$program %in% named_set
df$nm_lab <- ifelse(df$is_named, paste0(df$new_id, " ", str_trunc(df$name_short, 22, ellipsis = "…")), NA_character_)
class_med_pts <- df %>% group_by(big_class) %>% summarise(med = median(h_mac_cosine), .groups = "drop") %>% mutate(yy = as.numeric(big_class))
xmin <- 0
xmax <- 0.75
pt_to_mm <- 0.352777
ax_col <- "#3A3A3A"
th <- theme_minimal(base_size = 6.5, base_family = "Helvetica") + theme(plot.margin = margin(0, 0, 0, 0), panel.grid.major.x = element_line(colour = "#E8E8E8", linewidth = 0.25 * pt_to_mm), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = ax_col, linewidth = 0.4 * pt_to_mm), axis.ticks.x = element_line(colour = ax_col, linewidth = 0.4 * pt_to_mm), axis.ticks.length = unit(1.4, "pt"), axis.ticks.y = element_blank(), axis.title.x = element_text(size = 6.5, 
    colour = ax_col, margin = margin(t = 1.5)), axis.title.y = element_blank(), axis.text.x = element_text(size = 5.5, colour = ax_col), axis.text.y = element_text(size = 5.5, colour = ax_col, hjust = 1), legend.position = "none", plot.title = element_blank(), plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent", colour = NA))
ny <- length(class_levels)
p <- ggplot(df, aes(x = h_mac_cosine, y = yj)) + geom_hline(yintercept = seq(1.5, ny - 0.5, by = 1), colour = "#ECECEC", linewidth = 0.25 * pt_to_mm) + geom_point(aes(colour = big_class), size = 1.8 * pt_to_mm * 2.83, alpha = 0.9, stroke = 0) + geom_point(data = class_med_pts, aes(x = med, y = yy), shape = 23, size = 2.2 * pt_to_mm * 2.83, fill = "white", colour = ax_col, stroke = 0.4 * pt_to_mm) + geom_text_repel(data = df %>% filter(is_named), aes(label = nm_lab, colour = big_class), size = 5 * 
    pt_to_mm, family = "Helvetica", segment.size = 0.3 * pt_to_mm, segment.colour = "#B5B5B5", box.padding = 0.16, point.padding = 0.1, min.segment.length = 0, max.overlaps = Inf, hjust = 1, nudge_x = -0.012, direction = "both", force = 2.2, force_pull = 0.6, seed = 11, show.legend = FALSE) + scale_colour_manual(values = class_cols) + scale_x_continuous(limits = c(xmin, xmax), breaks = c(0, 0.25, 0.5, 0.75), expand = expansion(mult = c(0.02, 0.04))) + scale_y_continuous(breaks = seq_len(ny), labels = class_lab[levels(df$big_class)], 
    limits = c(0.4, ny + 0.6), expand = expansion(mult = c(0, 0))) + labs(x = "human→macaque cosine") + th
H_in <- 2.35
W_in <- H_in * 1.3
svglite(file.path(out_dir, "b.svg"), width = W_in, height = H_in, bg = "transparent")
invisible(dev.off())

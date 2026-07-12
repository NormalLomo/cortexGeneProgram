#!/usr/bin/env Rscript
# Fig. 6 panel g program-by-laminar heatmap.
# Functional-class blocks and a compact top annotation strip organize the panel.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(svglite)
  library(readr)
})

base_dir <- "CORTEX_PROGRAM_ROOT/results/xspecies_humanmap_v1/figures/Fig10_v2"
csv  <- file.path(base_dir, "data", "panelG_program_laminar.csv")
func_csv <- file.path(base_dir, "data", "func_class_curated.csv")
renum_p <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1/program_renumber_map.tsv"
out  <- Sys.getenv("FIG10_PANELG_OUT", file.path(base_dir, "svg_panels", "g.svg"))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

df <- read.csv(csv, stringsAsFactors = FALSE)
func_df <- read.csv(func_csv, stringsAsFactors = FALSE, check.names = FALSE)
renum <- read_tsv(renum_p, show_col_types = FALSE)

o2n <- setNames(
  paste0("P", renum$new_P[renum$status == "kept"]),
  paste0("P", renum$old_P[renum$status == "kept"])
)
df$new_id <- o2n[df$program]
stopifnot(!any(is.na(df$new_id)))

collapse_func_class <- function(fc) {
  if (fc %in% c("glia_astrocyte", "glia_oligo_myelin")) return("Glia")
  if (fc == "glia_microglia_immune") return("Immune")
  if (fc == "vascular") return("Vascular")
  if (fc == "cytoskeleton") return("Cyto")
  if (fc == "neuron_axon_guidance") return("Axon")
  if (fc %in% c("neuron_synapse", "neuron_ion_channel")) return("Syn")
  if (fc %in% c("neuron_neuropeptide", "neuron_activity_IEG")) return("NP/IEG")
  "Other"
}

func_df$broad <- vapply(func_df$function_class, collapse_func_class, character(1))
func_df$new_id <- o2n[func_df$program]
func_df <- func_df[!is.na(func_df$new_id), , drop = FALSE]
func_df$new_num <- as.integer(sub("^P", "", func_df$new_id))

group_order <- c("Glia", "Immune", "Vascular", "Cyto", "Axon", "Syn", "NP/IEG", "Other")
func_df$group_rank <- match(func_df$broad, group_order)
func_df <- func_df[order(func_df$group_rank, -func_df$h_mac_cosine, func_df$new_num), , drop = FALSE]
prog_levels <- unique(func_df$new_id)
nprog <- length(prog_levels)

df$x <- match(df$new_id, prog_levels)

depth_bottom <- c("WM", "L6", "L5", "L4", "L3", "L2", "L1", "ARACHNOID")
depth_labels_bottom <- c("WM", "L6", "L5", "L4", "L3", "L2", "L1", "AR")
df$y <- match(df$majorDomain, depth_bottom)

strip_cols <- c(
  Glia = "#B4684C",
  Immune = "#4F8F73",
  Vascular = "#B28A46",
  Cyto = "#7B6B43",
  Axon = "#4B6A8C",
  Syn = "#4F76B7",
  "NP/IEG" = "#8B6EBA",
  Other = "#888888"
)

group_runs <- rle(func_df$broad)
group_ends <- cumsum(group_runs$lengths)
group_starts <- c(1, head(group_ends + 1, -1))
group_centers <- (group_starts + group_ends) / 2
annot_df <- data.frame(
  group = group_runs$values,
  xmin = group_starts - 0.5,
  xmax = group_ends + 0.5,
  x = group_centers,
  width = group_runs$lengths,
  stringsAsFactors = FALSE
)
label_map <- c(
  Glia = "Glia", Immune = "Imm", Vascular = "Vasc", Cyto = "Cyto",
  Axon = "Axon", Syn = "Syn", "NP/IEG" = "NP/IEG", Other = "Other"
)
annot_df$label <- unname(label_map[annot_df$group])
vline_pos <- head(group_ends + 0.5, -1)

tick_breaks <- sort(unique(c(seq(1, nprog, by = 5), group_starts)))
tick_labels <- prog_levels[tick_breaks]

lim <- 0.85
div_cols <- c(
  "#1B3A4B", "#2D5A6E", "#6E9BAD", "#DCE6E8",
  "#F2EDE4", "#E8B58C", "#C8743C", "#9A4E22"
)

p <- ggplot(df, aes(x = x, y = y, fill = score)) +
  annotate(
    "rect",
    xmin = annot_df$xmin,
    xmax = annot_df$xmax,
    ymin = 8.48,
    ymax = 8.82,
    fill = strip_cols[annot_df$group],
    colour = NA
  ) +
  annotate(
    "text",
    x = annot_df$x[annot_df$width >= 3],
    y = 9.12,
    label = annot_df$label[annot_df$width >= 3],
    family = "sans",
    fontface = "bold",
    size = 1.7,
    colour = strip_cols[annot_df$group[annot_df$width >= 3]]
  ) +
  geom_vline(
    xintercept = vline_pos,
    colour = "#6F6F6F",
    linewidth = 0.55 / .pt
  ) +
  geom_tile(colour = "white", linewidth = 0.30 / .pt) +
  scale_fill_gradientn(
    colours = div_cols,
    limits = c(-lim, lim),
    oob = scales::squish,
    breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
    name = NULL,
    guide = guide_colourbar(
      barwidth = unit(16, "mm"),
      barheight = unit(2.1, "mm"),
      frame.colour = "#3A3A3A",
      frame.linewidth = 0.38 / .pt,
      ticks.colour = "#3A3A3A",
      ticks.linewidth = 0.38 / .pt,
      direction = "horizontal"
    )
  ) +
  scale_x_continuous(
    breaks = tick_breaks,
    labels = tick_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(depth_labels_bottom),
    labels = depth_labels_bottom,
    expand = c(0, 0),
    limits = c(0.5, 9.24)
  ) +
  coord_equal(expand = FALSE, clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 6) +
  theme(
    text = element_text(family = "sans", colour = "#3A3A3A"),
    axis.text.x.bottom = element_text(
      size = 4.8, angle = 90, hjust = 1, vjust = 0.5,
      margin = margin(t = 0.45, unit = "mm")
    ),
    axis.text.y = element_text(
      size = 5.4, hjust = 1,
      margin = margin(r = 0.55, unit = "mm")
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 4.7),
    legend.key.size = unit(2.7, "mm"),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.margin = margin(0, 0, 0.4, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.position = "top",
    legend.justification = "center",
    plot.margin = margin(0, 0, 0, 0)
  )

cell_mm <- 1.50
w_in <- (nprog * cell_mm + 4.0) / 25.4
h_in <- (8 * cell_mm + 4.4 + 2.9) / 25.4

svglite(out, width = w_in, height = h_in, bg = "transparent")
print(p)
invisible(dev.off())

cat("WROTE", out, "\n")
cat("grid", nprog, "cols x 8 rows\n")
cat("group blocks:", paste(annot_df$group, annot_df$width, sep = ":", collapse = " | "), "\n")
cat("tick labels:", paste(tick_labels, collapse = ","), "\n")
cat("aspect_w_in", round(w_in, 3), "h_in", round(h_in, 3),
    "ratio_w/h", round(w_in / h_in, 3), "\n")

#!/usr/bin/env Rscript
# Fig. 6 panel f program-by-subclass enrichment heatmap.
# Explicit subclass grouping, visible boundaries, and a compact legend organize
# the retained-program values.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(svglite)
})

base_dir <- "CORTEX_PROGRAM_ROOT/results/xspecies_humanmap_v1/figures/Fig10_v2"
csv      <- file.path(base_dir, "data", "panelF_program_subclass.csv")
func_csv <- file.path(base_dir, "data", "func_class_curated.csv")
out_svg  <- Sys.getenv("FIG10_PANELF_OUT", file.path(base_dir, "svg_panels", "f.svg"))
dir.create(dirname(out_svg), showWarnings = FALSE, recursive = TRUE)

df <- read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
func_df <- read.csv(func_csv, stringsAsFactors = FALSE, check.names = FALSE)

prog12 <- c("P54","P37","P26","P13","P40","P49","P56","P7","P24","P15","P1","P8")
df <- df[df$program %in% prog12, , drop = FALSE]

old2new_f <- c(
  P54 = "P49", P37 = "P33", P26 = "P23", P13 = "P12",
  P40 = "P36", P49 = "P45", P56 = "P51", P7 = "P7",
  P24 = "P21", P15 = "P14", P1 = "P1", P8 = "P8"
)

collapse_func_class <- function(fc) {
  if (fc %in% c("glia_astrocyte", "glia_oligo_myelin")) return("Glia")
  if (fc == "glia_microglia_immune") return("Immune")
  if (fc == "vascular") return("Vascular")
  if (fc == "cytoskeleton") return("Cytoskeleton")
  if (fc == "neuron_axon_guidance") return("Axon")
  if (fc %in% c("neuron_synapse", "neuron_ion_channel")) return("Synaptic")
  if (fc %in% c("neuron_neuropeptide", "neuron_activity_IEG")) return("NP/IEG")
  if (fc == "metabolic_housekeeping") return("Housekeeping")
  "Other"
}

func_df$broad_class <- vapply(func_df$function_class, collapse_func_class, character(1))
row_class_map <- setNames(func_df$broad_class, func_df$program)

subclass_info <- data.frame(
  subclass = c(
    "AST", "OLIGO", "OPC", "MICRO",
    "ENDO", "VLMC",
    "ET", "L2-L3 IT LINC00507", "L3-L4 IT RORB", "L4-L5 IT RORB",
    "L6 CAR3", "L6 CT", "L6 IT", "L6B", "NP",
    "CHANDELIER", "LAMP5", "NDNF", "PAX6", "PVALB", "SST", "VIP"
  ),
  group = c(
    rep("Glia", 3), "Immune",
    rep("Vascular", 2),
    rep("Exc", 9),
    rep("Inh", 7)
  ),
  stringsAsFactors = FALSE
)

group_cols <- c(
  Glia = "#B4684C",
  Immune = "#4F8F73",
  Vascular = "#B28A46",
  Exc = "#5A82BD",
  Inh = "#8B6EBA"
)
row_text_cols <- c(
  Glia = "#7C4535",
  Immune = "#2E6F57",
  Vascular = "#7A6431",
  Cytoskeleton = "#7B6B43",
  Axon = "#4B6A8C",
  Synaptic = "#345EA5",
  "NP/IEG" = "#7355A3",
  Housekeeping = "#6A6A6A",
  Other = "#7A7A7A"
)

subclass_order <- subclass_info$subclass
prog_top <- prog12
prog_bottom <- rev(prog_top)

df <- df[match(prog_bottom, df$program, nomatch = 0) > 0 | TRUE, , drop = FALSE]
df$x <- match(df$subclass, subclass_order)
df$y <- match(df$program, prog_bottom)
df <- df[!is.na(df$x) & !is.na(df$y), , drop = FALSE]

ylabs_new <- unname(old2new_f[prog_bottom])
y_text_cols <- unname(row_text_cols[row_class_map[prog_bottom]])
y_text_cols[is.na(y_text_cols)] <- row_text_cols["Other"]

group_runs <- rle(subclass_info$group)
group_ends <- cumsum(group_runs$lengths)
group_starts <- c(1, head(group_ends + 1, -1))
group_centers <- (group_starts + group_ends) / 2
annot_df <- data.frame(
  group = group_runs$values,
  xmin = group_starts - 0.5,
  xmax = group_ends + 0.5,
  x = group_centers,
  stringsAsFactors = FALSE
)
annot_df$label <- c("Glia", "Imm", "Vasc", "Exc", "Inh")
vline_pos <- head(group_ends + 0.5, -1)

div_cols <- c("#1B4F6A", "#7FA6BA", "#F4F1EA", "#E0A98C", "#B0502E")
ff <- "Liberation Sans"
col_axis <- "#3A3A3A"
nrow_plot <- length(prog_bottom)

p <- ggplot(df, aes(x = x, y = y, fill = z)) +
  annotate(
    "rect",
    xmin = annot_df$xmin,
    xmax = annot_df$xmax,
    ymin = nrow_plot + 0.58,
    ymax = nrow_plot + 0.90,
    fill = group_cols[annot_df$group],
    colour = NA
  ) +
  annotate(
    "text",
    x = annot_df$x,
    y = nrow_plot + 1.20,
    label = annot_df$label,
    family = ff,
    fontface = "bold",
    size = 1.8,
    colour = unname(group_cols[annot_df$group])
  ) +
  geom_vline(
    xintercept = vline_pos,
    colour = "#6F6F6F",
    linewidth = 0.55 / .pt
  ) +
  geom_tile(colour = "white", linewidth = 0.32 / .pt) +
  scale_fill_gradientn(
    colours = div_cols,
    limits = c(-2, 2),
    breaks = c(-2, -1, 0, 1, 2),
    oob = scales::squish,
    name = "Enrichment z",
    guide = guide_colourbar(
      direction = "horizontal",
      barwidth = unit(17, "mm"),
      barheight = unit(2.4, "mm"),
      frame.colour = col_axis,
      frame.linewidth = 0.38 / .pt,
      ticks.colour = col_axis,
      ticks.linewidth = 0.38 / .pt,
      title.position = "top"
    )
  ) +
  scale_x_continuous(
    breaks = seq_along(subclass_order),
    labels = subclass_order,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_len(nrow_plot),
    labels = ylabs_new,
    expand = c(0, 0),
    limits = c(0.5, nrow_plot + 1.34)
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = ff) +
  theme(
    text = element_text(family = ff, colour = col_axis),
    axis.text.x = element_text(
      size = 4.8, angle = 52, hjust = 1, vjust = 1,
      colour = col_axis, lineheight = 0.92,
      margin = margin(t = 0.5, unit = "mm")
    ),
    axis.text.y = element_text(size = 5.3, colour = y_text_cols),
    axis.ticks = element_line(colour = col_axis, linewidth = 0.38 / .pt),
    axis.ticks.length = unit(0.55, "mm"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 5.3, face = "bold", colour = col_axis),
    legend.text = element_text(size = 4.8, colour = col_axis),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.5, "mm"),
    plot.margin = margin(3.0, 0, 0, 0, unit = "mm")
  )

mm <- function(x) x / 25.4
ggsave(
  out_svg, p,
  device = svglite::svglite,
  width = mm(118),
  height = mm(61),
  units = "in",
  bg = "transparent"
)

cat(sprintf("WROTE %s\n", out_svg))
cat(sprintf("rows(programs)=%d cols(subclass)=%d cells=%d\n",
            nrow_plot, length(subclass_order), nrow(df)))
cat(sprintf("column groups: %s\n", paste(annot_df$group, collapse = " | ")))
cat(sprintf("row order (top->bottom): %s\n",
            paste(old2new_f[prog_top], collapse = " ")))

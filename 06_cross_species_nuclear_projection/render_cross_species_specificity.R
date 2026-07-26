PROJECT_ROOT <- Sys.getenv("CORTEX_PROGRAM_ROOT", unset = "/DATA/cortex_nmf_program")
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(svglite)
    library(systemfonts)
})
ttc_reg <- "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
ttc_bold <- "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"
register_font(name = "FigSans", plain = ttc_reg, bold = ttc_bold, italic = ttc_reg, bolditalic = ttc_bold)
FONT <- "FigSans"
base_dir <- file.path(PROJECT_ROOT, "results/xspecies_humanmap_v1/figures/cross_species_nuclear_projection")
csv_path <- file.path(base_dir, "data", "panelE_significance.csv")
out_dir <- file.path(base_dir, "svg_panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_svg <- file.path(out_dir, "e.svg")
df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
df$h_mac_cosine <- as.numeric(df$h_mac_cosine)
df$h_mou_cosine <- as.numeric(df$h_mou_cosine)
df$new_id <- df$program
df <- df[order(-df$h_mac_cosine), ]
n <- nrow(df)
df$x <- seq_len(n)
long <- bind_rows(data.frame(x = df$x, y = 2, sp = "human → macaque", cos = df$h_mac_cosine, rank = df$h_mac_diag_rank, stringsAsFactors = FALSE), data.frame(x = df$x, y = 1, sp = "human → mouse", cos = df$h_mou_cosine, rank = df$h_mou_diag_rank, stringsAsFactors = FALSE))
long$best <- long$rank == 1
ramp_cols <- c("#E4E5DB", "#C8D6CA", "#ABC6B9", "#8EB6A8", "#70A697", "#509686", "#3A887B", "#327C76", "#2A7171", "#256268", "#21525C", "#1D4250")
df$highlight <- (df$x <= 6) | (df$x >= (n - 5))
lab_df <- df %>% transmute(x = x, lab = new_id, highlight = highlight)
tick_df <- lab_df
tri_df <- df %>% filter(highlight) %>% transmute(x = x)
hdr_df <- data.frame(y = c(2, 1), lab = c("human → macaque", "human → mouse"))
xlim <- c(-9, n + 0.5)
ylim <- c(0.3, 3)
p <- ggplot() + geom_tile(data = long, aes(x = x, y = y, fill = cos), color = "white", linewidth = 0.3, width = 0.96, height = 0.96) + geom_point(data = subset(long, best), aes(x = x, y = y), shape = 21, fill = "white", color = "#1B3A4B", size = 0.9, stroke = 0.45) + geom_point(data = tri_df, aes(x = x, y = 2.74), shape = 25, fill = "#1B1B1B", colour = "#1B1B1B", size = 1.05, stroke = 0) + geom_text(data = hdr_df, aes(x = 0.2, y = y, label = lab), hjust = 1, vjust = 0.5, size = 6/.pt, family = FONT, 
    fontface = "bold", colour = "#1B1B1B") + geom_text(data = subset(lab_df, !highlight), aes(x = x, y = 0.34, label = lab), angle = 90, hjust = 1, vjust = 0.5, size = 5/.pt, family = FONT, colour = "#9AA0A6") + geom_text(data = subset(lab_df, highlight), aes(x = x, y = 0.34, label = lab), angle = 90, hjust = 1, vjust = 0.5, size = 5.5/.pt, family = FONT, fontface = "bold", colour = "#1B1B1B") + scale_fill_gradientn(colours = ramp_cols, limits = c(0, 0.75), breaks = c(0, 0.25, 0.5, 0.75), name = NULL, 
    guide = guide_colourbar(barwidth = unit(26, "mm"), barheight = unit(2.2, "mm"), ticks.colour = "#3A3A3A", frame.colour = "#3A3A3A", frame.linewidth = 0.4, direction = "horizontal")) + coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, clip = "off", expand = FALSE) + theme_void(base_size = 5.5) + theme(text = element_text(family = FONT, colour = "#3A3A3A"), legend.title = element_blank(), legend.text = element_text(size = 5, colour = "#3A3A3A"), legend.position = "top", legend.justification = "center", 
    legend.background = element_rect(fill = NA, colour = NA), legend.margin = margin(0, 0, 1, 0, unit = "mm"), legend.box.margin = margin(0, 0, 0, 0), legend.key.size = unit(3, "mm"), plot.margin = margin(0, 0, 0, 0))
mm2in <- 1/25.4
svglite(out_svg, width = 246 * mm2in, height = 35 * mm2in, bg = "transparent")
invisible(dev.off())
svg_txt <- readLines(out_svg, warn = FALSE)
svg_txt <- gsub("font-family: *\"?Noto Sans CJK [A-Z]{2}\"?;", "font-family: Helvetica Neue, Arial, sans-serif;", svg_txt)
writeLines(svg_txt, out_svg)

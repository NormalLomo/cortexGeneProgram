#!/usr/bin/env Rscript
# Fig.10 panel e — HORIZONTAL conservation-specificity banner (zero-margin vector SVG)
# Current contract: a horizontal 2-row x 54-column conservation banner using
# retained program IDs ordered by decreasing human-to-macaque cosine.
#   The viridis-teal cosine scale uses limits c(0,0.75).
#     端點: cos0=#F1ECE3 (淺米) -> cos0.75=#1B3B4B (深 teal); 驗證 cos0.715=#1D4351
#   欄序 = panelE_significance.csv (54 row, mac-cosine 降序, 已剔 6 excluded).
#   Labels use retained IDs from program_renumber_map.tsv; source IDs are provenance only.
#   best-reciprocal 白圈 (diag_rank==1) 保留.
#   科學數值 (cosine/sig/rank) 一律不動, 只改排版方向 + 顏色映射(還原) + 編號.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(svglite)
  library(systemfonts)
})

# stable Latin-capable family so svglite writes a real font-family (no tofu)
ttc_reg  <- "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
ttc_bold <- "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"
register_font(name = "FigSans", plain = ttc_reg, bold = ttc_bold,
              italic = ttc_reg, bolditalic = ttc_bold)
FONT <- "FigSans"

base_dir <- "CORTEX_PROGRAM_ROOT/results/xspecies_humanmap_v1/figures/Fig10_v2"
csv_path <- file.path(base_dir, "data", "panelE_significance.csv")
renum_p  <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1/program_renumber_map.tsv"
out_dir  <- file.path(base_dir, "svg_panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_svg  <- file.path(out_dir, "e.svg")

# ---- read 54-row significance (handles quoted commas in func_short) ----
df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(df) == 54)                       # 60 - 6 cohort-technical excluded
df$h_mac_cosine <- as.numeric(df$h_mac_cosine)
df$h_mou_cosine <- as.numeric(df$h_mou_cosine)

# ---- source P-id (= csv program) -> retained P-id (kept rows) ----
renum <- read_tsv(renum_p, show_col_types = FALSE)
o2n <- setNames(paste0("P", renum$new_P[renum$status == "kept"]),
                paste0("P", renum$old_P[renum$status == "kept"]))
df$new_id <- o2n[df$program]
stopifnot(!any(is.na(df$new_id)))               # all 54 must map

# ---- column order: mac-cosine descending (as file) -> x = 1..54 left..right ----
df <- df[order(-df$h_mac_cosine), ]
n  <- nrow(df)
df$x <- seq_len(n)                              # 1 (highest mac) .. 54

# ---- long: 2 rows (macaque y=2 top, mouse y=1 bottom) ----
long <- bind_rows(
  data.frame(x = df$x, y = 2, sp = "human → macaque",
             cos = df$h_mac_cosine, rank = df$h_mac_diag_rank,
             stringsAsFactors = FALSE),
  data.frame(x = df$x, y = 1, sp = "human → mouse",
             cos = df$h_mou_cosine, rank = df$h_mou_diag_rank,
             stringsAsFactors = FALSE)
)
long$best <- long$rank == 1                     # diag_rank==1 = self best reciprocal

# ---- viridis-teal ramp sampled from master colorbar image (12 stops) ----
ramp_cols <- c("#E4E5DB","#C8D6CA","#ABC6B9","#8EB6A8","#70A697","#509686",
               "#3A887B","#327C76","#2A7171","#256268","#21525C","#1D4250")

# ---- column labels: ALL 54 programs labelled BELOW the band, in TWO TIERS
#      Highlighted = the 6 highest-mac
#      (head) + 6 lowest-mac (tail) columns -> bold black + a ▼ marker ABOVE the
#      band; the other 42 columns -> grey, no marker. Shown as NEW ids. ----
df$highlight <- (df$x <= 6) | (df$x >= (n - 5))   # head6 + tail6 (12 cols)
lab_df <- df %>% transmute(x = x, lab = new_id, highlight = highlight)
tick_df <- lab_df                                   # back-compat for the final cat()
# ▼ down-triangle markers (above the top band row), one per highlighted column
tri_df  <- df %>% filter(highlight) %>% transmute(x = x)

# Row header labels are right-aligned left of column 1.
hdr_df <- data.frame(y = c(2, 1), lab = c("human → macaque", "human → mouse"))

# ---- geometry: square cells; x in [1,54], y in {1,2}. left band for row
#      headers, top room for ▼ markers + colorbar, bottom room for the rotated
#      two-tier program labels. ----
xlim <- c(-9.0, n + 0.5)        # left band hosts the two right-aligned row headers
# tight band ylim (cells stay square via coord_fixed); ▼ markers just above the
# top row. Rotated labels (below) + colorbar (top legend strip) overflow this
# range via clip="off" and are given physical room by the svglite canvas height.
ylim <- c(0.3, 3.0)

p <- ggplot() +
  geom_tile(data = long, aes(x = x, y = y, fill = cos),
            color = "white", linewidth = 0.3, width = 0.96, height = 0.96) +
  geom_point(data = subset(long, best), aes(x = x, y = y),
             shape = 21, fill = "white", color = "#1B3A4B",
             size = 0.9, stroke = 0.45) +
  # Down-triangle markers above highlighted columns.
  geom_point(data = tri_df, aes(x = x, y = 2.74),
             shape = 25, fill = "#1B1B1B", colour = "#1B1B1B",
             size = 1.05, stroke = 0) +
  # Bold row headers, right-aligned just left of column 1.
  geom_text(data = hdr_df, aes(x = 0.2, y = y, label = lab),
            hjust = 1, vjust = 0.5, size = 6 / .pt, family = FONT,
            fontface = "bold", colour = "#1B1B1B") +
  # program labels (rotated 90deg, BELOW the bottom row), TWO TIERS, NEW ids:
  #   highlighted -> bold black; others -> grey.
  geom_text(data = subset(lab_df, !highlight), aes(x = x, y = 0.34, label = lab),
            angle = 90, hjust = 1, vjust = 0.5, size = 5 / .pt,
            family = FONT, colour = "#9AA0A6") +
  geom_text(data = subset(lab_df, highlight), aes(x = x, y = 0.34, label = lab),
            angle = 90, hjust = 1, vjust = 0.5, size = 5.5 / .pt,
            family = FONT, fontface = "bold", colour = "#1B1B1B") +
  scale_fill_gradientn(
    colours = ramp_cols, limits = c(0, 0.75),
    breaks = c(0, 0.25, 0.5, 0.75), name = NULL,
    guide = guide_colourbar(
      barwidth = unit(26, "mm"), barheight = unit(2.2, "mm"),
      ticks.colour = "#3A3A3A", frame.colour = "#3A3A3A",
      frame.linewidth = 0.4, direction = "horizontal"
    )
  ) +
  coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, clip = "off", expand = FALSE) +
  theme_void(base_size = 5.5) +
  theme(
    text = element_text(family = FONT, colour = "#3A3A3A"),
    # The colorbar uses ticks only.
    legend.title  = element_blank(),
    legend.text   = element_text(size = 5, colour = "#3A3A3A"),
    # Center the colorbar above the banner so it remains inside the compose box.
    legend.position   = "top",
    legend.justification = "center",
    legend.background = element_rect(fill = NA, colour = NA),
    legend.margin     = margin(0, 0, 1, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key.size   = unit(3, "mm"),
    plot.margin = margin(0, 0, 0, 0)
  )

# ---- export zero-margin vector SVG (wide short banner) ----
# 54 cols x ~ (2 rows + headers/ticks/colorbar). Wide aspect (~banner).
# Match the placed panel aspect while retaining square tiles with coord_fixed.
mm2in <- 1 / 25.4
# Canvas aspect matches the placed panel box so compose FILL_WIDTH keeps
# the height within the banner box. The layout contains
# stack = colorbar (top) + 2-row band + ▼ + two-tier labels (below).
svglite(out_svg, width = 246 * mm2in, height = 35 * mm2in, bg = "transparent")
print(p)
invisible(dev.off())

# normalise font-family (EN-only; robust Latin fallback chain)
svg_txt <- readLines(out_svg, warn = FALSE)
svg_txt <- gsub("font-family: *\"?Noto Sans CJK [A-Z]{2}\"?;",
                "font-family: Helvetica Neue, Arial, sans-serif;", svg_txt)
writeLines(svg_txt, out_svg)

cat("WROTE:", out_svg, "\n")
cat("cols:", n, "(mac-desc) | rows: 2 (human->macaque, human->mouse)\n")
cat("highlighted (head6+tail6, bold+▼):",
    paste(subset(lab_df, highlight)$lab, collapse = ","), "\n")
cat("n triangles:", nrow(tri_df), "| all labels below band, two-tier\n")
cat("ramp endpoints: cos0=", ramp_cols[1], " cos0.75=", tail(ramp_cols,1), "\n")

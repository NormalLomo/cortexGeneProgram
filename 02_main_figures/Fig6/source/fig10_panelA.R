#!/usr/bin/env Rscript
# Fig.10 panel a — 三物種程式保守度 violin (zero-margin vector SVG)
# 依 ART_SPEC §1.1/§1.4/§1.5/§2.a 重生成。鐵律：科學數據/措辭不動，只改美學。
suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
})

base   <- "CORTEX_PROGRAM_ROOT/results/xspecies_humanmap_v1/figures/Fig10_v2"
data_d <- file.path(base, "data")
out_d  <- file.path(base, "svg_panels")
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

# ---- 資料 ----
pairs <- read.csv(file.path(data_d, "panelA_pairs_long.csv"), stringsAsFactors = FALSE)
meds  <- read.csv(file.path(data_d, "panelA_medians.csv"),    stringsAsFactors = FALSE)

# pair 順序：由低到高中位 → human-mouse < human-macaque < mouse-macaque
ord <- c("human-mouse", "human-macaque", "mouse-macaque")
lab <- c("human-mouse" = "human–mouse",
         "human-macaque" = "human–macaque",
         "mouse-macaque" = "mouse–macaque")
pairs$pair <- factor(pairs$pair, levels = ord)
meds$pair  <- factor(meds$pair,  levels = ord)
meds <- meds[order(meds$pair), ]

# §1.1 pair 色：human-mouse=青瓷綠 / human-macaque=赭橙 / mouse-macaque=金棕
pair_cols <- c("human-mouse"   = "#3E8E7E",
               "human-macaque" = "#C8743C",
               "mouse-macaque" = "#9A7B4F")

# §1.4 字族（拉丁主字族；中文靠 Noto fallback，svglite 內嵌處理）
ff <- "Liberation Sans"

# 中位數標籤位置（box 上方留空）
meds$ylab <- meds$median
lab_txt <- sprintf("%.2f", meds$median)

# ---- 繪圖 ----
set.seed(11)
p <- ggplot(pairs, aes(x = pair, y = cosine, fill = pair, colour = pair)) +
  geom_violin(width = 0.85, alpha = 0.85, linewidth = 0.4,
              trim = FALSE, scale = "width", show.legend = FALSE) +
  geom_jitter(width = 0.07, height = 0, size = 0.55, alpha = 0.5,
              stroke = 0, show.legend = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", colour = "#3A3A3A",
               linewidth = 0.4, outlier.shape = NA, show.legend = FALSE) +
  # 中位數值標於 violin 上方
  geom_text(data = meds,
            aes(x = pair, y = 1.0, label = sprintf("%.2f", median)),
            inherit.aes = FALSE, family = ff, size = 5.5 / .pt,
            colour = "#3A3A3A", vjust = 0) +
  scale_fill_manual(values = pair_cols) +
  scale_colour_manual(values = pair_cols) +
  scale_x_discrete(labels = lab, expand = expansion(add = c(0.55, 0.55))) +
  scale_y_continuous(limits = c(0, 1.06),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     expand = expansion(mult = c(0.01, 0.0))) +
  labs(y = "Cosine similarity (conservation)") +
  theme_minimal(base_family = ff) +
  theme(
    text             = element_text(family = ff, colour = "#3A3A3A"),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(size = 6.5, colour = "#3A3A3A",
                                    margin = margin(r = 1)),
    axis.text.x      = element_text(size = 5.5, colour = "#3A3A3A"),
    axis.text.y      = element_text(size = 5.5, colour = "#3A3A3A"),
    axis.line        = element_line(linewidth = 0.4, colour = "#3A3A3A"),
    axis.ticks       = element_line(linewidth = 0.4, colour = "#3A3A3A"),
    axis.ticks.length = unit(1.2, "pt"),
    panel.grid.major.y = element_line(linewidth = 0.25, colour = "#E8E8E8"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    plot.background    = element_rect(fill = "transparent", colour = NA),
    panel.background   = element_rect(fill = "transparent", colour = NA),
    plot.margin        = margin(0, 0, 0, 0)
  )

# ---- zero-margin 向量 SVG ----
# 目標 aspect ≈ 1.27 (W:H) wider landscape，與 b(1.288) 近等 → a/b 並重頭條
W <- 2.35  # inch
H <- 1.85  # inch  → W/H ≈ 1.27
svglite(file.path(out_d, "a.svg"), width = W, height = H,
        bg = "transparent", system_fonts = list(sans = ff))
print(p)
invisible(dev.off())

cat("DONE panel a ->", file.path(out_d, "a.svg"), "\n")
cat(sprintf("aspect W:H = %.2f : %.2f  (H/W=%.3f)\n", W, H, H/W))
cat("medians (order low->high):\n")
print(meds[, c("pair", "median")])

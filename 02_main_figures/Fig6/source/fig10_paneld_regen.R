#!/usr/bin/env Rscript
# Fig10_v2 panel d 重生 — zero-margin 向量 SVG
# 依 ART_SPEC: 4 連橫排 gene-loading scatter (P8/P1/P15/P58)，人 K60 loading(x) vs 投影物種 refit loading(y)
# 物種色 macaque 赭橙 / mouse 青瓷；點 0.6pt alpha 0.3；回歸 1pt + CI；cosine 標角；功能名小標
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(patchwork); library(svglite); library(grid)
})

base_dir <- "CORTEX_PROGRAM_ROOT/results/xspecies_humanmap_v1/figures/Fig10_v2"
data_dir <- file.path(base_dir, "data")
out_dir  <- file.path(base_dir, "svg_panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

FONT <- "Liberation Sans"  # Arial/Helvetica 替代

# 物種色（ART_SPEC §1.1）
COL_MAC <- "#C8743C"; COL_MOU <- "#3E8E7E"
COL_AXIS <- "#3A3A3A"; COL_REF <- "#B8B8B8"; COL_NOTE <- "#5A5A5A"
sz <- function(pt) pt / .pt

# ---- 讀資料 ----
sc <- read_csv(file.path(data_dir, "panelD_gene_scatter.csv"), show_col_types = FALSE)
cc <- read_csv(file.path(data_dir, "panelD_case_cos.csv"),     show_col_types = FALSE)

prog_order <- c("P8", "P1", "P15", "P58")   # CSV source cNMF IDs
# Source-to-retained provenance mapping (program_renumber_map.tsv):
#   P8->P8, P1->P1, P15->P14, P58->P52 (cNMF-58 = new P52 Neuropeptide signaling).
# 第4 範例 cNMF-58(=舊P58) 是低保守對照, 新號 P52。
old2new_d <- c(P8 = "P8", P1 = "P1", P15 = "P14", P58 = "P52")
sc$program <- factor(sc$program, levels = prog_order)
short_map <- setNames(cc$func_short, cc$program)
cos_mac   <- setNames(cc$cos_mac,   cc$program)
cos_mou   <- setNames(cc$cos_mou,   cc$program)
sc$species <- factor(sc$species, levels = c("macaque", "mouse"))
sp_cols <- c(macaque = COL_MAC, mouse = COL_MOU)

xr <- range(sc$human_loading, na.rm = TRUE)
yr <- range(sc$species_refit, na.rm = TRUE)

make_panel <- function(pg, show_y = FALSE, show_leg = FALSE) {
  d  <- sc %>% filter(program == pg)
  cm <- cos_mac[[pg]]; co <- cos_mou[[pg]]
  fn <- short_map[[pg]]; fn <- sub("^P[0-9]+ ", "", fn); fn <- sub(",?\\s*$", "", fn)
  # 斷詞避免溢出：每行約 ≤13 字元（窄格，積極換行）
  wrap_w <- 13
  words <- strsplit(fn, " ")[[1]]; lines <- c(); cur <- ""
  for (w in words) { test <- if (cur=="") w else paste(cur, w)
    if (nchar(test) > wrap_w) { lines <- c(lines, cur); cur <- w } else cur <- test }
  lines <- c(lines, cur); fn <- paste(lines, collapse = "\n")
  # brain-weak programs would get a star suffix on the id (weak/non-significant
  # brain GO:BP enrichment; FDR high). All four exemplars here (P8/P1/P15/P58)
  # are brain-sig (confidence=brain-sig in program_names.tsv) -> no star.
  weak <- character(0)
  # Display the retained ID while preserving the functional label from the CSV.
  pg_new  <- old2new_d[[pg]]
  pg_disp <- if (pg %in% weak) paste0(pg_new, "*") else pg_new
  ttl <- paste0(pg_disp, " ", fn)

  p <- ggplot(d, aes(x = human_loading, y = species_refit, color = species)) +
    geom_hline(yintercept = 0, color = COL_REF, linewidth = 0.25) +
    geom_point(size = 0.6, alpha = 0.30, stroke = 0, shape = 16) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.0,
                aes(fill = species), alpha = 0.15) +
    scale_color_manual(values = sp_cols, name = NULL,
                       guide = if (show_leg) guide_legend(override.aes=list(size=1.6, alpha=1)) else "none") +
    scale_fill_manual(values = sp_cols, guide = "none") +
    scale_x_continuous(limits = xr, expand = expansion(mult = c(0.02, 0.02)),
                       breaks = scales::pretty_breaks(3)) +
    scale_y_continuous(limits = yr, expand = expansion(mult = c(0.02, 0.02)),
                       breaks = scales::pretty_breaks(4)) +
    annotate("text", x = xr[2], y = yr[2], hjust = 1, vjust = 1.1,
             label = sprintf("mac %.2f", cm), color = COL_MAC, size = sz(5.5), family = FONT) +
    annotate("text", x = xr[2], y = yr[2], hjust = 1, vjust = 2.5,
             label = sprintf("mou %.2f", co), color = COL_MOU, size = sz(5.5), family = FONT) +
    labs(title = ttl, x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_classic(base_family = FONT) +
    theme(
      plot.title = element_text(size = 8, face = "bold", family = FONT, hjust = 0.5,
                                margin = margin(0,0,1,0), lineheight = 0.95),
      axis.line = element_line(color = COL_AXIS, linewidth = 0.4),
      axis.ticks = element_line(color = COL_AXIS, linewidth = 0.4),
      axis.ticks.length = unit(1.2, "pt"),
      axis.text.x = element_text(size = 5.5, color = COL_AXIS, family = FONT),
      axis.text.y = if (show_y) element_text(size = 5.5, color = COL_AXIS, family = FONT) else element_blank(),
      axis.ticks.y = if (show_y) element_line(color = COL_AXIS, linewidth = 0.4) else element_blank(),
      panel.grid = element_blank(), panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = if (show_leg) c(0.5, -0.02) else "none",
      legend.direction = "horizontal",
      legend.text = element_text(size = 5, family = FONT),
      legend.key.size = unit(5, "pt"), legend.margin = margin(0,0,0,0),
      legend.background = element_blank(),
      aspect.ratio = 1,
      plot.margin = margin(0, 0, 0, 0)
    )
  p
}

p8  <- make_panel("P8",  show_y = TRUE,  show_leg = FALSE)
p1  <- make_panel("P1",  show_y = FALSE, show_leg = FALSE)
p15 <- make_panel("P15", show_y = FALSE, show_leg = FALSE)
p58 <- make_panel("P58", show_y = FALSE, show_leg = FALSE)

# 小格間距，避免標題相撞（仍緊湊）
row <- p8 + p1 + p15 + p58 + plot_layout(nrow = 1) &
  theme(plot.margin = margin(0, 2.5, 0, 2.5))  # 左右各 2.5pt gutter，給標題斷行空間

# 共享 y / x 軸標題（外層 grid 包，zero-margin）+ 底部水平 legend
W_mm <- 124; H_mm <- 72
yt <- textGrob("projected-species loading", rot = 90,
               gp = gpar(fontsize = 6.5, fontfamily = FONT, col = COL_AXIS))
xt <- textGrob("human K60 program loading",
               gp = gpar(fontsize = 6.5, fontfamily = FONT, col = COL_AXIS))

# 手繪緊湊水平 legend grob（避免脆弱 guide-box 抽取）
draw_legend <- function() {
  # 兩物種色點 + 標籤，水平排，置中
  vp <- viewport(layout = grid.layout(1, 4,
        widths = unit.c(unit(7,"pt"), unit(38,"pt"), unit(7,"pt"), unit(30,"pt"))))
  pushViewport(viewport(width = unit(82,"pt"), height = unit(1,"npc")))
  pushViewport(vp)
  pushViewport(viewport(layout.pos.col=1)); grid.points(0.5,0.5, pch=16, size=unit(4,"pt"),
        gp=gpar(col=COL_MAC)); popViewport()
  pushViewport(viewport(layout.pos.col=2)); grid.text("macaque", x=0.05, hjust=0,
        gp=gpar(fontsize=5, fontfamily=FONT, col=COL_AXIS)); popViewport()
  pushViewport(viewport(layout.pos.col=3)); grid.points(0.5,0.5, pch=16, size=unit(4,"pt"),
        gp=gpar(col=COL_MOU)); popViewport()
  pushViewport(viewport(layout.pos.col=4)); grid.text("mouse", x=0.05, hjust=0,
        gp=gpar(fontsize=5, fontfamily=FONT, col=COL_AXIS)); popViewport()
  popViewport(); popViewport()
}

svglite(file.path(out_dir, "d.svg"), width = W_mm/25.4, height = H_mm/25.4,
        bg = "transparent", fix_text_size = FALSE)
grid.newpage()
lay <- grid.layout(nrow = 2, ncol = 2,
                   widths  = unit.c(unit(4.5, "mm"), unit(1, "null")),
                   heights = unit.c(unit(1, "null"), unit(5, "mm")))
pushViewport(viewport(layout = lay, gp = gpar(fontfamily = FONT)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1)); grid.draw(yt); popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2)); print(row, newpage = FALSE); popViewport()
# 底列：左中 x 標題、右 legend
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  pushViewport(viewport(x = 0.38, width = 0.55)); grid.draw(xt); popViewport()
  pushViewport(viewport(x = 0.85, width = 0.28)); draw_legend(); popViewport()
popViewport()
popViewport()
invisible(dev.off())

cat("DONE panel d ->", file.path(out_dir, "d.svg"), "\n")
cat(sprintf("W=%dmm H=%dmm aspect=%.3f:1\n", W_mm, H_mm, W_mm/H_mm))

#!/usr/bin/env Rscript
# figure-fine Stage1: 重生 Fig10_v2 panel c -> zero-margin 向量 SVG
# 依 ART_SPEC.md (§1.1 物種色 / §1.3 / §1.4 字>=5pt / §1.5 線寬點徑 / panel c spec)
# 鐵律: 人 K60 program 投影鼠猴, 零 joint, 三物種對全留, 數字/措辭一律不動.
# 只改美學: 棄廉價紅藍 -> 物種色(鼠=青瓷綠 #3E8E7E / 猴=赭橙 #C8743C);
#           類中位 人-鼠 + 人-猴 各一菱形 + dumbbell 連線; zero-margin; 緊湊 legend.
suppressMessages({
  library(ggplot2); library(data.table); library(scales); library(svglite)
})

CJK  <- "Liberation Sans"   # Latin only (英文圖, 無 CJK)
BASE <- "CORTEX_PROGRAM_ROOT/results/xspecies_humanmap_v1"
D    <- file.path(BASE, "figures", "Fig10_v2", "data")
OUTD <- file.path(BASE, "figures", "Fig10_v2", "svg_panels")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

# ---- pt -> ggplot geom size 換算 (1pt 字 ~ size 0.3528 in geom_text; 線寬 pt 直接給 linewidth*?) ----
PT <- 1/2.845276  # mm per pt 的倒數不需; 用 size = pt * 0.3528 給 geom_text
gtsize <- function(pt) pt * 0.3528  # geom_text size from pt

# ---- ART_SPEC §1.1 物種色 ----
sp_cols <- c("mouse" = "#3E8E7E",   # mouse 青瓷綠
             "macaque" = "#C8743C") # macaque 赭橙
# 類中位菱形: 人-鼠 用鼠色深, 人-猴 用猴色深 (與點同語言, 描邊深灰)
med_cols <- c("mouse" = "#2F6B5E", "macaque" = "#A85C2A")
GREY_AX  <- "#3A3A3A"; GREY_NOTE <- "#5A5A5A"; GREY_REF <- "#E0E0E0"; GREY_DUMB <- "#B8B8B8"

# ---- theme: zero margin, sans, >=5pt (ART_SPEC §1.4 / §1.5) ----
zt <- theme_minimal(base_size = 6, base_family = CJK) +
  theme(
    text             = element_text(family = CJK, color = GREY_AX),
    plot.margin      = margin(0, 0, 0, 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.18, color = "#EDEDED"),
    panel.grid.major.y = element_blank(),
    axis.line.x      = element_line(linewidth = 0.4, color = GREY_AX),
    axis.ticks.x     = element_line(linewidth = 0.3, color = GREY_AX),
    axis.ticks.length= unit(1.2, "pt"),
    axis.text.x      = element_text(size = 5.5, color = "grey20"),
    axis.text.y      = element_text(size = 5.5, color = "grey15"),
    axis.title.x     = element_text(size = 6.5, color = GREY_AX),
    plot.title       = element_blank(),     # 標題交給 figure-fine compose 層
    plot.subtitle    = element_blank(),
    legend.position  = c(0.985, 0.14),
    legend.justification = c(1, 0),
    legend.title     = element_text(size = 5.5, face = "bold"),
    legend.text      = element_text(size = 5,   color = "grey20"),
    legend.key.size  = unit(4, "pt"),
    legend.spacing.y = unit(0.5, "pt"),
    legend.background= element_rect(fill = "white", color = NA),
    legend.margin    = margin(1, 2, 1, 2)
  )

# ===================== 讀資料 =====================
cc   <- fread(file.path(D, "panelC_v2_strat.csv"))      # program x 鼠/猴 x cosine, inner, low_conf
cmed <- fread(file.path(D, "panelC_v2_inner_med.csv"))  # 7 內類: mou_med, mac_med
cst  <- fread(file.path(D, "panelC_v2_stat.csv"))       # MWU p (mou_p, mac_p), n_neuron/n_nonneuron

inner_lv <- cmed$inner   # CSV 已按 mouse 中位升序 (least -> most conserved)
cc$inner    <- factor(cc$inner,   levels = inner_lv)
cmed$inner  <- factor(cmed$inner, levels = inner_lv)
# CSV species 值為 中文 鼠/猴 -> remap 英文 (科學數據不動, 僅標籤)
cc$species  <- ifelse(cc$species == "鼠", "mouse", "macaque")
cc$species  <- factor(cc$species, levels = c("mouse", "macaque"))

# y 數值座標 + 物種偏移 (鼠下 / 猴上), 兩物種同列可分辨
ybase     <- setNames(seq_along(inner_lv), inner_lv)
off       <- 0.17
cc$ynum   <- ybase[as.character(cc$inner)] + ifelse(cc$species == "mouse", -off, off)
cmed$ynum_m <- ybase[as.character(cmed$inner)] - off  # 鼠中位
cmed$ynum_q <- ybase[as.character(cmed$inner)] + off  # 猴中位

# 低證據類 (unresolved) 整體弱化
cc$is_unres <- cc$inner == "unresolved"

fmt_p <- function(p) if (p < 1e-3) sprintf("%.1e", p) else sprintf("%.3f", p)

# ===================== panel c =====================
pC <- ggplot() +
  # 橫向類別參考線 (極淡, 0.25pt)
  geom_hline(yintercept = ybase, linewidth = 0.22, color = GREY_REF) +
  # dumbbell: 連同類 鼠中位 <-> 猴中位
  geom_segment(data = cmed,
               aes(x = mou_med, xend = mac_med, y = ynum_m, yend = ynum_q),
               linewidth = 0.5, color = GREY_DUMB, lineend = "round") +
  # 個別 program 點 (高信賴)
  geom_point(data = cc[low_conf == FALSE & is_unres == FALSE],
             aes(cosine, ynum, color = species), size = 0.85, alpha = 0.80) +
  # 個別 program 點 (低信賴 unresolved 弱化)
  geom_point(data = cc[is_unres == TRUE],
             aes(cosine, ynum, color = species), size = 0.85, alpha = 0.40, shape = 17) +
  geom_point(data = cc[low_conf == TRUE & is_unres == FALSE],
             aes(cosine, ynum, color = species), size = 0.85, alpha = 0.42, shape = 17) +
  # 類中位菱形: 人-鼠
  geom_point(data = cmed, aes(mou_med, ynum_m, fill = "mouse"),
             shape = 23, size = 1.9, color = med_cols["mouse"], stroke = 0.4) +
  # 類中位菱形: 人-猴
  geom_point(data = cmed, aes(mac_med, ynum_q, fill = "macaque"),
             shape = 23, size = 1.9, color = med_cols["macaque"], stroke = 0.4) +
  scale_color_manual(values = sp_cols, name = "Species",
                     guide = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  scale_fill_manual(values = sp_cols, guide = "none") +
  scale_y_continuous(breaks = ybase, labels = inner_lv,
                     limits = c(0.45, length(inner_lv) + 0.62), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 0.83), breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                     expand = expansion(mult = c(0, 0.004))) +
  # 軸端極小灰註 (相對更保守), 不搶戲
  annotate("text", x = 0.82, y = length(inner_lv) + 0.42,
           label = "more conserved →", size = gtsize(5), color = "grey45",
           hjust = 1, fontface = "italic") +
  # MWU 統計註 (左下空白區: glia/neuronal/neuropeptide 左側 x<0.13 無點區, 弱化灰, 數字不動)
  annotate("text", x = 0.005, y = 3.4, hjust = 0, vjust = 0.5,
           size = gtsize(5.0), color = GREY_NOTE, lineheight = 1.05,
           label = sprintf("neuron vs non-neuron\nMann–Whitney\nhuman–mouse p=%s ***\nhuman–macaque p=%s (trend)",
                           fmt_p(cst$mou_p[1]), fmt_p(cst$mac_p[1]))) +
  labs(x = "Cross-species alignment cosine", y = NULL) +
  coord_cartesian(clip = "off") +
  zt

ggsave(file.path(OUTD, "c.svg"), pC, width = 2.55, height = 2.32,
       units = "in", device = svglite::svglite)

cat("DONE c.svg\n")
cat(sprintf("aspect W:H = %.3f : 1 (%.2f x %.2f in)\n", 2.55/2.32, 2.55, 2.32))

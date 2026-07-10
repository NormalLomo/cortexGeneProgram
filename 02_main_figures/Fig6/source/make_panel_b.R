#!/usr/bin/env Rscript
# Fig.10 panel b — REDESIGN per ART_SPEC §2.b
# zero-margin vector SVG via svglite
# Aesthetic-only redesign. Science/numbers untouched.
#   - 60 human-K60 programs ranked by human->macaque cosine (h_mac_cosine)
#   - grouped into func_class strips (ART_SPEC §1.3 七大類 aggregation)
#   - per-strip horizontal lanes sharing x = cosine; points = func-class colour
#   - name ONLY top-6 + bottom-6 (P8/P1/P15/P24/P7/P56 high; P18/P52/P35/P4/P58/P48 low)
#   - DISPLAY NAMES: single authority = results/crossregion_v1/program_names.tsv name_short
#   - STARS: 24 brain-weak whitelist (program_names.tsv brain_term == "brain-weak"), NOT per-FDR

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(readr)
  library(svglite)
  library(stringr)
})

base_dir <- "CORTEX_PROGRAM_ROOT/results/xspecies_humanmap_v1/figures/Fig10_v2"
data_dir <- file.path(base_dir, "data")
out_dir  <- file.path(base_dir, "svg_panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- single naming + renumber authority ----
# [2026-06-20 restore] 還原提交版排版 + 套 54-program 新編號。
# 兩張權威表 (crossregion_v1):
#   program_renumber_map.tsv : old_P -> new_P (status; excluded 6 = cNMF 9/18/19/35/52/57)
#   program_names.tsv        : new_P -> name_short / confidence(brain-weak whitelist)
# 各 data csv 的 program 欄 = 舊 P 號(=cNMF號)。策略: 按舊號 join, 顯示一律 new_P + 新名。
renum_authority <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1/program_renumber_map.tsv"
names_authority <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1/program_names.tsv"
renum <- read_tsv(renum_authority, show_col_types = FALSE)
nm    <- read_tsv(names_authority, show_col_types = FALSE)

# old P-id -> new P-id (kept rows only; excluded -> NA, will be filtered by the 54-row csv)
old2new <- renum %>%
  filter(status == "kept") %>%
  transmute(program = paste0("P", old_P), new_id = paste0("P", new_P))

# new P-id -> name_short + brain-weak flag (program_names.tsv keyed by new_P)
new2name <- nm %>%
  transmute(new_id     = new_P,
            name_short = name_short,
            is_weak    = (confidence == "brain-weak"))   # 24 brain-weak whitelist

# combined: old P-id -> new_id + name_short + is_weak
auth_map <- old2new %>% left_join(new2name, by = "new_id")

# ---------------------------------------------------------------- load
lol  <- read_csv(file.path(data_dir, "panelB_lollipop.csv"), show_col_types = FALSE)
# func_class_curated.csv used ONLY for the seven-class grouping enum (function_class),
# NOT for any display name (its name columns are a forbidden stale source).
curc <- read_csv(file.path(data_dir, "func_class_curated.csv"), show_col_types = FALSE)

df <- lol %>%
  select(program, h_mac_cosine, order) %>%
  left_join(curc %>% select(program, function_class), by = "program") %>%
  left_join(auth_map, by = "program")

# ---- aggregate fine func_class -> ART_SPEC §1.3 seven big classes ----
agg_map <- c(
  neuron_synapse         = "Synaptic",
  neuron_ion_channel     = "Synaptic",
  neuron_axon_guidance   = "Synaptic",
  neuron_neuropeptide    = "Neuropeptide/IEG",
  neuron_activity_IEG    = "Neuropeptide/IEG",
  glia_astrocyte         = "Glia (oligo/astro)",
  glia_oligo_myelin      = "Glia (oligo/astro)",
  glia_microglia_immune  = "Microglia/immune",
  vascular               = "Vascular",
  cytoskeleton           = "Cytoskeleton",
  metabolic_housekeeping = "Cytoskeleton",   # housekeeping/structural support -> grouped with cytoskeleton class slot
  other_unresolved       = "Unresolved"
)
df$big_class <- agg_map[df$function_class]
df$big_class[is.na(df$big_class)] <- "Unresolved"

# ART_SPEC §1.3 palette (seven big classes)
class_cols <- c(
  "Synaptic"           = "#3B6FB6",  # 靛藍
  "Neuropeptide/IEG"   = "#6B4E9E",  # 紫
  "Glia (oligo/astro)" = "#3E8E7E",  # 青瓷
  "Vascular"           = "#C8743C",  # 赭橙
  "Microglia/immune"   = "#B0506A",  # 覆盆莓玫
  "Cytoskeleton"       = "#8A7A52",  # 橄欖棕
  "Unresolved"         = "#9AA0A6"   # 中灰 (弱類)
)

# strip order: FIXED to the submitted master's row order (top -> bottom).
# [2026-06-20 criticfix] The data-driven median sort flipped the top two strips
# after the 6-program exclusion (Cytoskeleton median edged above Microglia/immune);
# the submitted figure has Microglia/immune on top. Pin the order to the master so
# the restored figure matches; per-class n is still computed from the (54-row) data.
class_order_top <- c("Microglia/immune", "Cytoskeleton", "Glia (oligo/astro)",
                     "Synaptic", "Vascular", "Neuropeptide/IEG", "Unresolved")
class_n <- df %>% group_by(big_class) %>% summarise(n = n(), .groups = "drop")
class_levels <- class_order_top[class_order_top %in% df$big_class]   # top -> bottom
df$big_class <- factor(df$big_class, levels = rev(class_levels))     # ggplot y bottom->top
n_by <- setNames(class_n$n, class_n$big_class)

# strip side labels: "Class (n)"
class_lab <- setNames(sprintf("%s  (n=%d)", class_levels, n_by[class_levels]), class_levels)

# vertical jitter within each strip lane (deterministic)
set.seed(11)
df <- df %>% group_by(big_class) %>%
  mutate(yj = as.numeric(big_class) + (runif(n()) - 0.5) * 0.52) %>%
  ungroup()

# ---- programs to name: top-6 + bottom-6 by cosine ----
# [2026-06-20 fig54甲] 改為資料驅動 (取代寫死的 60-母體名單; 與 panel e 同邏輯):
# 在 54 母體 (panelB_lollipop.csv 已剔 6) 按 h_mac_cosine 動態取 top6/bottom6.
# 舊寫死名單含 P18/P52/P35 (已剔除), 故必改; top6 仍為 P8/P1/P15/P24/P7/P56,
# bottom6 = 54 母體最低 6 個 (P4/P58/P55/P16/P34/P48 一類).
ord_b     <- df$program[order(df$h_mac_cosine)]
named_bot <- as.character(head(ord_b, 6))
named_top <- as.character(tail(ord_b, 6))
named_set <- c(named_top, named_bot)
df$is_named <- df$program %in% named_set
# display label = NEW P-id + authority name_short; star ONLY if brain-weak (24-whitelist)
# [restore] 顯示用 new_id (新 54 編號), 不用 df$program (舊 cNMF 號).
df$nm_lab <- ifelse(
  df$is_named,
  paste0(df$new_id, " ",
         str_trunc(df$name_short, 22, ellipsis = "…"),
         ifelse(df$is_weak, " *", "")),
  NA_character_
)

# class median markers (diamonds) on each strip
class_med_pts <- df %>% group_by(big_class) %>%
  summarise(med = median(h_mac_cosine), .groups = "drop") %>%
  mutate(yy = as.numeric(big_class))

# axis range
xmin <- 0; xmax <- 0.75

# ---------------------------------------------------------------- theme
pt_to_mm <- 0.352777
ax_col   <- "#3A3A3A"
note_col <- "#5A5A5A"

th <- theme_minimal(base_size = 6.5, base_family = "Helvetica") +
  theme(
    plot.margin       = margin(0, 0, 0, 0),
    panel.grid.major.x = element_line(colour = "#E8E8E8", linewidth = 0.25 * pt_to_mm),
    panel.grid.major.y = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line.x       = element_line(colour = ax_col, linewidth = 0.4 * pt_to_mm),
    axis.ticks.x      = element_line(colour = ax_col, linewidth = 0.4 * pt_to_mm),
    axis.ticks.length = unit(1.4, "pt"),
    axis.ticks.y      = element_blank(),
    axis.title.x      = element_text(size = 6.5, colour = ax_col, margin = margin(t = 1.5)),
    axis.title.y      = element_blank(),
    axis.text.x       = element_text(size = 5.5, colour = ax_col),
    axis.text.y       = element_text(size = 5.5, colour = ax_col, hjust = 1),
    legend.position   = "none",
    plot.title        = element_blank(),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.background  = element_rect(fill = "transparent", colour = NA)
  )

ny <- length(class_levels)

p <- ggplot(df, aes(x = h_mac_cosine, y = yj)) +
  # faint strip dividers
  geom_hline(yintercept = seq(1.5, ny - 0.5, by = 1),
             colour = "#ECECEC", linewidth = 0.25 * pt_to_mm) +
  geom_point(aes(colour = big_class), size = 1.8 * pt_to_mm * 2.83,
             alpha = 0.9, stroke = 0) +
  # class median diamonds
  geom_point(data = class_med_pts, aes(x = med, y = yy),
             shape = 23, size = 2.2 * pt_to_mm * 2.83,
             fill = "white", colour = ax_col, stroke = 0.4 * pt_to_mm) +
  # named labels (repel) — push to the LEFT of each point (right side is clipped),
  # right-justified so the label text grows away from the panel edge
  geom_text_repel(
    data = df %>% filter(is_named),
    aes(label = nm_lab, colour = big_class),
    size = 5 * pt_to_mm, family = "Helvetica",
    segment.size = 0.3 * pt_to_mm, segment.colour = "#B5B5B5",
    box.padding = 0.16, point.padding = 0.10,
    min.segment.length = 0, max.overlaps = Inf,
    hjust = 1, nudge_x = -0.012,
    direction = "both", force = 2.2, force_pull = 0.6, seed = 11,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = class_cols) +
  scale_x_continuous(limits = c(xmin, xmax),
                     breaks = c(0, 0.25, 0.5, 0.75),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(
    breaks = seq_len(ny),
    labels = class_lab[levels(df$big_class)],
    limits = c(0.4, ny + 0.6),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "human→macaque cosine") +
  th

# brain-weak star note inside panel (top-left, grey, 5.5pt)
p <- p + annotate("text", x = xmin + 0.005, y = ny + 0.5,
                  label = "* brain-weak program (21/54; weak brain GO:BP enrichment)",
                  hjust = 0, vjust = 1, size = 5.5 * pt_to_mm,
                  colour = note_col, family = "Helvetica")

# ---------------------------------------------------------------- export zero-margin
# target aspect ~1.3:1 (W:H). pick H by class count, W = 1.3*H
H_in <- 2.35
W_in <- H_in * 1.3
svglite(file.path(out_dir, "b.svg"), width = W_in, height = H_in, bg = "transparent")
print(p)
invisible(dev.off())

cat("PANEL_B_DONE\n")
cat(sprintf("n_programs=%d  n_classes=%d\n", nrow(df), ny))
cat(sprintf("class_order_top_to_bottom: %s\n", paste(class_levels, collapse=" > ")))
cat(sprintf("aspect_W:H = %.3f : %.3f  (%.3f)\n", W_in, H_in, W_in/H_in))
cat(sprintf("out=%s\n", file.path(out_dir, "b.svg")))

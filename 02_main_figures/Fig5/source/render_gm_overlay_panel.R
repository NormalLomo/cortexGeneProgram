#!/usr/bin/env Rscript
# =============================================================================
# render_gm_overlay_panel.R
# GM-only confound-resolution panel for microglia × synaptic hub.
#
# Builds two side-by-side panels:
#  (1) per-pair Z bar/lollipop comparing FULL-CHIP vs GM-only across the 9
#      prompt-spec hub pairs + complement-type control rows
#  (2) per-chip Z_i scatter (full vs GM-only) — 44 chip × 9 hub pair, with
#      identity line; points off-line = chips where GM-only changes the answer
#
# Inputs:
#   CORTEX_PROGRAM_ROOT/results/crossregion_v1/markcorr_betweenchip_v1/
#     betweenchip_microglia_syn_GM_only_compare.tsv  (from 03 extractor)
#     betweenchip_progprog_per_chip_Z.tsv            (full-chip per-chip)
#     betweenchip_progprog_per_chip_Z_GM.tsv         (GM-only per-chip)
#
# Output PDFs:
#   CORTEX_PROGRAM_ROOT/figures/figGM_supp/
#     fig_GM_overlay_lollipop.pdf
#     fig_GM_overlay_perchip_scatter.pdf
#     fig_GM_overlay_combined.pdf  (patchwork compose)
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table); library(scales)
})

ROOT  <- "CORTEX_PROGRAM_ROOT"
BTWN  <- file.path(ROOT, "results/crossregion_v1/markcorr_betweenchip_v1")
OUT   <- file.path(ROOT, "figures/figGM_supp")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

cmp  <- fread(file.path(BTWN, "betweenchip_microglia_syn_GM_only_compare.tsv"))
pcf  <- fread(file.path(BTWN, "betweenchip_progprog_per_chip_Z.tsv"))
pcg  <- fread(file.path(BTWN, "betweenchip_progprog_per_chip_Z_GM.tsv"))

# Limit to prompt_9 + complement controls (v2 schema uses A_newP / B_newP)
cmp_hub <- cmp[pair_set == "prompt_9"]
cmp_ctrl <- cmp[pair_set == "complement_ctrl_3"]
cmp_all <- rbind(cmp_hub, cmp_ctrl)
cmp_all[, pair_label := paste0(A_newP, " × ", B_newP)]
cmp_all[, group := ifelse(pair_set == "prompt_9", "Microglial (P36/P45/P54) × syn (hub)",
                          "Microglial complement-type P49 × syn (control)")]

# Long-format for lollipop
long <- melt(cmp_all[, .(pair_label, group, Z_full, Z_GM)],
             id.vars=c("pair_label","group"),
             variable.name="cond", value.name="Z")
long[, cond := factor(cond, levels=c("Z_full","Z_GM"),
                      labels=c("full-chip", "GM-only"))]

theme_b <- theme_minimal(base_size=7) + theme(
  text          = element_text(size=7, colour="grey15"),
  plot.title    = element_text(size=8, face="bold"),
  plot.subtitle = element_text(size=6.5, colour="grey40"),
  axis.title    = element_text(size=7),
  axis.text     = element_text(size=6.5, colour="grey25"),
  legend.title  = element_text(size=6.5),
  legend.text   = element_text(size=6),
  legend.position = "top",
  panel.grid.minor = element_blank()
)

# ===== Panel 1: lollipop full vs GM per pair =====
COL_FULL <- "#888888"; COL_GM <- "#8E44AD"
p1 <- ggplot(long, aes(x=Z, y=reorder(pair_label, Z), colour=cond)) +
  geom_segment(data=cmp_all, aes(x=Z_full, xend=Z_GM,
                                 y=pair_label, yend=pair_label),
               inherit.aes=FALSE, colour="grey75", linewidth=0.3) +
  geom_point(size=2.2) +
  scale_colour_manual(values=c("full-chip"=COL_FULL, "GM-only"=COL_GM),
                      name=NULL) +
  facet_grid(group ~ ., scales="free_y", space="free_y", switch="y") +
  labs(x="Stouffer Z (combined across 44 sections)",
       y=NULL,
       title="microglia–synaptic hub: full-chip vs GM-only g(r=25 µm)",
       subtitle="GM-only mask drops WM bins, isolates within-grey-matter co-organization") +
  theme_b +
  theme(strip.text.y.left = element_text(angle=0, size=6.5, face="bold"),
        strip.placement = "outside")

ggsave(file.path(OUT, "fig_GM_overlay_lollipop.pdf"), p1,
       width=4.5, height=3.5, device=cairo_pdf)

# ===== Panel 2: per-chip Z_i scatter (full vs GM) =====
# Limit to hub-prompt 9 pairs — TSV column program_N == new_P N (per BUG validation
# 2026-06-26): micro = P36/P45/P54, syn = P16/P26/P43.
hub_micro <- c(36, 45, 54); hub_syn <- c(16, 26, 43)
pcf[, A_newp := as.integer(sub("program_", "", A_name))]
pcf[, B_newp := as.integer(sub("program_", "", B_name))]
pcg[, A_newp := as.integer(sub("program_", "", A_name))]
pcg[, B_newp := as.integer(sub("program_", "", B_name))]
pcf_hub <- pcf[A_newp %in% hub_micro & B_newp %in% hub_syn]
pcg_hub <- pcg[A_newp %in% hub_micro & B_newp %in% hub_syn]

mg <- merge(pcf_hub[, .(A_newp, B_newp, chip_id, Z_full = Z_i)],
            pcg_hub[, .(A_newp, B_newp, chip_id, Z_GM = Z_i)],
            by=c("A_newp","B_newp","chip_id"))
mg[, pair_label := paste0("P", A_newp, "×P", B_newp)]
zmax <- max(abs(c(mg$Z_full, mg$Z_GM)), na.rm=TRUE)

p2 <- ggplot(mg, aes(x=Z_full, y=Z_GM)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey60", linewidth=0.3) +
  geom_hline(yintercept=0, colour="grey80", linewidth=0.2) +
  geom_vline(xintercept=0, colour="grey80", linewidth=0.2) +
  geom_point(alpha=0.55, size=1.3, colour="#8E44AD") +
  facet_wrap(~ pair_label, ncol=3) +
  coord_equal(xlim=c(-zmax, zmax), ylim=c(-zmax, zmax)) +
  labs(x="per-chip Z (full-chip)", y="per-chip Z (GM-only)",
       title="per-section Z_i preserved under GM-only mask",
       subtitle="44 chips × 9 hub pairs; identity line = perfect agreement") +
  theme_b +
  theme(strip.text = element_text(size=6, face="bold"),
        panel.grid.major = element_line(linewidth=0.15, colour="grey92"))

ggsave(file.path(OUT, "fig_GM_overlay_perchip_scatter.pdf"), p2,
       width=4.5, height=4.5, device=cairo_pdf)

# ===== Combined panel =====
combined <- p1 + p2 + plot_layout(widths=c(1, 1)) +
  plot_annotation(
    title = "GM-only confound check for the microglia–synaptic hub (v19 P1.1)",
    subtitle = paste("Per-pair Stouffer Z (left) and per-section Z_i (right)",
                     "for the 9 prompt-spec hub pairs under both nulls."),
    theme = theme(plot.title = element_text(size=9, face="bold"),
                  plot.subtitle = element_text(size=7, colour="grey40"))
  )
ggsave(file.path(OUT, "fig_GM_overlay_combined.pdf"), combined,
       width=9, height=4.6, device=cairo_pdf)

cat("[DONE] panels written to", OUT, "\n")

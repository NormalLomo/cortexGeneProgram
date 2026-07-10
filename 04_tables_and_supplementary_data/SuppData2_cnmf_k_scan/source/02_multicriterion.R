#!/usr/bin/env Rscript
# ED figure: Multi-criterion K selection (K=30-200)
# Defends fixed K=60 as a conservative choice in the conclusion-stable plateau.
# 5 panels, CNS-grade, square-ish panels, fonts >=5pt, vector PDF + PNG, 180mm full width.
# ALL data are read-only from existing result files; input files remain unchanged.

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
  library(ggrepel); library(scales)
})

BASE <- "CORTEX_PROGRAM_ROOT"
EXT  <- file.path(BASE, "figures/extended")
RES  <- file.path(BASE, "results")

# ---------------- palette ----------------
col_k60    <- "#D7263D"   # K=60 accent (red)
col_line   <- "#1B3A5C"   # primary dark navy
col_green  <- "#2E8B57"   # plateau green band
col_red    <- "#D7263D"   # warning red
gauss_col  <- "#2166AC"   # Gaussian IC
pois_col   <- "#B2182B"   # Poisson IC
aic_col    <- "#737373"   # AIC grey

# ---------------- shared theme ----------------
base_fs <- 7
th <- theme_classic(base_size = base_fs) +
  theme(
    plot.title    = element_text(size = base_fs + 1, face = "bold", hjust = 0,
                                 margin = margin(b = 2)),
    plot.subtitle = element_text(size = base_fs - 1.2, colour = "grey25",
                                 lineheight = 0.95, margin = margin(b = 3)),
    axis.title    = element_text(size = base_fs),
    axis.text     = element_text(size = base_fs - 1, colour = "black"),
    axis.line     = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks    = element_line(linewidth = 0.35, colour = "black"),
    legend.title  = element_text(size = base_fs - 1),
    legend.text   = element_text(size = base_fs - 1.5),
    legend.key.size = unit(7, "pt"),
    legend.background = element_blank(),
    legend.margin = margin(0,0,0,0),
    plot.tag      = element_text(size = base_fs + 3, face = "bold"),
    plot.margin   = margin(5, 6, 4, 5),
    aspect.ratio  = 0.85
  )

# ================= DATA =================
full <- read.delim(file.path(RES, "highk_mu_n100/extend_full_30_200.tsv"),
                   check.names = FALSE)
# fill K60 silhouette from builtin
k60 <- read.delim(file.path(RES, "gpu_k60_n100/kstats_row_K60_builtin.tsv"),
                  check.names = FALSE)
k60_sil <- as.numeric(k60$silhouette[1])           # 0.5797
full$silhouette[full$K == 60] <- k60_sil

ic <- read.delim(file.path(RES, "highk_mu_n100/extend_ic_full_20_200.tsv"),
                 check.names = FALSE)

conc <- read.delim(file.path(RES,
        "crossregion_v1/k_robustness_ext_30_80/concordance_ext_30_80.tsv"),
        check.names = FALSE)
drv  <- read.delim(file.path(RES,
        "crossregion_v1/k_robustness_ext_30_80/topdriver_recurrence_ext.tsv"),
        check.names = FALSE)
# L6 CT is #1 driver?  parse first subclass token
drv$first_sub <- sub(";.*$", "", drv$top_subclasses)
drv$l6ct_first <- drv$first_sub == "L6 CT"
conc <- merge(conc, drv[, c("K", "l6ct_first")], by = "K", all.x = TRUE)

nat <- read.delim(file.path(EXT, "ed_kselection_native_sil.tsv"),
                  check.names = FALSE)

# ================= PANEL a : silhouette vs K =================
pa <- ggplot(full, aes(K, silhouette)) +
  geom_line(colour = col_line, linewidth = 0.7) +
  geom_point(colour = col_line, size = 1.3) +
  geom_point(data = subset(full, K == 60), aes(K, silhouette),
             colour = col_k60, size = 2.6) +
  annotate("text", x = 60, y = full$silhouette[full$K==60] + 0.035,
           label = "K=60", colour = col_k60, fontface = "bold",
           size = (base_fs-1)/.pt) +
  scale_x_continuous(breaks = c(30,60,100,150,200)) +
  scale_y_continuous(limits = c(0.38, 0.70), breaks = seq(0.4,0.7,0.1)) +
  labs(title = "Stability (silhouette) vs K",
       subtitle = "Monotone decline favours small K - a known artifact\nof stability criteria (Kotliar et al.)",
       x = "Number of programs K", y = "Consensus silhouette") + th

# ================= PANEL b : distinct@0.85 vs K =================
pb <- ggplot(full, aes(K, `distinct_0.85`)) +
  geom_line(colour = col_line, linewidth = 0.7) +
  geom_point(colour = col_line, size = 1.3) +
  geom_point(data = subset(full, K == 60), aes(K, `distinct_0.85`),
             colour = col_k60, size = 2.6) +
  annotate("text", x = 60, y = full$`distinct_0.85`[full$K==60] - 9,
           label = "K=60", colour = col_k60, fontface = "bold",
           size = (base_fs-1)/.pt) +
  scale_x_continuous(breaks = c(30,60,100,150,200)) +
  labs(title = "Distinct programs (cos < 0.85) vs K",
       subtitle = "Fit / complexity rises without a ceiling\n-> cannot select K on its own",
       x = "Number of programs K", y = "N distinct programs") + th

# ================= PANEL c : IC vs K =================
ic_long <- ic %>%
  transmute(K,
            `BIC (Gaussian)` = BIC_gauss,
            `MDL (Gaussian)` = MDL_gauss,
            `ICL (Gaussian)` = ICL_gauss,
            `AIC (Gaussian)` = AIC_gauss,
            `BIC (Poisson)`  = BIC_pois,
            `MDL (Poisson)`  = MDL_pois,
            `ICL (Poisson)`  = ICL_pois) %>%
  pivot_longer(-K, names_to = "crit", values_to = "val") %>%
  mutate(family = ifelse(grepl("Poisson", crit), "Poisson", "Gaussian"),
         metric = sub(" .*", "", crit),
         metric = factor(metric, levels = c("BIC","MDL","ICL","AIC")))
# scale each metric to its own range for joint display (z within metric+family)
ic_long <- ic_long %>% group_by(crit) %>%
  mutate(scaled = (val - min(val)) / (max(val) - min(val))) %>% ungroup()

# finite optima markers
gauss_bic_opt <- ic$K[which.min(ic$BIC_gauss)]   # 80
pois_bic_opt  <- ic$K[which.min(ic$BIC_pois)]    # 30
opt_df <- bind_rows(
  ic_long %>% filter(crit == "BIC (Gaussian)", K == gauss_bic_opt),
  ic_long %>% filter(crit == "BIC (Poisson)",  K == pois_bic_opt)
)

metric_cols <- c(BIC = "#2166AC", MDL = "#762A83", ICL = "#1B7837", AIC = aic_col)
pc <- ggplot(ic_long, aes(K, scaled, colour = metric, linetype = family,
                          group = crit)) +
  geom_line(linewidth = 0.55) +
  geom_point(data = opt_df, aes(K, scaled), colour = "black", size = 2.6,
             inherit.aes = FALSE) +
  geom_point(data = opt_df, aes(K, scaled, colour = metric), size = 1.6,
             inherit.aes = FALSE, show.legend = FALSE) +
  annotate("text", x = 95, y = 0.18,
           label = "Gaussian BIC\noptimum K=80", size = (base_fs-2.2)/.pt,
           colour = gauss_col, hjust = 0, lineheight = 0.82) +
  annotate("text", x = 35, y = 0.30,
           label = "Poisson BIC\noptimum K=30", size = (base_fs-2.2)/.pt,
           colour = pois_col, hjust = 0, lineheight = 0.82) +
  annotate("segment", x = 150, xend = 198, y = 0.86, yend = 0.985,
           arrow = arrow(length = unit(3,"pt")), colour = aic_col,
           linewidth = 0.4) +
  annotate("text", x = 145, y = 0.84, label = "AIC -> edge",
           size = (base_fs-2.2)/.pt, colour = aic_col, hjust = 1) +
  scale_colour_manual(values = metric_cols, name = "Criterion") +
  scale_linetype_manual(values = c(Gaussian = "solid", Poisson = "22"),
                        name = "Likelihood") +
  scale_x_continuous(breaks = c(20,60,100,150,200)) +
  labs(title = "Information criteria vs K (penalised fit)",
       subtitle = "Principled IC have finite optima (30-80),\nrefuting an unbounded climb",
       x = "Number of programs K", y = "IC (min-max scaled per criterion)") +
  th + theme(legend.position = c(0.80, 0.40),
             legend.spacing.y = unit(0.5,"pt"),
             legend.spacing.x = unit(2,"pt"),
             legend.box = "horizontal",
             legend.key.size = unit(5.5,"pt"),
             legend.title = element_text(size = base_fs - 1.5),
             legend.text = element_text(size = base_fs - 2))

# ================= PANEL d : cross-K rank-concordance =================
plat <- subset(conc, K >= 40 & K <= 80)
dmin <- min(conc$spearman_all); dmax <- max(conc$spearman_all)
pd <- ggplot(conc, aes(K, spearman_all)) +
  # plateau green band
  annotate("rect", xmin = 38, xmax = 82, ymin = 0.955, ymax = 0.995,
           fill = col_green, alpha = 0.12) +
  annotate("text", x = 60, y = 0.998, label = "high plateau (rho 0.96-0.99)",
           colour = col_green, size = (base_fs-2)/.pt, fontface = "italic") +
  # red marker on K30 drop
  annotate("rect", xmin = 27, xmax = 33, ymin = 0.910, ymax = 0.930,
           fill = col_red, alpha = 0.12) +
  geom_line(colour = col_line, linewidth = 0.6, linetype = "31") +
  geom_point(aes(fill = l6ct_first), shape = 21, size = 2.6,
             colour = "black", stroke = 0.4) +
  geom_point(data = subset(conc, K == 60), shape = 21, size = 3.6,
             colour = col_k60, stroke = 1.0, aes(fill = l6ct_first)) +
  annotate("text", x = 30, y = 0.916, label = "K30 = 0.921",
           colour = col_red, size = (base_fs-2)/.pt, hjust = -0.08) +
  annotate("text", x = 60, y = 0.940, label = "K=60\n(plateau centre)",
           colour = col_k60, fontface = "bold", size = (base_fs-1.5)/.pt,
           lineheight = 0.85) +
  scale_fill_manual(values = c(`TRUE` = col_green, `FALSE` = "white"),
                    labels = c(`TRUE` = "L6 CT is #1 driver",
                               `FALSE` = "L6 CT not #1"),
                    name = NULL) +
  scale_x_continuous(breaks = c(30,40,50,60,65,70,80), limits = c(26,84)) +
  scale_y_continuous(limits = c(0.905, 1.002)) +
  labs(title = "Cross-K rank-concordance (Spearman rho)",
       subtitle = "K40-80 stable plateau; K30 collapses (0.921)\nand loses the L6 CT #1 driver",
       x = "Number of programs K", y = "Spearman rho vs K=60 reference") +
  th + theme(legend.position = c(0.70, 0.14),
             legend.key.size = unit(6,"pt"),
             legend.text = element_text(size = base_fs - 2))

# ================= PANEL e : caliber reconcile =================
ee <- bind_rows(
  data.frame(K = nat$K, val = nat$stability,
             series = "cNMF native k-selection stability"),
  data.frame(K = full$K[full$K <= 90], val = full$silhouette[full$K <= 90],
             series = "Extended consensus stability (dt=0.5)")
)
cal_cols <- c("cNMF native k-selection stability" = "#E69F00",
              "Extended consensus stability (dt=0.5)" = col_line)
pe <- ggplot(ee, aes(K, val, colour = series)) +
  annotate("rect", xmin = 48, xmax = 82, ymin = 0.745, ymax = 0.785,
           fill = "#E69F00", alpha = 0.10) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  geom_point(data = subset(ee, K == 60), size = 2.6, colour = col_k60) +
  annotate("text", x = 65, y = 0.79, label = "native plateau (K50-80)",
           colour = "#B8860B", size = (base_fs-2)/.pt, fontface = "italic") +
  annotate("text", x = 30, y = 0.83, label = "0.805", colour = "#B8860B",
           size = (base_fs-2)/.pt) +
  annotate("text", x = 30, y = 0.62, label = "0.654", colour = col_line,
           size = (base_fs-2)/.pt) +
  scale_colour_manual(values = cal_cols, name = NULL) +
  scale_x_continuous(breaks = c(30,40,50,60,70,80,90)) +
  labs(title = "Two calibers of stability reconciled",
       subtitle = "Different calibers, same direction:\nboth disfavour large K",
       x = "Number of programs K", y = "Stability score") +
  th + theme(legend.position = c(0.52, 0.20),
             legend.key.size = unit(6,"pt"))

# ================= assemble =================
title_str <- paste0(
  "Multi-criterion K selection: no single criterion is decisive\n",
  "K=60 is a conservative choice in the conclusion-stable plateau")

# summary text block to fill empty bottom-right slot
summ_txt <- paste(
  "Synthesis",
  "",
  "- Stability (a) and IC penalties favour small K;",
  "  fit/complexity (b) favours large K -- they oppose.",
  "- Information criteria (c) have finite optima",
  "  (Gaussian BIC K=80, Poisson BIC K=30): the",
  "  apparent unbounded climb is rejected.",
  "- Biological conclusions are stable across K40-80",
  "  (d: rho 0.96-0.99, L6 CT remains #1 driver);",
  "  only K30 collapses.",
  "- Both stability calibers (e) point the same way.",
  "",
  "=> K=60 sits in the centre of the conclusion-",
  "   stable plateau: a conservative, defensible choice.",
  sep = "\n")
ptxt <- ggplot() +
  annotate("text", x = 0, y = 1, label = summ_txt, hjust = 0, vjust = 1,
           size = (base_fs-1.5)/.pt, lineheight = 1.05, colour = "grey15") +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(6, 6, 6, 8))

design <- "
AABB
CCDD
EEFF
"
# tag only the 5 data panels a-e; text block (F) gets no letter
fig <- pa + pb + pc + pd + pe + ptxt +
  plot_layout(design = design) +
  plot_annotation(
    title = title_str,
    tag_levels = list(c("a","b","c","d","e","")),
    theme = theme(plot.title = element_text(size = base_fs + 2, face = "bold",
                  hjust = 0, lineheight = 1.05, margin = margin(b = 5)))
  ) &
  theme(plot.tag = element_text(size = base_fs + 3, face = "bold"),
        plot.tag.position = c(0.0, 1.0))

# 180mm full width; 3-row layout -> ~ 210mm tall
W <- 180/25.4; H <- 200/25.4
ggsave(file.path(EXT, "ed_kselection_multicrit.pdf"), fig,
       width = W, height = H, units = "in", device = cairo_pdf)
ggsave(file.path(EXT, "ed_kselection_multicrit.png"), fig,
       width = W, height = H, units = "in", dpi = 320, bg = "white")
cat("DONE: ed_kselection_multicrit.{pdf,png}\n")
cat(sprintf("K60 sil=%.4f  GaussBICopt=K%d  PoisBICopt=K%d\n",
            k60_sil, gauss_bic_opt, pois_bic_opt))

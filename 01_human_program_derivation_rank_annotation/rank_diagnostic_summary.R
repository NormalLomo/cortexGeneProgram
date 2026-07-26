PROJECT_ROOT <- Sys.getenv("CORTEX_PROGRAM_ROOT", unset = "/DATA/cortex_nmf_program")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(Sys.getenv("CORTEX_PROGRAM_ROOT", "/DATA/cortex_nmf_program"), mustWork = FALSE)
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(tidyr)
    library(ggrepel)
    library(scales)
})
BASE <- project_root
EXT <- file.path(BASE, "figures/extended")
RES <- file.path(BASE, "results")
col_k60 <- "#D7263D"
col_line <- "#1B3A5C"
col_green <- "#2E8B57"
col_red <- "#D7263D"
gauss_col <- "#2166AC"
pois_col <- "#B2182B"
aic_col <- "#737373"
base_fs <- 7
th <- theme_classic(base_size = base_fs) + theme(plot.title = element_text(size = base_fs + 1, face = "bold", hjust = 0, margin = margin(b = 2)), plot.subtitle = element_text(size = base_fs - 1.2, colour = "grey25", lineheight = 0.95, margin = margin(b = 3)), axis.title = element_text(size = base_fs), axis.text = element_text(size = base_fs - 1, colour = "black"), axis.line = element_line(linewidth = 0.35, colour = "black"), axis.ticks = element_line(linewidth = 0.35, colour = "black"), legend.title = element_text(size = base_fs - 
    1), legend.text = element_text(size = base_fs - 1.5), legend.key.size = unit(7, "pt"), legend.background = element_blank(), legend.margin = margin(0, 0, 0, 0), plot.tag = element_text(size = base_fs + 3, face = "bold"), plot.margin = margin(5, 6, 4, 5), aspect.ratio = 0.85)
full <- read.delim(file.path(RES, "highk_mu_n100/extend_full_30_200.tsv"), check.names = FALSE)
k60 <- read.delim(file.path(RES, "gpu_k60_n100/kstats_row_K60_builtin.tsv"), check.names = FALSE)
k60_sil <- as.numeric(k60$silhouette[1])
full$silhouette[full$K == 60] <- k60_sil
ic <- read.delim(file.path(RES, "highk_mu_n100/extend_ic_full_20_200.tsv"), check.names = FALSE)
conc <- read.delim(file.path(RES, "crossregion_v1/k_robustness_ext_30_80/concordance_ext_30_80.tsv"), check.names = FALSE)
drv <- read.delim(file.path(RES, "crossregion_v1/k_robustness_ext_30_80/topdriver_recurrence_ext.tsv"), check.names = FALSE)
drv$first_sub <- sub(";.*$", "", drv$top_subclasses)
drv$l6ct_first <- drv$first_sub == "L6 CT"
conc <- merge(conc, drv[, c("K", "l6ct_first")], by = "K", all.x = TRUE)
nat <- read.delim(file.path(EXT, "ed_kselection_native_sil.tsv"), check.names = FALSE)
pa <- ggplot(full, aes(K, silhouette)) + geom_line(colour = col_line, linewidth = 0.7) + geom_point(colour = col_line, size = 1.3) + geom_point(data = subset(full, K == 60), aes(K, silhouette), colour = col_k60, size = 2.6) + annotate("text", x = 60, y = full$silhouette[full$K == 60] + 0.035, label = "K=60", colour = col_k60, fontface = "bold", size = (base_fs - 1)/.pt) + scale_x_continuous(breaks = c(30, 60, 100, 150, 200)) + scale_y_continuous(limits = c(0.38, 0.7), breaks = seq(0.4, 0.7, 0.1)) + 
    labs(title = "Stability (silhouette) vs K", subtitle = "", x = "Number of programs K", y = "Consensus silhouette") + th
pb <- ggplot(full, aes(K, distinct_0.85)) + geom_line(colour = col_line, linewidth = 0.7) + geom_point(colour = col_line, size = 1.3) + geom_point(data = subset(full, K == 60), aes(K, distinct_0.85), colour = col_k60, size = 2.6) + annotate("text", x = 60, y = full$distinct_0.85[full$K == 60] - 9, label = "K=60", colour = col_k60, fontface = "bold", size = (base_fs - 1)/.pt) + scale_x_continuous(breaks = c(30, 60, 100, 150, 200)) + labs(title = "Distinct programs (cos < 0.85) vs K", subtitle = "Fit and complexity increase across the sampled ranks", 
    x = "Number of programs K", y = "N distinct programs") + th
ic_long <- ic %>% transmute(K, `BIC (Gaussian)` = BIC_gauss, `MDL (Gaussian)` = MDL_gauss, `ICL (Gaussian)` = ICL_gauss, `AIC (Gaussian)` = AIC_gauss, `BIC (Poisson)` = BIC_pois, `MDL (Poisson)` = MDL_pois, `ICL (Poisson)` = ICL_pois) %>% pivot_longer(-K, names_to = "crit", values_to = "val") %>% mutate(family = ifelse(grepl("Poisson", crit), "Poisson", "Gaussian"), metric = sub(" .*", "", crit), metric = factor(metric, levels = c("BIC", "MDL", "ICL", "AIC")))
ic_long <- ic_long %>% group_by(crit) %>% mutate(scaled = (val - min(val))/(max(val) - min(val))) %>% ungroup()
gauss_bic_opt <- ic$K[which.min(ic$BIC_gauss)]
pois_bic_opt <- ic$K[which.min(ic$BIC_pois)]
opt_df <- bind_rows(ic_long %>% filter(crit == "BIC (Gaussian)", K == gauss_bic_opt), ic_long %>% filter(crit == "BIC (Poisson)", K == pois_bic_opt))
metric_cols <- c(BIC = "#2166AC", MDL = "#762A83", ICL = "#1B7837", AIC = aic_col)
pc <- ggplot(ic_long, aes(K, scaled, colour = metric, linetype = family, group = crit)) + geom_line(linewidth = 0.55) + geom_point(data = opt_df, aes(K, scaled), colour = "black", size = 2.6, inherit.aes = FALSE) + geom_point(data = opt_df, aes(K, scaled, colour = metric), size = 1.6, inherit.aes = FALSE, show.legend = FALSE) + annotate("text", x = 95, y = 0.18, label = "Gaussian BIC\noptimum K=80", size = (base_fs - 2.2)/.pt, colour = gauss_col, hjust = 0, lineheight = 0.82) + annotate("text", x = 35, 
    y = 0.3, label = "Poisson BIC\noptimum K=30", size = (base_fs - 2.2)/.pt, colour = pois_col, hjust = 0, lineheight = 0.82) + annotate("segment", x = 150, xend = 198, y = 0.86, yend = 0.985, arrow = arrow(length = unit(3, "pt")), colour = aic_col, linewidth = 0.4) + annotate("text", x = 145, y = 0.84, label = "AIC -> edge", size = (base_fs - 2.2)/.pt, colour = aic_col, hjust = 1) + scale_colour_manual(values = metric_cols, name = "Criterion") + scale_linetype_manual(values = c(Gaussian = "solid", 
    Poisson = "22"), name = "Likelihood") + scale_x_continuous(breaks = c(20, 60, 100, 150, 200)) + labs(title = "Information criteria vs K (penalised fit)", subtitle = "IC curves across K", x = "Number of programs K", y = "IC (min-max scaled per criterion)") + th + theme(legend.position = c(0.8, 0.4), legend.spacing.y = unit(0.5, "pt"), legend.spacing.x = unit(2, "pt"), legend.box = "horizontal", legend.key.size = unit(5.5, "pt"), legend.title = element_text(size = base_fs - 1.5), legend.text = element_text(size = base_fs - 
    2))
plat <- subset(conc, K >= 40 & K <= 80)
dmin <- min(conc$spearman_all)
dmax <- max(conc$spearman_all)
pd <- ggplot(conc, aes(K, spearman_all)) + annotate("rect", xmin = 38, xmax = 82, ymin = 0.955, ymax = 0.995, fill = col_green, alpha = 0.12) + annotate("text", x = 60, y = 0.998, label = "high plateau (rho 0.96-0.99)", colour = col_green, size = (base_fs - 2)/.pt, fontface = "italic") + annotate("rect", xmin = 27, xmax = 33, ymin = 0.91, ymax = 0.93, fill = col_red, alpha = 0.12) + geom_line(colour = col_line, linewidth = 0.6, linetype = "31") + geom_point(aes(fill = l6ct_first), shape = 21, size = 2.6, 
    colour = "black", stroke = 0.4) + geom_point(data = subset(conc, K == 60), shape = 21, size = 3.6, colour = col_k60, stroke = 1, aes(fill = l6ct_first)) + annotate("text", x = 30, y = 0.916, label = "K30 = 0.921", colour = col_red, size = (base_fs - 2)/.pt, hjust = -0.08) + annotate("text", x = 60, y = 0.94, label = "K=60\n(plateau centre)", colour = col_k60, fontface = "bold", size = (base_fs - 1.5)/.pt, lineheight = 0.85) + scale_fill_manual(values = c(`TRUE` = col_green, `FALSE` = "white"), 
    labels = c(`TRUE` = "L6 CT first", `FALSE` = "L6 CT other"), name = NULL) + scale_x_continuous(breaks = c(30, 40, 50, 60, 65, 70, 80), limits = c(26, 84)) + scale_y_continuous(limits = c(0.905, 1.002)) + labs(title = "Cross-K rank-concordance (Spearman rho)", subtitle = "K40-80 plateau", x = "Number of programs K", y = "Spearman rho vs K=60 reference") + th + theme(legend.position = c(0.7, 0.14), legend.key.size = unit(6, "pt"), legend.text = element_text(size = base_fs - 2))
ee <- bind_rows(data.frame(K = nat$K, val = nat$stability, series = "cNMF native k-selection stability"), data.frame(K = full$K[full$K <= 90], val = full$silhouette[full$K <= 90], series = "Extended consensus stability (dt=0.5)"))
cal_cols <- c(`cNMF native k-selection stability` = "#E69F00", `Extended consensus stability (dt=0.5)` = col_line)
pe <- ggplot(ee, aes(K, val, colour = series)) + annotate("rect", xmin = 48, xmax = 82, ymin = 0.745, ymax = 0.785, fill = "#E69F00", alpha = 0.1) + geom_line(linewidth = 0.7) + geom_point(size = 1.4) + geom_point(data = subset(ee, K == 60), size = 2.6, colour = col_k60) + annotate("text", x = 65, y = 0.79, label = "native plateau (K50-80)", colour = "#B8860B", size = (base_fs - 2)/.pt, fontface = "italic") + annotate("text", x = 30, y = 0.83, label = "0.805", colour = "#B8860B", size = (base_fs - 
    2)/.pt) + annotate("text", x = 30, y = 0.62, label = "0.654", colour = col_line, size = (base_fs - 2)/.pt) + scale_colour_manual(values = cal_cols, name = NULL) + scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80, 90)) + labs(title = "Stability scores", subtitle = "Native and extended stability", x = "Number of programs K", y = "Stability score") + th + theme(legend.position = c(0.52, 0.2), legend.key.size = unit(6, "pt"))
title_str <- paste0("K selection diagnostics\n", "K=60 in the stable plateau")
summ_txt <- ""
ptxt <- ggplot() + annotate("text", x = 0, y = 1, label = summ_txt, hjust = 0, vjust = 1, size = (base_fs - 1.5)/.pt, lineheight = 1.05, colour = "grey15") + xlim(0, 1) + ylim(0, 1) + theme_void() + theme(plot.margin = margin(6, 6, 6, 8))
design <- "\nAABB\nCCDD\nEEFF\n"
fig <- pa + pb + pc + pd + pe + ptxt + plot_layout(design = design) + plot_annotation(title = title_str, tag_levels = list(c("a", "b", "c", "d", "e", "")), theme = theme(plot.title = element_text(size = base_fs + 2, face = "bold", hjust = 0, lineheight = 1.05, margin = margin(b = 5)))) & theme(plot.tag = element_text(size = base_fs + 3, face = "bold"), plot.tag.position = c(0, 1))
W <- 180/25.4
H <- 200/25.4
ggsave(file.path(EXT, "ed_kselection_multicrit.pdf"), fig, width = W, height = H, units = "in", device = cairo_pdf)
ggsave(file.path(EXT, "ed_kselection_multicrit.png"), fig, width = W, height = H, units = "in", dpi = 320, bg = "white")

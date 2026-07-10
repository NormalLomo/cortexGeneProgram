# Panel B: per-subclass eta2 distribution ridgeline (over the 60 programs)
source("CORTEX_PROGRAM_ROOT/figures/fig4/fig4_theme.R")
suppressMessages({library(ggridges)})

eta <- read.delim(file.path(RES, "within_subclass_region_eta2.tsv"))
drv <- read.delim(file.path(RES, "subclass_driver_rank.tsv"))
ord <- drv %>% arrange(median_eta2) %>% pull(subclass)   # low->high so high on top
eta$subclass <- factor(eta$subclass, levels = ord)
eta$Class <- sc_class(as.character(eta$subclass))
med <- drv %>% transmute(subclass, median_eta2)

p <- ggplot(eta, aes(x = eta2, y = subclass, fill = Class)) +
  geom_density_ridges(scale = 2.2, alpha = 0.78, linewidth = 0.25,
                      colour = "white", quantile_lines = TRUE, quantiles = 2,
                      vline_colour = "grey25", vline_size = 0.3,
                      jittered_points = TRUE, point_size = 0.25,
                      point_alpha = 0.25, position = position_raincloud(height = 0)) +
  scale_fill_manual(values = pal_class, name = "Class") +
  scale_x_continuous(expand = expansion(c(0.01, 0.05))) +
  labs(x = expression(paste("Within-subclass region ", eta^2)),
       y = NULL, title = "Region-driver ranking of cell subclasses",
       subtitle = "Distribution over 60 gene programs; subclasses ordered by median") +
  theme_fig4(8) +
  theme(legend.position = c(0.82, 0.25),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        panel.grid.major.y = element_blank())

ggsave(file.path(FIG, "fig4_b.pdf"), p, width = 4.2, height = 4.6, device = cairo_pdf)
cat("panel b done\n")

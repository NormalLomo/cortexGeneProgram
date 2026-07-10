# Panel C: compositional vs cell-autonomous variance partition (stacked horizontal bar)
source("CORTEX_PROGRAM_ROOT/figures/fig4/fig4_theme.R")

pc <- read.delim(file.path(RES, "panel_c_partition.tsv"))
pc$program <- factor(prog_label(pc$program),
                     levels = prog_label(pc$program[order(pc$eta2_region)]))
long <- pc %>%
  select(program, eta2_region, cell_autonomous_frac, compositional_frac) %>%
  pivot_longer(c(cell_autonomous_frac, compositional_frac),
               names_to = "component", values_to = "frac") %>%
  mutate(component = recode(component,
           cell_autonomous_frac = "Cell-autonomous",
           compositional_frac = "Compositional"),
         contrib = frac * eta2_region)   # scale by total region eta2

pal_comp <- c(`Cell-autonomous` = "#8E44AD", Compositional = "#F39C12")

p <- ggplot(long, aes(x = contrib, y = program, fill = component)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = pal_comp, name = NULL) +
  scale_x_continuous(expand = expansion(c(0, 0.04))) +
  labs(x = expression(paste("Region ", eta^2, "  (partitioned)")), y = NULL,
       title = "Sources of cross-region program variance",
       subtitle = "Top 15 region-variable programs") +
  theme_fig4(8) +
  theme(legend.position = "top",
        legend.justification = "left",
        panel.grid.major.y = element_blank())

ggsave(file.path(FIG, "fig4_c.pdf"), p, width = 3.8, height = 4.4, device = cairo_pdf)
cat("panel c done\n")

# Panel D: spotlight L3-L4 IT RORB - program x region mean activity heatmap, p14 highlighted
source("CORTEX_PROGRAM_ROOT/figures/fig4/fig4_theme.R")

m <- read.delim(file.path(RES, "region_subclass_program_mean.tsv"))
SC <- "L3-L4 IT RORB"
sub <- m %>% filter(subclass == SC)
sub$program <- as.integer(sub$program)

# pick programs that vary across regions for this subclass: top 20 by region range
rng <- sub %>% group_by(program) %>%
  summarise(rng = max(mean) - min(mean), .groups = "drop") %>%
  arrange(desc(rng)) %>% head(20)
sub2 <- sub %>% filter(program %in% rng$program)

# z-score per program across regions for visualization
sub2 <- sub2 %>% group_by(program) %>%
  mutate(z = (mean - mean(mean)) / (sd(mean) + 1e-9)) %>% ungroup()

# order programs by clustering on z (simple: by region of peak then range)
prog_ord <- rng$program
sub2$program_f <- factor(prog_label(sub2$program),
                         levels = prog_label(rev(prog_ord)))
p14_lab <- prog_label(14)   # full functional label for the spotlight program
# order regions by p14 activity gradient
reg_ord <- sub %>% filter(program == 14) %>% arrange(mean) %>% pull(region)
sub2$region <- factor(sub2$region, levels = reg_ord)

p <- ggplot(sub2, aes(x = region, y = program_f, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, name = "z-score") +
  # highlight p14 row
  annotate("rect", xmin = 0.5, xmax = length(reg_ord) + 0.5,
           ymin = which(levels(sub2$program_f) == p14_lab) - 0.5,
           ymax = which(levels(sub2$program_f) == p14_lab) + 0.5,
           fill = NA, colour = "#111111", linewidth = 0.8) +
  labs(x = paste0("Region (ordered by ", p14_lab, ")"), y = NULL,
       title = paste0("Spotlight: ", SC),
       subtitle = paste0("Program activity across 14 regions; ", p14_lab, " highlighted")) +
  theme_fig4(8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        panel.grid = element_blank())

ggsave(file.path(FIG, "fig4_d.pdf"), p, width = 3.9, height = 4.4, device = cairo_pdf)
cat("panel d done\n")

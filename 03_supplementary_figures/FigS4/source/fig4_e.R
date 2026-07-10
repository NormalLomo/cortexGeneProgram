# Panel E: top driver (subclass, program) pairs by eta2 - lollipop+bubble
source("CORTEX_PROGRAM_ROOT/figures/fig4/fig4_theme.R")
suppressMessages(library(ggrepel))

eta <- read.delim(file.path(RES, "within_subclass_region_eta2.tsv"))
eta$Class <- sc_class(eta$subclass)
top <- eta %>% filter(fdr < 0.05) %>% arrange(desc(eta2)) %>% head(25)
top$lab <- paste0(top$subclass, " · ", prog_label(top$program))
top$lab <- factor(top$lab, levels = rev(top$lab))

p <- ggplot(top, aes(x = eta2, y = lab, colour = Class)) +
  geom_segment(aes(x = 0, xend = eta2, yend = lab), linewidth = 0.4,
               colour = "grey75") +
  geom_point(aes(size = n_cells)) +
  scale_colour_manual(values = pal_class, name = "Class") +
  scale_size_continuous(range = c(1.2, 5.5), name = "n cells",
                        breaks = c(10000, 50000, 100000),
                        labels = c("10k", "50k", "100k")) +
  scale_x_continuous(expand = expansion(c(0, 0.06))) +
  labs(x = expression(paste("Within-subclass region ", eta^2)), y = NULL,
       title = "Top region-driver subclass × program pairs",
       subtitle = "Ranked by effect size (FDR < 0.05)") +
  theme_fig4(8) +
  theme(axis.text.y = element_text(size = 6.2),
        legend.position = "right",
        legend.box = "vertical",
        panel.grid.major.y = element_blank())

ggsave(file.path(FIG, "fig4_e.pdf"), p, width = 4.3, height = 4.6, device = cairo_pdf)
cat("panel e done\n")

# Panel G: eta2 robustness via cell-bootstrap (top driver subclasses)
source("CORTEX_PROGRAM_ROOT/figures/fig4/fig4_theme.R")

b <- read.delim(file.path(RES, "panel_g_bootstrap.tsv"))
b$Class <- sc_class(b$subclass)
# focus: show eta2 boxplots per (subclass) top program, faceted by subclass
# order subclasses by driver rank
drv <- read.delim(file.path(RES, "subclass_driver_rank.tsv"))
sc_levels <- intersect(drv$subclass, unique(b$subclass))
b$subclass <- factor(b$subclass, levels = sc_levels)

# within each subclass, order programs by full_eta2 (descending)
prog_order <- b %>% distinct(subclass, program, full_eta2) %>%
  arrange(subclass, desc(full_eta2)) %>%
  mutate(pl = prog_label(program))
b <- b %>% mutate(pl = prog_label(program))

# build per-subclass factor ordering via interaction key
b$key <- paste(b$subclass, b$pl, sep = "::")
key_levels <- prog_order %>% mutate(key = paste(subclass, pl, sep = "::")) %>% pull(key)
b$key <- factor(b$key, levels = rev(key_levels))

p <- ggplot(b, aes(x = eta2, y = key, fill = Class)) +
  geom_boxplot(outlier.size = 0.25, linewidth = 0.25, width = 0.65,
               outlier.alpha = 0.4) +
  geom_point(aes(x = full_eta2), shape = 23, size = 0.9, fill = "white",
             colour = "black", stroke = 0.3) +
  facet_wrap(~subclass, scales = "free", ncol = 4) +
  scale_fill_manual(values = pal_class, name = "Class") +
  scale_y_discrete(labels = function(x) sub(".*::", "", x)) +
  labs(x = expression(paste("Bootstrap ", eta^2, " (20 resamples)")), y = NULL,
       title = "Region-driver effect size is robust to cell resampling",
       subtitle = "Top 10 programs per driver subclass; diamond = full-data estimate (K=60)") +
  theme_fig4(7.5) +
  theme(axis.text.y = element_text(size = 5),
        strip.text = element_text(size = 7, face = "bold"),
        legend.position = "top", legend.justification = "left",
        panel.spacing = unit(2, "mm"))

ggsave(file.path(FIG, "fig4_g.pdf"), p, width = 8.0, height = 4.6, device = cairo_pdf)
cat("panel g done\n")

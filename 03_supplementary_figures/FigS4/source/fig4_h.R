# Panel H: tripartite driver network subclass -> program -> region (ggraph)
source("CORTEX_PROGRAM_ROOT/figures/fig4/fig4_theme.R")
suppressMessages({library(ggraph); library(tidygraph); library(igraph)})

eta <- read.delim(file.path(RES, "within_subclass_region_eta2.tsv"))
eta$program <- as.integer(eta$program)
m   <- read.delim(file.path(RES, "region_subclass_program_mean.tsv"))
m$program <- as.integer(m$program)

# 1) subclass -> program edges: top driver pairs (highest eta2)
sp <- eta %>% filter(fdr < 0.05) %>% arrange(desc(eta2)) %>% head(22)
prog_keep <- sort(unique(sp$program))
sc_keep   <- unique(sp$subclass)

# 2) program -> region edges: for each kept program, the region where it peaks
#    (averaged over the driver subclasses), keep top region per program (+1 runner-up)
pr <- m %>% filter(program %in% prog_keep) %>%
  group_by(program, region) %>% summarise(mean = mean(mean), .groups = "drop") %>%
  group_by(program) %>%
  mutate(spec = (mean - min(mean)) / (max(mean) - min(mean) + 1e-9)) %>%
  arrange(program, desc(mean)) %>% slice_head(n = 2) %>% ungroup()
reg_keep <- unique(pr$region)

# build node table with explicit columns to control layout
nodes <- bind_rows(
  data.frame(name = sc_keep,  type = "Subclass", layer = 1),
  data.frame(name = paste0("p", prog_keep), type = "Program", layer = 2),
  data.frame(name = reg_keep, type = "Region", layer = 3)
) %>% distinct(name, .keep_all = TRUE)
nodes$class <- ifelse(nodes$type == "Subclass", sc_class(nodes$name), nodes$type)
# display label: program nodes get the full functional name, others unchanged
nodes$disp <- nodes$name
.is_prog <- nodes$type == "Program"
nodes$disp[.is_prog] <- prog_label(as.integer(sub("^p", "", nodes$name[.is_prog])))

e1 <- sp %>% transmute(from = subclass, to = paste0("p", program),
                       weight = eta2, etype = "sc_prog")
e2 <- pr %>% transmute(from = paste0("p", program), to = region,
                       weight = spec, etype = "prog_reg")
edges <- bind_rows(e1, e2)

g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

# manual tripartite layout: x by layer, y spread within layer
set.seed(1)
lay <- nodes
lay$x <- lay$layer
# assign y within each layer evenly
lay <- lay %>% group_by(layer) %>%
  mutate(y = scales::rescale(rank(name, ties.method = "first"), to = c(0, 1))) %>%
  ungroup()
layout_xy <- data.frame(x = lay$x, y = lay$y)

pal_node <- c(pal_class, Program = "#7F8C8D", Region = "#34495E")

p <- ggraph(g, layout = "manual", x = layout_xy$x, y = layout_xy$y) +
  geom_edge_diagonal(aes(edge_width = weight, edge_colour = etype,
                         edge_alpha = weight)) +
  scale_edge_width(range = c(0.2, 1.8), guide = "none") +
  scale_edge_alpha(range = c(0.25, 0.9), guide = "none") +
  scale_edge_colour_manual(values = c(sc_prog = "#C0392B", prog_reg = "#2471A3"),
                           labels = c(sc_prog = "subclass→program (η²)",
                                      prog_reg = "program→region (specificity)"),
                           name = NULL) +
  geom_node_point(aes(colour = class, shape = type), size = 2.4) +
  scale_colour_manual(values = pal_node, name = NULL,
                      breaks = c("Excitatory","Inhibitory","Non-neuronal")) +
  scale_shape_manual(values = c(Subclass = 16, Program = 15, Region = 17),
                     name = NULL) +
  geom_node_text(aes(label = disp, hjust = ifelse(layer == 1, 1.12,
                     ifelse(layer == 3, -0.12, 0.5)),
                     vjust = ifelse(layer == 2, -0.9, 0.5)),
                 size = 1.9, colour = "grey15") +
  scale_x_continuous(expand = expansion(c(0.18, 0.18))) +
  scale_y_continuous(expand = expansion(c(0.04, 0.06))) +
  annotate("text", x = 1, y = 1.06, label = "Subclass", size = 2.6, fontface = "bold") +
  annotate("text", x = 2, y = 1.06, label = "Program",  size = 2.6, fontface = "bold") +
  annotate("text", x = 3, y = 1.06, label = "Region",   size = 2.6, fontface = "bold") +
  labs(title = "Cellular-system → program → region driver network",
       subtitle = "Top driver pairs; program→region by peak specificity") +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(size = 6), legend.key.size = unit(3.5, "mm"),
        legend.margin = margin(2, 2, 2, 2), legend.box.spacing = unit(2, "mm"),
        plot.margin = margin(6, 10, 9, 10))

ggsave(file.path(FIG, "fig4_h.pdf"), p, width = 5.2, height = 5.2, device = cairo_pdf)
cat("panel h done\n")

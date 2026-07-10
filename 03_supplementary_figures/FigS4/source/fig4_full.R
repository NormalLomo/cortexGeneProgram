# Assemble fig4_full.pdf : all 8 panels in one vector layout (patchwork + grid)
source("CORTEX_PROGRAM_ROOT/figures/fig4/fig4_theme.R")
suppressMessages({
  library(patchwork); library(grid); library(ggridges); library(ggrepel)
  library(ComplexHeatmap); library(circlize)
  library(ggraph); library(tidygraph); library(igraph)
})

# ---------------- helper builders (return ggplot/grob) ----------------

build_a <- function() {
  eta <- read.delim(file.path(RES, "within_subclass_region_eta2.tsv"))
  eta$program <- as.character(eta$program)
  M <- eta %>% select(subclass, program, eta2) %>%
    pivot_wider(names_from = program, values_from = eta2) %>% as.data.frame()
  rownames(M) <- M$subclass; M$subclass <- NULL
  M <- M[, order(as.integer(colnames(M)))]; M <- as.matrix(M)
  drv <- read.delim(file.path(RES, "subclass_driver_rank.tsv")); M <- M[drv$subclass, ]
  cls <- sc_class(rownames(M))
  col_fun <- colorRamp2(quantile(M, c(0,0.5,0.8,0.95,1)),
                        c("#FCFDBF","#FEC287","#F1605D","#8C2981","#2D1160"))
  ra <- rowAnnotation(Class = cls, col = list(Class = pal_class),
        annotation_name_gp = gpar(fontsize = 6),
        annotation_legend_param = list(Class = list(labels_gp = gpar(fontsize = 6),
            title_gp = gpar(fontsize = 7))), simple_anno_size = unit(2.5, "mm"))
  top <- eta %>% arrange(desc(eta2)) %>% head(12)
  mk <- matrix("", nrow(M), ncol(M), dimnames = dimnames(M))
  for (i in seq_len(nrow(top))) { s<-top$subclass[i]; p<-top$program[i]
    if (s %in% rownames(mk) && p %in% colnames(mk)) mk[s,p] <- "*" }
  hi_progs <- sort(unique(as.integer(top$program)))
  col_labs <- prog_label_selective(colnames(M), hi_progs)
  ht <- Heatmap(M, name="eta2", col=col_fun, cluster_rows=TRUE, cluster_columns=TRUE,
    show_row_dend=TRUE, show_column_dend=TRUE,
    row_dend_width=unit(6,"mm"), column_dend_height=unit(6,"mm"),
    row_names_gp=gpar(fontsize=5.5), column_labels=col_labs, column_names_gp=gpar(fontsize=5),
    column_title="Programs (K=60)", column_title_gp=gpar(fontsize=7),
    row_title="Subclass", row_title_gp=gpar(fontsize=7), left_annotation=ra,
    heatmap_legend_param=list(title=expression(eta^2), labels_gp=gpar(fontsize=6.5),
      title_gp=gpar(fontsize=8), legend_height=unit(20,"mm"), legend_width=unit(5,"mm"),
      grid_width=unit(4,"mm")),
    cell_fun=function(j,i,x,y,w,h,fill){ if(mk[i,j]=="*")
      grid.text("*",x,y,gp=gpar(fontsize=6,col="white",fontface="bold")) },
    rect_gp=gpar(col=NA))
  grid.grabExpr(draw(ht, heatmap_legend_side="right", annotation_legend_side="right",
                     padding=unit(c(1,1,1,1),"mm")))
}

build_b <- function() {
  eta <- read.delim(file.path(RES,"within_subclass_region_eta2.tsv"))
  drv <- read.delim(file.path(RES,"subclass_driver_rank.tsv"))
  ord <- drv %>% arrange(median_eta2) %>% pull(subclass)
  eta$subclass <- factor(eta$subclass, levels=ord); eta$Class <- sc_class(as.character(eta$subclass))
  ggplot(eta, aes(eta2, subclass, fill=Class)) +
    geom_density_ridges(scale=2.0, alpha=0.78, linewidth=0.22, colour="white",
        quantile_lines=TRUE, quantiles=2, vline_colour="grey25", vline_size=0.25) +
    scale_fill_manual(values=pal_class, name="Class") +
    scale_x_continuous(expand=expansion(c(0.01,0.05))) +
    labs(x=expression(paste("Within-subclass region ",eta^2)), y=NULL,
         title="Region-driver ranking", subtitle="over 60 programs") +
    theme_fig4(7) + theme(legend.position=c(0.82,0.22),
      legend.background=element_rect(fill=alpha("white",0.7),colour=NA),
      panel.grid.major.y=element_blank())
}

build_c <- function() {
  pc <- read.delim(file.path(RES,"panel_c_partition.tsv"))
  pc$program <- factor(prog_label(pc$program), levels=prog_label(pc$program[order(pc$eta2_region)]))
  long <- pc %>% select(program,eta2_region,cell_autonomous_frac,compositional_frac) %>%
    pivot_longer(c(cell_autonomous_frac,compositional_frac),names_to="component",values_to="frac") %>%
    mutate(component=recode(component, cell_autonomous_frac="Cell-autonomous",
           compositional_frac="Compositional"), contrib=frac*eta2_region)
  ggplot(long, aes(contrib, program, fill=component)) +
    geom_col(width=0.72, colour="white", linewidth=0.18) +
    scale_fill_manual(values=c(`Cell-autonomous`="#8E44AD",Compositional="#F39C12"),name=NULL)+
    scale_x_continuous(expand=expansion(c(0,0.04))) +
    labs(x=expression(paste("Region ",eta^2," partitioned")), y=NULL,
         title="Variance sources", subtitle="top 15 programs") +
    theme_fig4(7) + theme(legend.position="top", legend.justification="left",
      panel.grid.major.y=element_blank())
}

build_d <- function() {
  m <- read.delim(file.path(RES,"region_subclass_program_mean.tsv"))
  SC <- "L3-L4 IT RORB"; sub <- m %>% filter(subclass==SC); sub$program <- as.integer(sub$program)
  rng <- sub %>% group_by(program) %>% summarise(rng=max(mean)-min(mean),.groups="drop") %>%
    arrange(desc(rng)) %>% head(20)
  sub2 <- sub %>% filter(program %in% rng$program) %>% group_by(program) %>%
    mutate(z=(mean-mean(mean))/(sd(mean)+1e-9)) %>% ungroup()
  sub2$program_f <- factor(prog_label(sub2$program), levels=prog_label(rev(rng$program)))
  p14_lab <- prog_label(14)
  reg_ord <- sub %>% filter(program==14) %>% arrange(mean) %>% pull(region)
  sub2$region <- factor(sub2$region, levels=reg_ord)
  ggplot(sub2, aes(region, program_f, fill=z)) +
    geom_tile(colour="white", linewidth=0.25) +
    scale_fill_gradient2(low="#2166AC",mid="#F7F7F7",high="#B2182B",midpoint=0,name="z")+
    annotate("rect", xmin=0.5, xmax=length(reg_ord)+0.5,
      ymin=which(levels(sub2$program_f)==p14_lab)-0.5, ymax=which(levels(sub2$program_f)==p14_lab)+0.5,
      fill=NA, colour="#111111", linewidth=0.7) +
    labs(x=paste0("Region (by ",p14_lab,")"), y=NULL, title="Spotlight L3-L4 IT RORB", subtitle=paste0(p14_lab," highlighted"))+
    theme_fig4(7) + theme(axis.text.x=element_text(angle=45,hjust=1,size=5),
      axis.text.y=element_text(size=5), panel.grid=element_blank())
}

build_e <- function() {
  eta <- read.delim(file.path(RES,"within_subclass_region_eta2.tsv")); eta$Class <- sc_class(eta$subclass)
  top <- eta %>% filter(fdr<0.05) %>% arrange(desc(eta2)) %>% head(25)
  top$lab <- factor(paste0(top$subclass," · ",prog_label(top$program)),
                    levels=rev(paste0(top$subclass," · ",prog_label(top$program))))
  ggplot(top, aes(eta2, lab, colour=Class)) +
    geom_segment(aes(x=0,xend=eta2,yend=lab), linewidth=0.35, colour="grey75") +
    geom_point(aes(size=n_cells)) +
    scale_colour_manual(values=pal_class,name="Class") +
    scale_size_continuous(range=c(1,4.5),name="n cells",breaks=c(10000,50000,100000),
      labels=c("10k","50k","100k")) + scale_x_continuous(expand=expansion(c(0,0.06))) +
    labs(x=expression(paste("region ",eta^2)), y=NULL, title="Top driver pairs", subtitle="FDR<0.05")+
    theme_fig4(7) + theme(axis.text.y=element_text(size=5), legend.position="right",
      legend.box="vertical", panel.grid.major.y=element_blank())
}

build_f <- function() {  # ggplot twin of matplotlib UMAP (two sub-panels)
  d <- read.delim(file.path(RES,"panel_f_umap.tsv"))
  d <- d[sample(nrow(d)),]
  p_reg <- ggplot(d, aes(UMAP1, UMAP2, colour=region)) +
    geom_point(size=0.18, alpha=0.6) +
    guides(colour=guide_legend(override.aes=list(size=2.2), ncol=2)) +
    labs(title="L3-L4 IT RORB · region", x=NULL, y=NULL) +
    theme_void(base_size=7) + theme(legend.position="right",
      legend.text=element_text(size=5.5), legend.key.size=unit(3,"mm"),
      legend.title=element_blank(), plot.title=element_text(face="bold",size=8,hjust=0))
  p14_lab <- prog_label(14)
  p_p14 <- ggplot(d, aes(UMAP1, UMAP2, colour=p14)) +
    geom_point(size=0.18, alpha=0.7) +
    scale_colour_gradientn(colours=c("#000004","#51127C","#B63679","#FB8861","#FCFDBF"),
      name=p14_lab, limits=quantile(d$p14,c(.02,.98)), oob=scales::squish) +
    labs(title=paste0(p14_lab," activity"), x=NULL, y=NULL) +
    theme_void(base_size=7) + theme(legend.position="right",
      legend.key.height=unit(6,"mm"), legend.key.width=unit(2.4,"mm"),
      legend.text=element_text(size=5.5), legend.title=element_text(size=6.5),
      plot.title=element_text(face="bold",size=8,hjust=0))
  # tighten the inter-panel gap: the region legend (right of p_reg) created a
  # wide floating gap between the two UMAPs. Pull the panels closer and trim
  # the legend box / spacing so the gap closes.
  p_reg <- p_reg + theme(legend.box.spacing = unit(1, "mm"),
                         legend.margin = margin(0, 2, 0, 0),
                         plot.margin = margin(2, 2, 2, 2))
  p_p14 <- p_p14 + theme(legend.box.spacing = unit(1, "mm"),
                         legend.margin = margin(0, 2, 0, 0),
                         plot.margin = margin(2, 2, 2, 2))
  (p_reg + labs(tag="f") + theme(plot.tag=element_text(face="bold",size=12))) +
    p_p14 + plot_layout(widths = c(1, 1)) & theme(plot.margin = margin(2,2,2,2))
}

build_g <- function() {
  b <- read.delim(file.path(RES,"panel_g_bootstrap.tsv")); b$Class <- sc_class(b$subclass)
  drv <- read.delim(file.path(RES,"subclass_driver_rank.tsv"))
  sc_levels <- intersect(drv$subclass, unique(b$subclass)); b$subclass <- factor(b$subclass,levels=sc_levels)
  prog_order <- b %>% distinct(subclass,program,full_eta2) %>% arrange(subclass,desc(full_eta2)) %>%
    mutate(pl=prog_label(program))
  b$pl <- prog_label(b$program); b$key <- paste(b$subclass,b$pl,sep="::")
  key_levels <- prog_order %>% mutate(key=paste(subclass,pl,sep="::")) %>% pull(key)
  b$key <- factor(b$key, levels=rev(key_levels))
  ggplot(b, aes(eta2, key, fill=Class)) +
    geom_boxplot(outlier.size=0.15, linewidth=0.2, width=0.65, outlier.alpha=0.3) +
    geom_point(aes(x=full_eta2), shape=23, size=0.7, fill="white", colour="black", stroke=0.25) +
    facet_wrap(~subclass, scales="free", ncol=4) +
    scale_fill_manual(values=pal_class, name="Class") +
    scale_y_discrete(labels=function(x) sub(".*::","",x)) +
    labs(x=expression(paste("bootstrap ",eta^2)), y=NULL,
         title="Bootstrap robustness of region-driver effect (K=60)",
         subtitle="top 10 programs per subclass; diamond = full-data estimate") +
    theme_fig4(7) + theme(axis.text.y=element_text(size=4), strip.text=element_text(size=6,face="bold"),
      legend.position="top", legend.justification="left", panel.spacing=unit(1.5,"mm"))
}

build_h <- function() {
  eta <- read.delim(file.path(RES,"within_subclass_region_eta2.tsv")); eta$program <- as.integer(eta$program)
  m <- read.delim(file.path(RES,"region_subclass_program_mean.tsv")); m$program <- as.integer(m$program)
  sp <- eta %>% filter(fdr<0.05) %>% arrange(desc(eta2)) %>% head(22)
  prog_keep <- sort(unique(sp$program)); sc_keep <- unique(sp$subclass)
  pr <- m %>% filter(program %in% prog_keep) %>% group_by(program,region) %>%
    summarise(mean=mean(mean),.groups="drop") %>% group_by(program) %>%
    mutate(spec=(mean-min(mean))/(max(mean)-min(mean)+1e-9)) %>%
    arrange(program,desc(mean)) %>% slice_head(n=2) %>% ungroup()
  reg_keep <- unique(pr$region)
  nodes <- bind_rows(
    data.frame(name=sc_keep,type="Subclass",layer=1),
    data.frame(name=paste0("p",prog_keep),type="Program",layer=2),
    data.frame(name=reg_keep,type="Region",layer=3)) %>% distinct(name,.keep_all=TRUE)
  nodes$class <- ifelse(nodes$type=="Subclass", sc_class(nodes$name), nodes$type)
  nodes$disp <- nodes$name
  .is_prog <- nodes$type=="Program"
  nodes$disp[.is_prog] <- prog_label(as.integer(sub("^p","",nodes$name[.is_prog])))
  e1 <- sp %>% transmute(from=subclass,to=paste0("p",program),weight=eta2,etype="sc_prog")
  e2 <- pr %>% transmute(from=paste0("p",program),to=region,weight=spec,etype="prog_reg")
  edges <- bind_rows(e1,e2)
  g <- tbl_graph(nodes=nodes, edges=edges, directed=TRUE)
  lay <- nodes %>% group_by(layer) %>%
    mutate(y=scales::rescale(rank(name,ties.method="first"),to=c(0,1))) %>% ungroup()
  pal_node <- c(pal_class, Program="#7F8C8D", Region="#34495E")
  ggraph(g, layout="manual", x=lay$layer, y=lay$y) +
    geom_edge_diagonal(aes(edge_width=weight, edge_colour=etype, edge_alpha=weight)) +
    scale_edge_width(range=c(0.15,1.4), guide="none") +
    scale_edge_alpha(range=c(0.22,0.85), guide="none") +
    scale_edge_colour_manual(values=c(sc_prog="#C0392B",prog_reg="#2471A3"),
      labels=c(sc_prog="subclass→program",prog_reg="program→region"), name=NULL) +
    geom_node_point(aes(colour=class, shape=type), size=1.8) +
    scale_colour_manual(values=pal_node, name=NULL, breaks=c("Excitatory","Inhibitory","Non-neuronal"))+
    scale_shape_manual(values=c(Subclass=16,Program=15,Region=17), name=NULL) +
    geom_node_text(aes(label=disp, hjust=ifelse(layer==1,1.12,ifelse(layer==3,-0.12,0.5)),
      vjust=ifelse(layer==2,-0.9,0.5)), size=1.5, colour="grey15") +
    scale_x_continuous(expand=expansion(c(0.18,0.18))) +
    scale_y_continuous(expand=expansion(c(0.04,0.08))) +
    annotate("text",x=1,y=1.08,label="Subclass",size=2.2,fontface="bold") +
    annotate("text",x=2,y=1.08,label="Program",size=2.2,fontface="bold") +
    annotate("text",x=3,y=1.08,label="Region",size=2.2,fontface="bold") +
    labs(title="Driver network", subtitle="subclass → program → region") +
    theme_void(base_size=7) + theme(plot.title=element_text(face="bold",size=8,hjust=0),
      plot.subtitle=element_text(size=6.5,colour="grey35"), legend.position="bottom",
      legend.text=element_text(size=5.5), legend.key.size=unit(3,"mm"),
      legend.margin=margin(2,2,2,2), legend.box.spacing=unit(2,"mm"),
      plot.margin=margin(4,8,8,8))
}

# ---------------- assemble ----------------
cat("building panels for assembly...\n")
gA <- build_a(); pB <- build_b(); pC <- build_c(); pD <- build_d()
pE <- build_e(); pF <- build_f(); pG <- build_g(); pH <- build_h()

lab <- function(t) theme(plot.tag=element_text(face="bold",size=11))
tagify <- function(p, t) p + labs(tag=t) + theme(plot.tag=element_text(face="bold",size=12))

wA <- wrap_elements(full = gA) + labs(tag="a") +
      theme(plot.tag=element_text(face="bold",size=12), plot.tag.position=c(0.01,0.99))

# layout: 4 columns conceptually. Use patchwork design.
# Row1: a (wide, spans) | b
# Row2: c | d | e
# Row3: f (wide UMAP pair)
# Row4: g (wide) ; h
design <- "
AAAABB
AAAABB
CCDDEE
CCDDEE
FFFFFF
FFFFFF
GGGGHH
GGGGHH
"
combined <- wA + tagify(pB,"b") + tagify(pC,"c") + tagify(pD,"d") +
            tagify(pE,"e") + pF +
            tagify(pG,"g") + tagify(pH,"h") +
  # Row heights: panel a (rows 1-2) was under-filling vertically (whitespace
  # below the subclass×program heatmap). Shave its two rows so the heatmap
  # content sits at ~1.6:1 with no trailing whitespace; other row-groups keep
  # their (equal) weight.
  plot_layout(design = design,
              heights = c(0.72, 0.72, 1, 1, 1, 1, 1, 1)) +
  plot_annotation(title = "Fig. 4  Cellular system × region driver of cross-region gene-program differences in human cortex",
                  theme = theme(plot.title=element_text(face="bold", size=13)))
ggsave(file.path(FIG,"fig4_full.pdf"), combined, width=15, height=19, device=cairo_pdf, limitsize=FALSE)
ggsave(file.path(FIG,"fig4_full.png"), combined, width=15, height=19, dpi=200, bg="white", limitsize=FALSE)
file.copy(file.path(FIG,"fig4_full.png"), "/tmp/fig4_aspect.png", overwrite = TRUE)
cat("fig4_full done\n")

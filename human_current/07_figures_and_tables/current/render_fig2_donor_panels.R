if (!nzchar(Sys.getenv("NMF_WORK_ROOT"))) stop("Set NMF_WORK_ROOT; original source directories must remain read-only.")
#!/usr/bin/env Rscript
options(future.globals.maxSize = 50 * 1024^3)

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

root <- paste0(Sys.getenv("NMF_WORK_ROOT"), "")
input_root <- file.path(root, "inputs", "cortex_nmf_program")
out_dir <- file.path(root, "figures", "human_revision", "panels", "fig2_donor_evidence")
# The dedicated subclass branch never enters the retained A/D producer.
if ("--within-subclass-only" %in% commandArgs(trailingOnly = TRUE)) {
  subclass_dir <- file.path(input_root, "archived", "results", "v60_formal21_repro",
                            "results", "crossregion_v2")
  sources <- lapply(file.path(subclass_dir, c("subclass_nature_lmm.tsv",
                    "subclass_limma_dupcor.tsv", "subclass_dream.tsv")), function(path) {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
               quote = "", comment.char = "", na.strings = c("", "NA", "NaN"))
  })
  annotation_c <- read.delim(file.path(root, "tables", "TableS3_program_annotation.tsv"),
                            check.names = FALSE, stringsAsFactors = FALSE,
                            quote = "", comment.char = "")
  programs_c <- as.character(annotation_c$new_P[match(paste0("P", seq_len(54L)), annotation_c$new_P)])
  subclasses_c <- unique(as.character(sources[[1L]]$subclass))
  keys_c <- as.vector(outer(subclasses_c, programs_c, paste, sep = "\r"))
  stored_boolean <- function(x) {
    text <- tolower(trimws(as.character(x)))
    answer <- rep(NA, length(text))
    answer[!is.na(text) & text %in% c("true", "1")] <- TRUE
    answer[!is.na(text) & text %in% c("false", "0")] <- FALSE
    answer
  }
  calls_c <- lapply(sources, function(source) {
    index <- match(keys_c, paste(source$subclass, source$program, sep = "\r"))
    rows <- source[index, , drop = FALSE]
    q <- suppressWarnings(as.numeric(rows$BH_global_22x54))
    calls <- stored_boolean(rows$significant_BH_global_22x54)
    status <- tolower(as.character(rows$status))
    available <- !is.na(index) & is.finite(q) & !is.na(calls) &
      !is.na(status) & status %in% c("ok", "ok_singular")
    if ("estimable" %in% names(rows)) {
      estimable <- stored_boolean(rows$estimable)
      available <- available & !is.na(estimable) & estimable
    }
    calls[!available] <- NA
    calls
  })
  # Only count original calls for display; no new q or consensus subset.
  displayed_c <- matrix(Reduce(`+`, lapply(calls_c, as.integer)),
                        nrow = length(subclasses_c), ncol = length(programs_c),
                        dimnames = list(subclasses_c, programs_c))
  fills_c <- c("0" = "#F3F3F3", "1" = "#E7CE8E", "2" = "#55A397", "3" = "#51417C")
  missing_c <- "#A6A6A6"
  render_subclass_panel <- function() {
    W <- 531; H <- 190
    grid.newpage()
    pushViewport(viewport(width = unit(W, "pt"), height = unit(H, "pt")))
    label <- function(text, x, y, size = 6.5, just = "left", bold = FALSE, rot = 0) {
      grid.text(text, x = unit(x, "pt"), y = unit(H-y, "pt"),
                just = c(just, "centre"), rot = rot,
                gp = gpar(fontfamily = "sans", fontsize = size,
                          fontface = if (bold) "bold" else "plain", col = "#18232D"))
    }
    label("f", 3, 10, 13, bold = TRUE)
    label("Donor-aware regional effects within subclasses", 24, 9, 8, bold = TRUE)
    label("LMM / limma duplicateCorrelation / dream; count, not effect size", 24, 20, 6)
    label("Methods with q < 0.05; global BH: 22 x 54", 277, 8, 6.3)
    for (k in 0:3) {
      xx <- 279+k*28
      grid.rect(x = unit(xx+4, "pt"), y = unit(H-20, "pt"),
                width = unit(8, "pt"), height = unit(7, "pt"),
                gp = gpar(fill = fills_c[as.character(k)], col = "#CAD2D8", lwd = 0.4))
      label(as.character(k), xx+11, 20, 6.3)
    }
    grid.rect(x = unit(400, "pt"), y = unit(H-20, "pt"), width = unit(8, "pt"),
              height = unit(7, "pt"), gp = gpar(fill = missing_c, col = NA))
    label("NA: any method unavailable", 408, 20, 6.1)
    left <- 89; right <- 526; top <- 27; row_height <- 6.5
    column_width <- (right-left)/length(programs_c)
    bottom <- top+row_height*length(subclasses_c)
    for (i in seq_along(subclasses_c)) {
      cy <- top+(i-0.5)*row_height
      label(subclasses_c[[i]], left-5, cy, 6.5, "right")
      for (j in seq_along(programs_c)) {
        value <- displayed_c[i,j]
        fill <- if (is.na(value)) missing_c else fills_c[as.character(value)]
        grid.rect(x = unit(left+(j-0.5)*column_width, "pt"), y = unit(H-cy, "pt"),
                  width = unit(column_width, "pt"), height = unit(row_height, "pt"),
                  gp = gpar(fill = fill, col = "#D9E0E5", lwd = 0.22))
      }
    }
    for (j in seq_along(programs_c)) {
      label(programs_c[[j]], left+(j-0.5)*column_width, bottom+4, 6.5, "right", rot = 90)
    }
    upViewport()
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  c_pdf <- file.path(out_dir, "Fig2_within_subclass_regional_support.pdf")
  c_png <- file.path(out_dir, "Fig2_within_subclass_regional_support.png")
  cairo_pdf(c_pdf, width = 531/72, height = 190/72, family = "sans", onefile = FALSE)
  render_subclass_panel()
  invisible(dev.off())
  grDevices::png(c_png, width = round(531/72*300), height = round(190/72*300),
                units = "px", res = 300, type = "cairo", bg = "white")
  render_subclass_panel()
  invisible(dev.off())
  cat("Produced", c_pdf, c_png, "\n")
  quit(save = "no")
}

table_s4_path <- file.path(
  input_root, "archived", "revision_v58_regional_consensus", "data",
  "TableS4_donor_aware_regional_programs.tsv"
)
consensus_path <- file.path(
  input_root, "archived", "results", "v60_formal21_repro", "results",
  "crossregion_v2", "regional_consensus_21.tsv"
)
emm_path <- file.path(
  input_root, "archived", "results", "v60_formal21_repro", "revision_v59_fig2",
  "panel_egh", "source_data", "adjusted_area_emmeans_21x14.tsv"
)
annotation_path <- file.path(root, "tables", "TableS3_program_annotation.tsv")

s4 <- read.delim(table_s4_path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
consensus <- read.delim(consensus_path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
annotation <- read.delim(annotation_path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
emm <- read.delim(emm_path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")

canonical_programs <- paste0("P", seq_len(54L))
consensus_programs <- as.character(consensus$program)
program_order <- c(consensus_programs, canonical_programs[!canonical_programs %in% consensus_programs])

class_map <- c(
  exc = "Excitatory",
  inh = "Inhibitory",
  glia = "Glial",
  vascular = "Vascular",
  vasc = "Vascular",
  nonneuron = "Non-neuronal"
)
class_levels <- c("Excitatory", "Inhibitory", "Glial", "Vascular", "Non-neuronal")
class_colors <- c(
  "Excitatory" = "#3E78B2",
  "Inhibitory" = "#C57938",
  "Glial" = "#4C9A7F",
  "Vascular" = "#B45556",
  "Non-neuronal" = "#8C65B8"
)
confidence_levels <- c("Higher", "Lower")
confidence_colors <- c("Higher" = "#4B5563", "Lower" = "#C2413F")

ann_index <- match(program_order, annotation$new_P)
class_values <- unname(class_map[as.character(annotation$dominant_class[ann_index])])
confidence_values <- ifelse(
  grepl("lower", tolower(as.character(annotation$confidence[ann_index]))),
  "Lower", "Higher"
)
functional_names <- as.character(annotation$functional_name[ann_index])
if (anyNA(functional_names) || any(!nzchar(functional_names))) {
  stop("TableS3 functional_name is incomplete for the displayed programs")
}
display_names <- paste0(
  program_order, " ", functional_names,
  ifelse(confidence_values == "Lower", "*", "")
)
main_names <- display_names[match(consensus_programs, program_order)]
support_columns <- c(
  "nature_lmm_overall_significant_BH_54",
  "limma_dupcor_overall_significant_BH_54",
  "dream_overall_significant_BH_54"
)
method_short <- c("LMM", "limma", "dream")
support_source <- s4[match(program_order, s4$program), support_columns, drop = FALSE]
support_all <- as.matrix(support_source)
support_all <- matrix(
  ifelse(tolower(as.character(support_all)) %in% c("true", "1"), "Yes", "No"),
  nrow = nrow(support_all), ncol = ncol(support_all),
  dimnames = list(program_order, method_short)
)
support21 <- support_all[match(consensus_programs, rownames(support_all)), , drop = FALSE]
other_programs <- canonical_programs[!canonical_programs %in% consensus_programs]
other_names <- display_names[match(other_programs, program_order)]
support_other <- support_all[match(other_programs, rownames(support_all)), , drop = FALSE]
support_colors <- c("No" = "#EEF2F4", "Yes" = "#2F6DA5")

class_values <- class_values[match(program_order, program_order)]
confidence_values <- confidence_values[match(program_order, program_order)]
class_annotation_bottom <- rowAnnotation(
  "Class" = factor(class_values[seq_along(consensus_programs)], levels = class_levels),
  "Conf." = factor(confidence_values[seq_along(consensus_programs)], levels = confidence_levels),
  col = list("Class" = class_colors, "Conf." = confidence_colors),
  simple_anno_size = unit(2.6, "mm"),
  annotation_name_gp = gpar(fontfamily = "sans", fontsize = 4.6),
  show_annotation_name = FALSE,
  show_legend = c(FALSE, FALSE)
)

program_name_annotation <- rowAnnotation(
  Program = anno_text(
    main_names,
    just = "left",
    location = unit(0, "npc"),
    gp = gpar(fontfamily = "sans", fontsize = 6.5)
  ),
  width = unit(180, "pt"),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht_support21 <- Heatmap(
  support21,
  name = "Overall support",
  col = support_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  column_names_gp = gpar(fontfamily = "sans", fontsize = 6.5),
  column_names_rot = 90,
  column_title = NULL,
  column_title_gp = gpar(fontfamily = "sans", fontsize = 6.5),
  rect_gp = gpar(col = "white", lwd = 0.3),
  border = TRUE,
  width = unit(28, "pt"),
  height = unit(170, "pt"),
  show_heatmap_legend = FALSE
)

other_program_labels_1 <- other_programs[seq_len(17L)]
other_program_labels_2 <- other_programs[18:33]
ht_other1 <- Heatmap(
  support_other[seq_len(17L), , drop = FALSE],
  name = "Other support 1",
  col = support_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = gpar(fontfamily = "sans", fontsize = 6.5),
  column_names_gp = gpar(fontfamily = "sans", fontsize = 6.5),
  column_names_rot = 90,
  column_title = NULL,
  column_title_gp = gpar(fontfamily = "sans", fontsize = 6.5),
  rect_gp = gpar(col = "white", lwd = 0.25),
  border = TRUE,
  width = unit(28, "pt"),
  height = unit(110, "pt"),
  show_heatmap_legend = FALSE
)
ht_other2 <- Heatmap(
  support_other[18:33, , drop = FALSE],
  name = "Other support 2",
  col = support_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = gpar(fontfamily = "sans", fontsize = 6.5),
  column_names_gp = gpar(fontfamily = "sans", fontsize = 6.5),
  column_names_rot = 90,
  rect_gp = gpar(col = "white", lwd = 0.25),
  border = TRUE,
  width = unit(28, "pt"),
  height = unit(110, "pt"),
  show_heatmap_legend = FALSE
)

region_meta <- emm[order(as.numeric(emm$area_order)), c("region", "area_order", "cortical_division"), drop = FALSE]
region_meta <- region_meta[!duplicated(region_meta$region), , drop = FALSE]
regions <- as.character(region_meta$region)
region_split <- factor(
  as.character(region_meta$cortical_division),
  levels = unique(as.character(region_meta$cortical_division))
)

profile_matrix <- matrix(
  NA_real_, nrow = length(consensus_programs), ncol = length(regions),
  dimnames = list(consensus_programs, regions)
)
for (i in seq_along(consensus_programs)) {
  rows <- emm[emm$program == consensus_programs[[i]], , drop = FALSE]
  profile_matrix[i, match(rows$region, regions)] <- as.numeric(rows$profile_z)
}

finite_values <- as.numeric(profile_matrix[is.finite(profile_matrix)])
profile_range <- range(finite_values)
profile_ticks <- pretty(profile_range, n = 3)
if (profile_range[1] < 0 && profile_range[2] > 0) {
  profile_col_fun <- colorRamp2(
    c(profile_range[1], 0, profile_range[2]),
    c("#2166AC", "#FFFFFF", "#B2182B")
  )
} else {
  profile_col_fun <- colorRamp2(
    c(profile_range[1], profile_range[2]),
    c("#2166AC", "#B2182B")
  )
}

profile_peak <- matrix(FALSE, nrow = nrow(profile_matrix), ncol = ncol(profile_matrix))
for (i in seq_len(nrow(profile_matrix))) {
  row_values <- profile_matrix[i, ]
  row_max <- max(row_values, na.rm = TRUE)
  profile_peak[i, ] <- is.finite(row_values) & row_values == row_max
}

profile_layer_fun <- function(j, i, x, y, width, height, fill) {
  marked <- vapply(seq_along(j), function(k) profile_peak[i[[k]], j[[k]]], logical(1))
  if (any(marked)) {
    grid.rect(
      x = x[marked], y = y[marked],
      width = width[marked] * 0.84, height = height[marked] * 0.84,
      gp = gpar(col = "#1F2937", fill = NA, lwd = 1.0)
    )
  }
}

ht_profile <- Heatmap(
  profile_matrix,
  name = "profile_z",
  col = profile_col_fun,
  na_col = "#F7F7F7",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = region_split,
  column_gap = unit(1.0, "mm"),
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  left_annotation = class_annotation_bottom,
  row_names_side = "left",
  show_row_names = FALSE,
  row_names_gp = gpar(fontfamily = "sans", fontsize = 6.5),
  column_names_gp = gpar(fontfamily = "sans", fontsize = 6.5),
  column_names_rot = 55,
  column_title = NULL,
  column_title_gp = gpar(fontfamily = "sans", fontsize = 7.0),
  rect_gp = gpar(col = "white", lwd = 0.22),
  border = TRUE,
  width = unit(300, "pt"),
  height = unit(170, "pt"),
  layer_fun = profile_layer_fun,
  heatmap_legend_param = list(
    title = "profile_z",
    at = profile_ticks,
    labels_gp = gpar(fontsize = 6.5),
    title_gp = gpar(fontsize = 6.5)
  ),
  show_heatmap_legend = FALSE
)

profile_legend <- Legend(
  col_fun = profile_col_fun,
  title = "profile_z",
  at = profile_ticks,
  labels = format(profile_ticks, trim = TRUE),
  direction = "horizontal",
  legend_width = unit(65, "pt"),
  title_position = "topcenter",
  labels_gp = gpar(fontfamily = "sans", fontsize = 6.5),
  title_gp = gpar(fontfamily = "sans", fontsize = 6.5)
)
support_legend <- Legend(
  title = "BH FDR < 0.05",
  labels = c("No", "Yes"),
  legend_gp = gpar(fill = unname(support_colors[c("No", "Yes")]), col = NA),
  grid_width = unit(3.5, "mm"),
  grid_height = unit(3.5, "mm"),
  direction = "horizontal",
  nrow = 1,
  title_gp = gpar(fontfamily = "sans", fontsize = 6.5),
  labels_gp = gpar(fontfamily = "sans", fontsize = 6.5)
)
draw_other_support <- function(mat, names, x, y, width, height) {
  n_rows <- nrow(mat)
  label_width <- 113
  cell_x0 <- x + label_width + 2
  cell_width <- (width - label_width - 2) / ncol(mat)
  row_height <- height / n_rows
  for (i in seq_len(n_rows)) {
    yy <- y + height - (i - 0.5) * row_height
    grid.text(
      names[[i]], x = unit(x, "pt"), y = unit(yy, "pt"),
      just = c("left", "center"),
      gp = gpar(fontfamily = "sans", fontsize = 6.5, col = "#111827")
    )
    for (j in seq_len(ncol(mat))) {
      xx <- cell_x0 + (j - 0.5) * cell_width
      grid.rect(
        x = unit(xx, "pt"), y = unit(yy, "pt"),
        width = unit(max(cell_width - 0.6, 1.0), "pt"),
        height = unit(max(row_height - 0.6, 1.0), "pt"),
        gp = gpar(col = "white", lwd = 0.25,
                  fill = unname(support_colors[as.character(mat[i, j])]))
      )
    }
  }
  for (j in seq_len(ncol(mat))) {
    xx <- cell_x0 + (j - 0.5) * cell_width
    grid.text(
      method_short[[j]], x = unit(xx, "pt"), y = unit(y + height + 2, "pt"),
      rot = 90, just = c("left", "center"),
      gp = gpar(fontfamily = "sans", fontsize = 5.5, col = "#111827")
    )
  }
}

render_heatmap_panel <- function() {
  grid.newpage()
  W <- 531; H <- 240
  pushViewport(viewport(width=unit(W,"pt"), height=unit(H,"pt"),
                        gp=gpar(fontfamily="sans",fontsize=6.5), clip="off"))
  text_at <- function(label, x, y, size=6.5, just="left", rot=0, bold=FALSE, col="#1F2937") {
    grid.text(label, x=unit(x,"pt"), y=unit(H-y,"pt"), just=just, rot=rot,
              gp=gpar(fontfamily="sans",fontsize=size,fontface=if(bold) "bold" else "plain",col=col))
  }
  rect_at <- function(x,y,w,h,fill,border=NA,lwd=0.25) {
    grid.rect(x=unit(x,"pt"),y=unit(H-y,"pt"),width=unit(w,"pt"),height=unit(h,"pt"),
              just=c("left","top"),gp=gpar(fill=fill,col=border,lwd=lwd))
  }
  text_width <- function(z) convertWidth(grobWidth(textGrob(z,gp=gpar(fontsize=6.5,fontfamily="sans"))),"pt",valueOnly=TRUE)
  text_at("a",3,10,13,bold=TRUE,col="black")
  text_at("Adjusted regional profiles of 21 consensus programs",21,10,8.5,bold=TRUE)
  name_right <- 3+max(vapply(main_names,text_width,numeric(1)))
  x0 <- name_right+5; top <- 24; body_h <- 145; body_w <- 284-x0
  divisions <- as.character(region_meta$cortical_division)
  gap_before <- c(FALSE,divisions[-1] != divisions[-length(divisions)])
  gap <- 2
  cw <- (body_w-gap*sum(gap_before))/length(regions)
  rh <- body_h/length(consensus_programs)
  col_left <- x0+(seq_along(regions)-1)*cw+cumsum(gap_before)*gap
  for (i in seq_along(consensus_programs)) {
    yy <- top+(i-0.5)*rh
    text_at(main_names[[i]],name_right,yy,just="right")
    for (j in seq_along(regions)) {
      rect_at(col_left[[j]],top+(i-1)*rh,cw,rh,profile_col_fun(profile_matrix[i,j]),"white",0.22)
      if (profile_peak[i,j]) {
        rect_at(col_left[[j]]+cw*0.08,top+(i-1)*rh+rh*0.08,
                cw*0.84,rh*0.84,NA,"#1F2937",0.75)
      }
    }
    for (j in 1:3) rect_at(288+(j-1)*3,top+(i-1)*rh,3,rh,
                            unname(support_colors[support21[i,j]]),"white",0.2)
    rect_at(301,top+(i-1)*rh,3,rh,unname(class_colors[class_values[[i]]]),"white",0.2)
    rect_at(307,top+(i-1)*rh,3,rh,unname(confidence_colors[confidence_values[[i]]]),"white",0.2)
  }
  for (g in unique(divisions)) {
    ix <- which(divisions==g)
    rect_at(col_left[min(ix)],top,col_left[max(ix)]+cw-col_left[min(ix)],body_h,NA,"#343A40",0.45)
  }
  for (j in seq_along(regions)) text_at(regions[[j]],col_left[[j]]+cw/2,top+body_h+5,
                                        just="right",rot=55)
  text_at("Methods",292.5,top+body_h+4,just="right",rot=90)
  text_at("Class",302.5,top+body_h+4,just="right",rot=90)
  text_at("Conf.",308.5,top+body_h+4,just="right",rot=90)

  # Two close-set legend lines, authored in final-size points.
  text_at("profile_z",6,212)
  limits <- range(c(profile_range,profile_ticks))
  for (k in 0:79) rect_at(43+k*0.9,208,0.92,6,
                         profile_col_fun(limits[1]+(k+0.5)/80*diff(limits)))
  for (v in profile_ticks) {
    xx <- 43+(v-limits[1])/diff(limits)*72
    grid.lines(x=unit(c(xx,xx),"pt"),y=unit(H-c(214,216),"pt"),gp=gpar(col="#374151",lwd=0.4))
    text_at(format(v,trim=TRUE),xx,221,just="center")
  }
  text_at("BH FDR < 0.05",130,212)
  rect_at(183,209,5,6,unname(support_colors["No"]))
  text_at("No",191,212)
  rect_at(207,209,5,6,unname(support_colors["Yes"]))
  text_at("Yes",215,212)
  xx <- 240
  text_at("Class:",xx,212); xx <- xx+25
  for (label in class_levels) {
    rect_at(xx,209,4,6,unname(class_colors[label]))
    text_at(label,xx+7,212)
    xx <- xx+7+text_width(label)+7
  }
  text_at("Confidence:",6,232)
  rect_at(47,229,5,6,unname(confidence_colors["Higher"])); text_at("Higher",55,232)
  rect_at(81,229,5,6,unname(confidence_colors["Lower"])); text_at("Lower*",89,232)
  text_at("LMM / limma / dream (left to right)",123,232)
  text_at("limma = duplicateCorrelation",226,232)
  rect_at(327,229,5,6,NA,"#1F2937",0.75)
  text_at("Peak: within-program max; descriptive, not regional significance",338,232)
  upViewport()
}

example_programs <- c("P6", "P3", "P13", "P52", "P16", "P9")
example_rows <- emm[emm$program %in% example_programs, , drop = FALSE]
example_region_meta <- region_meta
example_regions <- regions
example_values <- matrix(
  NA_real_, nrow = length(example_programs), ncol = length(example_regions),
  dimnames = list(example_programs, example_regions)
)
for (i in seq_along(example_programs)) {
  rows <- example_rows[example_rows$program == example_programs[[i]], , drop = FALSE]
  example_values[i, match(rows$region, example_regions)] <- as.numeric(rows$profile_z)
}
example_division <- as.character(example_region_meta$cortical_division)
example_division_levels <- unique(example_division)
example_palette <- c(
  Frontal = "#C44E52", Cingulate = "#8064A2",
  "Parietal/somatosensory" = "#33A089",
  Temporal = "#E2A22C", Occipital = "#4C6EB1"
)
example_colors <- unname(example_palette[example_division])
example_names <- as.character(annotation$functional_name[match(example_programs, annotation$new_P)])
example_range <- range(as.numeric(example_values[is.finite(example_values)]))
example_pad <- max(diff(example_range) * 0.06, 0.15)
example_xlim <- c(example_range[1] - example_pad, example_range[2] + example_pad)
example_ticks <- pretty(example_xlim, n = 3)
example_ticks <- example_ticks[example_ticks >= example_xlim[1] & example_ticks <= example_xlim[2]]

render_examples_panel <- function() {
  grid.newpage()
  W <- 531; H <- 148
  pushViewport(viewport(width=unit(W,"pt"),height=unit(H,"pt"),clip="off"))
  text_at <- function(label,x,y,size=6.5,just="center",bold=FALSE,col="#1F2937") {
    grid.text(label,x=unit(x,"pt"),y=unit(H-y,"pt"),just=just,
              gp=gpar(fontfamily="sans",fontsize=size,fontface=if(bold) "bold" else "plain",col=col))
  }
  text_width <- function(label,size=8.5,bold=FALSE) {
    convertWidth(grobWidth(textGrob(label,gp=gpar(fontfamily="sans",fontsize=size,
                   fontface=if(bold) "bold" else "plain"))),"pt",valueOnly=TRUE)
  }
  title_width <- function(line,id) {
    if (line==id) return(text_width(id,9,TRUE))
    if (startsWith(line,paste0(id," "))) {
      return(text_width(id,9,TRUE)+text_width(" ")+text_width(substring(line,nchar(id)+2)))
    }
    text_width(line)
  }
  wrap_title <- function(id,label,width) {
    words <- strsplit(label," ",fixed=TRUE)[[1]]
    lines <- character(); line <- id
    for (word in words) {
      candidate <- paste(line,word)
      if (title_width(candidate,id)>width) {lines <- c(lines,line); line <- word} else line <- candidate
    }
    lines <- c(lines,line)
    if (length(lines)>2) stop(paste("D full title does not fit two lines at the retained font sizes:",id,label))
    lines
  }
  line_at <- function(x1,y1,x2,y2,col="#374151",lwd=0.5,lty=1) {
    grid.lines(x=unit(c(x1,x2),"pt"),y=unit(H-c(y1,y2),"pt"),gp=gpar(col=col,lwd=lwd,lty=lty))
  }
  text_at("d",3,8,13,"left",TRUE,"black")
  step <- (W-34)/6; plot_w <- step-6
  yy <- 32+(seq_along(example_regions)-1)*8.05
  for (j in seq_along(example_regions)) text_at(example_regions[[j]],27,yy[[j]],8,just="right")
  line_at(31,yy[1],31,yy[length(yy)])
  for (y in yy) line_at(28.5,y,31,y)
  for (k in seq_along(example_programs)) {
    left <- 33+(k-1)*step
    cx <- left+plot_w/2
    id <- example_programs[[k]]
    lines <- wrap_title(id,example_names[[k]],step-1)
    first_x <- cx-title_width(lines[[1]],id)/2
    text_at(id,first_x,4.5,9,"left",TRUE)
    rest <- if (lines[[1]]==id) "" else substring(lines[[1]],nchar(id)+2)
    if (nzchar(rest)) text_at(rest,first_x+text_width(id,9,TRUE)+text_width(" "),4.7,8.5,"left")
    if (length(lines)==2) text_at(lines[[2]],cx,14.5,8.5)
    xpos <- function(v) left+(v-example_xlim[1])/diff(example_xlim)*plot_w
    line_at(left,29.3,left+plot_w,29.3)
    for (tick in example_ticks) {
      x <- xpos(tick); line_at(x,28.4,x,29.3)
      text_at(format(tick,trim=TRUE),x,23.2,7.5)
    }
    zero <- xpos(0)
    line_at(zero,yy[1]-1,zero,yy[length(yy)]+1,"#A1A7AE",0.45,2)
    for (j in seq_along(example_regions)) {
      vx <- xpos(example_values[k,j])
      line_at(zero,yy[[j]],vx,yy[[j]],example_colors[[j]],0.8)
      grid.points(x=unit(vx,"pt"),y=unit(H-yy[[j]],"pt"),pch=16,
                  size=unit(2.6,"pt"),gp=gpar(col=example_colors[[j]]))
    }
  }
  text_at("Within-program standardized adjusted EMM (profile_z); descriptive, not region-level significance",W/2,143.4,7.5)
  upViewport()
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
heatmap_pdf <- file.path(out_dir, "Fig2_donor_adjusted_regional_heatmap.pdf")
heatmap_png <- file.path(out_dir, "Fig2_donor_adjusted_regional_heatmap.png")
examples_pdf <- file.path(out_dir, "Fig2_donor_adjusted_examples.pdf")
examples_png <- file.path(out_dir, "Fig2_donor_adjusted_examples.png")

examples_only <- "--examples-only" %in% commandArgs(trailingOnly = TRUE)
if (!examples_only) {
  cairo_pdf(heatmap_pdf, width = 531 / 72, height = 240 / 72, family = "sans", onefile = FALSE)
  render_heatmap_panel()
  invisible(dev.off())
  
  grDevices::png(
    heatmap_png, width = round(531 / 72 * 300), height = round(240 / 72 * 300),
    units = "px", res = 300, type = "cairo", bg = "white"
  )
  render_heatmap_panel()
  invisible(dev.off())
}

heatmap_only <- "--heatmap-only" %in% commandArgs(trailingOnly = TRUE)
if (!heatmap_only) {
  cairo_pdf(examples_pdf, width = 531 / 72, height = 148 / 72, family = "sans", onefile = FALSE)
  render_examples_panel()
  invisible(dev.off())

  grDevices::png(
    examples_png, width = round(531 / 72 * 300), height = round(148 / 72 * 300),
    units = "px", res = 300, type = "cairo", bg = "white"
  )
  render_examples_panel()
  invisible(dev.off())
}

if (examples_only) {
  cat("Produced", examples_pdf, examples_png, "\n")
} else if (heatmap_only) {
  cat("Produced", heatmap_pdf, heatmap_png, "\n")
} else {
  cat("Produced", heatmap_pdf, heatmap_png, examples_pdf, examples_png, "\n")
}

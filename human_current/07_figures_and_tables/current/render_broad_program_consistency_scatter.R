if (!nzchar(Sys.getenv("NMF_WORK_ROOT"))) stop("Set NMF_WORK_ROOT; original source directories must remain read-only.")
#!/usr/bin/env Rscript
options(future.globals.maxSize = 50 * 1024^3)

suppressPackageStartupMessages({
  library(grid)
})

SOURCE <- paste0(Sys.getenv("NMF_WORK_ROOT"), "/analysis/external_single_cell_program_scores/broad_program_consistency_scatter.tsv")
OUTDIR <- paste0(Sys.getenv("NMF_WORK_ROOT"), "/analysis/external_single_cell_program_scores")
PDF_PATH <- file.path(OUTDIR, "broad_program_consistency_scatter.pdf")
PNG_PATH <- file.path(OUTDIR, "broad_program_consistency_scatter.png")
SUPPLEMENT_PATH <- file.path(OUTDIR, "external_score_consistency_supplement.pdf")
BY_STUDY_CELLTYPE_PATH <- file.path(OUTDIR, "broad_program_consistency_by_study.tsv")
CELL_NEAREST_BOX_PDF <- file.path(OUTDIR, "cell_nearest_reference_cosine_boxplot.pdf")
CELL_NEAREST_BOX_PNG <- file.path(OUTDIR, "cell_nearest_reference_cosine_boxplot.png")
CELL_NEAREST_OVERALL_PDF <- file.path(OUTDIR, "cell_nearest_reference_cosine_overall_boxplot.pdf")
CELL_NEAREST_OVERALL_PNG <- file.path(OUTDIR, "cell_nearest_reference_cosine_overall_boxplot.png")
HEATMAP_DIR <- file.path(OUTDIR, "program_score_heatmaps_by_dataset")
HEATMAP_NAMES <- c(
  "Reference_1M",
  "Allen",
  "SEAAD",
  "Schirmer_2019_MS_UCSC",
  "SingleSoma_AD_PFC_CELLxGENE",
  "Tran_2021_reward_cortex_CELLxGENE",
  "GSE144136",
  "GSE174367",
  "GSE291605"
)

STUDY_LEVELS <- c(
  "Allen",
  "SEAAD",
  "Schirmer_2019_MS_UCSC",
  "SingleSoma_AD_PFC_CELLxGENE",
  "Tran_2021_reward_cortex_CELLxGENE",
  "GSE144136",
  "GSE174367",
  "GSE291605"
)
CLASS_LEVELS <- c("Exc", "Inh", "Ast", "Oligo", "OPC", "Micro", "Endo", "VLMC")
CLASS_COL <- c(
  Exc = "#1f77b4",
  Inh = "#d62728",
  Ast = "#ff7f0e",
  Oligo = "#2ca02c",
  OPC = "#9467bd",
  Micro = "#8c564b",
  Endo = "#17becf",
  VLMC = "#e377c2"
)


if ("--cell-nearest-reference-overall-8class" %in% commandArgs(trailingOnly = TRUE)) {
  nearest_paths <- setNames(
    file.path(OUTDIR, paste0("cell_nearest_reference_", STUDY_LEVELS, ".parquet")),
    STUDY_LEVELS
  )
  nearest_columns <- c("broad_class", "cosine_similarity")
  if (requireNamespace("arrow", quietly = TRUE)) {
    read_nearest_8class <- function(path) {
      arrow::read_parquet(path, col_select = nearest_columns, as_data_frame = TRUE)
    }
  } else if (requireNamespace("reticulate", quietly = TRUE)) {
    reticulate::use_python("/usr/bin/python3", required = TRUE)
    pyarrow_parquet <- reticulate::import("pyarrow.parquet", convert = FALSE)
    read_nearest_8class <- function(path) {
      table <- pyarrow_parquet$read_table(path, columns = reticulate::r_to_py(as.list(nearest_columns)))
      reticulate::py_to_r(table$to_pandas())
    }
  } else {
    stop("arrow or reticulate is required to read the approved cell-nearest Parquet outputs")
  }
  nearest_table <- do.call(
    rbind,
    lapply(STUDY_LEVELS, function(study_name) read_nearest_8class(nearest_paths[[study_name]]))
  )
  nearest_table$broad_class <- as.character(nearest_table$broad_class)
  nearest_table$cosine_similarity <- as.numeric(nearest_table$cosine_similarity)
  nearest_table <- nearest_table[nearest_table$broad_class %in% CLASS_LEVELS, , drop = FALSE]
  values_by_class <- lapply(CLASS_LEVELS, function(class_name) {
    values <- nearest_table$cosine_similarity[nearest_table$broad_class == class_name]
    values[is.finite(values)]
  })
  names(values_by_class) <- CLASS_LEVELS
  total_counts <- vapply(CLASS_LEVELS, function(class_name) {
    sum(nearest_table$broad_class == class_name)
  }, integer(1))
  valid_counts <- vapply(values_by_class, length, integer(1))
  na_counts <- total_counts - valid_counts
  values_horizontal <- rev(values_by_class)
  class_horizontal <- rev(CLASS_LEVELS)
  colors_horizontal <- rev(unname(CLASS_COL[CLASS_LEVELS]))

  draw_query_class_boxplot <- function() {
    op <- par(
      mar = c(4.1, 5.0, 2.9, 0.5),
      oma = c(1.5, 0.1, 0.1, 0.1),
      family = "sans",
      bg = "white",
      xpd = NA
    )
    on.exit(par(op), add = TRUE)
    y_positions <- seq_along(class_horizontal)
    plot(
      NA,
      type = "n",
      xlim = c(0.0, 1.0),
      ylim = c(0.5, length(y_positions) + 0.5),
      xaxt = "n",
      yaxt = "n",
      xlab = "Cosine similarity (0-1)",
      ylab = "",
      cex.lab = 0.56,
      bty = "n"
    )
    abline(v = seq(0.0, 1.0, by = 0.2), col = "#edf1f5", lwd = 0.5)

    # Every finite cell is sent through this deterministic y-jitter mapping;
    # x-values are never jittered and values below 0 remain in the source
    # vectors/statistics but fall outside the displayed window.
    set.seed(20260912)
    nx <- 900L
    ny <- 900L
    point_raster <- matrix(NA_character_, nrow = ny, ncol = nx)
    for (i in y_positions) {
      values <- values_horizontal[[i]]
      if (!length(values)) next
      jitter_y <- y_positions[[i]] + runif(length(values), min = -0.22, max = 0.22)
      keep <- is.finite(values) & values >= 0.0 & values <= 1.0
      if (!any(keep)) next
      x_idx <- 1L + floor((values[keep] - 0.0) / 1.0 * (nx - 1L))
      y_idx <- 1L + floor((jitter_y[keep] - 0.5) / 8.0 * (ny - 1L))
      inside <- x_idx >= 1L & x_idx <= nx & y_idx >= 1L & y_idx <= ny
      if (!any(inside)) next
      x_idx <- x_idx[inside]
      y_idx <- ny - y_idx[inside] + 1L
      point_raster[cbind(y_idx, x_idx)] <- grDevices::adjustcolor(
        colors_horizontal[[i]], alpha.f = 0.34
      )
    }
    rasterImage(
      point_raster,
      xleft = 0.0,
      ybottom = 0.5,
      xright = 1.0,
      ytop = 8.5,
      interpolate = FALSE
    )

    # Draw the original R boxplot.stats whisker/box definition on top of the
    # raster points; all finite values (including <0) were retained here.
    box_height <- 0.34
    cap_height <- 0.12
    for (i in y_positions) {
      values <- values_horizontal[[i]]
      if (!length(values)) next
      stats <- as.numeric(boxplot.stats(values, coef = 1.5)$stats)
      y <- y_positions[[i]]
      rect(
        stats[2], y - box_height / 2,
        stats[4], y + box_height / 2,
        col = grDevices::adjustcolor(colors_horizontal[[i]], alpha.f = 0.76),
        border = "#253447"
      )
      segments(stats[1], y, stats[2], y, col = "#253447")
      segments(stats[4], y, stats[5], y, col = "#253447")
      segments(stats[1], y - cap_height, stats[1], y + cap_height, col = "#253447")
      segments(stats[5], y - cap_height, stats[5], y + cap_height, col = "#253447")
      segments(stats[3], y - box_height / 2, stats[3], y + box_height / 2, col = "#253447", lwd = 1.1)
    }
    axis(
      1,
      at = seq(0.0, 1.0, by = 0.2),
      labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"),
      cex.axis = 0.50,
      gap.axis = 0
    )
    axis(
      2,
      at = y_positions,
      labels = class_horizontal,
      las = 1,
      cex.axis = 0.52,
      tick = FALSE,
      gap.axis = 0
    )
    mtext(
      "GEP score similarity",
      side = 3,
      line = 1.0,
      cex = 0.74,
      font = 2,
      col = "#111827"
    )
    mtext(
      sprintf("8 query classes; valid=%.2fM; NA=%s", sum(valid_counts) / 1e6, format(sum(na_counts), big.mark = ",")),
      side = 3,
      line = 0.05,
      cex = 0.42,
      col = "#56616f"
    )
    mtext(
      "All finite query cells shown with y-jitter; unrestricted 1M reference; not independent validation.",
      side = 1,
      outer = TRUE,
      line = 0.05,
      cex = 0.34,
      col = "#56616f"
    )
  }

  pdf(
    CELL_NEAREST_OVERALL_PDF,
    width = 175 / 72,
    height = 196 / 72,
    bg = "white",
    useDingbats = FALSE
  )
  draw_query_class_boxplot()
  dev.off()
  pdf_only <- "--pdf-only" %in% commandArgs(trailingOnly = TRUE)
  if (!pdf_only) {
    png(
      CELL_NEAREST_OVERALL_PNG,
      width = 175 * 8,
      height = 196 * 8,
      res = 576,
      type = "cairo-png",
      bg = "white"
    )
    draw_query_class_boxplot()
    dev.off()
  }
  cat(sprintf(
    "CELL_NEAREST_REFERENCE_OVERALL_8CLASS_COMPLETE\nvalid_total\t%d\nNA_total\t%d\noutput\t%s\n",
    sum(valid_counts),
    sum(na_counts),
    CELL_NEAREST_OVERALL_PDF
  ))
  if (!pdf_only) cat(sprintf("output\t%s\n", CELL_NEAREST_OVERALL_PNG))
  for (class_name in CLASS_LEVELS) {
    cat(sprintf(
      "class\t%s\tvalid\t%d\tNA\t%d\n",
      class_name,
      valid_counts[[class_name]],
      na_counts[[class_name]]
    ))
  }
  quit(save = "no", status = 0, runLast = FALSE)
}

if ("--cell-nearest-reference-overall-and-supplement" %in% commandArgs(trailingOnly = TRUE)) {
  nearest_paths <- setNames(
    file.path(OUTDIR, paste0("cell_nearest_reference_", STUDY_LEVELS, ".parquet")),
    STUDY_LEVELS
  )
  nearest_columns <- c("cosine_similarity")
  if (requireNamespace("arrow", quietly = TRUE)) {
    read_nearest_overall <- function(path) {
      arrow::read_parquet(path, col_select = nearest_columns, as_data_frame = TRUE)
    }
  } else if (requireNamespace("reticulate", quietly = TRUE)) {
    reticulate::use_python("/usr/bin/python3", required = TRUE)
    pyarrow_parquet <- reticulate::import("pyarrow.parquet", convert = FALSE)
    read_nearest_overall <- function(path) {
      table <- pyarrow_parquet$read_table(path, columns = reticulate::r_to_py(as.list(nearest_columns)))
      reticulate::py_to_r(table$to_pandas())
    }
  } else {
    stop("arrow or reticulate is required to read the approved cell-nearest Parquet outputs")
  }
  cosine_values <- unlist(
    lapply(STUDY_LEVELS, function(study_name) {
      frame <- read_nearest_overall(nearest_paths[[study_name]])
      as.numeric(frame$cosine_similarity)
    }),
    use.names = FALSE
  )
  valid_values <- cosine_values[is.finite(cosine_values)]
  n_total <- length(cosine_values)
  n_valid <- length(valid_values)
  n_na <- n_total - n_valid
  pdf(CELL_NEAREST_OVERALL_PDF, width = 4.8, height = 2.4, bg = "white", useDingbats = FALSE)
  par(mar = c(3.7, 4.5, 2.5, 0.6), family = "sans", bg = "white")
  boxplot(
    list(`All retained query cells` = valid_values),
    ylim = c(0, 1),
    col = "#2f5d8c",
    border = "#243447",
    outline = FALSE,
    ylab = "Nearest-reference cosine similarity",
    main = "Cross-dataset applicability of GEP scores",
    cex.main = 0.82,
    cex.lab = 0.78,
    cex.axis = 0.72,
    bty = "n"
  )
  mtext(sprintf("Cell-pooled valid=%d; NA=%d", n_valid, n_na), side = 3, line = 0.1, cex = 0.64, col = "#56616f")
  mtext("Unrestricted 1M reference; not independent validation.", side = 1, line = 2.5, cex = 0.62, col = "#56616f")
  dev.off()
  png(CELL_NEAREST_OVERALL_PNG, width = 1440, height = 720, res = 300, type = "cairo-png", bg = "white")
  par(mar = c(3.7, 4.5, 2.5, 0.6), family = "sans", bg = "white")
  boxplot(
    list(`All retained query cells` = valid_values),
    ylim = c(0, 1),
    col = "#2f5d8c",
    border = "#243447",
    outline = FALSE,
    ylab = "Nearest-reference cosine similarity",
    main = "Cross-dataset applicability of GEP scores",
    cex.main = 0.82,
    cex.lab = 0.78,
    cex.axis = 0.72,
    bty = "n"
  )
  mtext(sprintf("Cell-pooled valid=%d; NA=%d", n_valid, n_na), side = 3, line = 0.1, cex = 0.64, col = "#56616f")
  mtext("Unrestricted 1M reference; not independent validation.", side = 1, line = 2.5, cex = 0.62, col = "#56616f")
  dev.off()

  grouped_box_path <- CELL_NEAREST_BOX_PDF
  heatmap_pdfs <- file.path(HEATMAP_DIR, paste0(HEATMAP_NAMES, ".pdf"))
  compose_py <- paste(c(
    "import io",
    "import sys",
    "from pypdf import PdfReader, PdfWriter",
    "from reportlab.pdfgen import canvas",
    "output_path = sys.argv[1]",
    "anchor_path = sys.argv[2]",
    "heatmap_paths = sys.argv[3:]",
    "page_width = 19 * 72",
    "page_height = 21 * 72",
    "margin = 0.4 * 72",
    "column_gap = 0.2 * 72",
    "row_gap = 0.3 * 72",
    "column_width = (page_width - 2 * margin - 3 * column_gap) / 4.0",
    "row_height = (page_height - 2 * margin - 2 * row_gap) / 3.0",
    "row_tops = (page_height - margin, page_height - margin - row_height - row_gap, page_height - margin - 2 * (row_height + row_gap))",
    "first_x = [margin + i * (column_width + column_gap) for i in range(4)]",
    "heatmap_scale = 0.80",
    "letters = list('ABCDEFGHIJ')",
    "source_readers = [PdfReader(path) for path in [anchor_path] + heatmap_paths]",
    "heatmap_pages = [reader.pages[0] for reader in source_readers[1:]]",
    "def fit_page(source_page, cell_x, top_y, cell_width, max_height):",
    "    source_w = float(source_page.mediabox.width)",
    "    source_h = float(source_page.mediabox.height)",
    "    scale = min(cell_width / source_w, max_height / source_h)",
    "    placed_w = source_w * scale",
    "    placed_h = source_h * scale",
    "    placed_x = cell_x + (cell_width - placed_w) / 2.0",
    "    placed_y = top_y - placed_h",
    "    return source_page, placed_x, placed_y, placed_w, placed_h, scale",
    "def fixed_page(source_page, cell_x, top_y, cell_width, scale):",
    "    source_w = float(source_page.mediabox.width)",
    "    source_h = float(source_page.mediabox.height)",
    "    placed_w = source_w * scale",
    "    placed_h = source_h * scale",
    "    placed_x = cell_x + (cell_width - placed_w) / 2.0",
    "    placed_y = top_y - placed_h",
    "    return source_page, placed_x, placed_y, placed_w, placed_h, scale",
    "a_width = 3 * column_width + 2 * column_gap",
    "placements = []",
    "placements.append(fit_page(source_readers[0].pages[0], first_x[0], row_tops[0], a_width, row_height))",
    "placements.append(fixed_page(heatmap_pages[0], first_x[3], row_tops[0], column_width, heatmap_scale))",
    "for idx in range(4):",
    "    placements.append(fixed_page(heatmap_pages[idx + 1], first_x[idx], row_tops[1], column_width, heatmap_scale))",
    "for idx in range(4):",
    "    placements.append(fixed_page(heatmap_pages[idx + 5], first_x[idx], row_tops[2], column_width, heatmap_scale))",
    "overlay_buffer = io.BytesIO()",
    "overlay = canvas.Canvas(overlay_buffer, pagesize=(page_width, page_height))",
    "overlay.setFillColorRGB(0.08, 0.10, 0.14)",
    "overlay.setFont('Helvetica-Bold', 14)",
    "for letter, (_, placed_x, placed_y, placed_w, placed_h, _) in zip(letters, placements):",
    "    overlay.drawString(placed_x, placed_y + placed_h + 6.0, letter)",
    "overlay.save()",
    "overlay_buffer.seek(0)",
    "writer = PdfWriter()",
    "page = PdfReader(overlay_buffer).pages[0]",
    "for source_page, placed_x, placed_y, placed_w, placed_h, scale in placements:",
    "    page.merge_transformed_page(source_page, (scale, 0.0, 0.0, scale, placed_x, placed_y), over=True)",
    "writer.add_page(page)",
    "with open(output_path, 'wb') as handle:",
    "    writer.write(handle)"
  ), collapse = "\n")
  compose_status <- system2(
    "/usr/bin/python3",
    args = c("-c", shQuote(compose_py), SUPPLEMENT_PATH, grouped_box_path, heatmap_pdfs),
    stdout = TRUE,
    stderr = TRUE
  )
  if (!file.exists(SUPPLEMENT_PATH)) {
    cat(paste(compose_status, collapse = "\n"), "\n")
    stop("vector composition failed while creating the approved supplement")
  }
  cat(sprintf(
    "CELL_NEAREST_REFERENCE_OVERALL_COMPLETE\\nvalid\\t%d\\nNA\\t%d\\noutput\\t%s\\noutput\\t%s\\noutput\\t%s\\n",
    n_valid,
    n_na,
    CELL_NEAREST_OVERALL_PDF,
    CELL_NEAREST_OVERALL_PNG,
    SUPPLEMENT_PATH
  ))
  quit(save = "no", status = 0, runLast = FALSE)
}

if ("--cell-nearest-reference-boxplot" %in% commandArgs(trailingOnly = TRUE)) {
  nearest_paths <- setNames(
    file.path(OUTDIR, paste0("cell_nearest_reference_", STUDY_LEVELS, ".parquet")),
    STUDY_LEVELS
  )
  nearest_columns <- c("study", "broad_class", "cosine_similarity", "na_reason")
  if (requireNamespace("arrow", quietly = TRUE)) {
    read_nearest <- function(path) {
      arrow::read_parquet(path, col_select = nearest_columns, as_data_frame = TRUE)
    }
  } else if (requireNamespace("reticulate", quietly = TRUE)) {
    reticulate::use_python("/usr/bin/python3", required = TRUE)
    pyarrow_parquet <- reticulate::import("pyarrow.parquet", convert = FALSE)
    read_nearest <- function(path) {
      table <- pyarrow_parquet$read_table(path, columns = nearest_columns)
      reticulate::py_to_r(table$to_pandas())
    }
  } else {
    stop("arrow or reticulate is required to read the approved cell-nearest Parquet outputs")
  }
  nearest_table <- do.call(
    rbind,
    lapply(STUDY_LEVELS, function(study_name) read_nearest(nearest_paths[[study_name]]))
  )
  nearest_table$study <- factor(as.character(nearest_table$study), levels = STUDY_LEVELS)
  nearest_table$broad_class <- factor(as.character(nearest_table$broad_class), levels = CLASS_LEVELS)
  study_short <- c("Allen", "SEAAD", "Schirmer", "SingleSoma", "Tran", "GSE144", "GSE174", "GSE291")
  draw_nearest_boxplot <- function() {
    op <- par(
      mfrow = c(2, 4),
      mar = c(5.9, 3.8, 2.8, 0.6),
      oma = c(4.2, 0.2, 2.2, 0.2),
      family = "sans",
      bg = "white"
    )
    on.exit(par(op), add = TRUE)
    for (class_name in CLASS_LEVELS) {
      class_rows <- nearest_table[as.character(nearest_table$broad_class) == class_name, , drop = FALSE]
      values_by_study <- lapply(STUDY_LEVELS, function(study_name) {
        values <- class_rows$cosine_similarity[as.character(class_rows$study) == study_name]
        values[is.finite(values)]
      })
      valid_counts <- vapply(values_by_study, length, integer(1))
      total_counts <- vapply(STUDY_LEVELS, function(study_name) {
        sum(as.character(class_rows$study) == study_name)
      }, integer(1))
      present <- valid_counts > 0L
      if (any(present)) {
        boxplot(
          values_by_study[present],
          names = paste0(study_short[present], "\n", valid_counts[present], "/", total_counts[present]),
          ylim = c(0, 1),
          col = CLASS_COL[class_name],
          border = "#374151",
          outline = FALSE,
          main = class_name,
          col.main = CLASS_COL[class_name],
          cex.main = 0.95,
          cex.axis = 0.56,
          las = 2,
          ylab = "Cosine similarity",
          cex.lab = 0.76,
          bty = "n"
        )
      } else {
        plot.new()
        plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
        text(0.5, 0.54, class_name, cex = 0.95, font = 2, col = CLASS_COL[class_name])
        text(0.5, 0.44, "No retained cells", cex = 0.72, col = "#56616f")
      }
      mtext(
        sprintf("valid=%d; NA=%d", sum(valid_counts), sum(total_counts) - sum(valid_counts)),
        side = 3,
        line = 0.15,
        cex = 0.60,
        col = "#56616f"
      )
    }
    mtext(
      "Nearest-reference cosine similarity by study and broad class",
      side = 3,
      outer = TRUE,
      line = 0.55,
      cex = 1.08,
      font = 2,
      col = "#111827"
    )
    mtext(
      "Cell-level best-match distribution; unrestricted 1M reference; not independent validation.",
      side = 1,
      outer = TRUE,
      line = 1.05,
      cex = 0.76,
      col = "#56616f"
    )
  }
  pdf(CELL_NEAREST_BOX_PDF, width = 16, height = 10, bg = "white", useDingbats = FALSE)
  draw_nearest_boxplot()
  dev.off()
  png(CELL_NEAREST_BOX_PNG, width = 4800, height = 3000, res = 300, type = "cairo-png", bg = "white")
  draw_nearest_boxplot()
  dev.off()
  cat(sprintf(
    "CELL_NEAREST_REFERENCE_BOXPLOT_COMPLETE\nrows\t%d\noutput\t%s\noutput\t%s\n",
    nrow(nearest_table),
    CELL_NEAREST_BOX_PDF,
    CELL_NEAREST_BOX_PNG
  ))
  quit(save = "no", status = 0, runLast = FALSE)
}

source_table <- read.delim(
  SOURCE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = "NA"
)
source_table$x_reference_score <- as.numeric(source_table$x_reference_score)
source_table$y_external_score <- as.numeric(source_table$y_external_score)
source_table$broad_class <- factor(source_table$broad_class, levels = CLASS_LEVELS)
plot_table <- source_table[
  source_table$status == "OK" &
    is.finite(source_table$x_reference_score) &
    is.finite(source_table$y_external_score),
  ,
  drop = FALSE
]
if (nrow(plot_table) == 0) {
  stop("no finite broad-class/program score pairs are available")
}

if ("--by-study-celltype-only" %in% commandArgs(trailingOnly = TRUE)) {
  result_rows <- vector("list", length(STUDY_LEVELS) * length(CLASS_LEVELS))
  result_i <- 1L
  for (study_name in STUDY_LEVELS) {
    for (class_name in CLASS_LEVELS) {
      subset <- plot_table[
        plot_table$study == study_name &
          as.character(plot_table$broad_class) == class_name,
        ,
        drop = FALSE
      ]
      n_pairs <- nrow(subset)
      rho <- NA_real_
      na_reason <- ""
      cosine_similarity <- NA_real_
      cosine_na_reason <- ""
      if (n_pairs == 0L) {
        na_reason <- "no_celltype"
      } else if (n_pairs < 2L) {
        na_reason <- "insufficient_pairs"
      } else {
        reference_unique <- length(unique(subset$x_reference_score))
        external_unique <- length(unique(subset$y_external_score))
        if (reference_unique < 2L && external_unique < 2L) {
          na_reason <- "constant_reference_and_external"
        } else if (reference_unique < 2L) {
          na_reason <- "constant_reference"
        } else if (external_unique < 2L) {
          na_reason <- "constant_external"
        } else {
          rho <- suppressWarnings(cor(
            subset$x_reference_score,
            subset$y_external_score,
            method = "spearman"
          ))
          if (is.na(rho)) {
            na_reason <- "undefined"
          }
        }
      }
      if (n_pairs == 0L) {
        cosine_na_reason <- "no_celltype"
      } else {
        reference_values <- subset$x_reference_score
        external_values <- subset$y_external_score
        reference_scale <- max(abs(reference_values))
        external_scale <- max(abs(external_values))
        if (reference_scale == 0 && external_scale == 0) {
          cosine_na_reason <- "zero_reference_and_external_norm"
        } else if (reference_scale == 0) {
          cosine_na_reason <- "zero_reference_norm"
        } else if (external_scale == 0) {
          cosine_na_reason <- "zero_external_norm"
        } else {
          reference_scaled <- reference_values / reference_scale
          external_scaled <- external_values / external_scale
          denominator <- sqrt(sum(reference_scaled * reference_scaled)) * sqrt(sum(external_scaled * external_scaled))
          if (!is.finite(denominator) || denominator == 0) {
            cosine_na_reason <- "zero_norm"
          } else {
            cosine_similarity <- sum(reference_scaled * external_scaled) / denominator
            if (!is.finite(cosine_similarity)) {
              cosine_similarity <- NA_real_
              cosine_na_reason <- "undefined"
            }
          }
        }
      }
      result_rows[[result_i]] <- data.frame(
        study = study_name,
        broad_class = class_name,
        spearman_rho = rho,
        cosine_similarity = cosine_similarity,
        n_pairs = n_pairs,
        na_reason = na_reason,
        cosine_na_reason = cosine_na_reason,
        stringsAsFactors = FALSE
      )
      result_i <- result_i + 1L
    }
  }
  by_study_celltype <- do.call(rbind, result_rows)
  write.table(
    by_study_celltype,
    BY_STUDY_CELLTYPE_PATH,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
  cat(sprintf("BROAD_PROGRAM_CONSISTENCY_BY_STUDY_CELLTYPE_COMPLETE\nrows\t%d\noutput\t%s\n", nrow(by_study_celltype), BY_STUDY_CELLTYPE_PATH))
  for (i in seq_len(nrow(by_study_celltype))) {
    row <- by_study_celltype[i, , drop = FALSE]
    rho_text <- if (is.na(row$spearman_rho)) "NA" else sprintf("%.6f", row$spearman_rho)
    reason_text <- if (nzchar(row$na_reason)) row$na_reason else "OK"
    cosine_text <- if (is.na(row$cosine_similarity)) "NA" else sprintf("%.6f", row$cosine_similarity)
    cosine_reason_text <- if (nzchar(row$cosine_na_reason)) row$cosine_na_reason else "OK"
    cat(sprintf(
      "study_rho\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n",
      row$study,
      row$broad_class,
      rho_text,
      row$n_pairs,
      reason_text,
      cosine_text,
      cosine_reason_text
    ))
  }
  quit(save = "no", status = 0, runLast = FALSE)
}

x <- plot_table$x_reference_score
y <- plot_table$y_external_score
limits <- range(c(x, y), finite = TRUE)
span <- diff(limits)
pad <- if (span > 0) 0.05 * span else 0.05
limits <- c(limits[1] - pad, limits[2] + pad)

format_p <- function(value) {
  if (is.na(value)) return("NA")
  if (value < 1e-3) return(format(value, scientific = TRUE, digits = 2))
  sprintf("%.3f", value)
}
overall_test <- suppressWarnings(cor.test(
  x,
  y,
  method = "spearman",
  alternative = "two.sided",
  exact = FALSE
))
overall_rho <- unname(overall_test$estimate)
overall_p <- unname(overall_test$p.value)
draw_plot <- function() {
  op <- par(
    family = "sans",
    bg = "white"
  )
  on.exit({
    layout(1)
    par(op)
  }, add = TRUE)
  layout(
    matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
    widths = c(0.64, 0.36),
    heights = c(0.87, 0.13)
  )

  par(
    mar = c(4.6, 4.8, 3.0, 0.6) + 0.1,
    mgp = c(2.25, 0.70, 0),
    tcl = -0.25
  )
  plot(
    x,
    y,
    xlim = limits,
    ylim = limits,
    asp = 1,
    pch = 16,
    cex = 0.58,
    col = grDevices::adjustcolor(
      unname(CLASS_COL[as.character(plot_table$broad_class)]),
      alpha.f = 0.48
    ),
    axes = TRUE,
    xlab = "1M reference score (same broad class)",
    ylab = "External study score (same broad class)",
    main = "Cross-study GEP score consistency",
    cex.main = 0.98,
    cex.lab = 0.86,
    cex.axis = 0.78,
    bty = "n"
  )
  abline(a = 0, b = 1, lty = 2, lwd = 1.0, col = "#606060")
  mtext(
    "Each point = study × broad class × program; dashed line = y=x reference",
    side = 3,
    line = 0.12,
    cex = 0.66,
    col = "#56616f"
  )

  par(mar = c(0.2, 0.2, 0.2, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  text(0.03, 0.96, "Overall score association", adj = c(0, 1), cex = 0.88, font = 2, col = "#111827")
  text(
    0.03,
    0.86,
    sprintf("Overall Spearman rho = %.3f; P_nominal = %s", overall_rho, format_p(overall_p)),
    adj = c(0, 1),
    cex = 0.82,
    col = "#111827"
  )
  text(
    0.03,
    0.77,
    sprintf("All finite plotted points (n = %d)", nrow(plot_table)),
    adj = c(0, 1),
    cex = 0.76,
    col = "#56616f"
  )
  text(0.03, 0.63, "Broad-class colors", adj = c(0, 1), cex = 0.80, font = 2, col = "#111827")
  legend_x <- c(0.06, 0.30, 0.54, 0.78)
  legend_y <- c(0.50, 0.34)
  for (i in seq_along(CLASS_LEVELS)) {
    row_i <- ((i - 1) %/% 4) + 1
    col_i <- ((i - 1) %% 4) + 1
    x_i <- legend_x[col_i]
    y_i <- legend_y[row_i]
    points(x_i, y_i, pch = 16, cex = 1.05, col = CLASS_COL[CLASS_LEVELS[i]])
    text(x_i + 0.045, y_i, CLASS_LEVELS[i], adj = c(0, 0.5), cex = 0.78, col = "#111827")
  }

  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(
    0.01,
    0.68,
    "P_nominal: uncorrected two-sided Spearman test across all finite cross-study/class/program points; descriptive only.",
    adj = c(0, 0.5),
    cex = 0.76,
    col = "#56616f"
  )
  text(
    0.01,
    0.30,
    "Shared 1M reference and repeated study × program points are not independent biological replicates.",
    adj = c(0, 0.5),
    cex = 0.76,
    col = "#56616f"
  )
}

A_WIDTH_IN <- 13.6
A_HEIGHT_IN <- 6.5333333333
pdf(PDF_PATH, width = A_WIDTH_IN, height = A_HEIGHT_IN, bg = "white", useDingbats = FALSE)
draw_plot()
dev.off()
png(PNG_PATH, width = 4080, height = 1960, res = 300, type = "cairo-png", bg = "white")
draw_plot()
dev.off()

heatmap_pdfs <- file.path(HEATMAP_DIR, paste0(HEATMAP_NAMES, ".pdf"))
compose_py <- paste(c(
  "import io",
  "import sys",
  "from pypdf import PdfReader, PdfWriter",
  "from reportlab.pdfgen import canvas",
  "output_path = sys.argv[1]",
  "scatter_path = sys.argv[2]",
  "heatmap_paths = sys.argv[3:]",
  "page_width = 19 * 72",
  "page_height = 21 * 72",
  "margin = 0.4 * 72",
  "column_gap = 0.2 * 72",
  "row_gap = 0.3 * 72",
  "column_width = (page_width - 2 * margin - 3 * column_gap) / 4.0",
  "row_height = (page_height - 2 * margin - 2 * row_gap) / 3.0",
  "row_tops = (page_height - margin, page_height - margin - row_height - row_gap, page_height - margin - 2 * (row_height + row_gap))",
  "first_x = [margin + i * (column_width + column_gap) for i in range(4)]",
  "heatmap_scale = 0.80",
  "letters = list('ABCDEFGHIJ')",
  "source_readers = [PdfReader(path) for path in [scatter_path] + heatmap_paths]",
  "scatter_page = source_readers[0].pages[0]",
  "heatmap_pages = [reader.pages[0] for reader in source_readers[1:]]",
  "def fit_page(source_page, cell_x, top_y, cell_width, max_height):",
  "    source_w = float(source_page.mediabox.width)",
  "    source_h = float(source_page.mediabox.height)",
  "    scale = min(cell_width / source_w, max_height / source_h)",
  "    placed_w = source_w * scale",
  "    placed_h = source_h * scale",
  "    placed_x = cell_x + (cell_width - placed_w) / 2.0",
  "    placed_y = top_y - placed_h",
  "    return source_page, placed_x, placed_y, placed_w, placed_h, scale",
  "def fixed_page(source_page, cell_x, top_y, cell_width, scale):",
  "    source_w = float(source_page.mediabox.width)",
  "    source_h = float(source_page.mediabox.height)",
  "    placed_w = source_w * scale",
  "    placed_h = source_h * scale",
  "    placed_x = cell_x + (cell_width - placed_w) / 2.0",
  "    placed_y = top_y - placed_h",
  "    return source_page, placed_x, placed_y, placed_w, placed_h, scale",
  "a_width = 3 * column_width + 2 * column_gap",
  "placements = []",
  "placements.append(fit_page(scatter_page, first_x[0], row_tops[0], a_width, row_height))",
  "placements.append(fixed_page(heatmap_pages[0], first_x[3], row_tops[0], column_width, heatmap_scale))",
  "for idx in range(4):",
  "    placements.append(fixed_page(heatmap_pages[idx + 1], first_x[idx], row_tops[1], column_width, heatmap_scale))",
  "for idx in range(4):",
  "    placements.append(fixed_page(heatmap_pages[idx + 5], first_x[idx], row_tops[2], column_width, heatmap_scale))",
  "overlay_buffer = io.BytesIO()",
  "overlay = canvas.Canvas(overlay_buffer, pagesize=(page_width, page_height))",
  "overlay.setFillColorRGB(0.08, 0.10, 0.14)",
  "overlay.setFont('Helvetica-Bold', 14)",
  "for letter, (_, placed_x, placed_y, placed_w, placed_h, _) in zip(letters, placements):",
  "    overlay.drawString(placed_x, placed_y + placed_h + 6.0, letter)",
  "overlay.save()",
  "overlay_buffer.seek(0)",
  "writer = PdfWriter()",
  "page = PdfReader(overlay_buffer).pages[0]",
  "for source_page, placed_x, placed_y, placed_w, placed_h, scale in placements:",
  "    page.merge_transformed_page(source_page, (scale, 0.0, 0.0, scale, placed_x, placed_y), over=True)",
  "writer.add_page(page)",
  "with open(output_path, 'wb') as handle:",
  "    writer.write(handle)"
), collapse = "\n")
compose_status <- system2(
  "/usr/bin/python3",
  args = c("-c", shQuote(compose_py), SUPPLEMENT_PATH, PDF_PATH, heatmap_pdfs),
  stdout = TRUE,
  stderr = TRUE
)
if (!file.exists(SUPPLEMENT_PATH)) {
  cat(paste(compose_status, collapse = "\n"), "\n")
  stop("vector composition failed while creating the approved supplement")
}
cat(sprintf("BROAD_PROGRAM_CONSISTENCY_RENDER_COMPLETE\npoints\t%d\n", nrow(plot_table)))
cat(sprintf(
  "overall_stat\t%d\trho\t%.6f\tp_nominal\t%s\n",
  nrow(plot_table),
  overall_rho,
  format_p(overall_p)
))
cat(sprintf("supplement_pages\t1\noutput\t%s\noutput\t%s\noutput\t%s\n", PDF_PATH, PNG_PATH, SUPPLEMENT_PATH))

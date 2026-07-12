# Shared theme + palette for Fig4 panels
suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
})

FIG <- "CORTEX_PROGRAM_ROOT/figures/fig4"
RES <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1"

# Map raw cNMF indices to retained P1-P54 labels.
# program_renumber_map.tsv: old_P (int) -> new_P ("EXCLUDED" or "1".."54")
.rmap <- read.delim(file.path(RES, "program_renumber_map.tsv"), stringsAsFactors = FALSE)
.EXCLUDED_OLD <- .rmap$old_P[.rmap$new_P == "EXCLUDED"]  # 9,18,19,35,52,57
.old2new_int  <- setNames(as.integer(.rmap$new_P[.rmap$new_P != "EXCLUDED"]),
                           .rmap$old_P[.rmap$new_P != "EXCLUDED"])  # old_int -> new_int

# ---- program -> functional label lookup (GO:BP names) ----
# program_names.tsv uses cnmf_component (old int) as the canonical key.
# Build label = "P{new} {name_short}" (new P number), "*" if brain-weak.
.prog_names <- read.delim(file.path(RES, "program_names.tsv"),
                          stringsAsFactors = FALSE)
# Filter out excluded programs, map old -> new number
.prog_names <- .prog_names[!.prog_names$cnmf_component %in% .EXCLUDED_OLD, ]
.prog_names$new_prog_int <- .old2new_int[as.character(.prog_names$cnmf_component)]

.mk_label <- function(new_prog, name_short, confidence) {
  ns <- gsub("\\s+P[0-9]+\\s*$", "", trimws(as.character(name_short)))  # strip trailing " P{n}"
  star <- ifelse(confidence == "brain-weak", "*", "")
  paste0("P", new_prog, " ", ns, star)
}
.prog_names$label <- .mk_label(.prog_names$new_prog_int,
                               .prog_names$name_short,
                               .prog_names$confidence)
# PROG_LABEL: keyed by OLD integer (as char), but returns NEW P-labelled string
# This is intentional: callers pass old cNMF index (from data files),
# we return the new P-number label.
PROG_LABEL <- setNames(.prog_names$label, as.character(.prog_names$cnmf_component))

# Full functional label for a program id (old cNMF int, numeric or character).
# Returns "P{new} {name_short}" for kept programs, "P{new}" fallback for unknowns.
prog_label <- function(p) {
  key <- as.character(as.integer(p))
  out <- unname(PROG_LABEL[key])
  # For excluded programs, use new index (should be filtered before calling, but safe fallback)
  new_n <- .old2new_int[key]
  out[is.na(out)] <- ifelse(!is.na(new_n[is.na(out)]), paste0("P", new_n[is.na(out)]), paste0("P?[old", key[is.na(out)], "]"))
  out
}
# Bare "P{new}" label for old cNMF index.
prog_pid <- function(p) {
  new_n <- .old2new_int[as.character(as.integer(p))]
  ifelse(!is.na(new_n), paste0("P", new_n), paste0("P?[old", as.integer(p), "]"))
}
# Selective labeller: full name only for `highlight` ids (OLD cNMF ints), else bare "P{new}".
prog_label_selective <- function(p, highlight) {
  pi <- as.integer(p)
  ifelse(pi %in% as.integer(highlight), prog_label(pi), prog_pid(pi))
}

# neuron class assignment for the 22 subclasses
exc <- c("L6 CT","L6 IT","ET","NP","L6B","L6 CAR3","L3-L4 IT RORB",
         "L4-L5 IT RORB","L2-L3 IT LINC00507")
inh <- c("CHANDELIER","PVALB","SST","LAMP5","PAX6","NDNF","VIP")
nonneu <- c("AST","OPC","OLIGO","MICRO","VLMC","ENDO")
sc_class <- function(x) ifelse(x %in% exc, "Excitatory",
                        ifelse(x %in% inh, "Inhibitory", "Non-neuronal"))

# palette
pal_class <- c(Excitatory = "#C0392B", Inhibitory = "#2471A3", `Non-neuronal` = "#27AE60")
# eta2 continuous palette (viridis-like magma trimmed)
pal_eta_fun <- function() scale_fill_gradientn(
  colours = c("#FCFDBF","#FEC287","#F1605D","#B63679","#721F81","#2D1160"),
  name = expression(eta^2))

theme_fig4 <- function(base = 8) {
  theme_minimal(base_size = base, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
      axis.text = element_text(colour = "grey20"),
      axis.title = element_text(colour = "grey10", size = base + 1),
      plot.title = element_text(face = "bold", size = base + 3, hjust = 0),
      plot.subtitle = element_text(size = base, colour = "grey35"),
      legend.key.size = unit(3.2, "mm"),
      legend.text = element_text(size = base - 1),
      legend.title = element_text(size = base),
      plot.margin = margin(6, 8, 6, 8)
    )
}

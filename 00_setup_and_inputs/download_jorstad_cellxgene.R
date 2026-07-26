PROJECT_ROOT <- Sys.getenv("CORTEX_PROGRAM_ROOT", unset = "/DATA/cortex_nmf_program")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(Sys.getenv("CORTEX_PROGRAM_ROOT", "/DATA/cortex_nmf_program"), mustWork = FALSE)
suppressPackageStartupMessages({
    library(digest)
    library(httr)
    library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
    stop(call. = FALSE)
}
out_dir <- args[[1]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
collection_id <- "d17249d2-0e6e-4500-abb8-e6c93fa1ac6f"
collection_url <- sprintf("https://api.cellxgene.cziscience.com/curation/v1/collections/%s", collection_id)
response <- GET(collection_url)
stop_for_status(response)
collection <- fromJSON(content(response, as = "text", encoding = "UTF-8"), simplifyVector = FALSE)
datasets <- Filter(function(dataset) grepl("Supercluster", dataset$title, fixed = TRUE), collection$datasets)
if (length(datasets) != 5L) {
    stop(call. = FALSE)
}
datasets <- datasets[order(vapply(datasets, function(dataset) dataset$dataset_id, character(1)))]
records <- list()
for (dataset in datasets) {
    assets <- Filter(function(asset) identical(asset$filetype, "H5AD"), dataset$assets)
    if (length(assets) != 1L) {
        stop(call. = FALSE)
    }
    asset <- assets[[1]]
    destination <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9._-]", "_", dataset$title), ".h5ad"))
    file_response <- GET(asset$url, write_disk(destination, overwrite = TRUE), progress())
    stop_for_status(file_response, task = sprintf("download %s", dataset$title))
    bytes <- file.info(destination)$size
    if (is.na(bytes) || bytes <= 0L) {
        stop(call. = FALSE)
    }
    records[[length(records) + 1L]] <- data.frame(dataset_id = dataset$dataset_id, dataset_title = dataset$title, asset_url = asset$url, file_name = basename(destination), bytes = bytes, sha256 = digest(destination, algo = "sha256", file = TRUE), stringsAsFactors = FALSE)
}
source_list <- do.call(rbind, records)
write.table(source_list, file.path(out_dir, "SOURCE_CHECKSUMS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

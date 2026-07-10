#!/usr/bin/env Rscript
# Install the two spatial R packages from immutable Git commits after conda setup.
suppressPackageStartupMessages(library(remotes))

install_github("mojaveazure/seurat-disk", ref = "877d4e18ab38c686f5db54f8cd290274ccdbe295", upgrade = "never")
install_github("dmcable/spacexr", ref = "9f5dc33c8060f946c6072a138b70e189636e1435", upgrade = "never")

if (!requireNamespace("SeuratDisk", quietly = TRUE)) stop("SeuratDisk installation failed", call. = FALSE)
if (!requireNamespace("spacexr", quietly = TRUE)) stop("spacexr installation failed", call. = FALSE)
message("Installed SeuratDisk ", as.character(packageVersion("SeuratDisk")), " and spacexr ", as.character(packageVersion("spacexr")))

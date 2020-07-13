#!/usr/bin/env Rscript

# Code to prepare the reference dataset goes

#' Capitalize the first letter of a string
#'
#' @param x A character vector
#'
#' @return \code{x} with the first letter capitalized
#'
#' @examples
#' Capitalize("word")
#' Capitalize(c("hello", "world"))
#'
Capitalize <- function(x) {
  x <- paste0(
    toupper(x = substr(x = x, start = 1, stop = 1)),
    substr(x = x, start = 2, stop = nchar(x = x))
  )
  return(x)
}

# Check dependencies
message("Checking dependencies")
for (i in c("Seurat", "pkgload", "remotes")) {
  if (!requireNamespace(i, quietly = TRUE)) {
    install.packages(i)
  }
}

for (i in c('satijalab/seurat-data', 'mojaveazure/seurat-disk')) {
  pkg <- Capitalize(x = strsplit(x = basename(path = i), split = '-')[[1]])
  pkg <- paste(pkg, collapse = "")
  if (!requireNamespace(pkg, quietly = TRUE)) {
    remotes::install_github(repo = i)
  }
}

# Make the reference directory
refdir <- file.path(pkgload::package_file(), 'inst', 'references')
out <- file.path(refdir, 'pbmc3k_final.h5Seurat')
message("Creating references directory ", refdir)
dir.create(path = refdir, showWarnings = FALSE, recursive = TRUE)
if (file.exists(out)) {
  warning(
    "Output reference already exists, overwriting",
    call. = FALSE,
    immediate. = TRUE
  )
}

# Get pbmc3k.final
message("Fetching pbmc3k.final")
SeuratData::InstallData("pbmc3k")
data('pbmc3k.final', package = 'pbmc3k.SeuratData')
pbmc3k.final <- Seurat::UpdateSeuratObject(pbmc3k.final)

# Save pbmc3k.final
message("Saving pbmc3k.final to ", out)
SeuratDisk::SaveH5Seurat(
  object = pbmc3k.final,
  filename = out,
  overwrite = TRUE
)

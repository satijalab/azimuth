#!/usr/bin/env Rscript

requireNamespace('jsonlite')
requireNamespace('pkgload')

opts <- list(
  Azimuth.app.reference = 'https://seurat.nygenome.org/references/pbmc',
  welcomebox = "box(\n  h3(\"Please upload a dataset to map to the Multimodal PBMC reference\"),\n  \"Upload a counts matrix from an scRNA-seq dataset of human PBMC in one\n        of the following formats: hdf5, rds, h5ad, h5seurat. For testing, we\n        also provide a demo dataset of 11,769 human PBMC from 10x Genomics,\n        which is loaded automatically with the 'Load demo dataset' button\n        or available for download \",\n  a(\"here\",\n    href=\"https://www.dropbox.com/s/cmbvq2og93lnl9z/pbmc_10k_v3_filtered_feature_bc_matrix.h5?dl=0\",\n    target=\"_blank\"), # open in new browser tab\n  \".\",\n  width = 12\n)",
  Azimuth.app.pbcorthresh = 0.75,
  shiny.maxRequestSize = 524288000
)

jsonlite::write_json(
  x = opts,
  path = pkgload::package_file('inst', 'resources', 'config.json'),
  pretty = 4,
  auto_unbox = TRUE
)

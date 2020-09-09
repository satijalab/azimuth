#!/usr/bin/env Rscript

# Ensure Seurat v4.0 or higher is installed
if (packageVersion(pkg = "Seurat") < package_version(x = "3.9.9000")) {
  stop("Mapping datasets requires Seurat v4 or higher", call. = FALSE)
}

library(Seurat)

# Load helper functions from Azimuth
source("https://raw.githubusercontent.com/mojaveazure/seurat-mapper/master/R/helpers.R")

# Download the multimodal PBMC reference from [LINK]

# Load the reference file
reference <- LoadReference(path = "${ref}")

# Load the query object for mapping
# Change the file path based on where the query file is located on your system.
query <- LoadFileInput(path = "${path}")

# Calculate nCount_RNA and nFeature_RNA if the query does not
# contain them already
if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
    calcn <- as.data.frame(x = Seurat:::CalcN(object = query))
    colnames(x = calcn) <- paste(
      colnames(x = calcn),
      "RNA",
      sep = '_'
    )
    query <- AddMetaData(
      object = query,
      metadata = calcn
    )
    rm(calcn)
}

# Calculate percent mitochondrial genes if the query contains genes
# matching the regular expression "${mito.pattern}"
if (any(grepl(pattern = '${mito.pattern}', x = rownames(x = query)))) {
  query <- PercentageFeatureSet(
    object = query,
    pattern = '${mito.pattern}',
    col.name = '${mito.key}',
    assay = "RNA"
  )
}

# Filter cells based on the thresholds for nCount_RNA and nFeature_RNA
# you set in the app
cells.use <- query[["nCount_RNA", drop = TRUE]] <= ${ncount.max} &&
  query[["nCount_RNA", drop = TRUE]] >= ${ncount.min} &&
  query[["nFeature_RNA", drop = TRUE]] <= ${nfeature.max} &&
  query[["nFeature_RNA", drop = TRUE]] >= ${nfeature.min}

# If the query contains mitochondrial genes, filter cells based on the
# thresholds for ${mito.key} you set in the app
if ("${mito.key}" %in% c(colnames(x = query[[]]))) {
  cells.use <- query[["${mito.key}", drop = TRUE]] <= ${mito.max} &&
    query[["${mito.key}", drop = TRUE]] >= ${mito.min}
}

# Remove filtered cells from the query
query <- query[, cells.use]

# Preprocess with SCTransform
query <- SCTransform(
  object = query,
  assay = "RNA",
  residual.features = rownames(x = reference),
  ncells = ${sct.ncells},
  n_genes = ${sct.nfeats},
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

# Find anchors between query and reference that will be used for mapping
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  mapping = TRUE,
  reference.neighbors = "spca.annoy.neighbors",
  reference.assay = "RNA",
  query.assay = "SCT",
  reference.reduction = "spca",
  normalization.method = "SCT",
  features = rownames(x = reference),
  dims = 1:50,
  nn.method = "annoy",
  verbose = TRUE
)

# Map cells using the anchors just computed. Transfer cell type labels and
# impute protein expression.
mapped <- MapQueryData(
  reference = reference,
  query = query,
  dims = 1:50,
  anchorset = anchors,
  reference.neighbors = "spca.annoy.neighbors",
  transfer.labels = reference$id,
  transfer.expression = GetAssayData(
    object = reference[['ADT']],
    slot = 'data'
  )
)

# The reference used in the app is downsampled compared to the reference on which
# the UMAP model was computed. This step, using the helper function NNTransform,
# corrects the Neighbors of the "mapped" object to account for the downsampling.
mapped <- NNTransform(
  object = mapped,
  meta.data = reference[[]]
)

# Project the query to the reference UMAP.
query[["proj.umap"]] <- RunUMAP(
  object = mapped[["query_ref.nn"]],
  reduction.model = reference[["jumap"]],
  reduction.key = 'UMAP_'
)

# Add the predicted cell types to the query object as metadata
query[["predicted.id"]] <- Misc(object = mapped, slot = "predictions")[,"predicted.id"]

# Add the maximum prediction score (predicted.id.score) and the prediction scores
# for all cell types to the query as an Assay
query[["predictions"]] <- CreateAssayObject(data =  t(x = Misc(object = mapped, slot = "predictions")[, -1]))

# Add the imputed protein to the query object as an Assay
query[['${adt.key}']] <- CreateAssayObject(data = mapped[['transfer']][, colnames(x = query)])


# VISUALIZATIONS

# DimPlot of the reference
DimPlot(object = reference, reduction = "jumap", group.by = "id", label = TRUE) + NoLegend()

# DimPlot of the query, colored by predicted cell type
DimPlot(object = query, reduction = "proj.umap", group.by = "predicted.id", label = TRUE) + NoLegend()

# Plot the score for the predicted cell type of the query
FeaturePlot(object = query, features = "predicted.id.score", reduction = "proj.umap")
VlnPlot(object = query, features = "predicted.id.score", group.by = "predicted.id") + NoLegend()

# Plot the prediction score for the class CD16 Mono
FeaturePlot(object = query, features = "prediction.score.CD16.Mono", reduction = "proj.umap")
VlnPlot(object = query, features = "prediction.score.CD16.Mono", group.by = "predicted.id") + NoLegend()

# Plot an RNA feature
FeaturePlot(object = query, features = "${plotgene}", reduction = "proj.umap")
VlnPlot(object = query, features = "${plotgene}", group.by = "predicted.id") + NoLegend()

# Plot an imputed protein feature
FeaturePlot(object = query, features = "${plotadt}", reduction = "proj.umap")
VlnPlot(object = query, features = "${plotadt}", group.by = "predicted.id") + NoLegend()

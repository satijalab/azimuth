#!/usr/bin/env Rscript

# Ensure Seurat v4.0 or higher is installed
if (packageVersion(pkg = "Seurat") < package_version(x = "3.9.9.9002")) {
  stop("Mapping datasets requires Seurat v4 or higher", call. = FALSE)
}

# Ensure glmGamPoi is installed
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    BiocManager::install("glmGamPoi")
  }
}

library(Seurat)

# Load helper functions from Azimuth
source("https://raw.githubusercontent.com/satijalab/azimuth/master/R/helpers.R")

# Download the multimodal PBMC reference from [LINK] and extract the archive

# Load the reference
# Change the file path based on where the reference is located on your system.
reference <- LoadReference(path = "${ref.uri}")

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
cells.use <- query[["nCount_RNA", drop = TRUE]] <= ${ncount.max} &
  query[["nCount_RNA", drop = TRUE]] >= ${ncount.min} &
  query[["nFeature_RNA", drop = TRUE]] <= ${nfeature.max} &
  query[["nFeature_RNA", drop = TRUE]] >= ${nfeature.min}

# If the query contains mitochondrial genes, filter cells based on the
# thresholds for ${mito.key} you set in the app
if ("${mito.key}" %in% c(colnames(x = query[[]]))) {
  cells.use <- query[["${mito.key}", drop = TRUE]] <= ${mito.max} &
    query[["${mito.key}", drop = TRUE]] >= ${mito.min}
}

# Remove filtered cells from the query
query <- query[, cells.use]

# Preprocess with SCTransform
query <- SCTransform(
  object = query,
  assay = "RNA",
  residual.features = rownames(x = reference$map),
  method = 'glmGamPoi',
  ncells = ${sct.ncells},
  n_genes = ${sct.nfeats},
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

# Find anchors between query and reference
anchors <- FindTransferAnchors(
  reference = reference$map,
  query = query,
  k.filter = NA,
  reference.neighbors = "spca.annoy.neighbors",
  reference.assay = "SCT",
  query.assay = "SCT",
  reference.reduction = "spca",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
  dims = 1:50,
  n.trees = ${ntrees},
  mapping.score.k = 100
)

# Transfer cell type labels and impute protein expression
#
# Transferred labels are in a metadata column named "predicted.id"
# The maximum prediction score is in a metadata column named "predicted.id.score"
# The prediction scores for each class are in an assay named "prediction.score.id"
# The imputed assay is named "impADT"
query <- TransferData(
  reference = reference$map,
  query = query,
  dims = 1:50,
  anchorset = anchors,
  refdata = list(
    id = Idents(reference$map),
    impADT = GetAssayData(
    object = reference$map[['ADT']],
    slot = 'data'
  )),
  store.weights = TRUE
)

# Calculate the embeddings of the query data on the reference SPCA
query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference$map,
  query = query,
  reductions = "pcaproject",
  reuse.weights.matrix = TRUE
)

# Calculate the query neighbors in the reference
# with respect to the integrated embeddings
query[["query_ref.nn"]] <- FindNeighbors(
  object = Embeddings(reference$map[["spca"]]),
  query = Embeddings(query[["integrated_dr"]]),
  return.neighbor = TRUE,
  l2.norm = TRUE
)

# The reference used in the app is downsampled compared to the reference on which
# the UMAP model was computed. This step, using the helper function NNTransform,
# corrects the Neighbors to account for the downsampling.
query <- NNTransform(
  object = query,
  meta.data = reference$map[[]]
)

# Project the query to the reference UMAP.
query[["proj.umap"]] <- RunUMAP(
  object = query[["query_ref.nn"]],
  reduction.model = reference$map[["jumap"]],
  reduction.key = 'UMAP_'
)


# Calculate mapping score and add to metadata
query <- AddMetaData(
  object = query,
  metadata = MappingScore(anchors = anchors),
  col.name = "mapping.score"
)

# VISUALIZATIONS

# DimPlot of the reference
DimPlot(object = reference$plot, reduction = "umap", group.by = "id", label = TRUE) + NoLegend()

# DimPlot of the query, colored by predicted cell type
DimPlot(object = query, reduction = "proj.umap", group.by = "predicted.id", label = TRUE) + NoLegend()

# Plot the score for the predicted cell type of the query
FeaturePlot(object = query, features = "predicted.id.score", reduction = "proj.umap")
VlnPlot(object = query, features = "predicted.id.score", group.by = "predicted.id") + NoLegend()

# Plot the mapping score
FeaturePlot(object = query, features = "mapping.score", reduction = "proj.umap")
VlnPlot(object = query, features = "mapping.score", group.by = "predicted.id") + NoLegend()

# Plot the prediction score for the class CD16 Mono
FeaturePlot(object = query, features = "CD16 Mono", reduction = "proj.umap")
VlnPlot(object = query, features = "CD16 Mono", group.by = "predicted.id") + NoLegend()

# Plot an RNA feature
FeaturePlot(object = query, features = "${plotgene}", reduction = "proj.umap")
VlnPlot(object = query, features = "${plotgene}", group.by = "predicted.id") + NoLegend()

# Plot an imputed protein feature
FeaturePlot(object = query, features = "${plotadt}", reduction = "proj.umap")
VlnPlot(object = query, features = "${plotadt}", group.by = "predicted.id") + NoLegend()

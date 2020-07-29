#!/usr/bin/env Rscript

library(Seurat)
library(SeuratMapper)

# Read in the data and generate QC stats
refs <- SeuratMapper:::LoadReference(path = '${ref.uri}')
object <- SeuratMapper:::LoadFileInput(path = '${path}')
default.assay <- DefaultAssay(object = object)

ncount <- paste0('nCount_', default.assay)
nfeature <- paste0('nFeature_', default.assay)

qc <- c(ncount, nfeature)

if (any(grepl(pattern = '${mito.pattern}', x = rownames(x = object)))) {
  object <- PercentageFeatureSet(
    object = object,
    pattern = '${mito.pattern}',
    col.name = '${mito.key}',
    assay = default.assay
  )
  qc <- append(x = qc, values = '${mito.key}')
}

# VlnPlot(object = object, features = qc, ncol = length(x = qc))
# qc.table <- apply(X = object[[qc]], MARGIN = 2, FUN = quantile)
# qc.table <- as.data.frame(x = qc.table)
# colnames(x = qc.table) <- c(
#   'nUMI per cell',
#   'Genes detected per cell',
#   'Mitochondrial percentage per cell'
# )[seq_along(along.with = qc)]
# qc.table <- t(x = qc.table)
# qc.table

# Preprocess with SCTransform
cells.use <- object[[ncount, drop = TRUE]] <= ${ncount.max} &&
  object[[ncount, drop = TRUE]] >= ${ncount.min} &&
  object[[nfeature, drop = TRUE]] <= ${nfeature.max} &&
  object[[nfeature, drop = TRUE]] >= ${nfeature.min}

object <- object[, cells.use]
object <- SCTransform(
  object = object,
  assay = default.assay,
  residual.features = rownames(x = refs$map),
  ncells = ${sct.ncells},
  n_genes = ${sct.nfeats},
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

# Do the mapping
cells <- colnames(x = object)
anchors <- Seurat:::FindTransferAnchors_Fast(
  reference = refs$map,
  query = object,
  reference.nn = refs$index,
  reference.nnidx = refs$index$annoy_index,
  reference.assay = DefaultAssay(object = object),
  npcs = NULL,
  k.filter = NA,
  query.assay = 'SCT',
  reference.reduction = 'spca',
  normalization.method = 'SCT',
  features = rownames(x = refs$map),
  dims = 1:50
)
ingested <- Seurat:::IngestNewData_Fast(
  reference = refs$map,
  query = object,
  dims = 1:50,
  transfer.anchors = anchors,
  reference.nnidx = refs$index$annoy_index,
  transfer.labels = Idents(object = refs$map),
  transfer.expression = GetAssayData(
    object = refs$map[['ADT']],
    slot = 'data'
  )
)
slot(object = ingested, name = 'neighbors')[['query_ref.nn']] <- SeuratMapper:::NNTransform(
  neighbors = ingested[['query_ref.nn']],
  meta.data = refs$map[[]]
)
object <- SeuratMapper:::AddPredictions(
  object = object,
  preds = ingested$predicted.id,
  scores = ingested$predicted.id.score,
  preds.levels = levels(x = refs$map),
  preds.drop = TRUE
)
object <- RunUMAP(
  object = ingested[['query_ref.nn']],
  reduction.model = refs$map[['jumap']],
  reduction.key = 'ProjU_'
)
object[['${adt.key}']] <- CreateAssayObject(
  data = ingested[['transfer']][, cells]
)
object[['int']] <- CreateDimReducObject(
  embeddings = Embeddings(object = ingested[['int']])[cells, ],
  assay = default.assay
)
dsqr <- SeuratMapper:::QueryReference(
  reference = refs$map,
  query = object,
  assay.query = default.assay
)
object <- AddMetaData(
  object = object,
  metadata = SeuratMapper:::CalcMappingMetric(object = dsqr)
)
app.env$object <- app.env$object[, app.env$object$mapped]

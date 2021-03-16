query=heimmune
reference=pbmc
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

# Preprocess with SCTransform
query <- SCTransform(
  object = query,
  assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = reference),
  reference.SCT.model = reference[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

# Find anchors between query and reference
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference), VariableFeatures(object = query)),
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)

# Transfer cell type labels and impute protein expression
#
# Transferred labels are in metadata columns named "predicted.*"
# The maximum prediction score is in a metadata column named "predicted.*.score"
# The prediction scores for each class are in an assay named "prediction.score.*"
# The imputed assay is named "impADT" if computed

refdata <- lapply(X = c("celltype.l2", "celltype.l1", "celltype.l3"), function(x) {
  reference[[x, drop = TRUE]]
})
names(x = refdata) <- c("celltype.l2", "celltype.l1", "celltype.l3")
if (TRUE) {
  refdata[["impADT"]] <- GetAssayData(
    object = reference[['ADT']],
    slot = 'data'
  )
}
query <- TransferData(
  reference = reference,
  query = query,
  dims = 1:50,
  anchorset = anchors,
  refdata = refdata,
  n.trees = 20,
  store.weights = TRUE
)

# Calculate the embeddings of the query data on the reference SPCA
query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference,
  query = query,
  reductions = "pcaproject",
  reuse.weights.matrix = TRUE
)

# Calculate the query neighbors in the reference
# with respect to the integrated embeddings
query[["query_ref.nn"]] <- FindNeighbors(
  object = Embeddings(reference[["refDR"]]),
  query = Embeddings(query[["integrated_dr"]]),
  return.neighbor = TRUE,
  l2.norm = TRUE
)

# The reference used in the app is downsampled compared to the reference on which
# the UMAP model was computed. This step, using the helper function NNTransform,
# corrects the Neighbors to account for the downsampling.
query <- Azimuth:::NNTransform(
  object = query,
  meta.data = reference[[]]
)

# Project the query to the reference UMAP.
query[["proj.umap"]] <- RunUMAP(
  object = query[["query_ref.nn"]],
  reduction.model = reference[["refUMAP"]],
  reduction.key = 'UMAP_'
)

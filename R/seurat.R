#' @include zzz.R
#' @importFrom Seurat Embeddings
#'
NULL

#' Pseudobulk correlation test of query against reference
#'
#' Computes the correlation between a pseudobulk of the query object and the
#' reference dataset. The feature set is the intersection of the reference
#' variable features and all features present in the query. Correlation is
#' computed on log normalized expression values using spearman correlation.
#'
#' @param object Query object
#' @param ref Reference expression averages
#' @param min.features If fewer than min.features exist in the intersection,
#' return 0 as correlation.
#'
#' @return Returns a list with the following values
#' \describe{
#'  \item{\code{cor.res}}{correlation between query and reference}
#'  \item{\code{plot}}{
#'   scatterplot of average expression of query and reference common genes
#'  }
#' }
#'
#' @importFrom stats cor
#' @importFrom Seurat AverageExpression Idents<- NormalizeData
#' @importFrom ggplot2 ggplot aes_string geom_point xlab ylab ggtitle theme_bw
#'
#' @keywords internal
#'
PBCorTest <- function(object, ref, min.features = 250) {
  Idents(object = object) <- "PBTest"
  features <- intersect(rownames(x = object), rownames(x = ref))
  features <- features[features %in% rownames(x = object) & features %in% rownames(x = ref)]
  if (length(x = features) < 250) {
    return(0)
  }
  object <- NormalizeData(object = object, assay = "RNA", verbose = FALSE)
  avg <- AverageExpression(
    object = object,
    features = rownames(x = object),
    assays = "RNA",
    verbose = FALSE
  )[[1]]
  cor.res <- cor(
    x = log1p(avg[features, ]), y = log1p(ref[features, ]),
    method = "spearman"
  )
  cor.data <- data.frame(log1p(x = avg[features, ]), log1p(x = ref[features, ]))
  colnames(cor.data) <- c("q", "r")
  plot <- ggplot(cor.data, aes_string(x = "q", y = "r")) +
    geom_point() +
    xlab("Query average expression") +
    ylab("Reference average expression") +
    ggtitle(label = paste(
      "Pseudo-bulk correlation of",
      length(x = features),
      "common genes:",
      round(x = cor.res, digits = 2)
    )) +
    theme_bw()
  return(list(cor.res = cor.res, plot = plot))
}


# Post mapping QC metric
#
# @param query Query object
# @param ds.amount Amount to downsample query
# @return Returns
#
#' @importFrom SeuratObject Cells Idents Indices as.Neighbor
#' @importFrom Seurat RunPCA FindNeighbors FindClusters
#
#' @keywords internal
#
#
MappingQCMetric <- function(query, ds.amount = 5000){
  if (ncol(x = query) > ds.amount) {
    query <- subset(x = query, cells = sample(x = Cells(x = query), size = ds.amount))
  }
  query <- RunPCA(object = query, verbose = FALSE)
  query <- FindNeighbors(object = query, reduction = 'pca', dims = 1:30, graph.name = 'pca_snn')
  query[["orig_neighbors"]] <- as.Neighbor(x = query[["pca_snn"]])
  query <- FindClusters(object = query, resolution = 0.6, graph.name = 'pca_snn')
  query <- FindNeighbors(object = query, reduction = 'integrated_dr', dims = 1:50, return.neighbor = TRUE, graph.name = "integrated_neighbors")
  query <- RunUMAP(object = query, dims = 1:30, reduction = 'pca', reduction.key = 'qcumap1_', reduction.name = 'qc.orig.umap')
  query <- RunUMAP(object = query, dims = 1:30, reduction = 'integrated_dr', reduction.key = 'qcumap2_', reduction.name = 'qc.intdr.umap')
  ids <- Idents(object = query)
  proj_ent <- unlist(x = lapply(X = 1:length(x = Cells(x = query)), function(x) {
    neighbors <- Indices(object = query[["integrated_neighbors"]])[x, ]
    nn_ids <- ids[neighbors]
    p_x <- prop.table(x = table(nn_ids))
    nn_entropy <- sum(p_x * log(x = p_x), na.rm = TRUE)
    return(nn_entropy)
  }))
  names(x = proj_ent) <- Cells(x = query)
  orig_ent <- unlist(x = lapply(X = 1:length(x = Cells(x = query)), function(x) {
    neighbors <- Indices(object = query[["orig_neighbors"]])[x, ]
    nn_ids <- ids[neighbors]
    p_x <- prop.table(x = table(nn_ids))
    nn_entropy <- sum(p_x * log(x = p_x), na.rm = TRUE)
    return(nn_entropy)
  }))
  names(x = orig_ent) <- Cells(x = query)
  stat <- median(
    x = tapply(X = orig_ent, INDEX = ids, FUN = mean) -
      tapply(X = proj_ent, INDEX = ids, FUN = mean)
  )
  # TODO: slim down query
  query <- DietSeurat(
    object = query, assays = "refAssay",
    dimreducs =  c("qc.orig.umap", "qc.intdr.umap"), counts = FALSE
  )
  Misc(object = query, slot = "mapping.qc") <- list(
    stat = stat,
    proj_ent = proj_ent,
    orig_ent = orig_ent
  )
  return(query)
}

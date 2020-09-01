#' @include zzz.R
#' @importFrom Seurat Embeddings
#'
NULL

#' Add Predictions
#'
#' @param object A \code{Seurat} object
#' @param preds A named vector with the predicted IDs
#' @param scores A named vector with the scores for \code{preds}
#' @param preds.levels Levels for predicted IDs, useful for
#' ordering \code{preds}
#' @param preds.drop Drop unused levels from \code{preds}
#'
#' @return \code{object} with transformations applied
#'
#' @importFrom Seurat Idents<-
#'
#' @keywords internal
#'
AddPredictions <- function(
  object,
  preds,
  scores,
  preds.levels = unique(x = as.character(x = preds)),
  preds.drop = TRUE
) {
  on.exit(expr = gc(verbose = FALSE))
  cells <- colnames(x = object)
  preds <- factor(x = preds[cells], levels = preds.levels)
  if (isTRUE(x = preds.drop)) {
    preds <- droplevels(x = preds)
  }
  object$predicted.id <- preds
  object$predicted.id.score <- scores[cells]
  Idents(object = object) <- 'predicted.id'
  # object[[reduction.name]] <- RunUMAP(
  #   object = neighbors,
  #   reduction.model = reduction.model,
  #   reduction.key = reduction.key
  # )
  return(object)
}

#' Calculate a mapping metric
#'
#' @inheritParams QueryReference
#' @param object A \code{Seurat} object with reference and query cells; easily
#' generated with \code{\link{QueryReference}}
#' @param reduction Name of reduction to use
#'
#' @return Returns a data frame with the following columns:
#' \describe{
#'  \item{\code{mapping.score}}{Mapping score}
#'  \item{\code{mapping.score.smoothed}}{Smoothed mapping score}
#'  \item{\code{r.metric}}{R metric}
#'  \item{\code{mapped}}{
#'   A logical indiciating if the cell has been mapped or not
#'  }
#'  This data frame has the same number of cells marked as \dQuote{query} in
#'  \code{object$ingest}; rownames are the names of query cells
#' }
#'
#' @importFrom rdist cdist
#'
#' @keywords internal
#'
CalcMappingMetric <- function(object, reduction = 'int', dims = 1:50) {
  on.exit(expr = gc(verbose = FALSE))
  embeddings <- Embeddings(object = object[[reduction]])[, dims]
  # Get intial mapping score
  # TODO: export AnnoyNN
  # TODO: export L2Norm
  all.nn <- Seurat:::AnnoyNN(
    data = Seurat:::L2Norm(mat = embeddings),
    metric = 'euclidean',
    n.trees = 10,
    k = 21
  )
  query.cells <- which(x = object$ingest == 'query')
  query.thresh <- min(query.cells)
  query.nn <- all.nn$nn.idx[query.cells, ]
  query.nn[query.nn < query.thresh] <- 1
  query.nn[query.nn >= query.thresh] <- 0
  map.score <- rowSums(x = query.nn)
  # Smooth mapping score
  query.l2 <- Seurat:::L2Norm(mat = embeddings[query.cells, ])
  query.nn <- Seurat:::AnnoyNN(
    data = query.l2,
    metric = 'euclidean',
    n.trees = 10,
    k = 21
  )
  map.smoothed <- sapply(
    X = seq_along(along.with = map.score),
    FUN = function(i) {
      return(mean(x = map.score[query.nn$nn.idx[i, ]]))
    }
  )
  # Calcualte predicted from reference
  ref.l2 <- Seurat:::L2Norm(mat = embeddings[-query.cells, ])
  query.ref.nn <- Seurat:::AnnoyNN(
    data = ref.l2,
    query = query.l2,
    metric = 'euclidean',
    n.trees = 10,
    k = 21
  )
  query.predicted <- t(x = sapply(
    X = seq_along(along.with = query.cells),
    FUN = function(i) {
      return(colMeans(x = ref.l2[query.ref.nn$nn.idx[i, -1], ]))
    }
  ))
  query.dist <- sapply(
    X = seq_along(along.with = query.cells),
    FUN = function(i) {
      return(cdist(
        X = query.l2[i, , drop = FALSE],
        Y = query.predicted[i, , drop = FALSE]
      ))
    }
  )
  query.dist <- query.dist - query.nn$nn.dists[, 2]
  k.dist <- query.nn$nn.dists[, ncol(x = query.nn$nn.dists)] - query.nn$nn.dists[, 2]
  r.metric <- query.dist / k.dist
  mapping.metric <- data.frame(
    mapping.score = map.score,
    mapping.score.smoothed = map.smoothed,
    r.metric = r.metric,
    mapped = ifelse(
      test = map.score < 2 & r.metric > 1,
      yes = FALSE,
      no = TRUE
    ),
    row.names = names(x = query.cells)
  )
  return(mapping.metric)
}

#' Find a minimum number of matching cells
#'
#' @inheritParams base::set.seed
#' @param reference Reference \code{\link[Seurat]{Seurat}} object
#' @param query Query \code{\link[Seurat]{Seurat}} object
#' @param match Ensure that we only use identity classes present in both objects
#'
#' @return A vector of cells in \code{reference} in the identites of
#' \code{query}; for each identity, there will be either the same number of
#' cells for that identity as in \code{reference} or in \code{query}, whichever
#' is smaller
#'
#' @importFrom Seurat Idents CellsByIdentities
#'
#' @keywords internal
#'
MinimalMatchingCells <- function(reference, query, match = TRUE, seed = NULL) {
  on.exit(expr = gc(verbose = FALSE))
  set.seed(seed = seed)
  query.table <- table(Idents(object = query))
  if (isTRUE(x = match)) {
    query.table <- query.table[names(x = query.table) %in% levels(x = reference)]
  }
  if (!length(x = query.table)) {
    stop("No matching identities between the reference and query objects")
  }
  cbi <- CellsByIdentities(object = reference, idents = names(x = query.table))
  cells <- sapply(
    X = names(x = query.table),
    FUN = function(ident) {
      cells.id <- cbi[[ident]]
      if (length(x = cells.id)) {
        cells.id <- sample(
          x = cells.id,
          size = min(length(x = cells.id), query.table[[ident]])
        )
      }
      return(cells.id)
    }
  )
  return(cells)
}

#' Transform an NN index
#'
#' @param neighbors Neighbors
#' @param meta.data Metadata
#' @param key Column of metadata to use
#'
#' @return \code{neighbors} transfomed
#'
#' @keywords internal
#'
NNTransform <- function(neighbors, meta.data, key = 'ori.index') {
  on.exit(expr = gc(verbose = FALSE))
  ori.index <- t(x = sapply(
    X = 1:nrow(x = neighbors$nn.idx),
    FUN = function(i) {
      return(meta.data[neighbors$nn.idx[i, ], key])
    }
  ))
  rownames(x = ori.index) <- rownames(x = neighbors$nn.idx)
  neighbors$nn.idx <- ori.index
  return(neighbors)
}

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
#' @importFrom Seurat AverageExpression Idents<-
#'
#' @return Returns correlation between query and reference

PBCorTest <- function(object, ref, min.features = 250) {
  Idents(object = object) <- "PBTest"
  features <- intersect(rownames(x = object), rownames(x = ref))
  features <- features[features %in% rownames(x = object) & features %in% rownames(x = ref)]
  if(length(x = features) < 250) {
    return(0)
  }
  avg <- AverageExpression(object = object, features = rownames(object), assays = "RNA", verbose = FALSE)[[1]]
  cor.res <- cor(x = log1p(avg[features, ]), y = log1p(ref[features, ]), method = "spearman")
  return(cor.res)
}

#' Merge reference and query objects together
#'
# @inheritParams base::set.seed
#' @inheritParams Seurat::RunUMAP
#' @inheritParams MinimalMatchingCells
#' @param assay.query Name of assay in \code{query} to use
#' @param reduction.reference Name of reduction in \code{reference} to use
#' @param reduction.query Name of reduction in \code{query} to use
#' @param dims Which dimensions to use
#' @inheritDotParams MinimalMatchingCells
#'
#' @return A \code{\link[Seurat]{Seurat}} object with the merged reference and
#' query datasets
#'
#' @importFrom utils head
#' @importFrom Seurat DietSeurat DefaultAssay<- CreateDimReducObject
#' DefaultAssay Idents
#'
#' @keywords internal
#'
QueryReference <- function(
  reference,
  query,
  assay.query = 'RNA',
  reduction.reference = 'spca',
  reduction.query = 'int',
  reduction.name = 'int',
  dims = 1:50,
  ...
) {
  # Pull the cell embeddings
  reduction.reference <- Embeddings(object = reference[[reduction.reference]])[, dims]
  reduction.query <- Embeddings(object = query[[reduction.query]])[, dims]
  cells <- unlist(
    x = MinimalMatchingCells(reference = reference, query = query, ...),
    use.names = FALSE
  )
  # Diet the reference object
  reference <- subset(x = reference, cells = cells)
  features <- head(
    x = intersect(x = rownames(x = reference), y = rownames(x = query)),
    n = 20L
  )
  reference <- DietSeurat(
    object = reference,
    counts = FALSE,
    data = TRUE,
    features = features,
    assays = 'RNA'
  )
  reference$ingest <- 'reference'
  # Diet the query object
  DefaultAssay(object = query) <- assay.query
  query <- DietSeurat(
    object = query,
    counts = FALSE,
    data = TRUE,
    features = features,
    assays = assay.query
  )
  query$ingest <- 'query'
  # Merge the objects
  object <- merge(x = reference, y = query)
  merged.embed <- rbind(
    reduction.reference[colnames(x = reference), ],
    reduction.query
  )
  object[[reduction.name]] <- suppressWarnings(expr = CreateDimReducObject(
    embeddings = merged.embed,
    assay = DefaultAssay(object = object),
    key = reduction.name
  ))
  object$merged.id <- Idents(object = object)
  return(object)
}

#' @useDynLib SeuratMapper
#' @include zzz.R
#' @importFrom Seurat Embeddings
#'
NULL

#' Add Predictions
#'
#' @param object A \code{Seurat} object
#' @param preds dataframe: Rows are cells. Columns are predicted.id,
#' predicted.id.score, and prediction.score. [CLASSNAME] for all reference
#' classes
#' @param preds.levels Levels for predicted IDs, useful for
#' ordering \code{preds}
#' @param preds.drop Drop unused levels from \code{preds} and
#' \dQuote{predictions} assay.
#'
#' @return \code{object} with predicted.id and predicted.id.score added to
#' metadata, and prediction.score. [CLASSNAME] added as an assay named
#' \dQuote{predictions}
#'
#' @importFrom Seurat Idents<-
#'
#' @keywords internal
#'
AddPredictions <- function(
  object,
  preds,
  preds.levels = unique(x = as.character(x = preds)),
  preds.drop = TRUE
) {
  on.exit(expr = gc(verbose = FALSE))
  cells <- colnames(x = object)
  # make sure order is the same as in the object
  preds <- preds[cells, , drop = FALSE]
  object$predicted.id.score <- preds$predicted.id.score
  ids <- factor(x = preds$predicted.id, levels = preds.levels)
  if (isTRUE(x = preds.drop)) {
    ids <- droplevels(x = ids)
    # Keep only prediction scores for levels not dropped;
    # get rid of columns already added as metadata in the process
    # spaces in levels in predicted.id are replaced with dots in the column names
    preds <- subset(
      x = preds,
      select = gsub(
        pattern = " ",
        replacement = ".",
        x = paste0("prediction.score.", levels(x = ids))
      )
    )
  } else {
    # Keep all prediction scores, just get rid of columns already added as metadata
    preds <- subset(x = preds, select = -c("predicted.id", "predicted.id.score"))
  }
  object$predicted.id <- ids
  Idents(object = object) <- 'predicted.id'
  object[["predictions"]] <- CreateAssayObject(data = t(preds))
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

#' Metric for evaluating transfer quality
#'
#' @param anchorset Anchorset object returned from FindTransferAnchors
#' @param ref Reference object
#' @param query Query object
#' @param k Number of anchors to use in projection steps when computing weights
#' @param ndim Number of dimensions to use when working with low dimensional
#' projections of the data
#' @param ksmooth Number of cells to average over when computing transition
#' probabilities
#' @param ksnn Number of cells to average over when determining the kernel
#' bandwidth from the SNN graph
#' @param snn.prune Amount of pruning to apply to edges in SNN graph
#' @param verbose Display messages/progress
#'
#' @return Returns a vector of cell scores
#' @importFrom Seurat GetIntegrationData Embeddings Reductions
#' CreateDimReducObject Cells Distances Indices Index
#'
#' @export
#'
MappingScore <- function(
  anchorset,
  ref,
  query,
  ref.reduction = "spca",
  query.reduction = "spca",
  k = 50,
  ndim = 50,
  ksmooth = 100,
  ksnn = 20,
  snn.prune = 0,
  subtract.first.nn = TRUE,
  approx = FALSE,
  verbose = TRUE,
  debug = FALSE,
  query.weights = NULL
) {
  # Input checks
  start.time <- Sys.time()
  if (!query.reduction %in% Reductions(object = query)) {
    stop("Please provide a query with", query.reduction, "precomputed.")
  }
  if (!ref.reduction %in% Reductions(object = ref)) {
    stop("Please provide a reference with", ref.reduction, "precomputed.")
  }
  # Extract info from anchorset object
  combined.object <- slot(object = anchorset, name = "object.list")[[1]]
  anchors <- slot(object = anchorset, name = "anchors")
  reference.cells <- slot(object = anchorset, name = "reference.cells")
  query.cells <- slot(object = anchorset, name = "query.cells")
  query.neighbors <- slot(object = anchorset, name = "neighbors")[["query.neighbors"]]
  # Project reference values onto query
  if (verbose) {
    message("Projecting reference PCA onto query")
  }
  ## Need to set up an IntegrationData object to use FindWeights here
  int.mat <- matrix(data = NA, nrow = nrow(x = anchors), ncol = 0)
  rownames(x = int.mat) <- query.cells[anchors[, "cell2"]]
  slot(object = combined.object, name = 'tools')[["IT1"]] <- new(
    Class = "IntegrationData",
    anchors = anchors,
    neighbors = list(cells1 = reference.cells, cells2 = query.cells),
    integration.matrix = int.mat
  )
  ## Finding weights of anchors in query pca space
  ref.pca.orig <- Embeddings(object = ref[[ref.reduction]])[, 1:ndim]
  rownames(x = ref.pca.orig) <- paste0(rownames(x = ref.pca.orig), "_reference")
  query.pca.orig <- Embeddings(object = query[[query.reduction]])[, 1:ndim]
  rownames(x = query.pca.orig) <- paste0(rownames(x = query.pca.orig), "_query")

  dr.weights <- suppressWarnings(CreateDimReducObject(
    embeddings = rbind(query.pca.orig, ref.pca.orig)
  ))
  if (!is.null(x = query.weights)) {
    weights.matrix <- query.weights
  } else {
    combined.object <- Seurat:::FindWeights(
      object = combined.object,
      integration.name = "IT1",
      reduction = dr.weights,
      dims = 1:ncol(x = dr.weights),
      k = k,
      sd.weight = 1,
      eps = 0,
      nn.method = "annoy",
      cpp = TRUE,
      verbose = verbose
    )
    weights.matrix <- GetIntegrationData(
      object = combined.object,
      integration.name = "IT1",
      slot = "weights"
    )
  }
  ## Perform projection of ref pca values using weights matrix
  ref.pca <- Embeddings(object = ref[[ref.reduction]])[Cells(x = ref)[anchors[, 1]], 1:ndim]
  rownames(x = ref.pca) <- paste0(rownames(x = ref.pca), "_reference")
  query.cells.projected <- crossprod(x = ref.pca, y = as.matrix(x = weights.matrix))
  colnames(x = query.cells.projected) <- query.cells
  rownames(x = query.cells.projected) <- colnames(x = ref.pca)

  # Re-project the query cells back onto query
  if (verbose) {
    message("Projecting back the query cells into original PCA space")
  }
  ## Compute new weights
  dr.weights <- suppressWarnings(CreateDimReducObject(
    embeddings = rbind(
      t(x = as.matrix(x = query.cells.projected)),
      ref.pca.orig[reference.cells, ]
    ),
  ))
  combined.object <- Seurat:::FindWeights(
    object = combined.object,
    integration.name = "IT1",
    reduction = dr.weights,
    dims = 1:ndim,
    k = k,
    sd.weight = 1,
    eps = 0,
    nn.method = "annoy",
    reverse = TRUE,
    cpp = TRUE,
    verbose = verbose
  )
  weights.matrix <- GetIntegrationData(combined.object, integration.name = "IT1", slot = "weights")
  ## Project back onto query
  orig.pca <- Embeddings(object = query[[query.reduction]])[Cells(x = query)[anchors[, 2]], ]
  query.cells.back.corrected <- Matrix::t(x = Matrix::crossprod(x = as(object = orig.pca, Class = "dgCMatrix"), y = weights.matrix)[1:ndim, ])
  rownames(x = query.cells.back.corrected) <- query.cells
  query.cells.orig <- gsub(pattern = "_query", replacement = "", x = query.cells)
  query.cells.pca <- Embeddings(object = query[[query.reduction]])[query.cells.orig, 1:ndim]
  if (verbose) {
    message("Computing scores:")
    message("    Finding neighbors of original query cells")
  }
  ## Compute original neighborhood of query cells
  nn.method <- "rann"
  if (approx) {
    nn.method <- "annoy"
  }
  if (is.null(x = query.neighbors)) {
    query.neighbors <- Seurat:::NNHelper(
      data = query.cells.pca,
      query = query.cells.pca,
      k = max(ksmooth, ksnn),
      method = nn.method,
      cache.index = TRUE
    )
  }
  if (verbose) message("    Finding neighbors of transformed query cells")
  ## Compute new neighborhood of query cells after projections
  if (nn.method == "annoy") {
    if (is.null(x = Index(object = query.neighbors))) {
      corrected.neighbors <- Seurat:::NNHelper(
        data = query.cells.pca,
        query = query.cells.back.corrected,
        k = max(ksmooth, ksnn),
        method = nn.method,
        cache.index = TRUE
      )
    } else {
      corrected.neighbors <- Seurat:::AnnoySearch(
        index = Index(object = query.neighbors),
        query = query.cells.back.corrected,
        k = max(ksmooth, ksnn)
      )
      corrected.neighbors <- Seurat:::Neighbor(
        nn.idx = corrected.neighbors$nn.idx,
        nn.dist = corrected.neighbors$nn.dists
      )
    }
  }
  if (verbose) message("    Computing query SNN")
  snn <- Seurat:::ComputeSNN(nn_ranked = Indices(query.neighbors)[, 1:ksnn], prune = snn.prune)
  query.cells.pca <- t(x = query.cells.pca)
  if (verbose) message("    Determining bandwidth and computing transition probabilities")
  scores <- ScoreHelper(
    snn = snn,
    query_pca = query.cells.pca,
    query_dists = Distances(query.neighbors),
    corrected_nns = Indices(corrected.neighbors),
    k_snn = ksnn,
    subtract_first_nn = subtract.first.nn,
    display_progress = verbose
  )
  scores[scores > 1] <- 1
  query.names <- gsub(pattern = "_query", replacement = "", x = query.cells)
  names(x = scores) <- query.names
  end.time <- Sys.time()
  if (verbose) {
    message("Total elapsed time: ", end.time - start.time)
  }
  return(scores)
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
#' @param object Seurat object
#' @param meta.data Metadata
#' @param neighbor.slot Name of Neighbor slot
#' @param key Column of metadata to use
#'
#' @return \code{object} with transfomed neighbor.slot
#'
#' @importFrom Seurat Indices
#'
#' @keywords internal
#'
NNTransform <- function(
  object,
  meta.data,
  neighbor.slot = "query_ref.nn",
  key = 'ori.index'
) {
  on.exit(expr = gc(verbose = FALSE))
  ind <- Indices(object[[neighbor.slot]])
  ori.index <- t(x = sapply(
    X = 1:nrow(x = ind),
    FUN = function(i) {
      return(meta.data[ind[i, ], key])
    }
  ))
  rownames(x = ori.index) <- rownames(x = ind)
  slot(object = object[[neighbor.slot]], name = "nn.idx") <- ori.index
  return(object)
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
#' @return Returns a list with the following values
#' \describe{
#'  \item{\code{cor.res}}{correlation between query and reference}
#'  \item{\code{plot}}{
#'   scatterplot of average expression of query and reference common genes
#'  }
#' }
#'
#' @importFrom stats cor
#' @importFrom Seurat AverageExpression Idents<-
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab ggtitle theme_bw
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
  plot <- ggplot(cor.data, aes(q, r)) +
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

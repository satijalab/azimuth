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

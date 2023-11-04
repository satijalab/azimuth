#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Run Azimuth annotation
#'
#' @param query Seurat object or following type of path:
#' \itemize{
#'  \item A \code{.h5} matrix
#'  \item A \code{.rds} file containing a Seurat object
#'  \item A \code{.h5ad} anndata object
#'  \item A \code{.h5seurat} object
#' }
#' @return Returns a Seurat object containing celltype annotations
#'
#' @rdname RunAzimuth
#' @export RunAzimuth
#'
RunAzimuth <- function(query, ...) {
  UseMethod(generic = 'RunAzimuth', object = query)
}


#' Run Azimuth annotation for ATAC query
#'
#' @param query Seurat object or following type of path:
#' \itemize{
#'  \item A \code{.h5} matrix
#'  \item A \code{.rds} file containing a Seurat object
#'  \item A \code{.h5ad} anndata object
#'  \item A \code{.h5seurat} object
#' }
#' @return Returns a Seurat object containing celltype annotations
#'
#' @rdname RunAzimuthATAC
#' @export RunAzimuthATAC
#'
RunAzimuthATAC <- function(query, ...) {
  UseMethod(generic = 'RunAzimuthATAC', object = query)
}

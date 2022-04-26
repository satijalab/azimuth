#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Run Azimuth annotation
#'
#' @param query Seurat object or H5AD path
#' @return Returns a Seurat object containing celltype annotations
#'
#' @rdname RunAzimuth
#' @export RunAzimuth
#'
RunAzimuth <- function(query, ...) {
  UseMethod(generic = 'RunAzimuth', object = query)
}

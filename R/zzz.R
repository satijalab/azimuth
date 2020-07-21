#' @docType package
#' @name SeuratMapper-pacakge
#' @rdname SeuratMapper-pacakge
#'
#' @aliases SeuratMapper
#'
"_PACKAGE"

#' Create an annoy index
#'
#' @note Function exists because it's not exported from \pkg{uwot}
#'
#' @param name,ndim Paramters
#'
#' @return An nn index object
#'
#' @importFrom methods new
#' @importFrom RcppAnnoy AnnoyAngular AnnoyManhattan AnnoyEuclidean AnnoyHamming
#'
#' @keywords internal
#'
CreateAnn <- function(name, ndim) {
  return(switch(
    EXPR = name,
    cosine = new(Class = AnnoyAngular, ndim),
    manhattan = new(Class = AnnoyManhattan, ndim),
    euclidean = new(Class = AnnoyEuclidean, ndim),
    hamming = new(Class = AnnoyHamming, ndim),
    stop("BUG: unknown Annoy metric '", name, "'")
  ))
}

#' Sanitize feature names for \code{\link[shiny]{selectInput}}
#'
#' \code{\link[shiny]{selectInput}} has some limitiations with biological
#' feature names. This function sanitizes feature names according to the
#' following rules:
#' \itemize{
#'  \item Names matching \dQuote{\\.1} are \strong{removed}
#' }
#'
#' @param features A character vector of feature names
#'
#' @return \code{features}, but sanitized
#'
#' @keywords internal
#'
#' @seealso \code{\link[shiny]{selectInput}}
#'
FilterFeatures <- function(features) {
  return(sort(x = grep(
    pattern = '\\.1',
    x = features,
    value = TRUE,
    invert = TRUE
  )))
}

#' Load file input into a \code{Seurat} object
#'
#' @param path Path to input data
#'
#' @return A \code{Seurat} object
#'
#' @importFrom tools file_ext
#' @importFrom SeuratDisk LoadH5Seurat
#' @importFrom Seurat Read10X_h5 CreateSeuratObject
#'
#' @keywords internal
#'
LoadFileInput <- function(path) {
  type <- tolower(x = file_ext(x = path))
  return(switch(
    EXPR = type,
    'h5' = {
      mat <- Read10X_h5(filename = path)
      if (is.list(x = mat)) {
        mat <- mat[[1]]
      }
      CreateSeuratObject(counts = mat)
    },
    'rds' = readRDS(file = path),
    'h5seurat' = LoadH5Seurat(file = path, assays = 'counts'),
    stop("Unknown file type: ", type, call. = FALSE)
  ))
}

#' Load the reference RDS files
#'
#' Read in a reference \code{\link[Seurat]{Seurat}} object and annoy index. This
#' function can read either from URLs or a file path. In order to read properly,
#' there must be the following files:
#' \itemize{
#'  \item \dQuote{ref.Rds} for the reference \code{Seurat} object
#'  \item \dQuote{idx.Rds} for the R-based index object
#'  \item \dQuote{idx.annoy} for the
#'  \link{RcppAnnoy::AnnoyIndex}{RcppAnnoy-based} index object
#' }
#'
#' @param path Path or URL to the two RDS files
#' @param seconds Timeout to check for URLs in seconds
#'
#' @return A list with two entries:
#' \describe{
#'  \item{\code{reference}}{The reference \code{\link[Seurat]{Seurat}} object}
#'  \item{\code{index}}{
#'    A list with index information, includes an
#'    \code{\link[RcppAnnoy]{AnnoyIndex}} object
#'  }
#' }
#'
#' @importFrom httr build_url parse_url modify_url status_code GET timeout
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Load from a URL
#' ref <- LoadReference("http://saucyx220.nygenome.org")
#' # Load from a directory
#' ref2 <- LoadReference("/var/www/html")
#' }
#'
LoadReference <- function(path, seconds = 10L) {
  ref.name <- 'ref.Rds'
  idx.name <- 'idx.Rds'
  ann.name <- 'idx.annoy'
  uri <- build_url(url = parse_url(url = path))
  if (grepl(pattern = '^://', x = uri)) {
    if (!dir.exists(paths = path)) {
      stop("Cannot find directory ", path, call. = FALSE)
    }
    reference <- file.path(path, ref.name)
    index <- file.path(path, idx.name)
    annoy <- file.path(path, ann.name)
    if (!file.exists(reference)) {
      stop("Cannot find reference RDS in ", path, call. = FALSE)
    } else if (!file.exists(index)) {
      stop("Cannot find index RDS in ", path, call. = FALSE)
    } else if (!file.exists(annoy)) {
      stop("Cannot find index annoy in ", path, call. = FALSE)
    }
  } else {
    Online <- function(url) {
      return(identical(
        x = status_code(x = GET(url = url, timeout(seconds = seconds))),
        y = 200L
      ))
    }
    uri <- parse_url(url = uri)
    ref.uri <- modify_url(url = uri, path = ref.name)
    idx.uri <- modify_url(url = uri, path = idx.name)
    ann.uri <- modify_url(url = uri, path = ann.name)
    if (!Online(url = ref.uri)) {
      stop("Cannot find reference RDS object at ", path, call. = FALSE)
    } else if (!Online(url = idx.uri)) {
      stop("Cannot find index RDS object at ", path, call = FALSE)
    } else if (!Online(url = ann.uri)) {
      stop("Cannot find index annoy object at ", path, call. = FALSE)
    }
    reference <- url(description = ref.uri)
    index <- url(description = idx.uri)
    annoy <- tempfile()
    download.file(url = ann.uri, destfile = annoy, quiet = TRUE)
    on.exit(expr = {
      close(con = reference)
      close(con = index)
      unlink(x = annoy)
    })
  }
  nn <- readRDS(file = index)
  annoy.index <- CreateAnn(name = nn$metric, ndim = nn$ndim)
  annoy.index$load(file = annoy)
  nn$annoy_index <- annoy.index
  return(list(
    reference = readRDS(file = reference),
    index = nn
  ))
}

GetDefaultArguments <- function(f, method = NULL) {
  UseMethod(generic = 'GetDefaultArguments', object = f)
}

#' @method GetDefaultArguments character
#'
GetDefaultArguments.character <- function(f, method = NULL) {
  ''
}

#' @method GetDefaultArguments function
#'
GetDefaultArguments.function <- function(f, method = NULL) {
  browser()
  return(GetDefaultArguments(
    f = as.character(x = substitute(expr = f)),
    method = method
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load Hooks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {
  op <- options()
  # TODO: replace this
  options(shiny.maxRequestSize = 100 * (1024 ^ 2))
}

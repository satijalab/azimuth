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
  # TODO: add support for H5AD and loom files
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
#'  \item \dQuote{ref.Rds} for the downsampled reference \code{Seurat}
#'  object (for mapping)
#'  \item \dQuote{plotref.Rds} for the reference \code{Seurat}
#'  object (for plotting)
#'  \item \dQuote{fullref.Rds} for the full reference \code{Seurat}
#'  object (for density estimation)
#'  \item \dQuote{adtref.Rds} for the ADT data for the reference \code{Seurat}
#'  object (for feature imputation)
#'  \item \dQuote{idx.Rds} for the nearest-neighbor index object
#' }
#'
#' @param path Path or URL to the two RDS files
#' @param seconds Timeout to check for URLs in seconds
#'
#' @return A list with four entries:
#' \describe{
#'  \item{\code{map}}{
#'    The downsampled reference \code{\link[Seurat]{Seurat}}
#'    object (for mapping)
#'  }
#'  \item{\code{plot}}{The reference \code{Seurat} object (for plotting)}
#'  \item{\code{full}}{
#'    The full reference \code{Seurat} object (for density estimation)
#'  }
#'  \item{\code{index}}{
#'    A list with nearest-neighbor index information, includes an
#'    \code{\link[RcppAnnoy]{AnnoyIndex}} object
#'  }
#' }
#'
#' @importFrom httr build_url parse_url status_code GET timeout
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
  ref.names <- list(
    map = 'ref.Rds',
    plt = 'plotref.Rds',
    ref = 'fullref.Rds',
    adt = 'adtref.Rds',
    idx = 'idx.Rds',
    ann = 'idx.annoy'
  )
  if (substr(x = path, start = nchar(x = path), stop = nchar(x = path)) == '/') {
    path <- substr(x = path, start = 1, stop = nchar(x = path) - 1)
  }
  uri <- build_url(url = parse_url(url = path))
  if (grepl(pattern = '^://', x = uri)) {
    if (!dir.exists(paths = path)) {
      stop("Cannot find directory ", path, call. = FALSE)
    }
    mapref <- file.path(path, ref.names$map)
    pltref <- file.path(path, ref.names$plt)
    fllref <- file.path(path, ref.names$ref)
    adtref <- file.path(path, ref.names$adt)
    idxref <- file.path(path, ref.names$idx)
    annref <- file.path(path, ref.names$ann)
    exists <- file.exists(c(mapref, pltref, fllref, adtref, idxref, annref))
    if (!all(exists)) {
      stop(
        "Missing the following files from the directory provided: ",
        Oxford(unlist(x = ref.names)[!exists], join = 'and')
      )
    }
  } else {
    Online <- function(url, strict = FALSE) {
      if (isTRUE(x = strict)) {
        code <- 200L
        comp <- identical
      } else {
        code <- 404L
        comp <- Negate(f = identical)
      }
      return(comp(
        x = status_code(x = GET(url = url, timeout(seconds = seconds))),
        y = code
      ))
    }
    ref.uris <- paste(uri, ref.names, sep = '/')
    names(x = ref.uris) <- names(x = ref.names)
    online <- vapply(
      X = ref.uris,
      FUN = Online,
      FUN.VALUE = logical(length = 1L),
      USE.NAMES = FALSE
    )
    if (!all(online)) {
      stop(
        "Cannot find the following files at the site given: ",
        Oxford(unlist(x = ref.names)[!online], join = 'and')
      )
    }
    mapref <- url(description = ref.uris[['map']])
    pltref <- url(description = ref.uris[['plt']])
    fllref <- url(description = ref.uris[['ref']])
    adtref <- url(description = ref.uris[['adt']])
    idxref <- url(description = ref.uris[['idx']])
    # annref <- url(description = ref.uris[6])
    annref <- tempfile()
    download.file(url = ref.uris[['ann']], destfile = annref, quiet = TRUE)
    on.exit(expr = {
      close(con = mapref)
      close(con = pltref)
      close(con = fllref)
      close(con = adtref)
      close(con = idxref)
      unlink(x = annref)
    })
  }
  # Load the map reference and ADT values for imputation
  map <- readRDS(file = mapref)
  map[['ADT']] <- readRDS(file = adtref)[['ADT']]
  # Load the annoy index
  nn <- readRDS(file = idxref)
  annoy.index <- CreateAnn(name = nn$metric, ndim = nn$ndim)
  annoy.index$load(file = annref)
  nn$annoy_index <- annoy.index
  gc(verbose = FALSE)
  return(list(
    map = map,
    plot = readRDS(file = pltref),
    full = readRDS(file = fllref),
    index = nn
  ))
}

#' Make An English List
#'
#' Joins items together to make an English list; uses the Oxford comma for the
#' last item in the list.
#'
#' @inheritParams base::paste
#' @param join either \dQuote{and} or \dQuote{or}
#'
#' @return A character vector of the values, joined together with commas and
#' \code{join}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' Oxford("red")
#' Oxford("red", "blue")
#' Oxford("red", "green", "blue")
#' }
#'
Oxford <- function(..., join = c('and', 'or')) {
  join <- match.arg(arg = join)
  args <- as.character(x = c(...))
  args <- Filter(f = nchar, x = args)
  if (length(x = args) == 1) {
    return(args)
  } else if (length(x = args) == 2) {
    return(paste(args, collapse = paste0(' ', join, ' ')))
  }
  return(paste0(
    paste(args[1:(length(x = args) - 1)], collapse = ', '),
    paste0(', ', join, ' '),
    args[length(x = args)]
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

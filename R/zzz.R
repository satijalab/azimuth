#' @docType package
#'
#' @section Package options:
#'
#' SeuratMapper uses the following options to control the behaviour of the app,
#' users can configure these with \code{\link[base]{options}}:
#'
#' \describe{
#'  \item{\code{Azimuth.de.mincells}}{
#'   Minimum number of cells per cluster for differential expression; defaults
#'   to \code{15}
#'  }
#'  \item{\code{Azimuth.map.pcthresh}}{
#'   Only show mapped plot if the percentage of cells mapped meets or
#'   exceeds this threshold; defaults to \code{60}
#'  }
#'  \item{\code{Azimuth.sct.ncells}, \code{Azimuth.sct.nfeats}}{
#'   Number of cells and features to use for
#'   \code{\link[Seurat]{SCTransform}}, respectively. Defaults to \code{1000}
#'   for each
#'  }
#' }
#'
#' @inheritSection AzimuthApp App options
#'
#'
#' @aliases SeuratMapper
#'
"_PACKAGE"

default.options <- list(
  Azimuth.de.mincells = 15L,
  Azimuth.map.pcthresh = 60L,
  Azimuth.sct.ncells = 1000L,
  Azimuth.sct.nfeats = 1000L
)

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
#' Take a file and load it into a \code{\link[Seurat]{Seurat}} object. Supports
#' a variety of file types and always returns a \code{Seurat} object
#'
#' @param path Path to input data
#'
#' @return A \code{\link[Seurat]{Seurat}} object
#'
#' @importFrom tools file_ext
#' @importFrom SeuratDisk Connect LoadH5Seurat
#' @importFrom Seurat Read10X_h5 CreateSeuratObject as.sparse
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
    'rds' = {
      object <- readRDS(file = path)
      if (inherits(x = object, what = c('Matrix', 'matrix', 'data.frame'))) {
        object <- CreateSeuratObject(counts = as.sparse(x = object))
      }
      if (!inherits(x = object, what = 'Seurat')) {
        stop("The RDS file must be a Seurat object", call. = )
      }
    },
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
#'   The downsampled reference \code{\link[Seurat]{Seurat}}
#'   object (for mapping)
#'  }
#'  \item{\code{plot}}{The reference \code{Seurat} object (for plotting)}
#'  \item{\code{full}}{
#'   The full reference \code{Seurat} object (for density estimation)
#'  }
#'  \item{\code{index}}{
#'   A list with nearest-neighbor index information, includes an
#'   \code{\link[RcppAnnoy]{AnnoyIndex}} object
#'  }
#' }
#'
#' @importFrom Seurat Idents<-
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
  # Load the other references
  plot <- readRDS(file = pltref)
  full <- readRDS(file = fllref)
  id.check <- vapply(
    X = c(map, plot, full),
    FUN = function(x) {
      return('id' %in% colnames(x = x[[]]))
    },
    FUN.VALUE = logical(length = 1L)
  )
  if (all(id.check)) {
    Idents(object = map) <-
      Idents(object = plot) <-
      Idents(object = full) <- 'id'
  }
  gc(verbose = FALSE)
  return(list(
    map = map,
    plot = plot,
    full = full,
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

#' Prepare differential expression results for rendering
#'
#' @param diff.exp A dataframe with differential expression results from
#' \code{\link[presto:wilcoxauc]{presto::wilcoxauc}}
#' @param groups.use Names of groups to filter \code{diff.exp} to; groups must
#' be found in \code{diff.exp$group}
#' @param n Number of feature to filter \code{diff.exp} to per group
#' @param logfc.thresh logFC threshold
#'
#' @return \code{diff.exp}, ordered by adjusted p-value, filtered to \code{n}
#' features per group in \code{group.use}
#'
#' @importFrom rlang %||%
#' @importFrom utils head
#'
#' @seealso \code{\link[presto]{wilcoxauc}}
#'
#' @keywords internal
#'
RenderDiffExp <- function(
  diff.exp,
  groups.use = NULL,
  n = 10L,
  logfc.thresh = 0L
) {
  cols.remove <- c('feature', 'logFC', 'auc')
  groups.use <- groups.use %||% unique(x = as.character(x = diff.exp$group))
  diff.exp <- lapply(
    X = groups.use,
    FUN = function(group) {
      group.de <- diff.exp[diff.exp$group == group, , drop = FALSE]
      group.de <- group.de[group.de$logFC > logfc.thresh, , drop = FALSE]
      group.de <- group.de[order(group.de$padj, -group.de$auc), , drop = FALSE]
      return(head(x = group.de, n = n))
    }
  )
  diff.exp <- do.call(what = 'rbind', diff.exp)
  rownames(x = diff.exp) <- make.unique(names = diff.exp$feature)
  diff.exp <- diff.exp[, !colnames(x = diff.exp) %in% cols.remove, drop = FALSE]
  return(diff.exp)
}

#' Manage tabs with \pkg{shinyjs}
#'
#' Quickly generate JavaScript IDs for Shiny tab panels. Also build a function
#' to hide Shiny tab panels in JavaScript using \pkg{shinyjs}
#'
#' @param id ID of a \code{\link[shiny]{tabsetPanel}}
#' @param values One or more values of a \code{\link[shiny]{tabPanel}}
#' (see the \code{value} parameter)
#' @param fxn Name of JavaScript call function
#'
#' @return \code{TabJSHide}: a string with a JavaScript function to hide a set
#' of tabs
#'
#' @name TabJS
#' @rdname TabJS
#'
#' @keywords internal
#'
#' @note These functions are designed to run custom JavaScript code using
#' \code{\link[shinyjs]{extendShinyJS}}; use of custom JavaScript code requires
#' the \pkg{V8} package. \pkg{V8} requires a local install of either the
#' \href{https://chromium.googlesource.com/v8/v8}{V8} JavaScript Engine or
#' \href{https://nodejs.org/en/}{Node.js}
#'
#' @seealso \code{\link[shinyjs]{extendShinyJS}} \code{\link[shiny]{tabPanel}}
#'
#' @examples
#' \dontrun{
#' TabJSHide('tabs', values = c('mapped', 'fexplorer'))
#' }
#'
TabJSHide <- function(id, values, fxn = 'hide') {
  keys <- paste0(
    '$(',
    sQuote(x = TabJSKey(id = id, values = values), q = FALSE),
    ').hide()',
    collapse = '; '
  )
  return(paste0(
    'shinyjs.',
    fxn,
    ' = function() {',
    keys,
    ';}'
  ))
}

#' @rdname TabJS
#'
#' @return \code{TabJSKey}: a string with the JavaScript ID for a given
#' set of tabs
#'
#' @examples
#' \dontrun{
#' TabJSKey('tabs', values = 'mapped')
#' }
#'

TabJSKey <- function(id, values) {
  return(paste0(
    "#",
    id,
    " li a[data-value=",
    dQuote(x = values, q = FALSE),
    "]"
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load Hooks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {
  op <- options()
  # TODO: replace this
  options(shiny.maxRequestSize = 100 * (1024 ^ 2))
  # Set some default options
  toset <- !names(x = default.options) %in% names(x = op)
  if (any(toset)) {
    options(default.options[toset])
  }
}

#' @docType package
#'
#' @section Package options:
#'
#' SeuratMapper uses the following options to control the behaviour of the app,
#' users can configure these with \code{\link[base]{options}}:
#'
#' \describe{
#' \item{\code{Azimuth.map.ncells}}{
#'   Minimum number of cells required to accept uploaded file.
#'   Defaults to \code{100}
#'  }
#'  \item{\code{Azimuth.map.ngenes}}{
#'   Minimum number of genes in common with reference to accept uploaded file.
#'   Defaults to \code{1000}
#'  }
#'  \item{\code{Azimuth.map.ngenes}}{
#'   Minimum number of anchors that must be found to complete mapping.
#'   Defaults to \code{50}
#'  }
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
  Azimuth.map.ncells = 100L,
  Azimuth.map.ngenes = 1000L,
  Azimuth.map.nanchors = 50L,
  Azimuth.map.pcthresh = 60L,
  Azimuth.sct.ncells = 1000L,
  Azimuth.sct.nfeats = 1000L
)

#' Attach dependent packages
#'
#' Attaches the following packages
#' \itemize{
#'  \item shinyBS
#' }
#'
#' @return Attaches the required packages and returns invisible \code{NULL}
#'
#' @keywords internal
#'
AttachDeps <- function() {
  deps <- c(
    'shinyBS'
  )
  for (d in deps) {
    packageStartupMessage("Attaching ", d)
    attachNamespace(ns = d)
  }
}

#' Returns a dataframe of the frequency or percentage of levels of category.2
#' (column) for object split by each level of category.1 (row)
#'
#' @param object a Seurat object
#' @param category.1 a metadata field in the object
#' @param category.2 another metadata field in the object
#' @param percentage if TRUE, returns percentages; otherwise, counts
#'
#' @importFrom Seurat FetchData
#'
#' @keywords internal
#'
CategoryTable <- function(
  object,
  category.1,
  category.2,
  percentage = FALSE
) {
  data <- FetchData(object, c(category.1, category.2))
  data[, category.1] <- droplevels(factor(x = data[, category.1]))
  data[, category.2] <- droplevels(factor(x = data[, category.2]))
  tbl <- table(data[, category.1], data[, category.2])
  if (percentage) {
    tbl <- t(apply(
      X = tbl,
      MARGIN = 1,
      FUN = function(x) round(100 * (x/sum(x)), digits = 1))
    )
  }
  return(as.data.frame.matrix(tbl))
}

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
#' @details \code{LoadFileInput} supports several file types to be read in as
#' \code{Seurat} objects. File type is determined by extension, matched in a
#' case-insensitive manner See sections below for details about supported
#' filtypes, required extension, and specifics for how data is loaded
#'
#' @param path Path to input data
#'
#' @return A \code{\link[Seurat]{Seurat}} object
#'
#' @importFrom tools file_ext
#' @importFrom SeuratDisk LoadH5Seurat
#' @importFrom Seurat Read10X_h5 CreateSeuratObject as.sparse DietSeurat
#' DefaultAssay
#'
#' @keywords internal
#'
#' @section 10X H5 File (extension \code{h5}):
#' 10X HDF5 files are supported for all versions of Cell Ranger; data is read
#' in using \code{\link[Seurat]{Read10X_h5}}. \strong{Note}: for multi-modal
#' 10X HDF5 files, only the \emph{first} matrix is read in
#'
#' @section Rds File (extension \code{rds}):
#' Rds files are supported as long as they contain one of the following data
#' types:
#' \itemize{
#'  \item A \code{\link[Seurat]{Seurat}} V3 object
#'  \item An S4 \code{\link[Matrix]{Matrix}} object
#'  \item An S3 \code{\link[base]{matrix}} object
#'  \item A \code{\link[base]{data.frame}} object
#' }
#' For S4 \code{Matrix}, S3 \code{matrix}, and \code{data.frame} objects, a
#' \code{Seurat} object will be made with
#' \code{\link[Seurat]{CreateSeuratObject}} using the default arguments
#'
#' @section h5Seurat File (extension \code{h5seurat}):
#' h5Seurat files and all of their features are fully supported. They are read
#' in via \code{\link[SeuratDisk]{LoadH5Seurat}}. \strong{Note}: only the
#' \dQuote{counts} matrices are read in and only the default assay is kept
#'
#' @inheritSection LoadH5AD AnnData H5AD File (extension \code{h5ad})
#'
LoadFileInput <- function(path) {
  # TODO: add support for loom files
  on.exit(expr = gc(verbose = FALSE))
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
      } else if (inherits(x = object, what = 'Seurat')) {
        object <- DietSeurat(
          object = object,
          assays = DefaultAssay(object = object)
        )
      } else {
        stop("The RDS file must be a Seurat object", call. = FALSE)
      }
      object
    },
    'h5ad' = LoadH5AD(path = path),
    'h5seurat' = {
      object <- LoadH5Seurat(file = path, assays = 'counts')
      object <- DietSeurat(
        object = object,
        assays = DefaultAssay(object = object)
      )
      object
    },
    stop("Unknown file type: ", type, call. = FALSE)
  ))
}

#' Load a diet H5AD file
#'
#' Read in only the counts matrix and (if present) metadata of an H5AD file and
#' return a \code{Seurat} object
#'
#' @inheritParams LoadFileInput
#'
#' @return A \code{Seurat} object
#'
#' @importFrom hdf5r h5attr
#' @importFrom SeuratDisk Connect
#' @importFrom Seurat AddMetaData CreateSeuratObject
#'
#' @keywords internal
#'
#' @section AnnData H5AD File (extension \code{h5ad}):
#' Only H5AD files from AnnData v0.7 or higher are supported. Data is read from
#' the H5AD file in the following manner
#' \itemize{
#'  \item The counts matrix is read from \dQuote{/raw/X}; if \dQuote{/raw/X} is
#'  not present, the matrix is read from \dQuote{/X}
#'  \item Feature names are read from feature-level metadata. Feature level
#'  metadata must be an HDF5 group, HDF5 compound datasets are \strong{not}
#'  supported. If counts are read from \code{/raw/X}, features names are looked
#'  for in \dQuote{/raw/var}; if counts are read from \dQuote{/X}, features
#'  names are looked for in \dQuote{/var}. In both cases, feature names are read
#'  from the dataset specified by the \dQuote{_index} attribute, \dQuote{_index}
#'  dataset, or \dQuote{index} dataset, in that order
#'  \item Cell names are read from cell-level metadata. Cell-level metadata must
#'  be an HDF5 group, HDF5 compound datasets are \strong{not} supported.
#'  Cell-level metadata is read from \dQuote{/obs}. Cell names are read from the
#'  dataset specified by the \dQuote{_index} attribute, \dQuote{_index} dataset,
#'  or \dQuote{index} dataset, in that order
#'  \item Cell-level metadata is read from the \dQuote{/obs} dataset. Columns
#'  will be returned in the same order as in the \dQuote{column-order}, if
#'  present, or in alphabetical order. If a dataset named \dQuote{__categories}
#'  is present, then all datasets in \dQuote{__categories} will serve as factor
#'  levels for datasets present in \dQuote{/obs} with the same name (eg. a
#'  dataset named \dQuote{/obs/__categories/leiden} will serve as the levels for
#'  \dQuote{/obs/leiden}). Row names will be set as cell names as described
#'  above. All datasets in \dQuote{/obs} will be loaded except for
#'  \dQuote{__categories} and the cell names dataset
#' }
#'
LoadH5AD <- function(path) {
  adata <- suppressWarnings(expr = Connect(filename = path, force = TRUE))
  on.exit(expr = adata$close_all())
  IsMetaData <- function(md) {
    check <- SeuratDisk:::Exists(x = adata, name = md)
    if (isTRUE(x = check)) {
      check <- inherits(x = adata[[md]], what = 'H5Group')
    }
    return(check)
  }
  GetIndex <- function(md) {
    return(
      if (adata[[md]]$attr_exists(attr_name = '_index')) {
        h5attr(x = adata[[md]], which = '_index')
      } else if (adata[[md]]$exists(name = '_index')) {
        '_index'
      } else if (adata[[md]]$exists(name = 'index')) {
        'index'
      } else {
        stop("Cannot find the rownames for '", md, "'", call. = FALSE)
      }
    )
  }
  GetRowNames <- function(md) {
    return(adata[[md]][[GetIndex(md = md)]][])
  }
  LoadMetadata <- function(md) {
    factor.cols <- if (adata[[md]]$exists(name = '__categories')) {
      names(x = adata[[md]][['__categories']])
    } else {
      NULL
    }
    index <- GetIndex(md = md)
    col.names <- names(x = adata[[md]])
    if (adata[[md]]$attr_exists(attr_name = 'column-order')) {
      tryCatch(
        expr = {
          col.order <- h5attr(x = adata[[md]], which = 'column-order')
          col.names <- c(
            intersect(x = col.order, y = col.names),
            setdiff(x = col.names, y = col.order)
          )
        },
        error = function(...) {
          return(invisible(x = NULL))
        }
      )
    }
    col.names <- col.names[!col.names %in% c('__categories', index)]
    df <- sapply(
      X = col.names,
      FUN = function(i) {
        x <- adata[[md]][[i]][]
        if (i %in% factor.cols) {
          x <- factor(x = x, levels = adata[[md]][['__categories']][[i]][])
        }
        return(x)
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    return(as.data.frame(x = df, row.names = GetRowNames(md = md)))
  }
  if (SeuratDisk:::Exists(x = adata, name = 'raw/X')) {
    md <- 'raw/var'
    counts <- as.matrix(x = adata[['raw/X']])
  } else if (SeuratDisk:::Exists(x = adata, name = 'X')) {
    md <- 'var'
    counts <- as.matrix(x = adata[['X']])
  } else {
    stop("Cannot find counts matrix")
  }
  if (!IsMetaData(md = md)) {
    stop("Cannot find feature-level metadata", call. = FALSE)
  } else if (!IsMetaData(md = 'obs')) {
    stop("Cannot find cell-level metadata", call. = FALSE)
  }
  metadata <- LoadMetadata(md = 'obs')
  rownames(x = counts) <- GetRowNames(md = md)
  colnames(x = counts) <- rownames(x = metadata)
  object <- CreateSeuratObject(counts = counts)
  if (ncol(x = metadata)) {
    object <- AddMetaData(object = object, metadata = metadata)
  }
  return(object)
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
#' @importFrom utils download.file
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

#' Return names of metadata columns in a Seurat object that have an
#' appropriate number of levels for plotting when converted to a factor
#'
#' @param object a Seurat object
#' @param min.levels minimum number of levels in a metadata factor to include
#' @param max.levels maximum number of levels in a metadata factor to include
#'
#' @keywords internal
#'
PlottableMetadataNames <- function(
  object,
  min.levels = 2,
  max.levels = 20
) {
  column.status <- sapply(
    object[[]],
    FUN = function(column) {
      length(levels(droplevels(as.factor(column)))) >= min.levels &&
        length(levels(droplevels(as.factor(column)))) <= max.levels
    }
  ) & (colnames(object[[]]) != "mapping.score") & (colnames(object[[]]) != "predicted.id")
  return(colnames(object[[]])[column.status])
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
  # cols.remove <- c('feature', 'logFC', 'auc')
  cols.keep <- c('avgExpr', 'auc', 'padj', 'pct_in', 'pct_out')
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
  # diff.exp <- diff.exp[, !colnames(x = diff.exp) %in% cols.remove, drop = FALSE]
  diff.exp <- diff.exp[, cols.keep, drop = FALSE]
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
  # Attach deps
  AttachDeps()
  op <- options()
  # TODO: replace this
  options(shiny.maxRequestSize = 100 * (1024 ^ 2))
  # Set some default options
  toset <- !names(x = default.options) %in% names(x = op)
  if (any(toset)) {
    options(default.options[toset])
  }
}

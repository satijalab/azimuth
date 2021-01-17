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
#' @importFrom SeuratDisk Connect
#' @importFrom Seurat Read10X_h5 CreateSeuratObject as.sparse Assays
#' GetAssayData DefaultAssay<- DietSeurat as.Seurat
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
  type <- tolower(x = tools::file_ext(x = path))
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
        if (!'RNA' %in% Assays(object = object)) {
          stop("No RNA assay provided", call. = FALSE)
        } else if (Seurat:::IsMatrixEmpty(x = GetAssayData(object = object, slot = 'counts', assay = 'RNA'))) {
          stop("No RNA counts matrix present", call. = )
        }
        DefaultAssay(object = object) <- "RNA"
        object <- DietSeurat(
          object = object,
          assays = "RNA"
        )
      } else {
        stop("The RDS file must be a Seurat object", call. = FALSE)
      }
      object
    },
    'h5ad' = LoadH5AD(path = path),
    'h5seurat' = {
      if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
        stop("Loading h5Seurat files requires SeuratDisk", call. = FALSE)
      }
      hfile <- suppressWarnings(expr = SeuratDisk::Connect(filename = path))
      on.exit(expr = hfile$close_all())
      if (!'RNA' %in% names(x = hfile[['assays']])) {
        stop("Cannot find the RNA assay in this h5Seurat file", call. = FALSE)
      } else if (!'counts' %in% names(x = hfile[['assays/RNA']])) {
        stop("No RNA counts matrix provided", call. = FALSE)
      }
      object <- as.Seurat(
        x = hfile,
        assays = list('RNA' = 'counts'),
        reductions = FALSE,
        graphs = FALSE,
        images = FALSE
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
#' @importFrom hdf5r H5File h5attr
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
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Loading H5AD files requires hdf5r", call. = FALSE)
  }
  adata <- hdf5r::H5File$new(filename = path, mode = 'r')
  on.exit(expr = adata$close_all())
  Exists <- function(name) {
    name <- unlist(x = strsplit(x = name[1], split = '/', fixed = TRUE))
    hpath <- character(length = 1L)
    exists <- TRUE
    for (i in seq_along(along.with = name)) {
      hpath <- paste(hpath, name[i], sep = '/')
      exists <- adata$exists(name = hpath)
      if (isFALSE(x = exists)) {
        break
      }
    }
    return(exists)
  }
  IsMetaData <- function(md) {
    return(Exists(name = md) && inherits(x = adata[[md]], what = 'H5Group'))
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
          col.order <- hdf5r::h5attr(x = adata[[md]], which = 'column-order')
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
  if (Exists(name = 'raw/X')) {
    md <- 'raw/var'
    counts <- as.matrix(x = adata[['raw/X']])
  } else if (Exists(name = 'X')) {
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
#'  \item \dQuote{idx.annoy} for the nearest-neighbor index object
#' }
#'
#' @param path Path or URL to the two RDS files
#' @param seconds Timeout to check for URLs in seconds
#'
#' @return A list with three entries:
#' \describe{
#'  \item{\code{map}}{
#'   The downsampled reference \code{\link[Seurat]{Seurat}}
#'   object (for mapping)
#'  }
#'  \item{\code{plot}}{The reference \code{Seurat} object (for plotting)}
#'  \item{\code{avgexp}}{Average expression (for pseudobulk check)}
#' }
#'
#' @importFrom Seurat Idents<- LoadAnnoyIndex
#' @importFrom httr build_url parse_url status_code GET timeout
#' @importFrom utils download.file
#' @importFrom Matrix sparseMatrix
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Load from a URL
#' ref <- LoadReference("https://seurat.nygenome.org/references/pbmc")
#' # Load from a directory
#' ref2 <- LoadReference("/var/www/html")
#' }
#'
LoadReference <- function(path, normalization.method, seconds = 10L) {
  ref.names <- list(
    map = 'ref.Rds',
    ann = 'idx.annoy'
  )
  if (substr(x = path, start = nchar(x = path), stop = nchar(x = path)) == '/') {
    path <- substr(x = path, start = 1, stop = nchar(x = path) - 1)
  }
  uri <- httr::build_url(url = httr::parse_url(url = path))
  if (grepl(pattern = '^://', x = uri)) {
    if (!dir.exists(paths = path)) {
      stop("Cannot find directory ", path, call. = FALSE)
    }
    mapref <- file.path(path, ref.names$map)
    annref <- file.path(path, ref.names$ann)
    exists <- file.exists(c(mapref, annref))
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
        x = httr::status_code(x = httr::GET(
          url = url,
          httr::timeout(seconds = seconds
        ))),
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
    annref <- tempfile()
    download.file(url = ref.uris[['ann']], destfile = annref, quiet = TRUE)
    on.exit(expr = {
      close(con = mapref)
      unlink(x = annref)
    })
  }
  # Load the map reference
  map <- readRDS(file = mapref)
  # Load the annoy index into the Neighbor object in the neighbors slot

  map[["refdr.annoy.neighbors"]] <- LoadAnnoyIndex(
    object = map[["refdr.annoy.neighbors"]],
    file = annref
  )
  if (normalization.method == 'SCT') {
    sct.model <- Misc(object = map[["SCT"]], slot = "vst.set")
    suppressWarnings(expr = Misc(object = map[["SCT"]], slot = "vst.set") <- list())
  }
  # Create plotref
  ad <- Tool(object = map, slot = "AzimuthReference")
  plotref.dr <- GetPlotRef(object = ad)
  cm <- sparseMatrix(
    i = 1, j = 1, x = 0, dims = c(1, nrow(x = plotref.dr)),
    dimnames = list("placeholder", Cells(x = plotref.dr))
  )
  plot <- CreateSeuratObject(
    counts = cm, names.field = NULL, names.delim = NULL # saves >10s on `feat/updateCSO` seurat branch
  )
  plot[["refUMAP"]] <- plotref.dr
  plot <- AddMetaData(object = plot, metadata = Misc(object = plotref.dr, slot = "plot.metadata"))
  avg <- GetAvgRef(object = ad)
  gc(verbose = FALSE)
  if (normalization.method == 'SCT') {
    return(list(
      map = map,
      sct.model = sct.model,
      plot = plot,
      avgexp = avg
    ))
  } else {
    return(list(
      map = map,
      plot = plot,
      avgexp = avg
    ))
  }
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

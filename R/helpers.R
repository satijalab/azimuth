#' Converts gene names of query to match type/species of reference names (human
#' or mouse).
#'
#' @param object Object to convert, must contain only RNA counts matrix
#' @param reference.names Gene names of reference
#' @param homolog.table Location of file (or URL) containing table with
#' human/mouse homologies
#' @return query object with converted feature names, likely subsetted
#'
#' @importFrom stringr str_match
#' @importFrom shiny isRunning
#' @useDynLib Azimuth
#' @export
#'
ConvertGeneNames <- function(object, reference.names, homolog.table) {
  uri <- httr::build_url(url = httr::parse_url(url = homolog.table))
  if (grepl(pattern = '^://', x = uri)) {
    if (!file.exists(homolog.table)) {
      stop("Homolog file doesn't exist at path provided: ", homolog.table)
    }
  } else {
    if (!Online(url = homolog.table)) {
      stop("Cannot find the homolog table at the URL given: ", homolog.table)
    }
    homolog.table <- url(description = homolog.table)
    on.exit(expr = {
      close(con = homolog.table)
    })
  }
  linked <- readRDS(file = homolog.table)
  query.names <- rownames(x = object)
  # remove version numbers
  ensembl.id <- '(?:ENSG|ENSMUS)'
  capture.nonversion <- paste0('(', ensembl.id, '.*)\\.(?:.*)')
  pattern <- paste0('(?:', capture.nonversion,')|(.*)')
  query.names <- Filter(
    f = function(x) !is.na(x = x),
    x = as.vector(x = t(x = str_match(string = query.names, pattern = pattern)[, 2:3]))
  )
  # determine idtype and species of query
  # 5000 because sometimes ensembl IDs are mixed with regular gene names
  query.names.sub = sample(
    x = query.names,
    size = min(length(x = query.names), 5000)
  )
  idtype <- names(x = which.max(x = apply(
    X = linked,
    MARGIN = 2,
    FUN = function(col) {
      length(x = intersect(x = col, y = query.names.sub))
    }
  )))
  species <- ifelse(test = length(x = grep(pattern = '\\.mouse', x = idtype)) > 0, 'mouse', 'human')
  idtype <- gsub(pattern = '\\.mouse|\\.human', replacement = '', x = idtype)
  message('detected inputs from ', toupper(x = species), ' with id type ', idtype)
  totalid <- paste0(idtype, '.', species)

  # determine idtype and species of ref
  reference.names.sub <- sample(x = reference.names, size = min(length(x = reference.names), 5000))
  idtype.ref <- names(x = which.max(x = apply(
    X = linked,
    MARGIN = 2,
    FUN = function(col) {
      length(x = intersect(x = col, y = reference.names))
    }
  )))
  species.ref <- ifelse(test = length(x = grep(pattern = '\\.mouse', x = idtype.ref)) > 0, 'mouse', 'human')
  idtype.ref <- gsub(pattern = '\\.mouse|\\.human', replacement = '', x = idtype.ref)
  message('reference rownames detected ', toupper(x = species.ref),' with id type ', idtype.ref)
  totalid.ref <- paste0(idtype.ref, '.', species.ref)

  if (totalid == totalid.ref) {
    return(object)
  } else {
    display.names <- setNames(
      c('ENSEMBL gene', 'ENSEMBL transcript', 'gene name', 'transcript name'),
      nm = c('Gene.stable.ID', 'Transcript.stable.ID','Gene.name','Transcript.name')
    )
    if (isRunning()) {
      showNotification(
        paste0(
          "Converted ",species," ",display.names[idtype]," IDs to ",
          species.ref," ",display.names[idtype.ref]," IDs"
        ),
        duration = 3,
        type = 'warning',
        closeButton = TRUE,
        id = 'no-progress-notification'
      )
    }
    # set up table indexed by query ids (totalid)
    linked.unique <- linked[!duplicated(x = linked[[totalid]]), ]
    new.indices <- which(query.names %in% linked.unique[[totalid]])
    message("Found ", length(x = new.indices), " out of ",
            length(x = query.names), " total inputs in conversion table")
    query.names <- query.names[new.indices]
    rownames(x = linked.unique) <- linked.unique[[totalid]]
    linked.unique <- linked.unique[query.names, ]
    # get converted totalid.ref
    new.names <- linked.unique[, totalid.ref]
    # remove duplicates
    notdup <- !duplicated(x = new.names)
    new.indices <- new.indices[notdup]
    new.names <- new.names[notdup]
    # subset/rename object accordingly
    counts <- GetAssayData(object = object[["RNA"]], slot = "counts")[rownames(x = object)[new.indices], ]
    rownames(x = counts) <- new.names
    object <- CreateSeuratObject(
      counts = counts,
      meta.data = object[[]]
    )
    return(object)
  }
}

# Return CSS styling for hover box on interactive plots
#
# @param x X hover position (hover$coords_css$x)
# @param y Y hover position (hover$coords_css$y)
#
# @return Returns a string of CSS to pass to style param
#
HoverBoxStyle <- function(x, y) {
  xpad <- 20 # important to avoid collisions between cursor and hover panel
  ypad <- 20
  paste0(
    "position:absolute; background-color:rgba(245, 245, 245, 0.85); ",
    "left:", (x + xpad), "px; top:",
    (y - ypad), "px;",
    "padding: 5px; margin-bottom: 0px;"
  )
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
#' @importFrom SeuratDisk Connect
#' @importFrom SeuratObject CreateSeuratObject Assays GetAssayData
#' DefaultAssay<-
#' @importFrom Seurat Read10X_h5 as.sparse Assays DietSeurat as.Seurat
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
#' @export
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
      CreateSeuratObject(counts = mat, min.cells = 1, min.features = 1)
    },
    'rds' = {
      object <- readRDS(file = path)
      if (inherits(x = object, what = c('Matrix', 'matrix', 'data.frame'))) {
        object <- CreateSeuratObject(counts = as.sparse(x = object), min.cells = 1, min.features = 1)
      } else if (inherits(x = object, what = 'Seurat')) {
        if (!'RNA' %in% Assays(object = object)) {
          stop("No RNA assay provided", call. = FALSE)
        } else if (Seurat:::IsMatrixEmpty(x = GetAssayData(object = object, slot = 'counts', assay = 'RNA'))) {
          stop("No RNA counts matrix present", call. = FALSE)
        }
        object <- CreateSeuratObject(
          counts = GetAssayData(object = object[["RNA"]], slot = "counts"),
          min.cells = 1,
          min.features = 1,
          meta.data = object[[]]
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
      object <- CreateSeuratObject(
        counts = GetAssayData(object = object[["RNA"]], slot = "counts"),
        min.cells = 1,
        min.features = 1,
        meta.data = object[[]]
      )
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
#' @importFrom SeuratObject AddMetaData CreateSeuratObject
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
    x <- adata[['raw/X']]
  } else if (Exists(name = 'X')) {
    md <- 'var'
    x <- adata[['X']]
  } else {
    stop("Cannot find counts matrix", call. = FALSE)
  }
  # check different possible attributes to try and get matrix shape
  if (isTRUE(x$attr_exists(attr_name = 'h5sparse_shape'))) {
    mtx.shape <- h5attr(x, 'h5sparse_shape')
  } else if (isTRUE(x$attr_exists(attr_name = 'shape'))) {
    mtx.shape <- h5attr(x, 'shape')
  } else {
    warning("Could not determine matrix shape")
  }
  # check different attributes to try and deduce matrix type
  if (isTRUE(x$attr_exists(attr_name = 'h5sparse_format'))) {
    mtx.type <- h5attr(x, 'h5sparse_format')
  } else if (isTRUE(x$attr_exists(attr_name = 'encoding-type'))) {
    mtx.type <- substr(h5attr(x, 'encoding-type'), 0, 3)
  } else {
    mtx.type <- 'csr' # assume matrix is csr
    warning("Could not determine matrix format")
  }
  if (mtx.type != 'csr') {
    p <- as.integer(x[['indptr']][])
    i <- as.integer(x[['indices']][])
    data <- as.double(x[['data']][])
    # csc -> csr
    converted.mtx <- csc_tocsr(
      n_row = as.integer(mtx.shape[1]),
      n_col = as.integer(mtx.shape[2]),
      Ap = p,
      Ai = i,
      Ax = data
    )
    # csr -> dgC
    counts <- new(
      Class = 'dgCMatrix',
      p = converted.mtx$p,
      i = converted.mtx$i,
      x = converted.mtx$x,
      Dim = c(mtx.shape[2], mtx.shape[1])
    )
  } else {
    # x must be a CSR matrix
    counts <- as.matrix(x = x)
  }
  metadata <- LoadMetadata(md = 'obs') # gather additional cell-level metadata
  rownames <- GetRowNames(md = md)
  colnames <- rownames(metadata)
  rownames(x = counts) <- rownames
  colnames(x = counts) <- colnames
  object <- CreateSeuratObject(counts = counts)
  if (ncol(x = metadata)) {
    object <- AddMetaData(object = object, metadata = metadata)
  }
  object <- subset(object, subset = nCount_RNA > 0)
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
#' @return A list with two entries:
#' \describe{
#'  \item{\code{map}}{
#'   The downsampled reference \code{\link[Seurat]{Seurat}}
#'   object (for mapping)
#'  }
#'  \item{\code{plot}}{The reference \code{Seurat} object (for plotting)}
#' }
#'
#' @importFrom SeuratObject Idents<-
#' @importFrom Seurat LoadAnnoyIndex
#' @importFrom httr build_url parse_url status_code GET timeout
#' @importFrom utils download.file
#' @importFrom Matrix sparseMatrix
#'
#' @export
#' @examples
#' \dontrun{
#' # Load from a URL
#' ref <- LoadReference("https://seurat.nygenome.org/references/pbmc")
#' # Load from a directory
#' ref2 <- LoadReference("/var/www/html")
#' }
#'
LoadReference <- function(path, seconds = 10L) {
  ref.names <- list(
    map = 'ref.Rds',
    ann = 'idx.annoy'
  )
  if (substr(x = path, start = nchar(x = path), stop = nchar(x = path)) == '/') {
    path <- substr(x = path, start = 1, stop = nchar(x = path) - 1)
  }
  uri <- httr::build_url(url = httr::parse_url(url = path))
  if (grepl(pattern = '^://', x = uri) | grepl(pattern = '^[a-zA-Z]{1}://', x = uri)) {
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

  # handle new parameters in uwot models beginning in v0.1.13
  if (!"num_precomputed_nns" %in% names(Misc(map[["refUMAP"]])$model)) {
    Misc(map[["refUMAP"]], slot="model")$num_precomputed_nns <- 1
  }

  # Load the annoy index into the Neighbor object in the neighbors slot
  map[["refdr.annoy.neighbors"]] <- LoadAnnoyIndex(
    object = map[["refdr.annoy.neighbors"]],
    file = annref
  )
  # Validate that reference contains required dims
  if (ncol(x = map[["refDR"]]) < getOption(x = "Azimuth.map.ndims")) {
    stop("Provided reference doesn't contain at least ",
         getOption(x = "Azimuth.map.ndims"), " dimensions. Please either
         regenerate reference with requested dimensionality or adjust ",
         "the Azimuth.map.ndims option.")
  }
  # Create plotref
  ad <- Tool(object = map, slot = "AzimuthReference")
  plotref.dr <- GetPlotRef(object = ad)
  cm <- sparseMatrix(
    i = 1, j = 1, x = 0, dims = c(1, nrow(x = plotref.dr)),
    dimnames = list("placeholder", Cells(x = plotref.dr))
  )
  plot <- CreateSeuratObject(
    counts = cm
  )
  plot[["refUMAP"]] <- plotref.dr
  plot <- AddMetaData(object = plot, metadata = Misc(object = plotref.dr, slot = "plot.metadata"))
  gc(verbose = FALSE)
  return(list(
    map = map,
    plot = plot
  ))
}

#' Load the bridge, reference, and reference extension RDS files
#'
#' Read in a reference \code{\link[Seurat]{Seurat}} object, the multiomic bridge, and the precomputed extended reference. This
#' function can read either from URLs or a file path. In order to read properly,
#' there must be the following files:
#' \itemize{
#'  \item \dQuote{ref.Rds} for the downsampled reference \code{Seurat}
#'  object (for mapping)
#'  \item \dQuote{bridge.Rds} for the multiomic bridge \code{Seurat} object 
#'  \item \dQuote{ext.Rds} for the extended reference \code{Seurat} object 
#' }
#'
#' @param path Path or URL to the three RDS files
#' @param seconds Timeout to check for URLs in seconds
#'
#' @return A list with four entries:
#' \describe{
#'  \item{\code{map}}{
#'   The downsampled reference \code{\link[Seurat]{Seurat}}
#'   object (for mapping)
#'  }
#'  \item{\code{plot}}{The reference \code{Seurat} object (for plotting)}
#'  \item{\code{nridhe}}{
#'   The multiomic bridge \code{\link[Seurat]{Seurat}}
#'   object 
#'  }
#'  \item{\code{ext}}{
#'   The extended reference \code{\link[Seurat]{Seurat}}
#'   object
#'  }
#' }
#'
#' @importFrom SeuratObject Idents<-
#' @importFrom Seurat LoadAnnoyIndex
#' @importFrom httr build_url parse_url status_code GET timeout
#' @importFrom utils download.file
#' @importFrom Matrix sparseMatrix
#'
#' @export
#' @examples
#' \dontrun{
#' # Load from a URL
#' ref <- LoadBridgeReference("https://seurat.nygenome.org/references/pbmc")
#' # Load from a directory
#' ref2 <- LoadBridgeReference("/var/www/html")
#' }
#'

LoadBridgeReference<- function(path, seconds = 10L) {
  ref.names <- list(
    bridge = 'bridge.Rds',
    map = 'ref.Rds',
    ext = 'ext.Rds'
  )
  if (substr(x = path, start = nchar(x = path), stop = nchar(x = path)) == '/') {
    path <- substr(x = path, start = 1, stop = nchar(x = path) - 1)
  }
  uri <- httr::build_url(url = httr::parse_url(url = path))
  if (grepl(pattern = '^://', x = uri) | grepl(pattern = '^[a-zA-Z]{1}://', x = uri)) {
    if (!dir.exists(paths = path)) {
      stop("Cannot find directory ", path, call. = FALSE)
    }
    bridgeref <- file.path(path, ref.names$bridge)
    mapref <- file.path(path, ref.names$map)
    extref <- file.path(path, ref.names$ext)
    exists <- file.exists(c(bridgeref, mapref, extref))
    if (!all(exists)) {
      stop(
        "Missing the following files from the directory provided: ",
        Oxford(unlist(x = ref.names)[!exists], join = 'and')
      )
    }
  } else {
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
    #mapref <- url(description = ref.uris[['map']])
    #annref <- tempfile()
    #download.file(url = ref.uris[['ann']], destfile = annref, quiet = TRUE)
    #on.exit(expr = {
    #  close(con = mapref)
    #  unlink(x = annref)
  }
  # Load the map reference
  bridge <- readRDS(file = bridgeref)
  map <- readRDS(file = mapref)
  ext <- readRDS(file = extref)
  
  # handle new parameters in uwot models beginning in v0.1.13
  if (!"num_precomputed_nns" %in% names(Misc(map[["refUMAP"]])$model)) {
    Misc(map[["refUMAP"]], slot="model")$num_precomputed_nns <- 1
  }
  
  # Load the annoy index into the Neighbor object in the neighbors slot
  #map[["refdr.annoy.neighbors"]] <- LoadAnnoyIndex(
  # object = map[["refdr.annoy.neighbors"]],
  #file = annref
  #)
  # Validate that reference contains required dims
  print(ncol(x = map[["refDR"]]))
  print(getOption(x = "Azimuth.map.ndims"))
  if (ncol(x = map[["refDR"]]) < getOption(x = "Azimuth.map.ndims")) {
    stop("Provided reference doesn't contain at least ",
         getOption(x = "Azimuth.map.ndims"), " dimensions. Please either
         regenerate reference with requested dimensionality or adjust ",
         "the Azimuth.map.ndims option.")
  }
  # Create plotref
  ad <- Tool(object = map, slot = "AzimuthReference")
  plotref.dr <- GetPlotRef(object = ad)
  cm <- sparseMatrix(
    i = 1, j = 1, x = 0, dims = c(1, nrow(x = plotref.dr)),
    dimnames = list("placeholder", Cells(x = plotref.dr))
  )
  plot <- CreateSeuratObject(
    counts = cm
  )
  plot[["refUMAP"]] <- plotref.dr
  plot <- AddMetaData(object = plot, metadata = Misc(object = plotref.dr, slot = "plot.metadata"))
  gc(verbose = FALSE)
  return(list(
    map = map,
    plot = plot,
    bridge = bridge,
    ext = ext
  ))
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
#' @importFrom SeuratObject Indices
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

# Check if file is available at given URL
#
# @param url URL to file
# @param strict Ensure http code is 200
#
# @return Boolean indicating file presence
#
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
      httr::timeout(seconds = getOption('timeout')
      ))),
    y = code
  ))
}

# Make An English List
#
# Joins items together to make an English list; uses the Oxford comma for the
# last item in the list.
#
#' @inheritParams base::paste
# @param join either \dQuote{and} or \dQuote{or}
#
# @return A character vector of the values, joined together with commas and
# \code{join}
#
#' @keywords internal
#
# @examples
# \donttest{
# Oxford("red")
# Oxford("red", "blue")
# Oxford("red", "green", "blue")
# }
#
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

# Determine if there are a prohibitive # of annotations for legend
OversizedLegend <- function(annotation.list) {
  return(length(x = unique(x = as.vector(x = annotation.list))) > 50)
}

# Toggle demo button enable/disable
#
# @param action Whether to enable or disable the buttons
# @param demos data.frame containing demo dataset name and file paht
#
# @return No return value
#
ToggleDemos <- function(action = c("enable", "disable"), demos = NULL) {
  if (!is.null(x = demos)) {
    if (action == "enable") {
      for (i in 1:nrow(x = demos)) {
        enable(id = paste0("triggerdemo", i))
      }
    }
   if (action == "disable") {
     for (i in 1:nrow(x = demos)) {
       disable(id = paste0("triggerdemo", i))
     }
   }
  }
}

# Theme for the plot on welcome page
#
WelcomePlot <- function(...) {
  welcomeplot.theme <- theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    legend.position = "none", plot.title = element_blank(),
    # if we want backgroundless... (also replace 'box' with 'fluidRow' in UI)
    panel.background = element_rect(color = '#ecf0f5', fill = '#ecf0f5'),
    plot.background = element_rect(color = '#ecf0f5', fill = '#ecf0f5')
  )
  return(welcomeplot.theme)
}



############### Overlap Functionality ###############
# QC For Overlap between multiomic bridge and atac query 
#
# @param query_assay 
# @param multiome_atac 
#
# @return No return value
#

bridge_qc <- function(query_assay, mutliome_atac){
  o_hits <- GetOverlaps(query_assay, obj.multi[["ATAC"]])
  
  atac_peaks <- OverlapQC(o_hits, query_assay, obj.multi)
  # Calculate the percent of overlap for each overlapping peak 
  # Right now this function is only working on peaks that are definitely overlapping. To generalize it you'd have to add a part that first checks if the ranges are within each other to confirm if theyre overlapping, if not, it should return 0  
  
  #### Density Plot
  d <- density(atac_peaks$perc_overlap) 
  plot(d, xlab='Percentage of Overlap', main = 'Distribution of Overlap Percentages')
  
  #Calculating the percent of query that is covered by the overlapping portions 
  # Sum of all query widths

  o_total <- OverlapTotal(atac_peaks)
  # print this 
}

OverlapDistPlot <- function(query_assay, multiome){
  o_hits <- findOverlaps(query_assay, multiome[["ATAC"]])
  atac_peaks <- OverlapQC(o_hits, query_assay, multiome)
  d <- density(atac_peaks$perc_overlap)
  plot(d, xlab='Percentage of Overlap', main = 'Distribution of Overlap Percentages')
}

# Calculate Overlap percentage for dataframe of overlap info per peak 
#
# @param atac_peaks dataframe with coordinates for each peak and overlap 
#
# @return Percentage of Overlap 
#
PercOverlap <- function(atac_peaks){
  # if no overlap  
  len_overlap <- atac_peaks$o_end - atac_peaks$o_start
  len_q <- atac_peaks$end - atac_peaks$start
  perc <- (len_overlap/len_q)
  return(perc)
}

# Calculate Create dataframe with info for each peaks's overlap in multiome 
#
# @param o_hits Iranges object of overlapping hits 
# @param atac
# @param multi
#
# @return Percentage of Overlap 
#
OverlapQC <- function(o_hits, atac, multi) {
  atac_inds <- queryHits(o_hits)
  multi_inds <- subjectHits(o_hits)
  atac_peaks <- as.data.table(GetAssayData(atac, slot = "ranges")[atac_inds,])
  multi_peaks <- as.data.table(GetAssayData(multi, assay = "ATAC", slot = "ranges")[multi_inds,])
  atac_peaks$o_start <- mapply(max, atac_peaks$start, multi_peaks$start)
  atac_peaks$o_end <- mapply(min, atac_peaks$end, multi_peaks$end)
  atac_peaks$perc_overlap <- PercOverlap(atac_peaks)
  atac_peaks
}

OverlapTotal <- function(atac_peaks){ # from overlap qc 
  q_width <- sum(atac_peaks$width) # but there will be repeats in this 
  o_width <- sum(atac_peaks$o_end - atac_peaks$o_start)
  amount_covered <- (o_width/q_width) * 100
  return(amount_covered)
}

# Requantify atac peaks to either multiomic peaks or to genes 
#
# @param o_hits Iranges object of overlapping hits (Should use same assay as assay for requantification)
# @param ATAC chromatin assay or Seurat Object
# @param subject ATAC assay from Bridge or Transcripts dataframe 
# @param assay assay to use in requantifying peaks to genes (original peaks "peak.orig" or requantified peaks "ATAC")
# @param verbose
#
# @return Percentage of Overlap 
#
RequantifyPeaks <- function(
    #o_hits, 
    atac, 
    subject,
    assay = "peak.orig",
    verbose = TRUE){
  # Query peaks that have overlap w/ multiome peaks
  if (inherits(x = atac, what = "ChromatinAssay")){
    o_hits <- findOverlaps(atac, subject[["ATAC"]])
    atac_inds <- queryHits(o_hits)
    atac_final <- atac[atac_inds, ]
    new_names <- rownames(subject[["ATAC"]][subjectHits(o_hits)]) 
    if (verbose){
      message("Requantifying query peaks to match multiome")
    }
  } else if (inherits(x = atac, what = "Seurat")){ 
    o_hits <- suppressWarnings(findOverlaps(atac[[assay]], transcripts))
    atac_inds <- queryHits(o_hits)
    DefaultAssay(atac) <- assay
    print(atac)
    atac_data <- GetAssayData(atac, assay = assay, slot = "counts")
    atac_final <- atac_data[atac_inds, ]
    new_names <- GRangesToString(subject[subjectHits(o_hits)])
    if (verbose){
      message("Requantifying query peaks to genes")
    }
  } else{
    stop("Incorrect object type ")
  }
  # Reassign query row names
  rownames(atac_final) <- new_names
  # Merge duplicates
  atac_final <- rowsum(atac_final, row.names(atac_final), reorder=FALSE)  
  atac_final <- Matrix::Matrix(atac_final, sparse = TRUE) 
  ##### code from signac 
  if (inherits(x = subject, what = "GRanges")){
    gene.key <- subject$gene_name
    names(x = gene.key) <- GRangesToString(grange = subject)
    rownames(x = atac_final) <- as.vector(x = gene.key[rownames(x = atac_final)])
    atac_final <- atac_final[rownames(x = atac_final) != "", ]
  }
  return(atac_final)
}


############## Gene Activity Scores: FUNCTIONS FROM SIGNAC ############################


CollapseToLongestTranscript <- function(ranges) { ########## from signac 
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

GetTranscripts <- function( ##########From signac
  object,
  assay = NULL,
  features = NULL,
  extend.upstream = 2000,
  extend.downstream = 0,
  biotypes = "protein_coding",
  max.width = 500000,
  process_n = 2000,
  gene.id = FALSE,
  verbose = TRUE
) {
  if (!is.null(x = features)) {
    if (length(x = features) == 0) {
      stop("Empty list of features provided")
    }
  }
  # collapse to longest protein coding transcript
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay.")
  }
  annotation <- Annotation(object = object[[assay]])
  # replace NA names with gene IDD
  annotation$gene_name <- ifelse(
    test = is.na(x = annotation$gene_name) | (annotation$gene_name == ""),
    yes = annotation$gene_id,
    no = annotation$gene_name
  )
  if (length(x = annotation) == 0) {
    stop("No gene annotations present in object")
  }
  if (verbose) {
    message("Extracting gene coordinates")
  }
  transcripts <- CollapseToLongestTranscript(ranges = annotation)
  if (gene.id) {
    transcripts$gene_name <- transcripts$gene_id
  }
  if (!is.null(x = biotypes)) {
    transcripts <- transcripts[transcripts$gene_biotype %in% biotypes]
    if (length(x = transcripts) == 0) {
      stop("No genes remaining after filtering for requested biotypes")
    }
  }
  
  # filter genes if provided
  if (!is.null(x = features)) {
    transcripts <- transcripts[transcripts$gene_name %in% features]
    if (length(x = transcripts) == 0) {
      stop("None of the requested genes were found in the gene annotation")
    }
  }
  if (!is.null(x = max.width)) {
    transcript.keep <- which(x = width(x = transcripts) < max.width)
    transcripts <- transcripts[transcript.keep]
    if (length(x = transcripts) == 0) {
      stop("No genes remaining after filtering for max.width")
    }
  }
  
  # extend to include promoters
  transcripts <- Extend(
    x = transcripts,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  
}


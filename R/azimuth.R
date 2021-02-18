#' @include zzz.R
#' @include helpers.R
#' @include ui.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Launch Command
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Launch the mapping app
#'
#' @param config Path to JSON-formatted configuration file specifying options;
#' for an example config file, see
#' \code{system.file("resources", "config.json", package = "Azimuth")}
#' @param ... Options to set, see \code{?`\link{Azimuth-package}`} for details
#' on \pkg{Azimuth}-provided options
#'
#' @section Specifying options:
#' R options can be provided as named arguments to \code{AzimuthApp} through
#' dots (...), set in a config file, or set globally. Arguments provided to
#' \code{AzimuthApp} through dots take precedence if the same option is provided
#' in a config file. Options provided through dots or a config file take
#' precedence if the same option was set globally.
#'
#' Options in the \code{\link[Azimuth:Azimuth-package]{Azimuth.app}} namespace
#' can be specified using a shorthand notation in both the config file and as
#' arguments to \code{AzimuthApp}. For example, the option
#' \code{Azimuth.app.reference} can be shortened to \code{reference} in the
#' config file or as an argument to \code{AzimuthApp}
#'
#' @return None, launches the mapping Shiny app
#'
#' @importFrom shiny runApp shinyApp
#' @importFrom withr with_options
#' @importFrom jsonlite read_json
#'
#' @export
#'
#' @seealso \code{\link{Azimuth-package}}
#'
#' @examples
#' if (interactive()) {
#'   AzimuthApp(system.file("resources", "config.json", package = "Azimuth"))
#' }
#'
AzimuthApp <- function(config = NULL, ...) {
  useShinyjs()
  # If multiple items have the same name in the named list, with_options sets
  # the option to the last entry with that name in the list. Therefore, putting
  # the config file options first, followed by options set in dots, followed by
  # hardcoded options, achieves the desired precedence.
  opts <- list()
  # Add options set through config file
  if (!is.null(x = config)) {
    opts <- c(opts, read_json(path = config, simplifyVector = TRUE))
  }
  # Add options set through named arguments in dots
  args <- list(...)
  if (length(x = args) && !is.null(x = names(x = args))) {
    # only add named elements
    opts <- c(opts, args[names(x = args) != ""])
  }
  # if any arguments from dots or config file have no "." character,
  # prepend the "Azimuth.app" namespace
  for (i in seq_along(along.with = opts)) {
    if (!grepl(pattern = '\\.', x = names(x = opts)[i])) {
      names(x = opts)[i] <- paste0('Azimuth.app.', names(x = opts)[i])
    }
  }
  # Add sensible defaults
  # Shiny doesn't set shiny.maxRequestSize on load
  if (!'shiny.maxRequestSize' %in% names(x = opts) && is.null(x = getOption(x = 'shiny.maxRequestSize'))) {
    opts$shiny.maxRequestSize <- 500 * (1024 ^ 2)
  }
  # Add pageLength to jQuery DataTables options
  opts$DT.options <- as.list(x = c(
    opts$DT.options,
    getOption(x = 'DT.options')
  ))
  if (!'pageLength' %in% names(x = opts$DT.options)) {
    opts$DT.options$pageLength <- 10L
  }
  # Set future.globals.maxSize; this is not user-configurable
  maxcells <- with_options(
    new = opts,
    code = getOption(x = 'Azimuth.app.max_cells')
  )
  opts$future.globals.maxSize <- maxcells * 320000
  # Launch the app
  with_options(
    new = opts,
    code = runApp(appDir = shinyApp(ui = AzimuthUI, server = AzimuthServer))
  )
  return(invisible(x = NULL))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' AzimuthData
#'
#' The AzimuthData class is used to store reference info needed for Azimuth
#'
#' @slot plotref DimReduc object containing UMAP for plotting and projection.
#' This should also contain the cell IDs in the misc slot
#' @slot colormap Vector of id-color mapping for specifying the plots.
#' @slot seurat.version Version of Seurat used in reference construction
#' @slot azimuth.version Version of Azimuth used in reference construction
#' @slot reference.version Version of the Azimuth reference
#'
#' @name AzimuthData-class
#' @rdname AzimuthData-class
#' @exportClass AzimuthData
#'
AzimuthData <- setClass(
  Class = 'AzimuthData',
  slots = c(
    plotref = 'DimReduc',
    colormap = 'list',
    seurat.version = 'package_version',
    azimuth.version = 'package_version',
    reference.version = 'character'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Get Azimuth color mapping
#'
#' Pull ID-color mapping for Azimuth plotting
#'
#' @param object An object
#' @param slot Name of tool
#' @param ... Arguments passed to other methods
#'
#' @return A named vector specifying the colors for all reference IDs
#'
#' @rdname GetColorMap
#' @export GetColorMap
#'
GetColorMap <- function(object, ...) {
  UseMethod(generic = 'GetColorMap', object = object)
}

#' Get Azimuth plotref
#'
#' Pull DimReduc used in Azimuth plotting/projection
#'
#' @inheritParams GetColorMap
#'
#' @return A DimReduc object
#'
#' @rdname GetPlotRef
#' @export GetPlotRef
#'
GetPlotRef <- function(object, ...) {
  UseMethod(generic = 'GetPlotRef', object = object)
}

#' Get Azimuth reference version number
#'
#' Pull the reference version information
#'
#' @return A character string specifying the reference version
#'
#' @rdname ReferenceVersion
#' @export ReferenceVersion
#'
ReferenceVersion <- function(object, ...) {
  UseMethod(generic = 'ReferenceVersion', object = object)
}

#' Set Azimuth color mapping
#'
#' Set ID-color mapping for Azimuth plotting
#'
#' @inheritParams GetColorMap
#'
#' @return An object with the colormap slot set
#'
#' @rdname SetColorMap
#' @export SetColorMap
#'
SetColorMap <- function(object, ...) {
  UseMethod(generic = 'SetColorMap', object = object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname GetColorMap
#' @export
#' @method GetColorMap AzimuthData
#'
GetColorMap.AzimuthData <- function(object, ...) {
  return(slot(object = object, name = "colormap"))
}

#' @rdname GetColorMap
#' @export
#' @method GetColorMap Seurat
#'
GetColorMap.Seurat <- function(object, slot = "AzimuthReference", ...) {
  return(GetColorMap(object = Tool(object = object, slot = slot)))
}

#' @rdname GetPlotRef
#' @export
#' @method GetPlotRef AzimuthData
#'
GetPlotRef.AzimuthData <- function(object, ...) {
  return(slot(object = object, name = "plotref"))
}

#' @rdname GetPlotRef
#' @export
#' @method GetPlotRef Seurat
#'
GetPlotRef.Seurat <- function(object, slot = "AzimuthReference", ...) {
  return(GetPlotRef(object = Tool(object = object, slot = slot)))
}

#' @param object Seurat or AzimuthData object
#' @param slot Name of the version to pull. Can be "seurat.version",
#' "azimuth.version", or "reference.version".
#' @param ... Not used
#' @rdname ReferenceVersion
#' @export
#' @method ReferenceVersion AzimuthData
#'
ReferenceVersion.AzimuthData <- function(object, ...) {
  return(slot(object = object, name = "reference.version"))
}

#' @rdname ReferenceVersion
#' @export
#' @method ReferenceVersion Seurat
ReferenceVersion.Seurat <- function(object, slot = "AzimuthReference", ...) {
  return(ReferenceVersion(object = Tool(object = object, slot = slot)))
}

#' @rdname SetColorMap
#' @param value New colormap to assign
#' @export
#' @method SetColorMap AzimuthData
#'
SetColorMap.AzimuthData <- function(object, value, ...) {
  return(slot(object = object, name = "colormap") <- value)
}

#' @rdname SetColorMap
#' @export
#' @method SetColorMap Seurat
#'
SetColorMap.Seurat <- function(object, slot = "AzimuthReference", value, ...) {
  return(SetColorMap(object = Tool(object = object, slot = slot), value = value, ...))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# AzimuthData Helper Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a Seurat object compatible with Azimuth.
#'
#' @inheritParams CreateAzimuthData
#' @param refUMAP Name of UMAP in reference to use for mapping
#' @param refDR Name of DimReduc in reference to use for mapping
#' @param refAssay Name of SCTAssay to use in reference
#' @param dims Dimensions to use in reference neighbor finding
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param ori.index Index of the cells used in mapping in the original object on
#' which UMAP was run. Only need to provide if UMAP was run on different set of
#' cells.
#' @param assays Assays to retain for transfer
#' @param metadata Metadata to retain for transfer
#' @param verbose Display progress/messages
#'
#' @return Returns a Seurat object with AzimuthData stored in the tools slot for
#' use with Azimuth.
#'
#' @importFrom SeuratObject Reductions Misc Misc<- Assays Cells Loadings Idents
#' DefaultAssay Tool<-
#' @importFrom Seurat FindNeighbors NormalizeData AverageExpression DietSeurat
#' @importFrom methods as
#'
#' @export
#'
AzimuthReference <- function(
  object,
  refUMAP = "umap",
  refDR = "spca",
  refAssay = "SCT",
  dims = 1:50,
  k.param = 31,
  plotref = "umap",
  plot.metadata = NULL,
  ori.index = NULL,
  colormap = NULL,
  assays = NULL,
  metadata = NULL,
  reference.version = "0.0.0",
  verbose = FALSE
) {
  # Parameter validation
  if (!refUMAP %in% Reductions(object = object)) {
    stop("refUMAP (", refUMAP, ") not found in Seurat object provided")
  }
  if (is.null(x = Misc(object = object[[refUMAP]], slot = "model"))) {
    stop("refUMAP (", refUMAP, ") does not have the umap model info stored. ",
         "Please rerun RunUMAP with return.model = TRUE.")
  }
  if (!refDR %in% Reductions(object = object)) {
    stop("refDR (", refDR, ") not found in Seurat object provided")
  }
  if (is.null(x = metadata)) {
    stop("Please specify at least one metadata field (for transfer and plotting).")
  }
  for(i in metadata) {
    if (! i %in% colnames(x = object[[]])) {
      warning(i, " not found in Seurat object metadata")
      next
    }
    if (! is.factor(x = object[[i, drop = TRUE]])) {
      warning(i, " is not a factor. Converting to factor with alphabetical ",
              "levels.", call. = FALSE)
      object[[i, drop = TRUE]] <- factor(x = object[[i, drop = TRUE]], levels = sort(x = unique(object[[i, drop = TRUE]])))
    }
  }
  if (!refAssay %in% Assays(object = object)) {
    stop("Seurat object provided must have the SCT Assay stored.")
  }
  if (!inherits(x = object[[refAssay]], what = "SCTAssay")) {
    stop("refAssay (", refAssay, ") is not an SCTAssay.")
  }
  if (length(x = levels(x = object[[refAssay]])) != 1) {
    stop("refAssay (", refAssay, ") should contain a single SCT model.")
  }

  suppressWarnings(expr = object[["refUMAP"]] <- object[[refUMAP]])
  suppressWarnings(expr = object[["refDR"]] <- object[[refDR]])

  # Calculate the Neighbors
  object <- FindNeighbors(
    object = object,
    reduction = "refDR",
    dims = dims,
    graph.name = "refdr.annoy.neighbors",
    k.param = k.param,
    nn.method = "annoy",
    annoy.metric = "cosine",
    cache.index = TRUE,
    return.neighbor = TRUE,
    l2.norm = FALSE,
    verbose = verbose
  )
  if (verbose) {
    message("Computing pseudobulk averages")
  }
  features <- rownames(x = Loadings(object = object[['refDR']]))
  plot.metadata <- plot.metadata %||% object[[metadata]]
  if (inherits(x = plotref, what = "DimReduc")) {
    plot.metadata <- plot.metadata[Cells(x = plotref), ]
  }
  ad <- CreateAzimuthData(
    object = object,
    plotref = plotref,
    plot.metadata  = plot.metadata,
    colormap = colormap,
    reference.version = reference.version
  )
  # Add the "ori.index" column.
  ori.index <- ori.index %||% match(Cells(x = object), Cells(x = object[["refUMAP"]]))
  object$ori.index <- ori.index

  # Subset the features of the RNA assay
  DefaultAssay(object = object) <- refAssay
  object[[refAssay]] <- subset(x = object[[refAssay]], features = features)
  # Preserves DR after DietSeurat
  DefaultAssay(object = object[["refDR"]]) <- refAssay
  object <- DietSeurat(
    object = object,
    counts = FALSE,
    assays = c(refAssay, assays),
    dimreducs = c("refDR","refUMAP")
  )
  metadata <- c(metadata, "ori.index")
  for (i in colnames(x = object[[]])) {
    if (!i %in% metadata){
      object[[i]] <- NULL
    }
  }
  sct.model <- slot(object = object[[refAssay]], name = "SCTModel.list")[[1]]
  object[["refAssay"]] <- as(object = suppressWarnings(Seurat:::CreateDummyAssay(assay = object[[refAssay]])), Class = "SCTAssay")
  slot(object = object[["refAssay"]], name = "SCTModel.list") <- list(refmodel = sct.model)
  DefaultAssay(object = object) <- "refAssay"
  DefaultAssay(object = object[["refDR"]]) <- "refAssay"
  Tool(object = object) <- ad
  object <- DietSeurat(
    object = object,
    counts = FALSE,
    assays = c("refAssay", assays),
    dimreducs = c("refDR","refUMAP")
  )
  ValidateAzimuthReference(object = object)
  return(object)
}

#' Create an \code{\link{AzimuthData}} object
#'
#' Create an auxiliary \code{\link{AzimuthData}} object for storing necessary
#' info when generating an Azimuth reference.
#'
#' @param object Seurat object
#' @param plotref Either the name of the DimReduc in the provided Seurat object
#' to use for the plotting reference or the DimReduc object itself.
#' @param plot.metadata A data.frame of discrete metadata fields for the cells
#' in the plotref.
#' @param colormap A list of named and ordered vectors specifying the colors and levels
#' for the metadata. See \code{\link{CreateColorMap}} for help
#' generating your own.
#' @param reference.version Version of the Azimuth reference
#'
#' @return Returns an \code{\link{AzimuthData}} object
#'
#' @importFrom SeuratObject Reductions Misc<-
#'
#' @export
#'
CreateAzimuthData <- function(
  object,
  plotref = "umap",
  plot.metadata = NULL,
  colormap = NULL,
  reference.version = '0.0.0'
) {
  if (inherits(x = plotref, what = "character")) {
    if (plotref %in% Reductions(object = object)) {
      plotref <- object[[plotref]]
    } else {
      stop("The DimReduc ", plotref, " was not found in the provided object.")
    }
  }
  plot.metadata <- plot.metadata %||% data.frame(id = Idents(object = object))
  if (is.null(x = colormap)) {
    colormap <- lapply(X = colnames(x = plot.metadata), FUN = function(x) {
      if (is.factor(x = plot.metadata[, x])) {
        return(CreateColorMap(ids = levels(x = plot.metadata[, x])))
      } else {
        CreateColorMap(ids = sort(x = unique(x = plot.metadata[, x])))
      }
    })
    names(x = colormap) <- colnames(x = plot.metadata)
  }
  for (i in colnames(x = plot.metadata)) {
    if (! i %in% names(x = colormap)) {
      colormap[i] <- list(CreateColorMap(ids = unique(x = plot.metadata[[i]])))
    }
    plot.metadata[[i]] <- factor(x = plot.metadata[[i]], levels = names(x = colormap[[i]]))
  }
  slot(object = plotref, name = 'misc')[["plot.metadata"]] <- plot.metadata
  colormap <- colormap[colnames(x = plot.metadata)]
  ad <- new(
    Class = "AzimuthData",
    plotref = plotref,
    colormap = colormap,
    seurat.version = packageVersion("Seurat"),
    azimuth.version = packageVersion("Azimuth"),
    reference.version = reference.version
  )
  return(ad)
}

#' Create A Color Map
#'
#' Create mapping between IDs and colors to use with reference plotting in
#' Azimuth
#'
#' @param object Seurat object
#' @param ids Vector of IDs to link to colors
#' @param colors Vector of colors to use
#' @param seed Set to randomly shuffle color assignments
#'
#' @return A named vector of colors
#'
#' @importFrom scales hue_pal
#'
#' @export
#'
CreateColorMap <- function(object, ids = NULL, colors = NULL, seed = NULL) {
  ids <- ids %||% levels(x = object)
  colors <- colors %||% hue_pal()(n = length(x = ids))
  if (length(x = ids) != length(x = colors)) {
    stop("Please provide equal length vectors for ids and colors.")
  }
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
    colors <- sample(x = colors)
  }
  names(x = colors) <- ids
  return(colors)
}

# Cluster preservation score
#
# @param query Query object
# @param ds.amount Amount to downsample query
# @return Returns
#
#' @importFrom SeuratObject Cells Idents Indices as.Neighbor
#' @importFrom Seurat RunPCA FindNeighbors FindClusters MinMax
#
#' @keywords internal
#
#
ClusterPreservationScore <- function(query, ds.amount) {
  query <- DietSeurat(object = query, assays = "refAssay", scale.data = TRUE, counts = FALSE, dimreducs = "integrated_dr")
  if (ncol(x = query) > ds.amount) {
    query <- subset(x = query, cells = sample(x = Cells(x = query), size = ds.amount))
  }
  query <- RunPCA(object = query, verbose = FALSE)
  query <- FindNeighbors(object = query, reduction = 'pca', dims = 1:50, graph.name = 'pca_snn')
  query[["orig_neighbors"]] <- as.Neighbor(x = query[["pca_snn"]])
  query <- FindClusters(object = query, resolution = 0.6, graph.name = 'pca_snn')
  query <- FindNeighbors(object = query, reduction = 'integrated_dr', dims = 1:50, return.neighbor = TRUE, graph.name = "integrated_neighbors")
  ids <- Idents(object = query)
  integrated.neighbor.indices <- Indices(object = query[["integrated_neighbors"]])
  proj_ent <- unlist(x = lapply(X = 1:length(x = Cells(x = query)), function(x) {
    neighbors <- integrated.neighbor.indices[x, ]
    nn_ids <- ids[neighbors]
    p_x <- prop.table(x = table(nn_ids))
    nn_entropy <- sum(p_x * log(x = p_x), na.rm = TRUE)
    return(nn_entropy)
  }))
  names(x = proj_ent) <- Cells(x = query)
  orig.neighbor.indices <- Indices(object = query[["orig_neighbors"]])
  orig_ent <- unlist(x = lapply(X = 1:length(x = Cells(x = query)), function(x) {
    neighbors <- orig.neighbor.indices[x, ]
    nn_ids <- ids[neighbors]
    p_x <- prop.table(x = table(nn_ids))
    nn_entropy <- sum(p_x * log(x = p_x), na.rm = TRUE)
    return(nn_entropy)
  }))
  names(x = orig_ent) <- Cells(x = query)
  stat <- median(
    x = tapply(X = orig_ent, INDEX = ids, FUN = mean) -
      tapply(X = proj_ent, INDEX = ids, FUN = mean)
  )
  if (stat <= 0) {
    stat <- 5.00
  } else {
    stat <- -1 * log2(x = stat)
    stat <- MinMax(data = stat, min = 0.00, max = 5.00)
  }
  return(stat)
}

#' Validate References for Azimuth
#'
#' Validate aspects of a Seurat object to be used as an Azimuth reference
#'
#' @param object Seurat object
#' @param ad.name Name in the tools slot containing the AzimuthData object.
#'
#' @return No return value
#'
#' @importFrom SeuratObject Tool Misc Reductions
#'
#' @export
#'
ValidateAzimuthReference <- function(object, ad.name = "AzimuthReference") {
  if (!inherits(x = Tool(object = object, slot = ad.name), what = "AzimuthData")) {
    stop ("Reference must contain an AzimuthData object in the tools slot.")
  }
  plotref <- GetPlotRef(object = object, slot = ad.name)
  colormap <- GetColorMap(object = object, slot = ad.name)
  # plotref needs to have IDs in misc
  plotids <- Misc(object = plotref, slot = "plot.metadata")
  if (is.null(x = plotids)) {
    stop("plotref in AzimuthData object must contain plot.metadata in the misc slot.")
  } else {
    for (id in colnames(x = plotids)) {
      if (length(x = plotids[[id]]) != nrow(x = plotref)) {
        stop(
          "Length of ", id, " in plotref in the AzimuthData object is not equal to ",
          "the number of cells in plotref."
        )
      }
      if (!all(sort(x = as.character(unique(x = plotids[[id]]))) == sort(x = names(x = colormap[[id]])))) {
        stop(
          "The colormap stored in the AzimuthData object must contain a ",
          "color-id mapping for every unique id present in the plotting ids."
        )
      }
    }
  }
  if (is.null(x = Misc(object = plotref, slot = "model"))) {
    stop("plotref must contain the umap model.")
  }
  if (!"refDR" %in% Reductions(object = object)) {
    stop("Object must contain a DimReduc called refDR to use in transfer/projection.")
  }
  if (!"ori.index" %in% colnames(x = object[[]])){
    stop(
      "Seurat object metadata must contain 'ori.index' field, storing the ",
      "mapping between the index of the cells in the object UMAP was run on ",
      "and the cell indices in the object here."
    )
  }
  if (!"refAssay" %in% Assays(object = object)) {
    stop("Must contain assay called 'refAssay'.")
  }
  if (!inherits(x = object[["refAssay"]], what = "SCTAssay")) {
    stop("refAssay must be an SCTAssay object.")
  }
  if (!"refmodel" %in% levels(x = object[["refAssay"]])) {
    stop("refAssay must contain the SCTModel called refmodel.")
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMethod(
  f = 'show',
  signature = 'AzimuthData',
  definition = function(object) {
    cat('An AzimuthData object - reference version:', slot(object = object, name = "reference.version"),
        '\nContains', length(x = GetColorMap(object = object)), 'meta.data field(s) to transfer.')
  }
)

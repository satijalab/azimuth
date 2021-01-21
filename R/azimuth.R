#' @include zzz.R
#' @include seurat.R
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
#' @slot avgref Average RNA expression for pseudobulk correlation tests
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
    avgref = 'matrix',
    colormap = 'list',
    seurat.version = 'package_version',
    azimuth.version = 'package_version',
    reference.version = 'character'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get Azimuth average expression
#'
#' Pull reference average RNA expression matrix
#'
#' @param object An object
#' @param slot Name of tool
#' @param ... Arguments passed to other methods
#'
#' @return A feature x 1 matrix of average RNA expression values
#'
#' @rdname GetAvgRef
#' @export GetAvgRef
#'
GetAvgRef <- function(object, ...) {
  UseMethod(generic = 'GetAvgRef', object = object)
}

#' Get Azimuth color mapping
#'
#' Pull ID-color mapping for Azimuth plotting
#'
#' @inheritParams GetAvgRef
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
#' @inheritParams GetAvgRef
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
#' @inheritParams GetAvgRef
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

#' @rdname GetAvgRef
#' @export
#' @method GetAvgRef AzimuthData
#'
GetAvgRef.AzimuthData <- function(object, ...) {
  return(slot(object = object, name = "avgref"))
}

#' @rdname GetAvgRef
#' @export
#' @method GetAvgRef Seurat
#'
GetAvgRef.Seurat <- function(object, slot = "AzimuthReference", ...) {
  return(GetAvgRef(object = Tool(object = object, slot = slot)))
}

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
#'
#' @export
#'
AzimuthReference <- function(
  object,
  refUMAP = "umap",
  refDR = "spca",
  dims = 1:50,
  k.param = 31,
  plotref = "umap",
  plot.metadata = NULL,
  ori.index = NULL,
  avgref = NULL,
  colormap = NULL,
  assays = NULL,
  metadata = NULL,
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
  if (!"SCT" %in% Assays(object = object)) {
    stop("Seurat object provided must have the SCT Assay stored.")
  }
  if (is.null(x = avgref) && !"RNA" %in% Assays(object = object)) {
    stop("Please provide either the RNA assay or avgref.")
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
  if (is.null(x = avgref)) {
    random.name <- "allcells"
    while (random.name %in% colnames(x = object[[]])) {
      random.name <- paste0(sample(letters, size = 10), collapse = "")
    }
    Idents(object = object) <- random.name
    object <- NormalizeData(object = object, assay = "RNA", verbose = verbose)
    avgref <- suppressMessages(expr = AverageExpression(
      object = object,
      assays = "RNA",
      verbose = verbose
    ))[[1]][features, , drop = FALSE]
    Idents(object = object) <- "id"
  }
  plot.metadata <- plot.metadata %||% object[[metadata]]
  if (inherits(x = plotref, what = "DimReduc")) {
    plot.metadata <- plot.metadata[Cells(x = plotref), ]
  }
  ad <- CreateAzimuthData(
    object = object,
    plotref = plotref,
    plot.metadata  = plot.metadata,
    avgref = avgref,
    colormap = colormap
  )
  # Add the "ori.index" column.
  ori.index <- ori.index %||% match(Cells(x = object), Cells(x = object[["refUMAP"]]))
  object$ori.index <- ori.index

  # Subset the features of the RNA assay
  DefaultAssay(object = object) <- "SCT"
  object[["SCT"]] <- subset(x = object[["SCT"]], features = features)
  # Preserves DR after DietSeurat
  DefaultAssay(object = object[["refDR"]]) <- "SCT"
  object <- DietSeurat(
    object = object,
    counts = FALSE,
    assays = c("SCT", assays),
    dimreducs = c("refDR","refUMAP")
  )
  metadata <- c(metadata, "ori.index")
  for (i in colnames(x = object[[]])) {
    if (!i %in% metadata){
      object[[i]] <- NULL
    }
  }
  object[["SCT"]] <- as(object = suppressWarnings(Seurat:::CreateDummyAssay(assay = object[["SCT"]])), Class = "SCTAssay")
  Misc(object = object[["SCT"]], slot = "vst.set") <- list()
  Tool(object = object) <- ad
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
#' @param avgref Matrix containing the average RNA expression for the reference,
#' used in the pseudobulk correlation test.
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
  avgref = NULL,
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
  avgref <- as.matrix(x = avgref)
  ad <- new(
    Class = "AzimuthData",
    plotref = plotref,
    avgref = avgref,
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
  avgref <- GetAvgRef(object = object, slot = ad.name)
  if (!all(rownames(x = Loadings(object = object[['refDR']])) %in% rownames(x = avgref))){
    stop(
      "avgref must contain average expression values for all features used ",
      "when computing refDR."
    )
  }
  if (!"ori.index" %in% colnames(x = object[[]])){
    stop(
      "Seurat object metadata must contain 'ori.index' field, storing the ",
      "mapping between the index of the cells in the object UMAP was run on ",
      "and the cell indices in the object here."
    )
  }
}

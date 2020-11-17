#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The AzimuthData class is used to store reference info needed for Azimuth
#'
#' @slot plotref DimReduc object containing UMAP for plotting and projection.
#' This should also contain the cell IDs in the misc slot
#' @slot avgref Average RNA expression for pseudobulk correlation tests
#' @slot colormap Vector of id-color mapping for specifying the plots.
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
    colormap = 'vector'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get Azimuth plotref
#'
#' Pull DimReduc used in Azimuth plotting/projection
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return A DimReduc object
#'
#' @rdname GetPlotRef
#' @export GetPlotRef
#'
GetPlotRef <- function(object, ...) {
  UseMethod(generic = 'GetPlotRef', object = object)
}

#' @rdname GetPlotRef
#' @export
#' @method GetPlotRef AzimuthData
#'
GetPlotRef.AzimuthData <- function(object) {
  return(slot(object = object, name = "plotref"))
}

#' @param slot Name of tool
#'
#' @rdname GetPlotRef
#' @export
#' @method GetPlotRef Seurat
#'
GetPlotRef.Seurat <- function(object, slot = "Azimuth") {
  return(GetPlotRef(object = Tool(object = object, slot = slot)))
}

#' Get Azimuth average expression
#'
#' Pull reference average RNA expression matrix
#'
#' @param object An object
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

#' @rdname GetAvgRef
#' @export
#' @method GetAvgRef AzimuthData
#'
GetAvgRef.AzimuthData <- function(object) {
  return(slot(object = object, name = "avgref"))
}

#' @param slot Name of tool
#'
#' @rdname GetAvgRef
#' @export
#' @method GetAvgRef Seurat
#'
GetAvgRef.Seurat <- function(object, slot = "Azimuth") {
  return(GetAvgRef(object = Tool(object = object, slot = slot)))
}

#' Get Azimuth color mapping
#'
#' Pull ID-color mapping for Azimuth plotting
#'
#' @param object An object
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

#' @rdname GetColorMap
#' @export
#' @method GetColorMap AzimuthData
#'
GetColorMap.AzimuthData <- function(object) {
  return(slot(object = object, name = "colormap"))
}

#' @param slot Name of tool
#'
#' @rdname GetColorMap
#' @export
#' @method GetColorMap Seurat
#'
GetColorMap.Seurat <- function(object, slot = "Azimuth") {
  return(GetColorMap(object = Tool(object = object, slot = slot)))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# AzimuthData helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create auxiliary AzimuthData object for storing necessary info when generating
#' an Azimuth reference.
#'
#' @param object Seurat object
#' @param plotref Either the name of the DimReduc in the provided Seurat object
#' to use for the plotting reference or the DimReduc object itself.
#' @param plotids Vector of IDs specifying the identities of the cells in the
#' plotref.
#' @param avgref Matrix containing the average RNA expression for the reference,
#' used in the pseudobulk correlation test.
#' @param colormap A named and ordered vector specifying the colors and levels
#' for the IDs being plotted. See \link{\code{CreateColorMap}} for help
#' generating your own.
#' @param ad.name Name in the tools slot to store the AzimuthData object.
#' Azimuth is expecting the name "Azimuth", so generally don't modify this.
#'
#' @return Returns a Seurat object with the AzimuthData object stored in the
#' tools slot
#'
#' @importFrom Seurat Reductions Misc<-
#' @export
#'
CreateAzimuthData <- function(
  object,
  plotref = "umap",
  plotids = NULL,
  avgref = NULL,
  colormap = NULL,
  ad.name = "Azimuth"
) {
  if (inherits(x = plotref, what = "character")) {
    if (plotref %in% Reductions(object = object)) {
      plotref <- object[[plotref]]
    } else {
      stop("The DimReduc ", plotref, " was not found in the provided object.")
    }
  }
  plotids <- plotids %||% Idents(object = object)
  Misc(object = plotref, slot = "ids") <- plotids
  colormap <- colormap %||% CreateColorMap(object = object)
  avgref <- as.matrix(x = avgref)
  slot(object = object, name = "tools")[[ad.name]] <- AzimuthData(
    plotref = plotref,
    avgref = avgref,
    colormap = colormap
  )
  ValidateAzimuthReference(object = object, ad.name = ad.name)
  return(object)
}

#' Create mapping between IDs and colors to use with reference plotting in
#' Azimuth
#'
#' @param object Seurat object
#' @param ids Vector of IDs to link to colors
#' @param colors Vector of colors to use
#' @param seed Set to randomly shuffle color assignments
#'
#' @export
#'
CreateColorMap <- function(object, ids = NULL, colors = NULL, seed = NULL) {
  ids <- ids %||% levels(x = object)
  colors <- colors %||% gg_color_hue(n = length(x = ids))
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

# From: https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


ValidateAzimuthReference <- function(object, ad.name = "Azimuth") {
  if (!inherits(x = Tool(object = object, slot = ad.name), what = "AzimuthData")) {
    stop ("Reference must contain an AzimuthData object in the tools slot.")
  }
  plotref <- GetPlotRef(object = object, slot = ad.name)
  colormap <- GetColorMap(object = object, slot = ad.name)
  # plotref needs to have IDs in misc
  plotids <- Misc(object = plotref, slot = "ids")
  if (is.null(x = plotids)) {
    stop("plotref in AzimuthData object must contain ids in the misc slot.")
  } else {
    if (length(x = plotids) != nrow(x = plotref)) {
      stop("Length of ids in plotref in the AzimuthData object is not equal to ",
           "the number of cells in plotref.")
    }
    if (!all(sort(x = as.character(unique(x = plotids))) == sort(x = names(x = colormap)))) {
      stop("The colormap stored in the AzimuthData object must contain a ",
           "color-id mapping for every unique id present in the plotting ids.")
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
    stop("avgref must contain average expression values for all features used ",
         "when computing refDR.")
  }
}

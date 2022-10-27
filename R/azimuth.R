#' @include zzz.R
#' @include helpers.R
#' @include ui.R
#'
NULL

#' @inheritParams RunAzimuth
#' @param reference Name of reference to map to or a path to a directory containing ref.Rds and idx.annoy
#' @param annotation.levels list of annotation levels to map. If not specified, all will be mapped.
#' @param umap.name name of umap reduction in the returned object
#' @param do.adt transfer ADT assay
#' @param assay query assay name
#'
#' @return Seurat object with reference reductions and annotations
#'
#' @importFrom SeuratData InstallData InstalledData LoadData AvailableData
#'
#' @export
#' @method RunAzimuth Seurat
#' @rdname RunAzimuth
#'
RunAzimuth.Seurat <- function(
  query,
  reference,
  annotation.levels = NULL,
  umap.name = "ref.umap",
  do.adt = FALSE,
  verbose = TRUE,
  assay = "RNA",
  k.weight = 50,
  n.trees = 20,
  mapping.score.k = 100
) {
  if (dir.exists(reference)) {
    reference <- LoadReference(reference)$map
  } else {
    reference <- tolower(reference)
    if (reference %in% InstalledData()$Dataset) {
      # only get the `map` object since no plotting is performed
      reference <- LoadData(reference, type = "azimuth")$map
    } else if (reference %in% AvailableData()$Dataset) {
      InstallData(reference)
      # only get the `map` object since no plotting is performed
      reference <- LoadData(reference, type = "azimuth")$map
    } else {
      possible.references <- AvailableData()$Dataset[grepl("*ref", AvailableData()$Dataset)]
      print("Choose one of:")
      print(possible.references)
      stop(paste("Could not find a reference for", reference))
    }
    # handle expected new parameters in uwot models beginning in v0.1.13
    if (!"num_precomputed_nns" %in% names(Misc(reference[["refUMAP"]])$model)) {
      Misc(reference[["refUMAP"]], slot="model")$num_precomputed_nns <- 1
    }
  }
  dims <- as.double(length(slot(reference, "reductions")$refDR))
  if (isTRUE(do.adt) && !("ADT" %in% Assays(reference))) {
    warning("Cannot impute an ADT assay because the reference does not have antibody data")
    do.adt = FALSE
  }
  reference.version <- ReferenceVersion(reference)
  azimuth.version <- as.character(packageVersion(pkg = "Azimuth"))
  seurat.version <- as.character(packageVersion(pkg = "Seurat"))
  meta.data <- names(slot(reference, "meta.data"))

  # is annotation levels are not specify, gather all levels of annotation
  if (is.null(annotation.levels)) {
    annotation.levels <- names(slot(object = reference, name = "meta.data"))
    annotation.levels <- annotation.levels[!grepl(pattern = "^nCount", x = annotation.levels)]
    annotation.levels <- annotation.levels[!grepl(pattern = "^nFeature", x = annotation.levels)]
    annotation.levels <- annotation.levels[!grepl(pattern = "^ori", x = annotation.levels)]
  }

  # Change the file path based on where the query file is located on your system.
  query <- ConvertGeneNames(
    object = query,
    reference.names = rownames(x = reference),
    homolog.table = 'https://seurat.nygenome.org/azimuth/references/homologs.rds'
  )

  # Calculate nCount_RNA and nFeature_RNA if the query does not
  # contain them already
  if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
      calcn <- as.data.frame(x = Seurat:::CalcN(object = query))
      colnames(x = calcn) <- paste(
        colnames(x = calcn),
        "RNA",
        sep = '_'
      )
      query <- AddMetaData(
        object = query,
        metadata = calcn
      )
      rm(calcn)
  }

  # Calculate percent mitochondrial genes if the query contains genes
  # matching the regular expression "^MT-"
  if (any(grepl(pattern = '^MT-', x = rownames(x = query)))) {
    query <- PercentageFeatureSet(
      object = query,
      pattern = '^MT-',
      col.name = 'percent.mt',
      assay = assay
    )
  }

  # Preprocess with SCTransform
  query <- SCTransform(
    object = query,
    assay = assay,
    new.assay.name = "refAssay",
    residual.features = rownames(x = reference),
    reference.SCT.model = reference[["refAssay"]]@SCTModel.list$refmodel,
    method = 'glmGamPoi',
    ncells = 2000,
    n_genes = 2000,
    do.correct.umi = FALSE,
    do.scale = FALSE,
    do.center = TRUE,
    verbose = verbose
  )
  # Find anchors between query and reference
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    k.filter = NA,
    reference.neighbors = "refdr.annoy.neighbors",
    reference.assay = "refAssay",
    query.assay = "refAssay",
    reference.reduction = "refDR",
    normalization.method = "SCT",
    features = intersect(rownames(x = reference), VariableFeatures(object = query)),
    dims = 1:dims,
    n.trees = n.trees,
    mapping.score.k = mapping.score.k,
    verbose = verbose
  )
  # Transferred labels are in metadata columns named "predicted.*"
  # The maximum prediction score is in a metadata column named "predicted.*.score"
  # The prediction scores for each class are in an assay named "prediction.score.*"
  # The imputed assay is named "impADT" if computed
  refdata <- lapply(X = annotation.levels, function(x) {
    reference[[x, drop = TRUE]]
  })
  names(x = refdata) <- annotation.levels

  if (isTRUE(do.adt)) {
    refdata[["impADT"]] <- GetAssayData(
      object = reference[["ADT"]],
      slot = "data"
    )
  }

  query <- TransferData(
    reference = reference,
    query = query,
    dims = 1:dims,
    anchorset = anchors,
    refdata = refdata,
    n.trees = 20,
    store.weights = TRUE,
    k.weight = k.weight,
    verbose = verbose
  )
  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = reference,
    query = query,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE,
    verbose = verbose
  )
  # Calculate the query neighbors in the reference
  # with respect to the integrated embeddings
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(reference[["refDR"]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE,
    verbose = verbose
  )
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  query <- NNTransform(
    object = query,
    meta.data = reference[[]]
  )
  # Project the query to the reference UMAP.
  query[[umap.name]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = reference[["refUMAP"]],
    reduction.key = 'UMAP_',
    verbose = verbose
  )
  # Calculate mapping score and add to metadata
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors, ndim = dims),
    col.name = "mapping.score"
  )
  return(query)
}


#' @inheritParams RunAzimuth
#' @param reference Name of reference to map to or a path to a directory containing ref.Rds bridge.Rds and ext.Rds
#' @param annotation.levels list of annotation levels to map. If not specified, all will be mapped.
#' @param umap.name name of umap reduction in the returned object
#' @param do.adt transfer ADT assay
#' @param assay query assay name
#'
#' @return Seurat object with reference reductions and annotations
#'
#' @importFrom SeuratData InstallData InstalledData LoadData AvailableData
#' @importFrom Signac CollapseToLongestTranscript GetTranscripts GetGRangesFromEnsDb RunTFIDF # issue is i think the transcript ones are internal functions so would this work
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom IRanges findOverlaps
#' @importFrom Seurat FindBridgeTransferAnchors MapQuery NormalizeData
#' @importFrom data.table as.data.table
#' @export
#' @method RunAzimuth Bridge
#' @rdname RunAzimuth.Bridge
#'
RunAzimuth.Bridge <- function(
    query,
    reference,
    annotation.levels = NULL,
    umap.name = "ref.umap",
    do.adt = FALSE,
    verbose = TRUE,
    assay = "RNA",
    k.weight = 50,
    n.trees = 20,
    mapping.score.k = 100,
    dims.atac = 2:50, 
    dims.rna = 1:50
) {
  if (dir.exists(reference)) { 
    reference_all <- LoadBridgeReference(reference)
    reference <- reference_all$map
    multiome <- reference_all$bridge
    obj.rna.ext <- reference_all$ext
  } else {
    stop("Can't find path to reference")
  }
  #reference.version <- ReferenceVersion(reference)
  azimuth.version <- as.character(packageVersion(pkg = "Azimuth"))
  seurat.version <- as.character(packageVersion(pkg = "Seurat"))
  meta.data <- names(slot(reference, "meta.data"))
  
  # is annotation levels are not specify, gather all levels of annotation
  if (is.null(annotation.levels)) {
    annotation.levels <- names(slot(object = reference, name = "meta.data"))
    annotation.levels <- annotation.levels[!grepl(pattern = "^nCount", x = annotation.levels)]
    annotation.levels <- annotation.levels[!grepl(pattern = "^nFeature", x = annotation.levels)]
    annotation.levels <- annotation.levels[!grepl(pattern = "^ori", x = annotation.levels)]
  }
  if (file.exists(query)) {
    query_counts <- LoadFileInput(path = query) 
  } else {
    stop("Can't find path to query")
  }
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # this could also be saved as an object if we wanna save 30 s 
  seqlevelsStyle(annotation) <- 'UCSC'
  print("got annotations")
  query_assay <- CreateChromatinAssay(
    counts = query_counts[["RNA"]]@counts,
    sep = c(":", "-"),
    annotation = annotation
  )
  print("made chromatin assay")
  print(query_assay)
  o_hits <- findOverlaps(query_assay, multiome[["ATAC"]])
  print("found overlaps")
  query_requantified  <- RequantifyPeaks(o_hits, query_assay, multiome)
  print("requantified peaks")
  # Create assay with requantified ATAC data
  ATAC_assay <- CreateChromatinAssay(
    counts = query_requantified,
    sep = c(":", "-"),
    annotation = annotation
  )
  
  # Create Seurat Object
  obj.atac <- CreateSeuratObject(counts = ATAC_assay, assay = 'ATAC')
  obj.atac[['peak.orig']] <- query_assay
  obj.atac <- subset(obj.atac, subset = nCount_ATAC < 7e4 & nCount_ATAC > 2000)
  
  # normalize query
  obj.atac <- RunTFIDF(obj.atac)
  
  # Find anchors between query and reference # deleted find transfer anchors 
  bridge.anchor <- FindBridgeTransferAnchors(extended.reference = obj.rna.ext,
                                             query = obj.atac,
                                             reduction = "lsiproject",
                                             dims = dims.atac)
  # Transferred labels are in metadata columns named "predicted.*"
  # The maximum prediction score is in a metadata column named "predicted.*.score"
  # The prediction scores for each class are in an assay named "prediction.score.*"
  # The imputed assay is named "impADT" if computed
  refdata <- lapply(X = annotation.levels, function(x) {
    reference[[x, drop = TRUE]]
  })
  names(x = refdata) <- annotation.levels
  
  #if (isTRUE(do.adt)) {
  #  refdata[["impADT"]] <- GetAssayData(
  #    object = reference[["ADT"]],
  #    slot = "data"
  #  )
  #}
  
  obj.atac <- MapQuery(anchorset = bridge.anchor,  # deleted transfer data 
                       reference = reference, 
                       query = obj.atac, 
                       refdata = refdata,
                       reduction.model = "refUMAP")
  # Get Gene Activities 
  transcripts <- GetTranscripts(obj.atac)
  o_hits <- findOverlaps(obj.atac[["ATAC"]], transcripts)
  atac_final <- RequantifyPeaks(o_hits, obj.atac, transcripts)
  #dd feature matrix to Chromatin Assay 
  obj.atac[['RNA']] <- CreateAssayObject(counts = atac_final)
  
  #Normalize the feature data
  obj.atac <- NormalizeData(
    object = obj.atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(obj.atac$nCount_RNA)
  )
  return(obj.atac)
}


#' @inheritParams RunAzimuth
#' @export
#' @method RunAzimuth character
#' @rdname RunAzimuth
#'
RunAzimuth.character <- function(
  query,
  ...
) {
  obj <- LoadFileInput(path = query)
  return(RunAzimuth(obj, ...))
}

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
AzimuthApp <- function(config = NULL, type = "standard", ...) {
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
  # Launch the app
  if (type == "standard"){
    with_options(
      new = opts,
      code = runApp(appDir = shinyApp(ui = AzimuthUI, server = AzimuthServer))
    )
  } else{
    with_options(
      new = opts,
      code = runApp(appDir = shinyApp(ui = AzimuthBridgeUI, server = AzimuthBridgeServer))
    )
  }
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
# @param type Standard or Bridge
# @return Returns
#
#' @importFrom SeuratObject Cells Idents Indices as.Neighbor
#' @importFrom Seurat RunPCA FindNeighbors FindClusters MinMax
#' @importFrom Signac RunSVD
#
#' @keywords internal
#
#
ClusterPreservationScore <- function(query, ds.amount, type = "standard") {
  dims <- min(50, getOption(x = "Azimuth.map.ndims"))
  if(type == "standard"){
    query <- DietSeurat(object = query, assays = "refAssay", scale.data = TRUE, counts = FALSE, dimreducs = "integrated_dr")
    if (ncol(x = query) > ds.amount) {
      query <- subset(x = query, cells = sample(x = Cells(x = query), size = ds.amount))
    }
    query <- RunPCA(object = query, npcs = dims, verbose = FALSE)
    query <- FindNeighbors(
      object = query,
      reduction = 'pca',
      dims = 1:dims,
      graph.name = paste0("pca_", c("nn", "snn"))
    )
    query[["orig_neighbors"]] <- as.Neighbor(x = query[["pca_nn"]])
    query <- FindClusters(object = query, resolution = 0.6, graph.name = 'pca_snn')
    query <- FindNeighbors(
      object = query,
      reduction = 'integrated_dr',
      dims = 1:dims,
      return.neighbor = TRUE,
      graph.name ="integrated_neighbors_nn"
    )
  } else if(type == "bridge") {
    query <- DietSeurat(object = query, assays = "refAssay", scale.data = TRUE, counts = FALSE, dimreducs = "ref.Bridge.reduc")
    if (ncol(x = query) > ds.amount) {
      query <- subset(x = query, cells = sample(x = Cells(x = query), size = ds.amount))
    }
    query <- RunSVD(object = query, verbose = FALSE)
    query <- FindNeighbors(
      object = query,
      reduction = 'lsi',
      dims = 1:dims,
      graph.name = paste0("lsi_", c("nn", "snn"))
      )
    query[["orig_neighbors"]] <- as.Neighbor(x = query[["lsi_nn"]])
    query <- FindClusters(object = query, resolution = 0.6, graph.name = 'lsi_snn')
    query <- FindNeighbors(
      object = query,
      reduction = 'ref.Bridge.reduc',
      dims = 1:dims,
      return.neighbor = TRUE,
      graph.name ="integrated_neighbors_nn")
    } else{
      print("Incorrect type: Must be either 'standard' or 'bridge'")
    }
  ids <- Idents(object = query)
  integrated.neighbor.indices <- Indices(object = query[["integrated_neighbors_nn"]])
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

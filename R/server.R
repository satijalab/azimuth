#' @include zzz.R
#' @include helpers.R
NULL

#' Server function for the mapping app
#'
#' @param input,output,session Required Shiny app server parameters
#'
#' @return The shiny server logic
#'
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom data.table as.data.table
#' @importFrom DT dataTableProxy renderDT selectRows
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom future future plan resolved value
#' @importFrom GenomeInfoDb seqlevelsStyle StandardChromosomes
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom ggplot2 annotate geom_hline ggtitle scale_colour_hue
#' theme_void xlab layer_scales xlim ylim ggplot aes geom_point theme
#' element_blank element_rect labs
#' @importFrom googlesheets4 gs4_auth gs4_get sheet_append
#' @importFrom IRanges findOverlaps
#' @importFrom JASPAR2020 JASPAR2020
#' @importFrom methods slot slot<- new
#' @importFrom presto wilcoxauc
#' @importFrom SeuratObject AddMetaData Assays Cells DefaultAssay Embeddings
#' GetAssayData Idents Idents<- Key RenameCells Reductions Tool SetAssayData
#' VariableFeatures
#' @importFrom Seurat DimPlot FeaturePlot FindNeighbors FindTransferAnchors
#' IntegrateEmbeddings MappingScore NoLegend PercentageFeatureSet
#' RunUMAP TransferData SCTransform VlnPlot LabelClusters
#' FindBridgeTransferAnchors MapQuery NormalizeData
#' @importFrom Signac CollapseToLongestTranscript GetTranscripts GetGRangesFromEnsDb RunTFIDF 
#' Runmotif AddMotifs
#' @importFrom shiny downloadHandler observeEvent isolate Progress
#' reactiveValues renderPlot renderTable renderText removeUI setProgress
#' safeError updateNumericInput updateSelectizeInput updateCheckboxInput updateTextAreaInput
#' withProgress renderUI onStop showNotification wellPanel nearPoints insertUI
#' modalDialog showModal getDefaultReactiveDomain
#' @importFrom shinydashboard menuItem renderMenu renderValueBox
#' sidebarMenu valueBox
#' @importFrom shinyjs addClass enable disable hide removeClass show onclick
#' disable
#' @importFrom stringr str_interp str_trim str_split
#' @importFrom patchwork wrap_plots
#' @importFrom stats na.omit quantile setNames median
#' @importFrom TFBSTools getMatrixSet
#' @importFrom utils write.table packageVersion
#' @importFrom plotly plotlyOutput renderPlotly toWebGL ggplotly plot_ly
#'
#' @keywords internal
#'
AzimuthServer <- function(input, output, session) {
  hide(id = "legend")
  disable(id = 'metacolor.ref')
  # hide demo dataset button if required
  if (is.null(x = getOption(x = 'Azimuth.app.demodataset'))) {
    hide(id = "demobuttons")
  }
  do.adt <- isTRUE(x = as.logical(getOption(x = 'Azimuth.app.do_adt'), default = TRUE))
  adt.key <- 'impADT'
  # Do Bridge Integration Workflow for ATAC query
  do.bridge <- isTRUE(x = as.logical(getOption(x = 'Azimuth.app.do_bridge'), default = FALSE))
  mt.key <- 'percent.mt'
  mito.pattern <- getOption(x = 'Azimuth.app.mito', default = '^MT-')
  
  n.trees <- getOption(x = "Azimuth.map.ntrees")
  app.env <- reactiveValues(
    adt.features = character(length = 0L),
    anchors = NULL,
    annotations = NULL,
    bridge_anchors = FALSE,
    chromatin_assay_1 = NULL,
    chromatin_assay_2 = NULL,
    motif.diff.expr = list(),
    motif.feature = "", 
    motif.features = character(length = 0L),
    clusterpreservationqc = NULL,
    counts = FALSE,
    demo = FALSE,
    demo.inputs = NULL,
    demo.tracker = NULL,
    demo.files = NULL,
    default.assay = NULL,
    default.motif.feature = NULL,
    default.feature = NULL,
    default.metadata = NULL,
    diff.exp = list(),
    feature = '',
    features = character(length = 0L),
    mapping.score = NULL,
    messages = 'Upload a file',
    nanchors = 0L,
    ncellsupload = 0L,
    ncellspreproc = 0L,
    object = NULL,
    metadata.cont = character(length = 0L),
    scorefeatures = character(length = 0L),
    plot.ranges = list(),
    plots.refdim_df = NULL,
    plots.refdim_intro_df = NULL,
    plots.objdim_df = NULL,
    plots.querydim_df = NULL,
    fresh.plot = TRUE,
    singlepred = NULL,
    emptyref = NULL,
    merged = NULL,
    metadata.discrete = NULL,
    metadata.notransfer = NULL,
    requantified_multiome = NULL, 
    requantified_genes = NULL,
    disable = FALSE,
    query.names = character()
  )
  react.env <- reactiveValues(
    no = FALSE,
    anchors = FALSE,
    annotations = FALSE,
    biomarkers = FALSE,
    bridge = FALSE,
    bridge.query = FALSE,
    bridge_anchors = FALSE,
    motif = FALSE, 
    motif.features = FALSE,
    chromatin_assay_1 = FALSE,
    cluster.score = FALSE,
    features = FALSE,
    get.motif.feature = FALSE,
    get.feature = FALSE,
    map = FALSE,
    markers = FALSE,
    metadata = FALSE,
    mt = NULL,
    xferopts = FALSE,
    path = NULL,
    progress = NULL,
    plot.qc = FALSE,
    qc = FALSE,
    requantify_multiome = FALSE, 
    requantify_genes = FALSE,
    score = FALSE,
    sctransform = FALSE,
    start = numeric(length = 0L),
    transform = FALSE
  )
  if (isTRUE(x = do.adt)) {
    output$imputedlabel <- renderUI(expr = h3('Imputed protein biomarkers'))
  } else {
    for (id in c('imputedinput', 'imputedtable', 'imputeddl')) {
      removeUI(selector = paste0('#', id), immediate = TRUE)
    }
    for (id in c('featureinput', 'scoreinput')) {
      removeClass(id = id, class = 'thirds')
      addClass(id = id, class = 'halves')
    }
    removeClass(id = 'biotable', class = 'halves')
    addClass(id = 'biotable', class = 'fulls')
  }
  if (!isTRUE(x = do.bridge)) {
    print("removing ui elements")
    for (id in c('dist.qc', 'q4', 'valuebox.overlap', 'valuebox.jaccard', 'motifinput', 'continput.motif', 'metagroup.motif', 'motifvln', 'markerclustersgroupinput.motif', 'motiftable', 'overlap_box')) {
      removeUI(selector = paste0('#', id), immediate = TRUE)
    }
  }
  ResetEnv <- function() {
    print('resetting...')
    app.env$disable <- TRUE
    output$menu2 <- NULL
    react.env$plot.qc <- FALSE
    app.env$messages <- NULL
    output$valubox.jaccard <- NULL
    output$valubox.upload <- NULL
    output$valuebox.preproc <- NULL
    output$valuebox.mapped <- NULL
    output$valubox.overlap <- NULL
    output$valuebox_panchors <- NULL
    output$valuebox_mappingqcstat <- NULL
    app.env$emptyref <- NULL
    app.env$merged <- NULL
    app.env$metadata.discrete <- NULL
    disable(id = 'map')
    hide(selector = '.rowhide')
  }
  motif.proxy <- dataTableProxy(outputId = "motifs")
  rna.proxy <- dataTableProxy(outputId = 'biomarkers')
  adt.proxy <- dataTableProxy(outputId = 'adtbio')
  logging <- all(vapply(
    X = paste0('Azimuth.app.google', c('sheet', 'token', 'tokenemail')),
    FUN = function(x) {
      return(!is.null(x = getOption(x = x)))
    },
    FUN.VALUE = logical(length = 1L)
  ))
  googlesheet <- NULL
  if (logging) {
    try(
      expr = {
        gs4_auth(
          email = getOption(x = "Azimuth.app.googletokenemail"),
          cache = getOption(x = "Azimuth.app.googletoken")
        )
        googlesheet <- gs4_get(ss = getOption(x = "Azimuth.app.googlesheet"))
        app_start_time <- Sys.time()
        app_session_id <- paste0(Sys.info()[["nodename"]], as.numeric(Sys.time()))
        onStop(
          fun = function() {
            try(expr = sheet_append(
              ss = googlesheet,
              data = data.frame(
                "SESSIONLENGTH",
                app_session_id,
                as.numeric(x = Sys.time() - app_start_time, units = "mins")
              )
            ))
          }
        )
      }
    )
  }
  if (!is.null(x = googlesheet)) {
    try(expr = sheet_append(
      ss = googlesheet,
      data = data.frame(
        "STARTUPTIME",
        app_session_id,
        Sys.time()
      )
    ))
    output$menu3 <- renderMenu(expr = {
      sidebarMenu(menuItem(
        text = 'Feedback',
        tabName = 'tab_feedback',
        icon = icon(name = 'comments'),
        selected = FALSE
      ))
    })
  }
  demos <- getOption("Azimuth.app.demodataset")
  if (!is.null(x = demos)) {
    if (!inherits(x = demos, what = "data.frame")) {
      if (is.null(x = names(x = demos))) {
        if (length(x = demos) > 1) {
          demo.names <- paste0("Demo", 1:length(x = demos))
        } else {
          demo.names <- "Load demo dataset"
        }
      } else {
        demo.names <- names(x = demos)
      }
      demos <- data.frame(name = demo.names, file = demos)
    }
    app.env$demo.files <- demos$file
    app.env$demo.inputs <- paste0("triggerdemo", 1:nrow(x = demos))
    app.env$demo.tracker <- rep(x = 0, times = nrow(x = demos))
    for (i in 1:nrow(x = demos)) {
      insertUI(
        selector = '#demobuttons',
        where = 'beforeEnd',
        immediate = TRUE,
        ui = actionButton(
          inputId = paste0('triggerdemo', i),
          label = demos$name[i],
          width = '85%'
        )
      )
    }
  }
  if (getOption(x = "Azimuth.app.metatableheatmap")) {
    insertUI(
      selector = '#tablemetadata',
      where = 'beforeEnd',
      immediate = TRUE,
      ui = plotlyOutput(outputId = 'metadata.heatmap')
    )
  } else {
    insertUI(
      selector = '#tablemetadata',
      where = 'beforeEnd',
      immediate = TRUE,
      ui = tableOutput(outputId = 'metadata.table')
    )
  }
  if (getOption(x = "Azimuth.app.overlayedreference")) {
    insertUI(
      selector = "#topdim",
      where = "beforeEnd",
      immediate = TRUE,
      ui = box(
        title = 'Mapped Query',
        checkboxInput(inputId = 'legend', label = 'Show legend'),
        checkboxGroupInput(
          inputId = "label.opts", label = NULL,
          choiceNames = c("Show labels", "Filter cluster labels (size <2%)"),
          choiceValues = c("labels", "filterlabels"),
          selected = c("labels", "filterlabels"), inline = TRUE),
        checkboxInput(inputId = 'showrefonly', label = 'View reference only'),
        selectizeInput(
          inputId = 'metacolor.ref',
          label = 'Reference metadata to color by',
          choices = '',
          multiple = TRUE,
        ),
        selectizeInput(
          inputId = 'metacolor.query',
          label = 'Query metadata to color by',
          choices = '',
          multiple = TRUE,
        ),
        div(
          style = "position:relative",
          plotOutput(
            outputId = 'objdim',
            hover = hoverOpts(
              id = "objdim_hover_location",
              delay = 5,
              delayType = "debounce",
              nullOutside = TRUE
            ),
            height='750px'
          ),
          uiOutput("objdim_hover_box")
        ),
        width = 12,
        height = 'auto'
      )
    )
  } else {
    insertUI(
      selector = "#topdim",
      where = "beforeEnd",
      immediate = TRUE,
      ui = box(
        title = 'Reference',
        checkboxGroupInput(inputId = "dimplot.opts", label = NULL, choiceNames = c("Show labels", "Show legend"), choiceValues = c("labels", "legend"), selected = "legend", inline = TRUE),
        selectizeInput(
          inputId = 'metacolor.ref',
          label = 'Metadata to color by',
          choices = '',
          multiple = TRUE,
        ),
        div(
          style = "position:relative",
          plotOutput(
            outputId = 'refdim',
            hover = hoverOpts(
              id = "refdim_hover_location",
              delay = 5,
              delayType = "debounce",
              nullOutside = TRUE
            )
          ),
          uiOutput("refdim_hover_box")
        ),
        width = 12
      )
    )
    insertUI(
      selector = "#bottomdim",
      where = "beforeEnd",
      immediate = TRUE,
      ui = box(
        title = 'Query',
        selectizeInput(
          inputId = 'metacolor.query',
          label = 'Metadata to color by',
          choices = '',
          multiple = TRUE,
        ),
        div(
          style = "position:relative",
          plotOutput(
            outputId = 'querydim',
            hover = hoverOpts(
              id = "querydim_hover_location",
              delay = 5,
              delayType = "debounce",
              nullOutside = TRUE
            )
          ),
          uiOutput("querydim_hover_box")
        ),
        width = 12
      )
    )
  }
  if (isTRUE(x = do.bridge)) {
    withProgress(
      message = "Loading bridge and reference",
      expr = {
        print("bridge version")
        disable(id = 'file')
        ToggleDemos(action = "disable", demos = demos)
        setProgress(value = 0.2)
        refs <- LoadBridgeReference(
          path = getOption(
            x = 'Azimuth.app.reference',
            default = stop(safeError(error = "No reference provided"))
          )
        )
        setProgress(value = 1)
        enable(id = 'file')
        ToggleDemos(action = "enable", demos = demos)
        react.env$bridge <- TRUE
      }
    )
  } else {
    withProgress(
      message = "Loading reference",
      expr = {
        print("standard version")
        disable(id = 'file')
        ToggleDemos(action = "disable", demos = demos)
        setProgress(value = 0)
        refs <- LoadReference(
          path = getOption(
            x = 'Azimuth.app.reference',
            default = stop(safeError(error = "No reference provided"))
          )
        )
        setProgress(value = 1)
        enable(id = 'file')
        ToggleDemos(action = "enable", demos = demos)
        react.env$standard = TRUE
      }
    )
  }
  if (!is.null(x = googlesheet)) {
    try(
      expr = sheet_append(
        ss = googlesheet,
        data = data.frame(
          c('REFERENCE_NAME', "REFERENCE_VERSION"),
          c(app_session_id, app_session_id),
          c(basename(getOption(x = 'Azimuth.app.reference')),  ReferenceVersion(object = refs$map))
        )
      ),
      silent = TRUE
    )
  }
  plotseed <- getOption(x = "Azimuth.app.plotseed")
  if (!is.null(x = plotseed)) {
    set.seed(seed = plotseed)
    colormap <- GetColorMap(object = refs$map)
    for (i in names(x = colormap)) {
      names(x = colormap[[i]]) <- sample(x = names(x = colormap[[i]]))
    }
    refs$map <- SetColorMap(object = refs$map, value = colormap)
  }
  metadata.annotate <- names(x = GetColorMap(object = refs$map))
  if (!is.null(x = getOption(x = "Azimuth.app.metadata_notransfer", default = NULL))) {
    metadata.notransfer <- str_trim(
      str_split(
        getOption(x = "Azimuth.app.metadata_notransfer", default = NULL),
        ','
      )[[1]]
    )
    possible.metadata.transfer <- setdiff(x = metadata.annotate, y = metadata.notransfer)
  } else {
    possible.metadata.transfer <- metadata.annotate
  }
  if (length(x = possible.metadata.transfer) > 1) {
    react.env$xferopts <- TRUE
  }
  default_xfer <- getOption(x = "Azimuth.app.default_metadata", default = possible.metadata.transfer[1])
  if (!default_xfer %in% possible.metadata.transfer) {
    default_xfer <- possible.metadata.transfer[1]
  }
  # React to events
  # Load the data and prepare for QC
  observeEvent(
    eventExpr = input$file,
    handlerExpr = {
      ResetEnv()
      if (nchar(x = input$file$datapath)) {
        react.env$path <- input$file$datapath
      }
    }
  )
  observeEvent(
    eventExpr = sapply(X = app.env$demo.inputs, FUN = function(x) input[[x]]),
    handlerExpr = {
      if (isTRUE(x = !all(sapply(X = app.env$demo.inputs, FUN = is.null)))) {
        ResetEnv()
        for (i in 1:length(x = app.env$demo.inputs)) {
          if (isTRUE(x = input[[app.env$demo.inputs[i]]] != app.env$demo.tracker[i])) {
            app.env$demo.tracker[i] <- app.env$demo.tracker[i] + 1
            react.env$path <- app.env$demo.files[i]
          }
        }
      }
    },
    ignoreInit = TRUE
  )
  observeEvent(
    eventExpr = list(react.env$path, react.env$standard),
    handlerExpr = {
      if (!is.null(x = react.env$path) && nchar(x = react.env$path)) {
        if (isTRUE(react.env$standard)) {
          withProgress(
            message = 'Reading Input',
            expr = {
              setProgress(value = 0)
              tryCatch(
                expr = {
                  print("Doing standard path")
                  app.env$object <- LoadFileInput(path = react.env$path)
                  app.env$object <- DietSeurat(
                    app.env$object,
                    assays = "RNA"
                  )
                  app.env$object <- ConvertGeneNames(
                    object = app.env$object,
                    reference.names = rownames(x = refs$map),
                    homolog.table = getOption(x = 'Azimuth.app.homologs')
                  )
                  if (react.env$path %in% app.env$demo.files) {
                    app.env$demo <- TRUE
                  } else {
                    app.env$demo <- FALSE
                  }
                  app.env$object$query <- 'query'
                  Idents(object = app.env$object) <- 'query'
                  # check that no names overlap with reference
                  query.cell.names <- paste0("query", 1:ncol(x = app.env$object))
                  while (any(query.cell.names %in% Cells(x = refs$map))) {
                    query.cell.names <- paste0(query.cell.names, "x")
                  }
                  app.env$query.names <- Cells(x = app.env$object)
                  app.env$object <- RenameCells(object = app.env$object, new.names = query.cell.names)
                  
                  app.env$default.assay <- DefaultAssay(object = app.env$object)
                  new.mt <- any(grepl(
                    pattern = mito.pattern,
                    x = rownames(x = app.env$object)
                  ))
                  if (isFALSE(x = new.mt) & !isFALSE(x = react.env$mt)) {
                    removeUI(selector = '#pctmt', immediate = TRUE)
                  } else if (!isFALSE(x = new.mt) & isFALSE(x = react.env$mt)) {
                    insertUI(
                      selector = '#nfeature',
                      where = 'afterEnd',
                      immediate = TRUE,
                      ui = div(
                        id = 'pctmt',
                        numericInput(
                          inputId = 'minmt',
                          label = NULL,
                          value = 0,
                          width = '90%'
                        ),
                        numericInput(
                          'maxmt',
                          label = NULL,
                          value = 0,
                          width = '90%'
                        )
                      )
                    )
                  }
                  react.env$mt <- new.mt
                  common.features <- intersect(
                    x = rownames(x = app.env$object),
                    y = rownames(x = refs$map)
                  )
                  reject <- c(
                    length(x = common.features) < getOption(x = 'Azimuth.map.ngenes'),
                    length(x = Cells(x = app.env$object)) > getOption(x = 'Azimuth.app.max_cells')
                  )
                  if (any(reject)) {
                    app.env$object <- NULL
                    gc(verbose = FALSE)
                    reject <- min(which(x = reject))
                    print("GOING TO REJECT")
                    app.env$messages <- paste(
                      c(
                        'Not enough genes in common with reference.',
                        'Too many cells.'
                      ),
                      'Try another dataset.'
                    )[reject]
                  }
                  if (isFALSE(x = react.env$xferopts)) {
                    removeUI(selector = '#xferopts', immediate = TRUE)
                  }
                  react.env$qc <- !any(reject)
                  react.env$path <- NULL
                },
                error = function(e) {
                  app.env$messages <- e$message
                  showNotification(
                    e$message,
                    duration = 10,
                    type = 'error',
                    closeButton = TRUE,
                    id = 'no-progress-notification'
                  )
                  app.env$object <- NULL
                  gc(verbose = FALSE)
                  react.env$path <- NULL
                }
              )
              setProgress(value = 1)
            }
          )
        }   
      }
    }
  )
  observeEvent(
    eventExpr = list(react.env$path, react.env$bridge),
    handlerExpr = {
      if (!is.null(x = react.env$path) && nchar(x = react.env$path)) {
        if (isTRUE(react.env$bridge)) {
          withProgress(
            message = 'Reading ATAC Peaks',
            expr = {
              setProgress(value = 0)
              tryCatch(
                expr = {
                  app.env$counts <- LoadFileInput(path = react.env$path, 
                                                  bridge = TRUE)
                  
                  app.env$counts <- DietSeurat(
                    app.env$counts,
                    assays = "RNA"
                  )
                  # app.env$object <- ConvertGeneNames(
                  #   object = app.env$object,
                  #   reference.names = rownames(x = refs$map),
                  #   homolog.table = getOption(x = 'Azimuth.app.homologs')
                  # )
                  if (react.env$path %in% app.env$demo.files) {
                    app.env$demo <- TRUE
                  } else {
                    app.env$demo <- FALSE
                  }
                  app.env$counts$query <- 'query'
                  react.env$chromatin_assay_1 <- TRUE
                  react.env$path <- NULL
                },
                error = function(e) {
                  app.env$messages <- e$message
                  showNotification(
                    e$message,
                    duration = 10,
                    type = 'error',
                    closeButton = TRUE,
                    id = 'no-progress-notification'
                  )
                  app.env$object <- NULL
                  gc(verbose = FALSE)
                  react.env$path <- NULL
                }
              )
              setProgress(value = 1)
            }
          )
        }
      }
    }
  )
  observeEvent(
    eventExpr = react.env$chromatin_assay_1, 
    handlerExpr = {
      if (isTRUE(x = react.env$chromatin_assay_1)) {
        withProgress(message = "Making Chromatin Assay", expr = {
          setProgress(value = 0.2)
          tryCatch(expr = {
            print("about to make chromatin assay")
            app.env$annotations <- refs$annotation
            app.env$chromatin_assay_1 <- CreateChromatinAssay(
              counts = app.env$counts[["RNA"]]@counts, # this should probably be clearer 
              sep = c(":", "-"),
              annotation = app.env$annotations
            )
            print(app.env$chromatin_assay_1)
            print(refs$bridge)
            #app.env$o_hits <- findOverlaps(app.env$chromatin_assay_1, refs$bridge[["ATAC"]])
            #qc_table <- OverlapQC(app.env$chromatin_assay_1, refs$bridge)
            perc_overlap <- round(x = OverlapTotal(app.env$chromatin_assay_1, refs$bridge[["ATAC"]]), digits = 4)
            print("PERC OVERLAP")
            print(perc_overlap)
            if (perc_overlap >= 80) {
              output$valuebox.overlap <- renderValueBox(expr = {
                valueBox(value = perc_overlap, subtitle = "Overlap Percentage",
                         icon = icon(name = "check"), color = "green")
              })
            }
            else if (perc_overlap < 80 & perc_overlap > 60) {
              output$valuebox.overlap<- renderValueBox(expr = {
                valueBox(value = perc_overlap, subtitle = "Overlap Percentage",
                         icon = icon(name = "exclamation-circle"), color = "yellow")
              })
              
            }
            else {
              output$valuebox.overlap <- renderValueBox(expr = {
                valueBox(value = perc_overlap, subtitle = "Overlap Percentage Too Low",
                         icon = icon(name = "exclamation-circle"), color = "red")
              })
            }
            jaccard <- round(x = PeakJaccard(app.env$chromatin_assay_1, refs$bridge[["ATAC"]]), digits = 4)
            print("JACCARD SIMILARITY")
            print(jaccard)
            if (jaccard >= 50) {
              output$valuebox.jaccard <- renderValueBox(expr = {
                valueBox(value = jaccard, subtitle = "Jaccard Similarity",
                         icon = icon(name = "check"), color = "green")
              })
            }
            else if (perc_overlap < 50 & perc_overlap > 20) {
              output$valuebox.jaccard<- renderValueBox(expr = {
                valueBox(value = jaccard, subtitle = "Jaccard Similarity",
                         icon = icon(name = "exclamation-circle"), color = "yellow")
              })
              
            }
            else {
              output$valuebox.jaccard <- renderValueBox(expr = {
                valueBox(value = jaccard, subtitle = "Jaccard Similarity is Low",
                         icon = icon(name = "exclamation-circle"), color = "red")
              })
            }
            print("made chromatin assay ")
            query.cell.names <- paste0("query", 1:ncol(x = app.env$chromatin_assay_1))
            print("got query cell names")
            head(query.cell.names)
            while (any(query.cell.names %in% Cells(x = refs$map))) {
              query.cell.names <- paste0(query.cell.names, 
                                         "x")
            }
            print("assessed query cell names")
            app.env$query.names <- Cells(x = app.env$chromatin_assay_1)
            print("QUERY CELLS")
            print(length(app.env$query.names))
            print("got cells ")
            app.env$chromatin_assay_1 <- RenameCells(object = app.env$chromatin_assay_1, 
                                                     new.names = query.cell.names)
            print("renamed cells")
            
            # remove this because we don't have mitochondrial genes, just peaks 
            removeUI(selector = '#pctmt', immediate = TRUE)
            react.env$mt <- FALSE
            
            react.env$requantify_multiome <- TRUE
            react.env$chromatin_assay_1 <- FALSE
          }, error = function(e) {
            app.env$messages <- e$message
            showNotification(e$message, duration = 10, 
                             type = "error", closeButton = TRUE, id = "no-progress-notification")
            app.env$chromatin_assay_1 <- NULL
            gc(verbose = FALSE)
            react.env$chromatin_assay_1 <- NULL
          })
          setProgress(value = 0.3)
        }
        )
      }
    })
  observeEvent(
    eventExpr = react.env$requantify_multiome, 
    handlerExpr = {
      if (isTRUE(x = react.env$requantify_multiome)) {
        withProgress(message = "Requantifying Peaks to Match Bridge", expr = {
          print("DOING REQUANTIFICATION")
          setProgress(value = 0.6)
          tryCatch(expr = {
            app.env$requantified_multiome <- RequantifyPeaks(app.env$chromatin_assay_1, refs$bridge)
            app.env$chromatin_assay_2 <- CreateChromatinAssay(
              counts = app.env$requantified_multiome,
              sep = c(":", "-"),
              annotation = app.env$annotations
            )
            app.env$object <- CreateSeuratObject(counts = app.env$chromatin_assay_2, assay = 'ATAC')
            app.env$object[['peak.orig']] <- app.env$chromatin_assay_1
            app.env$object$query <- "query"
            app.env$default.assay <- DefaultAssay(app.env$object)
            
            common.features <- intersect(
              x = rownames(x = app.env$object),
              y = rownames(x = refs$bridge[["ATAC"]])
            )
            print(length(rownames(x = app.env$object)))
            head(row.names(app.env$object))
            print(length(rownames(x = refs$bridge[["ATAC"]])))
            head(rownames(refs$bridge[["ATAC"]]))
            print(length(common.features))
            reject_peaks <- c(
              length(x = common.features) < getOption(x = 'Azimuth.map.ngenes'),
              length(x = Cells(x = app.env$object)) > getOption(x = 'Azimuth.app.max_cells')
            )
            if (any(reject_peaks)) {
              print("GOING TO REJECT 2")
              app.env$object <- NULL
              gc(verbose = FALSE)
              reject_peaks <- min(which(x = reject_peaks))
              app.env$messages <- paste(
                c(
                  'Not enough peaks in common with reference.',
                  'Too many cells.'
                ),
                'Try another dataset.'
              )[reject_peaks]
            }
            if (isFALSE(x = react.env$xferopts)) {
              removeUI(selector = '#xferopts', immediate = TRUE)
            }
            
            react.env$qc <- !any(reject_peaks)
            react.env$requantify_multiome <- FALSE
          }, error = function(e) {
            app.env$messages <- e$message
            showNotification(e$message, duration = 10, 
                             type = "error", closeButton = TRUE, id = "no-progress-notification")
            app.env$chromatin_assay_2 <- NULL
            gc(verbose = FALSE)
            react.env$requantify_multiome <- NULL
          })
          setProgress(value = 1)
        }
        )
      }
    })
  observeEvent(
    eventExpr = react.env$qc,
    handlerExpr = {
      if (isTRUE(x = react.env$qc)) {
        for (id in qc.ids) {
          try(expr = enable(id = id), silent = TRUE)
        }
        ncount <- paste0('nCount_', app.env$default.assay)
        nfeature <- paste0('nFeature_', app.env$default.assay)
        if (!all(c(ncount, nfeature) %in% colnames(x = app.env$object[[]]))) {
          withProgress(
            message = 'Calculating nCount and nFeature',
            expr = {
              setProgress(value = 0)
              calcn <- as.data.frame(x = Seurat:::CalcN(object = app.env$object))
              colnames(x = calcn) <- paste(
                colnames(x = calcn),
                app.env$default.assay,
                sep = '_'
              )
              app.env$object <- AddMetaData(
                object = app.env$object,
                metadata = calcn
              )
              rm(calcn)
              gc(verbose = FALSE)
              setProgress(value = 1)
            }
          )
        }
        ncount.val <- range(app.env$object[[ncount, drop = TRUE]])
        ncount.val <- c(
          floor(x = min(ncount.val)),
          ceiling(x = max(ncount.val))
        )
        ncount.min <- if (is.null(getOption(x = "Azimuth.app.ncount_min"))) {
          ncount.val[1]
        } else {
          max(ncount.val[1], getOption(x = "Azimuth.app.ncount_min"))
        }
        ncount.max <- if (is.null(getOption(x = "Azimuth.app.ncount_max"))) {
          ncount.val[2]
        } else {
          min(ncount.val[2], getOption(x = "Azimuth.app.ncount_max"))
        }
        updateNumericInput(
          session = session,
          inputId = 'num.ncountmin',
          label = paste('min', ncount),
          value = ncount.min,
          min = ncount.val[1],
          max = ncount.val[2]
        )
        updateNumericInput(
          session = session,
          inputId = 'num.ncountmax',
          label = paste('max', ncount),
          value = ncount.max,
          min = ncount.val[1],
          max = ncount.val[2]
        )
        nfeature.val <- range(app.env$object[[nfeature, drop = TRUE]])
        nfeature.val <- c(
          floor(x = min(nfeature.val)),
          ceiling(x = max(nfeature.val))
        )
        nfeature.min <- if (is.null(getOption(x = "Azimuth.app.nfeature_min"))) {
          nfeature.val[1]
        } else {
          max(nfeature.val[1], getOption(x = "Azimuth.app.nfeature_min"))
        }
        nfeature.max <- if (is.null(getOption(x = "Azimuth.app.nfeature_max"))) {
          nfeature.val[2]
        } else {
          min(nfeature.val[2], getOption(x = "Azimuth.app.nfeature_max"))
        }
        updateNumericInput(
          session = session,
          inputId = 'num.nfeaturemin',
          label = paste('min', nfeature),
          value = nfeature.min,
          min = nfeature.val[1],
          max = nfeature.val[2]
        )
        updateNumericInput(
          session = session,
          inputId = 'num.nfeaturemax',
          label = paste('max', nfeature),
          value = nfeature.max,
          min = nfeature.val[1],
          max = nfeature.val[2]
        )
        if (isTRUE(x = react.env$mt)) {
          app.env$object <- PercentageFeatureSet(
            object = app.env$object,
            pattern = mito.pattern,
            col.name = mt.key,
            assay = app.env$default.assay
          )
          mito.val <- range(app.env$object[[mt.key, drop = TRUE]])
          mito.val <- c(
            floor(x = min(mito.val)),
            ceiling(x = max(mito.val))
          )
          mito.min <- if (is.null(getOption(x = "Azimuth.app.pctmt_min"))) {
            mito.val[1]
          } else {
            max(mito.val[1], getOption(x = "Azimuth.app.pctmt_min"))
          }
          mito.max <- if (is.null(getOption(x = "Azimuth.app.pctmt_max"))) {
            mito.val[2]
          } else {
            min(mito.val[2], getOption(x = "Azimuth.app.pctmt_max"))
          }
          updateNumericInput(
            session = session,
            inputId = 'minmt',
            label = paste('min', mt.key),
            value = mito.min,
            min = mito.val[1],
            max = mito.val[2]
          )
          updateNumericInput(
            session = session,
            inputId = 'maxmt',
            label = paste('max', mt.key),
            value = mito.max,
            min = mito.val[1],
            max = mito.val[2]
          )
        }
        output$menu1 <- renderMenu(expr = {
          sidebarMenu(menuItem(
            text = 'Preprocessing',
            tabName = 'tab_preproc',
            icon = icon(name = 'filter'),
            selected = TRUE
          ))
        })
        ncellsupload <- length(x = colnames(x = app.env$object))
        app.env$ncellsupload <- ncellsupload
        app.env$messages <- paste(ncellsupload, 'cells uploaded')
        if (ncellsupload < getOption(x = 'Azimuth.map.ncells')) {
          output$valuebox.upload <- renderValueBox(expr = {
            valueBox(
              value = ncellsupload,
              subtitle = paste0(
                'cells uploaded - ',
                getOption(x = 'Azimuth.map.ncells'), ' required'
              ),
              icon = icon(name = 'times'),
              color = 'red'
            )
          })
        } else {
          output$valuebox.upload <- renderValueBox(expr = {
            valueBox(
              value = ncellsupload,
              subtitle = 'cells uploaded',
              icon = icon(name = 'check'),
              color = 'green'
            )
          })
          if (!is.null(x = googlesheet)) {
            try(
              expr = sheet_append(
                ss = googlesheet,
                data = data.frame(
                  'CELLSUPLOAD',
                  app_session_id,
                  ncellsupload
                )
              ),
              silent = TRUE
            )
          }
        }
        if (!is.null(x = react.env$progress)) {
          react.env$progress$close()
          enable(id = 'file')
          ToggleDemos(action = "enable", demos = demos)
          react.env$progress <- NULL
        }
        updateSelectizeInput(
          session = session,
          inputId = 'metadataxfer',
          choices = possible.metadata.transfer,
          selected = default_xfer,
          server = TRUE,
          options = selectize.opts[-which(x = names(x = selectize.opts) == 'maxItems')]
        )
        react.env$qc <- FALSE
        react.env$plot.qc <- TRUE
        if (isTRUE(x = do.bridge)) {
          react.env$dist.qc <- TRUE
        }
      }
    }
  )
  observeEvent(
    eventExpr = input$metadataxfer,
    handlerExpr = {
      if (length(x = input$metadataxfer) == 0) {
        disable(id = 'map')
      } else {
        enable(id = 'map')
      }
    },
    ignoreNULL = FALSE
  )
  # Filter and process the data
  observeEvent(
    eventExpr = input$map,
    handlerExpr = {
      react.env$start <- Sys.time()
      disable(id = 'file')
      ToggleDemos(action = "disable", demos = demos)
      for (id in qc.ids) {
        try(expr = disable(id = id), silent = TRUE)
      }
      react.env$progress <- Progress$new(style = 'notification')
      react.env$progress$set(
        value = 0,
        message = 'Filtering based on nCount and nFeature'
      )
      ncount <- paste0('nCount_', DefaultAssay(object = app.env$object))
      nfeature <- paste0('nFeature_', DefaultAssay(object = app.env$object))
      cells.use <- app.env$object[[ncount, drop = TRUE]] >= input$num.ncountmin &
        app.env$object[[ncount, drop = TRUE]] <= input$num.ncountmax &
        app.env$object[[nfeature, drop = TRUE]] >= input$num.nfeaturemin &
        app.env$object[[nfeature, drop = TRUE]] <= input$num.nfeaturemax
      if (isTRUE(x = react.env$mt)) {
        cells.use <- cells.use &
          app.env$object[[mt.key, drop = TRUE]] >= input$minmt &
          app.env$object[[mt.key, drop = TRUE]] <= input$maxmt
      }
      ncellspreproc <- sum(cells.use)
      app.env$ncellspreproc <- ncellspreproc
      # not enough cells available after filtering: reset filter elements
      if (ncellspreproc < getOption(x = "Azimuth.map.ncells")) {
        output$valuebox.preproc <- renderValueBox(expr = valueBox(
          value = ncellspreproc,
          subtitle = paste0(
            'cells after filtering - ',
            getOption(x = 'Azimuth.map.ncells'), ' required'
          ),
          icon = icon("times"),
          color = "red"
        ))
        react.env$qc <- TRUE
      } else {
        output$valuebox.preproc <- renderValueBox(expr = valueBox(
          value = ncellspreproc,
          subtitle = "cells after filtering",
          icon = icon("check"),
          color = "green"
        ))
        if (!is.null(googlesheet)) {
          try(sheet_append(
            ss = googlesheet,
            data = data.frame(
              "CELLSPREPROC",
              app_session_id,
              ncellspreproc
            )
          ))
        }
        app.env$object <- app.env$object[, cells.use]
        app.env$query.names <- app.env$query.names[cells.use]
        if (isTRUE(x = do.bridge)) {
          react.env$tfidf <- TRUE
        } else {
          react.env$sctransform <- TRUE
        }
      }
    }
  )
  observeEvent(
    eventExpr = react.env$sctransform,
    handlerExpr = {
      if (isTRUE(x = react.env$sctransform)) {
        react.env$progress$set(
          value = 0.2,
          message = 'Normalizing with SCTransform'
        )
        tryCatch(
          expr = {
            app.env$object <- suppressWarnings(expr = SCTransform(
              object = app.env$object,
              residual.features = rownames(x = refs$map),
              reference.SCT.model = slot(object = refs$map[["refAssay"]], name = "SCTModel.list")[["refmodel"]],
              method = "glmGamPoi",
              do.correct.umi = FALSE,
              do.scale = FALSE,
              do.center = TRUE,
              new.assay.name = "refAssay"
            ))
          },
          error = function(e) {
            app.env$object <- suppressWarnings(expr = SCTransform(
              object = app.env$object,
              residual.features = rownames(x = refs$map),
              reference.SCT.model = slot(object = refs$map[["refAssay"]], name = "SCTModel.list")[["refmodel"]],
              method = "poisson",
              do.correct.umi = FALSE,
              do.scale = FALSE,
              do.center = TRUE,
              new.assay.name = "refAssay"
            ))
          }
        )
        app.env$messages <- c(
          app.env$messages,
          paste(ncol(x = app.env$object), "cells preprocessed")
        )
        react.env$anchors <- TRUE
        react.env$sctransform <- FALSE
      }
    }
  )
  observeEvent(
    eventExpr = react.env$tfidf, 
    handlerExpr = {
      if (isTRUE(x = react.env$tfidf)) {
        react.env$progress$set(
          value = 0.2, 
          message = "Normalizing with TFIDF"
        )
        tryCatch(
          expr = {
            app.env$object <- suppressWarnings(expr = RunTFIDF(object = app.env$object,
                                                               method = 1))
          }, error = function(e) {
            app.env$object <- suppressWarnings(expr = RunTFIDF(object = app.env$object,
                                                               method = 1))
          })
        app.env$messages <- c(
          app.env$messages, 
          paste(ncol(x = app.env$object), "cells preprocessed")
        )
        react.env$bridge_anchors <- TRUE
        react.env$tfidf <- FALSE
      }
    }
  )
  observeEvent(
    eventExpr = react.env$anchors,
    handlerExpr = {
      if (isTRUE(x = react.env$anchors)) {
        react.env$progress$set(value = 0.3, message = 'Finding anchors')
        app.env$anchors <- FindTransferAnchors(
          reference = refs$map,
          query = app.env$object,
          k.filter = NA,
          reference.neighbors = "refdr.annoy.neighbors",
          reference.assay = "refAssay",
          query.assay = 'refAssay',
          reference.reduction = 'refDR',
          normalization.method = 'SCT',
          features = intersect(
            x = rownames(x = refs$map),
            y = VariableFeatures(object = app.env$object)
          ),
          dims = 1:getOption(x = "Azimuth.map.ndims"),
          n.trees = n.trees,
          verbose = TRUE,
          mapping.score.k = 100
        )
        nanchors <- nrow(x = slot(object = app.env$anchors, name = "anchors"))
        app.env$nanchors <- nanchors
        if (!is.null(googlesheet)) {
          try(sheet_append(
            ss = googlesheet,
            data = data.frame(
              "NANCHORS",
              app_session_id,
              nanchors
            )
          ))
        }
        if (nanchors < getOption(x = 'Azimuth.map.nanchors') |
            length(x = unique(x = slot(
              object = app.env$anchors, name = "anchors")[, 2]
            )) < 50
        ) {
          output$valuebox.mapped <- renderValueBox(expr = {
            valueBox(
              value = 'Failure',
              subtitle = paste0('Too few anchors identified (', nanchors, ')'),
              icon = icon(name = 'times'),
              color = 'red',
              width = 6
            )
          })
          app.env$object <- NULL
          app.env$anchors <- NULL
          react.env$progress$close()
          enable(id = 'file')
          ToggleDemos(action = "enable", demos = demos)
          gc(verbose = FALSE)
        } else {
          query.unique <- length(x = unique(x = slot(object = app.env$anchors, name = "anchors")[, "cell2"]))
          percent.anchors <- round(x = query.unique / ncol(x = app.env$object) * 100, digits = 2)
          if (percent.anchors <  getOption(x = "Azimuth.map.panchorscolors")[1]) {
            print("rendering value box 1")
            output$valuebox_panchors <- renderValueBox(expr = {
              valueBox(
                value = paste0(percent.anchors, "%"),
                subtitle = "% of query cells with anchors",
                color = 'red',
                icon = icon(name = 'times')
              )
            })
          } else if (percent.anchors <  getOption(x = "Azimuth.map.panchorscolors")[2]) {
            print("rendering value box 2")
            output$valuebox_panchors <- renderValueBox(expr = {
              valueBox(
                value = paste0(percent.anchors, "%"),
                subtitle = "% of query cells with anchors",
                color = 'yellow',
                icon = icon(name = 'exclamation-circle')
              )
            })
          } else {
            output$valuebox_panchors <- renderValueBox(expr = {
              valueBox(
                value = paste0(percent.anchors, "%"),
                subtitle = "% of query cells with anchors",
                color = 'green',
                icon = icon(name = 'check')
              )
            })
          }
          react.env$map <- TRUE
        }
        react.env$anchors <- FALSE
      }
    }
  )
  observeEvent(
    eventExpr = react.env$bridge_anchors, 
    handlerExpr = {
      if (isTRUE(x = react.env$bridge_anchors)) {
        react.env$progress$set(value = 0.3, message = "Finding anchors")
        app.env$anchors <- FindBridgeTransferAnchors(extended.reference = refs$ext,
                                                     query = app.env$object,
                                                     reduction = "lsiproject",
                                                     scale = FALSE,
                                                     dims = 2:50) # making this a default
        print("found anchors ")
        nanchors <- nrow(x = slot(object = app.env$anchors, 
                                  name = "anchors"))
        print("was able to get nanchors")
        print(nanchors)
        app.env$nanchors <- nanchors
        print("going onto googlesheet")
        if (!is.null(googlesheet)) {
          print("doing google sheet thing")
          try(sheet_append(ss = googlesheet, data = data.frame("NANCHORS", 
                                                               app_session_id, nanchors)))
        }
        if (nanchors < getOption(x = "Azimuth.map.nanchors") | 
            length(x = unique(x = slot(object = app.env$anchors, 
                                       name = "anchors")[, 2])) < 50) {
          print("doing this length if ")
          output$valuebox.mapped <- renderValueBox(expr = {
            valueBox(value = "Failure", subtitle = paste0("Too few anchors identified (", 
                                                          nanchors, ")"), icon = icon(name = "times"), 
                     color = "red", width = 6)
          })
          print("checking")
          app.env$object <- NULL
          app.env$anchors <- NULL
          react.env$progress$close()
          enable(id = "file")
          ToggleDemos(action = "enable", demos = demos)
          
          gc(verbose = FALSE)
        }
        else {
          print("doing the else")
          query.unique <- length(x = unique(x = slot(object = app.env$anchors, 
                                                     name = "anchors")[, "cell2"]))
          percent.anchors <- round(x = query.unique/ncol(x = app.env$object) * 
                                     100, digits = 2)
          if (percent.anchors < getOption(x = "Azimuth.map.panchorscolors")[1]) {
            print("rendering value box 1")
            output$valuebox_panchors <- renderValueBox(expr = {
              valueBox(value = paste0(percent.anchors, 
                                      "%"), subtitle = "% of query cells with anchors", 
                       color = "red", icon = icon(name = "times"))
            })
          }
          else if (percent.anchors < getOption(x = "Azimuth.map.panchorscolors")[2]) {
            print("rendering value box 2")
            output$valuebox_panchors <- renderValueBox(expr = {
              valueBox(value = paste0(percent.anchors, 
                                      "%"), subtitle = "% of query cells with anchors", 
                       color = "yellow", icon = icon(name = "exclamation-circle"))
            })
          }
          else {
            print("rendering value box 4")
            output$valuebox_panchors <- renderValueBox(expr = {
              valueBox(value = paste0(percent.anchors, 
                                      "%"), subtitle = "% of query cells with anchors", 
                       color = "green", icon = icon(name = "check"))
            })
          }
          print("finished anchors")
          react.env$mapquery <- TRUE
        }
        print("making bridge anchors false")
        react.env$bridge_anchors <- FALSE
      }
    }
  )
  observeEvent(
    eventExpr = list(react.env$map, input$metadataxfer),
    handlerExpr = {
      if (isTRUE(x = react.env$map)) {
        if (is.null(x = input$metadataxfer)) {
          app.env$metadataxfer <- names(x = GetColorMap(object = refs$map))
        } else {
          app.env$metadataxfer <- input$metadataxfer
        }
        react.env$progress$set(value = 0.5, message = 'Mapping cells')
        refdata <- lapply(X = app.env$metadataxfer, function(x) {
          refs$map[[x, drop = TRUE]]
        })
        names(x = refdata) <- app.env$metadataxfer
        if (do.adt) {
          refdata[["impADT"]] <- GetAssayData(
            object = refs$map[['ADT']],
            slot = 'data'
          )
        }
        app.env$object <- TransferData(
          reference = refs$map,
          query = app.env$object,
          dims = 1:getOption(x = "Azimuth.map.ndims"),
          anchorset = app.env$anchors,
          refdata = refdata,
          n.trees = n.trees,
          store.weights = TRUE
        )
        app.env$singlepred <- NULL
        for(i in app.env$metadataxfer) {
          app.env$singlepred <- c(app.env$singlepred, length(x = unique(x = as.vector(x = app.env$object[[paste0("predicted.", i), drop = TRUE]]))) == 1)
          app.env$object[[paste0("predicted.", i), drop = TRUE]] <- factor(
            x = app.env$object[[paste0("predicted.", i), drop = TRUE]],
            levels = levels(x = refs$map[[i, drop = TRUE]])
          )
        }
        singlepred <- all(app.env$singlepred)
        if (singlepred & (length(x = setdiff(possible.metadata.transfer, app.env$metadataxfer)) > 0)) {
          showNotification(
            paste0("Only one predicted class. Re-running with all metadata."),
            duration = 5,
            type = 'warning',
            closeButton = TRUE,
            id = 'no-progress-notification'
          )
          updateSelectizeInput(
            session = getDefaultReactiveDomain(),
            inputId = 'metadataxfer',
            choices = possible.metadata.transfer,
            selected = possible.metadata.transfer,
          )
          app.env$metadataxfer <- input$metadataxfer
        } else if (singlepred) {
          showNotification(
            paste0(
              "Only one predicted class: ",
              app.env$object[[paste0("predicted.", app.env$metadataxfer[1]), drop = TRUE]][1]
            ),
            duration = 5,
            type = 'warning',
            closeButton = TRUE,
            id = 'no-progress-notification'
          )
          app.env$object <- NULL
          app.env$anchors <- NULL
          react.env$path <- NULL
          react.env$map <- FALSE
          react.env$progress$close()
          enable(id = 'file')
          ToggleDemos(action = "enable", demos = demos)
          gc(verbose = FALSE)
        } else {
          app.env$object <- IntegrateEmbeddings(
            anchorset = app.env$anchors,
            reference = refs$map,
            query = app.env$object,
            reductions = "pcaproject",
            reuse.weights.matrix = TRUE
          )
          if (is.null(x = getOption(x = "Azimuth.app.default_metadata"))) {
            app.env$default.metadata <- names(x = refdata)[1]
          } else {
            if (getOption(x = "Azimuth.app.default_metadata") %in% names(x = refdata)) {
              app.env$default.metadata <- getOption(x = "Azimuth.app.default_metadata")
            } else {
              app.env$default.metadata <- names(x = refdata)[1]
            }
          }
          react.env$score <- TRUE
          # react.env$cluster.score <- TRUE
          react.env$map <- FALSE
        }
      }
    }
  )
  observeEvent(
    eventExpr = list(react.env$mapquery, input$metadataxfer), 
    handlerExpr = {
      if (isTRUE(x = react.env$mapquery)) {
        print("doing this mapquery")
        if (is.null(x = input$metadataxfer)) {
          app.env$metadataxfer <- names(x = GetColorMap(object = refs$map))
        }
        else {
          app.env$metadataxfer <- input$metadataxfer
        }
        print("mapping cells ")
        react.env$progress$set(value = 0.5, message = "Mapping cells")
        refdata <- lapply(X = app.env$metadataxfer, function(x) { 
          refs$map[[x, drop = TRUE]]
        })
        names(x = refdata) <- app.env$metadataxfer
        if (do.adt) {
          print("trying to do adt")
          refdata[["impADT"]] <- GetAssayData(object = refs$map[["ADT"]], 
                                              slot = "data")
        }
        app.env$object <-  MapQuery(anchorset = app.env$anchors,  # deleted transfer data 
                                    reference = refs$map, 
                                    query = app.env$object, 
                                    refdata = refdata,
                                    reduction.model = "refUMAP")
        print("finished mapquery")
        print(app.env$metadataxfer)
        app.env$singlepred <- NULL
        for (i in app.env$metadataxfer) { 
          app.env$singlepred <- c(app.env$singlepred, 
                                  length(x = unique(x = as.vector(x = app.env$object[[paste0("predicted.", 
                                                                                             i), drop = TRUE]]))) == 1)
          app.env$object[[paste0("predicted.", i), drop = TRUE]] <- factor(x = app.env$object[[paste0("predicted.", 
                                                                                                      i), drop = TRUE]], levels = levels(x = refs$map[[i, 
                                                                                                                                                       drop = TRUE]]))
        }
        singlepred <- all(app.env$singlepred)
        if (singlepred & (length(x = setdiff(possible.metadata.transfer, 
                                             app.env$metadataxfer)) > 0)) {
          showNotification(paste0("Only one predicted class. Re-running with all metadata."), 
                           duration = 5, type = "warning", closeButton = TRUE, 
                           id = "no-progress-notification")
          updateSelectizeInput(session = getDefaultReactiveDomain(), 
                               inputId = "metadataxfer", choices = possible.metadata.transfer, 
                               selected = possible.metadata.transfer, )
          app.env$metadataxfer <- input$metadataxfer
        }
        else if (singlepred) {
          showNotification(paste0("Only one predicted class: ", 
                                  app.env$object[[paste0("predicted.", app.env$metadataxfer[1]), 
                                                  drop = TRUE]][1]), duration = 5, type = "warning", 
                           closeButton = TRUE, id = "no-progress-notification")
          app.env$object <- NULL
          app.env$bridge_anchors <- NULL
          react.env$path <- NULL
          react.env$mapquery <- FALSE
          react.env$progress$close()
          enable(id = "file")
          ToggleDemos(action = "enable", demos = demos)
          gc(verbose = FALSE)
        }
        else {
          if (is.null(x = getOption(x = "Azimuth.app.default_metadata"))) {
            app.env$default.metadata <- names(x = refdata)[1]
          }
          else {
            if (getOption(x = "Azimuth.app.default_metadata") %in% 
                names(x = refdata)) {
              app.env$default.metadata <- getOption(x = "Azimuth.app.default_metadata")
            }
            else {
              app.env$default.metadata <- names(x = refdata)[1]
            }
          }
          #react.env$score <- TRUE - ill do this after getting gene activity scores 
          react.env$gene_activity <- TRUE
          react.env$mapquery <- FALSE
        }
      }
    }
  )
  observeEvent(
    eventExpr = react.env$gene_activity, 
    handlerExpr = {
      if (isTRUE(react.env$gene_activity)) {
        # Use original peaks 
        print("GENE ACTIVITY ")
        DefaultAssay(app.env$object) <- "peak.orig"
        print(app.env$object)
        startTime <- Sys.time()
        app.env$transcripts <- GetTranscripts(app.env$object)
        endTime <- Sys.time()
        print("TOTAL TIME TO GET TRANSCRIPTS")
        print(endTime - startTime)
        temp <- RequantifyPeaks(app.env$object, app.env$transcripts)
        #add feature matrix to Chromatin Assay 
        app.env$object[['RNA']] <- CreateAssayObject(counts = temp)
        
        #o_hits <- findOverlaps(app.env$object[["ATAC"]], app.env$transcripts)
        #temp <- RequantifyPeaks(o_hits, app.env$object, app.env$transcripts)
        #dd feature matrix to Chromatin Assay 
        #app.env$object[['RNA']] <- CreateAssayObject(counts = temp)
        
        #Normalize the feature data
        app.env$object <- NormalizeData(
          object = app.env$object,
          assay = 'RNA',
          normalization.method = 'LogNormalize',
          scale.factor = median(app.env$object$nCount_RNA)
        )
        print("feature data normalized")
        
        react.env$gene_activity <- FALSE
        react.env$score <- TRUE
      }
    }
  )
  observeEvent(
    eventExpr = react.env$cluster.score,
    handlerExpr = {
      if (isTRUE(react.env$cluster.score)) {
        # post mapping QC
        if (isTRUE(x = do.bridge)){
          qc.stat <- round(
            x = ClusterPreservationScore(
              query = app.env$object,
              ds.amount = getOption(x = "Azimuth.map.postmapqcds"),
              type = "bridge"
            ),
            digits = 2
          )
        } else {
          qc.stat <- round(
            x = ClusterPreservationScore(
              query = app.env$object,
              ds.amount = getOption(x = "Azimuth.map.postmapqcds"),
              type = "standard"
            ),
            digits = 2
          )
        }
        if (!is.null(googlesheet)) {
          try(sheet_append(
            ss = googlesheet,
            data = data.frame(
              "CLUSTERPRESERVATIONQC",
              app_session_id,
              qc.stat
            )
          ))
        }
        app.env$clusterpreservationqc <- qc.stat
        if (qc.stat <  getOption(x = "Azimuth.map.postmapqccolors")[1]) {
          output$valuebox_mappingqcstat <- renderValueBox(expr = {
            valueBox(
              value = paste0(qc.stat, "/5"),
              subtitle = "cluster preservation score",
              color = 'red',
              icon = icon(name = 'times')
            )
          })
        } else if (qc.stat <  getOption(x = "Azimuth.map.postmapqccolors")[2]) {
          output$valuebox_mappingqcstat <- renderValueBox(expr = {
            valueBox(
              value = paste0(qc.stat, "/5"),
              subtitle = "cluster preservation score",
              color = 'yellow',
              icon = icon(name = 'exclamation-circle')
            )
          })
        } else {
          output$valuebox_mappingqcstat <- renderValueBox(expr = {
            valueBox(
              value = paste0(qc.stat, "/5"),
              subtitle = "cluster preservation score",
              color = 'green',
              icon = icon(name = 'check')
            )
          })
        }
        react.env$cluster.score <- FALSE
        react.env$transform <- TRUE
      }
    }
  )
  observeEvent(
    eventExpr = react.env$score,
    handlerExpr = {
      if (isTRUE(x = react.env$score)) {
        react.env$progress$set(
          value = 0.7,
          message = 'Calculating mapping score'
        )
        # post mapping QC
        print("scoring")
        if (isTRUE(x = do.bridge)){
          print("doing bridge cluster preservation score")
          app.env$object[['refAssay']] <- app.env$object[['ATAC']]
          DefaultAssay(app.env$object) <- 'refAssay'
          DefaultAssay(app.env$object[["ref.Bridge.reduc"]]) <- 'refAssay'
          app.env$object <- FindTopFeatures(app.env$object,
                                            min.cutoff = "q0")
          qc.stat <- round(
            x = ClusterPreservationScore(
              query = app.env$object,
              ds.amount = getOption(x = "Azimuth.map.postmapqcds"),
              type = "bridge"
            ),
            digits = 2
          )
        } else {
          print("doing standard scoring")
          qc.stat <- round(
            x = ClusterPreservationScore(
              query = app.env$object,
              ds.amount = getOption(x = "Azimuth.map.postmapqcds"),
              type = "standard"
            ),
            digits = 2
          )
        }
        if (!is.null(googlesheet)) {
          try(sheet_append(
            ss = googlesheet,
            data = data.frame(
              "CLUSTERPRESERVATIONQC",
              app_session_id,
              qc.stat
            )
          ))
        }
        app.env$clusterpreservationqc <- qc.stat
        if (qc.stat <  getOption(x = "Azimuth.map.postmapqccolors")[1]) {
          output$valuebox_mappingqcstat <- renderValueBox(expr = {
            valueBox(
              value = paste0(qc.stat, "/5"),
              subtitle = "cluster preservation score",
              color = 'red',
              icon = icon(name = 'times')
            )
          })
        } else if (qc.stat <  getOption(x = "Azimuth.map.postmapqccolors")[2]) {
          output$valuebox_mappingqcstat <- renderValueBox(expr = {
            valueBox(
              value = paste0(qc.stat, "/5"),
              subtitle = "cluster preservation score",
              color = 'yellow',
              icon = icon(name = 'exclamation-circle')
            )
          })
        } else {
          output$valuebox_mappingqcstat <- renderValueBox(expr = {
            valueBox(
              value = paste0(qc.stat, "/5"),
              subtitle = "cluster preservation score",
              color = 'green',
              icon = icon(name = 'check')
            )
          })
        }
        if (isTRUE(x = do.bridge)){
          print("doing bridge and gonna subset")
          refdr <- subset(
            x = app.env$anchors@object.list[[1]][["Bridge.reduc"]], # im gonna try calling this Bridge.Reduc
            cells = paste0(Cells(x = app.env$object), "_query")
          )
          refdr <- RenameCells(
            object = refdr, 
            new.names = Cells(x = app.env$object)
          )
          print("BRIDGE REFDR")
          print(refdr)
          refdr.ref <- subset(
            x = app.env$anchors@object.list[[1]][["Bridge.reduc"]], 
            cells = paste0(Cells(x = refs$map), "_reference")
          )
          refdr.ref <- RenameCells(
            object = refdr.ref, 
            new.names = Cells(x = refs$map)
          )
        } else {
          print("doing the standard version")
          refdr <- subset(
            x = app.env$anchors@object.list[[1]][["pcaproject.l2"]],
            cells = paste0(Cells(x = app.env$object), "_query")
          )
          refdr <- RenameCells(
            object = refdr,
            new.names = Cells(x = app.env$object)
          )
          refdr.ref <- subset(
            x = app.env$anchors@object.list[[1]][["pcaproject.l2"]],
            cells = paste0(Cells(x = refs$map), "_reference")
          )
          refdr.ref <- RenameCells(
            object = refdr.ref,
            new.names = Cells(x = refs$map)
          )
        }
        if (Sys.getenv("RSTUDIO") == "1") {
          plan("sequential")
        }
        # reduce size of object in anchorset
        print("doing diet seurat")
        app.env$anchors@object.list[[1]] <- DietSeurat(
          object = app.env$anchors@object.list[[1]]
        )
        print("subsetting")
        app.env$anchors@object.list[[1]] <- subset(
          x = app.env$anchors@object.list[[1]],
          features = c(rownames(x = app.env$anchors@object.list[[1]])[1])
        )
        app.env$anchors@object.list[[1]] <- RenameCells(
          object = app.env$anchors@object.list[[1]],
          new.names = unname(obj = sapply(
            X = Cells(x = app.env$anchors@object.list[[1]]),
            FUN = function(x) {
              return(gsub(pattern = "_reference", replacement = "", x = x))
            }
          )))
        app.env$anchors@object.list[[1]] <- RenameCells(
          object = app.env$anchors@object.list[[1]],
          new.names = sapply(
            X = Cells(x = app.env$anchors@object.list[[1]]),
            FUN = function(x) {
              return(gsub(pattern = "_query", replacement = "", x = x))
            },
            USE.NAMES = FALSE
          )
        )
        app.env$anchors@object.list[[1]]@meta.data <- data.frame()
        app.env$anchors@object.list[[1]]@active.ident <- factor()
        mapping.score.k <- min(c(
          50,
          length(x = unique(x = app.env$anchors@anchors[, 1])),
          length(x = unique(x = app.env$anchors@anchors[, 2])))
        )
        app.env$mapping.score <- future(
          expr = {
            MappingScore(
              anchors = app.env$anchors@anchors,
              combined.object = app.env$anchors@object.list[[1]],
              query.neighbors =  slot(object = app.env$anchors, name = "neighbors")[["query.neighbors"]],
              query.weights = Tool(object = app.env$object, slot = "TransferData")$weights.matrix,
              query.embeddings = Embeddings(object = refdr),
              ref.embeddings = Embeddings(object = refdr.ref),
              nn.method = "annoy",
              n.trees = n.trees,
              ndim = getOption(x = "Azimuth.map.ndims"),
              kanchors = mapping.score.k
            )
          }
        )
        app.env$object <- AddMetaData(
          object = app.env$object,
          metadata = rep(x = 0, times = ncol(x = app.env$object)),
          col.name = "mapping.score"
        )
        app.env$anchors <- NULL
        rm(refdr)
        gc(verbose = FALSE)
        # react.env$transform <- TRUE
        react.env$cluster.score <- TRUE
        react.env$score <- FALSE
      }
    }
  )
  observeEvent(
    eventExpr = react.env$transform,
    handlerExpr = {
      if (isTRUE(x = react.env$transform)) {
        if (isTRUE(x = do.bridge)) {
          react.env$progress$set(value = 0.8)
          suppressWarnings(expr = app.env$object[["umap.proj"]] <- app.env$object[["ref.umap"]])
        }
        else {
          react.env$progress$set(value = 0.8, message = 'Running UMAP transform')
          app.env$object[["query_ref.nn"]] <- FindNeighbors(
            object = Embeddings(refs$map[["refDR"]])[, 1:getOption("Azimuth.map.ndims")],
            query = Embeddings(app.env$object[["integrated_dr"]]),
            return.neighbor = TRUE,
            l2.norm = TRUE,
            n.trees = n.trees
          )
          app.env$object <- NNTransform(
            object = app.env$object,
            meta.data = refs$map[[]]
          )
          app.env$object[['umap.proj']] <- RunUMAP(
            object = app.env$object[['query_ref.nn']],
            reduction.model = refs$map[['refUMAP']],
            reduction.key = 'UMAP_'
          )
          app.env$object <- SetAssayData(
            object = app.env$object,
            assay = 'refAssay',
            slot = 'scale.data',
            new.data = new(Class = 'matrix')
          )
        }
        gc(verbose = FALSE)
        app.env$messages <- c(
          app.env$messages,
          paste(ncol(x = app.env$object), "cells mapped")
        )
        react.env$biomarkers <- TRUE
        if (isTRUE(x = do.bridge)) {
          react.env$motif <- TRUE
        }
        react.env$transform <- FALSE
      }
    }
  )
  observeEvent(
    eventExpr = react.env$biomarkers,
    handlerExpr = {
      if (isTRUE(x = react.env$biomarkers)) {
        react.env$progress$set(
          value = 0.95,
          message = 'Running differential expression'
        )
        app.env$gene.assay <- "RNA"
        for (i in app.env$metadataxfer[!app.env$singlepred]) {
          app.env$diff.expr[[paste(app.env$gene.assay, i, sep = "_")]] <- wilcoxauc(
            X = app.env$object,
            group_by = paste0("predicted.", i),
            assay = 'data',
            seurat_assay = app.env$gene.assay
          )
          if (isTRUE(x = do.adt)) {
            app.env$diff.expr[[paste(adt.key, i, sep = "_")]] <- wilcoxauc(
              X = app.env$object,
              group_by = paste0("predicted.", i),
              assay = 'data',
              seurat_assay = adt.key
            )
          }
        }
        
        # Finalize the log
        mapping.time <- difftime(
          time1 = Sys.time(),
          time2 = react.env$start,
          units = 'secs'
        )
        time.fmt <- FormatDiffTime(dt = mapping.time)
        app.env$messages <- c(
          app.env$messages,
          time.fmt
        )
        if (!is.null(x = googlesheet)) {
          try(expr = sheet_append(
            ss = googlesheet,
            data = data.frame(
              "MAPPINGTIME",
              app_session_id,
              as.numeric(x = mapping.time)
            )
          ))
        }
        if (!is.null(x = googlesheet)) {
          try(
            expr = sheet_append(
              ss = googlesheet,
              data = data.frame(
                "SUMMARY",
                app_session_id,
                basename(getOption(x = 'Azimuth.app.reference')),
                ReferenceVersion(object = refs$map),
                app.env$demo,
                app.env$ncellsupload,
                app.env$ncellspreproc,
                as.numeric(x = mapping.time),
                Sys.Date(),
                app.env$nanchors,
                app.env$clusterpreservationqc
              )
            ),
            silent = TRUE
          )
        }
        if (isTRUE(x = do.bridge)) {
          output$menu2 <- renderMenu(expr = {
            sidebarMenu(
              menuItem(
                text = "Cell Plots",
                tabName = "tab_cell",
                icon = icon("chart-area")
              ),
              menuItem(
                text = "Feature Plots",
                tabName = "tab_feature",
                icon = icon("chart-area")
              ),
              menuItem(
                text = "Motifs",
                tabName = "tab_motif",
                icon = icon("chart-area")
              ),
              menuItem(
                text = "Download Results",
                tabName = "tab_download",
                icon = icon("file-download")
              )
            )
          })
        } else {
          output$menu2 <- renderMenu(expr = {
            sidebarMenu(
              menuItem(
                text = "Cell Plots",
                tabName = "tab_cell",
                icon = icon("chart-area")
              ),
              menuItem(
                text = "Feature Plots",
                tabName = "tab_feature",
                icon = icon("chart-area")
              ),
              menuItem(
                text = "Download Results",
                tabName = "tab_download",
                icon = icon("file-download")
              )
            )
          })
        }
        app.env$object <- RenameCells(object = app.env$object, new.names = app.env$query.names)
        if (!isTRUE(x = do.bridge)) {
          react.env$progress$close()
        }
        enable(id = 'file')
        ToggleDemos(action = "enable", demos = demos)
        react.env$metadata <- TRUE
        react.env$biomarkers <- FALSE
      }
    }
  )
  observeEvent(eventExpr = react.env$motif, handlerExpr = {
    if (isTRUE(x = react.env$motif)) {
      react.env$progress$set(value = 0.98, message = "Running Motif Analysis")
      DefaultAssay(app.env$object) <- "ATAC"
      # Remove peaks on scaffolds 
      main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
      keep.peaks <- which(as.character(seqnames(granges(app.env$object))) %in% main.chroms)
      app.env$object[["ATAC"]] <- subset(app.env$object[["ATAC"]], features = rownames(app.env$object[["ATAC"]])[keep.peaks])
      
      pfm <- getMatrixSet(
        x = JASPAR2020,
        opts = list(species = 9606, all_versions = FALSE)
      )
      print("adding motifs")
      # startTime <- Sys.time()
      # app.env$object <- AddMotifs(
      #   object = app.env$object,
      #   genome = BSgenome.Hsapiens.UCSC.hg38,
      #   pfm = pfm
      # )
      # endTime <- Sys.time()
      # print("TOTAL TIME TO ADD MOTIFS")
      # print(endTime - startTime)
      # print("calculating motif")
      # library(BiocParallel)
      # register(MulticoreParam(3))
      # startTime <- Sys.time()
      # app.env$object <- Runmotif(
      #   object = app.env$object,
      #   genome = BSgenome.Hsapiens.UCSC.hg38
      # )
      # endTime <- Sys.time()
      # print("TOTAL TIME TO RUN motif")
      # print(endTime - startTime)
      # # Rename motifs from ids
      # motif_name <- ConvertMotifID(app.env$object[["peak.orig"]]@motifs, id = rownames(app.env$object[[app.env$default.assay]]@data))
      # rownames(app.env$object[[app.env$default.assay]]@data) <- motif_name
      
      # print(head(row.names(app.env$object[[app.env$default.assay]]@data)))
      # for (i in app.env$metadataxfer[!app.env$singlepred]) {
      #   print("setting motif.diff.expr")
      #   Idents(app.env$object) <- paste0("predicted.", i)
      #   app.env$default.assay <- "motif"
      #   app.env$motif.diff.expr[[paste(app.env$default.assay, # changed all of these to default.assay
      #                                     i, sep = "_")]] <- FindAllMarkers(object = app.env$object, assay = app.env$default.assay, slot = "data", 
      #                                                                       only.pos = T, mean.fcn = rowMeans, fc.name = "avg_diff")
      #   motif_ids <- ConvertMotifID(app.env$object[["peak.orig"]]@motifs, name = app.env$motif.diff.expr[[paste(app.env$default.assay, i, sep = "_")]]$gene)
      #   
      #   print("MOTIF IDS")
      #   print(head(motif_ids))
      #   app.env$motif.diff.expr[[paste(app.env$default.assay, i, sep = "_")]]$motif_id <- motif_ids
      #   print("column names of differential expression")
      #   print(colnames(app.env$motif.diff.expr[[paste(app.env$default.assay, i, sep = "_")]]))
      #   
      # }
      
      
      # FindMotif version 
      print("finding motifs")
      for (i in app.env$metadataxfer[!app.env$singlepred]) {
        app.env$peaks.diff.expr[[paste(app.env$default.assay, i, sep = "_")]] <- wilcoxauc(X = app.env$object,
                                                                                           group_by = paste0("predicted.", i),
                                                                                           assay = "data", 
                                                                                           seurat_assay = app.env$default.assay)
        peaks.list <- split(app.env$peaks.diff.expr[[paste(app.env$default.assay, i, sep = "_")]], 
                            f = app.env$peaks.diff.expr[[paste(app.env$default.assay, i, sep = "_")]]$group)
        motif.list <- list()
        for (num in 1:length(peaks.list)){
          print("starting to find motifs")
          if (nrow(peaks.list[[num]]) > 0){
            peaks.list[[num]] <- peaks.list[[num]][order(peaks.list[[num]]$logFC, decreasing = TRUE), ]
            if (nrow(peaks.list[[num]]) > 1000) {
              print("over 1000 peaks")
              top.da.peak <- peaks.list[[num]][1:1000,]$feature   #[peaks.list[[num]]$logFC > 0.5, ]$feature
            } else {
              print("smaller set of peaks")
              top.da.peak <- peaks.list[[num]][peaks.list[[num]]$pval < 0.05, ]$feature
            }
            print(head(top.da.peak))
            enriched.motifs <- FindMotifs( 
              object = refs$bridge[["ATAC"]],
              features = top.da.peak)
            enriched.motifs$group <- names(peaks.list[num])
            motif.list[[num]] <- enriched.motifs
            print(head(enriched.motifs))
          }  
        }
        print(head(dplyr::bind_rows(motif.list)))
        app.env$motif.diff.expr[[paste(app.env$default.assay, i, sep = "_")]] <- dplyr::bind_rows(motif.list)
        print("MOTIF DIFF EXPR ")
        print(head(app.env$motif.diff.expr))
        print(head(app.env$motif.diff.expr[[paste(app.env$default.assay, i, sep = "_")]]))
        
      }
      print("about to close progress")
      react.env$progress$close()
      react.env$motif <- FALSE
    }
  })
  # Update input controls
  observeEvent(
    eventExpr = react.env$metadata,
    handlerExpr = {
      if (isTRUE(x = react.env$metadata)) {
        print("at metadata")
        #  Add the discrete metadata dropdowns
        metadata.discrete <- sort(
          x = PlottableMetadataNames(
            object = app.env$object,
            exceptions = app.env$metadataxfer,
            min.levels = 1,
            max.levels = 50
          )
        )
        app.env$metadata.discrete <- metadata.discrete
        for (id in c('metarow', 'metacol', 'metagroup', 'metagroup.motif')) {
          if (id == 'metarow') {
            show.metadata <- 'query'
          } else {
            show.metadata <- paste0("predicted.", app.env$default.metadata)
          }
          updateSelectizeInput(
            session = session,
            inputId = id,
            choices = metadata.discrete,
            selected = show.metadata,
            server = TRUE,
            options = selectize.opts
          )
        }
        updateSelectizeInput(
          session = session,
          inputId = 'metacolor.query',
          choices = c(grep(pattern = '^predicted.', x = metadata.discrete, value = TRUE),
                      grep(pattern = '^predicted.', x = metadata.discrete, value = TRUE, invert = TRUE)),
          selected = paste0("predicted.", app.env$default.metadata),
          server = TRUE,
          options = selectize.opts[-which(x = names(x = selectize.opts) == 'maxItems')]
        )
        # Add the continuous metadata dropdown
        metadata.cont <- sort(x = setdiff(
          x = colnames(x = app.env$object[[]]),
          y = metadata.discrete
        ))
        metadata.cont <- Filter(
          f = function(x) {
            return(is.numeric(x = app.env$object[[x, drop = TRUE]]))
          },
          x = metadata.cont
        )
        # Add prediction scores for all classes to continuous metadata
        metadata.cont <- sort(x = metadata.cont)
        if (any(grepl(pattern = "predicted.*.score", x = metadata.cont))) {
          metadata.cont <- metadata.cont[-grep(pattern = "predicted.*.score", x = metadata.cont)]
        }
        if (any(grepl(pattern = "*_refAssay", x = metadata.cont))) {
          metadata.cont <- metadata.cont[-grep(pattern = "*_refAssay", x = metadata.cont)]
        }
        max.predictions <- paste0("predicted.", app.env$metadataxfer, ".score")
        names(x = max.predictions) <- app.env$metadataxfer
        max.predictions <- as.list(x = max.predictions)
        prediction.score.names <-
          lapply(X = app.env$metadataxfer, FUN = function(x) {
            key <- Key(object = app.env$object[[paste0("prediction.score.", x)]])
            ids <- paste0(rownames(x = app.env$object[[paste0("prediction.score.", x)]]))
            values <- paste0(key, ids)
            names(x = values) <- ids
            return(values)
          })
        names(x = prediction.score.names) <- paste0("Prediction scores - ", app.env$metadataxfer)
        metadata.cont <- c(
          list("Max prediction scores" = max.predictions),
          prediction.score.names,
          list("Other Metadata" = metadata.cont)
        )
        app.env$metadata.cont <- metadata.cont
        updateSelectizeInput(
          session = session,
          inputId = 'metadata.cont',
          choices = app.env$metadata.cont,
          selected = '',
          server = TRUE,
          options = selectize.opts
        )
        updateSelectizeInput(
          session = session,
          inputId = 'metadata.cont.motif',
          choices = app.env$metadata.cont,
          selected = '',
          server = TRUE,
          options = selectize.opts
        )
        updateSelectizeInput(
          session = session,
          inputId = 'metacolor.ref',
          choices = c(grep(pattern = '^predicted.', x = app.env$metadataxfer, value = TRUE), # re-ordering not working...
                      grep(pattern = '^predicted.', x = app.env$metadataxfer, value = TRUE, invert = TRUE)),
          selected = app.env$default.metadata,
          server = TRUE,
          options = selectize.opts[-which(x = names(x = selectize.opts) == 'maxItems')]
        )
        react.env$features <- TRUE
        if (isTRUE(x = do.bridge)) {
          react.env$motif.features <- TRUE
        } 
        react.env$metadata <- FALSE
      }
    }
  )
  observeEvent(
    eventExpr = react.env$features,
    handlerExpr = {
      if (isTRUE(x = react.env$features)) {
        print("doing features")
        DefaultAssay(app.env$object) <- "RNA"
        app.env$default.feature <- ifelse(
          test = getOption(x = 'Azimuth.app.default_gene') %in% rownames(x = app.env$object),
          yes = getOption(x = 'Azimuth.app.default_gene'),
          no = VariableFeatures(object = app.env$object)[1]
        )
        print("DEFAULT FEATURE")
        print(app.env$default.feature)
        app.env$features <- unique(x = c(
          FilterFeatures(
            features = VariableFeatures(object = app.env$object)[1:selectize.opts$maxOptions]
          ),
          FilterFeatures(features = rownames(x = app.env$object))
        ))
        print(head(app.env$features))
        updateSelectizeInput(
          session = session,
          inputId = 'feature',
          label = 'Feature',
          choices = app.env$features,
          selected = app.env$default.feature,
          server = TRUE,
          options = selectize.opts
        )
        print('should have made input feature')
        print("featire ")
        if (isTRUE(x = do.adt)) {
          # app.env$adt.features <- sort(x = FilterFeatures(features = rownames(
          #   x = app.env$object[[adt.key]]
          # )))
          app.env$adt.features <- sort(x = rownames(
            x = app.env$object[[adt.key]]
          ))
          updateSelectizeInput(
            session = session,
            inputId = 'adtfeature',
            choices = app.env$adt.features,
            selected = '',
            server = TRUE,
            options = selectize.opts
          )
        }
        react.env$features <- FALSE
        if (!isTRUE(x = do.bridge)){
          react.env$markers <- TRUE
        }
      }
    }
  )
  observeEvent(
    eventExpr = react.env$motif.features, 
    handlerExpr = {
      if (isTRUE(x = react.env$motif.features)) {
        print("doing motif features")
        print("printing")
        DefaultAssay(app.env$object) <- app.env$default.assay
        print("APP ENV motif FEATURE DEFAULT")
        print(head(row.names(app.env$object[[app.env$default.assay]]@data)))
        app.env$default.motif.feature <- ifelse(test = getOption(x = 'Azimuth.app.default_motif') %in% 
                                                  row.names(x = app.env$object[[app.env$default.assay]]@data), 
                                                yes = getOption(x = 'Azimuth.app.default_motif'), 
                                                no = row.names(x = app.env$object[[app.env$default.assay]]@data)[1])
        print(app.env$default.motif.feature)
        app.env$motif.features <- unique(x = row.names(x = app.env$object[[app.env$default.assay]]@data)) # c(FilterFeatures(features =
        print(head(app.env$motif.features))
        print(app.env$default.motif.feature %in% app.env$motif.features)
        updateSelectizeInput(session = session, inputId = "motif.feature", 
                             label = "Motif", choices = app.env$motif.features, 
                             selected = app.env$default.motif.feature, server = TRUE, 
                             options = selectize.opts)
        print("should have made motif feature")
        
        if (isTRUE(x = do.adt)) {
          app.env$adt.features <- sort(x = rownames(x = app.env$object[[adt.key]]))
          updateSelectizeInput(session = session, inputId = "adtfeature", 
                               choices = app.env$adt.features, selected = "", 
                               server = TRUE, options = selectize.opts)
        }
        react.env$motif.features <- FALSE
        react.env$markers <- TRUE
        print("finished motif features")
      }
    }
  )
  observeEvent(
    eventExpr = react.env$markers,
    handlerExpr = {
      if (isTRUE(x = react.env$markers)) {
        print("doing markers")
        allowed.clusters <- names(x = which(
          x = table(app.env$object[[paste0("predicted.", app.env$default.metadata)]]) > getOption(x = 'Azimuth.de.mincells')
        ))
        allowed.clusters <- factor(
          x = allowed.clusters,
          levels = unique(x = app.env$object[[paste0("predicted.", app.env$default.metadata), drop = TRUE]])
        )
        allowed.clusters <- sort(x = levels(x = droplevels(x = na.omit(
          object = allowed.clusters
        ))))
        print(allowed.clusters)
        # updateSelectizeInput(
        #   session = session,
        #   inputId = 'select.prediction',
        #   choices = allowed.clusters,
        #   selected = allowed.clusters[1],
        #   server = TRUE,
        #   options = selectize.opts
        # )
        
        updateSelectizeInput(
          session = session,
          inputId = 'markerclusters',
          choices = allowed.clusters,
          selected = allowed.clusters[1],
          server = TRUE,
          options = selectize.opts
        )
        
        updateSelectizeInput(
          session = session,
          inputId = 'markerclustersgroup',
          choices = app.env$metadataxfer[!app.env$singlepred],
          selected = app.env$default.metadata,
          server = TRUE,
          options = selectize.opts
        )
        
        updateSelectizeInput(
          session = session,
          inputId = 'markerclustersgroup.motif',
          choices = app.env$metadataxfer[!app.env$singlepred],
          selected = app.env$default.metadata,
          server = TRUE,
          options = selectize.opts
        )
        
        react.env$markers <- FALSE
        app.env$disable <- FALSE
        react.env$get.feature <- TRUE
        if (isTRUE(x = do.bridge)){
          react.env$get.motif.feature <- TRUE
        }
      }
    }
  )
  observeEvent(
    eventExpr = react.env$no,
    handlerExpr = {
      if (FALSE) {
        # Enable the feature explorer
        
        # Add the predicted ID and score to the plots
        
        
        
        # Enable downloads
        
        
        react.env$no <- FALSE
      }
    }
  )
  # Handle input changes
  observeEvent( # RNA feature
    eventExpr = input$feature,
    handlerExpr = {
      if (nchar(x = input$feature)) {
        if (nchar(x = input$markerclustersgroup)) {
          print(paste0("FEATURE IN NEW BLOCK: ", input$feature))
          if (isTRUE(x = do.bridge)) {
            app.env$feature <- ifelse(
              test = input$feature %in% rownames(x = app.env$object[[app.env$gene.assay]]),
              yes = paste0(
                Key(object = app.env$object[[app.env$gene.assay]]),
                input$feature
              ),
              no = input$feature
            )
          } else {
            app.env$feature <- ifelse(
              test = input$feature %in% rownames(x = app.env$object[["refAssay"]]),
              yes = paste0(
                Key(object = app.env$object[["refAssay"]]),
                input$feature
              ),
              no = input$feature
            )
          }
          for (f in c('adtfeature', 'metadata.cont')) {
            updateSelectizeInput(
              session = session,
              inputId = f,
              choices = list(
                'adtfeature' = app.env$adt.features,
                'metadata.cont' = app.env$metadata.cont
              )[[f]],
              selected = '',
              server = TRUE,
              options = selectize.opts
            )
          }
          table.check <- input$feature %in% rownames(x = RenderDiffExp(
            diff.exp = app.env$diff.expr[[paste(app.env$gene.assay, input$markerclustersgroup, sep = "_")]],
            groups.use = input$markerclusters,
            n = Inf
          ))
          print("DID TABLE CHECK FOR GET FEATURE")
          tables.clear <- list(adt.proxy, rna.proxy)[c(TRUE, !table.check)]
          for (tab in tables.clear) {
            selectRows(proxy = tab, selected = NULL)
          }
        }
      }
    }
  )
  
  observeEvent( # motif feature
    eventExpr = input$motif.feature,
    handlerExpr = {
      if (nchar(x = input$motif.feature)) {
        if (nchar(x = input$markerclustersgroup.motif)) {
          print(paste0("motif FEATURE IN NEW BLOCK: ", input$motif.feature))
          app.env$motif.feature <- ifelse(
            test = input$motif.feature %in% rownames(x = app.env$object),
            yes = paste0(
              Key(object = app.env$object[[app.env$default.assay]]),
              input$motif.feature
            ),
            no = input$motif.feature
          )
          table.check <- input$motif.feature %in% rownames(x = RenderDiffMotifExp(
            diff.exp = app.env$motif.diff.expr[[paste(app.env$default.assay, input$markerclustersgroup.motif, sep = "_")]],
            groups.use = input$markerclusters.motif,
            n = Inf
          ))
          print("DID TABLE CHECK FOR GET MOTIF FEATURE")
          print(table.check)
          tables.clear <- list(adt.proxy, motif.proxy)[c(TRUE, !table.check)]
          for (tab in tables.clear) {
            selectRows(proxy = tab, selected = NULL)
          }
        }
      }
    }
  )
  
  observeEvent( # Protein feature
    eventExpr = input$adtfeature,
    handlerExpr = {
      if (nchar(x = input$adtfeature)) {
        app.env$feature <- paste0(
          Key(object = app.env$object[[adt.key]]),
          input$adtfeature
        )
        for (f in c('feature', 'metadata.cont')) {
          updateSelectizeInput(
            session = session,
            inputId = f,
            choices = list(
              'feature' = app.env$features,
              'metadata.cont' = app.env$metadata.cont
            )[[f]],
            selected = '',
            server = TRUE,
            options = selectize.opts
          )
        }
        table.check <- input$adtfeature %in% rownames(x = RenderDiffExp(
          diff.exp = app.env$diff.expr[[paste(adt.key, input$markerclustersgroup, sep = "_")]],
          groups.use = input$markerclusters,
          n = Inf
        ))
        tables.clear <- list(rna.proxy, adt.proxy)[c(TRUE, !table.check)]
        for (tab in tables.clear) {
          selectRows(proxy = tab, selected = NULL)
        }
      }
    }
  )
  # observeEvent(
  #   eventExpr = react.env$get.feature, 
  #   handlerExpr = {
  #     if (isTRUE(react.env$get.feature)) {
  #       print("react.env$get.feature")
  #       print(app.env$feature)
  #       #req(input$markerclustersgroup)
  #       #req(input$feature)
  #       print(paste0("input feature:", input$feature))
  #       print("MARKER CLUSTERS GROUP")
  #       print(input$markerclustersgroup)
  #       table.check <- input$feature %in% rownames(x = RenderDiffExp(
  #         diff.exp = app.env$diff.expr[[paste(app.env$gene.assay, input$markerclustersgroup, sep = "_")]],
  #         groups.use = input$markerclusters,
  #         n = Inf
  #       ))
  #       print("DID TABLE CHECK FOR GET FEATURE")
  #       tables.clear <- list(adt.proxy, rna.proxy)[c(TRUE, !table.check)]
  #       for (tab in tables.clear) {
  #         selectRows(proxy = tab, selected = NULL)
  #       }
  #     }
  #     react.env$get.feature <- FALSE
  #   }
  # )
  # observeEvent(
  #   eventExpr = react.env$get.motif.feature, 
  #   handlerExpr = {
  #     if (isTRUE(react.env$get.motif.feature)) {
  #       # req(input$markerclusters)
  #       # req(input$markerclustersgroup.motif)
  #       # req(input$motif.feature)
  #       print(paste0("input motif feature:", input$motif.feature))
  #       print("MARKER CLUSTERS GROUP")
  #       print(input$markerclustersgroup.motif)
  #       table.check <- input$motif.feature %in% rownames(x = RenderDiffMotifExp(
  #         diff.exp = app.env$motif.diff.expr[[paste(app.env$default.assay, input$markerclustersgroup.motif, sep = "_")]],
  #         groups.use = input$markerclusters,
  #         n = Inf
  #       ))
  #       print(table.check)
  #       tables.clear <- list(adt.proxy, motif.proxy)[c(TRUE, !table.check)]
  #       for (tab in tables.clear) {
  #         selectRows(proxy = tab, selected = NULL)
  #       }
  #     }
  #   react.env$get.motif.feature <- FALSE
  #   }
  # )
  observeEvent( # Continuous Metadata
    eventExpr = input$metadata.cont,
    handlerExpr = {
      if (nchar(x = input$metadata.cont)) {
        if (input$metadata.cont == "mapping.score") {
          if (resolved(x = app.env$mapping.score)) {
            app.env$object$mapping.score <- value(app.env$mapping.score)
          }
        }
        app.env$feature <- input$metadata.cont
        for (f in c('feature', 'adtfeature')) {
          updateSelectizeInput(
            session = session,
            inputId = f,
            choices = list(
              'feature' = app.env$features,
              'adtfeature' = app.env$adt.features
            )[[f]],
            selected = '',
            server = TRUE,
            options = selectize.opts
          )
        }
        for (tab in list(rna.proxy, adt.proxy)) {
          selectRows(proxy = tab, selected = NULL)
        }
      }
    }
  )
  observeEvent( # Continuous Metadata
    eventExpr = input$metadata.cont.motif,
    handlerExpr = {
      if (nchar(x = input$metadata.cont.motif)) {
        if (input$metadata.cont.motif == "mapping.score") {
          if (resolved(x = app.env$mapping.score)) {
            app.env$object$mapping.score <- value(app.env$mapping.score)
          }
        }
        app.env$feature <- input$metadata.cont.motif
        updateSelectizeInput(
          session = session,
          inputId = "motif.feature",
          choices = app.env$motiffeatures,
          selected = '',
          server = TRUE,
          options = selectize.opts
        )
        for (tab in list(rna.proxy, adt.proxy)) {
          selectRows(proxy = tab, selected = NULL)
        }
      }
    }
  )
  observeEvent( # Marker clusters group
    eventExpr = input$markerclustersgroup,
    handlerExpr = {
      if (nchar(x = input$markerclustersgroup)) {
        allowed.clusters <- names(x = which(
          x = table(app.env$object[[paste0("predicted.", input$markerclustersgroup)]]) > getOption(x = 'Azimuth.de.mincells')
        ))
        allowed.clusters <- factor(
          x = allowed.clusters,
          levels = unique(x = app.env$object[[paste0("predicted.", input$markerclustersgroup), drop = TRUE]])
        )
        allowed.clusters <- sort(x = levels(x = droplevels(x = na.omit(
          object = allowed.clusters
        ))))
        app.env$allowedclusters <- allowed.clusters
        updateSelectizeInput(
          session = session,
          inputId = "markerclusters",
          choices = app.env$allowedclusters,
          selected = app.env$allowedclusters[1],
          server = TRUE,
          options = selectize.opts
        )
      }
    }
  )
  observeEvent( # Marker clusters group motif
    eventExpr = input$markerclustersgroup.motif,
    handlerExpr = {
      if (nchar(x = input$markerclustersgroup.motif)) {
        allowed.clusters <- names(x = which(
          x = table(app.env$object[[paste0("predicted.", input$markerclustersgroup.motif)]]) > getOption(x = 'Azimuth.de.mincells')
        ))
        allowed.clusters <- factor(
          x = allowed.clusters,
          levels = unique(x = app.env$object[[paste0("predicted.", input$markerclustersgroup.motif), drop = TRUE]])
        )
        allowed.clusters <- sort(x = levels(x = droplevels(x = na.omit(
          object = allowed.clusters
        ))))
        app.env$allowedclusters <- allowed.clusters
        updateSelectizeInput(
          session = session,
          inputId = "markerclusters.motif",
          choices = app.env$allowedclusters,
          selected = app.env$allowedclusters[1],
          server = TRUE,
          options = selectize.opts
        )
      }
    }
  )
  observeEvent( # Select from biomarkers table
    eventExpr = input$biomarkers_rows_selected,
    handlerExpr = {
      if (length(x = input$biomarkers_rows_selected)) {
        updateSelectizeInput(
          session = session,
          inputId = 'feature',
          choices = app.env$features,
          selected = rownames(x = RenderDiffExp(
            diff.exp = app.env$diff.expr[[paste(app.env$default.assay, input$markerclustersgroup, sep = "_")]],
            groups.use = input$markerclusters,
            n = Inf
          ))[input$biomarkers_rows_selected],
          server = TRUE,
          options = selectize.opts
        )
      }
    }
  )
  observeEvent( # Select from adtbio table
    eventExpr = input$adtbio_rows_selected,
    handlerExpr = {
      if (length(x = input$adtbio_rows_selected)) {
        updateSelectizeInput(
          session = session,
          inputId = 'adtfeature',
          choices = app.env$adt.features,
          selected = rownames(x = RenderDiffExp(
            diff.exp = app.env$diff.expr[[paste(adt.key, input$markerclustersgroup, sep = "_")]],
            groups.use = input$markerclusters,
            n = Inf
          ))[input$adtbio_rows_selected],
          server = TRUE,
          options = selectize.opts
        )
      }
    }
  )
  observeEvent( # Select from motif table
    eventExpr = input$motif_rows_selected,
    handlerExpr = {
      if (length(x = input$motif_rows_selected)) {
        print("selecting motif rows ")
        updateSelectizeInput(
          session = session,
          inputId = 'motif',
          choices = app.env$motif.features,
          selected = rownames(x = RenderDiffMotifExp(
            diff.exp = app.env$motif.diff.expr[[paste(app.env$default.assay, input$markerclustersgroup.motif, sep = "_")]],
            groups.use = input$markerclusters.motif,
            n = Inf
          ))[input$motif_rows_selected],
          server = TRUE,
          options = selectize.opts
        )
      }
    }
  )
  observeEvent( # Once reference metadata is initialized, set the default legend/labels based on ref only
    eventExpr = input$metacolor.ref,
    handlerExpr = {
      if (length(x = unique(x = as.vector(x = refs$plot[[input$metacolor.ref[1], drop = TRUE]]))) >= 30) {
        if (isTRUE(app.env$fresh.plot)) {
          updateCheckboxInput(
            session = session,
            inputId = 'labels',
            value = TRUE
          )
          updateCheckboxInput(
            session = session,
            inputId = 'legend',
            value = FALSE
          )
          app.env$fresh.plot <- FALSE
        }
      }
    }
  )
  observeEvent( # Record feedback and update UI if feedback submitted
    eventExpr = input$submit_feedback,
    handlerExpr = {
      if (!is.null(x = googlesheet)) {
        try(expr = sheet_append(
          ss = googlesheet,
          data = data.frame(
            "FEEDBACK",
            app_session_id,
            paste0('feedback: \"', input$feedback, '\"')
          )
        ))
      }
      updateTextAreaInput(
        session = session,
        inputId = 'feedback',
        label = NULL,
        value = 'Thank you for your feedback!')
    }
  )
  observeEvent( # Change metadata appropriately
    eventExpr = input$showrefonly,
    handlerExpr = {
      if (!is.null(app.env$metadata.discrete)) {
        if (input$showrefonly) {
          # change to appropriate input$metacolor.ref if its an option
          disable(id = 'metacolor.query')
          enable(id = 'metacolor.ref')
        } else {
          # change to appropriate input$metacolor.query
          disable(id = 'metacolor.ref')
          enable(id = 'metacolor.query')
        }
      }
    }
  )
  # Plots
  output$plot.qc <- renderPlot(expr = {
    if (!is.null(x = isolate(expr = app.env$object)) & isTRUE(x = react.env$plot.qc)) {
      # all(paste0(c('nCount_', 'nFeature_'), app.env$default.assay) %in% colnames(app.env$object@meta.data)) ) {
      qc <- paste0(c('nCount_', 'nFeature_'), app.env$default.assay)
      if (isTRUE(x = react.env$mt)) {
        qc <- c(qc, mt.key)
      }
      check.qcpoints <- "qcpoints" %in% input$check.qc
      check.qcscale <- "qcscale" %in% input$check.qc
      vlnlist <- VlnPlot(
        object = isolate(app.env$object),
        features = qc,
        group.by = 'query',
        combine = FALSE,
        pt.size = ifelse(
          test = check.qcpoints,
          yes = 0,
          no = Seurat:::AutoPointSize(data = isolate(app.env$object))
        ),
        log = check.qcscale
      )
      # nCount
      vlnlist[[1]] <- vlnlist[[1]] +
        geom_hline(yintercept = input$num.ncountmin) +
        geom_hline(yintercept = input$num.ncountmax) +
        annotate(
          geom = "rect",
          alpha = 0.2,
          fill = "red",
          ymin = input$num.ncountmax,
          ymax = Inf,
          xmin = 0.5,
          xmax = 1.5
        ) +
        annotate(
          geom = "rect",
          alpha = 0.2,
          fill = "red",
          ymin = ifelse(test = check.qcscale, yes = 0, no = -Inf),
          ymax = input$num.ncountmin,
          xmin = 0.5,
          xmax = 1.5
        ) +
        NoLegend() +
        xlab("")
      # nFeature
      vlnlist[[2]] <- vlnlist[[2]] +
        geom_hline(yintercept = input$num.nfeaturemin) +
        geom_hline(yintercept = input$num.nfeaturemax) +
        annotate(
          geom = "rect",
          alpha = 0.2,
          fill = "red",
          ymin = input$num.nfeaturemax,
          ymax = Inf,
          xmin = 0.5,
          xmax = 1.5
        ) +
        annotate(
          geom = "rect",
          alpha = 0.2,
          fill = "red",
          ymin = ifelse(test = check.qcscale, yes = 0, no = -Inf),
          ymax = input$num.nfeaturemin,
          xmin = 0.5,
          xmax = 1.5
        ) +
        NoLegend() +
        xlab("")
      if (react.env$mt) {
        vlnlist[[3]] <- vlnlist[[3]] +
          geom_hline(yintercept = input$minmt) +
          geom_hline(yintercept = input$maxmt) +
          annotate(
            geom = "rect",
            alpha = 0.2,
            fill = "red",
            ymin = input$maxmt,
            ymax = Inf,
            xmin = 0.5,
            xmax = 1.5
          ) +
          annotate(
            geom = "rect",
            alpha = 0.2,
            fill = "red",
            ymin = ifelse(test = check.qcscale, yes = 0, no = -Inf),
            ymax = input$minmt,
            xmin = 0.5,
            xmax = 1.5
          ) +
          NoLegend() +
          xlab("")
      }
      wrap_plots(vlnlist, ncol = length(x = vlnlist))
    }
  })
  output$overlap_box <- renderUI(
    box(
      title = p(
        'Overlap QC',
        bsButton(
          inputId = 'q4',
          label = '',
          icon = icon(name = 'question'),
          style = 'info',
          size = 'extra-small'
        )
      ),
      bsPopover(
        id = 'q4',
        title = 'Overlap QC',
        content = paste(
          'The distribution of overlap percentages for each peak. A strongly left-skewed ',
          'distribution means that most of the peaks have ~100% overlap to the corresponding multiome peak', 
          'and thus the requantified peaks will (maintain) the data from the original peaks. Also, note the ', 
          'total overlap percentage for a summary of this information.'
        ),
        placement = 'right',
        trigger = 'focus',
        options = list(container = 'body')
      ),
      plotOutput(outputId = 'dist.qc'),
      width = 4
    )
  )
  output$dist.qc <- renderPlot(expr = {
    if (!is.null(x = isolate(expr = app.env$chromatin_assay_1)) & isTRUE(x = react.env$dist.qc)) {
      print("making dist plots")
      dist <- OverlapDistPlot(query_assay = isolate(app.env$chromatin_assay_1),
                              multiome = refs$bridge[["ATAC"]])
    }
  })
  output$refdim_intro <- renderPlot(expr = {
    # save plot dataframe to minimize on-hover computation
    app.env$plots.refdim_intro_df <- cbind(
      as.data.frame(x = Embeddings(object = refs$plot[['refUMAP']])),
      refs$plot[[]]
    )
    p <- DimPlot(
      object = refs$plot,
      combine = FALSE,
      group.by = default_xfer,
      cols = GetColorMap(object = refs$map)[[default_xfer]],
      repel = TRUE,
      label = TRUE,
      raster = FALSE
    )[[1]]
    # for later use by query plot:
    app.env$plot.ranges <- list(
      layer_scales(p)$x$range$range,
      layer_scales(p)$y$range$range
    )
    # strip down the intro plot-- no title, legend, or axes
    p + WelcomePlot()
  })
  output$refdim_intro_hover_box <- renderUI({
    hover <- input$refdim_intro_hover_location
    df <- app.env$plots.refdim_intro_df
    if (!is.null(x = hover)){
      hover[['mapping']] <- setNames(object = as.list(x = colnames(x = app.env$plots.refdim_intro_df)[1:2]), nm = c('x', 'y'))
    }
    point <- nearPoints(
      df = df,
      coordinfo = hover,
      threshold = 10,
      maxpoints = 1,
      addDist = TRUE
    )
    if (nrow(x = point) == 0) {
      return(NULL)
    }
    hovertext <- do.call(
      what = paste0,
      args = lapply(X = metadata.annotate, FUN = function(md) {
        paste0("<span>", md, "</span>: <i>", point[[md]], "</i><br>")
      })
    )
    wellPanel(
      style = HoverBoxStyle(x = hover$coords_css$x, y = hover$coords_css$y),
      p(HTML(text = hovertext))
    )
  })
  
  output$refdim <- renderPlot(expr = {
    if (!is.null(x = input$metacolor.ref)) {
      colormaps <- GetColorMap(object = refs$map)[input$metacolor.ref]
      # no interactivity if multiple plots per row (less useful in this case)
      if (length(x = colormaps) == 1) {
        ## already stored reference dataframe in app.env
        app.env$plots.refdim_df <- app.env$plots.refdim_intro_df
        DimPlot(
          object = refs$plot,
          label = isTRUE("labels" %in% input$dimplot.opts),
          group.by = input$metacolor.ref,
          cols = colormaps[[1]],
          repel = TRUE,
          raster = FALSE
        )[[1]] +
          labs(x = "UMAP 1", y = "UMAP 2") +
          if (isFALSE(x = "legend" %in% input$dimplot.opts) | OversizedLegend(refs$plot[[input$metacolor.ref, drop = TRUE]])) NoLegend()
      } else {
        app.env$plots.refdim_df <- NULL
        plots <- list()
        for (i in 1:length(x = colormaps)) {
          plots[[i]] <- DimPlot(
            object = refs$plot,
            label = isTRUE("labels" %in% input$dimplot.opts),
            group.by = input$metacolor.ref[i],
            cols = colormaps[[i]],
            repel = TRUE,
            raster = FALSE
          ) + labs(x = "UMAP 1", y = "UMAP 2") +
            if (isFALSE(x = "legend" %in% input$dimplot.opts) | OversizedLegend(refs$plot[[input$metacolor.ref[i], drop = TRUE]])) NoLegend()
        }
        wrap_plots(plots, nrow = 1)
      }
    }
  })
  
  output$refdim_hover_box <- renderUI({
    if (!is.null(x = app.env$plots.refdim_df)) {
      hover <- input$refdim_hover_location
      df <- app.env$plots.refdim_df
      if (!is.null(x = hover)){
        hover[['mapping']] <- setNames(object = as.list(x = colnames(x = app.env$plots.refdim_intro_df)[1:2]), nm = c('x', 'y'))
      }
      point <- nearPoints(
        df = df,
        coordinfo = hover,
        threshold = 10,
        maxpoints = 1,
        addDist = TRUE
      )
      if (nrow(x = point) == 0) {
        return(NULL)
      }
      hovertext <- do.call(
        what = paste0,
        args = as.list(c(
          paste0("<b>", point[[input$metacolor.ref]], "</b><br>"),
          sapply(X = setdiff(possible.metadata.transfer, input$metacolor.ref), FUN = function(md) {
            paste0("<span>", md, "</span>: <i>", point[[md]], "</i><br>")
          })
        ))
      )
      wellPanel(
        style = HoverBoxStyle(x = hover$coords_css$x, y = hover$coords_css$y),
        p(HTML(text = hovertext))
      )
    }
  })
  
  output$objdim <- renderPlot(expr = {
    if (!is.null(x = app.env$object) && app.env$disable == FALSE) {
      # create empty ref
      if (is.null(x = app.env$emptyref) | is.null(x = app.env$merged)) {
        app.env$emptyref <- refs$plot
        Idents(object = app.env$emptyref) <- '.'
        for (md in colnames(x = app.env$emptyref[[]])) {
          app.env$emptyref[[md]] <- '.'
        }
        for (md in setdiff(
          x = colnames(x = app.env$object[[]]),
          y = colnames(x = app.env$emptyref[[]])
        )) {
          app.env$emptyref[[md]] <- '.'
        }
        app.env$object[['refUMAP']] <- app.env$object[['umap.proj']]
        app.env$merged <- merge(app.env$emptyref, app.env$object, merge.dr = 'refUMAP')
      }
      
      if (isFALSE(x = input$showrefonly) &
          length(x = Reductions(object = app.env$object)) &
          !is.null(x = input$metacolor.query)) { # SHOW OVERLAY
        app.env$plots.refdim_df <- NULL # hide reference hover box
        if (length(x = input$metacolor.query) == 1) {
          # get colormap if avail
          group.var <- gsub(pattern = "^predicted.", replacement = "", x = input$metacolor.query)
          colormap <- GetColorMap(object = refs$map)[[group.var]]
          if (!grepl(pattern = "^predicted.", x = input$metacolor.query)) {
            colormap <- CreateColorMap(ids=unique(as.vector(app.env$object[[input$metacolor.query,drop=T]])))
          }
          colormap['.'] <- '#F1F1F1'
          
          # make dataframe so don't need to recompute during hover- QUERY only!
          app.env$plots.objdim_df <- cbind(
            as.data.frame(x = Embeddings(object = app.env$object[['umap.proj']])),
            app.env$object[[]]
          )
          p <- DimPlot(
            object = app.env$merged,
            group.by = input$metacolor.query,
            label = FALSE,
            cols = colormap[names(x = colormap) %in% c(
              '.', unique(x = as.vector(x = app.env$object[[input$metacolor.query, drop = TRUE]])))],
            repel = TRUE,
            raster = FALSE,
            reduction = "refUMAP"
          )[[1]] +
            xlim(app.env$plot.ranges[[1]]) +
            ylim(app.env$plot.ranges[[2]]) +
            labs(x = "UMAP 1", y = "UMAP 2") +
            if (isFALSE(x = input$legend)) NoLegend()
          if (isTRUE(x = 'labels' %in% input$label.opts)) {
            keep <- if (isTRUE(x = 'filterlabels' %in% input$label.opts)) {
              t <- table(as.vector(x = app.env$object[[input$metacolor.query, drop = TRUE]]))
              names(x = t)[which(x = t > 0.02 * ncol(x = app.env$object))]
            } else NULL
            return(LabelClusters(
              plot = p,
              id = input$metacolor.query,
              clusters = keep
            ))
          }
          return(p)
        } else {
          app.env$plots.objdim_df <- NULL # no interactivity
          plots <- list()
          for (i in 1:length(x = input$metacolor.query)) {
            group.var <- gsub(pattern = "^predicted.", replacement = "", x = input$metacolor.query[i])
            colormap <- GetColorMap(object = refs$map)[[group.var]]
            if (!grepl(pattern = "^predicted.", x = input$metacolor.query[i])) {
              colormap <- CreateColorMap(
                ids = unique(x = as.vector(x = app.env$object[[input$metacolor.query[i], drop = TRUE]]))
              )
            }
            colormap['.'] <- '#F1F1F1'
            p <- DimPlot(
              object = app.env$merged,
              group.by = input$metacolor.query[i],
              cols = colormap[names(x = colormap) %in% c(
                '.', unique(x = as.vector(x = app.env$object[[input$metacolor.query[i], drop = TRUE]])))],
              repel = TRUE,
              raster = FALSE,
              reduction = "refUMAP"
            )[[1]] + xlim(app.env$plot.ranges[[1]]) +
              ylim(app.env$plot.ranges[[2]]) +
              labs(x = "UMAP 1", y = "UMAP 2") +
              if (isFALSE(x = input$legend) | OversizedLegend(annotation.list = app.env$object[[input$metacolor.query[i], drop = TRUE]])) NoLegend()
            if (isTRUE('labels' %in% input$label.opts)) {
              keep <- if (isTRUE(x = 'filterlabels' %in% input$label.opts)) {
                t <- table(as.vector(x = app.env$object[[input$metacolor.query[i], drop = TRUE]]))
                print(names(x = t)[which(t > 0.02 * ncol(x = app.env$object))])
                names(x = t)[which(t > 0.02 * ncol(x = app.env$object))]
              } else NULL
              plots[[i]] <- LabelClusters(
                plot = p,
                id = input$metacolor.query[i],
                clusters = keep
              )
            } else {
              plots[[i]]<-p
            }
          }
          wrap_plots(plots, nrow = 1)
        }
      } else { # SHOW REFERENCE ONLY
        app.env$plots.objdim_df <- NULL # hide query hover box
        if (!is.null(x = input$metacolor.ref)) {
          colormaps <- GetColorMap(object = refs$map)[input$metacolor.ref]
          if (length(x = colormaps) == 1) {
            app.env$plots.refdim_df <- app.env$plots.refdim_intro_df
            DimPlot(
              object = refs$plot,
              label = isTRUE('labels' %in% input$label.opts),
              group.by = input$metacolor.ref,
              cols = colormaps[[1]],
              repel = TRUE,
              raster = FALSE
            )[[1]] +
              labs(x = "UMAP 1", y = "UMAP 2") +
              if (isFALSE(input$legend) | OversizedLegend(annotation.list = refs$plot[[input$metacolor.ref, drop = TRUE]])) NoLegend()
          } else {
            app.env$plots.refdim_df <- NULL
            plots <- list()
            for (i in 1:length(x = colormaps)) {
              plots[[i]] <- DimPlot(
                object = refs$plot,
                label = isTRUE('labels' %in% input$label.opts),
                group.by = input$metacolor.ref[i],
                cols = colormaps[[i]],
                repel = TRUE,
                raster = FALSE
              ) + labs(x = "UMAP 1", y = "UMAP 2") +
                if (isFALSE(x = input$legend) | OversizedLegend(annotation.list = refs$plot[[input$metacolor.ref[i], drop = TRUE]])) NoLegend()
            }
            wrap_plots(plots, nrow = 1)
          }
        }
      }
    }
  })
  output$querydim <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      if (length(x = Reductions(object = app.env$object)) & !is.null(x = input$metacolor.query)) {
        if (length(x = input$metacolor.query) == 1) {
          # get colormap if avail
          group.var <- gsub(pattern = "^predicted.", replacement = "", x = input$metacolor.query)
          colormap <- GetColorMap(object = refs$map)[[group.var]]
          if (!grepl(pattern = "^predicted.", x = input$metacolor.query)) {
            colormap <- NULL
          }
          # make dataframe so don't need to recompute during hover
          app.env$plots.querydim_df <- cbind(
            as.data.frame(x = Embeddings(object = app.env$object[['umap.proj']])),
            app.env$object[[]]
          )
          DimPlot(
            object = app.env$object,
            group.by = input$metacolor.query,
            label = isTRUE('labels' %in% input$dimplot.opts),
            cols = colormap[names(x = colormap) %in% unique(x = app.env$object[[input$metacolor.query, drop = TRUE]])],
            repel = TRUE,
            reduction = "umap.proj"
          )[[1]] +
            xlim(app.env$plot.ranges[[1]]) +
            ylim(app.env$plot.ranges[[2]]) +
            labs(x = "UMAP 1", y = "UMAP 2") +
            if (isFALSE(x = "legend" %in% input$dimplot.opts) | OversizedLegend(app.env$object[[input$metacolor.query, drop = TRUE]])) NoLegend()
        } else {
          app.env$plots.querydim_df <- NULL
          plots <- list()
          for (i in 1:length(x = input$metacolor.query)) {
            group.var <- gsub(pattern = "^predicted.", replacement = "", x = input$metacolor.query[i])
            colormap <- GetColorMap(object = refs$map)[[group.var]]
            if (!grepl(pattern = "^predicted.", x = input$metacolor.query[i])) {
              colormap <- NULL
            }
            plots[[i]] <- DimPlot(
              object = app.env$object,
              group.by = input$metacolor.query[i],
              label = isTRUE('labels' %in% input$dimplot.opts),
              cols = colormap[names(x = colormap) %in% unique(x = app.env$object[[input$metacolor.query[i], drop = TRUE]])],
              repel = TRUE,
              reduction = "umap.proj"
            ) + xlim(app.env$plot.ranges[[1]]) +
              ylim(app.env$plot.ranges[[2]]) +
              labs(x = "UMAP 1", y = "UMAP 2") +
              if (isFALSE(x = "legend" %in% input$dimplot.opts) | OversizedLegend(app.env$object[[input$metacolor.query[i], drop = TRUE]])) NoLegend()
          }
          wrap_plots(plots, nrow = 1)
        }
      }
    }
  })
  
  output$objdim_hover_box <- renderUI({
    if (!is.null(x = app.env$plots.objdim_df)) {
      hover <- input$objdim_hover_location
      df <- app.env$plots.objdim_df
      if (!is.null(x = hover)){
        hover[['mapping']] <- setNames(object = as.list(x = colnames(x = app.env$plots.objdim_df)[1:2]), nm = c('x', 'y'))
      }
      point <- nearPoints(
        df = df,
        coordinfo = hover,
        threshold = 10,
        maxpoints = 1,
        addDist = TRUE
      )
      if (nrow(x = point) == 0) {
        return(NULL)
      }
      hovertext <- do.call(
        what = paste0,
        args = as.list(c(
          paste0("<b>", point[[input$metacolor.query]], "</b><br>"),
          if (grepl(pattern = "^predicted.", x = input$metacolor.query)) {
            paste0(
              "<i>prediction score</i>: <span>",
              format(
                x = round(x = point[[paste0(input$metacolor.query,'.score')]], digits = 2),
                nsmall = 2
              ),
              "</span><br>"
            )
          }
        ))
      )
      wellPanel(
        style = HoverBoxStyle(x = hover$coords_css$x, y = hover$coords_css$y),
        p(HTML(text = hovertext))
      )
    } else if (!is.null(x = app.env$plots.refdim_df)) {
      hover <- input$objdim_hover_location
      df <- app.env$plots.refdim_df
      if (!is.null(x = hover)){
        hover[['mapping']] <- setNames(object = as.list(x = colnames(x = app.env$plots.refdim_intro_df)[1:2]), nm = c('x', 'y'))
      }
      point <- nearPoints(
        df = df,
        coordinfo = hover,
        threshold = 10,
        maxpoints = 1,
        addDist = TRUE
      )
      if (nrow(x = point) == 0) {
        return(NULL)
      }
      hovertext <- do.call(
        what = paste0,
        args = as.list(c(
          paste0("<b>", point[[input$metacolor.ref]], "</b><br>"),
          sapply(X = setdiff(metadata.annotate, input$metacolor.ref), FUN = function(md) {
            paste0("<span>", md, "</span>: <i>", point[[md]], "</i><br>")
          })
        ))
      )
      wellPanel(
        style = HoverBoxStyle(x = hover$coords_css$x, y = hover$coords_css$y),
        p(HTML(text = hovertext))
      )
    }
  })
  output$querydim_hover_box <- renderUI({
    if (!is.null(x = app.env$plots.querydim_df)) {
      hover <- input$querydim_hover_location
      df <- app.env$plots.querydim_df
      if (!is.null(x = hover)){
        hover[['mapping']] <- setNames(object = as.list(x = colnames(x = app.env$plots.querydim_df)[1:2]), nm = c('x', 'y'))
      }
      point <- nearPoints(
        df = df,
        coordinfo = hover,
        threshold = 10,
        maxpoints = 1,
        addDist = TRUE
      )
      if (nrow(x = point) == 0) {
        return(NULL)
      }
      hovertext <- do.call(
        what = paste0,
        args = as.list(c(
          paste0("<b>", point[[input$metacolor.query]], "</b><br>"),
          if (grepl(pattern = "^predicted.", x = input$metacolor.query)) {
            paste0(
              "<i>prediction score</i>: <span>",
              format(
                x = round(x = point[[paste0(input$metacolor.query,'.score')]], digits = 2),
                nsmall = 2
              ),
              "</span><br>"
            )
          }
        ))
      )
      wellPanel(
        style = HoverBoxStyle(x = hover$coords_css$x, y = hover$coords_css$y),
        p(HTML(text = hovertext))
      )
    }
  })
  output$evln <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      if (isTRUE(x = do.bridge)){
        avail <- c(
          paste0(
            Key(object = app.env$object[[app.env$gene.assay]]),
            rownames(x = app.env$object[[app.env$gene.assay]])
          ),
          colnames(x = app.env$object[[]])
        )
      } else {
        DefaultAssay(app.env$object) <- "refAssay"
        avail <- c(
          paste0(
            Key(object = app.env$object[["refAssay"]]),
            rownames(x = app.env$object)
          ),
          colnames(x = app.env$object[[]])
        )
      }
      # prediction assays
      prediction.names <- unlist(x = lapply(
        X = app.env$metadataxfer,
        FUN = function(x) {
          assay <- paste0("prediction.score.", x)
          pred <- rep(x = x, times = nrow(x = app.env$object[[assay]]))
          names(x = pred) <- paste0(
            Key(object = app.env$object[[assay]]),
            rownames(x = app.env$object[[assay]])
          )
          return(pred)
        })
      )
      max.pred.names <- paste0("predicted.", app.env$metadataxfer, ".score")
      avail <- c(avail, names(x = prediction.names))
      print("AVAIL")
      print(head(avail))
      print(app.env$feature)
      print(app.env$feature %in% avail)
      if (do.adt) {
        avail <- c(
          avail,
          paste0(
            Key(object = app.env$object[[adt.key]]),
            rownames(x = app.env$object[[adt.key]])
          )
        )
      }
      
      if (app.env$feature %in% avail) {
        if (app.env$feature == "mapping.score" && !resolved(x = app.env$mapping.score)) {
          ggplot() +
            annotate("text", x = 4, y = 25, size=8, label = "Mapping score still computing ... ") +
            theme_void()
        } else {
          title <- ifelse(
            test = grepl(pattern = '^refassay_', x = app.env$feature),
            yes = gsub(pattern = '^refassay_', replacement = '', x = app.env$feature),
            no = app.env$feature
          )
          if (app.env$feature %in% names(x = prediction.names)) {
            pred <- strsplit(x = app.env$feature, split = "_")[[1]][2]
            group <- prediction.names[app.env$feature]
            title <- paste0("Prediction Score (", group, ") ", pred)
          }
          if (app.env$feature %in% max.pred.names) {
            pred <- gsub(pattern = "predicted.", replacement = "", x = app.env$feature)
            pred <- gsub(pattern = ".score", replacement = "", x = pred)
            title <- paste0("Max Prediction Score - ", pred)
          }
          print("FEATURE FOR FEATURE PLOT")
          print(app.env$feature)
          VlnPlot(
            object = app.env$object,
            features = app.env$feature,
            group.by = input$metagroup,
            pt.size = ifelse(
              test = input$check.featpoints,
              yes = 0,
              no = Seurat:::AutoPointSize(data = app.env$object)
            )
          ) +
            ggtitle(label = title) +
            NoLegend()
        }
      }
    }
  })
  output$edim <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      print("RENDERING DIM PLOT")
      palettes <- list(
        c("lightgrey", "blue"),
        c('lightgrey', 'darkred')
      )
      if (isTRUE(x = do.bridge)){
        names(x = palettes) <- c(
          Key(object = app.env$object[[app.env$gene.assay]]),
          'md_'
        )
      } else{
        DefaultAssay(app.env$object) <- "refAssay"
        names(x = palettes) <- c(
          Key(object = app.env$object[["refAssay"]]),
          'md_'
        )
      }
      if (do.adt) {
        palettes[[Key(object = app.env$object[[adt.key]])]] <-  c('lightgrey', 'darkgreen')
      }
      # prediction assays
      prediction.names <- unlist(x = lapply(
        X = app.env$metadataxfer,
        FUN = function(x) {
          assay <- paste0("prediction.score.", x)
          pred <- rep(x = x, times = nrow(x = app.env$object[[assay]]))
          names(x = pred) <- paste0(
            Key(object = app.env$object[[assay]]),
            rownames(x = app.env$object[[assay]])
          )
          return(pred)
        })
      )
      max.pred.names <- paste0("predicted.", app.env$metadataxfer, ".score")
      md <- c(colnames(x = app.env$object[[]]), names(x = prediction.names))
      feature.key <- if (app.env$feature %in% md) {
        'md_'
      } else {
        paste0(
          unlist(x = strsplit(x = app.env$feature, split = '_'))[1],
          '_'
        )
      }
      pal.use <- palettes[[feature.key]]
      if (!is.null(x = pal.use)) {
        if (app.env$feature == "mapping.score" && !resolved(x = app.env$mapping.score)) {
          ggplot() +
            annotate("text", x = 4, y = 25, size=8, label = "Mapping score still computing ... ") +
            theme_void()
        } else {
          title <- ifelse(
            test = grepl(pattern = '^refassay_', x = app.env$feature),
            yes = gsub(pattern = '^refassay_', replacement = '', x = app.env$feature),
            no = app.env$feature
          )
          if (app.env$feature %in% names(x = prediction.names)) {
            pred <- strsplit(x = app.env$feature, split = "_")[[1]][2]
            group <- prediction.names[app.env$feature]
            title <- paste0("Prediction Score (", group, ") ", pred)
          }
          if (app.env$feature %in% max.pred.names) {
            pred <- gsub(pattern = "predicted.", replacement = "", x = app.env$feature)
            pred <- gsub(pattern = ".score", replacement = "", x = pred)
            title <- paste0("Max Prediction Score - ", pred)
          }
          print("FEATURE FOR FEATURE PLOT")
          print(app.env$feature)
          suppressWarnings(expr = FeaturePlot(
            object = app.env$object,
            features = app.env$feature,
            cols = pal.use,
            reduction = "umap.proj"
          )) + xlim(app.env$plot.ranges[[1]]) +
            ylim(app.env$plot.ranges[[2]]) +
            ggtitle(label = title)
        }
      }
    }
  })
  
  # output$motifvln <- renderPlot(expr = {
  #   if (!is.null(x = app.env$object)) {
  #     print("RENDERING VLN PLOT FOR MOTIF")
  #     avail <- c(paste0(Key(object = app.env$object[[app.env$default.assay]]), 
  #                       rownames(x = app.env$object[[app.env$default.assay]])), colnames(x = app.env$object[[]]))
  #     prediction.names <- unlist(x = lapply(X = app.env$metadataxfer, 
  #                                           FUN = function(x) {
  #                                             assay <- paste0("prediction.score.", x)
  #                                             pred <- rep(x = x, times = nrow(x = app.env$object[[assay]]))
  #                                             names(x = pred) <- paste0(Key(object = app.env$object[[assay]]), 
  #                                                                       rownames(x = app.env$object[[assay]]))
  #                                             return(pred)
  #                                           }))
  #     print("PREDICTION NAMES")
  #     print(head(prediction.names))
  #     max.pred.names <- paste0("predicted.", app.env$metadataxfer, 
  #                              ".score")
  #     avail <- c(avail, names(x = prediction.names))
  #     print("AVAIL")
  #     print(head(avail))
  #     print("second if statement")
  #     print(app.env$motif.feature)
  #     print(app.env$motif.feature %in% avail)
  #     if (app.env$motif.feature %in% avail) {
  #       print("went into the second if statement")
  #       if (app.env$motif.feature == "mapping.score" && !resolved(x = app.env$mapping.score)) {
  #         ggplot() + annotate("text", x = 4, y = 25, 
  #                             size = 8, label = "Mapping score still computing ... ") + 
  #           theme_void()
  #       }
  #       else {
  #         title <- ifelse(test = grepl(pattern = "^motif_", 
  #                                      x = app.env$motif.feature), yes = gsub(pattern = "^motif_", 
  #                                                                                replacement = "", x = app.env$motif.feature), no = app.env$motif.feature)
  #         if (app.env$motif.feature %in% names(x = prediction.names)) {
  #           print("went into the third if statement")
  #           pred <- strsplit(x = app.env$motif.feature, split = "_")[[1]][2]
  #           group <- prediction.names[app.env$motif.feature]
  #           title <- paste0("Prediction Score (", group, 
  #                           ") ", pred)
  #         }
  #         if (app.env$motif.feature %in% max.pred.names) {
  #           print("went into the fourth print statement")
  #           pred <- gsub(pattern = "predicted.", replacement = "", 
  #                        x = app.env$motif.feature)
  #           pred <- gsub(pattern = ".score", replacement = "", 
  #                        x = pred)
  #           title <- paste0("Max Prediction Score - ", 
  #                           pred)
  #         }
  #         print("making vln plot")
  #         print("MOTIF FOR VLN PLOT")
  #         print(app.env$motif.feature)
  #         VlnPlot(object = app.env$object, features = app.env$motif.feature, 
  #                 group.by = input$metagroup.motif, pt.size = ifelse(test = input$check.featpoints, 
  #                                                                    yes = 0, no = Seurat:::AutoPointSize(data = app.env$object))) + 
  #           ggtitle(label = title) + NoLegend()
  #       }
  #     }
  #   }
  # })
  
  output$motifdim <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      print("MAKING MOTIF DIM PLOT")
      palettes <- list(c("lightgrey", "blue"), c("lightgrey", 
                                                 "darkred"))
      names(x = palettes) <- c(Key(object = app.env$object[[app.env$default.assay]]), 
                               "md_")
      prediction.names <- unlist(x = lapply(X = app.env$metadataxfer, 
                                            FUN = function(x) {
                                              assay <- paste0("prediction.score.", x)
                                              pred <- rep(x = x, times = nrow(x = app.env$object[[assay]]))
                                              names(x = pred) <- paste0(Key(object = app.env$object[[assay]]), 
                                                                        rownames(x = app.env$object[[assay]]))
                                              return(pred)
                                            }))
      max.pred.names <- paste0("predicted.", app.env$metadataxfer, 
                               ".score")
      md <- c(colnames(x = app.env$object[[]]), names(x = prediction.names))
      print("MD")
      print(head(md))
      feature.key <- if (app.env$motif.feature %in% md) {
        "md_"
      }
      else {
        paste0(unlist(x = strsplit(x = app.env$motif.feature, 
                                   split = "_"))[1], "_")
      }
      print("FEATURE KEY:")
      print(head(feature.key))
      pal.use <- palettes[[feature.key]]
      print("PAL.USE:")
      print(head(pal.use))
      if (!is.null(x = pal.use)) {
        if (app.env$motif.feature == "mapping.score" && !resolved(x = app.env$mapping.score)) {
          ggplot() + annotate("text", x = 4, y = 25, 
                              size = 8, label = "Mapping score still computing ... ") + 
            theme_void()
        }
        else {
          title <- ifelse(test = grepl(pattern = "^motif_", 
                                       x = app.env$feature), yes = gsub(pattern = "^motif_", 
                                                                        replacement = "", x = app.env$motif.feature), no = app.env$motif.feature)
          if (app.env$motif.feature %in% names(x = prediction.names)) {
            pred <- strsplit(x = app.env$motif.feature, split = "_")[[1]][2]
            group <- prediction.names[app.env$motif.feature]
            title <- paste0("Prediction Score (", group, 
                            ") ", pred)
          }
          if (app.env$motif.feature %in% max.pred.names) {
            pred <- gsub(pattern = "predicted.", replacement = "", 
                         x = app.env$motif.feature)
            pred <- gsub(pattern = ".score", replacement = "", 
                         x = pred)
            title <- paste0("Max Prediction Score - ", 
                            pred)
          }
          suppressWarnings(expr = FeaturePlot(object = app.env$object, 
                                              features = app.env$motif.feature, cols = pal.use, 
                                              min.cutoff = 'q10', max.cutoff = 'q90',  
                                              reduction = "umap.proj")) + xlim(app.env$plot.ranges[[1]]) + 
            ylim(app.env$plot.ranges[[2]]) + ggtitle(label = title)
        }
      }
    }
  })
  # Messages
  output$message <- renderUI(expr = {
    p(HTML(text = paste(app.env$messages, collapse = "<br />")))
  })
  output$containerid <- renderUI(expr = {
    p(HTML(text = paste(
      paste("debug ID:", Sys.info()[["nodename"]]),
      paste('Azimuth version:', packageVersion(pkg = 'Azimuth')),
      paste('Seurat version:', packageVersion(pkg = 'Seurat')),
      paste('Reference version:', ReferenceVersion(object = refs$map)),
      sep = "<br />"
    )))
  })
  output$text.cellsremain <- renderText(expr = {
    if (!is.null(x = isolate(app.env$object))) {
      ncount <- paste0('nCount_', DefaultAssay(object = isolate(app.env$object)))
      nfeature <- paste0('nFeature_', DefaultAssay(object = isolate(app.env$object)))
      cells.use <- isolate(app.env$object)[[ncount, drop = TRUE]] >= input$num.ncountmin &
        isolate(app.env$object)[[ncount, drop = TRUE]] <= input$num.ncountmax &
        isolate(app.env$object)[[nfeature, drop = TRUE]] >= input$num.nfeaturemin &
        isolate(app.env$object)[[nfeature, drop = TRUE]] <= input$num.nfeaturemax
      if (any(grepl(pattern = mito.pattern, x = rownames(x = isolate(app.env$object))))) {
        cells.use <- cells.use &
          isolate(app.env$object)[[mt.key, drop = TRUE]] >= input$minmt &
          isolate(app.env$object)[[mt.key, drop = TRUE]] <= input$maxmt
      }
      paste(sum(cells.use), "cells remain after current filters")
    }
  })
  output$text.dladt <- renderText(
    expr = {
      c(
        "imputed.assay <- readRDS('azimuth_impADT.Rds')",
        "object <- object[, Cells(imputed.assay)]",
        "object[['impADT']] <- imputed.assay"
      )
    },
    sep = "\n"
  )
  output$text.dlumap <- renderText(
    expr = {
      c(
        "projected.umap <- readRDS('azimuth_umap.Rds')",
        "object <- object[, Cells(projected.umap)]",
        "object[['umap.proj']] <- projected.umap"
      )
    },
    sep = "\n"
  )
  output$text.dlpred <- renderText(
    expr = {
      c(
        "predictions <- read.delim('azimuth_pred.tsv', row.names = 1)",
        "object <- AddMetaData(",
        "\tobject = object,",
        "\tmetadata = predictions)"
      )
    },
    sep = "\n"
  )
  output$text.dlall <- renderText(
    expr = {
      c(
        "object <- AddAzimuthResults(object, filename = 'azimuth_results.Rds')"
      )
    },
    sep = "\n"
  )
  # Tables
  output$table.qc <- renderTable(
    expr = {
      if (!is.null(x = isolate(app.env$object))) {
        qc <- paste0(c('nCount_', 'nFeature_'), app.env$default.assay)
        tbl <- apply(X = isolate(app.env$object)[[qc]], MARGIN = 2, FUN = quantile)
        tbl <- as.data.frame(x = tbl)
        if (isTRUE(x = do.bridge)){
          colnames(x = tbl) <- c('Fragments per cell', 'Peaks detected per cell')
        } else{
          colnames(x = tbl) <- c('nUMI per cell', 'Genes detected per cell')
        }
        if (mt.key %in% colnames(x = isolate(app.env$object)[[]])) {
          tbl[, 3] <- quantile(x = isolate(app.env$object)[[mt.key, drop = TRUE]])
          colnames(x = tbl)[3] <- 'Mitochondrial percentage per cell'
        }
        t(x = tbl)
      }
    },
    rownames = TRUE
  )
  output$biomarkers <- renderDT(
    expr = {
      if (!is.null(x = app.env$diff.expr[[paste(app.env$gene.assay, input$markerclustersgroup, sep ="_")]])) {
        RenderDiffExp(
          diff.exp =  app.env$diff.expr[[paste(app.env$gene.assay, input$markerclustersgroup, sep ="_")]],
          groups.use = input$markerclusters,
          n = Inf
        )
      }
    },
    selection = 'single',
    options = list(dom = 't')
  )
  output$adtbio <- renderDT(
    expr = {
      if (!is.null(x = app.env$diff.expr[[paste(adt.key, input$markerclustersgroup, sep = "_")]])) {
        RenderDiffExp(
          diff.exp = app.env$diff.expr[[paste(adt.key, input$markerclustersgroup, sep = "_")]],
          groups.use = input$markerclusters,
          n = Inf
        )
      }
    },
    selection = 'single',
    options = list(dom = 't')
  )
  output$motifs <- renderDT(
    expr = {
      print("RENDERING MOTIF DATA TABLE")
      print("INPUT MARKER CLUSTERS: ")
      print(input$markerclustersgroup)
      print(paste(app.env$default.assay, input$markerclustersgroup, sep ="_"))
      print(head(app.env$motif.diff.expr[[paste(app.env$default.assay, input$markerclustersgroup.motif, sep ="_")]]))
      print(input$markerclusters)
      if (!is.null(x = app.env$motif.diff.expr[[paste(app.env$default.assay, input$markerclustersgroup.motif, sep ="_")]])) {
        RenderDiffMotifExp(
          diff.exp =  app.env$motif.diff.expr[[paste(app.env$default.assay, input$markerclustersgroup.motif, sep ="_")]],
          groups.use = input$markerclusters.motif,
          n = Inf
        )
      }
    },
    selection = 'single',
    options = list(dom = 't')
  )
  output$metadata.table <- renderTable(
    expr = {
      if (!is.null(x = app.env$object)) {
        CategoryTable(
          object = app.env$object,
          category.1 = input$metarow,
          category.2 = input$metacol,
          percentage = (input$radio.pct == "Percentage")
        )
      }
    },
    rownames = TRUE
  )
  output$metadata.heatmap <- renderPlotly({
    table <- CategoryTable(
      object = app.env$object,
      category.1 = input$metarow,
      category.2 = input$metacol,
      percentage = (input$radio.pct == "Percentage")
    )
    table <- as.matrix(table)
    plot_ly(x = colnames(table), y = rownames(table), z = table, type = 'heatmap',
            height='1000px')
  })
  # Downloads
  output$dlumap <- downloadHandler(
    filename = paste0(tolower(x = app.title), '_umap.Rds'),
    content = function(file) {
      if (!is.null(x = app.env$object)) {
        if ('umap.proj' %in% Reductions(object = app.env$object)) {
          saveRDS(object = app.env$object[['umap.proj']], file = file)
        }
      }
    }
  )
  output$dladt <- downloadHandler(
    filename = paste0(tolower(x = app.title), '_impADT.Rds'),
    content = function(file) {
      if (!is.null(x = app.env$object)) {
        if ('impADT' %in% Assays(object = app.env$object)) {
          saveRDS(object = app.env$object[['impADT']], file = file)
        }
      }
    }
  )
  output$dlpred <- downloadHandler(
    filename = paste0(tolower(x = app.title), '_pred.tsv'),
    content = function(file) {
      req <- paste0("predicted.", c(app.env$metadataxfer, paste0(app.env$metadataxfer, ".score")))
      if (resolved(x = app.env$mapping.score)) {
        req <- c(req, 'mapping.score')
      }
      if (all(req %in% colnames(x = app.env$object[[]]))) {
        pred.df <- app.env$object[[req]]
        if (resolved(x = app.env$mapping.score)) {
          pred.df$mapping.score <- value(app.env$mapping.score)
        }
        pred.df <- cbind(cell = rownames(x = pred.df), pred.df)
        write.table(
          x = pred.df,
          file = file,
          quote = FALSE,
          row.names = FALSE,
          col.names = TRUE,
          sep = '\t'
        )
      }
    }
  )
  output$dlall <- downloadHandler(
    filename = paste0(tolower(x = app.title), '_results.Rds'),
    content = function(file) {
      results <- list()
      if (!is.null(x = app.env$object)) {
        if ('impADT' %in% Assays(object = app.env$object)) {
          results$impADT <- app.env$object[['impADT']]
        }
        if ('umap.proj' %in% Reductions(object = app.env$object)) {
          results$umap <- app.env$object[['umap.proj']]
        }
      }
      
      req <- paste0("predicted.", c(app.env$metadataxfer, paste0(app.env$metadataxfer, ".score")))
      if (resolved(x = app.env$mapping.score)) {
        req <- c(req, 'mapping.score')
      }
      if (all(req %in% colnames(x = app.env$object[[]]))) {
        pred.df <- app.env$object[[req]]
        if (resolved(x = app.env$mapping.score)) {
          pred.df$mapping.score <- value(app.env$mapping.score)
        }
        pred.df <- cbind(cell = rownames(x = pred.df), pred.df)
        results$pred.df <- pred.df
      }
      
      saveRDS(results, file = file)
    }
  )
  output$dlscript <- downloadHandler(
    filename = paste0(tolower(x = app.title), '_analysis.R'),
    content = function(file) {
      template <- readLines(con = system.file(
        file.path('resources', 'template.R'),
        package = 'Azimuth'
      ))
      template <- paste(template, collapse = '\n')
      e <- new.env()
      # e$ref.uri <- 'https://seurat.nygenome.org/references/pbmc/'
      e$ref.uri <- getOption(
        x = 'Azimuth.app.refuri',
        default = getOption(x = 'Azimuth.app.reference')
      )
      e$path <- input$file$name
      e$mito.pattern <- getOption(x = 'Azimuth.app.mito', default = '^MT-')
      e$mito.key <- mt.key
      e$ncount.max <- input$num.ncountmax
      e$ncount.min <- input$num.ncountmin
      e$nfeature.max <- input$num.nfeaturemax
      e$nfeature.min <- input$num.nfeaturemin
      e$mito.max <- input$maxmt
      e$mito.min <- input$minmt
      e$sct.ncells <- getOption(x = 'Azimuth.sct.ncells')
      e$sct.nfeats <- getOption(x = 'Azimuth.sct.nfeats')
      e$ntrees <- getOption(x = 'Azimuth.map.ntrees')
      e$ndims <- getOption(x = "Azimuth.map.ndims")
      e$adt.key <- adt.key
      e$do.adt <- do.adt
      e$metadataxfer <- app.env$metadataxfer
      if (length(x = e$metadataxfer == 1)) {
        e$metadataxfer <- paste0("\"", e$metadataxfer, "\"")
      }
      e$plotgene <- getOption(x = 'Azimuth.app.default_gene')
      e$plotadt <- getOption(x = 'Azimuth.app.default_adt')
      writeLines(text = str_interp(string = template, env = e), con = file)
    }
  )
  # render UI elements that depend on arguments
  output$refdescriptor <- renderText(
    expr = eval(expr = HTML(getOption(x = "Azimuth.app.refdescriptor")))
  )
  output$welcomebox <- renderUI(
    expr = eval(expr = parse(text = getOption(x = "Azimuth.app.welcomebox")))
  )
  
  # render popup UI elements
  onclick('overlap_popup', showModal(modalDialog(
    title = "Overlap QC",
    div(
      paste(
        "In order to conduct bridge integration for ATAC data without uploading a large ", 
        "fragment file, we requantify the ATAC query peaks to match the multiomic bridge ", 
        "based on the overlap between each query peak to a bridge peak and rename the query ", 
        "peak to the bridge peak with highest overlap. The box color corresponds to the following bins: "
      ),
      tags$ul(list(
        tags$li(paste0("0% to ", getOption(x = "Azimuth.map.panchorscolors")[1], "%: Likely problematic (red)")),
        tags$li(paste0(getOption(x = "Azimuth.map.panchorscolors")[1], "% to ", getOption(x = "Azimuth.map.panchorscolors")[2], "%: Possibly problematic (yellow)")),
        tags$li(paste0(getOption(x = "Azimuth.map.panchorscolors")[2], "% to 100%: Likely successful (green)"))
      )),
      tags$h4("Caveats"),
      paste0(
        "A high percentage of overlap is expected if the query ATAC data and bridge ATAC data ", 
        "were processed with the same versions of Cell Ranger and means that there will ", 
        "likely be little loss of information by using this overlap renaming process. ", 
        "The mapping can still be sucessesful if this value has a low percentage, but downstream gene expression", 
        "calculations may be innacurate as this again uses another overlap process to requantify peaks to genes."
      )
    )
  )))
  onclick('panchors_popup', showModal(modalDialog(
    title = "Anchor QC",
    div(
      paste(
        "The Azimuth reference-mapping procedure first identifies a set of 'anchors', ",
        "or pairwise correspondences between cells predicted to be in a similar biological state, ",
        "between query and reference datasets. Here we report the percentage of query cells ",
        "participating in an anchor correspondence. The box color corresponds to the following bins: "
      ),
      tags$ul(list(
        tags$li(paste0("0% to ", getOption(x = "Azimuth.map.panchorscolors")[1], "%: Likely problematic (red)")),
        tags$li(paste0(getOption(x = "Azimuth.map.panchorscolors")[1], "% to ", getOption(x = "Azimuth.map.panchorscolors")[2], "%: Possibly problematic (yellow)")),
        tags$li(paste0(getOption(x = "Azimuth.map.panchorscolors")[2], "% to 100%: Likely successful (green)"))
      )),
      tags$h4("Caveats"),
      paste0(
        "If the query dataset consists of a homogeneous group of cells, or if the ",
        "query dataset contains cells from multiple batches (which would be corrected ",
        "by Azimuth), this metric may return a low value even in cases where mapping is ",
        "successful. Users in these cases should check results carefully. In particular, ",
        "we encourage users to verify identified differentially expressed marker genes for annotated cell types."
      )
    )
  )))
  onclick('mappingqcstat_popup', showModal(modalDialog(
    title = "Cluster Preservation",
    div(
      tags$h4("Overview"),
      paste0(
        "For each query dataset, we downsample to at most 5,000 cells, and perform an ",
        "unsupervised clustering. This score reflects the preservation of the unsupervised ",
        "cluster structure, and is based on the entropy of unsupervised cluster labels in ",
        "each query cell's local neighborhood after mapping. Scores are scaled from 0 (poor) to 5 (best)"
      ),
      tags$ul(list(
        tags$li(paste0("0 to ", getOption(x = "Azimuth.map.postmapqccolors")[1], ": Likely problematic (red)")),
        tags$li(paste0(getOption(x = "Azimuth.map.postmapqccolors")[1], " to ", getOption(x = "Azimuth.map.postmapqccolors")[2], ": Possibly problematic (yellow)")),
        tags$li(paste0(getOption(x = "Azimuth.map.postmapqccolors")[2], " to 5: Likely successful (green)"))
      )),
      tags$h4("Caveats"),
      paste0(
        "This metric relies on the unsupervised clustering representing corresponding to ",
        "biologically distinct cell states. If the query dataset consists of a homogeneous ",
        "group of cells, or if the query dataset contains cells from multiple batches ",
        "(which would be corrected by Azimuth), this metric may return a low value even ",
        "in cases where mapping is successful. Users in these cases should check results ",
        "carefully. In particular, we encourage users to verify identified differentially ",
        "expressed marker genes for annotated cell types."
      ),
      # tags$h4("Details"),
      # paste0(
      #   "To compute the mapping statistic, we first randomly downsample the ",
      #   "query to ", getOption(x = "Azimuth.map.postmapqcds"), " cells for ",
      #   "computational efficiency. We then compute an independent unsupervised ",
      #   "clustering on the query. Using these cluster IDs, we then examine the ",
      #   "neighborhoods of each cell in query PCA space and also in the mapped ",
      #   "(projected) space. We compute an entropy of cluster labels and then ",
      #   "take the mean entropy averaged over each cluster in both spaces. For ",
      #   "each cluster we take the difference and report a single statistic as ",
      #   "the median -log2 of these values, clipped to range between 0 and 5.",
      #   "For the exact implementation details, please see the ",
      #   "ClusterPreservationScore function in the azimuth github repo."
      # ),
      
    )
  )))
}










#' Server function for the mapping app for bridge integration 
#'
#' @param input,output,session Required Shiny app server parameters
#'
#' @return The shiny server logic
#'
#' @importFrom DT dataTableProxy renderDT selectRows
#' @importFrom future future plan resolved value
#' @importFrom ggplot2 annotate geom_hline ggtitle scale_colour_hue
#' theme_void xlab layer_scales xlim ylim ggplot aes geom_point theme
#' element_blank element_rect labs
#' @importFrom googlesheets4 gs4_auth gs4_get sheet_append
#' @importFrom methods slot slot<- new
#' @importFrom presto wilcoxauc
#' @importFrom SeuratObject AddMetaData Assays Cells DefaultAssay Embeddings
#' GetAssayData Idents Idents<- Key RenameCells Reductions Tool SetAssayData
#' VariableFeatures
#' @importFrom Seurat DimPlot FeaturePlot FindNeighbors FindTransferAnchors
#' IntegrateEmbeddings MappingScore NoLegend PercentageFeatureSet
#' RunUMAP TransferData SCTransform VlnPlot LabelClusters
#' FindBridgeTransferAnchors MapQuery NormalizeData
#' @importFrom shiny downloadHandler observeEvent isolate Progress
#' reactiveValues renderPlot renderTable renderText removeUI setProgress
#' safeError updateNumericInput updateSelectizeInput updateCheckboxInput updateTextAreaInput
#' withProgress renderUI onStop showNotification wellPanel nearPoints insertUI
#' modalDialog showModal getDefaultReactiveDomain
#' @importFrom shinydashboard menuItem renderMenu renderValueBox
#' sidebarMenu valueBox
#' @importFrom shinyjs addClass enable disable hide removeClass show onclick
#' disable
#' @importFrom stringr str_interp str_trim str_split
#' @importFrom patchwork wrap_plots
#' @importFrom stats na.omit quantile setNames median
#' @importFrom utils write.table packageVersion
#' @importFrom plotly plotlyOutput renderPlotly toWebGL ggplotly plot_ly
#' @importFrom Signac CollapseToLongestTranscript GetTranscripts GetGRangesFromEnsDb RunTFIDF 
#' RunChromVar AddMotifs
#' @importFrom TFBSTools getMatrixSet
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom IRanges findOverlaps
#' @importFrom data.table as.data.table
#'
#' @keywords internal
#'
AzimuthBridgeServer <- function(input, output, session) {
  hide(id = "legend")
  disable(id = "metacolor.ref")
  if (is.null(x = getOption(x = "Azimuth.app.demodataset"))) {
    hide(id = "demobuttons")
  }
  mt.key <- "percent.mt"
  mito.pattern <- getOption(x = "Azimuth.app.mito", default = "^MT-")
  do.adt <- isTRUE(x = getOption(x = "Azimuth.app.do_adt", 
                                 default = TRUE))
  adt.key <- "impADT"
  n.trees <- getOption(x = "Azimuth.map.ntrees")
  app.env <- reactiveValues(adt.features = character(length = 0L), 
                            anchors = NULL, annotations = NULL, bridge_anchors = FALSE, 
                            clusterpreservationqc = NULL, counts = NULL, chromatin_assay_1 = NULL,
                            chromatin_assay_2 = NULL, default.chromvar.feature = NULL, chromvar.feature = "", chromvar.features = character(length = 0L), 
                            chromvar.diff.expr = list(), demo = FALSE, 
                            demo.inputs = NULL, demo.tracker = NULL, demo.files = NULL, 
                            default.assay = NULL, default.feature = NULL, default.metadata = NULL, 
                            diff.exp = list(), gene.assay = NULL, feature = "", features = character(length = 0L), 
                            mapping.score = NULL, messages = "Upload a file", nanchors = 0L, 
                            ncellsupload = 0L, ncellspreproc = 0L, object = NULL, 
                            metadata.cont = character(length = 0L), o_hits = NULL, scorefeatures = character(length = 0L), 
                            plot.ranges = list(), plots.refdim_df = NULL, plots.refdim_intro_df = NULL, 
                            plots.objdim_df = NULL, plots.querydim_df = NULL, requantified_multiome = NULL, 
                            requantified_genes = NULL, fresh.plot = TRUE, 
                            singlepred = NULL, emptyref = NULL, merged = NULL, metadata.discrete = NULL, 
                            metadata.notransfer = NULL, disable = FALSE, query.names = character())
  react.env <- reactiveValues(no = FALSE, anchors = FALSE, annotations = FALSE,
                              biomarkers = FALSE, bridge_anchors = FALSE, cluster.score = FALSE, chromvar = FALSE, chromvar.features = FALSE,
                              chromatin_assay_1 = FALSE, features = FALSE, get.feature = FALSE, get.chromvar.feature = FALSE, 
                              map = FALSE, mapquery = FALSE, markers = FALSE, metadata = FALSE, mt = NULL, 
                              xferopts = FALSE, path = NULL, progress = NULL, plot.qc = FALSE, requantify_multiome = FALSE, 
                              requantify_genes = FALSE, qc = FALSE, score = FALSE, sctransform = FALSE, start = numeric(length = 0L), 
                              transform = FALSE)
  if (isTRUE(x = do.adt)) { # will i need this
    output$imputedlabel <- renderUI(expr = h3("Imputed protein biomarkers"))
  }
  else {
    for (id in c("imputedinput", "imputedtable", "imputeddl")) { 
      removeUI(selector = paste0("#", id), immediate = TRUE)
    }
    for (id in c("featureinput", "scoreinput")) {
      removeClass(id = id, class = "thirds")
      addClass(id = id, class = "halves")
    }
    removeClass(id = "biotable", class = "halves")
    addClass(id = "biotable", class = "fulls")
  }
  ResetEnv <- function() {
    print("resetting...")
    app.env$disable <- TRUE
    output$menu2 <- NULL
    react.env$plot.qc <- FALSE
    app.env$messages <- NULL
    output$valubox.upload <- NULL
    output$valuebox.preproc <- NULL
    output$valuebox.mapped <- NULL
    output$valuebox_panchors <- NULL
    output$valuebox_mappingqcstat <- NULL
    app.env$emptyref <- NULL
    app.env$merged <- NULL
    app.env$metadata.discrete <- NULL
    disable(id = "map")
    hide(selector = ".rowhide")
  }
  motif.proxy <- dataTableProxy(outputId = "motifs")
  rna.proxy <- dataTableProxy(outputId = "biomarkers")
  adt.proxy <- dataTableProxy(outputId = "adtbio")
  logging <- all(vapply(X = paste0("Azimuth.app.google", c("sheet", 
                                                           "token", "tokenemail")), FUN = function(x) {
                                                             return(!is.null(x = getOption(x = x)))
                                                           }, FUN.VALUE = logical(length = 1L)))
  googlesheet <- NULL
  if (logging) {
    try(expr = {
      gs4_auth(email = getOption(x = "Azimuth.app.googletokenemail"), 
               cache = getOption(x = "Azimuth.app.googletoken"))
      googlesheet <- gs4_get(ss = getOption(x = "Azimuth.app.googlesheet"))
      app_start_time <- Sys.time()
      app_session_id <- paste0(Sys.info()[["nodename"]], 
                               as.numeric(Sys.time()))
      onStop(fun = function() {
        try(expr = sheet_append(ss = googlesheet, data = data.frame("SESSIONLENGTH", 
                                                                    app_session_id, as.numeric(x = Sys.time() - 
                                                                                                 app_start_time, units = "mins"))))
      })
    })
  }
  if (!is.null(x = googlesheet)) {
    try(expr = sheet_append(ss = googlesheet, data = data.frame("STARTUPTIME", 
                                                                app_session_id, Sys.time())))
    output$menu3 <- renderMenu(expr = {
      sidebarMenu(menuItem(text = "Feedback", tabName = "tab_feedback", 
                           icon = icon(name = "comments"), selected = FALSE))
    })
  }
  demos <- getOption("Azimuth.app.demodataset")
  if (!is.null(x = demos)) {
    if (!inherits(x = demos, what = "data.frame")) {
      if (is.null(x = names(x = demos))) {
        if (length(x = demos) > 1) {
          demo.names <- paste0("Demo", 1:length(x = demos))
        }
        else {
          demo.names <- "Load demo dataset"
        }
      }
      else {
        demo.names <- names(x = demos)
      }
      demos <- data.frame(name = demo.names, file = demos)
    }
    app.env$demo.files <- demos$file
    app.env$demo.inputs <- paste0("triggerdemo", 1:nrow(x = demos))
    app.env$demo.tracker <- rep(x = 0, times = nrow(x = demos))
    for (i in 1:nrow(x = demos)) {
      insertUI(selector = "#demobuttons", where = "beforeEnd", 
               immediate = TRUE, ui = actionButton(inputId = paste0("triggerdemo", 
                                                                    i), label = demos$name[i], width = "85%"))
    }
  }
  if (getOption(x = "Azimuth.app.metatableheatmap")) {
    insertUI(selector = "#tablemetadata", where = "beforeEnd", 
             immediate = TRUE, ui = plotlyOutput(outputId = "metadata.heatmap"))
  }
  else {
    insertUI(selector = "#tablemetadata", where = "beforeEnd", 
             immediate = TRUE, ui = tableOutput(outputId = "metadata.table"))
  }
  if (getOption(x = "Azimuth.app.overlayedreference")) {
    insertUI(selector = "#topdim", where = "beforeEnd", immediate = TRUE, 
             ui = box(title = "Mapped Query", checkboxInput(inputId = "legend", 
                                                            label = "Show legend"), checkboxGroupInput(inputId = "label.opts", 
                                                                                                       label = NULL, choiceNames = c("Show labels", 
                                                                                                                                     "Filter cluster labels (size <2%)"), choiceValues = c("labels", 
                                                                                                                                                                                           "filterlabels"), selected = c("labels", "filterlabels"), 
                                                                                                       inline = TRUE), checkboxInput(inputId = "showrefonly", 
                                                                                                                                     label = "View reference only"), selectizeInput(inputId = "metacolor.ref", 
                                                                                                                                                                                    label = "Reference metadata to color by", choices = "", 
                                                                                                                                                                                    multiple = TRUE, ), selectizeInput(inputId = "metacolor.query", 
                                                                                                                                                                                                                       label = "Query metadata to color by", choices = "", 
                                                                                                                                                                                                                       multiple = TRUE, ), div(style = "position:relative", 
                                                                                                                                                                                                                                               plotOutput(outputId = "objdim", hover = hoverOpts(id = "objdim_hover_location", 
                                                                                                                                                                                                                                                                                                 delay = 5, delayType = "debounce", nullOutside = TRUE), 
                                                                                                                                                                                                                                                          height = "750px"), uiOutput("objdim_hover_box")), 
                      width = 12, height = "auto"))
  }
  else {
    insertUI(selector = "#topdim", where = "beforeEnd", immediate = TRUE, 
             ui = box(title = "Reference", checkboxGroupInput(inputId = "dimplot.opts", 
                                                              label = NULL, choiceNames = c("Show labels", 
                                                                                            "Show legend"), choiceValues = c("labels", 
                                                                                                                             "legend"), selected = "legend", inline = TRUE), 
                      selectizeInput(inputId = "metacolor.ref", label = "Metadata to color by", 
                                     choices = "", multiple = TRUE, ), div(style = "position:relative", 
                                                                           plotOutput(outputId = "refdim", hover = hoverOpts(id = "refdim_hover_location", 
                                                                                                                             delay = 5, delayType = "debounce", nullOutside = TRUE)), 
                                                                           uiOutput("refdim_hover_box")), width = 12))
    insertUI(selector = "#bottomdim", where = "beforeEnd", 
             immediate = TRUE, ui = box(title = "Query", selectizeInput(inputId = "metacolor.query", 
                                                                        label = "Metadata to color by", choices = "", 
                                                                        multiple = TRUE, ), div(style = "position:relative", 
                                                                                                plotOutput(outputId = "querydim", hover = hoverOpts(id = "querydim_hover_location", 
                                                                                                                                                    delay = 5, delayType = "debounce", nullOutside = TRUE)), 
                                                                                                uiOutput("querydim_hover_box")), width = 12))
  }
  withProgress(message = "Loading bridge and reference", expr = {
    disable(id = "file")
    ToggleDemos(action = "disable", demos = demos)
    setProgress(value = 0.2)
    refs <- LoadBridgeReference(path = getOption(x = "Azimuth.app.reference", 
                                           default = stop(safeError(error = "No reference provided"))))
    setProgress(value = 1)
    enable(id = "file")
    ToggleDemos(action = "enable", demos = demos)
  })
  if (!is.null(x = googlesheet)) {
    try(expr = sheet_append(ss = googlesheet, data = data.frame(c("REFERENCE_NAME", 
                                                                  "REFERENCE_VERSION"), c(app_session_id, app_session_id), 
                                                                c(basename(getOption(x = "Azimuth.app.reference")), 
                                                                  ReferenceVersion(object = refs$map)))), silent = TRUE)
  }
  plotseed <- getOption(x = "Azimuth.app.plotseed")
  if (!is.null(x = plotseed)) {
    set.seed(seed = plotseed)
    colormap <- GetColorMap(object = refs$map)
    for (i in names(x = colormap)) {
      names(x = colormap[[i]]) <- sample(x = names(x = colormap[[i]]))
    }
    refs$map <- SetColorMap(object = refs$map, value = colormap)
  }
  metadata.annotate <- names(x = GetColorMap(object = refs$map))
  if (!is.null(x = getOption(x = "Azimuth.app.metadata_notransfer", 
                             default = NULL))) {
    metadata.notransfer <- str_trim(str_split(getOption(x = "Azimuth.app.metadata_notransfer", 
                                                        default = NULL), ",")[[1]])
    possible.metadata.transfer <- setdiff(x = metadata.annotate, 
                                          y = metadata.notransfer)
  }
  else {
    possible.metadata.transfer <- metadata.annotate
  }
  if (length(x = possible.metadata.transfer) > 1) {
    react.env$xferopts <- TRUE
  }
  default_xfer <- getOption(x = "Azimuth.app.default_metadata", 
                            default = possible.metadata.transfer[1])
  if (!default_xfer %in% possible.metadata.transfer) {
    default_xfer <- possible.metadata.transfer[1]
  }
  observeEvent(eventExpr = input$file, handlerExpr = {
    ResetEnv()
    if (nchar(x = input$file$datapath)) {
      react.env$path <- input$file$datapath
    }
  })
  observeEvent(eventExpr = sapply(X = app.env$demo.inputs, 
                                  FUN = function(x) input[[x]]), handlerExpr = {
                                    if (isTRUE(x = !all(sapply(X = app.env$demo.inputs, FUN = is.null)))) {
                                      ResetEnv()
                                      for (i in 1:length(x = app.env$demo.inputs)) {
                                        if (isTRUE(x = input[[app.env$demo.inputs[i]]] != 
                                                   app.env$demo.tracker[i])) {
                                          app.env$demo.tracker[i] <- app.env$demo.tracker[i] + 
                                            1
                                          react.env$path <- app.env$demo.files[i]
                                        }
                                      }
                                    }
                                  }, ignoreInit = TRUE)
  observeEvent(eventExpr = react.env$path, handlerExpr = {
    if (!is.null(x = react.env$path) && nchar(x = react.env$path)) {
      withProgress(message = "Reading ATAC Peaks", expr = {
        setProgress(value = 0)
        tryCatch(expr = {
          app.env$counts <- LoadFileInput(path = react.env$path)
          if (react.env$path %in% app.env$demo.files) {
            app.env$demo <- TRUE
          }
          else {
            app.env$demo <- FALSE
          }
          app.env$counts$query <- "query"
          react.env$path <- NULL
          react.env$chromatin_assay_1 <- TRUE
        }, error = function(e) {
          app.env$messages <- e$message
          showNotification(e$message, duration = 10, 
                           type = "error", closeButton = TRUE, id = "no-progress-notification")
          app.env$counts <- NULL
          gc(verbose = FALSE)
          react.env$path <- NULL
        })
        setProgress(value = 0.2)
      })
    }
  })
  # observeEvent(eventExpr = react.env$annotations, handlerExpr = { # LOADING THIS IN NOW 
  #   if (isTRUE(x = react.env$annotations)){
  #     withProgress(message = "Loading Annotations", expr = {
  #       setProgress(value = 0.1)
  #       tryCatch(expr = {
  #         app.env$annotations <- suppressWarnings(GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)) # this could also be saved as an object if we wanna save 30 s 
  #         print("getting seq level styles")
  #         seqlevelsStyle(app.env$annotations) <- 'UCSC'
  #         print("got seq level styles")
  #         react.env$annotations <- NULL
  #         react.env$chromatin_assay_1 <- TRUE
  #       }, error = function(e) {
  #         app.env$messages <- e$message
  #         showNotification(e$message, duration = 10, 
  #                          type = "error", closeButton = TRUE, id = "no-progress-notification")
  #         app.env$annotations <- NULL
  #         gc(verbose = FALSE)
  #         react.env$annotations <- NULL
  #       }
  #       )
  #       setProgress(value = 0.15)
  #     })
  #   }
  # })
  observeEvent(eventExpr = react.env$chromatin_assay_1, handlerExpr = {
    if (isTRUE(x = react.env$chromatin_assay_1)) {
      withProgress(message = "Making Chromatin Assay", expr = {
        setProgress(value = 0.2)
        tryCatch(expr = {
          print("about to make chromatin assay")
          app.env$annotations <- refs$annotation
          app.env$chromatin_assay_1 <- CreateChromatinAssay(
            counts = app.env$counts[["RNA"]]@counts, # this should probably be clearer 
            sep = c(":", "-"),
            annotation = app.env$annotations
          )
          print(app.env$chromatin_assay_1)
          print(refs$bridge)
          #app.env$o_hits <- findOverlaps(app.env$chromatin_assay_1, refs$bridge[["ATAC"]])
          #qc_table <- OverlapQC(app.env$chromatin_assay_1, refs$bridge)
          perc_overlap <- round(x = OverlapTotal(app.env$chromatin_assay_1, refs$bridge[["ATAC"]]), digits = 4)
          print("PERC OVERLAP")
          print(perc_overlap)
          if (perc_overlap >= 80) {
            output$valuebox.overlap <- renderValueBox(expr = {
              valueBox(value = perc_overlap, subtitle = "Overlap Percentage",
                       icon = icon(name = "check"), color = "green")
            })
          }
          else if (perc_overlap < 80 & perc_overlap > 60) {
            output$valuebox.overlap<- renderValueBox(expr = {
              valueBox(value = perc_overlap, subtitle = "Overlap Percentage",
                       icon = icon(name = "exclamation-circle"), color = "yellow")
            })
            
          }
          else {
            output$valuebox.overlap <- renderValueBox(expr = {
              valueBox(value = perc_overlap, subtitle = "Overlap Percentage Too Low",
                       icon = icon(name = "exclamation-circle"), color = "red")
            })
          }
          jaccard <- round(x = PeakJaccard(app.env$chromatin_assay_1, refs$bridge[["ATAC"]]), digits = 4)
          print("JACCARD SIMILARITY")
          print(jaccard)
          if (jaccard >= 50) {
            output$valuebox.jaccard <- renderValueBox(expr = {
              valueBox(value = jaccard, subtitle = "Jaccard Similarity",
                       icon = icon(name = "check"), color = "green")
            })
          }
          else if (perc_overlap < 50 & perc_overlap > 20) {
            output$valuebox.jaccard<- renderValueBox(expr = {
              valueBox(value = jaccard, subtitle = "Jaccard Similarity",
                       icon = icon(name = "exclamation-circle"), color = "yellow")
            })
            
          }
          else {
            output$valuebox.jaccard <- renderValueBox(expr = {
              valueBox(value = jaccard, subtitle = "Jaccard Similarity is Low",
                       icon = icon(name = "exclamation-circle"), color = "red")
            })
          }
          print("made chromatin assay ")
          query.cell.names <- paste0("query", 1:ncol(x = app.env$chromatin_assay_1))
          print("got query cell names")
          head(query.cell.names)
          while (any(query.cell.names %in% Cells(x = refs$map))) {
            query.cell.names <- paste0(query.cell.names, 
                                       "x")
          }
          print("assessed query cell names")
          app.env$query.names <- Cells(x = app.env$chromatin_assay_1)
          print("QUERY CELLS")
          print(length(app.env$query.names))
          print("got cells ")
          app.env$chromatin_assay_1 <- RenameCells(object = app.env$chromatin_assay_1, 
                                                   new.names = query.cell.names)
          print("renamed cells")
          
          # remove this because we don't have mitochondrial genes, just peaks 
          removeUI(selector = '#pctmt', immediate = TRUE)
          react.env$mt <- FALSE
          
          
          react.env$requantify_multiome <- TRUE
          react.env$chromatin_assay_1 <- FALSE
        }, error = function(e) {
          app.env$messages <- e$message
          showNotification(e$message, duration = 10, 
                           type = "error", closeButton = TRUE, id = "no-progress-notification")
          app.env$chromatin_assay_1 <- NULL
          gc(verbose = FALSE)
          react.env$chromatin_assay_1 <- NULL
        })
        setProgress(value = 0.3)
      })
    }
  })
  observeEvent(eventExpr = react.env$requantify_multiome, handlerExpr = {
    if (isTRUE(x = react.env$requantify_multiome)) {
      withProgress(message = "Requantifying Peaks to Match Bridge", expr = {
        setProgress(value = 0.6)
        tryCatch(expr = {
          app.env$requantified_multiome <- RequantifyPeaks(app.env$chromatin_assay_1, refs$bridge)
          app.env$chromatin_assay_2 <- CreateChromatinAssay(
            counts = app.env$requantified_multiome,
            sep = c(":", "-"),
            annotation = app.env$annotations
          )
          app.env$object <- CreateSeuratObject(counts = app.env$chromatin_assay_2, assay = 'ATAC')
          app.env$object[['peak.orig']] <- app.env$chromatin_assay_1
          app.env$object$query <- "query"
          app.env$default.assay <- DefaultAssay(app.env$object)
          
          common.features <- intersect(
            x = rownames(x = app.env$object),
            y = rownames(x = refs$map)
          )
          reject <- c(
            length(x = common.features) < getOption(x = 'Azimuth.map.ngenes'),
            length(x = Cells(x = app.env$object)) > getOption(x = 'Azimuth.app.max_cells')
          )
          if (any(reject)) {
            app.env$object <- NULL
            gc(verbose = FALSE)
            reject <- min(which(x = reject))
            app.env$messages <- paste(
              c(
                'Not enough genes in common with reference.',
                'Too many cells.'
              ),
              'Try another dataset.'
            )[reject]
          }
          if (isFALSE(x = react.env$xferopts)) {
            removeUI(selector = '#xferopts', immediate = TRUE)
          }
          
          
          react.env$qc <- !any(reject)
          react.env$requantify_multiome <- FALSE
        }, error = function(e) {
          app.env$messages <- e$message
          showNotification(e$message, duration = 10, 
                           type = "error", closeButton = TRUE, id = "no-progress-notification")
          app.env$chromatin_assay_2 <- NULL
          gc(verbose = FALSE)
          react.env$requantify_multiome <- NULL
        })
        setProgress(value = 1)
      })
    }
  })
  
  observeEvent(eventExpr = react.env$qc, handlerExpr = {
    if (isTRUE(x = react.env$qc)) {
      for (id in qc.ids) {
        try(expr = enable(id = id), silent = TRUE)
      }
      ncount <- paste0("nCount_", app.env$default.assay)
      print(ncount)
      nfeature <- paste0("nFeature_", app.env$default.assay)
      print(nfeature)
      if (!all(c(ncount, nfeature) %in% colnames(x = app.env$object[[]]))) {
        withProgress(message = "Calculating nCount and nFeature", 
                     expr = {
                       setProgress(value = 0)
                       calcn <- as.data.frame(x = Seurat:::CalcN(object = app.env$object))
                       colnames(x = calcn) <- paste(colnames(x = calcn), 
                                                    app.env$default.assay, sep = "_")
                       app.env$object <- AddMetaData(object = app.env$object, 
                                                     metadata = calcn)
                       rm(calcn)
                       gc(verbose = FALSE)
                       setProgress(value = 1)
                     })
      }
      ncount.val <- range(app.env$object[[ncount, drop = TRUE]])
      ncount.val <- c(floor(x = min(ncount.val)), ceiling(x = max(ncount.val)))
      ncount.min <- if (is.null(getOption(x = "Azimuth.app.ncount_min"))) {
        ncount.val[1]
      }
      else {
        max(ncount.val[1], getOption(x = "Azimuth.app.ncount_min"))
      }
      ncount.max <- if (is.null(getOption(x = "Azimuth.app.ncount_max"))) {
        ncount.val[2]
      }
      else {
        min(ncount.val[2], getOption(x = "Azimuth.app.ncount_max"))
      }
      print(ncount.min)
      print(ncount.max)
      updateNumericInput(session = session, inputId = "num.ncountmin", 
                         label = paste("min", ncount), value = ncount.min, 
                         min = ncount.val[1], max = ncount.val[2])
      updateNumericInput(session = session, inputId = "num.ncountmax", 
                         label = paste("max", ncount), value = ncount.max, 
                         min = ncount.val[1], max = ncount.val[2])
      nfeature.val <- range(app.env$object[[nfeature, drop = TRUE]])
      nfeature.val <- c(floor(x = min(nfeature.val)), ceiling(x = max(nfeature.val)))
      nfeature.min <- if (is.null(getOption(x = "Azimuth.app.nfeature_min"))) {
        nfeature.val[1]
      }
      else {
        max(nfeature.val[1], getOption(x = "Azimuth.app.nfeature_min"))
      }
      nfeature.max <- if (is.null(getOption(x = "Azimuth.app.nfeature_max"))) {
        nfeature.val[2]
      }
      else {
        min(nfeature.val[2], getOption(x = "Azimuth.app.nfeature_max"))
      }
      print(nfeature.min)
      print(nfeature.max)
      updateNumericInput(session = session, inputId = "num.nfeaturemin", 
                         label = paste("min", nfeature), value = nfeature.min, 
                         min = nfeature.val[1], max = nfeature.val[2])
      updateNumericInput(session = session, inputId = "num.nfeaturemax", 
                         label = paste("max", nfeature), value = nfeature.max, 
                         min = nfeature.val[1], max = nfeature.val[2])
      if (isTRUE(x = react.env$mt)) {
        app.env$object <- PercentageFeatureSet(object = app.env$object, 
                                               pattern = mito.pattern, col.name = mt.key, 
                                               assay = app.env$default.assay)
        mito.val <- range(app.env$object[[mt.key, drop = TRUE]])
        mito.val <- c(floor(x = min(mito.val)), ceiling(x = max(mito.val)))
        mito.min <- if (is.null(getOption(x = "Azimuth.app.pctmt_min"))) {
          mito.val[1]
        }
        else {
          max(mito.val[1], getOption(x = "Azimuth.app.pctmt_min"))
        }
        mito.max <- if (is.null(getOption(x = "Azimuth.app.pctmt_max"))) {
          mito.val[2]
        }
        else {
          min(mito.val[2], getOption(x = "Azimuth.app.pctmt_max"))
        }
        updateNumericInput(session = session, inputId = "minmt", 
                           label = paste("min", mt.key), value = mito.min, 
                           min = mito.val[1], max = mito.val[2])
        updateNumericInput(session = session, inputId = "maxmt", 
                           label = paste("max", mt.key), value = mito.max, 
                           min = mito.val[1], max = mito.val[2])
      }
      print("about to render menu tab")
      output$menu1 <- renderMenu(expr = {
        sidebarMenu(menuItem(text = "Preprocessing", 
                             tabName = "tab_preproc", icon = icon(name = "filter"), 
                             selected = TRUE))
      })
      ncellsupload <- length(x = colnames(x = app.env$object))
      app.env$ncellsupload <- ncellsupload
      app.env$messages <- paste(ncellsupload, "cells uploaded")
      print("about to make value box")
      if (ncellsupload < getOption(x = "Azimuth.map.ncells")) {
        output$valuebox.upload <- renderValueBox(expr = {
          valueBox(value = ncellsupload, subtitle = paste0("cells uploaded - ", 
                                                           getOption(x = "Azimuth.map.ncells"), " required"), 
                   icon = icon(name = "times"), color = "red")
        })
      }
      else {
        output$valuebox.upload <- renderValueBox(expr = {
          valueBox(value = ncellsupload, subtitle = "cells uploaded", 
                   icon = icon(name = "check"), color = "green")
        })
        if (!is.null(x = googlesheet)) {
          try(expr = sheet_append(ss = googlesheet, data = data.frame("CELLSUPLOAD", 
                                                                      app_session_id, ncellsupload)), silent = TRUE)
        }
      }
      if (!is.null(x = react.env$progress)) {
        react.env$progress$close()
        enable(id = "file")
        ToggleDemos(action = "enable", demos = demos)
        react.env$progress <- NULL
      }
      updateSelectizeInput(session = session, inputId = "metadataxfer", 
                           choices = possible.metadata.transfer, selected = default_xfer, 
                           server = TRUE, options = selectize.opts[-which(x = names(x = selectize.opts) == 
                                                                            "maxItems")])
      react.env$qc <- FALSE
      react.env$plot.qc <- TRUE
      react.env$dist.qc <- TRUE
      print("about to make qc plots")
    }
  })
  observeEvent(eventExpr = input$metadataxfer, handlerExpr = {
    if (length(x = input$metadataxfer) == 0) {
      disable(id = "map")
    }
    else {
      enable(id = "map")
    }
  }, ignoreNULL = FALSE)
  observeEvent(eventExpr = input$map, handlerExpr = {
    react.env$start <- Sys.time()
    disable(id = "file")
    ToggleDemos(action = "disable", demos = demos)
    for (id in qc.ids) {
      try(expr = disable(id = id), silent = TRUE)
    }
    react.env$progress <- Progress$new(style = "notification")
    react.env$progress$set(value = 0, message = "Filtering based on nCount and nFeature")
    ncount <- paste0("nCount_", DefaultAssay(object = app.env$object))
    nfeature <- paste0("nFeature_", DefaultAssay(object = app.env$object))
    cells.use <- app.env$object[[ncount, drop = TRUE]] >= 
      input$num.ncountmin & app.env$object[[ncount, drop = TRUE]] <= 
      input$num.ncountmax & app.env$object[[nfeature, drop = TRUE]] >= 
      input$num.nfeaturemin & app.env$object[[nfeature, 
                                              drop = TRUE]] <= input$num.nfeaturemax
    # if (isTRUE(x = react.env$mt)) { # deleting this bc i think its leaving the mt bar thing 
    #   cells.use <- cells.use & app.env$object[[mt.key, 
    #                                            drop = TRUE]] >= input$minmt & app.env$object[[mt.key, 
    #                                                                                           drop = TRUE]] <= input$maxmt
    # }
    ncellspreproc <- sum(cells.use)
    app.env$ncellspreproc <- ncellspreproc
    if (ncellspreproc < getOption(x = "Azimuth.map.ncells")) {
      output$valuebox.preproc <- renderValueBox(expr = valueBox(value = ncellspreproc, 
                                                                subtitle = paste0("cells after filtering - ", 
                                                                                  getOption(x = "Azimuth.map.ncells"), " required"), 
                                                                icon = icon("times"), color = "red"))
      react.env$qc <- TRUE
    }
    else {
      output$valuebox.preproc <- renderValueBox(expr = valueBox(value = ncellspreproc, 
                                                                subtitle = "cells after filtering", icon = icon("check"), 
                                                                color = "green"))
      if (!is.null(googlesheet)) {
        try(sheet_append(ss = googlesheet, data = data.frame("CELLSPREPROC", 
                                                             app_session_id, ncellspreproc)))
      }
      app.env$object <- app.env$object[, cells.use]
      print("LENGTH OF CELLS USE: length(cells.use")
      print(length(cells.use))
      print(app.env$object)
      app.env$query.names <- app.env$query.names[cells.use]
      #react.env$sctransform <- TRUE replacing
      react.env$tfidf <- TRUE
    }
  })
  observeEvent(eventExpr = react.env$tfidf, handlerExpr = {
    if (isTRUE(x = react.env$tfidf)) {
      react.env$progress$set(value = 0.2, message = "Normalizing with TFIDF")
      tryCatch(expr = {
        app.env$object <- suppressWarnings(expr = RunTFIDF(object = app.env$object,
                                                           method = 1))
      }, error = function(e) {
        app.env$object <- suppressWarnings(expr = RunTFIDF(object = app.env$object,
                                                           method = 1))
      })
      app.env$messages <- c(app.env$messages, paste(ncol(x = app.env$object), 
                                                    "cells preprocessed"))
      #react.env$anchors <- TRUE
      react.env$bridge_anchors <- TRUE
      #react.env$sctransform <- FALSE
      react.env$tfidf <- FALSE
    }
  })
  observeEvent(eventExpr = react.env$sctransform, handlerExpr = {
    if (isTRUE(x = react.env$sctransform)) {
      react.env$progress$set(value = 0.2, message = "Normalizing with SCTransform")
      tryCatch(expr = {
        app.env$object <- suppressWarnings(expr = SCTransform(object = app.env$object, 
                                                              residual.features = rownames(x = refs$map), 
                                                              reference.SCT.model = slot(object = refs$map[["refAssay"]], 
                                                                                         name = "SCTModel.list")[["refmodel"]], method = "glmGamPoi", 
                                                              do.correct.umi = FALSE, do.scale = FALSE, do.center = TRUE, 
                                                              new.assay.name = "refAssay"))
      }, error = function(e) {
        app.env$object <- suppressWarnings(expr = SCTransform(object = app.env$object, 
                                                              residual.features = rownames(x = refs$map), 
                                                              reference.SCT.model = slot(object = refs$map[["refAssay"]], 
                                                                                         name = "SCTModel.list")[["refmodel"]], method = "poisson", 
                                                              do.correct.umi = FALSE, do.scale = FALSE, do.center = TRUE, 
                                                              new.assay.name = "refAssay"))
      })
      app.env$messages <- c(app.env$messages, paste(ncol(x = app.env$object), 
                                                    "cells preprocessed"))
      react.env$anchors <- TRUE
      react.env$sctransform <- FALSE
    }
  })
  observeEvent(eventExpr = react.env$anchors, handlerExpr = {
    if (isTRUE(x = react.env$anchors)) {
      react.env$progress$set(value = 0.3, message = "Finding anchors")
      app.env$anchors <- FindTransferAnchors(reference = refs$map, 
                                             query = app.env$object, k.filter = NA, reference.neighbors = "refdr.annoy.neighbors", 
                                             reference.assay = "refAssay", query.assay = "refAssay", 
                                             reference.reduction = "refDR", normalization.method = "SCT", 
                                             features = intersect(x = rownames(x = refs$map), 
                                                                  y = VariableFeatures(object = app.env$object)), 
                                             dims = 1:getOption(x = "Azimuth.map.ndims"), 
                                             n.trees = n.trees, verbose = TRUE, mapping.score.k = 100)
      nanchors <- nrow(x = slot(object = app.env$anchors, 
                                name = "anchors"))
      app.env$nanchors <- nanchors
      if (!is.null(googlesheet)) {
        try(sheet_append(ss = googlesheet, data = data.frame("NANCHORS", 
                                                             app_session_id, nanchors)))
      }
      if (nanchors < getOption(x = "Azimuth.map.nanchors") | 
          length(x = unique(x = slot(object = app.env$anchors, 
                                     name = "anchors")[, 2])) < 50) {
        output$valuebox.mapped <- renderValueBox(expr = {
          valueBox(value = "Failure", subtitle = paste0("Too few anchors identified (", 
                                                        nanchors, ")"), icon = icon(name = "times"), 
                   color = "red", width = 6)
        })
        app.env$object <- NULL
        app.env$anchors <- NULL
        react.env$progress$close()
        enable(id = "file")
        ToggleDemos(action = "enable", demos = demos)
        gc(verbose = FALSE)
      }
      else {
        query.unique <- length(x = unique(x = slot(object = app.env$anchors, 
                                                   name = "anchors")[, "cell2"]))
        percent.anchors <- round(x = query.unique/ncol(x = app.env$object) * 
                                   100, digits = 2)
        if (percent.anchors < getOption(x = "Azimuth.map.panchorscolors")[1]) {
          output$valuebox_panchors <- renderValueBox(expr = {
            valueBox(value = paste0(percent.anchors, 
                                    "%"), subtitle = "% of query cells with anchors", 
                     color = "red", icon = icon(name = "times"))
          })
        }
        else if (percent.anchors < getOption(x = "Azimuth.map.panchorscolors")[2]) {
          output$valuebox_panchors <- renderValueBox(expr = {
            valueBox(value = paste0(percent.anchors, 
                                    "%"), subtitle = "% of query cells with anchors", 
                     color = "yellow", icon = icon(name = "exclamation-circle"))
          })
        }
        else {
          output$valuebox_panchors <- renderValueBox(expr = {
            valueBox(value = paste0(percent.anchors, 
                                    "%"), subtitle = "% of query cells with anchors", 
                     color = "green", icon = icon(name = "check"))
          })
        }
        react.env$map <- TRUE
      }
      react.env$anchors <- FALSE
    }
  })
  observeEvent(eventExpr = react.env$bridge_anchors, handlerExpr = {
    if (isTRUE(x = react.env$bridge_anchors)) {
      react.env$progress$set(value = 0.3, message = "Finding anchors")
      app.env$bridge_anchors <- FindBridgeTransferAnchors(extended.reference = refs$ext,
                                                          query = app.env$object,
                                                          reduction = "lsiproject",
                                                          scale = FALSE,
                                                          dims = 2:50) # making this a default
      print("found anchors ")
      nanchors <- nrow(x = slot(object = app.env$bridge_anchors, 
                                name = "anchors"))
      print(nanchors)
      app.env$nanchors <- nanchors
      if (!is.null(googlesheet)) {
        print("doing google sheet thing")
        try(sheet_append(ss = googlesheet, data = data.frame("NANCHORS", 
                                                             app_session_id, nanchors)))
      }
      if (nanchors < getOption(x = "Azimuth.map.nanchors") | 
          length(x = unique(x = slot(object = app.env$bridge_anchors, 
                                     name = "anchors")[, 2])) < 50) {
        print("doing this length if ")
        output$valuebox.mapped <- renderValueBox(expr = {
          valueBox(value = "Failure", subtitle = paste0("Too few anchors identified (", 
                                                        nanchors, ")"), icon = icon(name = "times"), 
                   color = "red", width = 6)
        })
        app.env$object <- NULL
        app.env$bridge_anchors <- NULL
        react.env$progress$close()
        enable(id = "file")
        ToggleDemos(action = "enable", demos = demos)
        
        gc(verbose = FALSE)
      }
      else {
        print("doing the else")
        query.unique <- length(x = unique(x = slot(object = app.env$bridge_anchors, 
                                                   name = "anchors")[, "cell2"]))
        percent.anchors <- round(x = query.unique/ncol(x = app.env$object) * 
                                   100, digits = 2)
        if (percent.anchors < getOption(x = "Azimuth.map.panchorscolors")[1]) {
          output$valuebox_panchors <- renderValueBox(expr = {
            valueBox(value = paste0(percent.anchors, 
                                    "%"), subtitle = "% of query cells with anchors", 
                     color = "red", icon = icon(name = "times"))
          })
        }
        else if (percent.anchors < getOption(x = "Azimuth.map.panchorscolors")[2]) {
          print("doing the else if ")
          output$valuebox_panchors <- renderValueBox(expr = {
            valueBox(value = paste0(percent.anchors, 
                                    "%"), subtitle = "% of query cells with anchors", 
                     color = "yellow", icon = icon(name = "exclamation-circle"))
          })
        }
        else {
          print("doing the second else ")
          output$valuebox_panchors <- renderValueBox(expr = {
            valueBox(value = paste0(percent.anchors, 
                                    "%"), subtitle = "% of query cells with anchors", 
                     color = "green", icon = icon(name = "check"))
          })
        }
        print("finished anchors")
        react.env$mapquery <- TRUE
      }
      react.env$bridge_anchors <- FALSE
    }
  })
  observeEvent(eventExpr = list(react.env$map, input$metadataxfer), 
               handlerExpr = {
                 if (isTRUE(x = react.env$map)) {
                   if (is.null(x = input$metadataxfer)) {
                     app.env$metadataxfer <- names(x = GetColorMap(object = refs$map))
                   }
                   else {
                     app.env$metadataxfer <- input$metadataxfer
                   }
                   react.env$progress$set(value = 0.5, message = "Mapping cells")
                   refdata <- lapply(X = app.env$metadataxfer, function(x) {
                     refs$map[[x, drop = TRUE]]
                   })
                   names(x = refdata) <- app.env$metadataxfer
                   if (do.adt) {
                     refdata[["impADT"]] <- GetAssayData(object = refs$map[["ADT"]], 
                                                         slot = "data")
                   }
                   app.env$object <- TransferData(reference = refs$map, 
                                                  query = app.env$object, dims = 1:getOption(x = "Azimuth.map.ndims"), 
                                                  anchorset = app.env$anchors, refdata = refdata, 
                                                  n.trees = n.trees, store.weights = TRUE)
                   app.env$singlepred <- NULL
                   for (i in app.env$metadataxfer) {
                     app.env$singlepred <- c(app.env$singlepred, 
                                             length(x = unique(x = as.vector(x = app.env$object[[paste0("predicted.", 
                                                                                                        i), drop = TRUE]]))) == 1)
                     app.env$object[[paste0("predicted.", i), drop = TRUE]] <- factor(x = app.env$object[[paste0("predicted.", 
                                                                                                                 i), drop = TRUE]], levels = levels(x = refs$map[[i, 
                                                                                                                                                                  drop = TRUE]]))
                   }
                   singlepred <- all(app.env$singlepred)
                   if (singlepred & (length(x = setdiff(possible.metadata.transfer, 
                                                        app.env$metadataxfer)) > 0)) {
                     showNotification(paste0("Only one predicted class. Re-running with all metadata."), 
                                      duration = 5, type = "warning", closeButton = TRUE, 
                                      id = "no-progress-notification")
                     updateSelectizeInput(session = getDefaultReactiveDomain(), 
                                          inputId = "metadataxfer", choices = possible.metadata.transfer, 
                                          selected = possible.metadata.transfer, )
                     app.env$metadataxfer <- input$metadataxfer
                   }
                   else if (singlepred) {
                     showNotification(paste0("Only one predicted class: ", 
                                             app.env$object[[paste0("predicted.", app.env$metadataxfer[1]), 
                                                             drop = TRUE]][1]), duration = 5, type = "warning", 
                                      closeButton = TRUE, id = "no-progress-notification")
                     app.env$object <- NULL
                     app.env$anchors <- NULL
                     react.env$path <- NULL
                     react.env$map <- FALSE
                     react.env$progress$close()
                     enable(id = "file")
                     ToggleDemos(action = "enable", demos = demos)
                     gc(verbose = FALSE)
                   }
                   else {
                     app.env$object <- IntegrateEmbeddings(anchorset = app.env$anchors, 
                                                           reference = refs$map, query = app.env$object, 
                                                           reductions = "pcaproject", reuse.weights.matrix = TRUE)
                     if (is.null(x = getOption(x = "Azimuth.app.default_metadata"))) {
                       app.env$default.metadata <- names(x = refdata)[1]
                     }
                     else {
                       if (getOption(x = "Azimuth.app.default_metadata") %in% 
                           names(x = refdata)) {
                         app.env$default.metadata <- getOption(x = "Azimuth.app.default_metadata")
                       }
                       else {
                         app.env$default.metadata <- names(x = refdata)[1]
                       }
                     }
                     react.env$score <- TRUE
                     react.env$map <- FALSE
                   }
                 }
               })
  observeEvent(eventExpr = list(react.env$mapquery, input$metadataxfer), 
               handlerExpr = {
                 if (isTRUE(x = react.env$mapquery)) {
                   print("doing this mapquery")
                   if (is.null(x = input$metadataxfer)) {
                     app.env$metadataxfer <- names(x = GetColorMap(object = refs$map))
                   }
                   else {
                     app.env$metadataxfer <- input$metadataxfer
                   }
                   print("mapping cells ")
                   react.env$progress$set(value = 0.5, message = "Mapping cells")
                   refdata <- lapply(X = app.env$metadataxfer, function(x) { 
                     refs$map[[x, drop = TRUE]]
                   })
                   names(x = refdata) <- app.env$metadataxfer
                   if (do.adt) {
                     print("trying to do adt")
                     refdata[["impADT"]] <- GetAssayData(object = refs$map[["ADT"]], 
                                                         slot = "data")
                   }
                   app.env$object <-  MapQuery(anchorset = app.env$bridge_anchors,  # deleted transfer data 
                                               reference = refs$map, 
                                               query = app.env$object, 
                                               refdata = refdata,
                                               reduction.model = "refUMAP")
                   print("finished mapquery")
                   print(app.env$metadataxfer)
                   app.env$singlepred <- NULL
                   for (i in app.env$metadataxfer) { 
                     app.env$singlepred <- c(app.env$singlepred, 
                                             length(x = unique(x = as.vector(x = app.env$object[[paste0("predicted.", 
                                                                                                        i), drop = TRUE]]))) == 1)
                     app.env$object[[paste0("predicted.", i), drop = TRUE]] <- factor(x = app.env$object[[paste0("predicted.", 
                                                                                                                 i), drop = TRUE]], levels = levels(x = refs$map[[i, 
                                                                                                                                                                  drop = TRUE]]))
                   }
                   singlepred <- all(app.env$singlepred)
                   if (singlepred & (length(x = setdiff(possible.metadata.transfer, 
                                                        app.env$metadataxfer)) > 0)) {
                     showNotification(paste0("Only one predicted class. Re-running with all metadata."), 
                                      duration = 5, type = "warning", closeButton = TRUE, 
                                      id = "no-progress-notification")
                     updateSelectizeInput(session = getDefaultReactiveDomain(), 
                                          inputId = "metadataxfer", choices = possible.metadata.transfer, 
                                          selected = possible.metadata.transfer, )
                     app.env$metadataxfer <- input$metadataxfer
                   }
                   else if (singlepred) {
                     showNotification(paste0("Only one predicted class: ", 
                                             app.env$object[[paste0("predicted.", app.env$metadataxfer[1]), 
                                                             drop = TRUE]][1]), duration = 5, type = "warning", 
                                      closeButton = TRUE, id = "no-progress-notification")
                     app.env$object <- NULL
                     app.env$bridge_anchors <- NULL
                     react.env$path <- NULL
                     react.env$mapquery <- FALSE
                     react.env$progress$close()
                     enable(id = "file")
                     ToggleDemos(action = "enable", demos = demos)
                     gc(verbose = FALSE)
                   }
                   else {
                     if (is.null(x = getOption(x = "Azimuth.app.default_metadata"))) {
                       app.env$default.metadata <- names(x = refdata)[1]
                     }
                     else {
                       if (getOption(x = "Azimuth.app.default_metadata") %in% 
                           names(x = refdata)) {
                         app.env$default.metadata <- getOption(x = "Azimuth.app.default_metadata")
                       }
                       else {
                         app.env$default.metadata <- names(x = refdata)[1]
                       }
                     }
                     #react.env$score <- TRUE - ill do this after getting gene activity scores 
                     react.env$gene_activity <- TRUE
                     react.env$mapquery <- FALSE
                   }
                 }
               })
  observeEvent(eventExpr = react.env$gene_activity, handlerExpr = {
    if (isTRUE(react.env$gene_activity)) {
      # Use original peaks 
      print("GENE ACTIVITY ")
      DefaultAssay(app.env$object) <- "peak.orig"
      print(app.env$object)
      startTime <- Sys.time()
      app.env$transcripts <- GetTranscripts(app.env$object)
      endTime <- Sys.time()
      print("TOTAL TIME TO GET TRANSCRIPTS")
      print(endTime - startTime)
      temp <- RequantifyPeaks(app.env$object, app.env$transcripts)
      #add feature matrix to Chromatin Assay 
      app.env$object[['RNA']] <- CreateAssayObject(counts = temp)
      
      #o_hits <- findOverlaps(app.env$object[["ATAC"]], app.env$transcripts)
      #temp <- RequantifyPeaks(o_hits, app.env$object, app.env$transcripts)
      #dd feature matrix to Chromatin Assay 
      #app.env$object[['RNA']] <- CreateAssayObject(counts = temp)
      
      #Normalize the feature data
      app.env$object <- NormalizeData(
        object = app.env$object,
        assay = 'RNA',
        normalization.method = 'LogNormalize',
        scale.factor = median(app.env$object$nCount_RNA)
      )
      print("feature data normalized")
      
      react.env$gene_activity <- FALSE
      react.env$score <- TRUE
    }
  })
  observeEvent(eventExpr = react.env$cluster.score, handlerExpr = {
    if (isTRUE(react.env$cluster.score)) {
      qc.stat <- round(x = ClusterPreservationScore(query = app.env$object, 
                                                    ds.amount = getOption(x = "Azimuth.map.postmapqcds"), type = "bridge"), 
                       digits = 2)
      if (!is.null(googlesheet)) {
        try(sheet_append(ss = googlesheet, data = data.frame("CLUSTERPRESERVATIONQC", 
                                                             app_session_id, qc.stat)))
      }
      app.env$clusterpreservationqc <- qc.stat
      if (qc.stat < getOption(x = "Azimuth.map.postmapqccolors")[1]) {
        output$valuebox_mappingqcstat <- renderValueBox(expr = {
          valueBox(value = paste0(qc.stat, "/5"), subtitle = "cluster preservation score", 
                   color = "red", icon = icon(name = "times"))
        })
      }
      else if (qc.stat < getOption(x = "Azimuth.map.postmapqccolors")[2]) {
        output$valuebox_mappingqcstat <- renderValueBox(expr = {
          valueBox(value = paste0(qc.stat, "/5"), subtitle = "cluster preservation score", 
                   color = "yellow", icon = icon(name = "exclamation-circle"))
        })
      }
      else {
        output$valuebox_mappingqcstat <- renderValueBox(expr = {
          valueBox(value = paste0(qc.stat, "/5"), subtitle = "cluster preservation score", 
                   color = "green", icon = icon(name = "check"))
        })
      }
      react.env$cluster.score <- FALSE
      react.env$transform <- TRUE
    }
  })
  observeEvent(eventExpr = react.env$score, handlerExpr = {
    if (isTRUE(x = react.env$score)) {
      react.env$progress$set(value = 0.7, message = "Calculating mapping score")
      print("qc.stat")
      print(app.env$object@reductions)
      app.env$object[['refAssay']] <- app.env$object[['ATAC']]
      DefaultAssay(app.env$object) <- 'refAssay'
      DefaultAssay(app.env$object[["ref.Bridge.reduc"]]) <- 'refAssay'
      app.env$object <- FindTopFeatures(app.env$object,
                                        min.cutoff = "q0")
      qc.stat <- round(x = ClusterPreservationScore(query = app.env$object, 
                                                    ds.amount = getOption(x = "Azimuth.map.postmapqcds"), type = "bridge"), 
                       digits = 2)
      print("got qc.stat")
      if (!is.null(googlesheet)) {
        try(sheet_append(ss = googlesheet, data = data.frame("CLUSTERPRESERVATIONQC", 
                                                             app_session_id, qc.stat)))
      }
      print("clusterpreservationqx")
      app.env$clusterpreservationqc <- qc.stat
      if (qc.stat < getOption(x = "Azimuth.map.postmapqccolors")[1]) {
        output$valuebox_mappingqcstat <- renderValueBox(expr = {
          valueBox(value = paste0(qc.stat, "/5"), subtitle = "cluster preservation score", 
                   color = "red", icon = icon(name = "times"))
        })
      }
      else if (qc.stat < getOption(x = "Azimuth.map.postmapqccolors")[2]) {
        output$valuebox_mappingqcstat <- renderValueBox(expr = {
          valueBox(value = paste0(qc.stat, "/5"), subtitle = "cluster preservation score", 
                   color = "yellow", icon = icon(name = "exclamation-circle"))
        })
      }
      else {
        output$valuebox_mappingqcstat <- renderValueBox(expr = {
          valueBox(value = paste0(qc.stat, "/5"), subtitle = "cluster preservation score", 
                   color = "green", icon = icon(name = "check"))
        })
      }
      print("getting refdr")
      refdr <- subset(x = app.env$bridge_anchors@object.list[[1]][["Bridge.reduc"]], # im gonna try calling this Bridge.Reduc
                      cells = paste0(Cells(x = app.env$object), "_query"))
      refdr <- RenameCells(object = refdr, new.names = Cells(x = app.env$object))
      print('got refdr')
      refdr.ref <- subset(x = app.env$bridge_anchors@object.list[[1]][["Bridge.reduc"]], 
                          cells = paste0(Cells(x = refs$map), "_reference"))
      refdr.ref <- RenameCells(object = refdr.ref, new.names = Cells(x = refs$map))
      print("got refdr.ref")
      if (Sys.getenv("RSTUDIO") == "1") {
        plan("sequential")
      }
      print("about to diet seurat") # have to change all of the anchors to bridge_anchors
      app.env$bridge_anchors@object.list[[1]] <- DietSeurat(object = app.env$bridge_anchors@object.list[[1]])
      app.env$bridge_anchors@object.list[[1]] <- subset(x = app.env$bridge_anchors@object.list[[1]], 
                                                        features = c(rownames(x = app.env$bridge_anchors@object.list[[1]])[1]))
      print("doing app.env bridge anchors rename")
      app.env$bridge_anchors@object.list[[1]] <- RenameCells(object = app.env$bridge_anchors@object.list[[1]], 
                                                             new.names = unname(obj = sapply(X = Cells(x = app.env$bridge_anchors@object.list[[1]]), 
                                                                                             FUN = function(x) {
                                                                                               return(gsub(pattern = "_reference", replacement = "", 
                                                                                                           x = x))
                                                                                             })))
      print("doing 2")
      app.env$bridge_anchors@object.list[[1]] <- RenameCells(object = app.env$bridge_anchors@object.list[[1]], 
                                                             new.names = sapply(X = Cells(x = app.env$bridge_anchors@object.list[[1]]), 
                                                                                FUN = function(x) {
                                                                                  return(gsub(pattern = "_query", replacement = "", 
                                                                                              x = x))
                                                                                }, USE.NAMES = FALSE))
      app.env$bridge_anchors@object.list[[1]]@meta.data <- data.frame()
      app.env$bridge_anchors@object.list[[1]]@active.ident <- factor()
      mapping.score.k <- min(c(50, length(x = unique(x = app.env$bridge_anchors@anchors[, 
                                                                                        1])), length(x = unique(x = app.env$bridge_anchors@anchors[, 
                                                                                                                                                   2]))))
      app.env$mapping.score <- future(expr = {
        MappingScore(anchors = app.env$bridge_anchors@anchors, 
                     combined.object = app.env$bridge_anchors@object.list[[1]], 
                     query.neighbors = slot(object = app.env$bridge_anchors, 
                                            name = "neighbors")[["query.neighbors"]], 
                     query.weights = Tool(object = app.env$object, 
                                          slot = "TransferData")$weights.matrix, query.embeddings = Embeddings(object = refdr), 
                     ref.embeddings = Embeddings(object = refdr.ref), 
                     nn.method = "annoy", n.trees = n.trees, ndim = getOption(x = "Azimuth.map.ndims"), 
                     kanchors = mapping.score.k)
      })
      app.env$object <- AddMetaData(object = app.env$object, 
                                    metadata = rep(x = 0, times = ncol(x = app.env$object)), 
                                    col.name = "mapping.score")
      app.env$bridge_anchors <- NULL
      rm(refdr)
      gc(verbose = FALSE)
      react.env$cluster.score <- TRUE
      react.env$score <- FALSE
    }
  })
  observeEvent(eventExpr = react.env$transform, handlerExpr = {
    if (isTRUE(x = react.env$transform)) {
      react.env$progress$set(value = 0.8)
      #app.env$object[["query_ref.nn"]] <- FindNeighbors(object = Embeddings(refs$map[["refDR"]])[,
      #  1:getOption("Azimuth.map.ndims")], query = Embeddings(app.env$object[["integrated_dr"]]),
      #return.neighbor = TRUE, l2.norm = TRUE, n.trees = n.trees)
      #app.env$object <- NNTransform(object = app.env$object,
      #               meta.data = refs$map[[]])
      suppressWarnings(expr = app.env$object[["umap.proj"]] <- app.env$object[["ref.umap"]])
      # <- RunUMAP(object = app.env$object[["query_ref.nn"]],
      # reduction.model = refs$map[["refUMAP"]], reduction.key = "UMAP_")
      # app.env$object <- SetAssayData(object = app.env$object,
      #                                assay = "refAssay", slot = "scale.data", new.data = new(Class = "matrix"))
      gc(verbose = FALSE)
      app.env$messages <- c(app.env$messages, paste(ncol(x = app.env$object),
                                                    "cells mapped"))
      react.env$biomarkers <- TRUE
      react.env$chromvar <- TRUE
      react.env$transform <- FALSE
    }
  })
  observeEvent(eventExpr = react.env$chromvar, handlerExpr = {
    if (isTRUE(x = react.env$chromvar)) {
      react.env$progress$set(value = 0.98, message = "Running Motif Analysis")
      DefaultAssay(app.env$object) <- "peak.orig"
      # Remove peaks on scaffolds 
      main.chroms <- GenomeInfoDb:::standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
      keep.peaks <- which(as.character(seqnames(granges(app.env$object))) %in% main.chroms)
      app.env$object[["peak.orig"]] <- subset(app.env$object[["peak.orig"]], features = rownames(app.env$object[["peak.orig"]])[keep.peaks])
      
      pfm <- getMatrixSet(
        x = JASPAR2020,
        opts = list(species = 9606, all_versions = FALSE)
      )
      print("adding motifs")
      startTime <- Sys.time()
      app.env$object <- AddMotifs(
        object = app.env$object,
        genome = BSgenome.Hsapiens.UCSC.hg38,
        pfm = pfm
      )
      endTime <- Sys.time()
      print("TOTAL TIME TO ADD MOTIFS")
      print(endTime - startTime)
      print("calculating chromvar")
      library(BiocParallel)
      register(MulticoreParam(3))
      startTime <- Sys.time()
      app.env$object <- RunChromVAR(
        object = app.env$object,
        genome = BSgenome.Hsapiens.UCSC.hg38
      )
      endTime <- Sys.time()
      print("TOTAL TIME TO RUN CHROMVAR")
      print(endTime - startTime)
      # Rename motifs from ids
      motif_name <- ConvertMotifID(app.env$object[["peak.orig"]]@motifs, id = rownames(app.env$object[["chromvar"]]@data))
      rownames(app.env$object[["chromvar"]]@data) <- motif_name
      
      print(head(row.names(app.env$object[["chromvar"]]@data)))
      for (i in app.env$metadataxfer[!app.env$singlepred]) {
        print("setting chromvar.diff.expr")
        Idents(app.env$object) <- paste0("predicted.", i)
        app.env$chromvar.assay <- "chromvar"
        app.env$chromvar.diff.expr[[paste(app.env$chromvar.assay, # changed all of these to chromvar.assay
                                       i, sep = "_")]] <- FindAllMarkers(object = app.env$object, assay = app.env$chromvar.assay, slot = "data", 
                                                                         only.pos = T, mean.fcn = rowMeans, fc.name = "avg_diff")
        motif_ids <- ConvertMotifID(app.env$object[["peak.orig"]]@motifs, name = app.env$chromvar.diff.expr[[paste(app.env$chromvar.assay, i, sep = "_")]]$gene)
        
        print("MOTIF IDS")
        print(head(motif_ids))
        app.env$chromvar.diff.expr[[paste(app.env$chromvar.assay, i, sep = "_")]]$motif_id <- motif_ids
        print("column names of differential expression")
        print(colnames(app.env$chromvar.diff.expr[[paste(app.env$chromvar.assay, i, sep = "_")]]))
        
      }
      print(head(app.env$chromvar.diff.expr))
      print("about to close progress")
      react.env$progress$close()
      react.env$chromvar <- FALSE
    }
  })
  observeEvent(eventExpr = react.env$biomarkers, handlerExpr = {
    if (isTRUE(x = react.env$biomarkers)) {
      react.env$progress$set(value = 0.95, message = "Running differential expression")
      for (i in app.env$metadataxfer[!app.env$singlepred]) {
        app.env$gene.assay <- "RNA"
        print("setting diff.expr")
        app.env$diff.expr[[paste(app.env$gene.assay, # changed all of these to gene.assay
                                 i, sep = "_")]] <- wilcoxauc(X = app.env$object, 
                                                              group_by = paste0("predicted.", i), assay = "data", 
                                                              seurat_assay = app.env$gene.assay)
        print("DIFF EXPRESSION")
        print(head(app.env$diff.expr))
        if (isTRUE(x = do.adt)) {
          app.env$diff.expr[[paste(adt.key, i, sep = "_")]] <- wilcoxauc(X = app.env$object, 
                                                                         group_by = paste0("predicted.", i), assay = "data", 
                                                                         seurat_assay = adt.key)
        }
      }
      print("DIFF EXPRESSION")
      print(head(app.env$diff.expr))
      # print("calculating mapping time ")
      # mapping.time <- difftime(time1 = Sys.time(), time2 = react.env$start, 
      #                          units = "secs")
      # time.fmt <- FormatDiffTime(dt = mapping.time)
      #app.env$messages <- c(app.env$messages, time.fmt)
      if (!is.null(x = googlesheet)) {
        try(expr = sheet_append(ss = googlesheet, data = data.frame("MAPPINGTIME", 
                                                                    app_session_id, as.numeric(x = mapping.time))))
      }
      if (!is.null(x = googlesheet)) {
        try(expr = sheet_append(ss = googlesheet, data = data.frame("SUMMARY", 
                                                                    app_session_id, basename(getOption(x = "Azimuth.app.reference")), 
                                                                    ReferenceVersion(object = refs$map), app.env$demo, 
                                                                    app.env$ncellsupload, app.env$ncellspreproc, 
                                                                    as.numeric(x = mapping.time), Sys.Date(), app.env$nanchors, 
                                                                    app.env$clusterpreservationqc)), silent = TRUE)
      }
      print("rendering menu")
      output$menu2 <- renderMenu(expr = {
        sidebarMenu(menuItem(text = "Cell Plots", tabName = "tab_cell", 
                             icon = icon("chart-area")), menuItem(text = "Motif Plots", tabName = "tab_motif", 
                                                                  icon = icon("chart-area")), menuItem(text = "Feature Plots", 
                                                                  tabName = "tab_feature", icon = icon("chart-area")), 
                    menuItem(text = "Download Results", tabName = "tab_download", 
                             icon = icon("file-download")))
      })
      print("at the renaming step")
      print(app.env$object)
      print(length(app.env$query.names))
      app.env$object <- RenameCells(object = app.env$object, 
                                    new.names = app.env$query.names)
      enable(id = "file")
      ToggleDemos(action = "enable", demos = demos)
      react.env$metadata <- TRUE
      react.env$biomarkers <- FALSE
    }
  })
  observeEvent(eventExpr = react.env$metadata, handlerExpr = {
    if (isTRUE(x = react.env$metadata)) {
      print("at metadata")
      metadata.discrete <- sort(x = PlottableMetadataNames(object = app.env$object, 
                                                           exceptions = app.env$metadataxfer, min.levels = 1, 
                                                           max.levels = 50))
      app.env$metadata.discrete <- metadata.discrete
      print(app.env$metadata.discrete)
      for (id in c("metarow", "metacol", "metagroup", "metagroup.motif")) {
        if (id == "metarow") {
          show.metadata <- "query"
        }
        else {
          show.metadata <- paste0("predicted.", app.env$default.metadata)
        }
        updateSelectizeInput(session = session, inputId = id, 
                             choices = metadata.discrete, selected = show.metadata, 
                             server = TRUE, options = selectize.opts)
      }
      updateSelectizeInput(session = session, inputId = "metacolor.query", 
                           choices = c(grep(pattern = "^predicted.", x = metadata.discrete, 
                                            value = TRUE), grep(pattern = "^predicted.", 
                                                                x = metadata.discrete, value = TRUE, invert = TRUE)), 
                           selected = paste0("predicted.", app.env$default.metadata), 
                           server = TRUE, options = selectize.opts[-which(x = names(x = selectize.opts) == 
                                                                            "maxItems")])
      metadata.cont <- sort(x = setdiff(x = colnames(x = app.env$object[[]]), 
                                        y = metadata.discrete))
      print(metadata.cont)
      metadata.cont <- Filter(f = function(x) {
        return(is.numeric(x = app.env$object[[x, drop = TRUE]]))
      }, x = metadata.cont)
      metadata.cont <- sort(x = metadata.cont)
      print(metadata.cont)
      if (any(grepl(pattern = "predicted.*.score", x = metadata.cont))) {
        metadata.cont <- metadata.cont[-grep(pattern = "predicted.*.score", 
                                             x = metadata.cont)]
      }
      if (any(grepl(pattern = "*_refAssay", x = metadata.cont))) {
        metadata.cont <- metadata.cont[-grep(pattern = "*_refAssay", 
                                             x = metadata.cont)]
      }
      max.predictions <- paste0("predicted.", app.env$metadataxfer, 
                                ".score")
      names(x = max.predictions) <- app.env$metadataxfer
      max.predictions <- as.list(x = max.predictions)
      prediction.score.names <- lapply(X = app.env$metadataxfer, 
                                       FUN = function(x) {
                                         key <- Key(object = app.env$object[[paste0("prediction.score.", 
                                                                                    x)]])
                                         ids <- paste0(rownames(x = app.env$object[[paste0("prediction.score.", 
                                                                                           x)]]))
                                         values <- paste0(key, ids)
                                         names(x = values) <- ids
                                         return(values)
                                       })
      print(prediction.score.names)
      names(x = prediction.score.names) <- paste0("Prediction scores - ", 
                                                  app.env$metadataxfer)
      print(prediction.score.names)
      metadata.cont <- c(list(`Max prediction scores` = max.predictions), 
                         prediction.score.names, list(`Other Metadata` = metadata.cont))
      app.env$metadata.cont <- metadata.cont
      print("updating selective size input again")
      updateSelectizeInput(session = session, inputId = "metadata.cont", 
                           choices = app.env$metadata.cont, selected = "", 
                           server = TRUE, options = selectize.opts)
      updateSelectizeInput(session = session, inputId = "metadata.cont.motif", 
                           choices = app.env$metadata.cont, selected = "", 
                           server = TRUE, options = selectize.opts)
      updateSelectizeInput(session = session, inputId = "metacolor.ref", 
                           choices = c(grep(pattern = "^predicted.", x = app.env$metadataxfer, 
                                            value = TRUE), grep(pattern = "^predicted.", 
                                                                x = app.env$metadataxfer, value = TRUE, invert = TRUE)), 
                           selected = app.env$default.metadata, server = TRUE, 
                           options = selectize.opts[-which(x = names(x = selectize.opts) == 
                                                             "maxItems")])
      print("going onto features")
      print("calculating mapping time ")
      mapping.time <- difftime(time1 = Sys.time(), time2 = react.env$start,
                               units = "secs")
      time.fmt <- FormatDiffTime(dt = mapping.time)
      app.env$messages <- c(app.env$messages, time.fmt)
      react.env$chromvar.features <- TRUE
      react.env$features <- TRUE
      react.env$metadata <- FALSE
    }
  })
  observeEvent(eventExpr = react.env$chromvar.features, handlerExpr = {
    if (isTRUE(x = react.env$chromvar.features)) {
      print("doing chromvar features")
      print("printing")
      #DefaultAssay(app.env$object) <- app.env$chromvar.assay
      print("APP ENV CHROMVAR FEATURE DEFAULT")
      print(head(rownames(app.env$object)))
      print(head(row.names(app.env$object[["chromvar"]]@data)))
      app.env$default.chromvar.feature <- ifelse(test = "POU2F3" %in% 
                                                   row.names(x = app.env$object[["chromvar"]]@data), yes = "POU2F3", 
                                                 no = row.names(x = app.env$object[["chromvar"]]@data)[1])
      print(app.env$default.chromvar.feature)
      app.env$chromvar.features <- unique(x = row.names(x = app.env$object[["chromvar"]]@data)) # c(FilterFeatures(features =
      print(head(app.env$chromvar.features))
      print(app.env$default.chromvar.feature %in% app.env$chromvar.features)
      updateSelectizeInput(session = session, inputId = "chromvar.feature", 
                           label = "Motif", choices = app.env$chromvar.features, 
                           selected = app.env$default.chromvar.feature, server = TRUE, 
                           options = selectize.opts)

      if (isTRUE(x = do.adt)) {
        app.env$adt.features <- sort(x = rownames(x = app.env$object[[adt.key]]))
        updateSelectizeInput(session = session, inputId = "adtfeature", 
                             choices = app.env$adt.features, selected = "", 
                             server = TRUE, options = selectize.opts)
      }
      react.env$chromvar.features <- FALSE
      print("finished chromvar features")
    }
  })
  observeEvent(eventExpr = react.env$features, handlerExpr = {
    if (isTRUE(x = react.env$features)) {
      print("doing features")
      DefaultAssay(app.env$object) <- "RNA"
      app.env$default.feature <- ifelse(test = getOption(x = "Azimuth.app.default_gene") %in% 
                                          rownames(x = app.env$object), yes = getOption(x = "Azimuth.app.default_gene"), 
                                        no = VariableFeatures(object = app.env$object)[1])
      print(app.env$default.feature)
      app.env$features <- unique(x = c(FilterFeatures(features = VariableFeatures(object = app.env$object)[1:selectize.opts$maxOptions]), 
                                       FilterFeatures(features = rownames(x = app.env$object))))
      print(head(app.env$features))
      updateSelectizeInput(session = session, inputId = "feature", 
                           label = "Feature", choices = app.env$features, 
                           selected = app.env$default.feature, server = TRUE, 
                           options = selectize.opts)
 
      if (isTRUE(x = do.adt)) {
        app.env$adt.features <- sort(x = rownames(x = app.env$object[[adt.key]]))
        updateSelectizeInput(session = session, inputId = "adtfeature", 
                             choices = app.env$adt.features, selected = "", 
                             server = TRUE, options = selectize.opts)
      }
      react.env$markers <- TRUE
      react.env$features <- FALSE
    }
  })
  observeEvent(eventExpr = input$chromvar.feature, handlerExpr = {
    if (nchar(x = input$chromvar.feature)) {
      print(paste0("CHROMVAR FEATURE IN NEW BLOCK: ", input$chromvar.feature))
      DefaultAssay(app.env$object) <- app.env$chromvar.assay
      app.env$chromvar.feature <- ifelse(test = input$chromvar.feature %in%
                                           rownames(x = app.env$object), yes = paste0(Key(object = app.env$object[[app.env$chromvar.assay]]),
                                                                                      input$chromvar.feature), no = input$chromvar.feature)
      print("got app.env$chromvar.feature")
      print(app.env$chromvar.feature)
    }
  })
  observeEvent(eventExpr = input$feature, handlerExpr = {
    if (nchar(x = input$feature)) {
      print(paste0("FEATURE IN NEW BLOCK ", input$feature))
      app.env$feature <- ifelse(test = input$feature %in%
                                  rownames(x = app.env$object[[app.env$gene.assay]]), yes = paste0(Key(object = app.env$object[[app.env$gene.assay]]),
                                                                             input$feature), no = input$feature)
      print("got app.env$feature")
      print(app.env$feature)
    }
  })
  observeEvent(eventExpr = react.env$markers, handlerExpr = {
    if (isTRUE(x = react.env$markers)) {
      print("doing markers")
      DefaultAssay(app.env$object) <- "RNA"
      allowed.clusters <- names(x = which(x = table(app.env$object[[paste0("predicted.", 
                                                                           app.env$default.metadata)]]) > getOption(x = "Azimuth.de.mincells")))
      print(allowed.clusters)
      allowed.clusters <- factor(x = allowed.clusters, 
                                 levels = unique(x = app.env$object[[paste0("predicted.", 
                                                                            app.env$default.metadata), drop = TRUE]]))
      allowed.clusters <- sort(x = levels(x = droplevels(x = na.omit(object = allowed.clusters))))
      print(allowed.clusters)
      updateSelectizeInput(session = session, inputId = "markerclusters", 
                           choices = allowed.clusters, selected = allowed.clusters[1], 
                           server = TRUE, options = selectize.opts)
      updateSelectizeInput(session = session, inputId = "markerclustersgroup", 
                           choices = app.env$metadataxfer[!app.env$singlepred], 
                           selected = app.env$default.metadata, server = TRUE, 
                           options = selectize.opts)
      updateSelectizeInput(session = session, inputId = "markerclustersgroup.motif", 
                           choices = app.env$metadataxfer[!app.env$singlepred], 
                           selected = app.env$default.metadata, server = TRUE, 
                           options = selectize.opts)
      react.env$get.feature <- TRUE
      react.env$get.chromvar.feature <- TRUE
      react.env$markers <- FALSE
      app.env$disable <- FALSE
    }
  })
  observeEvent(eventExpr = react.env$no, handlerExpr = {
    if (FALSE) {
      react.env$no <- FALSE
    }
  })
  # observeEvent(eventExpr = list(react.env$get.chromvar.feature, input$chromvar.feature), handlerExpr = {
  #   if (isTRUE(x = react.env$get.chromvar.feature)) {
  #     browser()
  #     print("doing get chromvar features")
  #     if (nchar(x = input$chromvar.feature)) {
  #       print("doing input$chromvar.feature")
  #       head(input$chromvar.feature)
  #       DefaultAssay(app.env$object) <- app.env$chromvar.assay
  #       #app.env$chromvar.feature <- ifelse(test = input$chromvar.feature %in%
  #                                            #rownames(x = app.env$object), yes = paste0(Key(object = app.env$object[[app.env$chromvar.assay]]),
  #                                                                                       #input$chromvar.feature), no = input$chromvar.feature)
  #       print("got app.env$chrombar.feature")
  #       print(app.env$chromvar.feature)
  #       print("MARKER CLUSTERS GROUP INPUT: ")
  #       print(input$markerclustersgroup.motif)
  #       print("running render diff motif exp")
  #       table.check <- input$chromvar.feature %in% rownames(x = RenderDiffMotifExp(diff.exp = app.env$chromvar.diff.expr[[paste(app.env$chromvar.assay,
  #                                                                                                                               input$markerclustersgroup.motif, sep = "_")]], groups.use = input$markerclusters,
  #                                                                                  n = Inf))
  #       print(table.check)
  #       tables.clear <- list(adt.proxy, rna.proxy)[c(TRUE,
  #                                                    !table.check)]
  #       for (tab in tables.clear) {
  #         selectRows(proxy = tab, selected = NULL)
  #       }
  #     }
  #   }
  #   react.env$get.chromvar.feature <- FALSE
  # })
  observeEvent(eventExpr = react.env$get.chromvar.feature, handlerExpr = {
      if (isTRUE(x = react.env$get.chromvar.feature)) {
        #browser()
        req(input$markerclusters)
        req(input$markerclustersgroup.motif)
        req(input$chromvar.feature)
        print("doing react.env$get.chromvar.feature")
        print("got app.env$chrombar.feature")
        print(app.env$chromvar.feature)
        print("running render diff motif exp")
        print(paste(app.env$chromvar.assay,input$markerclustersgroup.motif, sep = "_"))
        table.check <- input$chromvar.feature %in% rownames(x = RenderDiffMotifExp(diff.exp = app.env$chromvar.diff.expr[[paste(app.env$chromvar.assay,
                                                                                                                                input$markerclustersgroup.motif, sep = "_")]], groups.use = input$markerclusters,
                                                                                   n = Inf))
        print(table.check)
        tables.clear <- list(adt.proxy, motif.proxy)[c(TRUE,
                                                     !table.check)]
        for (tab in tables.clear) {
          selectRows(proxy = tab, selected = NULL)
        }
      }
    react.env$get.chromvar.feature <- FALSE
  })
  # observeEvent(eventExpr = list(react.env$get.feature, input$feature), handlerExpr = {
  #   if (isTRUE(x = react.env$get.feature)) {
  #     browser()
  #     if (nchar(x = input$feature)) {
  #       print("doing input$feature")
  #       head(input$feature)
  #       app.env$feature <- ifelse(test = input$feature %in%
  #                                   rownames(x = app.env$object), yes = paste0(Key(object = app.env$object[[app.env$gene.assay]]),
  #                                                                              input$feature), no = input$feature)
  #       print("got app.env$feature")
  #       print(app.env$feature)
  #       for (f in c("adtfeature", "metadata.cont")) {
  #         print("updating selective size input")
  #         updateSelectizeInput(session = session, inputId = f,
  #                              choices = list(adtfeature = app.env$adt.features,
  #                                             metadata.cont = app.env$metadata.cont)[[f]],
  #                              selected = "", server = TRUE, options = selectize.opts)
  #       }
  #       for (f in c("feature", "metadata.cont")) {
  #         updateSelectizeInput(session = session, inputId = f, 
  #                              choices = list(feature = app.env$features, 
  #                                             metadata.cont = app.env$metadata.cont)[[f]], 
  #                              selected = "", server = TRUE, options = selectize.opts)
  #       }
  #       print("doing table check")
  #       print("MARKER CLUSTERS GROUP")
  #       print(input$markerclustersgroup)
  #       print("CHECKING FEATURE IN MARKER CLUSTERS GROUP")
  #       print(input$feature)
  #       table.check <- input$feature %in% rownames(x = RenderDiffExp(diff.exp = app.env$diff.expr[[paste(app.env$gene.assay,
  #                                                                                                        input$markerclustersgroup, sep = "_")]], groups.use = input$markerclusters,
  #                                                                    n = Inf))
  #       tables.clear <- list(adt.proxy, rna.proxy)[c(TRUE,
  #                                                    !table.check)]
  #       for (tab in tables.clear) {
  #         selectRows(proxy = tab, selected = NULL)
  #       }
  #     }
  #   }
  #   react.env$get.feature <- FALSE
  # })
  observeEvent(eventExpr = react.env$get.feature, handlerExpr = {
      if (isTRUE(react.env$get.feature)) {
        print("react.env$get.feature")
        print(app.env$feature)
        req(input$markerclustersgroup)
        req(input$feature)
        for (f in c("adtfeature", "metadata.cont")) {
          print("updating selective size input")
          updateSelectizeInput(session = session, inputId = f,
                               choices = list(adtfeature = app.env$adt.features,
                                              metadata.cont = app.env$metadata.cont)[[f]],
                               selected = "", server = TRUE, options = selectize.opts)
        }
        print("doing table check")
        print("MARKER CLUSTERS GROUP")
        print(input$markerclustersgroup)
        print("CHECKING FEATURE IN MARKER CLUSTERS GROUP")
        print(input$feature)
        table.check <- input$feature %in% rownames(x = RenderDiffExp(diff.exp = app.env$diff.expr[[paste(app.env$gene.assay,
                                                                                                         input$markerclustersgroup, sep = "_")]], groups.use = input$markerclusters,
                                                                     n = Inf))
        tables.clear <- list(adt.proxy, rna.proxy)[c(TRUE,
                                                     !table.check)]
        for (tab in tables.clear) {
          selectRows(proxy = tab, selected = NULL)
        }
      }
    react.env$get.feature <- FALSE
  })
  observeEvent(eventExpr = input$adtfeature, handlerExpr = {
    if (nchar(x = input$adtfeature)) {
      app.env$feature <- paste0(Key(object = app.env$object[[adt.key]]), 
                                input$adtfeature)
      for (f in c("feature", "metadata.cont")) {
        updateSelectizeInput(session = session, inputId = f, 
                             choices = list(feature = app.env$features, 
                                            metadata.cont = app.env$metadata.cont)[[f]], 
                             selected = "", server = TRUE, options = selectize.opts)
      }
      table.check <- input$adtfeature %in% rownames(x = RenderDiffExp(diff.exp = app.env$diff.expr[[paste(adt.key, 
                                                                                                          input$markerclustersgroup, sep = "_")]], groups.use = input$markerclusters, 
                                                                      n = Inf))
      tables.clear <- list(rna.proxy, adt.proxy)[c(TRUE, 
                                                   !table.check)]
      for (tab in tables.clear) {
        selectRows(proxy = tab, selected = NULL)
      }
    }
  })
  observeEvent(eventExpr = input$metadata.cont, handlerExpr = {
    if (nchar(x = input$metadata.cont)) {
      if (input$metadata.cont == "mapping.score") {
        if (resolved(x = app.env$mapping.score)) {
          app.env$object$mapping.score <- value(app.env$mapping.score)
        }
      }
      print('doing continued metadata')
      app.env$feature <- input$metadata.cont
      for (f in c("feature", "adtfeature")) {
        updateSelectizeInput(session = session, inputId = f, 
                             choices = list(feature = app.env$features, 
                                            adtfeature = app.env$adt.features)[[f]], 
                             selected = "", server = TRUE, options = selectize.opts)
      }
      for (tab in list(rna.proxy, adt.proxy)) {
        selectRows(proxy = tab, selected = NULL)
      }
    }
  })
  observeEvent(eventExpr = input$metadata.cont.motif, handlerExpr = {
    if (nchar(x = input$metadata.cont.motif)) {
      if (input$metadata.cont.motif == "mapping.score") {
        if (resolved(x = app.env$mapping.score)) {
          app.env$object$mapping.score <- value(app.env$mapping.score)
        }
      }
      print('doing continued metadata')
      app.env$chromvar.feature <- input$metadata.cont.motif
      updateSelectizeInput(session = session, inputId = "chromvar.feature", 
                           choices = app.env$chromvar.features, 
                           selected = "", server = TRUE, options = selectize.opts)
      for (tab in list(rna.proxy, adt.proxy)) { #changed
        selectRows(proxy = tab, selected = NULL)
      }
    }
  })
  observeEvent(eventExpr = input$markerclustersgroup, handlerExpr = {
    if (nchar(x = input$markerclustersgroup)) {
      print("doing markerclusters group")
      allowed.clusters <- names(x = which(x = table(app.env$object[[paste0("predicted.", 
                                                                           input$markerclustersgroup)]]) > getOption(x = "Azimuth.de.mincells")))
      allowed.clusters <- factor(x = allowed.clusters, 
                                 levels = unique(x = app.env$object[[paste0("predicted.", 
                                                                            input$markerclustersgroup), drop = TRUE]]))
      allowed.clusters <- sort(x = levels(x = droplevels(x = na.omit(object = allowed.clusters))))
      app.env$allowedclusters <- allowed.clusters
      updateSelectizeInput(session = session, inputId = "markerclusters", 
                           choices = app.env$allowedclusters, selected = app.env$allowedclusters[1], 
                           server = TRUE, options = selectize.opts)
    }
  })
  observeEvent(eventExpr = input$markerclustersgroup.motif, handlerExpr = {
    if (nchar(x = input$markerclustersgroup.motif)) {
      print("doing markerclusters group")
      allowed.clusters <- names(x = which(x = table(app.env$object[[paste0("predicted.", 
                                                                           input$markerclustersgroup.motif)]]) > getOption(x = "Azimuth.de.mincells")))
      allowed.clusters <- factor(x = allowed.clusters, 
                                 levels = unique(x = app.env$object[[paste0("predicted.", 
                                                                            input$markerclustersgroup.motif), drop = TRUE]]))
      allowed.clusters <- sort(x = levels(x = droplevels(x = na.omit(object = allowed.clusters))))
      app.env$allowedclusters <- allowed.clusters
      updateSelectizeInput(session = session, inputId = "markerclusters.motif", 
                           choices = app.env$allowedclusters, selected = app.env$allowedclusters[1], 
                           server = TRUE, options = selectize.opts)
    }
  })
  observeEvent(eventExpr = input$biomarkers_rows_selected, 
               handlerExpr = {
                 if (length(x = input$biomarkers_rows_selected)) {
                   print("doing biomarkers row selected")
                   updateSelectizeInput(session = session, inputId = "feature", 
                                        choices = app.env$features, selected = rownames(x = RenderDiffExp(diff.exp = app.env$diff.expr[[paste(app.env$gene.assay, 
                                                                                                                                              input$markerclustersgroup, sep = "_")]], 
                                                                                                          groups.use = input$markerclusters, n = Inf))[input$biomarkers_rows_selected], 
                                        server = TRUE, options = selectize.opts)
                 }
                 print("input biomarkers rows selected FEATURE:")
                 print(input$feature)
               })
  observeEvent(eventExpr = input$motif_rows_selected, 
               handlerExpr = {
                 if (length(x = input$motif_rows_selected)) {
                   print("doing motif row selected")
                   updateSelectizeInput(session = session, inputId = "motif", 
                                        choices = app.env$chromvar.features, selected = rownames(x = RenderDiffMotifExp(diff.exp = app.env$chromvar.diff.expr[[paste(app.env$chromvar.assay, 
                                                                                                                                                                  input$markerclustersgroup.motif, sep = "_")]], groups.use = input$markerclusters.motif, 
                                                                                                                        n = Inf))[input$motif_rows_selected], 
                                        server = TRUE, options = selectize.opts)
                 }
               })
  observeEvent(eventExpr = input$adtbio_rows_selected, handlerExpr = {
    if (length(x = input$adtbio_rows_selected)) {
      updateSelectizeInput(session = session, inputId = "adtfeature", 
                           choices = app.env$adt.features, selected = rownames(x = RenderDiffExp(diff.exp = app.env$diff.expr[[paste(adt.key, 
                                                                                                                                     input$markerclustersgroup, sep = "_")]], groups.use = input$markerclusters, 
                                                                                                 n = Inf))[input$adtbio_rows_selected], server = TRUE, 
                           options = selectize.opts)
    }
  })
  observeEvent(eventExpr = input$metacolor.ref, handlerExpr = {
    if (length(x = unique(x = as.vector(x = refs$plot[[input$metacolor.ref[1], 
                                                       drop = TRUE]]))) >= 30) {
      if (isTRUE(app.env$fresh.plot)) {
        updateCheckboxInput(session = session, inputId = "labels", 
                            value = TRUE)
        updateCheckboxInput(session = session, inputId = "legend", 
                            value = FALSE)
        app.env$fresh.plot <- FALSE
      }
    }
  })
  observeEvent(eventExpr = input$submit_feedback, handlerExpr = {
    if (!is.null(x = googlesheet)) {
      try(expr = sheet_append(ss = googlesheet, data = data.frame("FEEDBACK", 
                                                                  app_session_id, paste0("feedback: \"", input$feedback, 
                                                                                         "\""))))
    }
    updateTextAreaInput(session = session, inputId = "feedback", 
                        label = NULL, value = "Thank you for your feedback!")
  })
  observeEvent(eventExpr = input$showrefonly, handlerExpr = {
    if (!is.null(app.env$metadata.discrete)) {
      if (input$showrefonly) {
        disable(id = "metacolor.query")
        enable(id = "metacolor.ref")
      }
      else {
        disable(id = "metacolor.ref")
        enable(id = "metacolor.query")
      }
    }
  })
  output$plot.qc <- renderPlot(expr = {
    if (!is.null(x = isolate(expr = app.env$object)) & isTRUE(x = react.env$plot.qc)) {
      print("making some qc plots")
      qc <- paste0(c("nCount_", "nFeature_"), app.env$default.assay)
      # if (isTRUE(x = react.env$mt)) { # removing this for now bc i dont need it 
      #   qc <- c(qc, mt.key)
      # }
      print(qc)
      check.qcpoints <- "qcpoints" %in% input$check.qc
      check.qcscale <- "qcscale" %in% input$check.qc
      print("about to make vln plots ")
      vln_1 <- VlnPlot(object = isolate(app.env$object),
                       features = qc[1], group.by = "query", combine = FALSE,
                       pt.size = ifelse(test = check.qcpoints, yes = 0,
                                        no = Seurat:::AutoPointSize(data = isolate(app.env$object))),
                       log = check.qcscale)
      vln_2 <- VlnPlot(object = isolate(app.env$object),
                       features = qc[2], group.by = "query", combine = FALSE,
                       pt.size = ifelse(test = check.qcpoints, yes = 0,
                                        no = Seurat:::AutoPointSize(data = isolate(app.env$object))),
                       log = check.qcscale)
      #dist <- OverlapDistPlot(query_assay = isolate(app.env$object), 
      #multiome = refs$bridge)
      
      vlnlist <- c(vln_1, vln_2) # dist
      print(input$num.ncountmin)
      print(input$num.ncountmax)
      print(input$num.nfeaturemin)
      print(input$num.nfeaturemax)
      vlnlist[[1]] <- vlnlist[[1]] + geom_hline(yintercept = input$num.ncountmin) +
        geom_hline(yintercept = input$num.ncountmax) +
        annotate(geom = "rect", alpha = 0.2, fill = "red",
                 ymin = input$num.ncountmax, ymax = Inf, xmin = 0.5,
                 xmax = 1.5) + annotate(geom = "rect", alpha = 0.2,
                                        fill = "red", ymin = ifelse(test = check.qcscale,
                                                                    yes = 0, no = -Inf), ymax = input$num.ncountmin,
                                        xmin = 0.5, xmax = 1.5) + NoLegend() + xlab("")
      print("vlnlist1")
      vlnlist[[2]] <- vlnlist[[2]] + geom_hline(yintercept = input$num.nfeaturemin) +
        geom_hline(yintercept = input$num.nfeaturemax) +
        annotate(geom = "rect", alpha = 0.2, fill = "red",
                 ymin = input$num.nfeaturemax, ymax = Inf, xmin = 0.5,
                 xmax = 1.5) + annotate(geom = "rect", alpha = 0.2,
                                        fill = "red", ymin = ifelse(test = check.qcscale,
                                                                    yes = 0, no = -Inf), ymax = input$num.nfeaturemin,
                                        xmin = 0.5, xmax = 1.5) + NoLegend() + xlab("")
      # print("vlnlist2")
      # if (react.env$mt) { # trying deleting this bc i think its giving me issues and i dont need it 
      #   print("trying the mt")
      #   vlnlist[[3]] <- vlnlist[[3]] + geom_hline(yintercept = input$minmt) + 
      #     geom_hline(yintercept = input$maxmt) + annotate(geom = "rect", 
      #                                                     alpha = 0.2, fill = "red", ymin = input$maxmt, 
      #                                                     ymax = Inf, xmin = 0.5, xmax = 1.5) + annotate(geom = "rect", 
      #                                                                                                    alpha = 0.2, fill = "red", ymin = ifelse(test = check.qcscale, 
      #                                                                                                                                             yes = 0, no = -Inf), ymax = input$minmt, 
      #                                                                                                    xmin = 0.5, xmax = 1.5) + NoLegend() + xlab("")
      # 
      # }
      # plotlist[[3]] <- pl
      wrap_plots(vlnlist, ncol = length(x = vlnlist))
    }
  })
  output$dist.qc <- renderPlot(expr = {
    if (!is.null(x = isolate(expr = app.env$chromatin_assay_1)) & isTRUE(x = react.env$dist.qc)) {
      print("making dist plots")
      dist <- OverlapDistPlot(query_assay = isolate(app.env$chromatin_assay_1),
                              multiome = refs$bridge[["ATAC"]])
    }
  })
  
  output$refdim_intro <- renderPlot(expr = {
    app.env$plots.refdim_intro_df <- cbind(as.data.frame(x = Embeddings(object = refs$plot[["refUMAP"]])), 
                                           refs$plot[[]])
    
    p <- DimPlot(object = refs$plot, combine = FALSE, group.by = default_xfer, 
                 cols = GetColorMap(object = refs$map)[[default_xfer]], 
                 repel = TRUE, label = TRUE, raster = FALSE)[[1]]
    app.env$plot.ranges <- list(layer_scales(p)$x$range$range, 
                                layer_scales(p)$y$range$range)
    p + WelcomePlot()
  })
  output$refdim_intro_hover_box <- renderUI({
    hover <- input$refdim_intro_hover_location
    df <- app.env$plots.refdim_intro_df
    if (!is.null(x = hover)) {
      hover[["mapping"]] <- setNames(object = as.list(x = colnames(x = app.env$plots.refdim_intro_df)[1:2]), 
                                     nm = c("x", "y"))
    }
    point <- nearPoints(df = df, coordinfo = hover, threshold = 10, 
                        maxpoints = 1, addDist = TRUE)
    if (nrow(x = point) == 0) {
      return(NULL)
    }
    hovertext <- do.call(what = paste0, args = lapply(X = metadata.annotate, 
                                                      FUN = function(md) {
                                                        paste0("<span>", md, "</span>: <i>", point[[md]], 
                                                               "</i><br>")
                                                      }))
    wellPanel(style = HoverBoxStyle(x = hover$coords_css$x, 
                                    y = hover$coords_css$y), p(HTML(text = hovertext)))
  })
  output$refdim <- renderPlot(expr = {
    if (!is.null(x = input$metacolor.ref)) {
      colormaps <- GetColorMap(object = refs$map)[input$metacolor.ref]
      if (length(x = colormaps) == 1) {
        app.env$plots.refdim_df <- app.env$plots.refdim_intro_df
        DimPlot(object = refs$plot, label = isTRUE("labels" %in% 
                                                     input$dimplot.opts), group.by = input$metacolor.ref, 
                cols = colormaps[[1]], repel = TRUE, raster = FALSE)[[1]] + 
          labs(x = "UMAP 1", y = "UMAP 2") + if (isFALSE(x = "legend" %in% 
                                                         input$dimplot.opts) | OversizedLegend(refs$plot[[input$metacolor.ref, 
                                                                                                          drop = TRUE]])) 
            NoLegend()
      }
      else {
        app.env$plots.refdim_df <- NULL
        plots <- list()
        for (i in 1:length(x = colormaps)) {
          plots[[i]] <- DimPlot(object = refs$plot, label = isTRUE("labels" %in% 
                                                                     input$dimplot.opts), group.by = input$metacolor.ref[i], 
                                cols = colormaps[[i]], repel = TRUE, raster = FALSE) + 
            labs(x = "UMAP 1", y = "UMAP 2") + if (isFALSE(x = "legend" %in% 
                                                           input$dimplot.opts) | OversizedLegend(refs$plot[[input$metacolor.ref[i], 
                                                                                                            drop = TRUE]])) 
              NoLegend()
        }
        wrap_plots(plots, nrow = 1)
      }
    }
  })
  output$refdim_hover_box <- renderUI({
    if (!is.null(x = app.env$plots.refdim_df)) {
      hover <- input$refdim_hover_location
      df <- app.env$plots.refdim_df
      if (!is.null(x = hover)) {
        hover[["mapping"]] <- setNames(object = as.list(x = colnames(x = app.env$plots.refdim_intro_df)[1:2]), 
                                       nm = c("x", "y"))
      }
      point <- nearPoints(df = df, coordinfo = hover, threshold = 10, 
                          maxpoints = 1, addDist = TRUE)
      if (nrow(x = point) == 0) {
        return(NULL)
      }
      hovertext <- do.call(what = paste0, args = as.list(c(paste0("<b>", 
                                                                  point[[input$metacolor.ref]], "</b><br>"), sapply(X = setdiff(possible.metadata.transfer, 
                                                                                                                                input$metacolor.ref), FUN = function(md) {
                                                                                                                                  paste0("<span>", md, "</span>: <i>", point[[md]], 
                                                                                                                                         "</i><br>")
                                                                                                                                }))))
      wellPanel(style = HoverBoxStyle(x = hover$coords_css$x, 
                                      y = hover$coords_css$y), p(HTML(text = hovertext)))
    }
  })
  output$objdim <- renderPlot(expr = {
    if (!is.null(x = app.env$object) && app.env$disable == 
        FALSE) {
      if (is.null(x = app.env$emptyref) | is.null(x = app.env$merged)) {
        app.env$emptyref <- refs$plot
        Idents(object = app.env$emptyref) <- "."
        for (md in colnames(x = app.env$emptyref[[]])) {
          app.env$emptyref[[md]] <- "."
        }
        for (md in setdiff(x = colnames(x = app.env$object[[]]), 
                           y = colnames(x = app.env$emptyref[[]]))) {
          app.env$emptyref[[md]] <- "."
        }
        app.env$object[["refUMAP"]] <- app.env$object[["umap.proj"]]
        app.env$merged <- merge(app.env$emptyref, app.env$object, 
                                merge.dr = "refUMAP")
      }
      if (isFALSE(x = input$showrefonly) & length(x = Reductions(object = app.env$object)) & 
          !is.null(x = input$metacolor.query)) {
        app.env$plots.refdim_df <- NULL
        if (length(x = input$metacolor.query) == 1) {
          group.var <- gsub(pattern = "^predicted.", 
                            replacement = "", x = input$metacolor.query)
          colormap <- GetColorMap(object = refs$map)[[group.var]]
          if (!grepl(pattern = "^predicted.", x = input$metacolor.query)) {
            colormap <- CreateColorMap(ids = unique(as.vector(app.env$object[[input$metacolor.query, 
                                                                              drop = T]])))
          }
          colormap["."] <- "#F1F1F1"
          app.env$plots.objdim_df <- cbind(as.data.frame(x = Embeddings(object = app.env$object[["umap.proj"]])), 
                                           app.env$object[[]])
          p <- DimPlot(object = app.env$merged, group.by = input$metacolor.query, 
                       label = FALSE, cols = colormap[names(x = colormap) %in% 
                                                        c(".", unique(x = as.vector(x = app.env$object[[input$metacolor.query, 
                                                                                                        drop = TRUE]])))], repel = TRUE, raster = FALSE, 
                       reduction = "refUMAP")[[1]] + xlim(app.env$plot.ranges[[1]]) + 
            ylim(app.env$plot.ranges[[2]]) + labs(x = "UMAP 1", 
                                                  y = "UMAP 2") + if (isFALSE(x = input$legend)) 
                                                    NoLegend()
          if (isTRUE(x = "labels" %in% input$label.opts)) {
            keep <- if (isTRUE(x = "filterlabels" %in% 
                               input$label.opts)) {
              t <- table(as.vector(x = app.env$object[[input$metacolor.query, 
                                                       drop = TRUE]]))
              names(x = t)[which(x = t > 0.02 * ncol(x = app.env$object))]
            }
            else NULL
            return(LabelClusters(plot = p, id = input$metacolor.query, 
                                 clusters = keep))
          }
          return(p)
        }
        else {
          app.env$plots.objdim_df <- NULL
          plots <- list()
          for (i in 1:length(x = input$metacolor.query)) {
            group.var <- gsub(pattern = "^predicted.", 
                              replacement = "", x = input$metacolor.query[i])
            colormap <- GetColorMap(object = refs$map)[[group.var]]
            if (!grepl(pattern = "^predicted.", x = input$metacolor.query[i])) {
              colormap <- CreateColorMap(ids = unique(x = as.vector(x = app.env$object[[input$metacolor.query[i], 
                                                                                        drop = TRUE]])))
            }
            colormap["."] <- "#F1F1F1"
            p <- DimPlot(object = app.env$merged, group.by = input$metacolor.query[i], 
                         cols = colormap[names(x = colormap) %in% 
                                           c(".", unique(x = as.vector(x = app.env$object[[input$metacolor.query[i], 
                                                                                           drop = TRUE]])))], repel = TRUE, raster = FALSE, 
                         reduction = "refUMAP")[[1]] + xlim(app.env$plot.ranges[[1]]) + 
              ylim(app.env$plot.ranges[[2]]) + labs(x = "UMAP 1", 
                                                    y = "UMAP 2") + if (isFALSE(x = input$legend) | 
                                                                        OversizedLegend(annotation.list = app.env$object[[input$metacolor.query[i], 
                                                                                                                          drop = TRUE]])) 
                                                      NoLegend()
            if (isTRUE("labels" %in% input$label.opts)) {
              keep <- if (isTRUE(x = "filterlabels" %in% 
                                 input$label.opts)) {
                t <- table(as.vector(x = app.env$object[[input$metacolor.query[i], 
                                                         drop = TRUE]]))
                print(names(x = t)[which(t > 0.02 * ncol(x = app.env$object))])
                names(x = t)[which(t > 0.02 * ncol(x = app.env$object))]
              }
              else NULL
              plots[[i]] <- LabelClusters(plot = p, id = input$metacolor.query[i], 
                                          clusters = keep)
            }
            else {
              plots[[i]] <- p
            }
          }
          wrap_plots(plots, nrow = 1)
        }
      }
      else {
        app.env$plots.objdim_df <- NULL
        if (!is.null(x = input$metacolor.ref)) {
          colormaps <- GetColorMap(object = refs$map)[input$metacolor.ref]
          if (length(x = colormaps) == 1) {
            app.env$plots.refdim_df <- app.env$plots.refdim_intro_df
            DimPlot(object = refs$plot, label = isTRUE("labels" %in% 
                                                         input$label.opts), group.by = input$metacolor.ref, 
                    cols = colormaps[[1]], repel = TRUE, raster = FALSE)[[1]] + 
              labs(x = "UMAP 1", y = "UMAP 2") + if (isFALSE(input$legend) | 
                                                     OversizedLegend(annotation.list = refs$plot[[input$metacolor.ref, 
                                                                                                  drop = TRUE]])) 
                NoLegend()
          }
          else {
            app.env$plots.refdim_df <- NULL
            plots <- list()
            for (i in 1:length(x = colormaps)) {
              plots[[i]] <- DimPlot(object = refs$plot, 
                                    label = isTRUE("labels" %in% input$label.opts), 
                                    group.by = input$metacolor.ref[i], cols = colormaps[[i]], 
                                    repel = TRUE, raster = FALSE) + labs(x = "UMAP 1", 
                                                                         y = "UMAP 2") + if (isFALSE(x = input$legend) | 
                                                                                             OversizedLegend(annotation.list = refs$plot[[input$metacolor.ref[i], 
                                                                                                                                          drop = TRUE]])) 
                                                                           NoLegend()
            }
            wrap_plots(plots, nrow = 1)
          }
        }
      }
    }
  })
  output$querydim <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      if (length(x = Reductions(object = app.env$object)) & 
          !is.null(x = input$metacolor.query)) {
        if (length(x = input$metacolor.query) == 1) {
          group.var <- gsub(pattern = "^predicted.", 
                            replacement = "", x = input$metacolor.query)
          colormap <- GetColorMap(object = refs$map)[[group.var]]
          if (!grepl(pattern = "^predicted.", x = input$metacolor.query)) {
            colormap <- NULL
          }
          app.env$plots.querydim_df <- cbind(as.data.frame(x = Embeddings(object = app.env$object[["umap.proj"]])), 
                                             app.env$object[[]])
          DimPlot(object = app.env$object, group.by = input$metacolor.query, 
                  label = isTRUE("labels" %in% input$dimplot.opts), 
                  cols = colormap[names(x = colormap) %in% 
                                    unique(x = app.env$object[[input$metacolor.query, 
                                                               drop = TRUE]])], repel = TRUE, reduction = "umap.proj")[[1]] + 
            xlim(app.env$plot.ranges[[1]]) + ylim(app.env$plot.ranges[[2]]) + 
            labs(x = "UMAP 1", y = "UMAP 2") + if (isFALSE(x = "legend" %in% 
                                                           input$dimplot.opts) | OversizedLegend(app.env$object[[input$metacolor.query, 
                                                                                                                 drop = TRUE]])) 
              NoLegend()
        }
        else {
          app.env$plots.querydim_df <- NULL
          plots <- list()
          for (i in 1:length(x = input$metacolor.query)) {
            group.var <- gsub(pattern = "^predicted.", 
                              replacement = "", x = input$metacolor.query[i])
            colormap <- GetColorMap(object = refs$map)[[group.var]]
            if (!grepl(pattern = "^predicted.", x = input$metacolor.query[i])) {
              colormap <- NULL
            }
            plots[[i]] <- DimPlot(object = app.env$object, 
                                  group.by = input$metacolor.query[i], label = isTRUE("labels" %in% 
                                                                                        input$dimplot.opts), cols = colormap[names(x = colormap) %in% 
                                                                                                                               unique(x = app.env$object[[input$metacolor.query[i], 
                                                                                                                                                          drop = TRUE]])], repel = TRUE, reduction = "umap.proj") + 
              xlim(app.env$plot.ranges[[1]]) + ylim(app.env$plot.ranges[[2]]) + 
              labs(x = "UMAP 1", y = "UMAP 2") + if (isFALSE(x = "legend" %in% 
                                                             input$dimplot.opts) | OversizedLegend(app.env$object[[input$metacolor.query[i], 
                                                                                                                   drop = TRUE]])) 
                NoLegend()
          }
          wrap_plots(plots, nrow = 1)
        }
      }
    }
  })
  output$objdim_hover_box <- renderUI({
    if (!is.null(x = app.env$plots.objdim_df)) {
      hover <- input$objdim_hover_location
      df <- app.env$plots.objdim_df
      if (!is.null(x = hover)) {
        hover[["mapping"]] <- setNames(object = as.list(x = colnames(x = app.env$plots.objdim_df)[1:2]), 
                                       nm = c("x", "y"))
      }
      point <- nearPoints(df = df, coordinfo = hover, threshold = 10, 
                          maxpoints = 1, addDist = TRUE)
      if (nrow(x = point) == 0) {
        return(NULL)
      }
      hovertext <- do.call(what = paste0, args = as.list(c(paste0("<b>", 
                                                                  point[[input$metacolor.query]], "</b><br>"), 
                                                           if (grepl(pattern = "^predicted.", x = input$metacolor.query)) {
                                                             paste0("<i>prediction score</i>: <span>", format(x = round(x = point[[paste0(input$metacolor.query, 
                                                                                                                                          ".score")]], digits = 2), nsmall = 2), "</span><br>")
                                                           })))
      wellPanel(style = HoverBoxStyle(x = hover$coords_css$x, 
                                      y = hover$coords_css$y), p(HTML(text = hovertext)))
    }
    else if (!is.null(x = app.env$plots.refdim_df)) {
      hover <- input$objdim_hover_location
      df <- app.env$plots.refdim_df
      if (!is.null(x = hover)) {
        hover[["mapping"]] <- setNames(object = as.list(x = colnames(x = app.env$plots.refdim_intro_df)[1:2]), 
                                       nm = c("x", "y"))
      }
      point <- nearPoints(df = df, coordinfo = hover, threshold = 10, 
                          maxpoints = 1, addDist = TRUE)
      if (nrow(x = point) == 0) {
        return(NULL)
      }
      hovertext <- do.call(what = paste0, args = as.list(c(paste0("<b>", 
                                                                  point[[input$metacolor.ref]], "</b><br>"), sapply(X = setdiff(metadata.annotate, 
                                                                                                                                input$metacolor.ref), FUN = function(md) {
                                                                                                                                  paste0("<span>", md, "</span>: <i>", point[[md]], 
                                                                                                                                         "</i><br>")
                                                                                                                                }))))
      wellPanel(style = HoverBoxStyle(x = hover$coords_css$x, 
                                      y = hover$coords_css$y), p(HTML(text = hovertext)))
    }
  })
  output$querydim_hover_box <- renderUI({
    if (!is.null(x = app.env$plots.querydim_df)) {
      hover <- input$querydim_hover_location
      df <- app.env$plots.querydim_df
      if (!is.null(x = hover)) {
        hover[["mapping"]] <- setNames(object = as.list(x = colnames(x = app.env$plots.querydim_df)[1:2]), 
                                       nm = c("x", "y"))
      }
      point <- nearPoints(df = df, coordinfo = hover, threshold = 10, 
                          maxpoints = 1, addDist = TRUE)
      if (nrow(x = point) == 0) {
        return(NULL)
      }
      hovertext <- do.call(what = paste0, args = as.list(c(paste0("<b>", 
                                                                  point[[input$metacolor.query]], "</b><br>"), 
                                                           if (grepl(pattern = "^predicted.", x = input$metacolor.query)) {
                                                             paste0("<i>prediction score</i>: <span>", format(x = round(x = point[[paste0(input$metacolor.query, 
                                                                                                                                          ".score")]], digits = 2), nsmall = 2), "</span><br>")
                                                           })))
      wellPanel(style = HoverBoxStyle(x = hover$coords_css$x, 
                                      y = hover$coords_css$y), p(HTML(text = hovertext)))
    }
  })
  output$motifvln <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      print("RENDERING VLN PLOT FOR MOTIF")
      avail <- c(paste0(Key(object = app.env$object[[app.env$chromvar.assay]]), 
                        rownames(x = app.env$object[[app.env$chromvar.assay]])), colnames(x = app.env$object[[]]))
      prediction.names <- unlist(x = lapply(X = app.env$metadataxfer, 
                                            FUN = function(x) {
                                              assay <- paste0("prediction.score.", x)
                                              pred <- rep(x = x, times = nrow(x = app.env$object[[assay]]))
                                              names(x = pred) <- paste0(Key(object = app.env$object[[assay]]), 
                                                                        rownames(x = app.env$object[[assay]]))
                                              return(pred)
                                            }))
      print("PREDICTION NAMES")
      print(head(prediction.names))
      max.pred.names <- paste0("predicted.", app.env$metadataxfer, 
                               ".score")
      avail <- c(avail, names(x = prediction.names))
      print("AVAIL")
      print(head(avail))
      print("second if statement")
      print(app.env$chromvar.feature)
      print(app.env$chromvar.feature %in% avail)
      if (app.env$chromvar.feature %in% avail) {
        print("went into the second if statement")
        if (app.env$chromvar.feature == "mapping.score" && !resolved(x = app.env$mapping.score)) {
          ggplot() + annotate("text", x = 4, y = 25, 
                              size = 8, label = "Mapping score still computing ... ") + 
            theme_void()
        }
        else {
          title <- ifelse(test = grepl(pattern = "^chromvar_", 
                                       x = app.env$chromvar.feature), yes = gsub(pattern = "^chromvar_", 
                                                                        replacement = "", x = app.env$chromvar.feature), no = app.env$chromvar.feature)
          if (app.env$chromvar.feature %in% names(x = prediction.names)) {
            print("went into the third if statement")
            pred <- strsplit(x = app.env$chromvar.feature, split = "_")[[1]][2]
            group <- prediction.names[app.env$chromvar.feature]
            title <- paste0("Prediction Score (", group, 
                            ") ", pred)
          }
          if (app.env$chromvar.feature %in% max.pred.names) {
            print("went into the fourth print statement")
            pred <- gsub(pattern = "predicted.", replacement = "", 
                         x = app.env$chromvar.feature)
            pred <- gsub(pattern = ".score", replacement = "", 
                         x = pred)
            title <- paste0("Max Prediction Score - ", 
                            pred)
          }
          print("making vln plot")
          print("MOTIF FOR VLN PLOT")
          print(app.env$chromvar.feature)
          VlnPlot(object = app.env$object, features = app.env$chromvar.feature, 
                  group.by = input$metagroup.motif, pt.size = ifelse(test = input$check.featpoints, 
                                                               yes = 0, no = Seurat:::AutoPointSize(data = app.env$object))) + 
            ggtitle(label = title) + NoLegend()
        }
      }
    }
  })
  output$evln <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      print("EVLN")
      DefaultAssay(app.env$object) <- app.env$gene.assay
      print(head(rownames(x = app.env$object)))
      avail <- c(paste0(Key(object = app.env$object[[app.env$gene.assay]]), 
                        rownames(x = app.env$object)), colnames(x = app.env$object[[]]))
      prediction.names <- unlist(x = lapply(X = app.env$metadataxfer, 
                                            FUN = function(x) {
                                              assay <- paste0("prediction.score.", x)
                                              pred <- rep(x = x, times = nrow(x = app.env$object[[assay]]))
                                              names(x = pred) <- paste0(Key(object = app.env$object[[assay]]), 
                                                                        rownames(x = app.env$object[[assay]]))
                                              return(pred)
                                            }))
      max.pred.names <- paste0("predicted.", app.env$metadataxfer, 
                               ".score")
      avail <- c(avail, names(x = prediction.names))
      print(head(avail))
      if (do.adt) {
        avail <- c(avail, paste0(Key(object = app.env$object[[adt.key]]), 
                                 rownames(x = app.env$object[[adt.key]])))
      }
      print(app.env$feature %in% avail)
      print(app.env$feature)
      if (app.env$feature %in% avail) {
        if (app.env$feature == "mapping.score" && !resolved(x = app.env$mapping.score)) {
          ggplot() + annotate("text", x = 4, y = 25, 
                              size = 8, label = "Mapping score still computing ... ") + 
            theme_void()
        }
        else {
          title <- ifelse(test = grepl(pattern = "^RNA_", 
                                       x = app.env$feature), yes = gsub(pattern = "^RNA_", 
                                                                        replacement = "", x = app.env$feature), no = app.env$feature)
          if (app.env$feature %in% names(x = prediction.names)) {
            pred <- strsplit(x = app.env$feature, split = "_")[[1]][2]
            group <- prediction.names[app.env$feature]
            title <- paste0("Prediction Score (", group, 
                            ") ", pred)
          }
          if (app.env$feature %in% max.pred.names) {
            pred <- gsub(pattern = "predicted.", replacement = "", 
                         x = app.env$feature)
            pred <- gsub(pattern = ".score", replacement = "", 
                         x = pred)
            title <- paste0("Max Prediction Score - ", 
                            pred)
          }
          print("FEATURE FOR VLN PLOT")
          print(app.env$feature)
          print("GROUP BY")
          print(input$metagroup)
          VlnPlot(object = app.env$object, features = app.env$feature, 
                  group.by = input$metagroup, pt.size = ifelse(test = input$check.featpoints, 
                                                               yes = 0, no = Seurat:::AutoPointSize(data = app.env$object))) + 
            ggtitle(label = title) + NoLegend()
        }
      }
    }
  })
  output$motifdim <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      print("MAKING MOTIF DIM PLOT")
      palettes <- list(c("lightgrey", "blue"), c("lightgrey", 
                                                 "darkred"))
      names(x = palettes) <- c(Key(object = app.env$object[[app.env$chromvar.assay]]), 
                               "md_")
      prediction.names <- unlist(x = lapply(X = app.env$metadataxfer, 
                                            FUN = function(x) {
                                              assay <- paste0("prediction.score.", x)
                                              pred <- rep(x = x, times = nrow(x = app.env$object[[assay]]))
                                              names(x = pred) <- paste0(Key(object = app.env$object[[assay]]), 
                                                                        rownames(x = app.env$object[[assay]]))
                                              return(pred)
                                            }))
      max.pred.names <- paste0("predicted.", app.env$metadataxfer, 
                               ".score")
      md <- c(colnames(x = app.env$object[[]]), names(x = prediction.names))
      print("MD")
      print(head(md))
      feature.key <- if (app.env$chromvar.feature %in% md) {
        "md_"
      }
      else {
        paste0(unlist(x = strsplit(x = app.env$chromvar.feature, 
                                   split = "_"))[1], "_")
      }
      print("FEATURE KEY:")
      print(head(feature.key))
      pal.use <- palettes[[feature.key]]
      print("PAL.USE:")
      print(head(pal.use))
      if (!is.null(x = pal.use)) {
        if (app.env$chromvar.feature == "mapping.score" && !resolved(x = app.env$mapping.score)) {
          ggplot() + annotate("text", x = 4, y = 25, 
                              size = 8, label = "Mapping score still computing ... ") + 
            theme_void()
        }
        else {
          title <- ifelse(test = grepl(pattern = "^chromvar_", 
                                       x = app.env$feature), yes = gsub(pattern = "^chromvar_", 
                                                                        replacement = "", x = app.env$chromvar.feature), no = app.env$chromvar.feature)
          if (app.env$chromvar.feature %in% names(x = prediction.names)) {
            pred <- strsplit(x = app.env$chromvar.feature, split = "_")[[1]][2]
            group <- prediction.names[app.env$chromvar.feature]
            title <- paste0("Prediction Score (", group, 
                            ") ", pred)
          }
          if (app.env$chromvar.feature %in% max.pred.names) {
            pred <- gsub(pattern = "predicted.", replacement = "", 
                         x = app.env$chromvar.feature)
            pred <- gsub(pattern = ".score", replacement = "", 
                         x = pred)
            title <- paste0("Max Prediction Score - ", 
                            pred)
          }
          suppressWarnings(expr = FeaturePlot(object = app.env$object, 
                                              features = app.env$chromvar.feature, cols = pal.use, 
                                              min.cutoff = 'q10', max.cutoff = 'q90',  
                                              reduction = "umap.proj")) + xlim(app.env$plot.ranges[[1]]) + 
            ylim(app.env$plot.ranges[[2]]) + ggtitle(label = title)
        }
      }
    }
  })
  output$edim <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      palettes <- list(c("lightgrey", "blue"), c("lightgrey", 
                                                 "darkred"))
      names(x = palettes) <- c(Key(object = app.env$object[[app.env$gene.assay]]), 
                               "md_")
      if (do.adt) {
        palettes[[Key(object = app.env$object[[adt.key]])]] <- c("lightgrey", 
                                                                 "darkgreen")
      }
      prediction.names <- unlist(x = lapply(X = app.env$metadataxfer, 
                                            FUN = function(x) {
                                              assay <- paste0("prediction.score.", x)
                                              pred <- rep(x = x, times = nrow(x = app.env$object[[assay]]))
                                              names(x = pred) <- paste0(Key(object = app.env$object[[assay]]), 
                                                                        rownames(x = app.env$object[[assay]]))
                                              return(pred)
                                            }))
      max.pred.names <- paste0("predicted.", app.env$metadataxfer, 
                               ".score")
      md <- c(colnames(x = app.env$object[[]]), names(x = prediction.names))
      feature.key <- if (app.env$feature %in% md) {
        "md_"
      }
      else {
        paste0(unlist(x = strsplit(x = app.env$feature, 
                                   split = "_"))[1], "_")
      }
      pal.use <- palettes[[feature.key]]
      if (!is.null(x = pal.use)) {
        if (app.env$feature == "mapping.score" && !resolved(x = app.env$mapping.score)) {
          ggplot() + annotate("text", x = 4, y = 25, 
                              size = 8, label = "Mapping score still computing ... ") + 
            theme_void()
        }
        else {
          title <- ifelse(test = grepl(pattern = "^RNA_", 
                                       x = app.env$feature), yes = gsub(pattern = "^RNA_", 
                                                                        replacement = "", x = app.env$feature), no = app.env$feature)
          if (app.env$feature %in% names(x = prediction.names)) {
            pred <- strsplit(x = app.env$feature, split = "_")[[1]][2]
            group <- prediction.names[app.env$feature]
            title <- paste0("Prediction Score (", group, 
                            ") ", pred)
          }
          if (app.env$feature %in% max.pred.names) {
            pred <- gsub(pattern = "predicted.", replacement = "", 
                         x = app.env$feature)
            pred <- gsub(pattern = ".score", replacement = "", 
                         x = pred)
            title <- paste0("Max Prediction Score - ", 
                            pred)
          }
          suppressWarnings(expr = FeaturePlot(object = app.env$object, 
                                              features = app.env$feature, cols = pal.use, 
                                              reduction = "umap.proj")) + xlim(app.env$plot.ranges[[1]]) + 
            ylim(app.env$plot.ranges[[2]]) + ggtitle(label = title)
        }
      }
    }
  })
  output$message <- renderUI(expr = {
    p(HTML(text = paste(app.env$messages, collapse = "<br />")))
  })
  output$containerid <- renderUI(expr = {
    p(HTML(text = paste(paste("debug ID:", Sys.info()[["nodename"]]), 
                        paste("Azimuth version:", packageVersion(pkg = "Azimuth")), 
                        paste("Seurat version:", packageVersion(pkg = "Seurat")), 
                        paste("Reference version:", ReferenceVersion(object = refs$map)), 
                        sep = "<br />")))
  })
  output$text.cellsremain <- renderText(expr = {
    if (!is.null(x = isolate(app.env$object))) {
      ncount <- paste0("nCount_", DefaultAssay(object = isolate(app.env$object)))
      nfeature <- paste0("nFeature_", DefaultAssay(object = isolate(app.env$object)))
      cells.use <- isolate(app.env$object)[[ncount, drop = TRUE]] >= 
        input$num.ncountmin & isolate(app.env$object)[[ncount, 
                                                       drop = TRUE]] <= input$num.ncountmax & isolate(app.env$object)[[nfeature, 
                                                                                                                       drop = TRUE]] >= input$num.nfeaturemin & isolate(app.env$object)[[nfeature, 
                                                                                                                                                                                         drop = TRUE]] <= input$num.nfeaturemax
      if (any(grepl(pattern = mito.pattern, x = rownames(x = isolate(app.env$object))))) {
        cells.use <- cells.use & isolate(app.env$object)[[mt.key, 
                                                          drop = TRUE]] >= input$minmt & isolate(app.env$object)[[mt.key, 
                                                                                                                  drop = TRUE]] <= input$maxmt
      }
      paste(sum(cells.use), "cells remain after current filters")
    }
  })
  output$text.dladt <- renderText(expr = {
    c("imputed.assay <- readRDS('azimuth_impADT.Rds')", "object <- object[, Cells(imputed.assay)]", 
      "object[['impADT']] <- imputed.assay")
  }, sep = "\n")
  output$text.dlumap <- renderText(expr = {
    c("projected.umap <- readRDS('azimuth_umap.Rds')", "object <- object[, Cells(projected.umap)]", 
      "object[['umap.proj']] <- projected.umap")
  }, sep = "\n")
  output$text.dlpred <- renderText(expr = {
    c("predictions <- read.delim('azimuth_pred.tsv', row.names = 1)", 
      "object <- AddMetaData(", "\tobject = object,", "\tmetadata = predictions)")
  }, sep = "\n")
  output$text.dlall <- renderText(expr = {
    c("object <- AddAzimuthResults(object, filename = 'azimuth_results.Rds')")
  }, sep = "\n")
  output$table.qc <- renderTable(expr = {
    if (!is.null(x = isolate(app.env$object))) {
      qc <- paste0(c("nCount_", "nFeature_"), app.env$default.assay)
      tbl <- apply(X = isolate(app.env$object)[[qc]], MARGIN = 2, 
                   FUN = quantile, probs = c(0.0, 0.5, 1.0))
      tbl <- as.data.frame(x = tbl)
      colnames(x = tbl) <- c("Fragments per cell", "Peaks detected per cell")
      if (mt.key %in% colnames(x = isolate(app.env$object)[[]])) {
        tbl[, 3] <- quantile(x = isolate(app.env$object)[[mt.key, 
                                                          drop = TRUE]])
        colnames(x = tbl)[3] <- "Mitochondrial percentage per cell"
      }
      t(x = tbl)
    }
  }, rownames = TRUE)
  output$motifs <- renderDT(expr = {
    print("DIFF EXPRESSION FOR MOTIFS")
    print(head(app.env$chromvar.diff.expr))
    if (!is.null(x = app.env$chromvar.diff.expr[[paste(app.env$chromvar.assay, 
                                                       input$markerclustersgroup.motif, sep = "_")]])) {
      RenderDiffMotifExp(diff.exp = app.env$chromvar.diff.expr[[paste(app.env$chromvar.assay, 
                                                                      input$markerclustersgroup.motif, sep = "_")]], groups.use = input$markerclusters.motif, 
                         n = Inf)
    }
  }, selection = "single", options = list(dom = "t"))
  output$biomarkers <- renderDT(expr = {
    print("DIFF EXPRESSION FOR BIOMARKERS")
    print(head(app.env$diff.expr))
    if (!is.null(x = app.env$diff.expr[[paste(app.env$gene.assay, 
                                              input$markerclustersgroup, sep = "_")]])) {
      RenderDiffExp(diff.exp = app.env$diff.expr[[paste(app.env$gene.assay, 
                                                        input$markerclustersgroup, sep = "_")]], groups.use = input$markerclusters, 
                    n = Inf)
    }
  }, selection = "single", options = list(dom = "t"))
  output$adtbio <- renderDT(expr = {
    if (!is.null(x = app.env$diff.expr[[paste(adt.key, input$markerclustersgroup, 
                                              sep = "_")]])) {
      RenderDiffExp(diff.exp = app.env$diff.expr[[paste(adt.key, 
                                                        input$markerclustersgroup, sep = "_")]], groups.use = input$markerclusters, 
                    n = Inf)
    }
  }, selection = "single", options = list(dom = "t"))
  output$metadata.table <- renderTable(expr = {
    if (!is.null(x = app.env$object)) {
      CategoryTable(object = app.env$object, category.1 = input$metarow, 
                    category.2 = input$metacol, percentage = (input$radio.pct == 
                                                                "Percentage"))
    }
  }, rownames = TRUE)
  output$metadata.heatmap <- renderPlotly({
    table <- CategoryTable(object = app.env$object, category.1 = input$metarow, 
                           category.2 = input$metacol, percentage = (input$radio.pct == 
                                                                       "Percentage"))
    table <- as.matrix(table)
    plot_ly(x = colnames(table), y = rownames(table), z = table, 
            type = "heatmap", height = "1000px")
  })
  output$dlumap <- downloadHandler(filename = paste0(tolower(x = app.title), 
                                                     "_umap.Rds"), content = function(file) {
                                                       if (!is.null(x = app.env$object)) {
                                                         if ("umap.proj" %in% Reductions(object = app.env$object)) {
                                                           saveRDS(object = app.env$object[["umap.proj"]], 
                                                                   file = file)
                                                         }
                                                       }
                                                     })
  output$dladt <- downloadHandler(filename = paste0(tolower(x = app.title), 
                                                    "_impADT.Rds"), content = function(file) {
                                                      if (!is.null(x = app.env$object)) {
                                                        if ("impADT" %in% Assays(object = app.env$object)) {
                                                          saveRDS(object = app.env$object[["impADT"]], 
                                                                  file = file)
                                                        }
                                                      }
                                                    })
  output$dlpred <- downloadHandler(filename = paste0(tolower(x = app.title), 
                                                     "_pred.tsv"), content = function(file) {
                                                       req <- paste0("predicted.", c(app.env$metadataxfer, paste0(app.env$metadataxfer, 
                                                                                                                  ".score")))
                                                       if (resolved(x = app.env$mapping.score)) {
                                                         req <- c(req, "mapping.score")
                                                       }
                                                       if (all(req %in% colnames(x = app.env$object[[]]))) {
                                                         pred.df <- app.env$object[[req]]
                                                         if (resolved(x = app.env$mapping.score)) {
                                                           pred.df$mapping.score <- value(app.env$mapping.score)
                                                         }
                                                         pred.df <- cbind(cell = rownames(x = pred.df), pred.df)
                                                         write.table(x = pred.df, file = file, quote = FALSE, 
                                                                     row.names = FALSE, col.names = TRUE, sep = "\t")
                                                       }
                                                     })
  output$dlall <- downloadHandler(filename = paste0(tolower(x = app.title), 
                                                    "_results.Rds"), content = function(file) {
                                                      results <- list()
                                                      if (!is.null(x = app.env$object)) {
                                                        if ("impADT" %in% Assays(object = app.env$object)) {
                                                          results$impADT <- app.env$object[["impADT"]]
                                                        }
                                                        if ("umap.proj" %in% Reductions(object = app.env$object)) {
                                                          results$umap <- app.env$object[["umap.proj"]]
                                                        }
                                                      }
                                                      req <- paste0("predicted.", c(app.env$metadataxfer, paste0(app.env$metadataxfer, 
                                                                                                                 ".score")))
                                                      if (resolved(x = app.env$mapping.score)) {
                                                        req <- c(req, "mapping.score")
                                                      }
                                                      if (all(req %in% colnames(x = app.env$object[[]]))) {
                                                        pred.df <- app.env$object[[req]]
                                                        if (resolved(x = app.env$mapping.score)) {
                                                          pred.df$mapping.score <- value(app.env$mapping.score)
                                                        }
                                                        pred.df <- cbind(cell = rownames(x = pred.df), pred.df)
                                                        results$pred.df <- pred.df
                                                      }
                                                      saveRDS(results, file = file)
                                                    })
  output$dlscript <- downloadHandler(filename = paste0(tolower(x = app.title), 
                                                       "_analysis.R"), content = function(file) {
                                                         template <- readLines(con = system.file(file.path("resources", 
                                                                                                           "template.R"), package = "Azimuth"))
                                                         template <- paste(template, collapse = "\n")
                                                         e <- new.env()
                                                         e$ref.uri <- getOption(x = "Azimuth.app.refuri", default = getOption(x = "Azimuth.app.reference"))
                                                         e$path <- input$file$name
                                                         e$mito.pattern <- getOption(x = "Azimuth.app.mito", default = "^MT-")
                                                         e$mito.key <- mt.key
                                                         e$ncount.max <- input$num.ncountmax
                                                         e$ncount.min <- input$num.ncountmin
                                                         e$nfeature.max <- input$num.nfeaturemax
                                                         e$nfeature.min <- input$num.nfeaturemin
                                                         e$mito.max <- input$maxmt
                                                         e$mito.min <- input$minmt
                                                         e$sct.ncells <- getOption(x = "Azimuth.sct.ncells")
                                                         e$sct.nfeats <- getOption(x = "Azimuth.sct.nfeats")
                                                         e$ntrees <- getOption(x = "Azimuth.map.ntrees")
                                                         e$ndims <- getOption(x = "Azimuth.map.ndims")
                                                         e$adt.key <- adt.key
                                                         e$do.adt <- do.adt
                                                         e$metadataxfer <- app.env$metadataxfer
                                                         if (length(x = e$metadataxfer == 1)) {
                                                           e$metadataxfer <- paste0("\"", e$metadataxfer, "\"")
                                                         }
                                                         e$plotgene <- getOption(x = "Azimuth.app.default_gene")
                                                         e$plotadt <- getOption(x = "Azimuth.app.default_adt")
                                                         writeLines(text = str_interp(string = template, env = e), 
                                                                    con = file)
                                                       })
  output$refdescriptor <- renderText(expr = eval(expr = HTML(getOption(x = "Azimuth.app.refdescriptor"))))
  output$welcomebox <- renderUI(expr = eval(expr = parse(text = getOption(x = "Azimuth.app.welcomebox"))))
  onclick("panchors_overlap", showModal(modalDialog(title = "Overlap QC", 
                                                    div(paste("In order to conduct bridge integration for ATAC data without uploading a large ", 
                                                              "fragment file, we requantify the ATAC query peaks to match the multiomic bridge ", 
                                                              "based on the overlap between each query peak to a bridge peak and rename the query ", 
                                                              "peak to the bridge peak with highest overlap. The box color corresponds to the following bins: "), 
                                                        tags$ul(list(tags$li(paste0("0% to ", getOption(x = "Azimuth.map.panchorscolors")[1], 
                                                                                    "%: Likely problematic (red)")), tags$li(paste0(getOption(x = "Azimuth.map.panchorscolors")[1], 
                                                                                                                                    "% to ", getOption(x = "Azimuth.map.panchorscolors")[2], 
                                                                                                                                    "%: Possibly problematic (yellow)")), tags$li(paste0(getOption(x = "Azimuth.map.panchorscolors")[2], 
                                                                                                                                                                                         "% to 100%: Likely successful (green)")))), tags$h4("Caveats"), 
                                                        paste0("A high percentage of overlap is expected if the query ATAC data and bridge ATAC data ", 
                                                               "were processed with the same versions of Cell Ranger and means that there will ", 
                                                               "likely be little loss of information by using this overlap renaming process. ", 
                                                               "The mapping can still be sucessesful if this value has a low percentage, but downstream gene expression", 
                                                               "calculations may be innacurate as this again uses another overlap process to requantify peaks to genes.")))))
  onclick("panchors_popup", showModal(modalDialog(title = "Anchor QC", 
                                                  div(paste("The Azimuth reference-mapping procedure first identifies a set of 'anchors', ", 
                                                            "or pairwise correspondences between cells predicted to be in a similar biological state, ", 
                                                            "between query and reference datasets. Here we report the percentage of query cells ", 
                                                            "participating in an anchor correspondence. The box color corresponds to the following bins: "), 
                                                      tags$ul(list(tags$li(paste0("0% to ", getOption(x = "Azimuth.map.panchorscolors")[1], 
                                                                                  "%: Likely problematic (red)")), tags$li(paste0(getOption(x = "Azimuth.map.panchorscolors")[1], 
                                                                                                                                  "% to ", getOption(x = "Azimuth.map.panchorscolors")[2], 
                                                                                                                                  "%: Possibly problematic (yellow)")), tags$li(paste0(getOption(x = "Azimuth.map.panchorscolors")[2], 
                                                                                                                                                                                       "% to 100%: Likely successful (green)")))), tags$h4("Caveats"), 
                                                      paste0("If the query dataset consists of a homogeneous group of cells, or if the ", 
                                                             "query dataset contains cells from multiple batches (which would be corrected ", 
                                                             "by Azimuth), this metric may return a low value even in cases where mapping is ", 
                                                             "successful. Users in these cases should check results carefully. In particular, ", 
                                                             "we encourage users to verify identified differentially expressed marker genes for annotated cell types.")))))
  onclick("mappingqcstat_popup", showModal(modalDialog(title = "Cluster Preservation", 
                                                       div(tags$h4("Overview"), paste0("For each query dataset, we downsample to at most 5,000 cells, and perform an ", 
                                                                                       "unsupervised clustering. This score reflects the preservation of the unsupervised ", 
                                                                                       "cluster structure, and is based on the entropy of unsupervised cluster labels in ", 
                                                                                       "each query cell's local neighborhood after mapping. Scores are scaled from 0 (poor) to 5 (best)"), 
                                                           tags$ul(list(tags$li(paste0("0 to ", getOption(x = "Azimuth.map.postmapqccolors")[1], 
                                                                                       ": Likely problematic (red)")), tags$li(paste0(getOption(x = "Azimuth.map.postmapqccolors")[1], 
                                                                                                                                      " to ", getOption(x = "Azimuth.map.postmapqccolors")[2], 
                                                                                                                                      ": Possibly problematic (yellow)")), tags$li(paste0(getOption(x = "Azimuth.map.postmapqccolors")[2], 
                                                                                                                                                                                          " to 5: Likely successful (green)")))), tags$h4("Caveats"), 
                                                           paste0("This metric relies on the unsupervised clustering representing corresponding to ", 
                                                                  "biologically distinct cell states. If the query dataset consists of a homogeneous ", 
                                                                  "group of cells, or if the query dataset contains cells from multiple batches ", 
                                                                  "(which would be corrected by Azimuth), this metric may return a low value even ", 
                                                                  "in cases where mapping is successful. Users in these cases should check results ", 
                                                                  "carefully. In particular, we encourage users to verify identified differentially ", 
                                                                  "expressed marker genes for annotated cell types."), 
                                                       ))))
}



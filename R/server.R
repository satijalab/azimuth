#' @include zzz.R
#'
NULL

#' Server function for the mapping app
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
  mt.key <- 'percent.mt'
  mito.pattern <- getOption(x = 'Azimuth.app.mito', default = '^MT-')
  do.adt <- isTRUE(x = getOption(x = 'Azimuth.app.do_adt', default = TRUE))
  adt.key <- 'impADT'
  n.trees <- getOption(x = "Azimuth.map.ntrees")
  app.env <- reactiveValues(
    adt.features = character(length = 0L),
    anchors = NULL,
    clusterpreservationqc = NULL,
    demo = FALSE,
    demo.inputs = NULL,
    demo.tracker = NULL,
    demo.files = NULL,
    default.assay = NULL,
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
    disable = FALSE
  )
  react.env <- reactiveValues(
    no = FALSE,
    anchors = FALSE,
    biomarkers = FALSE,
    cluster.score = FALSE,
    features = FALSE,
    map = FALSE,
    markers = FALSE,
    metadata = FALSE,
    mt = NULL,
    xferopts = FALSE,
    path = NULL,
    progress = NULL,
    plot.qc = FALSE,
    qc = FALSE,
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
  ResetEnv <- function() {
    print('resetting...')
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
    disable(id = 'map')
    hide(selector = '.rowhide')
  }
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
  if (!inherits(x = demos, what = "data.frame") & !is.null(x = demos)) {
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
  withProgress(
    message = "Loading reference",
    expr = {
      disable(id = 'file')
      for (i in 1:nrow(x = demos)) {
        disable(id = paste0("triggerdemo", i))
      }
      setProgress(value = 0)
      refs <- LoadReference(
        path = getOption(
          x = 'Azimuth.app.reference',
          default = stop(safeError(error = "No reference provided"))
        )
      )
      setProgress(value = 1)
      enable(id = 'file')
      for (i in 1:nrow(x = demos)) {
        enable(id = paste0("triggerdemo", i))
      }
    }
  )
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
  # Load the data an prepare for QC
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
    eventExpr = react.env$path,
    handlerExpr = {
      if (!is.null(x = react.env$path) && nchar(x = react.env$path)) {
        withProgress(
          message = 'Reading Input',
          expr = {
            setProgress(value = 0)
            tryCatch(
              expr = {
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
  )
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
        updateNumericInput(
          session = session,
          inputId = 'num.ncountmin',
          label = paste('min', ncount),
          value = ncount.val[1],
          min = ncount.val[1],
          max = ncount.val[2]
        )
        updateNumericInput(
          session = session,
          inputId = 'num.ncountmax',
          label = paste('max', ncount),
          value = ncount.val[2],
          min = ncount.val[1],
          max = ncount.val[2]
        )
        nfeature.val <- range(app.env$object[[nfeature, drop = TRUE]])
        nfeature.val <- c(
          floor(x = min(nfeature.val)),
          ceiling(x = max(nfeature.val))
        )
        updateNumericInput(
          session = session,
          inputId = 'num.nfeaturemin',
          label = paste('min', nfeature),
          value = nfeature.val[1],
          min = nfeature.val[1],
          max = nfeature.val[2]
        )
        updateNumericInput(
          session = session,
          inputId = 'num.nfeaturemax',
          label = paste('max', nfeature),
          value = nfeature.val[2],
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
          updateNumericInput(
            session = session,
            inputId = 'minmt',
            label = paste('min', mt.key),
            value = mito.val[1],
            min = mito.val[1],
            max = mito.val[2]
          )
          updateNumericInput(
            session = session,
            inputId = 'maxmt',
            label = paste('max', mt.key),
            value = mito.val[2],
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
          for (demo.id in app.env$demo.inputs) {
            enable(id = demo.id)
          }
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
      for (demo.id in app.env$demo.inputs) {
        disable(id = demo.id)
      }
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
        react.env$sctransform <- TRUE
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
          for (demo.id in app.env$demo.inputs) {
            enable(id = demo.id)
          }
          gc(verbose = FALSE)
        } else {
          query.unique <- length(x = unique(x = slot(object = app.env$anchors, name = "anchors")[, "cell2"]))
          percent.anchors <- round(x = query.unique / ncol(x = app.env$object) * 100, digits = 2)
          if (percent.anchors <  getOption(x = "Azimuth.map.panchorscolors")[1]) {
            output$valuebox_panchors <- renderValueBox(expr = {
              valueBox(
                value = paste0(percent.anchors, "%"),
                subtitle = "% of query cells with anchors",
                color = 'red',
                icon = icon(name = 'times')
              )
            })
          } else if (percent.anchors <  getOption(x = "Azimuth.map.panchorscolors")[2]) {
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
          for (demo.id in app.env$demo.inputs) {
            enable(id = demo.id)
          }
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
    eventExpr = react.env$cluster.score,
    handlerExpr = {
      if (isTRUE(react.env$cluster.score)) {
        # post mapping QC
        qc.stat <- round(
          x = ClusterPreservationScore(
            query = app.env$object,
            ds.amount = getOption(x = "Azimuth.map.postmapqcds")
          ),
          digits = 2
        )
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
        qc.stat <- round(
          x = ClusterPreservationScore(
            query = app.env$object,
            ds.amount = getOption(x = "Azimuth.map.postmapqcds")
          ),
          digits = 2
        )
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
        if (Sys.getenv("RSTUDIO") == "1") {
          plan("sequential")
        }
        # reduce size of object in anchorset
        app.env$anchors@object.list[[1]] <- DietSeurat(
          object = app.env$anchors@object.list[[1]]
        )
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
        react.env$progress$set(value = 0.8, message = 'Running UMAP transform')
        app.env$object[["query_ref.nn"]] <- FindNeighbors(
          object = Embeddings(refs$map[["refDR"]]),
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
        gc(verbose = FALSE)
        app.env$messages <- c(
          app.env$messages,
          paste(ncol(x = app.env$object), "cells mapped")
        )
        react.env$biomarkers <- TRUE
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
        for (i in app.env$metadataxfer[!app.env$singlepred]) {
          app.env$diff.expr[[paste(app.env$default.assay, i, sep = "_")]] <- wilcoxauc(
            X = app.env$object,
            group_by = paste0("predicted.", i),
            assay = 'data',
            seurat_assay = app.env$default.assay
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
        react.env$progress$close()
        enable(id = 'file')
        for (demo.id in app.env$demo.inputs) {
          enable(id = demo.id)
        }
        react.env$metadata <- TRUE
        react.env$biomarkers <- FALSE
      }
    }
  )
  # Update input controls
  observeEvent(
    eventExpr = react.env$metadata,
    handlerExpr = {
      if (isTRUE(x = react.env$metadata)) {
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
        for (id in c('metarow', 'metacol', 'metagroup')) {
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
          inputId = 'metacolor.ref',
          choices = c(grep(pattern = '^predicted.', x = app.env$metadataxfer, value = TRUE), # re-ordering not working...
                      grep(pattern = '^predicted.', x = app.env$metadataxfer, value = TRUE, invert = TRUE)),
          selected = app.env$default.metadata,
          server = TRUE,
          options = selectize.opts[-which(x = names(x = selectize.opts) == 'maxItems')]
        )
        react.env$features <- TRUE
        react.env$metadata <- FALSE
      }
    }
  )
  observeEvent(
    eventExpr = react.env$features,
    handlerExpr = {
      if (isTRUE(x = react.env$features)) {
        app.env$default.feature <- ifelse(
          test = getOption(x = 'Azimuth.app.default_gene') %in% rownames(x = app.env$object),
          yes = getOption(x = 'Azimuth.app.default_gene'),
          no = VariableFeatures(object = app.env$object)[1]
        )
        app.env$features <- unique(x = c(
          FilterFeatures(
            features = VariableFeatures(object = app.env$object)[1:selectize.opts$maxOptions]
          ),
          FilterFeatures(features = rownames(x = app.env$object))
        ))
        updateSelectizeInput(
          session = session,
          inputId = 'feature',
          label = 'Feature',
          choices = app.env$features,
          selected = app.env$default.feature,
          server = TRUE,
          options = selectize.opts
        )
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
        react.env$markers <- TRUE
        react.env$features <- FALSE
      }
    }
  )

  observeEvent(
    eventExpr = react.env$markers,
    handlerExpr = {
      if (isTRUE(x = react.env$markers)) {
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

        react.env$markers <- FALSE
        app.env$disable <- FALSE
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
        app.env$feature <- ifelse(
          test = input$feature %in% rownames(x = app.env$object),
          yes = paste0(
            Key(object = app.env$object[["refAssay"]]),
            input$feature
          ),
          no = input$feature
        )
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
          diff.exp = app.env$diff.expr[[paste(app.env$default.assay, input$markerclustersgroup, sep = "_")]],
          groups.use = input$markerclusters,
          n = Inf
        ))
        tables.clear <- list(adt.proxy, rna.proxy)[c(TRUE, !table.check)]
        for (tab in tables.clear) {
          selectRows(proxy = tab, selected = NULL)
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
      avail <- c(
        paste0(
          Key(object = app.env$object[["refAssay"]]),
          rownames(x = app.env$object)
        ),
        colnames(x = app.env$object[[]])
      )
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
      palettes <- list(
        c("lightgrey", "blue"),
        c('lightgrey', 'darkred')
      )
      names(x = palettes) <- c(
        Key(object = app.env$object[["refAssay"]]),
        'md_'
      )
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
        "object <- AddAzimuthResults(object, azimuth_results = readRDS('azimuth_results.Rds'))"
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
        colnames(x = tbl) <- c('nUMI per cell', 'Genes detected per cell')
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
      if (!is.null(x = app.env$diff.expr[[paste(app.env$default.assay, input$markerclustersgroup, sep ="_")]])) {
        RenderDiffExp(
          diff.exp =  app.env$diff.expr[[paste(app.env$default.assay, input$markerclustersgroup, sep ="_")]],
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
    filename = paste0(tolower(x = app.title), '_all.Rds'),
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

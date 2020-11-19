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
#' theme_void xlab
#' @importFrom googlesheets4 gs4_auth gs4_get sheet_append
#' @importFrom methods slot slot<-
#' @importFrom presto wilcoxauc
#' @importFrom Seurat AddMetaData Assays Cells DimPlot DefaultAssay Embeddings
#' FeaturePlot FindNeighbors FindTransferAnchors GetAssayData Idents Idents<-
#' IntegrateEmbeddings Key MappingScore NoLegend PercentageFeatureSet
#' RenameCells Reductions RunUMAP Tool TransferData SetAssayData SCTransform
#' VariableFeatures VlnPlot
#' @importFrom shiny downloadHandler observeEvent isolate Progress
#' reactiveValues renderPlot renderTable renderText removeUI setProgress
#' safeError updateNumericInput updateSelectizeInput withProgress
#' @importFrom shinydashboard menuItem renderMenu renderValueBox
#' sidebarMenu valueBox
#' @importFrom shinyjs addClass enable disable hide removeClass show
#' @importFrom stringr str_interp
#' @importFrom patchwork wrap_plots
#' @importFrom stats na.omit quantile
#' @importFrom utils write.table
#'
#' @keywords internal
#'
AzimuthServer <- function(input, output, session) {
  # hide demo dataset button if required
  if (is.null(x = getOption(x = 'Azimuth.app.demodataset'))) {
    hide(id="triggerdemo")
  }
  mt.key <- 'percent.mt'
  mito.pattern <- getOption(x = 'Azimuth.app.mito', default = '^MT-')
  do.adt <- isTRUE(x = getOption(x = 'Azimuth.app.do_adt', default = TRUE))
  adt.key <- 'impADT'
  n.trees <- getOption(x = "Azimuth.map.ntrees")
  app.env <- reactiveValues(
    adt.features = character(length = 0L),
    anchors = NULL,
    default.assay = NULL,
    default.feature = NULL,
    diff.exp = list(),
    feature = '',
    features = character(length = 0L),
    mapping.score = NULL,
    messages = 'Upload a file',
    object = NULL,
    scorefeatures = character(length = 0L)
  )
  react.env <- reactiveValues(
    no = FALSE,
    anchors = FALSE,
    biomarkers = FALSE,
    features = FALSE,
    map = FALSE,
    markers = FALSE,
    metadata = FALSE,
    mt = FALSE,
    path = NULL,
    pbcor = FALSE,
    progress = NULL,
    qc = FALSE,
    score = FALSE,
    sctransform = FALSE,
    start = numeric(length = 0L),
    transform = FALSE
  )
  if (isTRUE(x = do.adt)) {
    output$imputedlabel <- renderUI(expr = h3('Imputed protein biomarkers'))
  } else {
    for (id in c('imputedinput', 'imputedtable')) {
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
    app.env$messages <- NULL
    output$valubox.upload <- NULL
    output$valuebox.preproc <- NULL
    output$valuebox.mapped <- NULL
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
  if (logging) {
    tryCatch(
      expr = {
        gs4_auth(
          email = getOption(x = "Azimuth.app.googletokenemail"),
          cache = getOption(x = "Azimuth.app.googletoken")
        )
        googlesheet <- gs4_get(ss = getOption(x = "Azimuth.app.googlesheet"))
        app_start_time <- Sys.time()
        onStop(
          fun = function() {
            try(expr = sheet_append(
              ss = googlesheet,
              data = data.frame(
                "SESSIONLENGTH",
                Sys.info()[["nodename"]],
                as.numeric(x = Sys.time() - app_start_time, units = "mins")
              )
            ))
          }
        )
      },
      error = function(e) {
        googlesheet <- NULL
      }
    )
  } else {
    googlesheet <- NULL
  }
  if (!is.null(x = googlesheet)) {
    try(expr = sheet_append(
      ss = googlesheet,
      data = data.frame(
        "STARTUPTIME",
        Sys.info()[["nodename"]],
        Sys.time()
      )
    ))
  }
  withProgress(
    message = "Loading reference",
    expr = {
      setProgress(value = 0)
      refs <- LoadReference(
        path = getOption(
          x = 'Azimuth.app.reference',
          default = stop(safeError(error = "No reference provided"))
        )
      )
      setProgress(value = 1)
    }
  )
  set.seed(seed = getOption(x = "Azimuth.app.plotseed"))
  plotlevels <- sample(x = levels(as.factor(refs$plot$id)), size = length(levels(as.factor(refs$plot$id))))
  # React to events
  # Load the data an prepare for QC
  observeEvent(
    eventExpr = input$file,
    handlerExpr = {
      if (nchar(x = input$file$datapath)) {
        react.env$path <- input$file$datapath
      }
    }
  )
  observeEvent(
    eventExpr = input$triggerdemo,
    handlerExpr = react.env$path <- getOption(x = 'Azimuth.app.demodataset')
  )
  observeEvent(
    eventExpr = react.env$path,
    handlerExpr = {
      if (!is.null(x = react.env$path) && nchar(x = react.env$path)) {
        ResetEnv()
        withProgress(
          message = 'Reading Input',
          expr = {
            setProgress(value = 0)
            tryCatch(
              expr = {
                app.env$object <- LoadFileInput(path = react.env$path)
                app.env$object$query <- 'query'
                Idents(object = app.env$object) <- 'query'
              },
              error = function(e) {
                app.env$messages <- e$message
                safeError(error = e$message)
              }
            )
            setProgress(value = 1)
          }
        )
        app.env$default.assay <- DefaultAssay(object = app.env$object)
        react.env$mt <- any(grepl(
          pattern = mito.pattern,
          x = rownames(x = app.env$object)
        ))
        if (isFALSE(x = react.env$mt)) {
          removeUI(selector = '#pctmt', immediate = TRUE)
        }
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
        react.env$qc <- !any(reject)
        react.env$path <- NULL
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
          label = paste('min', ncount),
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
          label = paste('min', nfeature),
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
        app.env$messages <- paste(ncellsupload, 'cells uploaded')
        if (ncellsupload < getOption(x = 'Azimuth.map.ncells')) {
          output$valuebox.upload <- renderValueBox(expr = {
            valueBox(
              value = ncellsupload,
              subtitle = 'cells uploaded',
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
          enable(id = 'map')
          if (!is.null(x = googlesheet)) {
            try(
              expr = sheet_append(
                ss = googlesheet,
                data = data.frame(
                  'CELLSUPLOAD',
                  Sys.info()[['nodename']],
                  ncellsupload
                )
              ),
              silent = TRUE
            )
          }
        }
        if (!is.null(x = react.env$progress)) {
          react.env$progress$close()
          react.env$progress <- NULL
        }
        react.env$qc <- FALSE
      }
    }
  )
  # Filter and process the data
  observeEvent(
    eventExpr = input$map,
    handlerExpr = {
      react.env$start <- Sys.time()
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
      # not enough cells available after filtering: reset filter elements
      if (ncellspreproc < getOption(x = "Azimuth.map.ncells")) {
        output$valuebox.preproc <- renderValueBox(expr = valueBox(
          value = ncellspreproc,
          subtitle = "cells after filtering",
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
              Sys.info()[["nodename"]],
              ncellspreproc
            )
          ))
        }
        app.env$object <- app.env$object[, cells.use]
        react.env$pbcor <- TRUE
      }
    }
  )
  observeEvent(
    eventExpr = react.env$pbcor,
    handlerExpr = {
      if (isTRUE(x = react.env$pbcor)) {
        react.env$progress$set(
          value = 0.1,
          message = 'Running pseudobulk correlation test'
        )
        pbcor <- PBCorTest(
          object = app.env$object,
          ref = refs$avg
        )
        if (!is.null(googlesheet)) {
          try(sheet_append(
            ss = googlesheet,
            data = data.frame(
              "PBCOR",
              Sys.info()[["nodename"]],
              pbcor[["cor.res"]]
            )
          ))
        }
        if (pbcor[["cor.res"]] < getOption(x = 'Azimuth.map.pbcorthresh')) {
          output$valuebox.mapped <- renderValueBox(expr = {
            valueBox(
              value = 'Failure',
              subtitle = 'Query is too dissimilar',
              icon = icon(name = 'times'),
              color = 'red',
              width = 6
            )
          })
          output$plot.pbcor <- renderPlot(expr = pbcor[['plot']])
          show(selector = ".rowhide")
          app.env$object <- NULL
          react.env$progress$close()
          gc(verbose = FALSE)
        } else {
          react.env$sctransform <- TRUE
        }
        react.env$pbcor <- FALSE
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
              method = "glmGamPoi",
              ncells = getOption(x = 'Azimuth.sct.ncells'),
              n_genes = getOption(x = 'Azimuth.sct.nfeats'),
              do.correct.umi = FALSE,
              do.scale = FALSE,
              do.center = TRUE
            ))
          },
          error = function(e) {
            app.env$object <- suppressWarnings(expr = SCTransform(
              object = app.env$object,
              residual.features = rownames(x = refs$map),
              method = "poisson",
              ncells = getOption(x = 'Azimuth.sct.ncells'),
              n_genes = getOption(x = 'Azimuth.sct.nfeats'),
              do.correct.umi = FALSE,
              do.scale = FALSE,
              do.center = TRUE
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
          reference.neighbors = "spca.annoy.neighbors",
          reference.assay = "SCT",
          query.assay = 'SCT',
          reference.reduction = 'spca',
          normalization.method = 'SCT',
          features = intersect(
            x = rownames(x = refs$map),
            y = VariableFeatures(object = app.env$object)
          ),
          dims = 1:50,
          n.trees = n.trees,
          verbose = TRUE,
          mapping.score.k = 100
        )
        # TODO: fail if not enough anchors (Azimuth.map.anchors)
        react.env$map <- TRUE
        react.env$anchors <- FALSE
      }
    }
  )
  observeEvent(
    eventExpr = react.env$map,
    handlerExpr = {
      if (isTRUE(x = react.env$map)) {
        react.env$progress$set(value = 0.5, message = 'Mapping cells')
        app.env$object <- TransferData(
          reference = refs$map,
          query = app.env$object,
          dims = 1:50,
          anchorset = app.env$anchors,
          refdata = list(
            id = Idents(object = refs$map),
            impADT = GetAssayData(
              object = refs$map[['ADT']],
              slot = 'data'
            )),
          n.trees = n.trees,
          store.weights = TRUE
        )
        app.env$object <- IntegrateEmbeddings(
          anchorset = app.env$anchors,
          reference = refs$map,
          query = app.env$object,
          reductions = "pcaproject",
          reuse.weights.matrix = TRUE
        )
        react.env$score <- TRUE
        react.env$map <- FALSE
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
        spca <- subset(
          x = app.env$anchors@object.list[[1]][["pcaproject.l2"]],
          cells = paste0(Cells(x = app.env$object), "_query")
        )
        spca <- RenameCells(
          object = spca,
          new.names = Cells(x = app.env$object)
        )
        spca.ref <- subset(
          x = app.env$anchors@object.list[[1]][["pcaproject.l2"]],
          cells = paste0(Cells(x = refs$map), "_reference")
        )
        spca.ref <- RenameCells(
          object = spca.ref,
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
        app.env$mapping.score <- future(
          expr = {
            MappingScore(
              anchors = app.env$anchors@anchors,
              combined.object = app.env$anchors@object.list[[1]],
              query.neighbors =  slot(object = app.env$anchors, name = "neighbors")[["query.neighbors"]],
              query.weights = Tool(object = app.env$object, slot = "TransferData")$weights.matrix,
              query.embeddings = Embeddings(object = spca),
              ref.embeddings = Embeddings(object = spca.ref),
              nn.method = "annoy",
              n.trees = n.trees
            )
          }
        )
        app.env$object <- AddMetaData(
          object = app.env$object,
          metadata = rep(x = 0, times = ncol(x = app.env$object)),
          col.name = "mapping.score"
        )
        app.env$anchors <- NULL
        rm(spca)
        gc(verbose = FALSE)
        react.env$transform <- TRUE
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
          object = Embeddings(refs$map[["spca"]]),
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
          reduction.model = refs$map[['jumap']],
          reduction.key = 'UMAP_'
        )
        app.env$object <- SetAssayData(
          object = app.env$object,
          assay = 'SCT',
          slot = 'scale.data',
          new.data = new(Class = 'matrix')
        )
        gc(verbose = FALSE)
        app.env$messages <- c(
          app.env$messages,
          paste(ncol(x = app.env$object), "cells mapped")
        )
        output$valuebox.mapped <- renderValueBox(expr = {
          valueBox(
            value = 'Success',
            subtitle = 'Mapping complete',
            icon = icon(name = 'check'),
            color = 'green'
          )
        })
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
        app.env$diff.expr[[app.env$default.assay]] <- wilcoxauc(
          X = app.env$object,
          group_by = 'predicted.id',
          assay = 'data',
          seurat_assay = app.env$default.assay
        )
        if (isTRUE(x = do.adt)) {
          app.env$diff.expr[[adt.key]] <- wilcoxauc(
            X = app.env$object,
            group_by = 'predicted.id',
            assay = 'data',
            seurat_assay = adt.key
          )
        }
        # Finalize the log
        time.fmt <- FormatDiffTime(dt = difftime(
          time1 = Sys.time(),
          time2 = react.env$start,
          units = 'secs'
        ))
        app.env$messages <- c(
          app.env$messages,
          time.fmt
        )
        if (!is.null(x = googlesheet)) {
          try(expr = sheet_append(
            ss = googlesheet,
            data = data.frame(
              "MAPPINGTIME",
              Sys.info()[["nodename"]],
              as.numeric(x = maptime.diff, units="secs")
            )
          ))
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
        metadata.discrete <- sort(x = c(
          "predicted.id",
          PlottableMetadataNames(
            object = app.env$object,
            min.levels = 1,
            max.levels = 50
          )
        ))
        for (id in c('metacolor', 'metarow', 'metacol', 'metagroup')) {
          updateSelectizeInput(
            session = session,
            inputId = id,
            choices = metadata.discrete,
            selected = 'predicted.id',
            server = TRUE,
            options = selectize.opts
          )
        }
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
        metadata.cont <- sort(x = c(
          metadata.cont,
          rownames(x = app.env$object[["prediction.score.id"]])
        ))
        app.env$scorefeatures <- metadata.cont
        updateSelectizeInput(
          session = session,
          inputId = 'scorefature',
          choices = app.env$scorefeatures,
          selected = '',
          server = TRUE,
          options = selectize.opts
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
        app.env$features <- FilterFeatures(
          features = rownames(x = app.env$object)
        )
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
          app.env$adt.features <- sort(x = FilterFeatures(features = rownames(
            x = app.env$object[[adt.key]]
          )))
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
          x = table(app.env$object$predicted.id) > getOption(x = 'Azimuth.de.mincells')
        ))
        allowed.clusters <- factor(
          x = allowed.clusters,
          levels = unique(x = as.vector(x = app.env$object$predicted.id))
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
        react.env$markers <- FALSE
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
            Key(object = app.env$object[["SCT"]]),
            input$feature
          ),
          no = input$feature
        )
        for (f in c('adtfeature', 'scorefeature')) {
          updateSelectizeInput(
            session = session,
            inputId = f,
            choices = list(
              'adtfeature' = app.env$adt.features,
              'scorefeature' = app.env$scorefeatures
            )[[f]],
            selected = '',
            server = TRUE,
            options = selectize.opts
          )
        }
        table.check <- input$feature %in% rownames(x = RenderDiffExp(
          diff.exp = app.env$diff.expr[[app.env$default.assay]],
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
        for (f in c('feature', 'scorefeature')) {
          updateSelectizeInput(
            session = session,
            inputId = f,
            choices = list(
              'feature' = app.env$features,
              'scorefeature' = app.env$scorefeatures
            )[[f]],
            selected = '',
            server = TRUE,
            options = selectize.opts
          )
        }
        table.check <- input$adtfeature %in% rownames(x = RenderDiffExp(
          diff.exp = app.env$diff.expr[[adt.key]],
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
  observeEvent( # Prediction score
    eventExpr = input$scorefeature,
    handlerExpr = {
      if (nchar(x = input$scorefeature)) {
        if (input$scorefeature == "mapping.score") {
          if (resolved(x = app.env$mapping.score)) {
            app.env$object$mapping.score <- value(app.env$mapping.score)
          }
        }
        app.env$feature <- input$scorefeature
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
  observeEvent( # Select from biomarkers table
    eventExpr = input$biomarkers_rows_selected,
    handlerExpr = {
      if (length(x = input$biomarkers_rows_selected)) {
        updateSelectizeInput(
          session = session,
          inputId = 'feature',
          choices = app.env$features,
          selected = rownames(x = RenderDiffExp(
            diff.exp = app.env$diff.expr[[app.env$default.assay]],
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
            diff.exp = app.env$diff.expr[[adt.key]],
            groups.use = input$markerclusters,
            n = Inf
          ))[input$adtbio_rows_selected],
          server = TRUE,
          options = selectize.opts
        )
      }
    }
  )
  # Plots
  output$plot.qc <- renderPlot(expr = {
    if (!is.null(x = isolate(app.env$object))) {
      qc <- paste0(c('nCount_', 'nFeature_'), app.env$default.assay)
      if (isTRUE(x = react.env$mt)) {
        qc <- c(qc, mt.key)
      }
      vlnlist <- VlnPlot(
        object = isolate(app.env$object),
        features = qc,
        group.by = 'query',
        combine = FALSE,
        pt.size = ifelse(
          test = input$check.qcpoints,
          yes = 0,
          no = Seurat:::AutoPointSize(data = isolate(app.env$object))
        ),
        log = input$check.qcscale
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
          ymin = ifelse(test = input$check.qcscale, yes = 0, no = -Inf),
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
          ymin = ifelse(test = input$check.qcscale, yes = 0, no = -Inf),
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
            ymin = ifelse(test = input$check.qcscale, yes = 0, no = -Inf),
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
  output$refdim <- renderPlot(expr = {
    DimPlot(
      object = refs$plot,
      label = input$labels,
      group.by = "id",
      repel = TRUE
    ) +
      scale_colour_hue(
        limits = plotlevels,
        breaks = sort(x = levels(x = as.factor(x = refs$plot$id))),
        drop = FALSE
      )
  })
  output$objdim <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      if (length(x = Reductions(object = app.env$object))) {
        if (input$metacolor == "predicted.id") {
          DimPlot(
            object = app.env$object,
            group.by = "predicted.id",
            label = input$labels,
            repel = TRUE,
            reduction = "umap.proj"
          ) + scale_colour_hue(
            limits = plotlevels,
            breaks = sort(x = levels(x = as.factor(x = refs$plot$id))),
            drop = FALSE
          )
        } else {
          DimPlot(
            object = app.env$object,
            group.by = input$metacolor,
            label = input$labels,
            repel = TRUE,
            reduction = "umap.proj"
          )
        }
      }
    }
  })
  output$evln <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      avail <- c(
        paste0(
          Key(object = app.env$object[["SCT"]]),
          rownames(x = app.env$object)
        ),
        paste0(
          Key(object = app.env$object[[adt.key]]),
          rownames(x = app.env$object[[adt.key]])
        ),
        colnames(x = app.env$object[[]]),
        rownames(x = app.env$object[["prediction.score.id"]])
      )
      if (app.env$feature %in% avail) {
        if (app.env$feature == "mapping.score" && !resolved(x = app.env$mapping.score)) {
          ggplot() +
            annotate("text", x = 4, y = 25, size=8, label = "Mapping score still computing ... ") +
            theme_void()
        } else {
          title <- ifelse(
            test = grepl(pattern = '^sct_', x = app.env$feature),
            yes = gsub(pattern = '^sct_', replacement = '', x = app.env$feature),
            no = app.env$feature
          )
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
        c('lightgrey', 'darkgreen'),
        c('lightgrey', 'darkred')
      )
      names(x = palettes) <- c(
        Key(object = app.env$object[["SCT"]]),
        Key(object = app.env$object[[adt.key]]),
        'md_'
      )
      md <- c(
        colnames(x = app.env$object[[]]),
        rownames(x = app.env$object[['prediction.score.id']])
      )
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
            test = grepl(pattern = '^sct_', x = app.env$feature),
            yes = gsub(pattern = '^sct_', replacement = '', x = app.env$feature),
            no = app.env$feature
          )
          suppressWarnings(expr = FeaturePlot(
            object = app.env$object,
            features = app.env$feature,
            cols = pal.use,
            reduction = "umap.proj"
          )) + ggtitle(label = title)
        }
      }
    }
  })
  # Messages
  output$message <- renderUI(expr = {
    p(HTML(text = paste(app.env$messages, collapse = "<br />")))
  })
  output$containerid <- renderText(c("debug ID: ", Sys.info()[["nodename"]]))
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
      if (!is.null(x = app.env$diff.expr[[app.env$default.assay]])) {
        RenderDiffExp(
          diff.exp = app.env$diff.expr[[app.env$default.assay]],
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
      if (!is.null(x = app.env$diff.expr[[adt.key]])) {
        RenderDiffExp(
          diff.exp = app.env$diff.expr[[adt.key]],
          groups.use = input$markerclusters,
          n = Inf
        )
      }
    },
    selection = 'single',
    options = list(dom = 't')
  )
  output$table.metadata <- renderTable(
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
      req <- c('predicted.id', 'predicted.id.score')
      if (resolved(x = app.env$mapping.score)) {
        req <- c(req, 'mapping.score')
      }
      if (all(req %in% colnames(x = app.env$object[[]]))) {
        pred.df <- data.frame(
          cell = colnames(x = app.env$object),
          predicted.id = app.env$object$predicted.id,
          predicted.score = app.env$object$predicted.id.score,
          stringsAsFactors = FALSE
        )
        if (resolved(x = app.env$mapping.score)) {
          pred.df$mapping.score <- value(app.env$mapping.score)
        }
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
  output$dlscript <- downloadHandler(
    filename = paste0(tolower(x = app.title), '_analysis.R'),
    content = function(file) {
      template <- readLines(con = system.file(
        file.path('resources', 'template.R'),
        package = 'Azimuth'
      ))
      template <- paste(template, collapse = '\n')
      e <- new.env()
      e$ref.uri <- 'https://seurat.nygenome.org/references/pbmc/'
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
      e$plotgene <- getOption(x = 'Azimuth.app.default_gene')
      e$plotadt <- getOption(x = 'Azimuth.app.default_adt')
      writeLines(text = str_interp(string = template, env = e), con = file)
    }
  )
  # render UI elements that depend on arguments
  output$welcomebox <- renderUI(
    expr = eval(expr = parse(text = getOption(x = "Azimuth.app.welcomebox")))
  )
}

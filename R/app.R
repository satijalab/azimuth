#' @include zzz.R
#' @include seurat.R
#' @import V8
#' @importFrom DT DTOutput
#' @importFrom htmltools tagList h4 hr h3 tags HTML p div
#' @importFrom shinyjs useShinyjs extendShinyjs disabled
#' @importFrom shiny fluidPage sidebarLayout sidebarPanel fileInput sliderInput
#' actionButton selectInput downloadButton mainPanel tabsetPanel tabPanel
#' plotOutput tableOutput verbatimTextOutput numericInput icon fluidRow
#' updateNumericInput radioButtons textOutput htmlOutput column
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#' dashboardBody menuItem tabItems tabItem sidebarMenu box valueBoxOutput
#' sidebarMenuOutput
#' @importFrom shinyBS bsPopover bsButton
#'
NULL

app.title <- 'Azimuth'

ui <- tagList(
  useShinyjs(),
  tags$style(type = "text/css", "#message {padding: 0px 10px;}"),
  tags$style(type = "text/css", "#containerid {position:fixed;bottom:0;right:0;left:0;padding: 0px 10px;}"),
  dashboardPage(
  dashboardHeader(title = app.title),
  dashboardSidebar(
    fileInput(
      inputId = "file",
      label = p(
        "File Upload",
        tags$style(type = "text/css", "#q1 {display: inline-block; vertical-align: middle;}"),
        bsButton("q1", label = "", icon = icon("question"), style = "info", size = "extra-small")
      ),
      accept = c('.h5', '.h5ad', '.h5seurat', '.rds')
    ),
    bsPopover(
      id = "q1",
      title = "Supported file types",
      content = "10x Genomics H5; Seurat object (RDS); H5AD; H5Seurat; Matrix/matrix/data.frame (RDS)",
      placement = "right",
      trigger = "focus",
      options = list(container = "body")
    ),
    htmlOutput(outputId = "message", inline = FALSE),
    sidebarMenu(
      menuItem(
        text = "Welcome",
        tabName = "tab_welcome"
      ),
      sidebarMenuOutput(outputId = "menu1"),
      sidebarMenuOutput(outputId = "menu2")
    ),
    htmlOutput(outputId = "containerid", inline = FALSE)
  ),
  dashboardBody(
    # fills background color to bottom of page when scrolling
    tags$head(tags$style(HTML('.content-wrapper { overflow: auto; }'))),
    tabItems(
    # Welcome tab
    tabItem(
      tabName = "tab_welcome",
      box(
        h3("Welcome to our app"),
        width = 12
      )
    ),
    # Preprocessing + QC Tab
    tabItem(
      tabName = "tab_preproc",
      fluidRow(
        box(
          title = p(
            "QC Filters",
            tags$style(type = "text/css", "#q2 {display: inline-block; vertical-align: middle;}"),
            bsButton(
              inputId = "q2",
              label = "",
              icon = icon("question"),
              style = "info",
              size = "extra-small"
            )
          ),
          column(
            disabled(numericInput(inputId = "num.ncountmin", label = NULL, value = 0)),
            disabled(numericInput(inputId = "num.nfeaturemin", label = NULL, value = 0)),
            disabled(numericInput(inputId = "num.mtmin", label = NULL, value = 0)),
            disabled(actionButton(inputId = "proc1", label = "Preprocess Input")),
            disabled(actionButton(inputId = "map", label = "Map cells to reference")),
            width = 6
          ),
          column(
            disabled(numericInput(inputId = "num.ncountmax", label = NULL, value = 0)),
            disabled(numericInput(inputId = "num.nfeaturemax", label = NULL, value = 0)),
            disabled(numericInput(inputId = "num.mtmax", label = NULL, value = 0)),
            width = 6
          ),
          width = 4
        ),
        bsPopover(
          id = "q2",
          title = "QC Filters",
          content = "Select a minimum and maximum value for nCount (number of molecules), nFeature (number of genes expressed), and mitochondrial percentage (if applicable)",
          placement = "right",
          trigger = "focus",
          options = list(container = "body")
        ),
        box(
          plotOutput(outputId = "plot.qc"),
          tableOutput(outputId = 'table.qc'),
          width = 8
        )
      ),
      fluidRow(
        valueBoxOutput(outputId = "valuebox.upload", width = 3),
        valueBoxOutput(outputId = "valuebox.preproc", width = 3),
        valueBoxOutput(outputId = "valuebox.mapped", width = 3)
      )
    ),
    # Cell tab
    tabItem(
      tabName = "tab_cell",
      box(
        title = "Reference",
        plotOutput(outputId = 'refdim'),
        width = 12
      ),
      box(
        title = "Query",
        disabled(selectInput(
          inputId = 'select.metadata',
          label = 'Metadata to color by',
          choices = '',
          selectize = FALSE,
          width = "25%"
        )),
        plotOutput(outputId = 'objdim'),
        width = 12
      ),
      box(
        title = "Metadata table",
        div(style="display: inline-block;vertical-align:top;width: 25%",
        disabled(selectInput(
          inputId = 'select.metadata1',
          label = 'Table rows',
          choices = '',
          selectize = FALSE
        ))),
        div(style="display: inline-block;vertical-align:top;width: 25%",
        disabled(selectInput(
          inputId = 'select.metadata2',
          label = 'Table columns',
          choices = '',
          selectize = FALSE
        ))),
        div(style="display: inline-block;vertical-align:top;width: 50%",
        disabled(radioButtons(
          inputId = 'radio.pct',
          label = NULL,
          choices = c('Percentage','Frequency'),
          inline = TRUE
        ))),
        tableOutput(outputId = 'table.metadata'),
        width = 12
      )
    ),
    # Feature tab
    tabItem(
      tabName = "tab_feature",
      box(
        title = "Feature Plots",
        div(style="display: inline-block;vertical-align:top;width: 33%",
        disabled(selectInput(
          inputId = "feature",
          label = "Feature",
          choices = '',
          selectize = FALSE
        ))),
        div(style="display: inline-block;vertical-align:top;width: 33%",
        disabled(selectInput(
          inputId = 'adtfeature',
          label = 'Imputed protein',
          choices = '',
          selectize = FALSE
        ))),
        div(style="display: inline-block;vertical-align:top;width: 33%",
        disabled(selectInput(
          inputId = 'scorefeature',
          label = "", #'Prediction Scores',
          choices = '',
          selectize = FALSE
        ))),
        plotOutput(outputId = 'edim'),
        plotOutput(outputId = 'evln'),
        width = 12
      ),
      box(
        title = p(
          "Predicted cell type cluster biomarkers",
          tags$style(type = "text/css", "#q3 {display: inline-block; vertical-align: middle;}"),
          bsButton(
            inputId = "q3",
            label = "",
            icon = icon("question"),
            style = "info",
            size = "extra-small"
          )
        ),
        bsPopover(
          id = "q3",
          title = "Biomarkers Table",
          content = "avgExpr: mean value of feature for cells in cluster; auc: area under ROC; padj: Benjamini-Hochberg adjusted p value; pct_in: percent of cells in the cluster with nonzero feature value; pct_out: percent of cells out of the cluster with nonzero feature value",
          placement = "right",
          trigger = "focus",
          options = list(container = "body")
        ),
        disabled(selectInput(
          inputId = 'select.biomarkers',
          label = 'Predicted cell type',
          choices = '',
          selectize = FALSE,
          width = "25%"
        )),
        column(
          h3("RNA biomarkers"),
          DTOutput(outputId = 'biomarkers'),
          width = 6
        ),
        column(
          h3("Imputed protein biomarkers"),
          DTOutput(outputId = 'adtbio'),
          width = 6
        ),
        width = 12
      )
    ),
    # Downloads tab
    tabItem(
      tabName = "tab_download",
      box(
        title = "Analysis script template ",
        verbatimTextOutput(outputId = "text.dlscript", placeholder = TRUE),
        disabled(downloadButton(
          outputId = 'dlscript',
          label = 'Download'
        )),
        width = 6
      ),
      box(
        title = "UMAP (Seurat Reduction RDS)",
        verbatimTextOutput(outputId = "text.dlumap"),
        disabled(downloadButton(
          outputId = 'dlumap',
          label = 'Download'
        )),
        width = 6
      ),
      box(
        title = "Imputed protein (Seurat Assay RDS)",
        verbatimTextOutput(outputId = "text.dladt"),
        disabled(downloadButton(
          outputId = 'dladt',
          label = 'Download'
        )),
        width = 6
      ),
      box(
       title = "Predicted cell types and scores (TSV)",
       verbatimTextOutput(outputId = "text.dlpred"),
       disabled(downloadButton(
         outputId = 'dlpred',
         label = 'Download'
       )),
       width = 6
      )
    )
  )
)))

#' Server function for the mapping app
#'
#' @param input,output,session Required Shiny app server parameters
#'
#' @return The shiny server logic
#'
#' @name AzimuthServer
#' @rdname AzimuthServer
#'
#' @importFrom methods slot<- slot
#' @importFrom ggplot2 ggtitle scale_colour_hue xlab geom_hline annotate
#' @importFrom presto wilcoxauc
#' @importFrom shinyjs show enable disable
#' @importFrom Seurat DefaultAssay PercentageFeatureSet SCTransform
#' VariableFeatures Idents GetAssayData RunUMAP CreateAssayObject
#' CreateDimReducObject Embeddings AddMetaData SetAssayData Key
#' VlnPlot DimPlot Reductions FeaturePlot Assays NoLegend Idents<- Cells
#' @importFrom shiny reactiveValues safeError appendTab observeEvent
#' withProgress setProgress updateSliderInput renderText updateSelectInput
#' updateTabsetPanel renderPlot renderTable downloadHandler renderUI
#' isolate
#' @importFrom shinydashboard renderValueBox valueBox renderMenu
#' @importFrom stats quantile
#' @importFrom utils write.table
#' @importFrom patchwork wrap_plots
#' @importFrom DT dataTableProxy selectRows renderDT
#'
#' @keywords internal
#'
server <- function(input, output, session) {
  mt.key <- 'percent.mt'
  mito.pattern <- getOption(x = 'Azimuth.app.mito', default = '^MT-')
  adt.key <- 'impADT'
  scores.key <- "scores"
  app.env <- reactiveValues(
    object = NULL,
    default.assay = NULL,
    default.feature = NULL,
    default.adt = NULL,
    feature = '',
    diff.exp = list(),
    messages = 'Upload a file'
  )
  rna.proxy <- dataTableProxy(outputId = "biomarkers")
  adt.proxy <- dataTableProxy(outputId = "adtbio")
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
  # React to events
  observeEvent( # Load the data
    eventExpr = input$file,
    handlerExpr = {
      # TODO disable file?
      withProgress(
        message = "Reading input",
        expr = {
          setProgress(value = 0)
          tryCatch(
            expr = {
              app.env$object <- LoadFileInput(path = input$file$datapath)
              app.env$object$query <- 'query'
              Idents(app.env$object) <- 'query'
            },
            error = function(e) { app.env$messages <- e$message }
          )
          setProgress(value = 1)
        }
      )
      if (!is.null(app.env$object)) {
        # Validate that there are genes in common with the reference
        genes.common <- intersect(
          y = rownames(x = refs$map),
          x = rownames(x = app.env$object)
        )
        if (length(x = genes.common) < getOption(x = "Azimuth.map.ngenes")) {
          app.env$object <- NULL
          gc(verbose = FALSE)
          app.env$messages <- "Not enough genes in common with reference. Try another dataset."
        }
        # Validate that there aren't too many cells
        else if (length(Cells(app.env$object)) > getOption(x = "Azimuth.app.max.cells")) {
          app.env$object <- NULL
          gc(verbose = FALSE)
          app.env$messages <- "Too many cells. Try another dataset."
        } else {
          app.env$default.assay <- DefaultAssay(object = app.env$object)
          ncount <- paste0('nCount_', app.env$default.assay)
          nfeature <- paste0('nFeature_', app.env$default.assay)
          if (!all(c(ncount, nfeature) %in% c(colnames(x = app.env$object[[]])))) {
            withProgress(
              message = "Calculating nCount and nFeature",
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
          updateNumericInput(
            session = session,
            inputId = 'num.ncountmin',
            label = paste("min", ncount),
            value = min(ncount.val),
            min = min(ncount.val),
            max = max(ncount.val),
            step = 1
          )
          enable(id = 'num.ncountmin')
          updateNumericInput(
            session = session,
            inputId = 'num.ncountmax',
            label = paste("max", ncount),
            value = max(ncount.val),
            min = min(ncount.val),
            max = max(ncount.val),
            step = 1
          )
          enable(id = 'num.ncountmax')
          nfeature.val <- range(app.env$object[[nfeature, drop = TRUE]])
          updateNumericInput(
            session = session,
            inputId = 'num.nfeaturemin',
            label = paste("min", nfeature),
            value = min(nfeature.val),
            min = min(nfeature.val),
            max = max(nfeature.val),
            step = 1
          )
          enable(id = 'num.nfeaturemin')
          updateNumericInput(
            session = session,
            inputId = 'num.nfeaturemax',
            label = paste("max", nfeature),
            value = max(nfeature.val),
            min = min(nfeature.val),
            max = max(nfeature.val),
            step = 1
          )
          enable(id = 'num.nfeaturemax')
          if (any(grepl(pattern = mito.pattern, x = rownames(x = app.env$object)))) {
            app.env$object <- PercentageFeatureSet(
              object = app.env$object,
              pattern = mito.pattern,
              col.name = mt.key,
              assay = app.env$default.assay
            )
            mito.val <- range(app.env$object[[mt.key, drop = TRUE]])
            updateNumericInput(
              session = session,
              inputId = 'num.mtmin',
              label = paste("min", mt.key),
              value = min(mito.val),
              min = min(mito.val),
              max = max(mito.val),
              step = 1
            )
            enable(id = 'num.mtmin')
            updateNumericInput(
              session = session,
              inputId = 'num.mtmax',
              label = paste("max", mt.key),
              value = max(mito.val),
              min = min(mito.val),
              max = max(mito.val),
              step = 1
            )
            enable(id = 'num.mtmax')
          }
          output$menu1 <- renderMenu(expr = {
            sidebarMenu(menuItem(
              text = "Preprocessing",
              tabName = "tab_preproc",
              icon = icon("filter"),
              selected = TRUE
            ))
          })
          ncellsupload <- length(x = colnames(x = app.env$object))
          app.env$messages <- paste(ncellsupload, "cells uploaded")
          if (ncellsupload < getOption(x = "Azimuth.map.ncells")) {
            output$valuebox.upload <- renderValueBox(expr = valueBox(
              value = ncellsupload,
              subtitle = "cells uploaded",
              icon = icon("times"),
              color = "red"
            ))
            # TODO what should happen when not enough cells are uploaded?
            # stop(safeError(error = "Not enough cells in uploaded dataset to proceed"))
          } else {
            output$valuebox.upload <- renderValueBox(expr = valueBox(
              value = ncellsupload,
              subtitle = "cells uploaded",
              icon = icon("check"),
              color = "green"
            ))
            enable(id = 'proc1')
          }
        }
      }
    }
  )
  observeEvent( # Process the user data
    eventExpr = input$proc1,
    handlerExpr = {
      disable(id = 'proc1')
      disable(id = 'num.ncountmin')
      disable(id = 'num.ncountmax')
      disable(id = 'num.nfeaturemin')
      disable(id = 'num.nfeaturemax')
      disable(id = 'num.mtmax')
      disable(id = 'num.mtmin')
      # Run SCTransform and enable mapping
      withProgress(
        message = "Normalizing with SCTransform",
        expr = {
          # output$sct <- renderText(expr = NULL)
          setProgress(
            value = 0,
            message = "Filtering based on nCount and nFeature"
          )
          ncount <- paste0('nCount_', DefaultAssay(object = app.env$object))
          nfeature <- paste0('nFeature_', DefaultAssay(object = app.env$object))
          cells.use <- app.env$object[[ncount, drop = TRUE]] >= input$num.ncountmin &
            app.env$object[[ncount, drop = TRUE]] <= input$num.ncountmax &
            app.env$object[[nfeature, drop = TRUE]] >= input$num.nfeaturemin &
            app.env$object[[nfeature, drop = TRUE]] <= input$num.nfeaturemax
          if (any(grepl(pattern = mito.pattern, x = rownames(x = app.env$object)))) {
            cells.use <- cells.use &
              app.env$object[[mt.key, drop = TRUE]] >= input$num.mtmin &
              app.env$object[[mt.key, drop = TRUE]] <= input$num.mtmax
          }
          ncellspreproc <- sum(cells.use)
          if (ncellspreproc < getOption(x = "Azimuth.map.ncells")) {
            output$valuebox.preproc <- renderValueBox(expr = valueBox(
              value = ncellspreproc,
              subtitle = "cells after filtering",
              icon = icon("times"),
              color = "red"
            ))
            # TODO what should happen when not enough cells are available after filtering?
            # stop(safeError(error = "Not enough cells after filtering to proceed"))
            # reset values and re-enable all filters; don't enable the mapping button; re-enable preproc button
            ncount.val <- range(app.env$object[[ncount, drop = TRUE]])
            updateNumericInput(
              session = session,
              inputId = 'num.ncountmin',
              label = paste("min", ncount),
              value = min(ncount.val),
              min = min(ncount.val),
              max = max(ncount.val),
              step = 1
            )
            enable(id = 'num.ncountmin')
            updateNumericInput(
              session = session,
              inputId = 'num.ncountmax',
              label = paste("max", ncount),
              value = max(ncount.val),
              min = min(ncount.val),
              max = max(ncount.val),
              step = 1
            )
            enable(id = 'num.ncountmax')
            nfeature.val <- range(app.env$object[[nfeature, drop = TRUE]])
            updateNumericInput(
              session = session,
              inputId = 'num.nfeaturemin',
              label = paste("min", nfeature),
              value = min(nfeature.val),
              min = min(nfeature.val),
              max = max(nfeature.val),
              step = 1
            )
            enable(id = 'num.nfeaturemin')
            updateNumericInput(
              session = session,
              inputId = 'num.nfeaturemax',
              label = paste("max", nfeature),
              value = max(nfeature.val),
              min = min(nfeature.val),
              max = max(nfeature.val),
              step = 1
            )
            enable(id = 'num.nfeaturemax')
            if (any(grepl(pattern = mito.pattern, x = rownames(x = app.env$object)))) {
              app.env$object <- PercentageFeatureSet(
                object = app.env$object,
                pattern = mito.pattern,
                col.name = mt.key,
                assay = app.env$default.assay
              )
              mito.val <- range(app.env$object[[mt.key, drop = TRUE]])
              updateNumericInput(
                session = session,
                inputId = 'num.mtmin',
                label = paste("min", mt.key),
                value = min(mito.val),
                min = min(mito.val),
                max = max(mito.val),
                step = 1
              )
              enable(id = 'num.mtmin')
              updateNumericInput(
                session = session,
                inputId = 'num.mtmax',
                label = paste("max", mt.key),
                value = max(mito.val),
                min = min(mito.val),
                max = max(mito.val),
                step = 1
              )
              enable(id = 'num.mtmax')
              enable(id = 'proc1')
            }
          } else {
            output$valuebox.preproc <- renderValueBox(expr = valueBox(
              value = ncellspreproc,
              subtitle = "cells after filtering",
              icon = icon("check"),
              color = "green"
            ))
            app.env$object <- app.env$object[, cells.use]
            setProgress(value = 0.2, message = "Normalizing with SCTransform")
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
            setProgress(value = 1)
            app.env$messages <- c(
              app.env$messages,
              paste(ncellspreproc, "cells preprocessed")
            )
            enable(id = "map")
          }
        }
      )
    }
  )
  observeEvent( # Map data
    eventExpr = input$map,
    handlerExpr = {
      disable(id = 'map')
      withProgress(
        message = 'Mapping data',
        expr = {
          setProgress(value = 0, message = "Finding anchors")
          cells <- colnames(x = app.env$object)
          # TODO: export FindTransferAnchors_Fast
          anchors <- Seurat:::FindTransferAnchors_Fast(
            reference = refs$map,
            query = app.env$object,
            reference.nn = refs$index,
            reference.nnidx = refs$index$annoy_index,
            reference.assay = DefaultAssay(object = refs$map),
            npcs = NULL,
            k.filter = NA,
            query.assay = 'SCT',
            reference.reduction = 'spca',
            normalization.method = 'SCT',
            features = rownames(x = refs$map),
            dims = 1:50
          )
          # TODO fail if not enough anchors (Azimuth.map.nanchors)
          setProgress(value = 0.6, message = 'Integrating data')
          # TODO: export IngestNewData_Fast
          ingested <- Seurat:::IngestNewData_Fast(
            reference = refs$map,
            query = app.env$object,
            dims = 1:50,
            transfer.anchors = anchors,
            reference.nnidx = refs$index$annoy_index,
            transfer.labels = Idents(object = refs$map),
            transfer.expression = GetAssayData(
              object = refs$map[['ADT']],
              slot = 'data'
            )
          )
          rm(anchors)
          gc(verbose = FALSE)
          setProgress(value = 0.8, message = "Running UMAP transform")
          slot(object = ingested, name = 'neighbors')[['query_ref.nn']] <- NNTransform(
            neighbors = ingested[['query_ref.nn']],
            meta.data = refs$map[[]]
          )
          app.env$object <- AddPredictions(
            object = app.env$object,
            preds = ingested$predicted.id,
            scores = ingested$predicted.id.score,
            preds.levels = levels(x = refs$map),
            preds.drop = TRUE
          )
          app.env$object[['umap.proj']] <- RunUMAP(
            object = ingested[['query_ref.nn']],
            reduction.model = refs$map[['jumap']],
            reduction.key = 'UMAP_'
          )
          suppressWarnings(expr = app.env$object[[adt.key]] <- CreateAssayObject(
            data = ingested[['transfer']][, cells]
          ))
          setProgress(value = 0.9, message = 'Calculating mapping metrics')
          app.env$object[['int']] <- CreateDimReducObject(
            embeddings = Embeddings(object = ingested[['int']])[cells, ],
            assay = app.env$default.assay
          )
          dsqr <- QueryReference(
            reference = refs$map,
            query = app.env$object,
            assay.query = app.env$default.assay,
            seed = 4
          )
          app.env$object <- AddMetaData(
            object = app.env$object,
            metadata = CalcMappingMetric(object = dsqr)
          )
          rm(dsqr, ingested)
          app.env$object <- SetAssayData(
            object = app.env$object,
            assay = 'SCT',
            slot = 'scale.data',
            new.data = new(Class = 'matrix')
          )
          gc(verbose = FALSE)
          setProgress(value = 1)
        }
      )
      mappingtext <- paste(sum(app.env$object$mapped)," cells mapped")
      mappingpct <- round(
        x = sum(app.env$object$mapped) / ncol(x = app.env$object) * 100,
        digits = 0
      )
      app.env$messages <- c(app.env$messages, mappingtext)
      if (mappingpct < getOption(x = "Azimuth.map.pcthresh")) {
        output$valuebox.mapped <- renderValueBox(expr = valueBox(
          value = paste0(mappingpct, "%"),
          subtitle = "cells mapped",
          icon = icon("times"),
          color = "red"
        ))
        # TODO what should happen when not enough cells map?
        # enable script download and download tab
        output$menu2 <- renderMenu(expr = {
          sidebarMenu(
            menuItem(
              text = "Download Results",
              tabName = "tab_download",
              icon = icon("file-download")
            ))}
        )
        # stop(safeError(error = "Query dataset could not be mapped to the reference"))
      } else {
        output$valuebox.mapped <- renderValueBox(expr = valueBox(
          value = paste0(mappingpct, "%"),
          subtitle = "cells mapped",
          icon = icon("check"),
          color = "green"
        ))
        # set unmapped cells predicted.id to NA
        app.env$object[["predicted.id"]][!app.env$object[["mapped"]]] <- NA
        # Enable the feature explorer
        enable(id = 'feature')
        app.env$default.feature <- ifelse(
          test = getOption(x = 'Azimuth.app.default.gene') %in% rownames(x = app.env$object),
          yes = getOption(x = 'Azimuth.app.default.gene'),
          no = VariableFeatures(object = app.env$object)[1]
        )
        updateSelectInput(
          session = session,
          inputId = 'feature',
          label = 'Feature',
          choices = FilterFeatures(features = rownames(x = app.env$object)),
          selected = app.env$default.feature
        )
        # Add the predicted ID and score to the plots
        enable(id = 'adtfeature')
        updateSelectInput(
          session = session,
          inputId = 'feature',
          label = 'Feature',
          choices = c(
            'predicted.id.score',
            FilterFeatures(features = rownames(x = app.env$object))
          ),
          selected = app.env$default.feature
        )
        adt.features <- sort(x = FilterFeatures(features = rownames(
          x = app.env$object[[adt.key]]
        )))
        app.env$default.adt <- ifelse(
          test = getOption(x = 'Azimuth.app.default.adt') %in% adt.features,
          yes = getOption(x = 'Azimuth.app.default.adt'),
          no = adt.features[1]
        )
        updateSelectInput(
          session = session,
          inputId = 'adtfeature',
          choices = adt.features,
          # selected = app.env$default.adt
          selected = ''
        )
        metadata.choices <- sort(x = c(
          "predicted.id",
          PlottableMetadataNames(
            object = app.env$object,
            min.levels = 1
          )
        ))
        updateSelectInput(
          session = session,
          inputId = 'select.metadata',
          choices = metadata.choices,
          selected = 'predicted.id'
        )
        updateSelectInput(
          session = session,
          inputId = 'select.metadata1',
          choices = metadata.choices,
          selected = 'predicted.id'
        )
        updateSelectInput(
          session = session,
          inputId = 'select.metadata2',
          choices = metadata.choices,
          selected = 'predicted.id'
        )
        enable(id = 'select.metadata')
        enable(id = 'select.metadata1')
        enable(id = 'select.metadata2')
        enable(id = 'radio.pct')
        # Compute biomarkers
        withProgress(
          message = "Running differential expression",
          expr = {
            setProgress(value = 0)
            app.env$diff.expr[[app.env$default.assay]] <- wilcoxauc(
              X = app.env$object,
              group_by = 'predicted.id',
              assay = 'data',
              seurat_assay = app.env$default.assay
            )
            setProgress(value = 0.6)
            app.env$diff.expr[[adt.key]] <- wilcoxauc(
              X = app.env$object,
              group_by = 'predicted.id',
              assay = 'data',
              seurat_assay = adt.key
            )
            setProgress(value = 1)
          }
        )
        allowed.clusters <- names(x = which(
          x = table(app.env$object$predicted.id) > getOption(x = 'Azimuth.de.mincells')
        ))
        allowed.clusters <- factor(
          x = allowed.clusters,
          levels = levels(x = app.env$object)
        )
        allowed.clusters <- levels(x = droplevels(x = allowed.clusters))
        enable(id = 'select.prediction')
        updateSelectInput(
          session = session,
          inputId = 'select.prediction',
          choices = allowed.clusters,
          selected = allowed.clusters[1]
        )
        enable(id = 'select.biomarkers')
        updateSelectInput(
          session = session,
          inputId = 'select.biomarkers',
          choices = allowed.clusters,
          selected = allowed.clusters[1]
        )
        # Enable downloads
        enable(id = 'dlumap')
        enable(id = 'dladt')
        enable(id = 'dlpred')
        output$text.dladt <- renderText(expr = {
          c("imputed.assay <- readRDS('azimuth_impADT.Rds')",
          "object <- object[, Cells(imputed.assay)]",
          "object[['impADT']] <- imputed.assay")
          },
          sep = "\n"
        )
        output$text.dlumap <- renderText(expr = {
          c("projected.umap <- readRDS('azimuth_umap.Rds')",
            "object <- object[, Cells(projected.umap)]",
            "object[['umap.proj']] <- projected.umap")
          },
          sep = "\n"
        )
        output$text.dlpred <- renderText(expr = {
          c("predictions <- read.delim('azimuth_pred.tsv', row.names = 1)",
            "object <- AddMetaData(",
            "\tobject = object,",
            "\tmetadata = predictions)")
          },
          sep = "\n"
        )
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
          ))}
        )
      }
    }
  )
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
          updateSelectInput(session = session, inputId = f, selected = '')
        }
        table.check <- input$feature %in% rownames(x = RenderDiffExp(
          diff.exp = app.env$diff.expr[[app.env$default.assay]],
          groups.use = input$select.biomarkers
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
          updateSelectInput(session = session, inputId = f, selected = '')
        }
        table.check <- input$adtfeature %in% rownames(x = RenderDiffExp(
          diff.exp = app.env$diff.expr[[adt.key]],
          groups.use = input$select.biomarkers
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
        app.env$feature <- paste0(
          Key(object = app.env$object[[scores.key]]),
          input$scorefeature
        )
        for (f in c('feature', 'adtfeature')) {
          updateSelectInput(session = session, inputId = f, selected = '')
        }
      }
    }
  )
  observeEvent( # Select from biomarkers table
    eventExpr = input$biomarkers_rows_selected,
    handlerExpr = {
      if (length(x = input$biomarkers_rows_selected)) {
        updateSelectInput(
          session = session,
          inputId = 'feature',
          selected = rownames(x = RenderDiffExp(
            diff.exp = app.env$diff.expr[[app.env$default.assay]],
            groups.use = input$select.biomarkers
          ))[input$biomarkers_rows_selected]
        )
      }
    }
  )
  observeEvent( # Select from adtbio table
    eventExpr = input$adtbio_rows_selected,
    handlerExpr = {
      if (length(x = input$adtbio_rows_selected)) {
        updateSelectInput(
          session = session,
          inputId = 'adtfeature',
          selected = rownames(x = RenderDiffExp(
            diff.exp = app.env$diff.expr[[adt.key]],
            groups.use = input$select.biomarkers
          ))[input$adtbio_rows_selected]
        )
      }
    }
  )
  # Plots
  output$plot.qc <- renderPlot(expr = {
  if (!is.null(x = isolate(app.env$object))) {
      qc <- paste0(c('nCount_', 'nFeature_'), app.env$default.assay)
      if (mt.key %in% colnames(x = isolate(app.env$object[[]]))) {
        qc <- c(qc, mt.key)
      }
      vlnlist <- VlnPlot(
        object = isolate(app.env$object),
        features = qc,
        group.by = 'query',
        combine = FALSE,
        pt.size = Seurat:::AutoPointSize(data = isolate(app.env$object))
      )
      # nCount
      vlnlist[[1]] <- vlnlist[[1]] +
        geom_hline(yintercept = input$num.ncountmin) +
        geom_hline(yintercept = input$num.ncountmax) +
        annotate(geom = "rect", alpha = 0.2, fill = "red",
                 ymin = input$num.ncountmax, ymax = Inf, xmin = 0.5, xmax = 1.5) +
        annotate(geom = "rect", alpha = 0.2, fill = "red",
                 ymin = -Inf, ymax = input$num.ncountmin, xmin = 0.5, xmax = 1.5) +
        NoLegend() + xlab("")
      # nFeature
      vlnlist[[2]] <- vlnlist[[2]] +
        geom_hline(yintercept = input$num.nfeaturemin) +
        geom_hline(yintercept = input$num.nfeaturemax) +
        annotate(geom = "rect", alpha = 0.2, fill = "red",
                 ymin = input$num.nfeaturemax, ymax = Inf, xmin = 0.5, xmax = 1.5) +
        annotate(geom = "rect", alpha = 0.2, fill = "red",
                 ymin = -Inf, ymax = input$num.nfeaturemin, xmin = 0.5, xmax = 1.5) +
        NoLegend() + xlab("")
      if (mt.key %in% colnames(x = isolate(app.env$object[[]]))) {
        vlnlist[[3]] <- vlnlist[[3]] +
          geom_hline(yintercept = input$num.mtmin) +
          geom_hline(yintercept = input$num.mtmax) +
          annotate(geom = "rect", alpha = 0.2, fill = "red",
                   ymin = input$num.mtmax, ymax = Inf, xmin = 0.5, xmax = 1.5) +
          annotate(geom = "rect", alpha = 0.2, fill = "red",
                   ymin = -Inf, ymax = input$num.mtmin, xmin = 0.5, xmax = 1.5) +
          NoLegend() + xlab("")
        wrap_plots(vlnlist, ncol = 3)
      } else {
        wrap_plots(vlnlist, ncol = 2)
      }
    }
  })
  output$refdim <- renderPlot(expr = {
    DimPlot(object = refs$plot)
  })
  output$objdim <- renderPlot(expr = {
  if (!is.null(x = app.env$object)) {
      if (length(x = Reductions(object = app.env$object))) {
        if (input$select.metadata == "predicted.id") {
          plotlevels <- c(levels(refs$plot$id)[levels(refs$plot$id) != "Doublet"], NA)
          DimPlot(
            object = app.env$object,
            group.by = "predicted.id") +
            scale_colour_hue(limits = plotlevels, drop = FALSE)
        } else {
          DimPlot(
            object = app.env$object,
            group.by = input$select.metadata
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
        colnames(x = app.env$object[[]])
      )
      if (app.env$feature %in% avail) {
        VlnPlot(object = app.env$object, features = app.env$feature) +
          NoLegend()
      }
    }
  })
  output$edim <- renderPlot({
    if (!is.null(x = app.env$object)) {
      palettes <- list(
        c("lightgrey", "blue"),
        c('lightgrey', 'darkgreen'),
        c('lightgrey', 'blue')
      )
      names(x = palettes) <- c(
        Key(object = app.env$object[["SCT"]]),
        Key(object = app.env$object[[adt.key]]),
        'md_'
      )
      feature.key <- paste0(
        unlist(x = strsplit(x = app.env$feature, split = '_'))[1],
        '_'
      )
      if (feature.key == paste0(app.env$feature, '_')) {
        feature.key <- 'md_'
      }
      pal.use <- palettes[[feature.key]]
      if (!is.null(x = pal.use)) {
        FeaturePlot(
          object = app.env$object,
          features = app.env$feature,
          cols = pal.use
        )
      }
    }
  })
  # Messages
  output$message <- renderUI(expr = {
    p(HTML(text = paste(app.env$messages, collapse = "<br />")))
  })
  output$containerid <- renderText(c("debug ID: ", Sys.info()[["nodename"]]))
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
          groups.use = input$select.biomarkers
        )
      }
    },
    selection = 'single',
    options = list(dom = 't', ordering = FALSE)
  )
  output$adtbio <- renderDT(
    expr = {
      if (!is.null(x = app.env$diff.expr[[adt.key]])) {
        RenderDiffExp(
          diff.exp = app.env$diff.expr[[adt.key]],
          groups.use = input$select.biomarkers
        )
      }
    },
    selection = 'single',
    options = list(dom = 't', ordering = FALSE)
  )
  output$table.metadata <- renderTable(
    expr = {
      if (!is.null(x = app.env$object)) {
        # TODO allow user to toggle between frequency and percentage
        CategoryTable(
          object = app.env$object,
          category.1 = input$select.metadata1,
          category.2 = input$select.metadata2,
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
  # TODO: Add data from CalcMappingMetric to _pred.tsv?
  output$dlpred <- downloadHandler(
    filename = paste0(tolower(x = app.title), '_pred.tsv'),
    content = function(file) {
      req <- c('predicted.id', 'predicted.id.score')
      if (all(req %in% colnames(x = app.env$object[[]]))) {
        write.table(
          x = data.frame(
            cell = colnames(x = app.env$object),
            predicted.id = app.env$object$predicted.id,
            predicted.score = app.env$object$predicted.id.score,
            stringsAsFactors = FALSE
          ),
          file = file,
          quote = FALSE,
          row.names = FALSE,
          col.names = TRUE,
          sep = '\t'
        )
      }
    }
  )
}

#' Launch the mapping app
#'
#' @param reference,mito,max.upload.mb,max.cells,default.gene,default.adt See \strong{App options} for more details
#'
#' @section App options:
#'
#' The following options are used for passing information into the app; users
#' can configure these either \link[base:options]{globally} or via arguments to
#' \code{\link{AzimuthApp}} (omitting the \dQuote{Azimuth.app} prefix):
#'
#' \describe{
#'  \item{\code{Azimuth.app.mito}}{
#'   Regular expression pattern indicating mitochondrial features in query object
#'  }
#'  \item{\code{Azimuth.app.reference}}{
#'   URL or directory path to reference dataset; see \code{\link{LoadReference}}
#'   for more details
#'  }
#'  \item{\code{Azimuth.app.max.upload.mb}}{
#'   Maximum file size (in MB) allowed to upload
#'  }
#'  \item{\code{Azimuth.app.max.cells}}{
#'   Maximum number of cells allowed to upload
#'  }
#'  \item{\code{Azimuth.app.default.gene}}{
#'   Gene to select by default in feature/violin plot
#'  }
#'  \item{\code{Azimuth.app.default.adt}}{
#'   ADT to select by default in feature/violin plot
#'  }
#' }
#'
#' @return None, launches the mapping Shiny app
#'
#' @importFrom shiny runApp shinyApp
#' @importFrom withr with_options
#'
#' @export
#'
#' @seealso \code{\link{SeuratMapper-package}}
#'
AzimuthApp <- function(
  mito = getOption(x = 'Azimuth.app.mito', default = '^MT-'),
  reference = getOption(
    x = 'Azimuth.app.reference',
    default = 'http://satijalab04.nygenome.org/pbmc'
  ),
  max.upload.mb = getOption(
    x = 'Azimuth.app.max.upload.mb',
    default = 500
  ),
  max.cells = getOption(
    x = 'Azimuth.app.max.cells',
    default = 50000
  ),
  default.gene = getOption(
    x = 'Azimuth.app.default.gene',
    default = "GNLY"
  ),
  default.adt = getOption(
    x = 'Azimuth.app.default.adt',
    default = "CD3-1"
  )
) {
  useShinyjs()
  opts <- list(
    shiny.maxRequestSize = max.upload.mb * (1024 ^ 2),
    Azimuth.app.mito = mito,
    Azimuth.app.reference = reference,
    Azimuth.app.max.cells = max.cells,
    Azimuth.app.default.gene = default.gene,
    Azimuth.app.default.adt = default.adt
  )
  with_options(
    new = opts,
    code = runApp(appDir = shinyApp(ui = ui, server = server))
  )
  return(invisible(x = NULL))
}

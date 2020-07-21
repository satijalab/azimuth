#' @include zzz.R
#' @importFrom htmltools tagList
#' @importFrom shinyjs useShinyjs disabled
#' @importFrom shiny fluidPage sidebarLayout sidebarPanel fileInput sliderInput
#' actionButton selectInput downloadButton mainPanel tabsetPanel tabPanel
#' plotOutput verbatimTextOutput shinyApp
#'
NULL

# app.title <- 'Azimuth'
app.title <- 'SuperCoolMappingNameTBD'

ui <- tagList(
  useShinyjs(),
  fluidPage(
    title = app.title,
    sidebarLayout(
      sidebarPanel = sidebarPanel(
        fileInput(inputId = "file", label = app.title),
        disabled(sliderInput(
          inputId = 'ncount',
          label = 'nCount',
          min = 0L,
          max = 1L,
          value = c(0L, 1L),
        )),
        disabled(sliderInput(
          inputId = 'nfeature',
          label = 'nFeature',
          min = 0L,
          max = 1L,
          value = c(0L, 1L),
        )),
        disabled(actionButton(inputId = "proc1", label = "Preprocess Input")),
        disabled(actionButton(inputId = "map", label = "Map cells to reference")),
        disabled(selectInput(
          inputId = 'feature',
          label = 'Feature',
          choices = '',
          selectize = FALSE,
          width = '100%'
        )),
        disabled(downloadButton(
          outputId = 'dlumap',
          label = 'Download UMAP RDS'
        )),
        disabled(downloadButton(
          outputId = 'dlpred',
          label = 'Download the predicted IDs and scores'
        ))
      ),
      mainPanel = mainPanel(tabsetPanel(
        id = "tabs",
        tabPanel(
          title = "Preprocessing and Quality Control",
          plotOutput(outputId = "qcvln"),
          verbatimTextOutput(outputId = "sct"),
          verbatimTextOutput(outputId = "mapping")
        )
      ))
    )
  )
)

#' Server function for the mapping app
#'
#' @param input,output,session Required Shiny app server parameters
#'
#' @return The shiny server logic
#'
#' @name AzimuthServer
#' @rdname AzimuthServer
#'
#' @importFrom shinyjs enable disable
#' @importFrom Seurat Idents<- DefaultAssay SCTransform VariableFeatures Idents
#' RunUMAP VlnPlot DimPlot Reductions FeaturePlot
#' @importFrom shiny reactiveValues appendTab tabPanel plotOutput observeEvent
#' withProgress setProgress updateSliderInput renderText updateSelectInput
#' updateTabsetPanel renderPlot downloadHandler
#'
#' @keywords internal
#'
server <- function(input, output, session) {
  app.env <- reactiveValues(object = NULL)
  refs <- LoadReference(path = "http://saucyx220.nygenome.org")
  Idents(object = refs$reference) <- 'id'
  appendTab(
    inputId = 'tabs',
    tab = tabPanel(
      title = "Mapped Data",
      plotOutput(outputId = 'refdim'),
      plotOutput(outputId = 'objdim')
    )
  )
  # React to events
  observeEvent( # Load the data
    eventExpr = input$file,
    handlerExpr = {
      withProgress(
        message = "Reading input",
        expr = {
          setProgress(value = 0)
          app.env$object <- LoadFileInput(path = input$file$datapath)
          setProgress(value = 1)
        }
      )
      enable(id = 'ncount')
      enable(id = 'nfeature')
      ncount <- paste0('nCount_', DefaultAssay(object = app.env$object))
      nfeature <- paste0('nFeature_', DefaultAssay(object = app.env$object))
      ncount.val <- range(app.env$object[[ncount, drop = TRUE]])
      nfeature.val <- range(app.env$object[[nfeature, drop = TRUE]])
      updateSliderInput(
        session = session,
        inputId = 'ncount',
        label = ncount,
        value = ncount.val,
        min = min(ncount.val),
        max = max(ncount.val)
      )
      updateSliderInput(
        session = session,
        inputId = 'nfeature',
        label = nfeature,
        value = nfeature.val,
        min = min(nfeature.val),
        max = max(nfeature.val)
      )
    }
  )
  observeEvent(eventExpr = input$file, handlerExpr = enable(id = 'proc1'))
  observeEvent( # Process the user data
    eventExpr = input$proc1,
    handlerExpr = {
      # Run SCTransform and enable mapping
      withProgress(
        message = "Normalizing with SCTransform",
        expr = {
          output$sct <- renderText(expr = NULL)
          setProgress(
            value = 0,
            message = "Filtering based on nCount and nFeature"
          )
          ncount <- paste0('nCount_', DefaultAssay(object = app.env$object))
          nfeature <- paste0('nFeature_', DefaultAssay(object = app.env$object))
          cells.use <- app.env$object[[ncount, drop = TRUE]] >= min(input$ncount) &
            app.env$object[[ncount, drop = TRUE]] <= max(input$ncount) &
            app.env$object[[nfeature, drop = TRUE]] >= min(input$nfeature) &
            app.env$object[[nfeature, drop = TRUE]] <= max(input$nfeature)
          app.env$object <- app.env$object[, cells.use]
          disable(id = 'ncount')
          disable(id = 'nfeature')
          setProgress(value = 0.2, message = "Normalizing with SCTransform")
          app.env$object <- suppressWarnings(expr = SCTransform(
            object = app.env$object,
            residual.features = rownames(x = refs$reference),
            ncells = min(3000, ncol(x = app.env$object))
          ))
          setProgress(value = 1)
          output$sct <- renderText(expr = "SCTransform complete")
        }
      )
      enable(id = "map")
      # Enable the feature explorer
      enable(id = 'feature')
      updateSelectInput(
        session = session,
        inputId = 'feature',
        label = 'Feature',
        choices = FilterFeatures(features = rownames(x = app.env$object)),
        selected = VariableFeatures(object = app.env$object)[1]
      )
      appendTab(
        inputId = 'tabs',
        tabPanel(
          title = "Feature Explorer",
          plotOutput(outputId = 'fvln'),
          plotOutput(outputId = 'fdim')
        )
      )
    }
  )
  observeEvent( # Map data
    eventExpr = input$map,
    handlerExpr = {
      withProgress(
        message = 'Mapping data',
        expr = {
          setProgress(value = 0, message = "Finding anchors")
          # TODO: export FindTransferAnchors_Fast
          anchors <- Seurat:::FindTransferAnchors_Fast(
            reference = refs$reference,
            query = app.env$object,
            reference.nn = refs$index,
            reference.nnidx = refs$index$annoy_index,
            reference.assay = DefaultAssay(object = refs$reference),
            npcs = NULL,
            k.filter = NA,
            query.assay = 'SCT',
            reference.reduction = 'spca',
            normalization.method = 'SCT',
            features = rownames(x = refs$reference),
            dims = 1:50
          )
          setProgress(value = 0.6, message = 'Integrating data')
          # TODO: export IngestNewData_Fase
          ingested <- Seurat:::IngestNewData_Fast(
            reference = refs$reference,
            query = app.env$object,
            dims = 1:50,
            transfer.anchors = anchors,
            reference.nnidx = refs$index$annoy_index,
            transfer.labels = Idents(object = refs$reference)
          )
          setProgress(value = 0.8, message = "Running UMAP transform")
          cells <- colnames(x = app.env$object)
          app.env$object$predicted.id <- ingested$predicted.id[cells]
          app.env$object$predicted.id.score <- ingested$predicted.id.score[cells]
          app.env$object[['umap.proj']] <- RunUMAP(
            object = ingested[['query_ref.nn']],
            reduction.model = refs$reference[['jumap']],
            reduction.key = 'ProjU_'
          )
          setProgress(value = 1)
        }
      )
      # Add the predicted ID and score to the plots
      Idents(object = app.env$object) <- 'predicted.id'
      updateSelectInput(
        session = session,
        inputId = 'feature',
        label = 'Feature',
        choices = c(
          'predicted.id.score',
          FilterFeatures(features = rownames(x = app.env$object))
        ),
        selected = 'predicted.id.score'
      )
      updateTabsetPanel(
        session = session,
        inputId = 'tabs',
        selected = 'Mapped Data'
      )
      # Enable downloads
      enable(id = 'dlumap')
      enable(id = 'dlpred')
    }
  )
  # Plots
  output$qcvln <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      qc <- paste0(
        c('nCount_', 'nFeature_'),
        DefaultAssay(object = app.env$object)
      )
      VlnPlot(object = app.env$object, features = qc)
    }
  })
  output$refdim <- renderPlot(expr = {
    DimPlot(object = refs$reference)
  })
  output$objdim <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      if (length(x = Reductions(object = app.env$object))) {
        DimPlot(object = app.env$object)
      }
    }
  })
  output$fvln <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      avail <- c(rownames(x = app.env$object), colnames(x = app.env$object[[]]))
      if (input$feature %in% avail) {
        VlnPlot(object = app.env$object, features = input$feature)
      }
    }
  })
  output$fdim <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      if (length(x = Reductions(object = app.env$object))) {
        FeaturePlot(object = app.env$object, features = input$feature)
      }
    }
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
#' @return None, launches the mapping Shiny app
#'
#' @importFrom shiny runApp
#'
#' @export
#'
AzimuthApp <- function() {
  runApp(appDir = shinyApp(ui = ui, server = server))
  return(invisible(x = NULL))
}

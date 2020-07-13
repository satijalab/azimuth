#' @importFrom shiny NS moduleServer reactive validate need
#'
NULL

#' Create a \code{\link[Seurat]{DimPlot}} Shiny Module
#'
#' @param id ...
#' @param object A \link[shiny]{reactive} \code{\link[Seurat]{Seurat}} object
#'
#' @return ...
#'
#' @importFrom Seurat DimPlot
#' @importFrom shiny renderPlot
#'
#' @export
#'
dimPlotServer <- function(id, object, ...) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      # object <- reactive(x = {
      #   validate(need(expr = object, message = FALSE))
      #   object
      # })
      output$dimplot <- renderPlot(
        expr = DimPlot(object = object())
      )
    }
  )
}

#' Load a 10X Dataset
#'
#' @param id ...
#'
#' @return A \link[shiny]{reactive} \code{\link[Seurat]{Seurat}} object
#'
#' @name tenxFileServer
#' @rdname tenxFileServer
#'
#' @importFrom Seurat Read10X_h5 CreateSeuratObject
#'
#' @export
#'
tenxH5Server <- function(id, ...) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      # userFile <- reactive(x = {
      #   validate(need(expr = input$file, message = FALSE))
      #   input$file
      # })
      counts <- shiny::eventReactive(
        eventExpr = input$file,
        valueExpr = {
          # message("Reading in data from ", userFile()$name)
          message("Reading in data from ", input$file()$name)
          mat <- Read10X_h5(filename = input$file()$datapath)
          if (is.list(x = mat)) {
            idx <- grep(pattern = 'Gene Expression', x = names(x = mat))
            if (!length(x = idx)) {
              idx <- 1
            }
            mat <- mat[[idx]]
          }
          mat
        }
      )
      return(counts)
      # counts <- reactive(
      #   x = {
      #
      #   },
      #   label = "Load Data"
      # )
      # shiny::observe(x = message("Uploaded ", userFile()$name))
      # return(counts)
    }
  )
}

#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA RunUMAP
stdProcessingServer <- function(id, object, ...) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      object <- reactive(x = {
        validate(need(expr = object, message = FALSE))
        object
      })
      observe("Processing data")
      object <- reactive(x = NormalizeData(object = object()))
      object <- reactive(x = FindVariableFeatures(object = object()))
      object <- reactive(x = ScaleData(object = object()))
      object <- reactive(x = RunPCA(object = object()))
      object <- reactive(x = RunUMAP(object = object()))
      return(object)
    }
  )
}

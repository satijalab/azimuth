#' @importFrom shiny NS plotOutput
#'
NULL

#' View a \code{\link[Seurat]{DimPlot}} in a Shiny App
#'
#' @param id ...
#'
#' @return ...
#'
#' @export
#'
dimPlotUI <- function(id) {
  ns <- NS(namespace = id)
  plotOutput(outputId = ns("dimplot"))
}

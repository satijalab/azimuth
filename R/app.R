#' @include ui.R
#' @include input.R
#' @include servers.R
#'
NULL

#' @importFrom shiny fluidPage sidebarLayout sidebarPanel mainPanel
#'
ui <- fluidPage(
  title = "SeuratMapper",
  sidebarLayout(
    sidebarPanel =  sidebarPanel(
      # tenxH5Input(id = "tenxh5")
      fileInput(inputId = "tenxh5", label = "10X HDF5 Input")
    ),
    mainPanel = mainPanel(
      # dimPlotUI(id = "reference"),
      # dimPlotUI(id = "user")
      shiny::plotOutput(outputId = "reference"),
      shiny::plotOutput(outputId = "user")
    )
  )
)

#' @importFrom SeuratDisk LoadH5Seurat
#'
server <- function(input, output, session) {
  # Load the reference object
  reference <- LoadH5Seurat(
    file = system.file(
      file.path('references', 'pbmc3k_final.h5Seurat'),
      package = 'SeuratMapper',
      mustWork = TRUE
    ),
    assays = "data",
    reductions = "umap",
    graphs = FALSE,
    images = FALSE,
    meta.data = FALSE,
    commands = FALSE
  )
  plot.env <- shiny::reactiveValues()
  # Load the user-supplied data
  observeEvent(
    eventExpr = input$tenxh5,
    handlerExpr = {
      # print(input$tenxh5)
      # print(class(input$tenxh5))
      mat <- Seurat::Read10X_h5(input$tenxh5$datapath)
      if (is.list(x = mat)) {
        mat <- mat[[1]]
      }
      object <- Seurat::CreateSeuratObject(counts = mat)
      object <- Seurat::NormalizeData(object)
      object <- Seurat::FindVariableFeatures(object)
      object <- Seurat::ScaleData(object)
      object <- Seurat::RunPCA(object)
      object <- Seurat::RunUMAP(object, dims = 1:30)
      plot.env$user <- object
    }
  )
  # counts <- tenxH5Server(id = "tenxh5")
  # object <- reactive(x = {
  #   obj <- CreateSeuratObject(counts = counts())
  #   message("Created Seurat object")
  #   obj
  # })
  # Create plots
  # dimPlotServer(id = "reference", object = reactive(x = reference))
  # dimPlotServer(id = "user", object = reactive(reference))
  output$reference <- shiny::renderPlot(
    Seurat::DimPlot(object = reference)
  )
  output$user <- shiny::renderPlot(expr = Seurat::DimPlot(plot.env$user))
}

app <- shinyApp(ui = ui, server = server)

# TODO:
# Download:
# umap coordinates
# number of cells per identity
# prediction score
# feature plot on query (separate tab?)

#' Launch the Mapper App
#'
#' @return None, launches the Shiny Mapper app
#'
#' @importFrom withr with_options
#'
#' @export
#'
Mapper <- function() {
  runApp(appDir = app)
  return(invisible(x = NULL))
}

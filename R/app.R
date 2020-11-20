#' @include zzz.R
#' @include seurat.R
#' @include helpers.R
#' @importFrom DT DTOutput
#' @importFrom htmltools tagList h4 hr h3 tags HTML p div a
#' @importFrom shinyjs useShinyjs disabled
#' @importFrom shiny fluidPage sidebarLayout sidebarPanel fileInput sliderInput
#' actionButton selectizeInput downloadButton mainPanel tabsetPanel tabPanel
#' plotOutput tableOutput verbatimTextOutput numericInput icon fluidRow
#' radioButtons textOutput htmlOutput column checkboxInput
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#' dashboardBody menuItem tabItems tabItem sidebarMenu box valueBoxOutput
#' sidebarMenuOutput
#' @importFrom shinyBS bsPopover bsButton
#'
NULL

app.title <- 'Azimuth'
selectize.opts <- list(
  maxOptions = 10000L,
  maxItems = 1L
)

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
    actionButton(inputId = "triggerdemo", label = "Load demo dataset"),
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
      htmlOutput(outputId = "welcomebox")
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
            textOutput(outputId = "text.cellsremain"),
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
          checkboxInput(inputId = 'check.qcscale', label = 'Log-scale Y-axis'),
          checkboxInput(inputId = 'check.qcpoints', label = 'Hide points'),
          plotOutput(outputId = "plot.qc"),
          tableOutput(outputId = 'table.qc'),
          width = 8
        )
      ),
      fluidRow(
        valueBoxOutput(outputId = "valuebox.upload", width = 3),
        valueBoxOutput(outputId = "valuebox.preproc", width = 3),
        valueBoxOutput(outputId = "valuebox.mapped", width = 3)
      ),
      fluidRow(
        class = "rowhide",
        box(
          plotOutput(outputId = "plot.pbcor"),
          width = 8
        )
      )
    ),
    # Cell tab
    tabItem(
      tabName = "tab_cell",
      box(
        title = "Reference",
        checkboxInput(inputId = 'labels', label = 'Show labels'),
        plotOutput(outputId = 'refdim'),
        width = 12
      ),
      box(
        title = "Query",
        disabled(selectizeInput(
          inputId = 'select.metadata',
          label = 'Metadata to color by',
          choices = '',
          width = "25%"
        )),
        plotOutput(outputId = 'objdim'),
        width = 12
      ),
      box(
        title = "Metadata table",
        div(
          style="display: inline-block;vertical-align:top;width: 25%",
          disabled(selectizeInput(
            inputId = 'select.metadata1',
            label = 'Table rows',
            choices = ''
          ))
        ),
        div(
          style="display: inline-block;vertical-align:top;width: 25%",
          disabled(selectizeInput(
            inputId = 'select.metadata2',
            label = 'Table columns',
            choices = ''
          ))
        ),
        div(
          style="display: inline-block;vertical-align:top;width: 50%",
          disabled(radioButtons(
            inputId = 'radio.pct',
            label = NULL,
            choices = c('Percentage','Frequency'),
            inline = TRUE
          ))
        ),
        tableOutput(outputId = 'table.metadata'),
        width = 12
      )
    ),
    # Feature tab
    tabItem(
      tabName = "tab_feature",
      box(
        title = "Feature Plots",
        div(
          style="display: inline-block;vertical-align:top;width: 33%",
          disabled(selectizeInput(
            inputId = 'feature',
            label = 'Feature',
            choices = ''
          ))
        ),
        div(
          style="display: inline-block;vertical-align:top;width: 33%",
          disabled(selectizeInput(
            inputId = 'adtfeature',
            label = 'Imputed protein',
            choices = ''
          ))
        ),
        div(
          style="display: inline-block;vertical-align:top;width: 33%",
          disabled(selectizeInput(
            inputId = 'scorefeature',
            label = "Prediction Scores and Continuous Metadata",
            choices = ''
          ))
        ),
        plotOutput(outputId = 'edim'),
        disabled(selectizeInput(
          inputId = 'groupfeature',
          label = 'Metadata to group by',
          choices = '',
          width = "25%"
        )),
        checkboxInput(inputId = 'check.featpoints', label = 'Hide points'),
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
          # content = "Only available for clusters with at least 15 cells. logFC: log fold-change between cells in the cluster specified and other cells; auc: area under ROC; padj: Benjamini-Hochberg adjusted p value; pct_in: percent of cells in the cluster with nonzero feature value; pct_out: percent of cells out of the cluster with nonzero feature value",
          content = "Only available for clusters with at least 15 cells. auc: area under ROC; padj: Benjamini-Hochberg adjusted p value; pct_in: percent of cells in the cluster with nonzero feature value; pct_out: percent of cells out of the cluster with nonzero feature value",
          placement = "right",
          trigger = "focus",
          options = list(container = "body")
        ),
        disabled(selectizeInput(
          inputId = 'select.biomarkers',
          label = 'Predicted cell type',
          choices = '',
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


#' Set up uploaded or demo file for QC
#'
#' @param app.env App environment
#' @param refs Reference object as from \code{\link{LoadReference}}
#' @param session Shiny app session
#' @param output Shiny app output
#' @param googlesheet Google sheet object as from \code{googlesheets4::gs4_get()}
#' @param mt.key String identifier for mitochondrial content
#' @param mito.pattern String pattern for identifying mitochondrial genes
#'
#' @importFrom Seurat DefaultAssay PercentageFeatureSet AddMetaData Cells
#' @importFrom shiny withProgress setProgress updateNumericInput
#' @importFrom shinydashboard renderValueBox valueBox renderMenu
#' @importFrom googlesheets4 sheet_append
#'
initQC <- function(app.env,
                   refs,
                   session,
                   output,
                   googlesheet,
                   mt.key,
                   mito.pattern){
  n.trees <- getOption(x = "Azimuth.map.ntrees")
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
    else if (length(Cells(app.env$object)) > getOption(x = "Azimuth.app.max_cells")) {
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
        value = floor(min(ncount.val)),
        min = floor(min(ncount.val)),
        max = ceiling(max(ncount.val))
      )
      enable(id = 'num.ncountmin')
      updateNumericInput(
        session = session,
        inputId = 'num.ncountmax',
        label = paste("max", ncount),
        value = ceiling(max(ncount.val)),
        min = floor(min(ncount.val)),
        max = ceiling(max(ncount.val))
      )
      enable(id = 'num.ncountmax')
      nfeature.val <- range(app.env$object[[nfeature, drop = TRUE]])
      updateNumericInput(
        session = session,
        inputId = 'num.nfeaturemin',
        label = paste("min", nfeature),
        value = floor(min(nfeature.val)),
        min = floor(min(nfeature.val)),
        max = ceiling(max(nfeature.val))
      )
      enable(id = 'num.nfeaturemin')
      updateNumericInput(
        session = session,
        inputId = 'num.nfeaturemax',
        label = paste("max", nfeature),
        value = ceiling(max(nfeature.val)),
        min = floor(min(nfeature.val)),
        max = ceiling(max(nfeature.val))
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
          value = floor(min(mito.val)),
          min = floor(min(mito.val)),
          max = ceiling(max(mito.val))
        )
        enable(id = 'num.mtmin')
        updateNumericInput(
          session = session,
          inputId = 'num.mtmax',
          label = paste("max", mt.key),
          value = ceiling(max(mito.val)),
          min = floor(min(mito.val)),
          max = ceiling(max(mito.val))
        )
        enable(id = 'num.mtmax')
        enable(id = 'check.qcscale')
        enable(id = 'check.qcpoints')
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
      # Not enough cells are uploaded
      if (ncellsupload < getOption(x = "Azimuth.map.ncells")) {
        output$valuebox.upload <- renderValueBox(expr = valueBox(
          value = ncellsupload,
          subtitle = "cells uploaded",
          icon = icon("times"),
          color = "red"
        ))
      } else {
        output$valuebox.upload <- renderValueBox(expr = valueBox(
          value = ncellsupload,
          subtitle = "cells uploaded",
          icon = icon("check"),
          color = "green"
        ))
        enable(id = 'map')
        if (!is.null(googlesheet)) {
          try(sheet_append(ss = googlesheet,
                           data = data.frame("CELLSUPLOAD", Sys.info()[["nodename"]], ncellsupload)))
        }
      }
    }
  }
}

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
#' theme_void
#' @importFrom presto wilcoxauc
#' @importFrom stringr str_interp
#' @importFrom shinyjs show hide enable disable
#' @importFrom Seurat DefaultAssay PercentageFeatureSet SCTransform
#' VariableFeatures Idents GetAssayData RunUMAP CreateAssayObject
#' CreateDimReducObject Embeddings AddMetaData SetAssayData Key
#' VlnPlot DimPlot Reductions FeaturePlot Assays NoLegend Idents<- Cells
#' FindTransferAnchors Misc Key<- RenameCells MappingScore
#' GetIntegrationData TransferData IntegrateEmbeddings FindNeighbors Tool
#' @importFrom shiny reactiveValues safeError appendTab observeEvent
#' withProgress setProgress updateNumericInput updateSliderInput renderText
#' updateSelectizeInput updateTabsetPanel renderPlot renderTable downloadHandler
#' renderUI isolate onStop
#' @importFrom shinydashboard renderValueBox valueBox renderMenu
#' @importFrom patchwork wrap_plots
#' @importFrom utils packageVersion write.table
#' @importFrom stats quantile na.omit
#' @importFrom DT dataTableProxy selectRows renderDT
#' @importFrom future future plan value resolved
#' @importFrom googlesheets4 gs4_auth gs4_get sheet_append
#'
#' @keywords internal
#'
server <- function(input, output, session) {
  # hide demo dataset button if required
  if (is.null(getOption(x = 'Azimuth.app.demodataset'))){
    hide(id="triggerdemo")
  }
  mt.key <- 'percent.mt'
  mito.pattern <- getOption(x = 'Azimuth.app.mito', default = '^MT-')
  adt.key <- 'impADT'
  n.trees <- getOption(x = "Azimuth.map.ntrees")
  app.env <- reactiveValues(
    object = NULL,
    default.assay = NULL,
    default.feature = NULL,
    default.adt = NULL,
    mapping.score = NULL,
    feature = '',
    diff.exp = list(),
    messages = 'Upload a file',
    features = character(length = 0L),
    adt.features = character(length = 0L),
    scorefeatures = character(length = 0L)
  )
  rna.proxy <- dataTableProxy(outputId = "biomarkers")
  adt.proxy <- dataTableProxy(outputId = "adtbio")
  if (!is.null(getOption(x = "Azimuth.app.googlesheet")) &&
      !is.null(getOption(x = "Azimuth.app.googletoken")) &&
      !is.null(getOption(x = "Azimuth.app.googletokenemail"))) {
    tryCatch(
      expr = {
        gs4_auth(email = getOption(x = "Azimuth.app.googletokenemail"),
                 cache = getOption(x = "Azimuth.app.googletoken"))
        googlesheet <- gs4_get(ss = getOption(x = "Azimuth.app.googlesheet"))
        app_start_time <- Sys.time()
        onStop(fun = function() try(sheet_append(ss = googlesheet, data = data.frame("SESSIONLENGTH", Sys.info()[["nodename"]], as.numeric(Sys.time() - app_start_time, units = "mins")))))
      },
      error = function(e) {
        googlesheet <- NULL
      }
    )
  } else {
    googlesheet <- NULL
  }
  if (!is.null(googlesheet)) {
    try(sheet_append(ss = googlesheet, data = data.frame("STARTUPTIME", Sys.info()[["nodename"]], Sys.time())))
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
  observeEvent( # Load the data
    eventExpr = input$file,
    handlerExpr = {
      # In case this is a re-upload, reset some elements
      app.env$messages <- ""
      output$valuebox.upload <- NULL
      output$valuebox.preproc <- NULL
      output$valuebox.mapped <- NULL
      output$menu2 <- NULL
      disable(id = "map")
      hide(selector = ".rowhide")
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
            error = function(e) {
              app.env$messages <- e$message
            }
          )
          setProgress(value = 1)
        }
      )
      initQC(app.env = app.env,
             refs = refs,
             session = session,
             output = output,
             googlesheet = googlesheet,
             mt.key = mt.key,
             mito.pattern = mito.pattern)
    }
  )
  observeEvent(eventExpr = input$triggerdemo,
               handlerExpr = {
                 # In case this is a re-upload, reset some elements
                 app.env$messages <- ""
                 output$valuebox.upload <- NULL
                 output$valuebox.preproc <- NULL
                 output$valuebox.mapped <- NULL
                 output$menu2 <- NULL
                 disable(id = 'map')
                 disable(id = 'triggerdemo')
                 hide(selector = ".rowhide")
                 withProgress(
                   message = "Reading input",
                   expr = {
                     setProgress(value = 0)
                     tryCatch(
                       expr = {
                         demodataset <- getOption(x = "Azimuth.app.demodataset")
                         app.env$object <- LoadFileInput(path = demodataset)
                         app.env$object$query <- 'query'
                         Idents(app.env$object) <- 'query'
                       },
                       error = function(e) {
                         app.env$messages <- e$message
                       }
                     )
                     setProgress(value = 1)
                   }
                 )
                 enable(id = 'triggerdemo')
                 initQC(app.env = app.env,
                        refs = refs,
                        session = session,
                        output = output,
                        googlesheet = googlesheet,
                        mt.key = mt.key,
                        mito.pattern = mito.pattern)
                 hide(selector = ".rowhide")
                 }
               )
  observeEvent( # Map data
    eventExpr = input$map,
    handlerExpr = {
      maptime.start <- Sys.time()

      disable(id = 'map')
      disable(id = 'num.ncountmin')
      disable(id = 'num.ncountmax')
      disable(id = 'num.nfeaturemin')
      disable(id = 'num.nfeaturemax')
      disable(id = 'num.mtmax')
      disable(id = 'num.mtmin')
      disable(id = 'check.qcscale')
      disable(id = 'check.qcpoints')
      # Run SCTransform and enable mapping
      withProgress(
        message = "Normalizing with SCTransform",
        expr = {
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
          # not enough cells available after filtering: reset filter elements
          if (ncellspreproc < getOption(x = "Azimuth.map.ncells")) {
            output$valuebox.preproc <- renderValueBox(expr = valueBox(
              value = ncellspreproc,
              subtitle = "cells after filtering",
              icon = icon("times"),
              color = "red"
            ))
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
                value = floor(min(mito.val)),
                min = floor(min(mito.val)),
                max = max(mito.val)
              )
              enable(id = 'num.mtmin')
              updateNumericInput(
                session = session,
                inputId = 'num.mtmax',
                label = paste("max", mt.key),
                value = ceiling(max(mito.val)),
                min = min(mito.val),
                max = ceiling(max(mito.val))
              )
              enable(id = 'num.mtmax')
              enable(id = 'check.qcscale')
              enable(id = 'check.qcpoints')
              enable(id = 'map')
            }
          } else {
            output$valuebox.preproc <- renderValueBox(expr = valueBox(
              value = ncellspreproc,
              subtitle = "cells after filtering",
              icon = icon("check"),
              color = "green"
            ))
            if (!is.null(googlesheet)) {
              try(sheet_append(ss = googlesheet, data = data.frame("CELLSPREPROC", Sys.info()[["nodename"]], ncellspreproc)))
            }
            app.env$object <- app.env$object[, cells.use]
            # Pseudobulk correlation test
            pbcor <- PBCorTest(
              object = app.env$object,
              ref = refs$avg
            )
            if (!is.null(googlesheet)) {
              try(sheet_append(ss = googlesheet, data = data.frame("PBCOR", Sys.info()[["nodename"]], pbcor[["cor.res"]])))
            }
            if (pbcor[["cor.res"]] < getOption(x = 'Azimuth.map.pbcorthresh')) {
              output$valuebox.mapped <- renderValueBox(expr = valueBox(
                value = "Failure",
                subtitle = "Query is too dissimilar",
                icon = icon("times"),
                color = "red", width = 6
              ))
              output$plot.pbcor <- renderPlot(pbcor[["plot"]])
              show(selector = ".rowhide")
              app.env$object <- NULL
              gc(verbose = FALSE)
            } else {
              setProgress(value = 0.1, message = "Normalizing with SCTransform")
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
                paste(ncellspreproc, "cells preprocessed")
              )
              setProgress(value = 0.3, message = "Finding anchors")
              cells <- colnames(x = app.env$object)
              anchors <- FindTransferAnchors(
                reference = refs$map,
                query = app.env$object,
                k.filter = NA,
                reference.neighbors = "spca.annoy.neighbors",
                reference.assay = "SCT",
                query.assay = 'SCT',
                reference.reduction = 'spca',
                normalization.method = 'SCT',
                features = intersect(rownames(x = refs$map), VariableFeatures(object = app.env$object)),
                dims = 1:50,
                n.trees = n.trees,
                verbose = TRUE,
                mapping.score.k = 100
              )
              # TODO fail if not enough anchors (Azimuth.map.nanchors)
              setProgress(value = 0.5, message = 'Mapping cells')
              app.env$object <- TransferData(
                reference = refs$map,
                query = app.env$object,
                dims = 1:50,
                anchorset = anchors,
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
                anchorset = anchors,
                reference = refs$map,
                query = app.env$object,
                reductions = "pcaproject",
                reuse.weights.matrix = TRUE
              )
              setProgress(value = 0.7, message = "Calculating mapping score")
              spca <- subset(
                x = anchors@object.list[[1]][["pcaproject.l2"]],
                cells = paste0(Cells(x = app.env$object), "_query")
              )
              spca <- RenameCells(object = spca, new.names = Cells(x = app.env$object))
              spca.ref <- subset(
                x = anchors@object.list[[1]][["pcaproject.l2"]],
                cells = paste0(Cells(x = refs$map), "_reference")
              )
              spca.ref <- RenameCells(object = spca.ref, new.names = Cells(x = refs$map))
              if (Sys.getenv("RSTUDIO") == "1") {
                plan("sequential")
              }
              # reduce size of object in anchorset
              anchors@object.list[[1]] <- DietSeurat(object = anchors@object.list[[1]])
              anchors@object.list[[1]] <- subset(
                x = anchors@object.list[[1]],
                features = c(rownames(x = anchors@object.list[[1]])[1])
              )
              anchors@object.list[[1]] <- RenameCells(
                object = anchors@object.list[[1]],
                new.names = unname(obj = sapply(
                  X = Cells(x = anchors@object.list[[1]]),
                  FUN = function(x) {
                    return(gsub(pattern = "_reference", replacement = "", x = x))
                  }
                )))
              anchors@object.list[[1]] <- RenameCells(
                object = anchors@object.list[[1]],
                new.names = unname(obj = sapply(
                  X = Cells(x = anchors@object.list[[1]]),
                  FUN = function(x) {
                    return(gsub(pattern = "_query", replacement = "", x = x))
                  }
                )))
              anchors@object.list[[1]]@meta.data <- data.frame()
              anchors@object.list[[1]]@active.ident <- factor()
              app.env$mapping.score <- future(
                expr = {
                  MappingScore(
                    anchors = anchors@anchors,
                    combined.object = anchors@object.list[[1]],
                    query.neighbors =  slot(object = anchors, name = "neighbors")[["query.neighbors"]],
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
              rm(anchors, spca)
              gc(verbose = FALSE)
              setProgress(value = 0.8, message = "Running UMAP transform")
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
                paste(ncellspreproc, "cells mapped")
              )
              output$valuebox.mapped <- renderValueBox(expr = valueBox(
                value = "Success",
                subtitle = "Mapping complete",
                icon = icon("check"),
                color = "green"
              ))
              # Enable the feature explorer
              enable(id = 'feature')
              app.env$default.feature <- ifelse(
                test = getOption(x = 'Azimuth.app.default_gene') %in% rownames(x = app.env$object),
                yes = getOption(x = 'Azimuth.app.default_gene'),
                no = VariableFeatures(object = app.env$object)[1]
              )
              app.env$features <- FilterFeatures(
                features = rownames(x = app.env$object)
              )
              shiny::updateSelectizeInput(
                session = session,
                inputId = 'feature',
                label = 'Feature',
                choices = app.env$features,
                selected = app.env$default.feature,
                server = TRUE,
                options = selectize.opts
              )
              # Add the predicted ID and score to the plots
              enable(id = 'adtfeature')
              app.env$adt.features <- sort(x = FilterFeatures(features = rownames(
                x = app.env$object[[adt.key]]
              )))
              app.env$default.adt <- ifelse(
                test = getOption(x = 'Azimuth.app.default_adt') %in% app.env$adt.features,
                yes = getOption(x = 'Azimuth.app.default_adt'),
                no = app.env$adt.features[1]
              )
              updateSelectizeInput(
                session = session,
                inputId = 'adtfeature',
                choices = app.env$adt.features,
                selected = '',
                server = TRUE,
                options = selectize.opts
              )
              metadata.choices <- sort(x = c(
                "predicted.id",
                PlottableMetadataNames(
                  object = app.env$object,
                  min.levels = 1,
                  max.levels = 50
                )
              ))
              updateSelectizeInput(
                session = session,
                inputId = 'select.metadata',
                choices = metadata.choices,
                selected = 'predicted.id',
                server = TRUE,
                options = selectize.opts
              )
              updateSelectizeInput(
                session = session,
                inputId = 'select.metadata1',
                choices = metadata.choices,
                selected = 'predicted.id',
                server = TRUE,
                options = selectize.opts
              )
              updateSelectizeInput(
                session = session,
                inputId = 'select.metadata2',
                choices = metadata.choices,
                selected = 'predicted.id',
                server = TRUE,
                options = selectize.opts
              )
              enable(id = 'select.metadata')
              enable(id = 'select.metadata1')
              enable(id = 'select.metadata2')
              enable(id = 'radio.pct')
              # Enable continuous metadata
              metadata.cont <- sort(x = setdiff(
                x = colnames(x = app.env$object[[]]),
                y = metadata.choices
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
              enable(id = 'scorefeature')
              updateSelectizeInput(
                session = session,
                inputId = 'groupfeature',
                choices = metadata.choices,
                selected = 'predicted.id',
                server = TRUE,
                options = selectize.opts
              )
              enable(id = 'groupfeature')
              # Compute biomarkers
              setProgress(
                value = 0.95,
              message = "Running differential expression"
              )
              app.env$diff.expr[[app.env$default.assay]] <- wilcoxauc(
                X = app.env$object,
                group_by = 'predicted.id',
                assay = 'data',
                seurat_assay = app.env$default.assay
              )
              app.env$diff.expr[[adt.key]] <- wilcoxauc(
                X = app.env$object,
                group_by = 'predicted.id',
                assay = 'data',
                seurat_assay = adt.key
              )
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
              enable(id = 'select.prediction')
              updateSelectizeInput(
                session = session,
                inputId = 'select.prediction',
                choices = allowed.clusters,
                selected = allowed.clusters[1],
                server = TRUE,
                options = selectize.opts
              )
              enable(id = 'select.biomarkers')
              updateSelectizeInput(
                session = session,
                inputId = 'select.biomarkers',
                choices = allowed.clusters,
                selected = allowed.clusters[1],
                server = TRUE,
                options = selectize.opts
              )
              # Enable downloads
              enable(id = 'dlumap')
              enable(id = 'dladt')
              enable(id = 'dlpred')
              enable(id = 'dlscript')
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
              maptime.diff <- difftime(
                time1 = Sys.time(),
                time2 = maptime.start,
                units = "secs"
              )
              if (maptime.diff < 60) {
                time.fmt <- gsub(
                  pattern = " 0",
                  replacement = " ",
                  x = format(x = .POSIXct(xx = maptime.diff), format = "in %S seconds"),
                  fixed = TRUE
                )
              } else {
                time.fmt <- gsub(
                  pattern = " 0",
                  replacement = " ",
                  x = format(
                    x = .POSIXct(xx = maptime.diff),
                    format = "in %M minutes %S seconds"
                  ),
                  fixed = TRUE
                )
              }
              app.env$messages <- c(
                app.env$messages,
                time.fmt
              )
              if (!is.null(googlesheet)) {
                try(sheet_append(ss = googlesheet, data = data.frame("MAPPINGTIME", Sys.info()[["nodename"]], as.numeric(maptime.diff, units="secs"))))
              }
            }
          }
          setProgress(value = 1)
        }
      )
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
          groups.use = input$select.biomarkers,
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
          groups.use = input$select.biomarkers,
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
            groups.use = input$select.biomarkers,
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
            groups.use = input$select.biomarkers,
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
      if (mt.key %in% colnames(x = isolate(app.env$object[[]]))) {
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
      if (mt.key %in% colnames(x = isolate(app.env$object[[]]))) {
        vlnlist[[3]] <- vlnlist[[3]] +
          geom_hline(yintercept = input$num.mtmin) +
          geom_hline(yintercept = input$num.mtmax) +
          annotate(
            geom = "rect",
            alpha = 0.2,
            fill = "red",
            ymin = input$num.mtmax,
            ymax = Inf,
            xmin = 0.5,
            xmax = 1.5
          ) +
          annotate(
            geom = "rect",
            alpha = 0.2,
            fill = "red",
            ymin = ifelse(test = input$check.qcscale, yes = 0, no = -Inf),
            ymax = input$num.mtmin,
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
        if (input$select.metadata == "predicted.id") {
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
            group.by = input$select.metadata,
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
            group.by = input$groupfeature,
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
  output$containerid <- renderUI(expr = {
    p(HTML(text = paste(
      paste("debug ID:", Sys.info()[["nodename"]]),
      paste('Azimuth version:', packageVersion(pkg = 'Azimuth')),
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
          isolate(app.env$object)[[mt.key, drop = TRUE]] >= input$num.mtmin &
          isolate(app.env$object)[[mt.key, drop = TRUE]] <= input$num.mtmax
      }
      paste(sum(cells.use), "cells remain after current filters")
    }
  })
  output$text.dladt <- renderText(
    expr = {
      c("imputed.assay <- readRDS('azimuth_impADT.Rds')",
        "object <- object[, Cells(imputed.assay)]",
        "object[['impADT']] <- imputed.assay")
    },
    sep = "\n"
  )
  output$text.dlumap <- renderText(
    expr = {
      c("projected.umap <- readRDS('azimuth_umap.Rds')",
        "object <- object[, Cells(projected.umap)]",
        "object[['umap.proj']] <- projected.umap")
    },
    sep = "\n"
  )
  output$text.dlpred <- renderText(
    expr = {
      c("predictions <- read.delim('azimuth_pred.tsv', row.names = 1)",
        "object <- AddMetaData(",
        "\tobject = object,",
        "\tmetadata = predictions)")
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
          groups.use = input$select.biomarkers,
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
          groups.use = input$select.biomarkers,
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
      e$mito.max <- input$num.mtmax
      e$mito.min <- input$num.mtmin
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
  output$welcomebox <- renderUI(expr = eval(parse(text = getOption(x = "Azimuth.app.welcomebox"))))
}

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
AzimuthApp <- function(config = NULL, ...) {
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
  opts$future.globals.maxSize <- maxcells * 320000
  # Launch the app
  with_options(
    new = opts,
    code = runApp(appDir = shinyApp(ui = ui, server = server))
  )
  return(invisible(x = NULL))
}

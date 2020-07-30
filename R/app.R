#' @include zzz.R
#' @include seurat.R
#' @import V8
#' @importFrom htmltools tagList h4 hr h3
#' @importFrom shinyjs useShinyjs extendShinyjs disabled
#' @importFrom shiny fluidPage sidebarLayout sidebarPanel fileInput sliderInput
#' actionButton selectInput downloadButton mainPanel tabsetPanel tabPanel
#' plotOutput tableOutput verbatimTextOutput
#'
NULL

app.title <- 'Azimuth'

ui <- tagList(
  useShinyjs(),
  extendShinyjs(
    text = TabJSHide(
      id = 'tabs',
      values = c('mapped', 'fexplorer', 'imputed', 'diffexp'),
      fxn = 'init'
    )
  ),
  fluidPage(
    title = app.title,
    sidebarLayout(
      sidebarPanel = sidebarPanel(
        # TODO: disable this initially? I don't think we can browse for data
        # while the reference is being loaded
        fileInput(
          inputId = "file",
          # TODO list supported filetypes in tooltip or helptext
          label = "File Upload (h5, h5seurat, or Seurat object as rds)",
          accept = c('.h5', '.h5seurat', '.rds')
        ),
        h4("Preprocessing Controls"),
        hr(),
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
        h4("Feature Selection"),
        hr(),
        disabled(selectInput(
          inputId = 'feature',
          label = 'Feature',
          choices = '',
          selectize = FALSE,
          width = '100%'
        )),
        disabled(selectInput(
          inputId = 'adtfeature',
          label = 'Imputed Protein',
          choices = '',
          selectize = FALSE,
          width = '100%'
        )),
        h4("Cluster Selection"),
        hr(),
        disabled(selectInput(
          inputId = 'declusters',
          label = 'Cell Type',
          choices = '',
          selectize = FALSE,
          width = '100%'
        )),
        h4("Downloads"),
        hr(),
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
          value = 'preprocessing',
          plotOutput(outputId = "qcvln"),
          tableOutput(outputId = 'qctbl'),
          verbatimTextOutput(outputId = "sct")
        ),
        tabPanel(
          title = 'Mapped Data',
          value = 'mapped',
          plotOutput(outputId = 'refdim'),
          plotOutput(outputId = 'objdim'),
          verbatimTextOutput(outputId = "mapping")
        ),
        tabPanel(
          title = 'Feature Explorer',
          value = 'fexplorer',
          plotOutput(outputId = 'fvln'),
          plotOutput(outputId = 'fdim'),
          hr(),
          h3("Imputed Proteins"),
          hr(),
          plotOutput(outputId = 'ivln'),
          plotOutput(outputId = 'idim')
        ),
        tabPanel(
          title = 'Biomarkers',
          value = 'diffexp',
          h3('Biomarkers'),
          tableOutput(outputId = 'biomarkers'),
          h3('ADT Biomarkers'),
          tableOutput(outputId = 'adtbio')
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
#' @importFrom methods slot<-
#' @importFrom ggplot2 ggtitle scale_colour_hue
#' @importFrom presto wilcoxauc
#' @importFrom shinyjs show enable disable
#' @importFrom Seurat DefaultAssay PercentageFeatureSet SCTransform
#' VariableFeatures Idents GetAssayData RunUMAP CreateAssayObject
#' CreateDimReducObject Embeddings AddMetaData SetAssayData Key
#' VlnPlot DimPlot Reductions FeaturePlot Assays NoLegend Idents<-
#' @importFrom shiny reactiveValues safeError appendTab observeEvent
#' withProgress setProgress updateSliderInput renderText updateSelectInput
#' updateTabsetPanel renderPlot renderTable downloadHandler
#'
#' @keywords internal
#'
server <- function(input, output, session) {
  mt.key <- 'percent.mt'
  adt.key <- 'impADT'
  app.env <- reactiveValues(
    object = NULL,
    default.assay = NULL,
    default.feature = NULL,
    default.adt = NULL,
    diff.exp = list()
  )
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
      # enable(id = 'file')
    }
  )
  shinyjs::show(selector = TabJSKey(id = 'tabs', values = 'mapped'))
  # React to events
  observeEvent( # Load the data
    eventExpr = input$file,
    handlerExpr = {
      withProgress(
        message = "Reading input",
        expr = {
          setProgress(value = 0)
          app.env$object <- LoadFileInput(path = input$file$datapath)
          app.env$object$query <- 'query'
          Idents(app.env$object) <- 'query'
          setProgress(value = 1)
        }
      )
      app.env$default.assay <- DefaultAssay(object = app.env$object)
      enable(id = 'ncount')
      enable(id = 'nfeature')
      ncount <- paste0('nCount_', app.env$default.assay)
      nfeature <- paste0('nFeature_', app.env$default.assay)
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
      mito.pattern <- getOption(x = 'Azimuth.app.mito', default = '^MT-')
      if (any(grepl(pattern = mito.pattern, x = rownames(x = app.env$object)))) {
        app.env$object <- PercentageFeatureSet(
          object = app.env$object,
          pattern = mito.pattern,
          col.name = mt.key,
          assay = app.env$default.assay
        )
      }
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
            residual.features = rownames(x = refs$map),
            ncells = getOption(x = 'Azimuth.sct.ncells'),
            n_genes = getOption(x = 'Azimuth.sct.nfeats'),
            do.correct.umi = FALSE,
            do.scale = FALSE,
            do.center = TRUE
          ))
          setProgress(value = 1)
          output$sct <- renderText(expr = "SCTransform complete")
        }
      )
      enable(id = "map")
      disable(id = 'proc1')
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
      shinyjs::show(selector = TabJSKey(id = 'tabs', values = 'fexplorer'))
    }
  )
  observeEvent( # Map data
    eventExpr = input$map,
    handlerExpr = {
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
            assay.query = app.env$default.assay
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
      if (sum(app.env$object$mapped) * 100 < getOption(x = "Azimuth.map.pcthresh")) {
        stop(safeError(error = "Query dataset could not be mapped to the reference"))
      }
      mappingtext <- paste0(
        sum(app.env$object$mapped),
        " cells mapped to reference (",
        round(
          x = sum(app.env$object$mapped) / ncol(x = app.env$object) * 100,
          digits = 2
        ),
        "%)"
      )
      output$mapping <- renderText(expr = mappingtext)
      app.env$object <- app.env$object[, app.env$object$mapped]
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
        selected = app.env$default.adt
      )
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
        x = table(app.env$object$predicted.id) > getOption(x = 'Azimuth.de.mincells'),
      ))
      allowed.clusters <- factor(
        x = allowed.clusters,
        levels = levels(x = app.env$object)
      )
      allowed.clusters <- levels(x = droplevels(x = allowed.clusters))
      enable(id = 'declusters')
      updateSelectInput(
        session = session,
        inputId = 'declusters',
        label = 'Cell Type',
        choices = allowed.clusters,
        selected = allowed.clusters[1]
      )
      shinyjs::show(selector = TabJSKey(id = 'tabs', 'diffexp'))
      # Enable downloads and downstream analyses
      updateTabsetPanel(
        session = session,
        inputId = 'tabs',
        selected = 'mapped'
      )
      enable(id = 'dlumap')
      enable(id = 'dlpred')
      disable(id = 'map')
    }
  )
  # Plots
  output$qcvln <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      qc <- paste0(c('nCount_', 'nFeature_'), app.env$default.assay)
      if (mt.key %in% colnames(x = app.env$object[[]])) {
        qc <- c(qc, mt.key)
      }
      VlnPlot(object = app.env$object, features = qc, group.by = 'query')
    }
  })
  output$refdim <- renderPlot(expr = {
    DimPlot(object = refs$plot) + ggtitle(label = 'Reference')
  })
  output$objdim <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      if (length(x = Reductions(object = app.env$object))) {
        DimPlot(object = app.env$object) +
          scale_colour_hue(limits = levels(refs$plot$id), drop = FALSE) +
          ggtitle(label = 'Query')
      }
    }
  })
  output$fvln <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      avail <- c(rownames(x = app.env$object), colnames(x = app.env$object[[]]))
      if (input$feature %in% avail) {
        VlnPlot(object = app.env$object, features = input$feature) +
          NoLegend()
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
  output$ivln <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      if (adt.key %in% Assays(object = app.env$object)) {
        VlnPlot(
          object = app.env$object,
          features = paste0(
            Key(object = app.env$object[[adt.key]]),
            input$adtfeature
          )
        ) +
          NoLegend()
      }
    }
  })
  output$idim <- renderPlot(expr = {
    if (!is.null(x = app.env$object)) {
      if (adt.key %in% Assays(object = app.env$object)) {
        FeaturePlot(
          object = app.env$object,
          features = paste0(
            Key(object = app.env$object[[adt.key]]),
            input$adtfeature
          ),
          cols = c('lightgrey', 'darkred'),
          min.cutoff = 'q10',
          max.cutoff = 'q99'
        )
      }
    }
  })
  # Tables
  output$qctbl <- renderTable(
    expr = {
      if (!is.null(x = app.env$object)) {
        qc <- paste0(c('nCount_', 'nFeature_'), app.env$default.assay)
        tbl <- apply(X = app.env$object[[qc]], MARGIN = 2, FUN = quantile)
        tbl <- as.data.frame(x = tbl)
        colnames(x = tbl) <- c('nUMI per cell', 'Genes detected per cell')
        if (mt.key %in% colnames(x = app.env$object[[]])) {
          tbl[, 3] <- quantile(x = app.env$object[[mt.key, drop = TRUE]])
          colnames(x = tbl)[3] <- 'Mitochondrial percentage per cell'
        }
        t(x = tbl)
      }
    },
    rownames = TRUE
  )
  output$biomarkers <- renderTable(
    expr = {
      if (!is.null(x = app.env$diff.expr[[app.env$default.assay]])) {
        RenderDiffExp(
          diff.exp = app.env$diff.expr[[app.env$default.assay]],
          groups.use = input$declusters
        )
      }
    },
    rownames = TRUE
  )
  output$adtbio <- renderTable(
    expr = {
      if (!is.null(x = app.env$diff.expr[[adt.key]])) {
        RenderDiffExp(
          diff.exp = app.env$diff.expr[[adt.key]],
          groups.use = input$declusters
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
#' @param reference,mito,max.upload.mb See \strong{App options} for more details
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
#'  \item{\code{Azimuth.app.default.gene}}{
#'   Gene to select by default in Feature Explorer
#'  }
#'  \item{\code{Azimuth.app.default.adt}}{
#'   ADT to select by default in Feature Explorer
#'  }
#' }
#'
#' @return None, launches the mapping Shiny app
#'
#' @importFrom shiny runApp shinyApp
#'
#' @export
#'
#' @seealso \code{\link{SeruatMapper-package}}
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
    Azimuth.app.default.gene = default.gene,
    Azimuth.app.default.adt = default.adt
  )
  withr::with_options(
    new = opts,
    code = runApp(appDir = shinyApp(ui = ui, server = server))
  )
  return(invisible(x = NULL))
}

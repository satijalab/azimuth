
library(devtools)
library(shinyjs)
library(shiny)
library(uwot)
library(sctransform)
library(ggplot2)
load_all("~/seurat/")
options(shiny.maxRequestSize = 100 * (1024 ^ 2))
# example app for prepending/appending a navbarMenu

reference <- readRDS("~/demo/reference_pbmc.rds")
reference.nn <- readAnnoyNN(file = "~/demo/reference.nn_pbmc.rds")
plotref <- readRDS("~/demo/plotref_pbmc.rds")
adtref <- readRDS("/home/satijar/demo/adtreference_pbmc.rds")
reference[["ADT"]] <- adtref[["ADT"]]

ui <-tagList(useShinyjs(),  fluidPage(title = "HuBMAPer",  
                                      sidebarLayout(sidebarPanel = sidebarPanel(
                                        fileInput(inputId = "tenxh5", label = "Input file"),
                                        selectInput(inputId = "objtype", label = "object type",selected='rds',
                                                    choices = c( "10X H5DF" = "10XH5", 
                                                                 "10X"="10X",
                                                                 "rds"= "rds",
                                                                 "h5Seurat" = "h5Seurat")),
                                        textInput(inputId = "name", label = "query object name", value = "query"), 
                                        disabled( actionButton("proc1", "Preprocess input") ),
                                        disabled(actionButton("map", "Map cells to reference"))
                                      ), 
                                      mainPanel = mainPanel( tabsetPanel(id = "tabs", 
                                                                         tabPanel(title = "File Upload", 
                                                                                  verbatimTextOutput(outputId = "sct"), 
                                                                                  verbatimTextOutput(outputId = "mapping")) 
                                      )
                                      ))))
server <- function(input, output, session) {
  plot.env <- shiny::reactiveValues()
  
  
  ## read in file
  observeEvent(
    eventExpr = input$tenxh5,
    handlerExpr = {
      enable("proc1")
    }
  )
  
  ## pre-processing
  observeEvent(
    eventExpr =  input$proc1 ,
    handlerExpr = {
      withProgress(message = 'Reading input', {
        setProgress(value = 0)
        setProgress(value = 0.2, message = 'Creating Seurat object')
        object <- switch(
          EXPR = input$objtype, 
          '10XH5' = { 
            mat <- Seurat::Read10X_h5(input$tenxh5$datapath)
            Seurat::CreateSeuratObject(counts = mat)
          }, 
          '10X' = {
            Read10X(input$tenxh5$datapath)
          }, 
          'rds' = {readRDS( input$tenxh5$datapath )}, 
          #'rds' = {readRDS("~/data/dropseq_pbmc.rds")}, 
          'h5Seurat'= { LoadH5Seurat(input$tenxh5$datapath  )}
        )
        object <- subset(object,nFeature_RNA>300)
        setProgress(value = 0.4,message = 'Normalization (SCTransform)')
        suppressWarnings(object <- SCTransform(object,
                                               residual.features  = rownames(reference),
                                               ncells = 1000,n_genes=1000))
        setProgress(value = 1,message = 'Finished pre-processing')
        
        output$qcViolin <- shiny::renderPlot(expr = Seurat::VlnPlot(object,'nFeature_RNA'))
        
        # return SCT information to main page
        output$sct <- renderText({"SCTransform is complete"})
        enable("map")
        disable("proc1")
        insertTab(inputId = "tabs", position = 'after',
                  tabPanel("QC", "We should add a table of QC stats here",
                           shiny::plotOutput(outputId = "qcViolin")),
                  target = "File Upload"
        )
      }) 
      plot.env$user <- object
    }
  )
  
  ##mapping
  observeEvent(
    eventExpr = input$map,
    handlerExpr = {
      withProgress(message = 'Mapping', {
        setProgress(value = 0,message = 'Mapping - Finding anchors')
        spca.anchor <- FindTransferAnchors_Fast(reference = reference, 
                                                query = plot.env$user,
                                                reference.nnidx = reference.nn$annoy_index,
                                                reference.nn = reference.nn,
                                                reference.assay = "RNA",
                                                npcs = NULL, 
                                                k.filter = NA,
                                                query.assay = "SCT", 
                                                reference.reduction = "spca", 
                                                normalization.method = "SCT",
                                                features = rownames(reference), 
                                                dims = 1:50)
        setProgress(value = 0.6,message = 'Integrating data')
        ingest.spca <- IngestNewData_Fast(reference = reference, 
                                          query = plot.env$user,
                                          dims = 1:50,k.weight = 25,
                                          transfer.anchors = spca.anchor,
                                          reference.nnidx = reference.nn$annoy_index,
                                          transfer.labels = reference$id,transfer.expression = reference[["ADT"]]@data)
        setProgress(value = 0.8, message = 'Running UMAP transform')
        
        plot.env$user$predicted.id <- ingest.spca$predicted.id[ Cells(plot.env$user) ]
        plot.env$user$predicted.id.score <- ingest.spca$predicted.id.score[ Cells(plot.env$user)  ]
        
        nn.idx <- ingest.spca@neighbors$query_ref.nn$nn.idx 
        ori.nn.idx<- t(sapply( 1:nrow(nn.idx) , function(x) reference@meta.data[ nn.idx[x,] , "ori.index"]))
        rownames(ori.nn.idx) <- rownames(nn.idx)
        ingest.spca@neighbors$query_ref.nn$nn.idx  <- ori.nn.idx
        
        
        plot.env$user[["umap"]] <- RunUMAP(object = ingest.spca[["query_ref.nn"]],
                                                reduction.model =  reference[["jumap"]], 
                                                reduction.key = "UMAP_")
      
        plot.env$user$predicted.id <- factor(plot.env$user$predicted.id,levels=levels(reference$id))
        output$refDimPlot <- shiny::renderPlot(expr = Seurat::DimPlot(plotref,group.by = 'id',label = F))
        mapPlot <- Seurat::DimPlot(plot.env$user,group.by = 'predicted.id',label = F)+scale_colour_hue(limits=levels(reference$id),drop=FALSE)
        output$mapDimPlot <- shiny::renderPlot(expr = mapPlot)
        
        plot.env$user[["ADT"]] <- CreateAssayObject( data = ingest.spca[["transfer"]]@data[,Cells(plot.env$user)])
        
        # provide download corrected UMAP and meta.data
        output$umap <- downloadHandler(filename = function(){
          paste0(input$name, "_umap.rds") } ,
          content = function(file){ 
            saveRDS(plot.env$user[["umap"]], file = file )
          } )
        output$metadata <- downloadHandler(filename =  function(){
          paste0(input$name, "_metadata.csv")} ,
          content = function(file){
            Seurat::write.table(x = plot.env$user@meta.data[,c("predicted.id", "predicted.id.score")], file = file, sep = "," )  })
        
        
        # return mapping information to main page
        output$mapping <- renderText({"Mapping is complete"})
        disable("map")
        setProgress(value = 1,message = 'Done mapping')
        insertTab(inputId = "tabs", position = 'after',
                  tabPanel("Mapped",
                           shiny::plotOutput(outputId = "refDimPlot"), 
                           shiny::plotOutput(outputId = "mapDimPlot"), 
                           downloadButton(outputId = "umap", label = "Download the umap"),
                           downloadButton(outputId = "metadata",
                                          label = "Download the predicted.id and predicted.score")
                  ),
                  target = "QC"
        )
        
        #grep removes useless features, otherwise too many and you get weird errors starting with T
        features <- sort(grep("\\.1",rownames(plot.env$user[["RNA"]]), value = T, invert = T))
        features_adt <- paste0('adt_',rownames(reference[["ADT"]]))
        insertTab(inputId = "tabs",position = 'after',
                  tabPanel("Gene Expression",  
                           sidebarLayout(
                             sidebarPanel(width = 3, 
                                          selectInput(
                                            inputId = 'feature',
                                            label = 'Gene',
                                            choices = c(sort(features),sort(features_adt),'predicted.id.score','nFeature_RNA'),
                                            selected = "GNLY",
                                            selectize = FALSE,
                                            width = '100%'
                                          )
                             ),
                             mainPanel(
                               shiny::plotOutput(outputId = "featureGene"),
                               shiny::plotOutput(outputId = "violinGene")
                             )
                           )
                  ),
                  target = "Mapped"
        )
      }
      ) 
    }
  )
  
  ##gene expression
  output$violinGene <- shiny::renderPlot(expr = Seurat::VlnPlot(plot.env$user, input$feature, group.by = "predicted.id")+NoLegend())
  output$featureGene <- shiny::renderPlot(expr = Seurat::FeaturePlot(plot.env$user,input$feature,pt.size = 0.1,reduction = 'umap'))
  
}
shinyApp(ui, server)

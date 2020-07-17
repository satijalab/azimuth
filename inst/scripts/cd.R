# TODO can you load code using load_all() within the shiny app?
# TODO should Seurat functions be prefixed with "Seurat::" as in the existing code in app.R?
# TODO should messages be printed with the "message" function?
# TODO do we want messages at all?
# TODO rendering DimPlot?

# Load code from seurat-private (branch feat/FastIngest)
library(devtools)
load_all("/home/darbyc/seurat-private") #TODO update path to Github repo
# Load UWOT for UMAP algorithm
load_all("/home/haoy/uwot/") #TODO update path to Github repo
# Load updated sctransform
load_all("/home/haoy/sctransform/") #TODO update path to Github repo

starttime.global <- Sys.time()
starttime.step <- starttime.global
message("Loading reference data from disk.")

# Load reference object from RDS
reference <- readRDS("/home/haoy/fast-mapping/pbmc_ref.diet.rds") #TODO update file path
# Load reference neighbors + Annoy index from RDS
reference.nnidx <- Seurat::readAnnoyNN("/home/haoy/fast-mapping/pbmc_ref.nn.rds") #TODO update file path

endtime.step <- Sys.time()
message("Reference load step: ", round(as.numeric(endtime.step - starttime.step, units="secs"), digits=1), "s")
message("Total time so far: ", round(as.numeric(endtime.step - starttime.global, units="secs"), digits=1), "s")


starttime.step <- Sys.time()
message("Loading query matrix.")

# Load query matrix (pbmc3k)
# mat <- Seurat::Read10X_h5(input$tenxh5$datapath) #TODO update path variable
# if (is.list(x = mat)) {
#   mat <- mat[[1]]
# }
# query <- Seurat::CreateSeuratObject(counts = mat)

# Alternate route from SeuratData
library(SeuratData)
data("pbmc3k")
DefaultAssay(pbmc3k) <- "RNA"
query <- pbmc3k

endtime.step <- Sys.time()
message("Query load step: ", round(as.numeric(endtime.step - starttime.step, units="secs"), digits=1), "s")
message("Total time so far: ", round(as.numeric(endtime.step - starttime.global, units="secs"), digits=1), "s")


starttime.step <- Sys.time()
message("Normalizing query matrix with SCTransform.")

# Normalize query matrix with SCTransform
query <- Seurat::SCTransform(
  query, 
  residual.features = rownames(reference))

endtime.step <- Sys.time()
message("Normalize step: ", round(as.numeric(endtime.step - starttime.step, units="secs"), digits=1), "s")
message("Total time so far: ", round(as.numeric(endtime.step - starttime.global, units="secs"), digits=1), "s")


starttime.step <- Sys.time()
message("Finding transfer anchors.")

# Find transfer anchors
anchors <- Seurat::FindTransferAnchors_Fast(
  reference = reference, 
  query = query, 
  reference.nn = reference.nnidx,
  reference.assay = "RNA",
  reference.nnidx = reference.nnidx$annoy_index,
  npcs = NULL, 
  k.filter = NA,
  query.assay = "SCT", 
  reference.reduction = "spca", 
  normalization.method = "SCT",
  features = rownames(reference[["spca"]]@feature.loadings), 
  dims = 1:50)

endtime.step <- Sys.time()
message("Anchor step: ", round(as.numeric(endtime.step - starttime.step, units="secs"), digits=1), "s")
message("Total time so far: ", round(as.numeric(endtime.step - starttime.global, units="secs"), digits=1), "s")


starttime.step <- Sys.time()
message("Correcting query embeddings and transferring features.")

# Ingest
ingest <- Seurat::IngestNewData_Fast(
  reference = reference, 
  query = query,
  dims = 1:50,
  transfer.anchors = anchors, 
  transfer.labels = reference$id,
  reference.nnidx = reference.nnidx$annoy_index)
query$predicted.id <- ingest$predicted.id[ Cells(query) ]
query$predicted.id.score <- ingest$predicted.id.score[ Cells(query)  ]

endtime.step <- Sys.time()
message("Ingest step: ", round(as.numeric(endtime.step - starttime.step, units="secs"), digits=1), "s")
message("Total time so far: ", round(as.numeric(endtime.step - starttime.global, units="secs"), digits=1), "s")


starttime.step <- Sys.time()
message("Computing query UMAP projection.")

# UMAP
query[["umap.proj"]] <- Seurat::RunUMAP(
  object = ingest[["query_ref.nn"]], 
  reduction.model = reference[["jumap"]], 
  reduction.key = "ProjU_")

endtime.step <- Sys.time()
message("UMAP step: ", round(as.numeric(endtime.step - starttime.step, units="secs"), digits=1), "s")
message("Total time so far: ", round(as.numeric(endtime.step - starttime.global, units="secs"), digits=1), "s")

# Plot UMAP
Seurat::DimPlot( #TODO need to change to render this in the app
  query, 
  reduction = "umap.proj", 
  group.by = "predicted.id", 
  label = T) + 
  NoLegend()


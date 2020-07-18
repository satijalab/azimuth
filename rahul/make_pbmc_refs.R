#make Reference and Annoy refs
object.ref.ori <- readRDS("/home/haoy/fast-mapping/pbmc_ref.diet.rds")
#  downsample and build downsample ref.nn
Idents(object.ref.ori) <- 'cluster'
object.ref <- subset(object.ref.ori,  downsample=500)
newcells <- unique(c(Cells(object.ref),WhichCells(object.ref.ori,expression = id%in%c('CD8 Memory','CD4 Memory'))))
object.ref <- subset(object.ref.ori,cells = newcells)
object.ref <- subset(object.ref, id=='Doublet',invert=TRUE)

object.ref$ori.index <-  match( Cells(object.ref), Cells(object.ref.ori) )
ref.nn <- AnnoyNN(data = object.ref[["spca"]]@cell.embeddings[,1:50], metric = "cosine", k = 31, return.annoy_index = T )

reference <- object.ref
reference.nn <- ref.nn
saveRDS(reference,file = "/home/satijar/demo/reference_pbmc.rds")
saveAnnoyNN(nn = reference.nn,file = "/home/satijar/demo/reference.nn_pbmc.rds")



#make plotting refs
Idents(object.ref.ori) <- 'cluster'
o1 <- subset(object.ref.ori,  downsample=500)
Idents(object.ref.ori) <- 'id'
o2 <- subset(object.ref.ori,  downsample=2000,idents = c("CD14 Mono","CD4 Naive"))
plotref <- subset(object.ref.ori,cells = unique(c(Cells(o1),Cells(o2))))
plotref <- subset(plotref, id=='Doublet',invert=TRUE)

plotref[["spca"]] <- NULL
plotref[["umap"]] <- CreateDimReducObject(embeddings = plotref[["jumap"]]@cell.embeddings,key = "UMAP_")
plotref[["jumap"]] <- NULL
saveRDS(plotref,file = "/home/satijar/demo/plotref_pbmc.rds")

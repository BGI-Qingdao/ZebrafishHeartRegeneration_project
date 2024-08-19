library(argparse)
parser <- ArgumentParser()
parser$add_argument("-i", "--input", dest = "input", type="character",  required=TRUE)
parser$add_argument("-o", "--output",dest = "output",type="character",  required=TRUE)
parser$add_argument("-m", "--trim",  dest = "trim",  type="double",default=0.1)
opts = parser$parse_args()

trim = opts$trim
prefix = opts$output

library(Seurat)
library(CellChat)
###############################
# read RDS now ...
###############################
seurat_objects <- readRDS(opts$input)
obj.list <- SplitObject(seurat_objects,split.by='orig.ident')

for ( slice in  names(obj.list) ) {
       print(slice)
       seurat_object = obj.list[[slice]]
       rdsname = paste0(prefix,".",slice,".cellchat.rds")
       if( file.exists(rdsname) ) {
           cellchat = readRDS(rdsname)
       } else {
           ###############################
           # prepare data now
           ###############################
           DefaultAssay(seurat_object) = 'Spatial'

           seurat_object<- NormalizeData(object = seurat_object,
                                         normalization.method = "LogNormalize",
                                         scale.factor = 1000)

           seurat_object$row = seurat_object$grid_x/70
           seurat_object$col = (seurat_object$grid_y*-1)/70

           data.input <- GetAssayData(seurat_object,
                                      assay = "Spatial",
                                      slot = "data") # normalized data matrix
           labels <- seurat_object$celltype
           meta <- data.frame(group = labels,
                              row.names = names(labels)) # create a dataframe of the cell labels
           spatial.loc <- seurat_object@meta.data[,c('col','row')]
           scale.factors <- list(spot.diameter = 35 , spot = 1)
           ###############################
           # create cellchat object
           ###############################
           cellchat <- createCellChat(object = data.input,
                                      meta = meta,
                                      group.by = "group",
                                      datatype = "spatial",
                                      coordinates = spatial.loc,
                                      scale.factors = scale.factors)
           #cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
           CellChatDB <- CellChatDB.zebrafish
           CellChatDB.use <- CellChatDB
           cellchat@DB <- CellChatDB.use

           future::plan("multiprocess", workers = 15)
           cellchat <- subsetData(cellchat)
           cellchat <- identifyOverExpressedGenes(cellchat)
           cellchat <- identifyOverExpressedInteractions(cellchat)
           cellchat <- computeCommunProb(cellchat,
                                         type = "truncatedMean",
                                         trim = trim,
                                         distance.use = TRUE,
                                         interaction.range = 250,
                                         scale.distance = 0.025)
           cellchat <- filterCommunication(cellchat,
                                           min.cells = 5)
           cellchat <- computeCommunProbPathway(cellchat)
           cellchat <- aggregateNet(cellchat)
           cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

           ###########################################################
           # save rds
           ###########################################################
           saveRDS(cellchat, file = rdsname)
       }
       ###########################################################
       # save pdf graphics
       ###########################################################

       groupSize <- as.numeric(table(cellchat@idents))
       pdf(paste0(prefix,".",slice,'.circle_all.pdf'))
       netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
       netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

       mat <- cellchat@net$weight
       for (i in 1:nrow(mat)) {
           mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
           mat2[i, ] <- mat[i, ]
           netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
       }
       dev.off()

       pdf(paste0(prefix,".",slice,'.out.pdf'))
       print(netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing"))
       dev.off()
       pdf(paste0(prefix,".",slice,'.in.pdf'))
       print(netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming"))
       dev.off()
       pdf(paste0(prefix,'.',slice,'.spatial.pathway.pdf'))
       for( path in cellchat@netP$pathways ) {
            tryCatch(
               {
                    print(netVisual_aggregate(cellchat, layout= "spatial",signaling =c(path), edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5))
               },
               error = function(e){
               },
               warning = function(w){
               },
               finally = {
               }
            )
       }
       dev.off()
}

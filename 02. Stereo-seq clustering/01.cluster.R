## Get the parameters
parser = argparse::ArgumentParser(description = "spatial")

parser$add_argument('-i', dest = 'input', help = 'input directory')
parser$add_argument('-o', dest = 'out', help = 'out directory')
parser$add_argument('-n', dest = 'name', help = 'sample name')

opts = parser$parse_args()

##lib
library(Seurat)
library(dplyr)
library(stringr)
library(harmony)

##Merge
setwd(opts$input)
rds_list <- list.files(opts$input, pattern = '\\.rds$')

vf <- vector()
merge.list <- list()

for (index in 1:length(rds_list)){
	
	obj <- readRDS(rds_list[index])
	
	obj <- subset(obj, subset = nCount_Spatial >= 100)

	obj@meta.data[,c("Time","Chip","Heart")] <- str_split_fixed(obj@meta.data$orig.ident,"_",3)[,c(1,2,3)]

	obj@meta.data$Sample <- paste(obj@meta.data$Time,obj@meta.data$Heart,sep="_")

	obj <- SCTransform(obj, assay = 'Spatial', variable.features.n = 2000, return.only.var.genes = FALSE, n_genes=NULL, min_cells=5, method='qpoisson')

	vf <- append(vf, VariableFeatures(obj))

	merge.list[[index]] <- obj
}

data_merge <- merge(x = merge.list[[1]], y = merge.list[2:length(merge.list)])

DefaultAssay(data_merge) <- "SCT"
VariableFeatures(data_merge) <- unique(vf)

##Cluster
data_merge <- RunPCA(data_merge)
data_merge <- RunHarmony(data_merge, group.by.vars = "orig.ident", reduction = "pca")
data_merge <- FindNeighbors(data_merge, dims = 1:30,reduction = "harmony")
data_merge <- FindClusters(data_merge, verbose = FALSE, resolution = c(seq(0.1,2,0.1)))
data_merge <- RunUMAP(data_merge, dims = 1:30, reduction = "harmony")

##Save
saveRDS(data_merge, file = paste0(opts$out, '/', opts$name, "_merge.rds"))
write.table(data_merge@meta.data,file=paste0(opts$out, '/', opts$name, "_merge_meta.txt"), sep ="\t", quote=FALSE, row.names =TRUE, col.names =TRUE)


library(Seurat)
library(anndata)
library(ggplot2)
library(cowplot)
library(patchwork)
library(RColorBrewer)

h5path <-'nice.h5ad'
data = read_h5ad(h5path)
mtx = as.matrix(data$X)
rownames(mtx) = data$obs_names
colnames(mtx) = data$var_names
data.merge <- CreateSeuratObject(counts=t(mtx),meta.data=data$obs,project='zebrafish')

scVi <- read.csv('scVI.csv',header=FALSE)
colnames(scVi) <- paste0("SCVI_", 1:10)
rownames(scVi) = data$obs_names
data.merge[["scvi"]] <- CreateDimReducObject(embeddings = as.matrix(scVi), key = "SCVI_", assay = DefaultAssay(data.merge))
data.merge[["scvi1"]] <- CreateDimReducObject(embeddings = as.matrix(scVi), key = "SCVI_", assay = DefaultAssay(data.merge))

data.merge <- RunUMAP(data.merge,reduction = "scvi1", dims = 1:10)
data.merge <- FindNeighbors(data.merge, reduction = "scvi1", dims = 1:10)
data.merge <- FindClusters(data.merge, resolution = 0.1)
data.merge <- FindClusters(data.merge, resolution = 0.3)
data.merge <- FindClusters(data.merge, resolution = 0.4)
data.merge <- FindClusters(data.merge, resolution = 0.5)
data.merge <- FindClusters(data.merge, resolution = 0.8)
data.merge <- FindClusters(data.merge, resolution = 0.2)
data.merge <- identity(data.merge)

saveRDS(data.merge,'scvi_reduction.integrate.rds')
write.table(data.merge@meta.data,"scvi_reduction.meta.txt",quote=F,sep="\t")

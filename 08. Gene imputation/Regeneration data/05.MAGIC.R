library(Rmagic)
library(Seurat)

obj <- readRDS("../../../../w03_rename/st_merge.rds")
DefaultAssay(obj) <- "Spatial"
obj <- NormalizeData(object = obj)
dat <- t(as.matrix(obj@assays$Spatial@data))

mk <- read.table("../../../01.Spatial/02_output/marker_uniq.lst",header=F,sep="\t")
mk <- mk$V1

dat_MAGIC <- magic(dat, genes=mk, n.jobs = 110)

seurat_obj <- CreateSeuratObject(t(dat_MAGIC$result), assay = "MAGIC", meta.data = obj@meta.data)
saveRDS(seurat_obj, './st_merge.rds')

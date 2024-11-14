##par
parser = argparse::ArgumentParser(description = "RCTD")
parser$add_argument('-i', '--input', help = 'input RDS')
parser$add_argument('-r', '--ref', help = 'reference RDS')
parser$add_argument('-o', '--out', help = 'out directory')
parser$add_argument('-n', '--name', dest = 'name', help = 'sample name')
opts = parser$parse_args()

##lib
library(spacexr)
library(Seurat)
library(tidyr)

##ref
scRNA <- readRDS(opts$ref)
scRNA <- subset(scRNA, subset = celltype != "Others")
ref_counts <- GetAssayData(scRNA, assay = "RNA", slot = "counts")
ref_celltypes <- as.factor(scRNA@meta.data$celltype)
names(ref_celltypes) <- rownames(scRNA@meta.data)
nUMI <- scRNA@meta.data$nCount_RNA
names(nUMI) <- rownames(scRNA@meta.data)
reference <- Reference(ref_counts, ref_celltypes, nUMI)

##query
spatial <- readRDS(opts$input)
spatial@meta.data$x <- spatial@meta.data$grid_x
spatial@meta.data$y <- spatial@meta.data$grid_y
spatial <- subset(spatial, subset = nCount_Spatial >= 100)
spatial_counts <- GetAssayData(spatial, assay = "Spatial", slot = "counts")
spatial_coords <- spatial@meta.data[,c("x","y")]
nUMI <- spatial@meta.data$nCount_Spatial
names(nUMI) <- rownames(spatial@meta.data)
query <- SpatialRNA(spatial_coords, spatial_counts, nUMI)

##rctd
myRCTD <- create.RCTD(query, reference, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
saveRDS(myRCTD, paste0(opts$out, '/', opts$name, "_RCTD.rds"))

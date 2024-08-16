library(spacexr)
library(anndata)
library(Seurat)
library(ggplot2)
library(ggsci)
library(tidyr)

scRNA <- readRDS('ref_nice.RCTD.rds')
ref_counts <- GetAssayData(scRNA, assay = "RNA", slot = "counts")
ref_celltypes <- scRNA@meta.data$cell_type
ref_celltypes <- as.factor(ref_celltypes)
names(ref_celltypes) <- rownames(scRNA@meta.data)
nUMI <- scRNA@meta.data$nCount_RNA
names(nUMI) <- rownames(scRNA@meta.data)
reference <- Reference(ref_counts, ref_celltypes, nUMI)

vlist = c('other.s100.h5ad',
          'other.s101.h5ad',
          'other.s102.h5ad',
          'other.s104.h5ad',
          'other.s105.h5ad',
          'other.s107.h5ad',
          'other.s108.h5ad',
          'other.s109.h5ad',
          'other.s110.h5ad',
          'other.s111.h5ad',
          'other.s113.h5ad',
          'other.s116.h5ad',
          'other.s118.h5ad',
          'other.s119.h5ad',
          'other.s120.h5ad',
          'other.s121.h5ad',
          'other.s122.h5ad',
          'other.s123.h5ad',
          'other.s125.h5ad',
          'other.s126.h5ad',
          'other.s127.h5ad',
          'other.s136.h5ad',
          'other.s140.h5ad',
          'other.s141.h5ad',
          'other.s148.h5ad',
          'other.s149.h5ad',
          'other.s150.h5ad',
          'other.s151.h5ad',
          'other.s155.h5ad',
          'other.s156.h5ad',
          'other.s159.h5ad',
          'other.s163.h5ad',
          'other.s22.h5ad',
          'other.s23.h5ad',
          'other.s24.h5ad',
          'other.s25.h5ad',
          'other.s26.h5ad',
          'other.s28.h5ad',
          'other.s30.h5ad',
          'other.s31.h5ad',
          'other.s32.h5ad',
          'other.s37.h5ad',
          'other.s44.h5ad',
          'other.s50.h5ad',
          'other.s52.h5ad',
          'other.s54.h5ad',
          'other.s56.h5ad',
          'other.s61.h5ad',
          'other.s69.h5ad',
          'other.s70.h5ad',
          'other.s73.h5ad',
          'other.s76.h5ad',
          'other.s77.h5ad',
          'other.s83.h5ad',
          'other.s85.h5ad',
          'other.s86.h5ad',
          'other.s87.h5ad',
          'other.s88.h5ad',
          'other.s89.h5ad',
          'other.s90.h5ad',
          'other.s91.h5ad',
          'other.s92.h5ad',
          'other.s93.h5ad',
          'other.s95.h5ad',
          'other.s96.h5ad',
          'other.s98.h5ad',
          'other.s99.h5ad')


for (index in 1:length(vlist)) {
    h5 <- vlist[index]
    h5path <- paste0("../", h5)
    print(h5path)
    data = read_h5ad(h5path)
    mtx = as.matrix(data$X)
    rownames(mtx) = data$obs_names
    colnames(mtx) = data$var_names
    spatial_counts = t(mtx)
    spatial_coords <- data$obs[,c("new_x","new_y")]
    nUMI <- colSums(spatial_counts)
    puck <- SpatialRNA(spatial_coords, spatial_counts, nUMI)
    barcodes <- colnames(puck@counts)
    print('create.RCTD begin')
    myRCTD <- create.RCTD(puck, reference, max_cores = 2, CELL_MIN_INSTANCE = 5, UMI_min = 10)
    print('create.RCTD done')
    myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
    saveRDS(myRCTD, paste0("./", h5, ".rds"))
    write.csv(myRCTD@results$results_df,paste0(h5,'.result.csv'))
}
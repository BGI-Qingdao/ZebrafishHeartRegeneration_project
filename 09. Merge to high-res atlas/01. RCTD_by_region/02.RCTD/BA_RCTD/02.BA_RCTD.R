library(spacexr)
library(anndata)
library(Seurat)
library(ggplot2)
library(ggsci)
library(tidyr)

scRNA <- readRDS('./scRNA.H_BA.rds')
scRNA <- subset(scRNA,subset =  celltype != "Endocardium (V)")
scRNA <- subset(scRNA,subset =  celltype != "Endocardium (A)")
scRNA <- subset(scRNA,subset =  celltype != "Cardiomyocytes (V)")
scRNA  <- subset(scRNA,subset =  celltype != "Endothelial cells")
scRNA  <- subset(scRNA,subset =  celltype != "Valves")

ref_counts <- GetAssayData(scRNA, assay = "RNA", slot = "counts")
ref_celltypes <- scRNA@meta.data$celltype
ref_celltypes <- as.factor(ref_celltypes)
names(ref_celltypes) <- rownames(scRNA@meta.data)
nUMI <- scRNA@meta.data$nCount_RNA
names(nUMI) <- rownames(scRNA@meta.data)
reference <- Reference(ref_counts, ref_celltypes, nUMI)

for (index in 1:200){
    h5path <- paste0("./00.data/s",index,".h5ad")# , h5)
    if ( file.exists(h5path) ) {
        print(h5path)
        rest = paste0("./s",index,'.result.csv')
        if ( file.exists(rest))
                next
        if (h5path ==  "./00.data/s51.h5ad")
                next
        if (h5path ==  "./00.data/s52.h5ad")
                next
        if (h5path ==  "./00.data/s61.h5ad")
            next
        if (h5path ==  "./00.data/s77.h5ad")
            next
        if (h5path ==  "./00.data/s89.h5ad")
            next
        if (h5path == "./00.data/s112.h5ad")
            next
        if (h5path == "./00.data/s123.h5ad")
            next
        if (h5path == "./00.data/s135.h5ad")
            next
        data = read_h5ad(h5path)
        mtx = as.matrix(data$X)
        rownames(mtx) = data$obs_names
        colnames(mtx) = data$var_names
        spatial_counts = t(mtx)
        spatial_coords <- data$obs[,c("x","y")]
        nUMI <- colSums(spatial_counts)
        puck <- SpatialRNA(spatial_coords, spatial_counts, nUMI)
        barcodes <- colnames(puck@counts)
        print('create.RCTD begin')
        myRCTD <- create.RCTD(puck, reference, max_cores = 2, CELL_MIN_INSTANCE = 5, UMI_min = 10)
        print('create.RCTD done')
        myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
        saveRDS(myRCTD, paste0("./s", index, ".rds"))
        write.csv(myRCTD@results$results_df,paste0("./s",index,'.result.csv'))
    }
}

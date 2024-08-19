library(Seurat)

obj1 <- readRDS('00.data/Zebrafish_scRNA.0dpa.rds')
obj_st <- readRDS('00.data/Zebrafish_stereoseq.0dpa.rds')
obj.list <- SplitObject(obj_st,split.by='orig.ident')
for ( slice in  names(obj.list) ) {
    print(slice)
    seurat_object = obj.list[[slice]]
    i2 <- FindTransferAnchors( reference = obj1,
                               query = seurat_object,
                               features = rownames(seurat_object),
                               reduction = 'cca',
                               reference.assay = 'RNA',
                               query.assay = 'Spatial')

    saveRDS(i2,paste0(slice,".Anchors.rds"))

    refdata <- GetAssayData( object = obj1,
                             assay = 'RNA',
                             slot = 'data')

    imputation <- TransferData( anchorset = i2,
                                query = seurat_object,
                                refdata = refdata,
                                weight.reduction = 'cca')

    saveRDS(imputation, paste0(slice,'.imputation.rds'))
}

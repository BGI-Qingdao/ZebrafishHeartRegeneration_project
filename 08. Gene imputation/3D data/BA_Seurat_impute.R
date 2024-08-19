library(Seurat)
library(future)
library(anndata)

plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100000 * 1024^3)

obj1 <- readRDS('scRNA.H_BA.rds')
data = read_h5ad('BA.3D.h5ad')
mtx = as.matrix(data$X)
rownames(mtx) = data$obs_names
colnames(mtx) = data$var_names
spa <- CreateSeuratObject(counts=t(mtx),meta.data=data$obs)
saveRDS(spa,'tmp.rds')

obj.list <- SplitObject(spa,split.by='sid')
for ( slice in  names(obj.list) ) {
    print(slice)
    seurat_object = obj.list[[slice]]
    if (dim(seurat_object@meta.data)[1]<=100){
        print('to small')
    }else{
        print('FindTransferAnchors')
        i2 <- FindTransferAnchors(  reference = obj1,
                                query = seurat_object,
                                features = rownames(seurat_object),
                                reduction = 'cca',
                                reference.assay = 'RNA',
                                query.assay = 'RNA')


        print('GetAssayData')
        refdata <- GetAssayData( object = obj1,
                                 assay = 'RNA',
                                 slot = 'data')

        print('TransferData')
        imputation <- TransferData( anchorset = i2,
                                    query =seurat_object ,
                                    refdata = refdata,
                                    weight.reduction = 'cca')
        saveRDS(imputation,paste0(slice,'imputation.rds'))
    }
}

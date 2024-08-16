library(Seurat)
library(SummarizedExperiment)
library(scuttle)

data<- readRDS("scvi_reduction.integrate.anno3.rds")
data <-subset(data, subset = cell_type != 'NA')
saveRDS(data,"ref_nice.RCTD.rds")

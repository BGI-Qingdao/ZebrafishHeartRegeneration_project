##Rscript for merge data zebrafish heart reg
setwd('/dellfsqd2/C_OCEAN/USERS/c-lumeina1/data/scrna/')
library(Seurat)
library(dplyr)
####H_0dpa data all 6 samples
dir_H_0dpa_1<-('Matrix/H-0d/0197-1-220119/')
count_H_0dpa_1<-Read10X(data.dir=dir_H_0dpa_1,gene.column=1)
H_0dpa_1<-CreateSeuratObject(counts=count_H_0dpa_1,project='H_0dpa',sample=1,min.cells=3,min.features=200)
H_0dpa_1@meta.data$sample<-c('0197-1-220119')

dir_H_0dpa_2<-('Matrix/H-0d/0197-2-220119/')
count_H_0dpa_2<-Read10X(data.dir=dir_H_0dpa_2,gene.column=1)
H_0dpa_2<-CreateSeuratObject(counts=count_H_0dpa_2,project='H_0dpa',min.cells=3,min.features=200)
H_0dpa_2@meta.data$sample<-c('0197-2-220119')

dir_H_0dpa_3<-('Matrix/H-0d/0197-3-220119/')
count_H_0dpa_3<-Read10X(data.dir=dir_H_0dpa_3,gene.column=1)
H_0dpa_3<-CreateSeuratObject(counts=count_H_0dpa_3,project='H_0dpa',min.cells=3,min.features=200)
H_0dpa_3@meta.data$sample<-c('0197-3-220119')

dir_H_0dpa_4<-('Matrix/H_0d/7057-1-220908/')
count_H_0dpa_4<-Read10X(data.dir=dir_H_0dpa_4,gene.column=1)
H_0dpa_4<-CreateSeuratObject(counts=count_H_0dpa_4,project='H_0dpa',min.cells=3,min.features=200)
H_0dpa_4@meta.data$sample<-c('7057-1-220908')

dir_H_0dpa_5<-('Matrix/H_0d/7057-2-220908/')
count_H_0dpa_5<-Read10X(data.dir=dir_H_0dpa_5,gene.column=1)
H_0dpa_5<-CreateSeuratObject(counts=count_H_0dpa_5,project='H_0dpa',min.cells=3,min.features=200)
H_0dpa_5@meta.data$sample<-c('7057-2-220908')

dir_H_0dpa_6<-('Matrix/H_0d/7057-3-220908/')
count_H_0dpa_6<-Read10X(data.dir=dir_H_0dpa_6,gene.column=1)
H_0dpa_6<-CreateSeuratObject(counts=count_H_0dpa_6,project='H_0dpa',min.cells=3,min.features=200)
H_0dpa_6@meta.data$sample<-c('7057-3-220908')

H_0dpa<-merge(H_0dpa_1,y=c(H_0dpa_2,H_0dpa_3,H_0dpa_4,H_0dpa_5,H_0dpa_6),add.cell.ids=c('0dpa_1','0dpa_2','0dpa_3','0dpa_4','0dpa_5','0dpa_6'),project='H_0dpa')

H_0dpa[["percent.mt"]] <- PercentageFeatureSet(H_0dpa, pattern = "^MT-|^mt-")
pdf("H_0dpa_QC.pdf",width=15,height=8)
VlnPlot(object = H_0dpa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(H_0dpa, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(H_0dpa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
H_0dpa<- subset(H_0dpa, subset = nFeature_RNA>200 & nFeature_RNA < 4000 & percent.mt < 40)
saveRDS(H_0dpa, file = "./H_0dpa.Rds")
rm(list = ls())

####H_atrium all 3 samples
dir_H_atrium_1<-('Matrix/H_atrium/4100-1-220526/')
count_H_atrium_1<-Read10X(data.dir=dir_H_atrium_1,gene.column=1)
H_atrium_1<-CreateSeuratObject(counts=count_H_atrium_1,project='H_atrium',min.cells=3,min.features=200)
H_atrium_1@meta.data$sample<-c('4100-1-220526')

dir_H_atrium_2<-('Matrix/H_atrium/4100-2-220526/')
count_H_atrium_2<-Read10X(data.dir=dir_H_atrium_2,gene.column=1)
H_atrium_2<-CreateSeuratObject(counts=count_H_atrium_2,project='H_atrium',min.cells=3,min.features=200)
H_atrium_2@meta.data$sample<-c('4100-2-220526')

dir_H_atrium_3<-('Matrix/H_atrium/4100-3-220526/')
count_H_atrium_3<-Read10X(data.dir=dir_H_atrium_3,gene.column=1)
H_atrium_3<-CreateSeuratObject(counts=count_H_atrium_3,project='H_atrium',min.cells=3,min.features=200)
H_atrium_3@meta.data$sample<-c('4100-3-220526')

H_atrium<-merge(H_atrium_1,y=c(H_atrium_2,H_atrium_3),add.cell.ids=c('atrium_1','atrium_2','atrium_3'),project='H_atrium')

H_atrium[["percent.mt"]] <- PercentageFeatureSet(H_atrium, pattern = "^MT-|^mt-")
pdf("H_atrium_QC.pdf",width=15,height=8)
VlnPlot(object = H_atrium, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(H_atrium, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(H_atrium, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
H_atrium<- subset(H_atrium, subset = nFeature_RNA>200 & nFeature_RNA < 4000 & percent.mt < 20)
saveRDS(H_atrium, file = "./H_atrium.Rds")
rm(list = ls())

####H_ventricle all 3 samples
dir_H_ventricle_1<-('Matrix/H_ventricle/4101-1-220526/')
count_H_ventricle_1<-Read10X(data.dir=dir_H_ventricle_1,gene.column=1)
H_ventricle_1<-CreateSeuratObject(counts=count_H_ventricle_1,project='H_ventricle',min.cells=3,min.features=200)
H_ventricle_1@meta.data$sample<-c('4101-1-220526')

dir_H_ventricle_2<-('Matrix/H_ventricle/4101-2-220526/')
count_H_ventricle_2<-Read10X(data.dir=dir_H_ventricle_2,gene.column=1)
H_ventricle_2<-CreateSeuratObject(counts=count_H_ventricle_2,project='H_ventricle',min.cells=3,min.features=200)
H_ventricle_2@meta.data$sample<-c('4101-2-220526')

dir_H_ventricle_3<-('Matrix/H_ventricle/4101-3-220526/')
count_H_ventricle_3<-Read10X(data.dir=dir_H_ventricle_3,gene.column=1)
H_ventricle_3<-CreateSeuratObject(counts=count_H_ventricle_3,project='H_ventricle',min.cells=3,min.features=200)
H_ventricle_3@meta.data$sample<-c('4101-3-220526')

H_ventricle<-merge(H_ventricle_1,y=c(H_ventricle_2,H_ventricle_3),add.cell.ids=c('ventricle_1','ventricle_2','ventricle_3'),project='H_ventricle')

H_ventricle[["percent.mt"]] <- PercentageFeatureSet(H_ventricle, pattern = "^MT-|^mt-")
pdf("H_ventricle_QC.pdf",width=15,height=8)
VlnPlot(object = H_ventricle, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(H_ventricle, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(H_ventricle, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
H_ventricle<- subset(H_ventricle, subset = nFeature_RNA>200 & nFeature_RNA < 4000 & percent.mt < 40)
saveRDS(H_ventricle, file = "./H_ventricle.Rds")
rm(list = ls())

####H_BA all 3 samples
dir_H_BA_1<-('Matrix/H_aortic_bulb/4099-1-220526/')
count_H_BA_1<-Read10X(data.dir=dir_H_BA_1,gene.column=1)
H_BA_1<-CreateSeuratObject(counts=count_H_BA_1,project='H_BA',min.cells=3,min.features=200)
H_BA_1@meta.data$sample<-c('4099-1-220526')

dir_H_BA_2<-('Matrix/H_aortic_bulb/4099-2-220526/')
count_H_BA_2<-Read10X(data.dir=dir_H_BA_2,gene.column=1)
H_BA_2<-CreateSeuratObject(counts=count_H_BA_2,project='H_BA',min.cells=3,min.features=200)
H_BA_2@meta.data$sample<-c('4099-2-220526')

dir_H_BA_3<-('Matrix/H_aortic_bulb/4099-3-220526/')
count_H_BA_3<-Read10X(data.dir=dir_H_BA_3,gene.column=1)
H_BA_3<-CreateSeuratObject(counts=count_H_BA_3,project='H_BA',min.cells=3,min.features=200)
H_BA_3@meta.data$sample<-c('4099-3-220526')

H_BA<-merge(H_BA_1,y=c(H_BA_2,H_BA_3),add.cell.ids=c('BA_1','BA_2','BA_3'),project='H_BA')

H_BA[["percent.mt"]] <- PercentageFeatureSet(H_BA, pattern = "^MT-|^mt-")
pdf("H_BA_QC.pdf",width=15,height=8)
VlnPlot(object = H_BA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(H_BA, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(H_BA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
H_BA<- subset(H_BA, subset = nFeature_RNA>200 & nFeature_RNA < 4000 & percent.mt < 5)
saveRDS(H_BA, file = "./H_BA.Rds")
rm(list = ls())

####H_6hpa all 3 samples
dir_H_6hpa_1<-('Matrix/H_6h/7058-1-220908/')
count_H_6hpa_1<-Read10X(data.dir=dir_H_6hpa_1,gene.column=1)
H_6hpa_1<-CreateSeuratObject(counts=count_H_6hpa_1,project='H_6hpa',min.cells=3,min.features=200)
H_6hpa_1@meta.data$sample<-c('7058-1-220908')

dir_H_6hpa_2<-('Matrix/H_6h/7058-2-220908/')
count_H_6hpa_2<-Read10X(data.dir=dir_H_6hpa_2,gene.column=1)
H_6hpa_2<-CreateSeuratObject(counts=count_H_6hpa_2,project='H_6hpa',min.cells=3,min.features=200)
H_6hpa_2@meta.data$sample<-c('7058-2-220908')

dir_H_6hpa_3<-('Matrix/H_6h/7058-3-220908/')
count_H_6hpa_3<-Read10X(data.dir=dir_H_6hpa_3,gene.column=1)
H_6hpa_3<-CreateSeuratObject(counts=count_H_6hpa_3,project='H_6hpa',min.cells=3,min.features=200)
H_6hpa_3@meta.data$sample<-c('7058-3-220908')

H_6hpa<-merge(H_6hpa_1,y=c(H_6hpa_2,H_6hpa_3),add.cell.ids=c('6hpa_1','6hpa_2','6hpa_3'),project='H_6hpa')

H_6hpa[["percent.mt"]] <- PercentageFeatureSet(H_6hpa, pattern = "^MT-|^mt-")
pdf("H_6hpa_QC.pdf",width=15,height=8)
VlnPlot(object = H_6hpa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(H_6hpa, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(H_6hpa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
H_6hpa<- subset(H_6hpa, subset = nFeature_RNA>200 & nFeature_RNA < 4000 & percent.mt < 20)
saveRDS(H_6hpa, file = "./H_6hpa.Rds")
rm(list = ls())

####H_12hpa all 3 samples
dir_H_12hpa_1<-('Matrix/H_12h/4165-1-220620/')
count_H_12hpa_1<-Read10X(data.dir=dir_H_12hpa_1,gene.column=1)
H_12hpa_1<-CreateSeuratObject(counts=count_H_12hpa_1,project='H_12hpa',min.cells=3,min.features=200)
H_12hpa_1@meta.data$sample<-c('4165-1-220620')

dir_H_12hpa_2<-('Matrix/H_12h/4165-2-220620/')
count_H_12hpa_2<-Read10X(data.dir=dir_H_12hpa_2,gene.column=1)
H_12hpa_2<-CreateSeuratObject(counts=count_H_12hpa_2,project='H_12hpa',min.cells=3,min.features=200)
H_12hpa_2@meta.data$sample<-c('4165-2-220620')

dir_H_12hpa_3<-('Matrix/H_12h/4165-3-220620/')
count_H_12hpa_3<-Read10X(data.dir=dir_H_12hpa_3,gene.column=1)
H_12hpa_3<-CreateSeuratObject(counts=count_H_12hpa_3,project='H_12hpa',min.cells=3,min.features=200)
H_12hpa_3@meta.data$sample<-c('4165-3-220620')

H_12hpa<-merge(H_12hpa_1,y=c(H_12hpa_2,H_12hpa_3),add.cell.ids=c('12hpa_1','12hpa_2','12hpa_3'),project='H_12hpa')

H_12hpa[["percent.mt"]] <- PercentageFeatureSet(H_12hpa, pattern = "^MT-|^mt-")
pdf("H_12hpa_QC.pdf",width=15,height=8)
VlnPlot(object = H_12hpa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(H_12hpa, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(H_12hpa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
H_12hpa<- subset(H_12hpa, subset = nFeature_RNA>200 & nFeature_RNA < 5000 & percent.mt < 20)
saveRDS(H_12hpa, file = "./H_12hpa.Rds")
rm(list = ls())


####H_1dpa all 3 samples
dir_H_1dpa_1<-('Matrix/H_1d/4166-1-220620/')
count_H_1dpa_1<-Read10X(data.dir=dir_H_1dpa_1,gene.column=1)
H_1dpa_1<-CreateSeuratObject(counts=count_H_1dpa_1,project='H_1dpa',min.cells=3,min.features=200)
H_1dpa_1@meta.data$sample<-c('4166-1-220620')

dir_H_1dpa_2<-('Matrix/H_1d/4166-2-220620/')
count_H_1dpa_2<-Read10X(data.dir=dir_H_1dpa_2,gene.column=1)
H_1dpa_2<-CreateSeuratObject(counts=count_H_1dpa_2,project='H_1dpa',min.cells=3,min.features=200)
H_1dpa_2@meta.data$sample<-c('4166-2-220620')

dir_H_1dpa_3<-('Matrix/H_1d/4166-3-220620/')
count_H_1dpa_3<-Read10X(data.dir=dir_H_1dpa_3,gene.column=1)
H_1dpa_3<-CreateSeuratObject(counts=count_H_1dpa_3,project='H_1dpa',min.cells=3,min.features=200)
H_1dpa_3@meta.data$sample<-c('4166-3-220620')

H_1dpa<-merge(H_1dpa_1,y=c(H_1dpa_2,H_1dpa_3),add.cell.ids=c('1dpa_1','1dpa_2','1dpa_3'),project='H_1dpa')

H_1dpa[["percent.mt"]] <- PercentageFeatureSet(H_1dpa, pattern = "^MT-|^mt-")
pdf("H_1dpa_QC.pdf",width=15,height=8)
VlnPlot(object = H_1dpa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(H_1dpa, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(H_1dpa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
H_1dpa<- subset(H_1dpa, subset = nFeature_RNA>200 & nFeature_RNA < 5000 & percent.mt < 20)
saveRDS(H_1dpa, file = "./H_1dpa.Rds")
rm(list = ls())

####H_3dpa all 3 samples
dir_H_3dpa_1<-('Matrix/H_3d/4167-1-220620/')
count_H_3dpa_1<-Read10X(data.dir=dir_H_3dpa_1,gene.column=1)
H_3dpa_1<-CreateSeuratObject(counts=count_H_3dpa_1,project='H_3dpa',min.cells=3,min.features=200)
H_3dpa_1@meta.data$sample<-c('4167-1-220620')

dir_H_3dpa_2<-('Matrix/H_3d/4167-2-220620/')
count_H_3dpa_2<-Read10X(data.dir=dir_H_3dpa_2,gene.column=1)
H_3dpa_2<-CreateSeuratObject(counts=count_H_3dpa_2,project='H_3dpa',min.cells=3,min.features=200)
H_3dpa_2@meta.data$sample<-c('4167-2-220620')

dir_H_3dpa_3<-('Matrix/H_3d/4167-3-220620/')
count_H_3dpa_3<-Read10X(data.dir=dir_H_3dpa_3,gene.column=1)
H_3dpa_3<-CreateSeuratObject(counts=count_H_3dpa_3,project='H_3dpa',min.cells=3,min.features=200)
H_3dpa_3@meta.data$sample<-c('4167-3-220620')

H_3dpa<-merge(H_3dpa_1,y=c(H_3dpa_2,H_3dpa_3),add.cell.ids=c('3dpa_1','3dpa_2','3dpa_3'),project='H_3dpa')

H_3dpa[["percent.mt"]] <- PercentageFeatureSet(H_3dpa, pattern = "^MT-|^mt-")
pdf("H_3dpa_QC.pdf",width=15,height=8)
VlnPlot(object = H_3dpa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(H_3dpa, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(H_3dpa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
H_3dpa<- subset(H_3dpa, subset = nFeature_RNA>200 & nFeature_RNA < 4000 & percent.mt < 20)
saveRDS(H_3dpa, file = "./H_3dpa.Rds")
rm(list = ls())

####H_7dpa all 3 samples
dir_H_7dpa_1<-('Matrix/xinzang_H_7d/5565-1-220728/')
count_H_7dpa_1<-Read10X(data.dir=dir_H_7dpa_1,gene.column=1)
H_7dpa_1<-CreateSeuratObject(counts=count_H_7dpa_1,project='H_7dpa',min.cells=3,min.features=200)
H_7dpa_1@meta.data$sample<-c('5565-1-220728')

dir_H_7dpa_2<-('Matrix/xinzang_H_7d/5565-2-220728/')
count_H_7dpa_2<-Read10X(data.dir=dir_H_7dpa_2,gene.column=1)
H_7dpa_2<-CreateSeuratObject(counts=count_H_7dpa_2,project='H_7dpa',min.cells=3,min.features=200)
H_7dpa_2@meta.data$sample<-c('5565-2-220728')

dir_H_7dpa_3<-('Matrix/xinzang_H_7d/5565-3-220728/')
count_H_7dpa_3<-Read10X(data.dir=dir_H_7dpa_3,gene.column=1)
H_7dpa_3<-CreateSeuratObject(counts=count_H_7dpa_3,project='H_7dpa',min.cells=3,min.features=200)
H_7dpa_3@meta.data$sample<-c('5565-3-220728')

H_7dpa<-merge(H_7dpa_1,y=c(H_7dpa_2,H_7dpa_3),add.cell.ids=c('7dpa_1','7dpa_2','7dpa_3'),project='H_7dpa')

H_7dpa[["percent.mt"]] <- PercentageFeatureSet(H_7dpa, pattern = "^MT-|^mt-")
pdf("H_7dpa_QC.pdf",width=15,height=8)
VlnPlot(object = H_7dpa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(H_7dpa, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(H_7dpa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
H_7dpa<- subset(H_7dpa, subset = nFeature_RNA>200 & nFeature_RNA < 4000 & percent.mt < 20)
saveRDS(H_7dpa, file = "./H_7dpa.Rds")
rm(list = ls())


####H_14dpa all 3 samples
dir_H_14dpa_1<-('Matrix/H_14d/7059-1-220908/')
count_H_14dpa_1<-Read10X(data.dir=dir_H_14dpa_1,gene.column=1)
H_14dpa_1<-CreateSeuratObject(counts=count_H_14dpa_1,project='H_14dpa',min.cells=3,min.features=200)
H_14dpa_1@meta.data$sample<-c('7059-1-220908')

dir_H_14dpa_2<-('Matrix/H_14d/7059-2-220908/')
count_H_14dpa_2<-Read10X(data.dir=dir_H_14dpa_2,gene.column=1)
H_14dpa_2<-CreateSeuratObject(counts=count_H_14dpa_2,project='H_14dpa',min.cells=3,min.features=200)
H_14dpa_2@meta.data$sample<-c('7059-2-220908')

dir_H_14dpa_3<-('Matrix/H_14d/7059-3-220908/')
count_H_14dpa_3<-Read10X(data.dir=dir_H_14dpa_3,gene.column=1)
H_14dpa_3<-CreateSeuratObject(counts=count_H_14dpa_3,project='H_14dpa',min.cells=3,min.features=200)
H_14dpa_3@meta.data$sample<-c('7059-3-220908')

H_14dpa<-merge(H_14dpa_1,y=c(H_14dpa_2,H_14dpa_3),add.cell.ids=c('14dpa_1','14dpa_2','14dpa_3'),project='H_14dpa')

H_14dpa[["percent.mt"]] <- PercentageFeatureSet(H_14dpa, pattern = "^MT-|^mt-")
pdf("H_14dpa_QC.pdf",width=15,height=8)
VlnPlot(object = H_14dpa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(H_14dpa, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(H_14dpa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
H_14dpa<- subset(H_14dpa, subset = nFeature_RNA>200 & nFeature_RNA < 4000 & percent.mt < 30)
saveRDS(H_14dpa, file = "./H_14dpa.Rds")
rm(list = ls())


####H_28dpa all 3 samples
dir_H_28dpa_1<-('Matrix/xinzang_H_28d/5566-1-220728/')
count_H_28dpa_1<-Read10X(data.dir=dir_H_28dpa_1,gene.column=1)
H_28dpa_1<-CreateSeuratObject(counts=count_H_28dpa_1,project='H_28dpa',min.cells=3,min.features=200)
H_28dpa_1@meta.data$sample<-c('5566-1-220728')

dir_H_28dpa_2<-('Matrix/xinzang_H_28d/5566-2-220728/')
count_H_28dpa_2<-Read10X(data.dir=dir_H_28dpa_2,gene.column=1)
H_28dpa_2<-CreateSeuratObject(counts=count_H_28dpa_2,project='H_28dpa',min.cells=3,min.features=200)
H_28dpa_2@meta.data$sample<-c('5566-2-220728')

dir_H_28dpa_3<-('Matrix/xinzang_H_28d/5566-3-220728/')
count_H_28dpa_3<-Read10X(data.dir=dir_H_28dpa_3,gene.column=1)
H_28dpa_3<-CreateSeuratObject(counts=count_H_28dpa_3,project='H_28dpa',min.cells=3,min.features=200)
H_28dpa_3@meta.data$sample<-c('5566-3-220728')

H_28dpa<-merge(H_28dpa_1,y=c(H_28dpa_2,H_28dpa_3),add.cell.ids=c('28dpa_1','28dpa_2','28dpa_3'),project='H_28dpa')

H_28dpa[["percent.mt"]] <- PercentageFeatureSet(H_28dpa, pattern = "^MT-|^mt-")
pdf("H_28dpa_QC.pdf",width=15,height=8)
VlnPlot(object = H_28dpa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(H_28dpa, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(H_28dpa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
H_28dpa<- subset(H_28dpa, subset = nFeature_RNA>200 & nFeature_RNA < 5000 & percent.mt < 25)
saveRDS(H_28dpa, file = "./H_28dpa.Rds")
rm(list = ls())

H_0dpa <- readRDS("H_0dpa.Rds")
H_1dpa <- readRDS("H_1dpa.Rds")
H_3dpa <- readRDS("H_3dpa.Rds")
H_BA <- readRDS("H_BA.Rds")
H_ventricle <- readRDS("H_ventricle.Rds")
H_atrium <- readRDS("H_atrium.Rds")
H_28dpa <- readRDS("H_28dpa.Rds")
H_14dpa <- readRDS("H_14dpa.Rds")
H_12hpa <- readRDS("H_12hpa.Rds")
H_7dpa <- readRDS("H_7dpa.Rds")
H_6hpa <- readRDS("H_6hpa.Rds")

Heart_merge<-merge(H_0dpa,y=c(H_1dpa,H_3dpa,H_BA,H_ventricle,H_atrium,H_28dpa,H_14dpa,H_12hpa,H_7dpa,H_6hpa),add.cell.ids=c('H_0dpa','H_1dpa','H_3dpa','H_BA','H_ventricle','H_atrium','H_28dpa','H_14dpa','H_12hpa','H_7dpa','H_6hpa'),project='Heart_zebrafish')
saveRDS(Heart_merge, file = "./Heart_merge.Rds")
pdf("Heart_merge_afterQC.pdf",width=15,height=8)
VlnPlot(object = Heart_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(Heart_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Heart_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

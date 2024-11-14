library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
#library(tidyverse)
library(patchwork)
library(harmony)
library(clustree)
library(ComplexHeatmap)
library(VennDiagram)
library(UpSetR)
###merge初步质控后的数据
H_0dpa <- readRDS("../H_0dpa/H_0dpa.Rds")
H_6hpa <- readRDS("../H_6hpa/H_6hpa.Rds")
H_12hpa <- readRDS("../H_12hpa/H_12hpa.Rds")
H_1dpa <- readRDS("../H_1dpa/H_1dpa.Rds")
H_3dpa <- readRDS("../H_3dpa/H_3dpa.Rds")
H_7dpa <- readRDS("../H_7dpa/H_7dpa.Rds")
H_14dpa <- readRDS("../H_14dpa/H_14dpa.Rds")
H_28dpa <- readRDS("../H_28dpa/H_28dpa.Rds")
H_atrium <- readRDS("../H_atrium/H_atrium.Rds")
H_ventricle <- readRDS("../H_ventricle/H_ventricle.Rds")
H_BA <- readRDS("../H_BA/H_BA.Rds")

list<-list(H_0dpa=H_0dpa,H_6hpa=H_6hpa,H_12hpa=H_12hpa,H_1dpa=H_1dpa,H_3dpa=H_3dpa,H_7dpa=H_7dpa,H_14dpa=H_14dpa,H_28dpa=H_28dpa,H_atrium=H_atrium,H_ventricle=H_ventricle,H_BA=H_BA)


##还要进一步质控红细胞和核糖体基因
for (i in 1:length(list)) {
    if(T){
        list[[i]][['percent.rp']]<-PercentageFeatureSet(list[[i]], pattern = "^rp[sl]")
    }
    if(T){
        HB.genes <- c("hbaa1",'hbae1.1','hbbe1.1',"hbbe3","hbae1.3","hbae3","hbba1","hbae5","hbba2","hbbe1.2","hbbe2","hbaa2")
        HB.genes <- CaseMatch(HB.genes, rownames(list[[i]])) 
        list[[i]][['percent.HB']]<-PercentageFeatureSet(list[[i]], features=HB.genes)
    }
}


###merge数据
Heart_all <- merge(list[[1]], list[2:length(list)])
saveRDS(Heart_all, file = "beforeqc_Heart_all.Rds")
#查看每个样本的细胞数
table(Heart_all@meta.data$orig.ident) 
Heart_all@meta.data$orig.ident <- factor( Heart_all@meta.data$orig.ident,levels=c("H_atrium","H_ventricle","H_BA","H_0dpa","H_6hpa","H_12hpa","H_1dpa","H_3dpa","H_7dpa","H_14dpa","H_28dpa"))

###画小提琴图
theme.set2=theme(axis.title.x=element_blank())
###设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.HB", "percent.rp")
group = "orig.ident"
###质控前小提琴图
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(Heart_all, group.by=group, pt.size = 0,features = plot.featrures[i])+ theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2) 
ggsave("before_QC.pdf", plot = violin,width = 20, height = 8)

qclist<-list(H_0dpa=H_0dpa,H_BA=H_BA,H_ventricle=H_ventricle,H_atrium=H_atrium)
###再次质控
Heart_all <- subset(Heart_all,subset = percent.HB < 100)
#查看每个样本的细胞数
table(Heart_all@meta.data$orig.ident) 

plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(Heart_all, group.by=group, pt.size = 0,features = plot.featrures[i])+ theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2) 
ggsave("after_QC.pdf", plot = violin,width = 20, height = 8)
saveRDS(Heart_all, file = "afterqc_Heart_all.Rds")

Heart_all<-readRDS('afterqc_Heart_all.Rds')
afterqc.list <- SplitObject(Heart_all, split.by = "orig.ident")
# 分别对每个样本数据进行标准化并筛选高变基因
for (i in 1:length(afterqc.list)) {
    afterqc.list[[i]] <- NormalizeData(afterqc.list[[i]], verbose = FALSE)
    afterqc.list[[i]] <- FindVariableFeatures(afterqc.list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}
# 提取每个数据集中筛选到的高变基因
hvgs_per_dataset <- lapply(afterqc.list, function(x) {
    x@assays$RNA@var.features
})
# 绘制venn图查看不同样本中高变基因的重叠情况
pdf('venn_hvgs.pdf',wi=15,he=10)
#venn::venn(hvgs_per_dataset, opacity = 0.5, zcolor = (scales::hue_pal())(9), cexsn = 1, 
#    cexil = 1, lwd = 1, col = "white", frame = F, borders = NA)
upset(fromList(hvgs_per_dataset),nsets = 11,main.bar.color = "black", sets = c("H_atrium","H_ventricle","H_BA","H_0dpa","H_6hpa","H_12hpa","H_1dpa","H_3dpa","H_7dpa","H_14dpa","H_28dpa"),point.size = 3)   
dev.off()
write.csv(hvgs_per_dataset,'hvgs_per_dataset.csv')

###merge后的数据降维聚类
Heart_all <- NormalizeData(Heart_all) %>%FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
###查看高可变基因
top20 <- head(VariableFeatures(Heart_all), 20)
plot1 <- VariableFeaturePlot(Heart_all)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
pdf("Top20_VariableFeaturePlot.pdf",width=10,height=8)
plot2
dev.off()

Heart_all<- RunPCA(Heart_all, verbose = F)
pdf('ElbowPlot.pdf')
ElbowPlot(Heart_all, ndims = 50)
dev.off()
Heart_all <- Heart_all %>% RunTSNE(dims=1:30) %>% RunUMAP(dims=1:30)
Heart_all <- FindNeighbors(Heart_all, dims=1:30) %>% FindClusters()

pdf('merge_pca_umap_tsne.pdf',wi=20,he=6)
plot_grid(ncol = 3,
  DimPlot(Heart_all, reduction = "pca", group.by = "orig.ident",pt.size = 1)+NoAxes()+ggtitle("PCA raw_data"),
  DimPlot(Heart_all, reduction = "tsne", group.by = "orig.ident",pt.size = 1)+NoAxes()+ggtitle("tSNE raw_data"),
  DimPlot(Heart_all, reduction = "umap", group.by = "orig.ident",pt.size = 1)+NoAxes()+ggtitle("UMAP raw_data"))
dev.off()
pdf('merge_everysample.pdf',wi=16,he=5)
DimPlot(Heart_all, group.by = "orig.ident",split.by = "orig.ident",pt.size = 1,order = T)
dev.off()

########################################################################################################################################################开始去批次,这时候要没有降维聚类的数据
###harmony整合
#library(harmony)
###整合方法1：单个样本间进行整合（推荐，效果更好）
# group.by.vars参数是设置按哪个分组来整合
#max.iter.harmony设置迭代次数，默认是10,运行RunHarmony结果会提示在迭代多少次后完成了收敛。
#Heart_all.harmony <- RunHarmony(Heart_all, group.by.vars = "orig.ident", reduction = "pca", 
#                                                     dims.use = 1:50, assay.use = "RNA")
#saveRDS(Heart_all.harmony, file = "harmony_Heart_all.Rds")


#pdf('harmony_ElbowPlot.pdf')
#ElbowPlot(Heart_all.harmony, ndims = 50)
#dev.off()

#Heart_all.harmony <- Heart_all.harmony %>% RunTSNE(dims=1:30) %>% RunUMAP(dims=1:30)

#pdf('compare_harmony_merge.pdf',wi=25,he=15)
#plot_grid(ncol = 3,
#  DimPlot(Heart_all, reduction = "pca", group.by = "orig.ident",pt.size = 1)+NoAxes()+ggtitle("PCA raw_data"),
#  DimPlot(Heart_all, reduction = "tsne", group.by = "orig.ident",pt.size = 1)+NoAxes()+ggtitle("tSNE raw_data"),
#  DimPlot(Heart_all, reduction = "umap", group.by = "orig.ident",pt.size = 1)+NoAxes()+ggtitle("UMAP raw_data"),#

#  DimPlot(Heart_all.harmony, reduction = "pca", group.by = "orig.ident",pt.size = 1)+NoAxes()+ggtitle("PCA harmony"),
#  DimPlot(Heart_all.harmony, reduction = "tsne", group.by = "orig.ident",pt.size = 1)+NoAxes()+ggtitle("tSNE harmony"),
#  DimPlot(Heart_all.harmony, reduction = "umap", group.by = "orig.ident",pt.size = 1)+NoAxes()+ggtitle("UMAP harmony"))
#dev.off()
#pdf('compare_harmony_merge_everysample.pdf',wi=16,he=5)
#DimPlot(Heart_all, group.by = "orig.ident",split.by = "orig.ident",pt.size = 1)
#DimPlot(Heart_all.harmony, group.by = "orig.ident",split.by = "orig.ident",pt.size = 1)
#dev.off()
#saveRDS(Heart_all.harmony, file = "harmony_Heart_all.Rds")
##################################################################################################################################################################################################

###harmony和单纯merge的结果没大有什么区别
###正式降维聚类
tree <- FindClusters(data,resolution = c(seq(0,2,0.1)))
clustree.out <- clustree(tree)
pdf('clustree.pdf',width=20,height=20)
print(clustree.out)
dev.off()

Heart_all<- FindClusters(Heart_all, resolution = 0.3)
write.table(table(Heart_all@active.ident),file="Heart_all_cellnumbers.xls")

pdf('umap_tsne_cluster.pdf',wi=25,he=20)
plot_grid(ncol = 1,
  DimPlot(Heart_all, reduction = "tsne", group.by = c('ident',"orig.ident"),pt.size = 0.5,label=T)+NoAxes()+ggtitle("tSNE Heart_all"),
  DimPlot(Heart_all, reduction = "umap", group.by = c('ident',"orig.ident"),pt.size = 0.5,label=T)+NoAxes()+ggtitle("UMAP Heart_all"))
dev.off()
saveRDS(Heart_all, file = "aftercluster.Rds")

test.use = "wilcox"
Heart_all.marker <- FindAllMarkers(Heart_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = test.use)    
write.csv(Heart_all.marker ,file='findmarkers.csv')
top50marker <- Heart_all.marker  %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50marker,file="top50marker.csv")
top20marker <- Heart_all.marker  %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20marker,file="top20marker.csv")
top10marker <- Heart_all.marker  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10marker,file="top10marker.csv")
top5marker <- Heart_all.marker  %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(top5marker,file="top5marker.csv")


subobj <- subset(Heart_all, downsample = 300)
pdf("top10marerheatmap.pdf",width=15,height=15)
DoHeatmap(subobj, features = top10marker$gene)+scale_fill_gradientn(colors= c("navy","white","firebrick3"))
dev.off()





mycolors <- c('#A6CEE3FF','#1F78B4FF', '#B2DF8AFF', '#33A02CFF' ,'#FB9A99FF' ,'#E31A1CFF' ,'#FDBF6FFF' ,'#FF7F00FF','#CAB2D6FF' ,'#6A3D9AFF' ,'#FFFF99FF' ,'#B15928FF','#803620FF','#E5D2DD',  '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658','#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')

pdf("top5_dotplot.pdf",width=30,height = 12)
DotPlot(Heart_all,features = unique(top5marker$gene),cols=mycolors) + RotatedAxis()
dev.off()


new.cluster.ids <- c("Endocardium(Ventricle)","Red blood cells","Smooth muscle cells","Endocardium(Atrium)","Macrophages1","Fibroblasts1","Cardiomyocytes(Ventricle)","T-cells", "Endocardium(frzb)", "Neuron cells", "Macrophages2", "Valve_Fibroblasts", "Epicardium", "Endothelial cells(plvapb)", "Coronary Ecs", "Cardiomyocytes(Atrium)", "Neutrophils", "Macrophages3", "Iymphatic ECs", "Endocardium", "Myelin cells",'Red blood cells2','Fibroblasts2','NA')

names(new.cluster.ids) <- levels(Heart_all) 
Heart_all <- RenameIdents(Heart_all, new.cluster.ids)      
Heart_all$celltype <- Idents(Heart_all)

saveRDS(Heart_all, file = "celltype.Rds")


pdf('marker_ident_umap.pdf',wi=15,he=10)
DimPlot(Heart_all,cols= mycolors , reduction = "umap", group.by = c('ident'),pt.size = 0.8,label=T)+ggtitle("Heart_all")
dev.off()
pdf('marker_ident_umap_sample.pdf',wi=160,he=10)
DimPlot(Heart_all,cols= mycolors , reduction = "umap", group.by = c('ident'),split.by=c('orig.ident'),pt.size = 0.8,label=T,order = T)+ggtitle("Heart_all")
dev.off()
###细胞比例
library(ggalluvial)
Cellratio <- prop.table(table(Idents(Heart_all), Heart_all$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)

pdf('cellratio.pdf',wi=15,he=6)
ggplot(Cellratio, aes(x =Var2, y= Freq, fill = Var1,
                  stratum=Var1, alluvium=Var1)) + 
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = mycolors)
dev.off()

subobj <- subset(Heart_all, downsample = 300)
pdf("marker_top5marerheatmap.pdf",width=15,height=15)
DoHeatmap(subobj, features = top5marker$gene)+scale_fill_gradientn(colors= c("navy","white","firebrick3"))
dev.off()


#####选择marker

celltype_marker=c(
"cdh5",'spock3','si:ch73-86n18.1','mb','igfbp6b',#Endocardium(Ventricle)
"vcam1b",'tmem88b','slc6a4a',#Endocardium(Atrium)
"vwf",'frzb','si:ch211-153b23.5',#Endocardium(frzb)

"kdrl",'plvapa','plvapb','pecam1','wu:fj16a03',#Endothelial cells(plvapb)
"cxcl12a","lyve1a","CU929150.1",#Endothelial cells(lyve1)

"acta2",'tagln','myh11a','rgs5a','tpm1',#Smooth muscle cells
"ckba","vim","fbln5",#Smooth muscle cells

"myl7","myh7",'tnnt2a','tnni4a','cmlc1',#cardiomyocytes(Ventricle)
"myh6","nppa","tnnc1b","tcap","tnni1b",#cardiomyocytes(Atrium)
"ttn.2","ttn.1",#cardiomyocytes

"tcf21","gstm.3","mmp2","tfa",#Epicardium

"fn1b",'col1a1a','mfap5','dcn','c4',#Fibroblasts
"zgc:153704","abi3bpb","angptl7","cyp26b1","fibinb",#Valve, Fibroblasts
"txn","serpine1","sele",#Endocardium(serpine1)

"cd74a","mpeg1.1","ccl35.1","c1qc",#Macrophages
"mfap4","c1qb","grn1","grn2","lygl1",#Macrophages
"sla2","ccl36.1","si:ch211-214p16.1",'ccl34b.4',#T-cells
"lyz","mpx","fcer1gl","mmp13a","lect2l",#Neutrophils

"spink2.2","spink2.4","spink2.5"#spinkgene
)

#############vlnplot
library(reshape2)
Heart_all$celltype<-factor(Heart_all$celltype,levels = sort(unique(Heart_all$celltype)))
Idents(Heart_all)="celltype"

vln.df<-as.data.frame(Heart_all[["RNA"]]@data[celltype_marker,])
vln.df$gene<-rownames(vln.df)
vln.df<-melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]<-c("CB","exp")

anno<-Heart_all@meta.data[c("celltype")]
anno$CB<-rownames(anno)
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = celltype_marker) #为了控制画图的基因顺序

pdf("heartmarker_vlnplot.pdf",width=10,height = 40)
vln.df%>%ggplot(aes(celltype,exp),cols=mycolors)+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
dev.off()

pdf("heartmarker_dotplot.pdf",width=20,height = 10)
DotPlot(object = Heart_all, 
        features=celltype_marker, 
        assay = "RNA") + theme(axis.text.x = element_text(angle = 45, 
                                                          vjust = 0.5, hjust=0.5))
dev.off()

library(FlexDotPlot)
dp = DotPlot(Heart_all, features = celltype_marker) + RotatedAxis()
pdf('newdot.pdf',wi=60,he=20)
dot_plot(dp$data[,c(3,4,1,2,5)],
         size_var = "pct.exp", 
         col_var = "avg.exp.scaled",
         size_legend = "Percent Expressed", 
         col_legend = "Average Expression",
         x.lab.pos = "bottom", 
         display_max_sizes = F, 
         hclust_method = "ward.D2",
         dend_x_var = c("pct.exp", "avg.exp.scaled"),
         dend_y_var = c("pct.exp", "avg.exp.scaled"),
         text.size = 2.5,
         cols.use = c('#330066','#336699','#66CC66','#FFCC33'),
         x.lab.size.factor = 2.5,
         y.lab.size.factor = 2.5)
dev.off()

dp = DotPlot(Heart_all, features = top5marker$gene) + RotatedAxis()
pdf('newdot_top5marker.pdf',wi=60,he=20)
dot_plot(dp$data[,c(3,4,1,2,5)],
         size_var = "pct.exp", 
         col_var = "avg.exp.scaled",
         size_legend = "Percent Expressed", 
         col_legend = "Average Expression",
         x.lab.pos = "bottom", 
         display_max_sizes = F, 
         hclust_method = "ward.D2",
         dend_x_var = c("pct.exp", "avg.exp.scaled"),
         dend_y_var = c("pct.exp", "avg.exp.scaled"),
         text.size = 2.5,
         cols.use = c('#330066','#336699','#66CC66','#FFCC33'),
         x.lab.size.factor = 2.5,
         y.lab.size.factor = 2.5)
dev.off()

#############Featureplot
library(scCustomize)
pdf('Featureplot.pdf',wi=50,he=320)
FeaturePlot_scCustom(Heart_all, features = celltype_marker, split.by = "orig.ident", num_columns = 11)
dev.off()

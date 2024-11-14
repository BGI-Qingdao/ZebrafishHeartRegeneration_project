library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(RColorBrewer)

library(patchwork)
library(harmony)
library(clustree)
library(ComplexHeatmap)

###merge初步质控后的数据
H_0dpa <- readRDS("H_0dpa.Rds")
H_BA <- readRDS(".H_BA.Rds")
H_ventricle <- readRDS("H_ventricle.Rds")
H_atrium <- readRDS("H_atrium.Rds")        
list<-list(H_0dpa=H_0dpa,H_BA=H_BA,H_ventricle=H_ventricle,H_atrium=H_atrium)


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
H_0dpa_atrium_ventricle_BA_merge <- merge(list[[1]], list[2:length(list)])
saveRDS(H_0dpa_atrium_ventricle_BA_merge, file = "beforeqc_H_0dpa_atrium_ventricle_BA_merge.Rds")
#查看每个样本的细胞数
table(H_0dpa_atrium_ventricle_BA_merge@meta.data$orig.ident) 
#H_0dpa    H_atrium        H_BA H_ventricle 
#55901       38333       26831       37874 
###画小提琴图
theme.set2=theme(axis.title.x=element_blank())
###设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.HB", "percent.rp")
group = "orig.ident"
###质控前小提琴图
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(H_0dpa_atrium_ventricle_BA_merge, group.by=group, pt.size = 0,features = plot.featrures[i])+ theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2) 
ggsave("before_QC.pdf", plot = violin,width = 9, height = 8)

qclist<-list(H_0dpa=H_0dpa,H_BA=H_BA,H_ventricle=H_ventricle,H_atrium=H_atrium)
###再次质控
#H_0dpa_atrium_ventricle_BA_merge <- subset(H_0dpa_atrium_ventricle_BA_merge,subset = percent.HB < 1)
#查看每个样本的细胞数
table(H_0dpa_atrium_ventricle_BA_merge@meta.data$orig.ident) 
#H_0dpa    H_atrium        H_BA H_ventricle 
#47776       15688       18139       36799 
###质控后小提琴图
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(H_0dpa_atrium_ventricle_BA_merge, group.by=group, pt.size = 0,features = plot.featrures[i])+ theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2) 
ggsave("after_QC.pdf", plot = violin,width = 9, height = 8)
saveRDS(H_0dpa_atrium_ventricle_BA_merge, file = "afterqc_H_0dpa_atrium_ventricle_BA_merge.Rds")

H_0dpa_atrium_ventricle_BA_merge<-readRDS('afterqc_H_0dpa_atrium_ventricle_BA_merge.Rds')
afterqc.list <- SplitObject(H_0dpa_atrium_ventricle_BA_merge, split.by = "orig.ident")
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
pdf('venn_hvgs.pdf')
venn::venn(hvgs_per_dataset, opacity = 0.5, zcolor = (scales::hue_pal())(9), cexsn = 1, 
    cexil = 1, lwd = 1, col = "white", frame = F, borders = NA)
dev.off()


###merge后的数据降维聚类
H_0dpa_atrium_ventricle_BA_merge <- NormalizeData(H_0dpa_atrium_ventricle_BA_merge) %>%FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
###查看高可变基因
top20 <- head(VariableFeatures(H_0dpa_atrium_ventricle_BA_merge), 20)
plot1 <- VariableFeaturePlot(H_0dpa_atrium_ventricle_BA_merge)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
pdf("Top20_VariableFeaturePlot.pdf",width=10,height=8)
plot2
dev.off()

H_0dpa_atrium_ventricle_BA_merge<- RunPCA(H_0dpa_atrium_ventricle_BA_merge, verbose = F)
pdf('ElbowPlot.pdf')
ElbowPlot(H_0dpa_atrium_ventricle_BA_merge,ndims=50)
dev.off()
H_0dpa_atrium_ventricle_BA_merge <- H_0dpa_atrium_ventricle_BA_merge %>% RunTSNE(dims=1:25) %>% RunUMAP(dims=1:25)
H_0dpa_atrium_ventricle_BA_merge <- FindNeighbors(H_0dpa_atrium_ventricle_BA_merge, dims=1:25) %>% FindClusters()

pdf('merge_pca_umap_tsne.pdf',wi=20,he=6)
plot_grid(ncol = 3,
  DimPlot(H_0dpa_atrium_ventricle_BA_merge, reduction = "pca", group.by = "orig.ident",pt.size = 1)+NoAxes()+ggtitle("PCA raw_data"),
  DimPlot(H_0dpa_atrium_ventricle_BA_merge, reduction = "tsne", group.by = "orig.ident",pt.size = 1)+NoAxes()+ggtitle("tSNE raw_data"),
  DimPlot(H_0dpa_atrium_ventricle_BA_merge, reduction = "umap", group.by = "orig.ident",pt.size = 1)+NoAxes()+ggtitle("UMAP raw_data"))
dev.off()
pdf('merge_everysample.pdf',wi=16,he=5)
DimPlot(H_0dpa_atrium_ventricle_BA_merge, group.by = "orig.ident",split.by = "orig.ident",pt.size = 1)
dev.off()

###正式降维聚类
tree <- FindClusters(H_0dpa_atrium_ventricle_BA_merge,resolution = c(seq(0,1,0.1)))
clustree.out <- clustree(tree)
pdf('clustree.pdf',width=20,height=20)
print(clustree.out)
dev.off()

H_0dpa_atrium_ventricle_BA_merge<- FindClusters(H_0dpa_atrium_ventricle_BA_merge, resolution = 0.2)
write.table(table(H_0dpa_atrium_ventricle_BA_merge@active.ident),file="H_0dpa_atrium_ventricle_BA_merge_cellnumbers.xls")
write.table(table(H_0dpa_atrium_ventricle_BA_merge@active.ident),file="cellnumbers.txt")

pdf('umap_tsne_cluster.pdf',wi=25,he=20)
plot_grid(ncol = 1,
  DimPlot(H_0dpa_atrium_ventricle_BA_merge, reduction = "tsne", group.by = c('ident',"orig.ident"),pt.size = 0.5,label=T)+NoAxes()+ggtitle("tSNE H_0dpa_atrium_ventricle_BA_merge"),
  DimPlot(H_0dpa_atrium_ventricle_BA_merge, reduction = "umap", group.by = c('ident',"orig.ident"),pt.size = 0.5,label=T)+NoAxes()+ggtitle("UMAP H_0dpa_atrium_ventricle_BA_merge"))
dev.off()
saveRDS(H_0dpa_atrium_ventricle_BA_merge, file = "aftercluster.Rds")


test.use = "wilcox"
H_0dpa_atrium_ventricle_BA_merge.marker <- FindAllMarkers(H_0dpa_atrium_ventricle_BA_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = test.use)    
write.csv(H_0dpa_atrium_ventricle_BA_merge.marker ,file='findmarkers.csv')
top20marker <- H_0dpa_atrium_ventricle_BA_merge.marker  %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20marker,file="top20marker.csv")
top10marker <- H_0dpa_atrium_ventricle_BA_merge.marker  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10marker,file="top10marker.csv")
top5marker <- H_0dpa_atrium_ventricle_BA_merge.marker  %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(top5marker,file="top5marker.csv")


subobj <- subset(H_0dpa_atrium_ventricle_BA_merge, downsample = 300)
pdf("top10marerheatmap.pdf",width=15,height=15)
DoHeatmap(subobj, features = top10marker$gene)+scale_fill_gradientn(colors= c("navy","white","firebrick3"))
dev.off()





mycolors <- c('#A6CEE3FF','#1F78B4FF', '#B2DF8AFF', '#33A02CFF' ,'#FB9A99FF' ,'#E31A1CFF' ,'#FDBF6FFF' ,'#FF7F00FF','#CAB2D6FF' ,'#6A3D9AFF' ,'#FFFF99FF' ,'#B15928FF','#803620FF','#E5D2DD',  '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658','#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')

pdf("top5_dotplot.pdf",width=30,height = 12)
DotPlot(H_0dpa_atrium_ventricle_BA_merge,features =unique(top5marker$gene),cols=mycolors) + RotatedAxis()
dev.off()


new.cluster.ids <- c("Endocardium(Ventricle)","Red cells","Endocardium(Atrium)","Smooth muscle cells","Cardiomyocytes(Ventricle)","Fibroblasts","T-cells","Endocardium(frzb)","Macrophages1","Cardiomyocytes(Atrium)","Endothelial cells(plvapb)", "Valve_Fibroblasts", "Macrophages2", "Macrophages3","Epicardium", "Endothelial cells(lyve1)", "Endocardium(serpine1)", "Endothelial cells(apnln)", "Neutrophils")

names(new.cluster.ids) <- levels(H_0dpa_atrium_ventricle_BA_merge) 
H_0dpa_atrium_ventricle_BA_merge <- RenameIdents(H_0dpa_atrium_ventricle_BA_merge, new.cluster.ids)      
H_0dpa_atrium_ventricle_BA_merge$celltype <- Idents(H_0dpa_atrium_ventricle_BA_merge)

saveRDS(H_0dpa_atrium_ventricle_BA_merge, file = "celltype.Rds")


pdf('marker_ident_umap.pdf',wi=18,he=10)
DimPlot(H_0dpa_atrium_ventricle_BA_merge,cols= mycolors , reduction = "umap", group.by = c('ident'),pt.size = 0.8,label=T)+ggtitle("H_0dpa_atrium_ventricle_BA_merge")
dev.off()
pdf('marker_ident_umap_sample.pdf',wi=60,he=15)
DimPlot(H_0dpa_atrium_ventricle_BA_merge,cols= mycolors , reduction = "umap", group.by = c('ident'),split.by=c('orig.ident'),pt.size = 0.8,label=T)+ggtitle("H_0dpa_atrium_ventricle_BA_merge")
dev.off()
###细胞比例
library(ggalluvial)
Cellratio <- prop.table(table(Idents(H_0dpa_atrium_ventricle_BA_merge), H_0dpa_atrium_ventricle_BA_merge$orig.ident), margin = 2)
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
#"ttn.2","ttn.1",#cardiomyocytes
"tcf21","gstm.3","mmp2","tfa",#Epicardium
"fn1b",'col1a1a','mfap5','dcn','c4',#Fibroblasts
"zgc:153704","abi3bpb","angptl7","cyp26b1","fibinb",#Valve, Fibroblasts
"txn","serpine1","sele",#valvelike
"cd74a","mpeg1.1","ccl35.1","c1qc",#Macrophages
"mfap4","c1qb","grn1","grn2","lygl1",#Macrophages
"sla2","ccl36.1","si:ch211-214p16.1",'ccl34b.4',#T-cells
"lyz","mpx","fcer1gl","mmp13a","lect2l"#Neutrophils
#"spink2.2","spink2.4","spink2.5"#spinkgene
)

#############vlnplot
library(reshape2)
H_0dpa_atrium_ventricle_BA_merge$celltype<-factor(H_0dpa_atrium_ventricle_BA_merge$celltype,levels = sort(unique(H_0dpa_atrium_ventricle_BA_merge$celltype)))
Idents(H_0dpa_atrium_ventricle_BA_merge)="celltype"

vln.df<-as.data.frame(H_0dpa_atrium_ventricle_BA_merge[["RNA"]]@data[celltype_marker,])
vln.df$gene<-rownames(vln.df)
vln.df<-melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]<-c("CB","exp")

anno<-H_0dpa_atrium_ventricle_BA_merge@meta.data[c("celltype")]
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
DotPlot(object = H_0dpa_atrium_ventricle_BA_merge, 
        features=celltype_marker, 
        assay = "RNA") + theme(axis.text.x = element_text(angle = 45, 
                                                          vjust = 0.5, hjust=0.5))
dev.off()

library(FlexDotPlot)
dp = DotPlot(H_0dpa_atrium_ventricle_BA_merge, features = celltype_marker) + RotatedAxis()
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
#############Featureplot
library(scCustomize)
pdf('Featureplot.pdf',wi=20,he=320)
FeaturePlot_scCustom(H_0dpa_atrium_ventricle_BA_merge, features = celltype_marker, split.by = "orig.ident", num_columns = 4)
dev.off()

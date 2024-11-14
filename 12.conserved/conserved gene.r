library(Seurat)
library(pheatmap)
library(ggplot2)

setwd('F:/02.Heart regeneration_BGI/00.NC_Review/06.Conserved gene/')
zeb<-readRDS('F:/02.Heart regeneration_BGI/01.生信分析/01.sc/00.final_version/00.regeneration/res1.0/01.final/Zebrafish_alltime_umap.rds')

gene<-read.table('Regeneration_gene.txt')
gene<-unique(gene)
colnames(gene)<-c('gene')
rownames(gene)<-gene$gene
allgeneexp <- AverageExpression(zeb,group.by = 'day',slot = 'data') 

geneexp<-merge(gene,allgeneexp,by='row.names')
geneexp<-geneexp[,-1]
rownames(geneexp)<-geneexp$gene
geneexp<-geneexp[,-1]
colnames(geneexp) <- gsub('RNA.g', 'Zeb_', colnames(geneexp))
write.csv(geneexp,'Zebrafish_reggene_day.csv')


#library(biomaRt)
#ensembl = useEnsembl(biomart="ensembl")
#listDatasets(ensembl)
#zebrafish <-  useMart("ensembl", dataset = "drerio_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
#human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
#mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 

#fish_list<-rownames(gene)
blast<-read.csv('addblass_v_lyx.csv')
three<-blast[,c('Hsa','Mmu','Dre')]
colnames(three)<-c('HGNC.symbol','MGI.symbol','gene')
genelist<-merge(gene,three,by='gene',all.x=T)
#fish2human = getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", 
#                 values = fish_list, 
#                 mart = zebrafish, 
##                 attributesL = c("hgnc_symbol"), 
#                 martL = human, uniqueRows=T)
#write.csv(fish2human,'fish2human.csv')

#fish2mouse = getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", 
#                 values = fish_list, 
#                 mart = zebrafish, 
#                 attributesL = c("mgi_symbol"), 
 #                martL = mouse, uniqueRows=T)
#write.csv(fish2mouse,'fish2mouse.csv')


#################################mouse exp
mmu<-readRDS('F:/02.Heart regeneration_BGI/00.NC_Review/02.lum/Mouse_sc_harmony.rds')

mmugene<-genelist$MGI.symbol
mmugene <- na.omit(mmugene)
mmugene <- unique(mmugene)
mmugene <- as.data.frame(mmugene)
rownames(mmugene)<-mmugene$mmugene
colnames(mmugene)<-c('gene')

mmu <- ScaleData(mmu,features=mmugene$gene)
mmuallgeneexp <- AverageExpression(mmu,group.by = 'orig.ident',slot = 'data') 

mmugeneexp<-merge(mmugene,mmuallgeneexp,by='row.names')
mmugeneexp<-mmugeneexp[,-1]
rownames(mmugeneexp)<-mmugeneexp$gene
mmugeneexp<-mmugeneexp[,-1]
colnames(mmugeneexp) <- gsub('RNA.mouse', 'Mmu_', colnames(mmugeneexp))
write.csv(mmugeneexp,'Mouse_reggene_day.csv')

##########################Human exp
Hum<-readRDS('F:/02.Heart regeneration_BGI/01.生信分析/05.HUMAN/02.SC/snRNA-seq-submission.rds')
humgene<-genelist$HGNC.symbol
humgene <- na.omit(humgene)
humgene <- unique(humgene)
humgene <- as.data.frame(humgene)
rownames(humgene)<-humgene$humgene
colnames(humgene)<-c('gene')

Hum <- ScaleData(Hum,features=humgene$gene)
humallgeneexp <- AverageExpression(Hum,group.by = 'major_labl',slot = 'data') 

humgeneexp<-merge(humgene,humallgeneexp,by='row.names')
humgeneexp<-humgeneexp[,-1]
rownames(humgeneexp)<-humgeneexp$gene
humgeneexp<-humgeneexp[,-1]
colnames(humgeneexp) <- gsub('RNA', 'Hsa_', colnames(humgeneexp))
write.csv(humgeneexp,'Human_reggene_day.csv')

##########################Fish heatmap
library(ClusterGVis)
getClusters(exp = geneexp)
cm<- clusterData(exp = geneexp,cluster.method = "mfuzz", cluster.num = 6, seed= 123)
table(cm$long.res$cluster)
p1 <- visCluster(object= cm,plot.type = "line", add.mline = TRUE)
p2 <- visCluster(object= cm,plot.type = "heatmap",column_names_rot= 60)
ggsave('mfuzz_line.png',p1,wi=10,he=6)
pdf('mfuzz_heatmap.pdf',wi=6,he=10)
p2 
dev.off()

mark <-read.table('top5marker.txt')
mark <- unique(mark)

pdf('mfuzz_heatmap_gene.pdf', wi=10,he=10)
visCluster(object = cm,
           plot.type = "heatmap",
           column_names_rot = 60,             
           markGenes = mark$V1,             
           genes.gp = c('italic',fontsize = 8)
           ) 
dev.off()

pdf("mfuzz_heatmap_line.pdf", width=12, height=10)
visCluster(object = cm
           ,plot.type = "both"
           ,column_names_rot = 60
           , line.side = "left"
           , add.mline = TRUE
           , markGenes = mark$V1
           ,genes.gp = c('italic',fontsize = 8)
           , show_row_dend = F)
dev.off()

library(clusterProfiler)
library(org.Dr.eg.db)

write.csv(cm$wide.res,'reggene_cluster.csv')
g<-cm$wide.res[,c(1,10)]
ENTREZID <- bitr(g$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = 'org.Dr.eg.db')

colnames(g)<-c('SYMBOL','cluster')
g1 <- merge(ENTREZID,g,by='SYMBOL')

enrich_go <- data.frame()
for (i in 1:6) {
  p<- enrichGO(gene= g1$ENTREZID[which(g1$cluster==i)]
              ,OrgDb= org.Dr.eg.db
              ,keyType= 'ENTREZID' 
              ,ont= "BP"
              ,pAdjustMethod = "BH"
              ,pvalueCutoff  = 0.05 
              ,qvalueCutoff  = 0.05
  )
  pp <- data.frame(id=paste0("C",i),term=p[,2],pval=p[,5])[c(1:5),]
  enrich_go <- rbind(enrich_go,pp)
}
write.csv(enrich_go,'reggene_cluster_enrich_go.csv')

pdf('mfuzz_heatmap_gene_GO.pdf',height = 13,width = 16)
visCluster(object = cm
           ,plot.type = "both"
           ,column_names_rot = 60
           , line.side = "left"
           , add.mline = TRUE
           , annoTerm.data = enrich_go
           , markGenes = mark$V1
           #,ctAnno.col = ggsci::pal_npg()(6)
           ,genes.gp = c('italic',fontsize = 8)
           , show_row_dend = F
           , go.col = rep(ggsci::pal_npg()(6),each = 5) #直接指定颜色
           , go.size = "pval"    #大小映射P值
           , add.bar = T, textbar.pos = c(0.8,0.2)
           , annoTerm.mside = "left"
            )
dev.off()

zeb<-readRDS('F:/02.Heart regeneration_BGI/00.NC_Review/07.Figure/单细胞注释更新/单细胞注释更新/st_merge.rds')
gene<-read.csv('../reggene_cluster.csv')
gene<-subset(gene,cluster == 6)
gene<-gene$gene
gene<-as.data.frame(gene)
colnames(gene)<-c('gene')
rownames(gene)<-gene$gene

allgeneexp <- AverageExpression(zeb,group.by = 'stage',slot = 'data') 

geneexp<-merge(gene,allgeneexp,by='row.names')
geneexp<-geneexp[,-1]
rownames(geneexp)<-geneexp$gene
geneexp<-geneexp[,-1]
colnames(geneexp) <- gsub('Spatial.', 'Zeb_', colnames(geneexp))
write.csv(geneexp,'Zebrafish_reggene_day_st.csv')

getClusters(exp = geneexp)
cm<- clusterData(exp = geneexp,cluster.method = "mfuzz", cluster.num = 2, seed= 123)
table(cm$long.res$cluster)
p1 <- visCluster(object= cm,plot.type = "line", add.mline = TRUE)
p2 <- visCluster(object= cm,plot.type = "heatmap",column_names_rot= 60)
ggsave('mfuzz_line_st.png',p1,wi=10,he=6)
pdf('mfuzz_heatmap_st.pdf',wi=6,he=10)
p2 
dev.off()
mark <-read.table('../top5marker.txt')
mark <- unique(mark)
pdf("mfuzz_heatmap_line_st.pdf", width=12, height=3)
visCluster(object = cm
           ,plot.type = "both"
           ,column_names_rot = 60
           , line.side = "left"
           , add.mline = TRUE
           , markGenes = mark$V1
           ,genes.gp = c('italic',fontsize = 8)
           , show_row_dend = F)
dev.off()

##########################Mouse heatmap
mmugeneexp <-read.csv('Mouse_reggene_day.csv',row.names=1)
mmugeneexp<-mmugeneexp[,c(2,1,4,3,6,5,8,7)]
getClusters(exp = mmugeneexp)

mcm<- clusterData(exp = mmugeneexp,cluster.method = "mfuzz", cluster.num = 10, seed= 123)
table(mcm$long.res$cluster)
p3 <- visCluster(object= mcm,plot.type = "line", add.mline = TRUE)
p4 <- visCluster(object= mcm,plot.type = "heatmap",column_names_rot= 60)
ggsave('mmu_mfuzz_line.png',p3,wi=15,he=10)
pdf('mmu_mfuzz_heatmap.pdf',wi=6,he=10)
p4
dev.off()
write.csv(mcm$wide.res,'mmu_reggene_cluster.csv')

colnames(mark) <- c('gene')
mmumark <- merge(mark,three,by='gene')

pdf("mmu_mfuzz_heatmap_line_3gene.pdf", width=12, height=10)
visCluster(object = mcm
           ,plot.type = "both"
           ,column_names_rot = 60
           , line.side = "left"
           , add.mline = TRUE
           , markGenes = mmumark
           ,genes.gp = c('italic',fontsize = 8)
           , show_row_dend = F)
dev.off()

library(org.Mm.eg.db)
g2<-mcm$wide.res[,c(1,10)]
ENTREZID2 <- bitr(g2$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = 'org.Mm.eg.db')

colnames(g2)<-c('SYMBOL','cluster')
g3 <- merge(ENTREZID2,g2,by='SYMBOL')

enrich_go <- data.frame()
for (i in 1:10) {
  p<- enrichGO(gene= g3$ENTREZID[which(g3$cluster==i)]
              ,OrgDb= org.Mm.eg.db
              ,keyType= 'ENTREZID' 
              ,ont= "BP"
              ,pAdjustMethod = "BH"
              ,pvalueCutoff  = 0.05 
              ,qvalueCutoff  = 0.05
  )
  pp <- data.frame(id=paste0("C",i),term=p[,2],pval=p[,5])[c(1:5),]
  enrich_go2 <- rbind(enrich_go,pp)
}
write.csv(enrich_go2,'mmu_reggene_cluster_enrich_go.csv')

#pdf('mmu_mfuzz_heatmap_GO.pdf',height = 13,width = 13)
#if (is.null(mcm) || nrow(enrich_go2) == 0) {
#    stop("请确保 mcm 和 enrich_go2 都不为空")
#}

#visCluster(object = mcm,
#           plot.type = "both",
#          column_names_rot = 60,
#           line.side = "left",
#           add.mline = TRUE,
#           annoTerm.data = enrich_go2,
#           show_row_dend = FALSE,
#           go.col = ggsci::pal_npg()(length(unique(enrich_go2$term))), # 确保颜色足够
#           go.size = "pval",    # 大小映射 P 值
#           add.bar = TRUE,
#           textbar.pos = c(0.8, 0.2),
#           annoTerm.mside = "left"
#)
#dev.off()

##########################Huamn heatmap
humgeneexp <-read.csv('Human_reggene_day.csv',row.names=1)
humgeneexp<-humgeneexp[,c(2,5,1,4,3)]
getClusters(exp = humgeneexp)

hcm<- clusterData(exp = humgeneexp,cluster.method = "mfuzz", cluster.num =4, seed= 123)
table(hcm$long.res$cluster)
p5 <- visCluster(object= hcm,plot.type = "line", add.mline = TRUE)
p6 <- visCluster(object= hcm,plot.type = "heatmap",column_names_rot= 60)
ggsave('hsa_mfuzz_line.png',p5,wi=15,he=10)
pdf('hsa_mfuzz_heatmap.pdf',wi=6,he=10)
p6
dev.off()
write.csv(hcm$wide.res,'hsa_reggene_cluster.csv')

pdf("hsa_mfuzz_heatmap_line.pdf", width=8, height=10)
visCluster(object = hcm
           ,plot.type = "both"
           ,column_names_rot = 60
           , line.side = "left"
           , add.mline = TRUE
           #, markGenes = mmumark$HGNC.symbol
           #,genes.gp = c('italic',fontsize = 8)
           , show_row_dend = F)
dev.off()

library(org.Hs.eg.db)
g3<-hcm$wide.res[,c(1,7)]
ENTREZID3 <- bitr(g3$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = 'org.Hs.eg.db')

colnames(g3)<-c('SYMBOL','cluster')
g4 <- merge(ENTREZID3,g3,by='SYMBOL')

enrich_go <- data.frame()
for (i in 1:6) {
  p<- enrichGO(gene= g4$ENTREZID[which(g4$cluster==i)]
              ,OrgDb= org.Hs.eg.db
              ,keyType= 'ENTREZID' 
              ,ont= "BP"
              ,pAdjustMethod = "BH"
              ,pvalueCutoff  = 0.05 
              ,qvalueCutoff  = 0.05
  )
  pp <- data.frame(id=paste0("C",i),term=p[,2],pval=p[,5])[c(1:5),]
  enrich_go3 <- rbind(enrich_go,pp)
}
write.csv(enrich_go3,'hsa_reggene_cluster_enrich_go.csv')



#########################Mouse heatmap
setwd('F:/02.Heart regeneration_BGI/00.NC_Review/06.Conserved gene/02.mouse/')
mmugeneexp <-read.csv('../Mouse_reggene_day.csv',row.names=1)
mmugeneexp<-mmugeneexp[,c(2,1,4,3,6,5,8,7)]

mmugeneexp1<-mmugeneexp[,c(1,2,5,6)]
getClusters(exp = mmugeneexp1)
mcm<- clusterData(exp = mmugeneexp1,cluster.method = "mfuzz", cluster.num = 6, seed= 123)
table(mcm$long.res$cluster)
p3 <- visCluster(object= mcm,plot.type = "line", add.mline = TRUE)
p4 <- visCluster(object= mcm,plot.type = "heatmap",column_names_rot= 60)
ggsave('mmu_mfuzz_line1.png',p3,wi=15,he=10)
pdf('mmu_mfuzz_heatmap1.pdf',wi=6,he=10)
p4
dev.off()
write.csv(mcm$wide.res,'mmu_reggene_cluster1.csv')

#colnames(mark) <- c('gene')
#mmumark <- merge(mark,three,by='gene')

pdf("mmu_mfuzz_heatmap_line_3gene1.pdf", width=6, height=10)
visCluster(object = mcm
           ,plot.type = "both"
           ,column_names_rot = 60
           , line.side = "left"
           , add.mline = TRUE
           #, markGenes = mmumark
           #,genes.gp = c('italic',fontsize = 8)
           , show_row_dend = F)
dev.off()


mmugeneexp8<-mmugeneexp[,c(3,4,7,8)]
getClusters(exp = mmugeneexp8)
mcm<- clusterData(exp = mmugeneexp8,cluster.method = "mfuzz", cluster.num = 7, seed= 123)
table(mcm$long.res$cluster)
p3 <- visCluster(object= mcm,plot.type = "line", add.mline = TRUE)
p4 <- visCluster(object= mcm,plot.type = "heatmap",column_names_rot= 60)
ggsave('mmu_mfuzz_line8.png',p3,wi=15,he=10)
pdf('mmu_mfuzz_heatmap8.pdf',wi=6,he=10)
p4
dev.off()
write.csv(mcm$wide.res,'mmu_reggene_cluster8.csv')

#colnames(mark) <- c('gene')
#mmumark <- merge(mark,three,by='gene')

pdf("mmu_mfuzz_heatmap_line_3gene8.pdf", width=6, height=10)
visCluster(object = mcm
           ,plot.type = "both"
           ,column_names_rot = 60
           , line.side = "left"
           , add.mline = TRUE
           #, markGenes = mmumark
           #,genes.gp = c('italic',fontsize = 8)
           , show_row_dend = F)
dev.off()


setwd('F:/02.Heart regeneration_BGI/00.NC_Review/06.Conserved gene/02.mouse/')
mmugeneexp <-read.csv('../Mouse_reggene_day.csv',row.names=1)
mmugeneexp<-mmugeneexp[,c(2,1,4,3,6,5,8,7)]

mmu1<-mmugeneexp[,c(1,2,5,6)]
#mmu1_sub <- subset(mmu1, mmu1[,1] < mmu1[,2] & (mmu1[,3] >= mmu1[,4] | mmu1[,3] >= 0.8 * mmu1[,4]))
mmu1_sub <- subset(mmu1, mmu1[,1] < 0.8 * mmu1[,2])
mmu1_sub2 <- subset(mmu1_sub, mmu1_sub[,3] > mmu1_sub[,4] | abs(mmu1_sub[,3] - mmu1_sub[,4]) < 0.15*mmu1_sub[,3] | abs(mmu1_sub[,3] - mmu1_sub[,4]) < 0.15*mmu1_sub[,4])
write.csv(mmu1_sub2,'mmu_1dpagene.csv')

mmu3<-mmugeneexp[,c(3,4,7,8)]
#mmu3_sub <- subset(mmu3,mmu3[,1] < mmu3[,2] & mmu3[,3] >= mmu3[,4])
#mmu3_sub <- subset(mmu3, mmu3[,1] < mmu3[,2] & (mmu3[,3] >= mmu3[,4] | mmu3[,3] >= 0.8 * mmu3[,4]))
mmu3_sub <- subset(mmu3, mmu3[,1] < 0.8*mmu3[,2])
mmu3_sub2 <- subset(mmu3_sub, mmu3_sub[,3] > mmu3_sub[,4]) 
#| abs(mmu3_sub[,3] - mmu3_sub[,4]) < 0.15*mmu3_sub[,3] | abs(mmu3_sub[,3] - mmu3_sub[,4]) < 0.15*mmu3_sub[,4])
write.csv(mmu3_sub2,'mmu_3dpagene.csv')

mmuall <- subset(mmugeneexp,0.8 *mmugeneexp[,1] > 0.8 *mmugeneexp[,5] & mmugeneexp[,2] > 0.8 *mmugeneexp[,6] & mmugeneexp[,3] > 0.8 *mmugeneexp[,7] & mmugeneexp[,4] > mmugeneexp[,8])

gene1 <- rownames(mmu1_sub2)
gene3 <- rownames(mmu3_sub2)
gene2 <- rownames(mmuall)
mousegene <- c(gene1,gene3,gene2)
mousegene <- unique(mousegene)
write.csv(mousegene,'mousegene_list.csv')

humangene<-read.csv('humangene.csv')
humangene<-humangene$X
Mousegene <- toupper(mousegene)  
intergene <- intersect(Mousegene,humangene)
write.csv(intergene,'intergene.csv')

############################################merged gene heatmap
setwd('F:/02.Heart regeneration_BGI/00.NC_Review/06.Conserved gene/03.gene/')
zeb<-readRDS('F:/02.Heart regeneration_BGI/00.NC_Review/07.Figure/单细胞注释更新/单细胞注释更新/sc_merge.rds')
fishgene<-read.csv('fishgene.txt')

zeb<-ScaleData(zeb,features=fishgene$gene)
gene_cell_exp <- AverageExpression(zeb,group.by = 'stage',slot = 'scale.data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene<-as.data.frame(fishgene)
rownames(gene)<-gene$gene
data<-merge(gene,gene_cell_exp,by='row.names')
data<-data[,-1]
rownames(data)<-data$gene
data<-data[,-1]

my_palette <- colorRampPalette(c('#adb5c5','#3e61ab','#6bc9e8','#f5e80b','#e3080b'))
p1<-pheatmap(data,color = my_palette(100),border = F,scale='row',cluster_cols = F,cutree_rows = 3)
ggsave('Fish_gene.pdf',p1,wi=5,he=7)

zeb<-readRDS('F:/02.Heart regeneration_BGI/00.NC_Review/07.Figure/单细胞注释更新/单细胞注释更新/st_merge.rds')
fishgene<-read.csv('fishgene.txt')

zeb<-ScaleData(zeb,features=fishgene$gene)
gene_cell_exp <- AverageExpression(zeb,group.by = 'stage',slot = 'scale.data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$Spatial)
gene<-as.data.frame(fishgene)
rownames(gene)<-gene$gene
data<-merge(gene,gene_cell_exp,by='row.names')
data<-data[,-1]
rownames(data)<-data$gene
data<-data[,-1]

my_palette <- colorRampPalette(c('#adb5c5','#3e61ab','#6bc9e8','#f5e80b','#e3080b'))
p1<-pheatmap(data,color = my_palette(100),border = F,scale='row',cluster_cols = F,cutree_rows = 3)
ggsave('Fish_st_gene.pdf',p1,wi=5,he=7)




mmu<-readRDS('F:/02.Heart regeneration_BGI/00.NC_Review/02.lum/Mouse_sc_harmony.rds')
mmugene<-read.csv('mousegene.txt')

mmu$day<-mmu$orig.ident
mmu$day[mmu$orig.ident == "mouse_P1_1MI"] <- "P1"
mmu$day[mmu$orig.ident == "mouse_P1_3MI"] <- "P1"
mmu$day[mmu$orig.ident == "mouse_P1_1Sham"] <- "P1"
mmu$day[mmu$orig.ident == "mouse_P1_3Sham"] <- "P1"

mmu$day[mmu$orig.ident == "mouse_P8_1MI"] <- "P8"
mmu$day[mmu$orig.ident == "mouse_P8_3MI"] <- "P8"
mmu$day[mmu$orig.ident == "mouse_P8_1Sham"] <- "P8"
mmu$day[mmu$orig.ident == "mouse_P8_3Sham"] <- "P8"

mmu$group<-mmu$orig.ident
mmu$group[mmu$orig.ident == "mouse_P1_1MI"] <- "P1_1MI"
mmu$group[mmu$orig.ident == "mouse_P1_3MI"] <- "P1_3MI"
mmu$group[mmu$orig.ident == "mouse_P1_1Sham"] <- "P1_1Sham"
mmu$group[mmu$orig.ident == "mouse_P1_3Sham"] <- "P1_3Sham"

mmu$group[mmu$orig.ident == "mouse_P8_1MI"] <- "P8_1MI"
mmu$group[mmu$orig.ident == "mouse_P8_3MI"] <- "P8_3MI"
mmu$group[mmu$orig.ident == "mouse_P8_1Sham"] <- "P8_1Sham"
mmu$group[mmu$orig.ident == "mouse_P8_3Sham"] <- "P8_3Sham"

mmu<-ScaleData(mmu,features=mmugene$gene)
gene_cell_exp <- AverageExpression(mmu,group.by = 'group',slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene<-as.data.frame(mmugene)
rownames(gene)<-gene$gene
data<-merge(gene,gene_cell_exp,by='row.names')
data<-data[,-1]
rownames(data)<-data$gene
data<-data[,-1]
data<-data[,c(2,1,6,5,4,3,8,7)]

my_palette <- colorRampPalette(c('#adb5c5','#3e61ab','#6bc9e8','#f5e80b','#e3080b'))
p1<-pheatmap(data[,1:4],color = my_palette(100),border = F,scale='row',cluster_cols = F,gaps_col = 4)
p2<-pheatmap(data[,5:8],color = my_palette(100),border = F,scale='row',cluster_cols = F,gaps_col = 4)
ggsave('Mouse_gene_1.pdf',p1,wi=4,he=7)
ggsave('Mouse_gene_3.pdf',p2,wi=4,he=7)

hum<-readRDS('F:/02.Heart regeneration_BGI/01.生信分析/05.HUMAN/02.SC/snRNA-seq-submission.rds')
humgene<-read.csv('humangene.txt')

hum<-ScaleData(hum,features=humgene$gene)
gene_cell_exp <- AverageExpression(hum,group.by = 'major_labl',slot = 'scale.data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene<-as.data.frame(humgene)
rownames(gene)<-gene$gene
data<-merge(gene,gene_cell_exp,by='row.names')
data<-data[,-1]
rownames(data)<-data$gene
data<-data[,-1]
data<-data[,c(2,5,1,4,3)]

my_palette <- colorRampPalette(c('#adb5c5','#3e61ab','#6bc9e8','#f5e80b','#e3080b'))
p1<-pheatmap(data,color = my_palette(100),border = F,scale='row',cluster_cols = F,gaps_col = 1)
ggsave('Human_gene.pdf',p1,wi=5,he=7)



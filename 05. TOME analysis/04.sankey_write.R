#TOME refer to
#Qiu, C., Cao, J., Martin, B.K. et al. Systematic reconstruction of cellular trajectories across mouse embryogenesis. Nat Genet 54, 328–341 (2022). https://doi.org/10.1038/s41588-022-01018-x
##lib
args = commandArgs(T)
library(dplyr); library(stringr); library(networkD3)

##
work_path = args[1]
setwd(work_path)
species = args[2]

##
TOME_sankey = function(myspecies, fontsize = 7, width = NULL, height = NULL){
  #list files
  myfile = list.files(paste0(work_path, '/output/knn_umap/')) %>% .[str_detect(string = ., 'longdata')] %>% data.frame(file = .)
  #read and rbind files
  mydf = paste0(work_path, '/output/knn_umap/', myfile$file) %>% lapply(function(file){read.table(file, sep ="\t", header=T)}) %>% data.table::rbindlist()
  colnames(mydf)[1:3] = c('from', 'to', 'value')
  write.table(mydf,paste0(work_path, "/output/trajectory/", myspecies, ".result.txt"),sep="\t",quote=F)
  mydf = subset(mydf,value > 0.25)
  write.table(mydf,paste0(work_path, "/output/trajectory/", myspecies, ".result_0.25.txt"),sep="\t",quote=F)

  Sankeynodes = mydf[,1:2] %>% c %>% unlist %>% unique %>% data.frame(name = .)
  Sankeynodes$index = 0:(nrow(Sankeynodes)-1)
  Sankeylinks = Sankeynodes %>% merge(mydf, ., by.x = 'from', by.y = 'name') %>% merge(., Sankeynodes, by.x = 'to', by.y = 'name')
  Sankeydata = Sankeylinks[,c('index.x','index.y','value')]
  #Sankeydata = Sankeydata %>% .[.$value!=0,]
  #Sankeydata = Sankeydata %>% .[.$value > 0.05,]

  p = sankeyNetwork(Links = Sankeydata, Nodes = Sankeynodes,
                    Source = 'index.x', Target = 'index.y', Value = 'value', 
                    NodeID = 'name', #NodeGroup = 'variable', 
                    fontSize = fontsize, sinksRight = F, width = width, height = height)
  saveNetwork(p, paste0(work_path, "/output/trajectory/", myspecies, "_Knn_umap_0.25.html"), selfcontained = T)
}

##
TOME_sankey(myspecies = species, fontsize = 10, width = 3000, height = 1000)

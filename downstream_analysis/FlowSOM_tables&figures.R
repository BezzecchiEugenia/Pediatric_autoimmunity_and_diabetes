#!/home/bezzecchi.eugenia/.conda/envs/eugenia_env/bin/Rscript

setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data/")


library(openCyto)
#BiocManager::install("CytoML")
library(CytoML)
library(flowWorkspace)
library(BiocManager)
library(devtools) #load it
library(reshape2)
library(stats)
library(cytolib)
library(flowCore)
library(CytoExploreR)
library(na.tools)
library(CytoML)
library(ggcyto)
library(data.table)
library(ggplot2)
library(ggridges)
library(dplyr)
library(CytoML)
library(CytoExploreR)
library("openCyto")
library(flowWorkspace)
library(flowCore)
library("CATALYST")
library(tidyverse)
library(data.table)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(data.table)
library(rstatix)
library(ggpubr)
library(moments)
library(rstatix)
library(ggpubr)
library("bestNormalize")
library(car)
library("multcomp")
library("gtools")

files=list.files(path,"sce",recursive = T)
tabella=as.data.frame(files)
tabella=as.data.frame(tabella)
colnames(tabella)[1]="files"
tabella$folder=tabella$files
tabella$folder=sub("/.*", "",tabella$folder)
tabella$folder=paste0(tabella$folder,"","/")
tabella$sce=sub('.*/', '', tabella$files)
tabella$prefix=c("Bcells","DCs","NKcells","Tcells","TREG")
tabella$clusters=c("meta20","meta20","meta40","meta20","meta15")


for (i in 1:5) {
  setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")
  folder=tabella$folder[i]
  sce=tabella$sce[i]
  setwd(folder)
  sce <- readRDS(sce)
  md=read.table("md.txt",sep='\t',header = T,row.names = 1)
  
  merging=list.files(path = ".", pattern =paste0("merging_",tabella$prefix[i],".txt"),recursive = F)
  merging_table1=read.delim(merging)
  clusters=tabella$clusters[i]
  
  sce <- mergeClusters(sce, k = clusters, 
                     table = merging_table1, id = "merging1",overwrite = T)

  
  p3 <- plotFreqHeatmap(sce, k = "merging1", perc = TRUE,
                        row_anno = FALSE, col_clust = FALSE,normalize = F)
  
  z<-p3@matrix
  z<-as.data.frame(p3@matrix)
  pop=as.factor(row.names(z))
  z<-data.frame(t(z))
  colnames(z)=pop
  z$sample_id=row.names(z)

  z=merge.data.table(z,md,by="sample_id",all.x = T)
  
  data<-z
  data$file_name=NULL
  data$condition=NULL
  data$age=NULL
  
  data_M<-melt(data,id.vars = "sample_id")
  data_M=merge(data_M,md,by="sample_id",all.x=T)
  fileName=paste0(tabella$prefix[i],"_","unsup_stats.txt")
  
  setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")
  write.table(data_M,fileName,sep="\t",quote = F,row.names = FALSE)
 
  dir.create("figures_paper")
  
  pdf("figures_paper/Tcells_Heatmap.pdf",12,3.4)
  plotExprHeatmap(sce, features = "type", scale="last",col_clust = F,
                  by = "cluster_id", k ="merging1", bin_anno = F,
                  bars = TRUE, perc = TRUE)
  dev.off()
  
  tiff("figures_paper/UMAP_Tcells.tiff",width = 12, height = 4,units = 'in',res=200)
  plotDR(sce, "UMAP", color_by ="merging1")
  dev.off()
  
  pdf("figures_paper/merging_strategy.pdf",8,9)
   plotExprHeatmap(sce, features = "type", scale="last",col_clust = F,
                by = "cluster_id", k =clusters,m= "merging1", bin_anno = F,
                bars = TRUE, perc = F)
   dev.off()
  
  
  }
  
 
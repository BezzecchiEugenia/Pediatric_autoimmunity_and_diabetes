library(openCyto)
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

setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")

pcx = 1
pcy = 2

DF=read.delim("merged_unsup_intersect.txt")

DF=DF[!grepl("Unclassified",row.names(DF)),]
DF=DF[!grepl("unclassified",row.names(DF)),]
DF=DF[!grepl("Double positive CD4 CD8 T cells",row.names(DF)),]

metadata=read.delim("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/metadata_samples.txt",sep="\t",header = T)
colnames(metadata)[1]="sample_id"


############
library(limma)
samp=as.data.frame(colnames(DF))
colnames(samp)[1]="sample_id"
samp$sample_id=gsub('\\.', '-',samp$sample_id)
tec=read.delim("batch_effect.txt")
samp=merge(samp,tec,by="sample_id",all.x=T)
samp=merge(samp,metadata,by="sample_id")
samp=samp[,-c(2,3,6)]
my_design = model.matrix(~samp$condition)
DF <- removeBatchEffect(DF, samp$batch,design=my_design)
write.table(DF,"merged_unsup_intersect_batch_removed.txt",sep="\t")
############

metadata=read.delim("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/metadata_samples.txt",sep="\t",header = T)
colnames(metadata)[1]="sample_id"

dir= 'PCA_paper'
dir.create(dir, recursive = TRUE)
setwd(dir)

metadata=metadata[,c(1,2,3)]
samples=as.data.frame(as.character(colnames(DF)))
colnames(samples)[1]="sample_id"
samples$sample_id=gsub("\\.","-", samples$sample_id)
samples=merge(samples,metadata,by="sample_id",all.x=T)
is.na(samples$condition)
centering = TRUE
scaling = TRUE

pca = prcomp(t(DF),center = centering, scale=scaling)
var = round(matrix(((pca$sdev^2)/(sum(pca$sdev^2))), ncol=1)*100,1)
score = as.data.frame(pca$x)


score$age=samples$age
score[]<-as.numeric(as.character(unlist(score, use.names = FALSE)))
score$Quartile<-ifelse(score$age<=quantile(score$age,c(0.5)),1,2)
score$age<-ifelse(score$Quartile==1,"young","old")
score$Quartile=NULL
score$factor=samples$condition

# plot paramters
xlab = paste("PC", pcx, " (",var[pcx],"%)", sep="")
ylab = paste("PC", pcy, " (",var[pcy],"%)", sep="")
cum = var[pcx]+var[pcy]
names = rownames(pca$x)

library("scales")
library(ggsci)
pal_jco(palette = c("default"), alpha = 1)
show_col(pal_jco("default", alpha = 0.6)(10))

score$factor=factor(score$factor,levels=c("HC","CD_THY","preT1D_LR","T1D"))
mycolors_c <- c("blue","forestgreen","grey44","red");names(mycolors_c) = levels(as.factor(score$factor))

  
pdf(paste("all_cond",'_PCA_',pcx,'_',pcy,'.pdf',sep=''),width=15, height=13)

ggplot(score, aes(x=score[,pcx], y=score[,pcy], shape=factor))+
    geom_point(size=6)+
  scale_shape_manual(values=c(0,4,2,20))+
    labs(x=xlab, y=ylab, title=paste("PC",pcx," vs PC",pcy," scoreplot",sep="")) + 
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 24),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
          axis.title.y = element_text(face = "bold", color = "black", size = 24),
          legend.text = element_text(face = "bold", color = "black", size = 16),legend.title = element_text(face = "bold", color = "black", size = 24),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))+
  scale_color_manual(values = mycolors_c)  

  dev.off()
  
setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/ULTIMATE_analysis/unsupervised_last/final_stats/")

pcx = 1
pcy = 2

DF=read.delim("merged_unsup_intersect_batch_removed.txt")
DF=DF[!grepl("Unclassified",row.names(DF)),]
DF=DF[!grepl("unclassified",row.names(DF)),]
DF=DF[!grepl("Double positive CD4 CD8 T cells",row.names(DF)),]

metadata=read.delim("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/metadata_samples.txt",sep="\t",header = T)
colnames(metadata)[1]="sample_id"

META=read.delim("clinical_var_v2.txt")
is.na(META) <- META==''


DF=DF[,!grepl("OSR",colnames(DF))]
DF=DF[,!grepl("preT1D",colnames(DF))]
DF=DF[,grepl("T1D",colnames(DF))]
sample_id=colnames(DF)
sample_id=as.data.frame(sample_id)
DF=as.data.frame(t(DF))
DF$sample_id=sample_id$sample_id

colnames(META)
META=META[,c(1,5,16,29)]

dir= 'PCA_paper'
dir.create(dir, recursive = TRUE)
setwd(dir)


pca_plot <- list()

i=0
for (col in 2:4) {
  z=colnames(META[col])
  dir.create(z)
  setwd(z)
  metadata=META[,c(1,col)]
  merged=merge(DF,metadata,by="sample_id")
  colnames(merged)
  merged$sample_id=NULL
  factor=as.data.frame(merged[,47])
  colnames(factor)="factor"
  factor$factor=as.character(factor$factor)
  merged=merged[,-47]
  merged=t(merged)
  merged=data.frame(merged)
  merged[]<-as.numeric(as.character(unlist(merged, use.names = FALSE)))
  
  
centering = TRUE
scaling = TRUE

# PCA

pca = prcomp(t(merged),center = centering, scale=scaling)
var = round(matrix(((pca$sdev^2)/(sum(pca$sdev^2))), ncol=1)*100,1)
score = as.data.frame(pca$x)

score$factor=factor$factor

# plot paramters
xlab = paste("PC", pcx, " (",var[pcx],"%)", sep="")
ylab = paste("PC", pcy, " (",var[pcy],"%)", sep="")
cum = var[pcx]+var[pcy]
names = rownames(pca$x)

score$factor=factor(score$factor,levels = c("low","high"))

 i=i+1
  pca_plot[[i]] <- ggplot(score, aes(x=score[,pcx], y=score[,pcy], shape=factor))+
    geom_point(size=7)+
  scale_shape_manual(values=c(4,20))+
    labs(x=xlab, y=ylab, title=paste("PC",pcx," vs PC",pcy," scoreplot",sep="")) + 
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
   labs(color=z) +
    theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 24),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
          axis.title.y = element_text(face = "bold", color = "black", size = 24),
          legend.text = element_text(face = "bold", color = "black", size = 16),legend.title = element_text(face = "bold", color = "black", size = 24),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) #+
pca_plot[[i]]=pca_plot[[i]]+guides(linetype = guide_legend(nrow = 2)) 
  
  pdf(paste(z,'_PCA_',pcx,'_',pcy,'.pdf',sep=''),width=10, height=7)
  plot(pca_plot[[i]])
  dev.off()
  setwd(dirname(getwd()))
}


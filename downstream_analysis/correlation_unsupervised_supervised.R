setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")


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

#import files with supervised and unsupervised statistics in a long format
load("all_panels_stats.RData")

vector=c("Tcells","TREG","DCs","NK","Bcells")
list=mget(vector, envir = globalenv())

for (i in 1:5) {
list[[i]]$panel=vector[i]
}

ALL_merged=rbindlist(list)

ALL_merged=ALL_merged[!grepl("Double positive CD4 CD8 T cells",ALL_merged$Population),]

merged=ALL_merged
corr.coef<-cor(dcast(merged,sample_id+Population~method,value.var = "prop")[,c("supervised","unsupervised")],method="spearman",use="na.or.complete")
et<-element_text(size=28)
et2<-element_text(size=20)

#devtools::install_github("kevinushey/Kmisc")
library(Kmisc)

p1<-ggplot(dcast(merged,sample_id+Population~method,value.var = "prop"))+
  geom_point(aes(x=supervised,y=unsupervised,color=basename(as.character(Population))))+
  scale_x_log10(wrap("Supervised gating",30))+
  scale_y_log10(wrap("Unsupervised - FlowSOM",30))+theme_bw()+geom_abline(lty=3)+
  theme(axis.text.x=et,axis.text.y=et,axis.title.x=et,axis.title.y=et,plot.title=et,legend.text=et2,legend.position="none")+
  scale_color_discrete("Population")+geom_text(aes(x=0,y=1,size=50,label=sprintf("rho=%s",signif(corr.coef[1,2],3))),vjust = "outward",hjust = "inward",
                                               data=data.frame(corr.coef))+ggtitle("All pop")

pdf("dotplot_sup_unsup.pdf",10,10)
plot(p1)
dev.off()

data=p1[["data"]]
r=data.frame(table(data$Population))
r$Var1=as.factor(r$Var1)
dim=dim(r)[1]


list=list()
for (i in 1:dim) {
data=p1[["data"]]
Pop=as.character(r$Var1[i])
pop2=gsub("\\+",'\\\\+',Pop)
pop2=gsub("\\-",'\\\\-',pop2)
data=data[grepl(pop2,data$Population),]
cor=data.frame(cor(data$unsupervised,data$supervised,method = "spearman"))
cor$pop=Pop
colnames(cor)[1]="Spearman_correlation"
list[[i]]=cor
}

corr=rbindlist(list)
data=p1[["data"]]
list=list()
list[[1]]=data

all=list()

for (i in 1:dim(r)[1]) {
  q=as.factor(r$Var1[i])
  all[[i]]=list[[1]][list[[1]]$Population==q,]
}


all_corr=list()

library(confintr)

for (i in 1:dim(r)[1]) {
  all_corr[[i]]=ci_cor(all[[i]]$supervised,all[[i]]$unsupervised, method = "spearman", type = "bootstrap", R = 30000, seed = 1)
  all_corr[[i]]=as.data.frame(rbind(all_corr[[i]]$interval,all_corr[[i]]$estimate))
  all_corr[[i]]$cor=all_corr[[i]]$V1[2]
  colnames(all_corr[[i]])=c("min","max","cor")
  all_corr[[i]]=all_corr[[i]][-2,]
  rownames(all_corr[[i]])=all[[i]]$Population[1]
  all_corr[[i]]$Population=rownames(all_corr[[i]])
}

all_corr=rbindlist(all_corr)



pdf("barplot_sup_unsup.pdf",18,12)
p=ggplot(all_corr) +
  geom_bar( aes(x=reorder(Population,cor), y=cor), stat="identity",fill="white",color="black",alpha=0.5) +
  geom_errorbar( aes(x=Population, ymin=min, ymax=max), width=0.4, colour="red", alpha=0.9, size=1.3)# +
p=p+theme_bw()
p=p+theme(plot.title = element_text(color="black", size=26, face="bold.italic"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5), 
          axis.title.x = element_blank(),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
          axis.title.y =  element_blank(),
         legend.text = element_text(face = "bold", color = "black", size = 16),legend.title = element_text(face = "bold", color = "black", size = 24))

p=p+rotate_x_text(angle = 90)
plot(p)
dev.off()
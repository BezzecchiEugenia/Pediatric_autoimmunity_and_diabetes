setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")

library(openCyto)
#BiocManager::install("CytoML")
library(CytoML)
library(flowWorkspace)
library(BiocManager)
library(devtools) #load it
library(reshape2)
library(stats)
#library(Kmisc)
#devtools::install_github("DillonHammill/CytoExploreR")
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
#library(CytoRSuite)
library(flowWorkspace)
library(flowCore)
library("CATALYST")
#library(here)
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

stats2=ALL_merged

stats2=stats2[grepl("unsupervised",stats2$method),]
stats2=stats2[!grepl("Double positive CD4 CD8 T cells",stats2$Population),]

  list=list()
  Boxplot_list=list()
  skewNESS=list()
  skewNESS_orig=list()
  
  
  
metadata=read.delim("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/metadata_samples.txt",sep="\t",header = T)
colnames(metadata)[1]="sample_id"
stats2=merge(stats2,metadata,by="sample_id",all.x=T)


tec=read.delim("../batch_effect.txt",row.names = 1)
tec$Tecnician=NULL
colnames(tec)[1]="sample_id"
stats2=merge.data.table(stats2,tec,by="sample_id",all.x = T)
stats2$correct_sample=NULL
stats2$to_rem=NULL

stats2$sex=ifelse(stats2$sex=="F",1,0)
stats2$sex=as.factor(stats2$sex)

stats2$batch=ifelse(stats2$batch=="batch2",1,0)
stats2$batch=as.factor(stats2$batch)

stats2$panel=factor(stats2$panel,levels=c("Tcells","NK","TREG","DCs","Bcells"))

r=data.frame(table(stats2$Population))
r$Var1=as.factor(r$Var1)
r$Population=r$Var1
stats=stats2[!duplicated(stats2$Population),][,c(2,5)]
r=merge(r,stats,by="Population",all.x = TRUE, all.y = FALSE)
r$panel=factor(r$panel,levels=c("Tcells","NK","TREG","DCs","Bcells"))
r=r[order(r$panel),]

stats2$parent=ifelse(stats2$panel=="Tcells","T",ifelse(stats2$panel=="Bcells","B",ifelse(stats2$panel=="DCs","HLADR+ LIN-",ifelse(stats2$panel=="TREG","CD4 T","Lymphocytes"))))
  
dim=dim(r)[1]  
 
  for (i in 1:dim) {
  #  i=4
  Pop=as.character(r$Var1[i])
  DF=stats2[which(stats2$Population==Pop),]
  
  
  DF$condition=factor(DF$condition,levels=c("HC","CD_THY","preT1D_LR","T1D"))
  norm <- bestNormalize(DF$prop,allow_lambert_s = TRUE)
  DF$prop<- predict(norm)
  DF$orig <- predict(norm, newdata =DF$prop, inverse = TRUE)
  
  skewNESS[[i]]=skewness(DF$prop)
  names(skewNESS)[[i]]=Pop
  
  stat.test <- shapiro_test(DF$prop)
  stat.test<-paste0("Shapiro-Wilk, W : ", stat.test$statistic," pval : ",stat.test$p.value)
  
  
  shapiro.p<-signif(shapiro.test(DF$prop)$p.value,1)
  
  hist(DF$prop,main=Pop)
  mtext(paste0(stat.test),
        side=3, cex=0.8)
  

  dev.new()
  parent=DF$Parent[i]
  DF$age=as.numeric(DF$age)
  res =lm(prop~condition+age+sex+batch,data=DF)
  summary(res)
  stat.test=glht(res, linfct =mcp(condition="Tukey"))
  summary(stat.test)
  pval=round(summary(stat.test)$test$pvalues,digits=4)
  
  a=summary(stat.test)$test
  b=names(a[["coefficients"]])
  pval=data.frame(t(rbind(pval,b)))
  pval= separate(data = pval, col = b, into = c("group1", "group2"), sep = "\\ - ")
  pval$p.signif=stars.pval(as.numeric(pval$pval))
  pval$p.signif=ifelse(pval$p.signif==" ","ns",pval$p.signif)
  pval$p.signif=ifelse(pval$p.signif==".","ns",pval$p.signif)
  list[[i]]=pval
  names(list)[[i]]=Pop
  list[[i]]$pop=Pop
  
  y_pos=get_y_position(
    DF,
    orig~condition,
    fun = "max",
    comparisons = NULL,
    step.increase = 0.16,
    y.trans = NULL,
    stack = FALSE,
    scales = c("free_y")
  )
  
  
  
  pval_filt=pval[!grepl("ns",pval$p.signif),]
  
  
  
    
    if (nrow(pval_filt)!= 0) {
      pval_filt$comparison=paste0(pval_filt$group1," | ",pval_filt$group2)
      y_pos$comparison=paste0(y_pos$group2," | ",y_pos$group1)
      y=merge(pval_filt,y_pos,by="comparison")
      
      
      
      
      DF$condition=factor(DF$condition,levels=c("HC","CD_THY","preT1D_LR","T1D"))
      
      Boxplot_list[[i]]=ggboxplot(DF, x = "condition", y = "orig", fill = "gray",show.legend = FALSE,
                                  add = "jitter", title=Pop,ylab = paste("cells fraction / ",DF$parent[i]),xlab="")+rotate_x_text(angle = 45)+stat_pvalue_manual(
                                    pval,size = 8, face = "bold",
                                    y.position = y$y.position, step.increase = 0.15,
                                    label = "p.signif",hide.ns = "T")+theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18),plot.title = element_text(size = 17, face = "bold",lineheight = .8))
      
    }else{
      
      Boxplot_list[[i]]=ggboxplot(DF, x = "condition", y = "orig", fill = "gray",show.legend = FALSE,
                                  add = "jitter",title=Pop, ylab = paste("cells fraction / ",DF$parent[i]),xlab="")+rotate_x_text(angle = 45)+theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18),plot.title = element_text(size = 17, face = "bold",lineheight = .8))
    }
    
  }
 
  
  library(gridExtra)
  n <- length(Boxplot_list)
  nCol <- floor(sqrt(n))
  pdf(paste0("unsupervised","_","_stats.pdf"),25,40)
  do.call("grid.arrange", c(Boxplot_list, ncol=5))
  dev.off()
  

  
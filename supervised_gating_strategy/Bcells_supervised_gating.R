setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data/")

library(openCyto)
library(CytoML)
library(flowWorkspace)
library(BiocManager)
library(devtools) #load it
library(ggcyto)
library(reshape2)
library(stats)
library(cytolib)
library(flowCore)
library(CytoExploreR)
library(na.tools)
library(CytoML)
library(ggcyto)
library(ggridges)
library(flowStats)
library(data.table)
library(bestNormalize)
library(ggpubr)
library(moments)
library(rstatix)
library(multcomp)
library(tidyr)
library(gtools)
library(ggpubr)


setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/ULTIMATE_analysis")
gs=load_gs("Bcells_gs/")


template=data.frame()

template <- rbind(template, add_pop(gs,
                                    parent = "Bcells",
                                    alias = "Naive_B",
                                    dims="APC-Cy7-A,FITC-A",
                                    gating_method ="boundary",gating_args ="min=c(-0.5,1.6),max=c(1.9,3.4)"))

ggcyto(gs[100:130],mapping = aes(x = "APC-Cy7-A",y = "FITC-A"),subset = "Bcells") +
geom_hex(bins = 100) +geom_gate("Naive_B")+axis_x_inverse_trans()+axis_y_inverse_trans()+geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.7,size=2.5)


template <-rbind(
  template, add_pop(gs,pop="-",
                    parent = "Naive_B",
                    alias = "IGD_only_naive",
                    dims="PE-A",
                    gating_method ="mindensity",gating_args ="min=2.4,max=2.8",collapseDataForGating = T))


template <- rbind(
  template,add_pop(gs,pop="+",
                   parent = "Naive_B",
                   alias = "Naive_Bcells",
                   dims="PE-A",
                    gating_method ="mindensity",gating_args ="min=2.4,max=2.8",collapseDataForGating = T))

ggcyto(gs[250:290],mapping = aes(x = "PE-A",y = "SSC-A"),subset = "Naive_B") +
geom_hex(bins=100) + geom_gate("IGD_only_naive")+ geom_gate("Naive_Bcells")+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.9,size=2.7)


template <- rbind(template, add_pop(gs,
                                    parent = "Bcells",
                                    alias = "IgM-IgD-",
                                    dims="PE-A,FITC-A",
                                    gating_method ="boundary",gating_args ="min=c(0,0),max=c(1.5,1.5)"))


template <- rbind(template, add_pop(gs,
                                    parent = "Bcells",
                                    alias = "IgM-IgD+",
                                    dims="PE-A,FITC-A",
                                    gating_method ="boundary",gating_args ="min=c(0,1.8),max=c(1.5,3.6)"))


template <- rbind(template, add_pop(gs,
                                    parent = "Bcells",
                                    alias = "IgM+IgD+",
                                    dims="PE-A,FITC-A",
                                    gating_method ="boundary",gating_args ="min=c(1.5,1.8),max=c(3.7,3.6)"))

ggcyto(gs[100:130],mapping = aes(x = "PE-A",y = "FITC-A"),subset = "Bcells") +
geom_hex(bins = 100) +geom_gate("IgM-IgD-")+geom_gate("IgM-IgD+")+geom_gate("IgM+IgD+")+axis_x_inverse_trans()+axis_y_inverse_trans()+geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.7,size=2.5)

template <- rbind(template, add_pop(gs,
                                    parent = "IgM+IgD+",
                                    alias = "CD27-",
                                    dims="PerCP-A,APC-Cy7-A",
                                    gating_method ="boundary",gating_args ="min=c(-0.5,-0.7),max=c(4,2)"))


ggcyto(gs[100:130],mapping = aes(x = "PerCP-A",y = "APC-Cy7-A"),subset = "IgM+IgD+") +
geom_hex(bins = 100) +geom_gate("CD27-")

template <- rbind(template, add_pop(gs,
                                    parent = "CD27-",
                                    alias = "transitional_B",
                                    dims="APC-A,PerCP-A",
                                    gating_method ="flowClust",gating_args ="K=1, transitional = T, translation = 0.8,transitional_angle = pi*.35"))

ggcyto(gs[210:240],mapping = aes(x = "APC-A",y = "PerCP-A"),subset = "CD27-") +
geom_hex(bins = 100) +geom_gate("transitional_B")+axis_x_inverse_trans()+axis_y_inverse_trans()+geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.7,size=2.5)


template <- rbind(
  template,add_pop(gs,pop="+",
                   parent = "IgM-IgD-",
                   alias = "Early_PB",
                   dims="PErCP-A,APC-Cy7-A",
                   gating_method ="boundary",gating_args ="min=c(2.1,-0.4),max=c(3.5,1.3)",collapseDataForGating = T))

ggcyto(gs[1:20],mapping = aes(x = "PErCP-A",y = "APC-Cy7-A"),subset = "IgM-IgD-") +
geom_hex(bins=80)+ geom_gate("Early_PB") 

template <- rbind(
  template,add_pop(gs,pop="+",
                   parent = "IgM-IgD-",
                   alias = "PB",
                   dims="PErCP-A,APC-Cy7-A",
                   gating_method ="boundary",gating_args ="min=c(2.1,1.3),max=c(3.5,4)",collapseDataForGating = T))

ggcyto(gs[1:20],mapping = aes(x = "PErCP-A",y = "APC-Cy7-A"),subset = "IgM-IgD-") +
geom_hex(bins=80)+ geom_gate("Early_PB") + geom_gate("PB") 

template <- rbind(
  template,add_pop(gs,pop="+",
                   parent = "Bcells",
                   alias = "DN",
                   dims="APC-Cy7-A,FITC-A",
                   gating_method ="boundary",gating_args ="min=c(-0.3,0),max=c(1,1.3)",collapseDataForGating = T))

ggcyto(gs[100:130],mapping = aes(x = "APC-Cy7-A",y = "FITC-A"),subset = "Bcells") +
geom_hex(bins=100)+ geom_gate("DN") 


template <- rbind(
  template,add_pop(gs,pop="+",
                   parent = "DN",
                   alias = "DN2",
                   dims="PE-Cy7-A,APC-A",
                   gating_method ="boundary",gating_args ="min=c(0,0),max=c(1.8,1.7)",collapseDataForGating = T))


template <- rbind(template, add_pop(gs,
                                    parent = "IgM+IgD+",
                                    alias = "Unswitched_Mem_B",
                                    dims="PerCP-A,APC-Cy7-A",
                                    gating_method ="boundary",gating_args ="min=c(-0.5,2),max=c(2.2,3.5)"))


ggcyto(gs[100:130],mapping = aes(x = "PerCP-A",y = "APC-Cy7-A"),subset = "IgM+IgD+") +
geom_hex(bins = 100)+geom_gate("Unswitched_Mem_B")+axis_x_inverse_trans()+axis_y_inverse_trans()+geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.7,size=2.5)

template <- rbind(
  template,add_pop(gs,pop="-",
                   parent = "Unswitched_Mem_B",
                   alias = "IgD_only_switched",
                   dims="PE-A",
                   gating_method ="mindensity",gating_args ="min=1.7",collapseDataForGating = T))

ggcyto(gs[120:140],mapping = aes(x = "PE-A",y = "SSC-A"),subset = "Unswitched_Mem_B") +
geom_hex(bins=40)+ geom_gate("IgD_only_switched") 

template <- rbind(template, add_pop(gs,
                                    parent = "IgM-IgD-",
                                    alias = "switched_Mem_B",
                                    dims="PerCP-A,APC-Cy7-A",
                                    gating_method ="boundary",gating_args ="min=c(0,1.3),max=c(2.1,4)"))

ggcyto(gs[1:30],mapping = aes(x = "PerCP-A",y = "APC-Cy7-A"),subset = "IgM-IgD-") +
geom_hex(bins = 100)+geom_gate("switched_Mem_B")+axis_x_inverse_trans()+axis_y_inverse_trans()+geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.7,size=2.5)


pdf("Bcells_norm/switched_Mem_B_BN.pdf")
ggcyto(gs_pop_get_data(gs[1:40], "switched_Mem_B"), aes(x="APC-Cy7-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()


gs <- normalize(gs, populations=c("switched_Mem_B"), bwFac = 2,
                dims=c("APC-Cy7-A"), minCountThreshold = 50)


pdf("Bcells_norm/switched_Mem_B_AN.pdf")
ggcyto(gs_pop_get_data(gs[1:40], "switched_Mem_B"), aes(x="APC-Cy7-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()


template <- rbind(template, add_pop(gs,
                                    parent = "IgM-IgD-",
                                    alias = "switched_Mem_B",
                                    dims="PerCP-A,APC-Cy7-A",
                                    gating_method ="boundary",gating_args ="min=c(0,2),max=c(2,3.5)"))

ggcyto(gs[100:130],mapping = aes(x = "PerCP-A",y = "APC-Cy7-A"),subset = "IgM-IgD-") +
geom_hex(bins = 100)+geom_gate("switched_Mem_B")+axis_x_inverse_trans()+axis_y_inverse_trans()+geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.7,size=2.5)



template <- rbind(template, add_pop(gs,
                                    parent = "IgM-IgD+",
                                    alias = "IgD_only_switched_mem",
                                    dims="PerCP-A,APC-Cy7-A",
                                    gating_method ="boundary",gating_args ="min=c(0,1.5),max=c(2,3)"))

ggcyto(gs[100:102],mapping = aes(x = "PerCP-A",y = "APC-Cy7-A"),subset = "IgM-IgD+") +
geom_hex(bins = 40) +geom_gate("IgD_only_switched_mem")+axis_x_inverse_trans()+axis_y_inverse_trans()+geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.7,size=2.5)


dir.create("Bcells_final_statistics")
setwd("Bcells_final_statistics/")

stats2<-getPopStats(gs,statistic="count")#extract stats
stats2[,prop := Count/ParentCount] #compute the cell proportions


md <- data.frame()[1:dim(stats2), ]
md$file_name<-stats2$name
md$sample_id<-gsub("Bcells_*",'',md$file_name)
md$sample_id=gsub("_001.fcs",".fcs",md$sample_id)
md$sample_id<-gsub("*.fcs",'',md$sample_id)


metadata=read.delim("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/metadata_samples.txt",sep="\t",header = T)
md=merge(md,metadata,by="sample_id",all.x = T)

table(md$condition)
condition=as.data.frame(cbind(name=md$file_name,condition=md$condition))
condition=condition[!duplicated(condition$name),]
stats2 = merge(stats2,condition,by="name",all.x=T)
data.table::setDT(stats2)

table(is.na(stats2$condition))
stats2=na.omit(stats2)

require(dplyr)
stats2=stats2 %>%
  mutate(Population = gsub(".*/", "",Population))

stats2=stats2 %>%
  mutate(Parent = gsub(".*/", "",Parent))


pdf("boxplots_prop2.pdf",15,20)
ggplot(stats2) + geom_boxplot(outlier.shape = NA) + aes(y = prop,x = Population) + geom_jitter() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) + facet_wrap( ~
                                                                          Parent,scales = "free")
dev.off()


pdf("boxplots_coeff_variation2.pdf",15,20)
ggplot(stats2[,.(CV = mad(prop) / median(prop)),.(Parent,condition,Population)]) +
  geom_bar(position = "dodge",stat = "identity") + aes(y = CV,x = Population,fill =
                                                         condition) + theme(axis.text.x = element_text(angle = 90,hjust = 1)) + facet_wrap( ~
                                                                                                                                              Parent,scales = "free_x")


dev.off()

pop=as.data.frame(table(stats2$Population))


norm <- bestNormalize(stats2$prop,allow_lambert_s = TRUE)
stats2$prop<- predict(norm)
stats2$orig <- predict(norm, newdata =stats2$prop, inverse = TRUE)

nm=list()
for (i in 1:dim(pop)[1]) {
  # i=1
  z=as.character(pop$Var1)[i]
  nm[[i]] = stats2[,.(name,rZ = (prop - median(prop)) / mad(prop)),.(Population)][,.(outlier = abs(rZ) >
                                                                                       3),.(name,Population)][outlier == TRUE][order(name)][Population %in% z]$name
  names(nm)[[i]]=z
}




library("sjmisc")
mylist2 = nm[-which(sapply(nm,is_empty))]
nm=mylist2

dir.create("outliers_last_norm")
for (j in 1:length(nm)){
  p=names(nm)[j]
  q=nm[[j]]
  pdf(paste("outliers_last_norm/",p,".pdf",sep=""))
  a=autoplot(gs[q],p,bins=200) + stats_null() + geom_stats(type = "percent",adjust = 0.8)
  plot(a)
  dev.off()
}


plot(gs)


setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")
save_gs(gs,path="Bcells_supervised_gates_gs")
write.table(stats2,"Bcells_supervised.txt",quote=F,col.names = NA)
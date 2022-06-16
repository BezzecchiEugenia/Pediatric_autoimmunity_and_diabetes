setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")

library(openCyto)
library(CytoML)
library(flowWorkspace)
library(BiocManager)
library(devtools) 
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


gs=load_gs("NK_gs/")

template=data.frame()


template <- rbind(template, add_pop(gs,
                                    parent = "Tcells",
                                    alias = "CD4-CD8-Tcells",
                                    dims="PE-Cy7-A,Pacific Blue-A",
                                    gating_method ="boundary",gating_args ="min=c(-0.7,-0.5),max=c(2.1,1.7),quantile=0.95"))


template <- rbind(template, add_pop(gs,
                                    parent = "Tcells",
                                    alias = "CD4+CD8+Tcells",
                                    dims="PE-Cy7-A,Pacific Blue-A",
                                    gating_method ="boundary",gating_args ="min=c(2.2,1.7),max=c(4.4,3.2),quantile=0.95"))

ggcyto(gs[1:20],mapping = aes(x = "PE-Cy7-A",y = "Pacific Blue-A"),subset = "Tcells") +
geom_hex(bins=200)+geom_gate("CD8_Tcells")+geom_gate("CD4_Tcells")+geom_gate("CD4-CD8-Tcells")+geom_gate("CD4+CD8+Tcells")


ggcyto(gs[1:24],mapping = aes(x = "PE-Cy7-A",y = "Pacific Blue-A"),subset = "Tcells") +
geom_hex(bins=200)+geom_gate("CD8_Tcells")+geom_gate("CD4_Tcells")+geom_gate("CD4-CD8-Tcells")+geom_gate("CD4+CD8+Tcells")


ggcyto(gs["NK cells_preT1D668-6F.fcs"],mapping = aes(x = "PE-Cy7-A",y = "Pacific Blue-A"),subset = "Tcells") +
geom_hex(bins=200)+geom_gate("CD8_Tcells")+geom_gate("CD4_Tcells")+geom_gate("CD4-CD8-Tcells")+geom_gate("CD4+CD8+Tcells")



template = rbind(
  template, add_pop(
    gs,alias ="CD4+CD27+CD28+",parent = "CD4_Tcells",dims = "APC-Cy7-A,FITC-A",gating_method = "boundary",gating_args ="min=c(1.6,1.2),max=c(3.5,3)"))



template = rbind(
  template, add_pop(
    gs,alias ="CD4+CD27-CD28+",parent = "CD4_Tcells",dims = "APC-Cy7-A,FITC-A",gating_method = "boundary",gating_args ="min=c(-0.2,1.2),max=c(1.6,3)"))



template = rbind(
  template, add_pop(
    gs,alias ="CD4+CD27-CD28-",parent = "CD4_Tcells",dims = "APC-Cy7-A,FITC-A",gating_method = "boundary",gating_args ="min=c(-0.2,0),max=c(1.6,1.2)"))


template = rbind(
  template, add_pop(
    gs,alias ="CD8+CD27+CD28+",parent = "CD8_Tcells",dims = "APC-Cy7-A,FITC-A",gating_method = "boundary",gating_args ="min=c(1.8,1.3),max=c(3.9,2.8)"))


ggcyto(gs[1:20],mapping = aes(x = "APC-Cy7-A",y = "FITC-A"),subset = "CD4_Tcells") +
geom_hex(bins=100) + geom_gate("CD8+CD27+CD28+")+ geom_gate("CD8+CD27-CD28+")


template = rbind(
  template, add_pop(
    gs,alias ="CD8+CD27+CD28+CD69+",parent = "CD8+CD27+CD28+",dims = "APC-A",gating_method = "mindensity",gating_args ="min=1.8,max=2.1"))

ggcyto(gs[1:20],mapping = aes(x = "APC-A",y = "SSC-A"),subset = "CD8+CD27+CD28+") +
geom_hex(bins=100) + geom_gate("CD8+CD27+CD28+CD69+")


template = rbind(
  template, add_pop(
    gs,alias ="CD8+CD27-CD28+",parent = "CD8_Tcells",dims = "APC-Cy7-A,FITC-A",gating_method = "boundary",gating_args ="min=c(-1,1.3),max=c(1.8,2.8)"))



template = rbind(
  template, add_pop(
    gs,alias ="CD8+CD27-CD28-",parent = "CD8_Tcells",dims = "APC-Cy7-A,FITC-A",gating_method = "boundary",gating_args ="min=c(-1,0),max=c(1.8,1.3)"))



template = rbind(
  template, add_pop(
    gs,alias ="CD8+CD27+CD28-",parent = "CD8_Tcells",dims = "APC-Cy7-A,FITC-A",gating_method = "boundary",gating_args ="min=c(1.8,0),max=c(3.9,1.3)"))


template = rbind(
  template, add_pop(
    gs,alias ="CD56+",pop="+",parent = "CD3-",dims = "PE-A",gating_method = "mindensity",gating_args = "min=1.1,max=1.3"))

 
template = rbind(
  template, add_pop(
    gs,alias ="CD56high",pop="+",parent = "CD56+",dims = "PE-A",gating_method = "mindensity",gating_args = "min=2.3,max=2.5"))


opt <- getOption("openCyto")
opt[["check.pop"]] <- FALSE
options(openCyto = opt)
template = rbind(
  template, add_pop(
    gs,alias ="CD3-CD56-",pop="-/-",parent = "CD3-", dims = "PerCP-A,PE-A",gating_method ="refGate",gating_args = "CD3-:CD56high"))



template = rbind(
  template, add_pop(
    gs,alias ="*",pop="-/+",parent = "CD3-CD56-",dims = "PE-A",gating_method = "mindensity",gating_args = "min=1.1,max=1.3"))

ggcyto(gs[1:15],mapping = aes(x = "PerCP-A",y = "PE-A"),subset = "CD3-CD56-") +
  geom_hex(bins=100)+geom_gate("PE-A+")+ geom_stats(type = c("count"))


opt <- getOption("openCyto")
opt[["check.pop"]] <- FALSE
options(openCyto = opt)
template = rbind(
  template, add_pop(
    gs,alias ="NK",pop="+",parent = "CD3-", dims = "PerCP-A,PE-A",gating_method ="refGate",gating_args = "PE-A+",collapseDataForGating = T))


template = rbind(
  template, add_pop(
    gs,alias ="NKact",pop="+",parent = "CD3-", dims = "PerCP-A,PE-A",gating_method ="refGate",gating_args = "CD56high",collapseDataForGating = T))


ggcyto(gs[16:40],mapping = aes(x = "PerCP-A",y = "PE-A"),subset = "CD3-") +
geom_hex(bins=100)+geom_gate("NK")+ geom_gate("NKact")+ geom_stats(type = c("percent"),digits = 1,label.padding = unit(0.01, "lines"),adjust=0.1,size=2.6)


template = rbind(
  template, add_pop(
    gs,alias="NK_CD69",parent = "NK",dims = "APC-A",gating_method = "mindensity",gating_args ="min=1.7,max=2.1"))

ggcyto(gs[100:130],mapping = aes(x = "APC-A",y = "SSC-A"),subset = "NK") +
geom_hex(bins=100) + geom_gate("NK_CD69")+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_stats(type = c("percent"),digits = 1,label.padding = unit(0.01, "lines"),adjust=0.8,size=2.6)

ggcyto(gs[1:30],mapping = aes(x = "APC-Cy7-A",y = "FITC-A"),subset = "CD4_Tcells") +
    geom_hex(bins=100) + geom_gate("CD4+CD27-CD28+")+ geom_gate("CD4+CD27+CD28+")+ geom_gate("CD4+CD27-CD28-")+ geom_gate("CD4+CD27+CD28-")+ggcyto_par_set(limits = list(x = c(-0.5, 4),y=c(0,3)))+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_stats(type = c("percent"),digits = 1,label.padding = unit(0.01, "lines"),adjust=0.8,size=2.6)
  
ggcyto(gs[1:30],mapping = aes(x = "APC-Cy7-A",y = "FITC-A"),subset = "CD8_Tcells") +
    geom_hex(bins=100) + geom_gate("CD8+CD27-CD28+")+ geom_gate("CD8+CD27+CD28+")+ geom_gate("CD8+CD27-CD28-")+ geom_gate("CD8+CD27+CD28-")+ggcyto_par_set(limits = list(x = c(-1, 4),y=c(0,3)))+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_stats(type = c("percent"),digits = 1,label.padding = unit(0.01, "lines"),adjust=0.8,size=2.6)



template = rbind(
  template, add_pop(
    gs,alias ="NKT",parent = "Tcells",dims ="PE-A",gating_method = "tailgate",gating_args ="adjust=0.9,tol=1e-2,min=1,max=1.3",collapseDataForGating = T))


ggcyto(gs[40:80],mapping = aes(x = "PerCP-A",y = "PE-A"),subset = "Tcells") +
  geom_hex(bins=100)+geom_gate("NKT")+ geom_stats(type = c("percent"),digits = 1,label.padding = unit(0.01, "lines"),adjust=0.1,size=2.6)


template = rbind(
  template, add_pop(
    gs,alias ="CD8+NKT",parent = "NKT",dims ="PE-Cy7-A,SSC-A",gating_method = "boundary",gating_args ="min=c(2.5,0),max=c(4.5,40000)",collapseDataForGating = T))

ggcyto(gs[100:120],mapping = aes(x = "PE-Cy7-A",y = "SSC-A"),subset = "NKT") +
  geom_hex(bins=100)+geom_gate("CD8+NKT")+ geom_stats(type = c("percent"),digits = 1,label.padding = unit(0.01, "lines"),adjust=0.1,size=2.6)


template=data.frame()
template = rbind(
  template, add_pop(
    gs,alias="NKact_CD69",parent = "NKact",dims = "APC-A",gating_method = "mindensity",gating_args ="min=0.8,max=1.1"))

ggcyto(gs[100:130],mapping = aes(x = "APC-A",y = "SSC-A"),subset = "NKact") +
geom_hex(bins=100) + geom_gate("NKact_CD69")+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_stats(type = c("percent"),digits = 1,label.padding = unit(0.01, "lines"),adjust=0.8,size=2.6)



nodesToHide <- c("Tcells_CD27+", "Tcells_CD28+", "CD3-CD56+", "CD3-CD56-", "PE-A+", "PE-A-","CD56high")
lapply(nodesToHide, function(thisNode) setNode(gs, thisNode, FALSE))


dir.create("NK_cells_statistics_last")
setwd("NK_cells_statistics_last/")


stats2<-getPopStats(gs,statistic="count")#extract stats
stats2[,prop := Count/ParentCount] #compute the cell proportions


md <- data.frame()[1:dim(stats2),]
md$file_name<-stats2$name
md$sample_id<-gsub("NKcells_*",'',md$file_name)
md$sample_id=gsub("_001.fcs",".fcs",md$sample_id)
md$sample_id<-gsub("*.fcs",'',md$sample_id)



metadata=read.delim("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/metadata_samples.txt",sep="\t",header = T)
md=merge(md,metadata,by="sample_id",all.x = T)

table(md$condition)
condition=as.data.frame(cbind(name=md$file_name,condition=md$condition))
condition=condition[!duplicated(condition$name),]
stats2 = merge(stats2,condition,by="name",all.x=T)
colnames(stats2)[7]="condition"
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

library(bestNormalize)
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

dir.create("outliers_last")
for (j in 1:length(nm)){
  p=names(nm)[j]
  q=nm[[j]]
  pdf(paste("outliers_last/",p,".pdf",sep=""))
  a=autoplot(gs[q],p,bins=200) + stats_null() + geom_stats(type = "percent",adjust = 0.8)
  plot(a)
  dev.off()
}

setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")
save_gs(gs,path="NK_supervised_gates_gs")
write.table(stats2,"NKcells_supervised.txt",quote=F,col.names = NA)

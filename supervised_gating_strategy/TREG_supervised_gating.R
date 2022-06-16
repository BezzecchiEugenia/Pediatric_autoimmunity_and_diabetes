setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")

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


gs=load_gs("TREG_gs/")
plot(gs)


template <- rbind(template, add_pop(gs,
                                    parent = "CD4_Tcells",
                                    alias = "CD25+CD127+",
                                    dims="APC-A,PE-Cy7-A",
                                    gating_method ="boundary",gating_args ="min=c(1.6,1.8),max=c(3,3)"))

template <- add_pop(gs,
                    parent = "CD4_Tcells",
                    alias = "CD127pos_CD45RAneg",
                    dims="PE-Cy7-A,PE-A",
                    gating_method ="boundary",gating_args = "min=c(1,0),max=c(3,1.5)")

template <- add_pop(gs,
                    parent = "CD4_Tcells",
                    alias = "CD127pos_CD45RApos",
                    dims="PE-Cy7-A,PE-A",
                    gating_method ="boundary",gating_args = "min=c(1,1.5),max=c(3,3.5)")


template <- add_pop(gs,
                    parent = "CD4_Tcells",
                    alias = "CD127neg_CD45RApos",
                    dims="PE-Cy7-A,PE-A",
                    gating_method ="boundary",gating_args = "min=c(-0.5,1.5),max=c(1,3.5)")

ggcyto(gs[1:20],mapping = aes(x = "PE-Cy7-A",y = "PE-A"),subset = "CD4_Tcells") +
  geom_hex(bins=100)+ geom_gate("CD127neg_CD45RApos") + geom_gate("CD127pos_CD45RAneg")+ geom_gate("CD127pos_CD45RApos")+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.9,size=2.7)


template <- rbind(template, add_pop(gs,
                                    parent = "CD4_Tcells",
                                    alias = "Foxp3low_CD45RA+",
                                    dims="FITC-A,PE-A",
                                    gating_method ="boundary",gating_args ="min=c(1.4,1.5),max=c(1.6,3.5)"))



template <- rbind(template, add_pop(gs,
                                    parent = "CD4_Tcells",
                                    alias = "Foxp3low_CD45RA-",
                                    dims="FITC-A,PE-A",
                                    gating_method ="boundary",gating_args ="min=c(1.45,-0.5),max=c(1.60,1.5)"))


template <- rbind(template, add_pop(gs,
                                    parent = "CD4_Tcells",
                                    alias = "Foxp3high_CD45RA-",
                                    dims="FITC-A,PE-A",
                                    gating_method ="boundary",gating_args ="min=c(1.60,-0.5),max=c(3,1.5)"))

ggcyto(gs[100:102],mapping = aes(x = "FITC-A",y = "PE-A"),subset = "CD4_Tcells") +
 geom_hex(bins = 80)+ggcyto_par_set(limits = list(x = c(-1, 5),y=c(-1,4)))+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_gate("Foxp3low_CD45RA+")+ geom_gate("Foxp3low_CD45RA-")+ geom_gate("Foxp3high_CD45RA-")+ geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.9,size=2.2)


ggcyto(gs[100:102],mapping = aes(x = "FITC-A",y = "SSC-A"),subset = "CD4_Tcells") +
 geom_hex(bins = 80)

ggcyto(gs[100:102],mapping = aes(x = "FITC-A",y = "PE-A"),subset = "CD4_Tcells") +
 geom_hex(bins = 80)+ggcyto_par_set(limits = list(x = c(-1, 5),y=c(-1,4)))+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_gate("Foxp3low_CD45RA-")


stats2<-getPopStats(gs,statistic="count")#extract stats
stats2[,prop := Count/ParentCount] #compute the cell proportions


dir.create("TREG_statistics_final")
setwd("TREG_statistics_final/")

md <- data.frame()[1:dim(stats2), ]
md$file_name<-stats2$name
md$sample_id<-gsub("TREG_*",'',md$file_name)
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

setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")
save_gs(gs,path="TREG_supervised_gates_gs")
write.table(stats2,"TREG_supervised.txt",quote=F,col.names = NA)



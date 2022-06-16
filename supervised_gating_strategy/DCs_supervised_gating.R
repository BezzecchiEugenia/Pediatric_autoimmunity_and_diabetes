
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
library(flowStats)
library(data.table)


gs=load_gs("DCs_gs")
plot(gs)

template=data.frame()


template <- rbind(template, add_pop(gs,
                                    parent = "mDCs",
                                    alias = "CD1c_mDCs",
                                    dims="Pacific Blue-A,PerCP-A",
                                    gating_method ="boundary",gating_args ="min=c(-0.2,1.5),max=c(1.5,3)"))

template <- rbind(template, add_pop(gs,
                                    parent = "mDCs",
                                    alias = "CD16_mDCs",
                                    dims="Pacific Blue-A,PerCP-A",
                                    gating_method ="boundary",gating_args ="min=c(2.3,-0.5),max=c(4.5,1.9)"))


template <- rbind(template, add_pop(gs,
                                    parent = "mDCs",
                                    alias = "CD1c_CD16_mDCs",
                                    dims="Pacific Blue-A,PerCP-A",
                                    gating_method ="boundary",gating_args ="min=c(1.5,1.9),max=c(3.5,3)"))


ggcyto(gs[100:110],mapping = aes(x = "Pacific Blue-A",y = "PerCP-A"),subset = "mDCs") +
    geom_hex(bins = 80) +geom_gate("CD1c_mDCs")+geom_gate("CD16_mDCs")+geom_gate("CD1c_CD16_mDCs")+ggcyto_par_set(limits = list(x = c(-1, 5),y=c(-1,4)))+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.01,size=3)
    

ggcyto(gs[200:240],mapping = aes(x = "Pacific Blue-A",y = "SSC-A"),subset = "HLADR+LIN-") +
    geom_hex(bins = 80) 

template <- rbind(template, add_pop(gs,
                                    parent = "HLADR+LIN-",
                                    alias = "classical_mono",
                                    dims="FITC-A,Pacific Blue-A",
                                    gating_method ="boundary",gating_args ="min=c(2,-0.5),max=c(3.4,2.5)"))

template <- rbind(template, add_pop(gs,
                                    parent = "HLADR+LIN-",
                                    alias = "intermediate_mono",
                                    dims="FITC-A,Pacific Blue-A",
                                    gating_method ="boundary",gating_args ="min=c(2.2,2.5),max=c(3.4,4)"))


template <- rbind(template, add_pop(gs,
                                    parent = "HLADR+LIN-",
                                    alias = "NC_mono",
                                    dims="FITC-A,Pacific Blue-A",
                                    gating_method ="boundary",gating_args ="min=c(1,2.5),max=c(2.2,4)"))


ggcyto(gs[1:30],mapping = aes(x = "FITC-A",y = "Pacific Blue-A"),subset = "HLADR+LIN-") +
  geom_hex(bins=200)+ geom_gate("classical_mono")+ geom_gate("intermediate_mono")+ geom_gate("NC_mono")+ggcyto_par_set(limits = list(x = c(0, 4),y=c(0,4.5)))+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.5,size=3)


ggcyto(gs[200:240],mapping = aes(x = "SSC-A",y = "Pacific Blue-A"),subset = "HLADR+LIN-") +
  geom_hex(bins=200)

dir.create("DCs_statistics")
setwd("DCs_statistics")

stats2<-getPopStats(gs,statistic="count")#extract stats
stats2[,prop := Count/ParentCount] #compute the cell proportions


md <- data.frame()[1:dim(stats2), ]
md$file_name<-stats2$name
md$sample_id<-gsub("DCs_*",'',md$file_name)
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


nm=mylist2[["HLADR+LIN-"]]
  
    n="DCs_T1D548.fcs"
    mygate = getGate(gs[[n]],"HLADR+LIN-") #extract a gates
    mygate@max["APC-Cy7-A"] = 4 #adjust the min value for CD8
    mygate@min["APC-Cy7-A"] = 2
    mygate@max["APC-A"] = 2.8 #adjust the min value for CD8
    mygate@min["APC-A"] = 0
    setGate(gs[[n]],"HLADR+LIN-",mygate) #set the gate for the cell population
    recompute(gs[[n]]) # recompute the statistics.
    
    
    n="DCs_preT1D1352.fcs"
    mygate = getGate(gs[[n]],"NC_mono") #extract a gates
    mygate@max["FITC-A"] = 2.2 #adjust the max value 
    mygate@min["FITC-A"] = 1
    mygate@max["Pacific Blue-A"] = 3 #adjust the min value
    mygate@min["Pacific Blue-A"] = 2
    setGate(gs[[n]],"NC_mono",mygate) #set the gate for the cell population
    recompute(gs[[n]]) # recompute the statistics.
    
    mygate = getGate(gs[[n]],"intermediate_mono") #extract a gates
    mygate@max["FITC-A"] = 3.3 #adjust the max value 
    mygate@min["FITC-A"] = 2.2
    mygate@max["Pacific Blue-A"] = 3 #adjust the min value 
    mygate@min["Pacific Blue-A"] = 2
    setGate(gs[[n]],"intermediate_mono",mygate) #set the gate for the cell population
    recompute(gs[[n]]) # recompute the statistics.
    
    
    mygate = getGate(gs[[n]],"classical_mono") #extract a gates
    mygate@max["FITC-A"] = 3.3 #adjust the max value
    mygate@min["FITC-A"] = 2
    mygate@max["Pacific Blue-A"] = 2 #adjust the min value 
    mygate@min["Pacific Blue-A"] = -1
    setGate(gs[[n]],"classical_mono",mygate) #set the gate for the cell population
    recompute(gs[[n]]) # recompute the statistics.
    
ggcyto(gs[1:40],mapping = aes(x = "Pacific Blue-A",y = "PerCP-A"),subset = "mDCs") +
    geom_hex(bins = 80) +geom_gate("CD1c_mDCs")+geom_gate("CD16_mDCs")+ggcyto_par_set(limits = list(x = c(-1, 5),y=c(-1,4)))+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.01,size=3)
 

setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")
save_gs(gs,path="DCs_supervised_gates_gs")
write.table(stats2,"DCs_supervised.txt",quote=F,col.names = NA)

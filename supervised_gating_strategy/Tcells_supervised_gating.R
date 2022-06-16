
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
library(data.table)


setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")
gs=load_gs("Tcells_gs")

dir.create("Tcells_norm")



template=data.frame()
template <- rbind(template, add_pop(gs,
                                    parent = "Tcells",
                                    alias = "CD4_Tcells",
                                    dims="APC-Cy7-A,APC-A",
                                    gating_method ="boundary",gating_args ="min=c(-1,2.4),max=c(2.5,4),quantile=0.95"))


template <- rbind(template, add_pop(gs,
                                    parent = "Tcells",
                                    alias = "CD8_Tcells",
                                    dims="APC-Cy7-A,APC-A",
                                    gating_method ="boundary",gating_args ="min=c(2.5,-1.25),max=c(4.4,2.4),quantile=0.95"))


gs <- normalize(gs, populations=c("CD4_Tcells", "CD8_Tcells"),
                dims=c( "APC-A","APC-Cy7-A"), minCountThreshold = 50)




template = rbind(
  template, add_pop(
    gs,pop="+",alias = "CD4_CD57+",parent = "CD4_Tcells",dims = "Pacific Blue-A",gating_method ="mindensity",gating_args = "min=0.5,max=1"
  )
)



pdf("Tcells_norm/CD4_CD57+.pdf")
ggcyto(gs[1:10],mapping = aes(x = "Pacific Blue-A",y = "SSC-A"),subset = "CD4_Tcells") +
  geom_hex(bins=200)+geom_gate("CD4_CD57+")
dev.off()

library(ggridges)
pdf("Tcells_norm/CD4_CD57+_BN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD4_CD57+"), aes(x="Pacific Blue-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()

gs <- normalize(gs, populations=c("CD4_CD57+"), bwFac = 3,
                dims=c("Pacific Blue-A"), minCountThreshold = 50)


pdf("Tcells_norm/CD4_CD57+_AN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD4_CD57+"), aes(x="Pacific Blue-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()


template = rbind(
  template, add_pop(
    gs,pop="-",alias = "CD8_CD57-",parent = "CD8_Tcells",dims = "Pacific Blue-A",gating_method ="mindensity",gating_args = "min=1.9,max=2"
  )
)

pdf("Tcells_norm/CD8_CD57-.pdf")
ggcyto(gs[1:10],mapping = aes(x = "Pacific Blue-A",y = "SSC-A"),subset = "CD8_Tcells") +
  geom_hex(bins=200)+geom_gate("CD8_CD57-")
dev.off()


pdf("Tcells_norm/CD8_CD57-_BN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD8_CD57-"), aes(x="Pacific Blue-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 5)+ facet_null()
dev.off()

gs <- normalize(gs, populations=c("CD8_CD57-"), bwFac =2,
                dims=c("Pacific Blue-A"), minCountThreshold = 50)


pdf("Tcells_norm/CD8_CD57-_AN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD8_CD57-"), aes(x="Pacific Blue-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()


template = rbind(
  template, add_pop(
    gs,pop="+",alias = "CD4_PD1+",parent = "CD4_Tcells",dims = "PE-Cy7-A",gating_method ="mindensity",gating_args = "min=1,max=1.5"
  )
)



pdf("Tcells_norm/CD4_PD1+.pdf")
ggcyto(gs[1:10],mapping = aes(x = "PE-Cy7-A",y = "SSC-A"),subset = "CD4_Tcells") +
  geom_hex(bins=200)+geom_gate("CD4_PD1+")
dev.off()


pdf("Tcells_norm/CD4_PD1+_BN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD4_PD1+"), aes(x="PE-Cy7-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()

gs <- normalize(gs, populations=c("CD4_PD1+"), bwFac = 3,
                dims=c("PE-Cy7-A"), minCountThreshold = 50)


pdf("Tcells_norm/CD4_PD1+_AN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD4_PD1+"), aes(x="PE-Cy7-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()


template = rbind(
  template, add_pop(
    gs,pop="+",alias = "CD8_PD1+",parent = "CD8_Tcells",dims = "PE-Cy7-A",gating_method ="mindensity",gating_args = "min=1,max=1.5"#gating_args = "max=0.5"
  )
)

pdf("Tcells_norm/CD8_PD1+.pdf")
ggcyto(gs[1:10],mapping = aes(x = "PE-Cy7-A",y = "SSC-A"),subset = "CD8_Tcells") +
  geom_hex(bins=200)+geom_gate("CD8_PD1+")
dev.off()


pdf("Tcells_norm/CD8_PD1+_BN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD8_PD1+"), aes(x="PE-Cy7-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()

gs <- normalize(gs, populations=c("CD8_PD1+"), bwFac = 3,
                dims=c("PE-Cy7-A"), minCountThreshold = 50)


pdf("Tcells_norm/CD8_PD1+_AN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD8_PD1+"), aes(x="PE-Cy7-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()


template = rbind(
  template, add_pop(
    gs,pop="+",alias = "CD4_CD45RA+",parent = "CD4_Tcells",dims = "FITC-A",gating_method ="mindensity",gating_args = "min=1,max=2"
  )
)


pdf("Tcells_norm/CD4_CD45RA+.pdf")
ggcyto(gs[1:10],mapping = aes(x = "FITC-A",y = "SSC-A"),subset = "CD4_Tcells") +
  geom_hex(bins=200)+geom_gate("CD4_CD45RA+")
dev.off()


pdf("Tcells_norm/CD4_CD45RA+_BN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD4_CD45RA+"), aes(x="FITC-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()

gs <- normalize(gs, populations=c("CD4_CD45RA+"), bwFac = 2,
                dims=c("FITC-A"), minCountThreshold = 50)


pdf("Tcells_norm/CD4_CD45RA+_AN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD4_CD45RA+"), aes(x="FITC-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()




template = rbind(
  template, add_pop(
    gs,pop="+",alias = "CD4_CCR7+",parent = "CD4_Tcells",dims = "PE-A",gating_method ="mindensity",gating_args = "min=1,max=2"
  )
)



pdf("Tcells_norm/CD4_CCR7+.pdf")
ggcyto(gs[1:10],mapping = aes(x = "PE-A",y = "SSC-A"),subset = "CD4_Tcells") +
  geom_hex(bins=200)+geom_gate("CD4_CCR7+")
dev.off()


pdf("Tcells_norm/CD4_CCR7+_BN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD4_CCR7+"), aes(x="PE-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()

gs <- normalize(gs, populations=c("CD4_CCR7+"), bwFac = 2,
                dims=c("PE-A"), minCountThreshold = 50)


pdf("Tcells_norm/CD4_CCR7+_AN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD4_CCR7+"), aes(x="PE-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()



template = rbind(
  template, add_pop(
    gs,pop="+",alias = "CD8_CD45RA+",parent = "CD8_Tcells",dims = "FITC-A",gating_method ="mindensity",gating_args = "min=2.4,max=3"
  )
)


pdf("Tcells_norm/CD8_CD45RA+.pdf")
ggcyto(gs[1:10],mapping = aes(x = "FITC-A",y = "SSC-A"),subset = "CD8_Tcells") +
  geom_hex(bins=200)+geom_gate("CD8_CD45RA+")
dev.off()


pdf("Tcells_norm/CD8_CD45RA+_BN.pdf")
ggcyto(gs_pop_get_data(gs[1:30], "CD8_CD45RA+"), aes(x="FITC-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()

gs <- normalize(gs, populations=c("CD8_CD45RA+"), bwFac = 1,
                dims=c("FITC-A"), minCountThreshold = 50)


pdf("Tcells_norm/CD8_CD45RA+_AN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD8_CD45RA+"), aes(x="FITC-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()


template = rbind(
  template, add_pop(
    gs,pop="+",alias = "CD8_CCR7+",parent = "CD8_Tcells",dims = "PE-A",gating_method ="mindensity"
  )
)

pdf("Tcells_norm/CD8_CCR7+.pdf")
ggcyto(gs[1:10],mapping = aes(x = "PE-A",y = "SSC-A"),subset = "CD8_Tcells") +
  geom_hex(bins=200)+geom_gate("CD8_CCR7+")
dev.off()

pdf("Tcells_norm/CD8_CCR7+_BN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD8_CCR7+"), aes(x="PE-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()

gs <- normalize(gs, populations=c("CD8_CCR7+"), bwFac = 2,
                dims=c("PE-A"), minCountThreshold = 50)


pdf("Tcells_norm/CD8_CCR7+_AN.pdf")
ggcyto(gs_pop_get_data(gs[1:10], "CD8_CCR7+"), aes(x="PE-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()


template <- rbind(template, add_pop(gs,
                                    parent = "Tcells",
                                    alias = "CD4_CD8_Tcells",
                                    dims="APC-Cy7-A,APC-A",
                                    gating_method ="boundary",gating_args ="min=c(2.5,2.4),max=c(4.4,4),quantile=0.95"))


ggcyto(gs[1:10],mapping = aes(x = "APC-Cy7-A",y = "APC-A"),subset = "Tcells") +
  geom_hex(bins=200)+geom_gate("CD4_Tcells")+geom_gate("CD8_Tcells")+geom_gate("CD4_CD8_Tcells")



template = rbind(
  template, add_pop(
    gs,alias ="Naive_Th",parent = "CD4_Tcells",dims = "FITC-A,PE-A",gating_method = "boundary",gating_args ="min=c(1,1.5),max=c(3.5,3)"))


template = rbind(
  template, add_pop(
    gs,alias ="Eff_Th",parent = "CD4_Tcells",dims = "FITC-A,PE-A",gating_method = "boundary",gating_args ="min=c(-0.5,-0.5),max=c(1,1.5)"))

template = rbind(
  template, add_pop(
    gs,alias ="Temra_Th",parent = "CD4_Tcells",dims = "PE-A,FITC-A",gating_method = "boundary",gating_args ="min=c(-0.5,1),max=c(1.5,3.5)"))


ggcyto(gs[1:10],mapping = aes(x = "PE-A",y = "FITC-A"),subset = "CD4_Tcells") +
  geom_hex(bins=200)+geom_gate("Eff_Th")+geom_gate("Naive_Th")+geom_gate("Temra_Th")


template = rbind(
  template, add_pop(
    gs,alias ="CM_Th",parent = "CD4_Tcells",dims = "PE-A,FITC-A",gating_method = "boundary",gating_args ="min=c(1.5,-0.5),max=c(3,1)"))

ggcyto(gs[1:10],mapping = aes(x = "PE-A",y = "FITC-A"),subset = "CD4_Tcells") +
  geom_hex(bins=200)+geom_gate("Eff_Th")+geom_gate("Naive_Th")+geom_gate("Temra_Th")+geom_gate("CM_Th")


template = rbind(
  template, gs_add_gating_method(
    gs,alias ="Temra_CD57_Th",parent = "Temra_Th",dims = "Pacific Blue-A",gating_method = "mindensity",gating_args ="min=1.9,max=2.1"))


ggcyto(gs[1:10],mapping = aes(x = "Pacific Blue-A",y = "SSC-A"),subset = "Temra_Th") +
  geom_hex(bins=200)+geom_gate("Temra_CD57_Th")


template = rbind(
  template, add_pop(
    gs,alias ="Eff_CD57_Th",parent = "Eff_Th",dims = "Pacific Blue-A",gating_method = "mindensity",gating_args ="min=1.9,max=2.1"))

ggcyto(gs[1:10],mapping = aes(x = "Pacific Blue-A",y = "SSC-A"),subset = "Eff_Th") +
  geom_hex(bins=200)+geom_gate("Eff_CD57_Th")

template <- rbind(
  template, add_pop(gs,pop="+",
                    parent = "Eff_Th",
                    alias = "CD4_EM_PE-Cy7-A+",
                    dims="PE-Cy7-A",
                    gating_method ="mindensity",gating_args ="min=-1,max=0.5",collapseDataForGating = T))


template <- rbind(
  template, add_pop(gs,pop="+",
                    parent = "Eff_Th",
                    alias = "EM_PD1+_Th",
                    dims="PE-Cy7-A",
                    gating_method ="mindensity",gating_args ="min=1.5,max=2",collapseDataForGating = T))

ggcyto(gs[1:10],mapping = aes(x = "PE-Cy7-A",y = "SSC-A"),subset = "Eff_Th") +
  geom_hex(bins=200)+geom_gate("EM_PD1+_Th")


opt <- getOption("openCyto")
opt[["check.pop"]] <- FALSE
options(openCyto = opt)

template = rbind(
  template, add_pop(
    gs,alias ="PD1_CD57_EM_CD4",pop="+/+",parent = "Eff_Th", dims = "PE-Cy7-A,Pacific Blue-A",gating_method ="refGate",gating_args = "EM_PD1+_Th:Eff_CD57_Th"))

ggcyto(gs[1:10],mapping = aes(x = "PE-Cy7-A",y = "Pacific Blue-A"),subset = "Eff_Th") +
  geom_hex(bins=200)+geom_gate("PD1_CD57_EM_CD4")


template = rbind(
  template, add_pop(
    gs,alias ="Naive_T_CD8",parent = "CD8_Tcells",dims = "FITC-A,PE-A",gating_method = "boundary",gating_args ="min=c(1.8,1.3),max=c(3.5,2.9)"))

template = rbind(
  template, add_pop(
    gs,alias ="Temra_CD8",parent = "CD8_Tcells",dims = "FITC-A,PE-A",gating_method = "boundary",gating_args ="min=c(1.8,-0.5),max=c(3.5,1.3)"))


template = rbind(
  template, add_pop(
    gs,alias ="Eff_CD8",parent = "CD8_Tcells",dims = "FITC-A,PE-A",gating_method = "boundary",gating_args ="min=c(-0.8,-0.5),max=c(1.8,1.3)"))

template = rbind(
  template, add_pop(
    gs,alias ="Temra_CD57_CD8",parent = "Temra_CD8",dims = "Pacific Blue-A",gating_method = "mindensity",gating_args ="min=1.9,max=2.1"))

ggcyto(gs[1:10],mapping = aes(x = "Pacific Blue-A",y = "SSC-A"),subset = "Temra_CD8") +
  geom_hex(bins=200)+geom_gate("Temra_CD57_CD8")

template = rbind(
  template, add_pop(
    gs,alias ="Eff_CD57_CD8",parent = "Eff_CD8",dims = "Pacific Blue-A",gating_method = "mindensity",gating_args ="min=1.9,max=2.1"))

ggcyto(gs[1:10],mapping = aes(x = "Pacific Blue-A",y = "SSC-A"),subset = "Eff_CD8") +
  geom_hex(bins=200)+geom_gate("Eff_CD57_CD8")

template <- rbind(
  template, add_pop(gs,pop="+",
                    parent = "Eff_CD8",
                    alias = "EM_PD1+_CD8",
                    dims="PE-Cy7-A",
                    gating_method ="gate_mindensity2",gating_args ="min=1.4,max=2"))

ggcyto(gs[1:10],mapping = aes(x = "PE-Cy7-A",y = "SSC-A"),subset = "Eff_CD8") +
  geom_hex(bins=70)+geom_gate("EM_PD1+_CD8")


opt <- getOption("openCyto")
opt[["check.pop"]] <- FALSE
options(openCyto = opt)
template = rbind(
  template, add_pop(
    gs,alias ="PD1_CD57_EM_CD8",pop="+/+",parent = "Eff_CD8", dims = "PE-Cy7-A,Pacific Blue-A",gating_method ="refGate",gating_args = "EM_PD1+_CD8:Eff_CD57_CD8"))

ggcyto(gs[1:10],mapping = aes(x = "PE-Cy7-A",y = "Pacific Blue-A"),subset = "Eff_CD8") +
  geom_hex(bins=200)+geom_gate("PD1_CD57_EM_CD8")


nodesToHide <- c( "CD4_CD45RA+", "CD4_CCR7+","CD4_PD1+", "CD4_CD57+","CD8_PD1+", "CD8_CD57-",
                  "CD8_CD45RA+", "CD8_CCR7+","CD4_EM_PE-Cy7-A+")
lapply(nodesToHide, function(thisNode) setNode(gs, thisNode, FALSE))

stats2<-getPopStats(gs,statistic="count")#extract stats
stats2[,prop := Count/ParentCount] #compute the cell proportions

dir.create("Tcells_statistics_final")
setwd("Tcells_statistics_final")

md <- data.frame()[1:dim(stats2), ]
md$file_name<-stats2$name
md$sample_id<-gsub("Tcells_*",'',md$file_name)
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
save_gs(gs,path="Tcells_supervised_gates_gs")
write.table(stats2,"Tcells_supervised.txt",quote=F,col.names = NA)


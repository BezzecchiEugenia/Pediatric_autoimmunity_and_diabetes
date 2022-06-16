
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
library(data.table)


path<-"/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/myFiles/"
files <- list.files(pattern = "Bcells",path, full =F,recursive = F)
files=files[1:15]


fs<-read.flowSet(files,path,transformation = F,truncate_max_range = F)
gs <- GatingSet(fs)


Comp=read.table("../spillover_autospill/autospill-spillover-B2.csv",sep=",",header=T,row.names = 1)
colnames(Comp)=row.names(Comp)


for (i in 1:length(gs)){
  gs[[i]]<-cyto_compensate(gs[[i]],Comp)
}

a=gs_get_compensations(gs)
gs<-cyto_transform(gs,type = "logicle")


setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data/")


polygonGate=function(fr=NA, pp_res=NULL, channels=NA, ...){
  flowCore:::polygonGate(fr, stains=channels, ...)}
registerPlugins(fun=polygonGate,methodName="polygonGate",dep="flowCore","gating")

template <- add_pop(gs,
                    parent = "root",
                    alias = "Singlets",
                    dims="FSC-A,FSC-H",
                    gating_method ="polygonGate",gating_args ="boundaries=structure(c(0,220000,220000,20000,0,125000,170000,50000),.Dim = c(4L, 2L), .Dimnames = list(NULL, c('FSC-A', 'FSC-H')))")


template <- rbind(
  template,add_pop(gs,
                   parent = "Singlets",
                   alias = "white cells",
                   dims="BV510-A",
                   gating_method ="mindensity",gating_args ="min=1.7,max=2"))


template = rbind(
  template, add_pop(
    gs,alias = "Lymphocytes",parent = "white cells",dims = "FSC-A,SSC-A",gating_method =
      "flowClust.2d",preprocessing_method = "prior_flowClust",preprocessing_args="K=3",collapseDataForGating = FALSE,gating_args =
      "K=3, target=c(1e5,2.5e4), quantile=0.98,min=c(50000,0),max=c(1.5e5,3e4)"
  )
)


nm=sampleNames(gs)

mat<- matrix(c(78346624,24069920,9969920,25426784), nrow = 2, dimnames = list(c("FSC-A","SSC-A"), c("FSC-A","SSC-A")))

for(n in nm){
  mygate = getGate(gs[[n]],"Lymphocytes") #extract a gates
  mygate@cov=mat
  setGate(gs[[n]],"Lymphocytes",mygate) #set the gate for the cell population
  recompute(gs[[n]]) # recompute the statistics.
}


template = rbind(
  template, add_pop(
    gs,alias = "Bcells",parent = "Lymphocytes",dims = "Pacific Blue-A",gating_method = "mindensity"
  )
)

template = rbind(
  template, add_pop(
    gs,pop="-",alias = "Bcells_CD27-",parent = "Bcells",dims = "APC-Cy7-A",gating_method ="mindensity",gating_args = "min=1.8,max=2"
  )
)

dir.create("Bcells_norm")

pdf("Bcells_norm/Bcells_CD27-.pdf",25,20)
ggcyto(gs[1:40],mapping = aes(x = "APC-Cy7-A",y = "SSC-A"),subset = "Bcells") +
  geom_hex(bins=200)+geom_gate("Bcells_CD27-")
dev.off()

library(ggridges)
pdf("Bcells_norm/Bcells_CD27-_BN.pdf")
ggcyto(gs_pop_get_data(gs[1:40], "Bcells_CD27-"), aes(x="APC-Cy7-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()

gs <- normalize(gs, populations=c("Bcells_CD27-"), bwFac = 3,
                dims=c("APC-Cy7-A"), minCountThreshold = 50)


pdf("Bcells_norm/Bcells_CD27-_AN.pdf")
ggcyto(gs_pop_get_data(gs[1:40], "Bcells_CD27-"), aes(x="APC-Cy7-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()

pdf("Bcells_norm/Bcells_CD27-_AN_2D.pdf",25,30)
ggcyto(gs[1:40],mapping = aes(x = "APC-Cy7-A",y = "SSC-A"),subset = "Bcells") +
  geom_hex(bins=200)+geom_gate("Bcells_CD27-")
dev.off()

template = rbind(
  template, add_pop(
    gs,pop="+",alias = "Bcells_IgD+",parent = "Bcells",dims = "FITC-A",gating_method ="mindensity",gating_args = "min=1,max=2"
  )
)

pdf("Bcells_norm/Bcells_IgD+.pdf")
ggcyto(gs[1:40],mapping = aes(x = "FITC-A",y = "SSC-A"),subset = "Bcells") +
  geom_hex(bins=200)+geom_gate("Bcells_IgD+")
dev.off()

library(ggridges)
pdf("Bcells_norm/Bcells_IgD+_BN.pdf")
ggcyto(gs_pop_get_data(gs[1:40], "Bcells_IgD+"), aes(x="FITC-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()

gs <- normalize(gs, populations=c("Bcells_IgD+"), bwFac = 1,
                dims=c("FITC-A"), minCountThreshold = 50)


pdf("Bcells_norm/Bcells_Bcells_IgD+_AN.pdf")
ggcyto(gs_pop_get_data(gs[1:40], "Bcells_IgD+"), aes(x="FITC-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()

pdf("Bcells_norm/Bcells_IgD+_AN_2D.pdf")
ggcyto(gs[1:40],mapping = aes(x = "FITC-A",y = "SSC-A"),subset = "Bcells") +
  geom_hex(bins=200)+geom_gate("Bcells_IgD+")
dev.off()


pdf("Bcells_norm/CD38.pdf")
ggcyto(gs[1:40],mapping = aes(x = "PerCP-A",y = "SSC-A"),subset = "Bcells") +
  geom_hex(bins=100)
dev.off()


pdf("Bcells_norm/CD21.pdf")
ggcyto(gs[1:40],mapping = aes(x = "Pe-Cy7-A",y = "SSC-A"),subset = "Bcells") +
  geom_hex(bins=100)
dev.off()


template <- rbind(template, add_pop(gs,
                                    parent = "Bcells",
                                    alias = "Naive_B",
                                    dims="APC-Cy7-A,FITC-A",
                                    gating_method ="boundary",gating_args ="min=c(-0.5,1.6),max=c(1.9,3.4)"))


template = rbind(
  template, add_pop(
    gs,pop="+",alias = "Naive_IgM+",parent = "Naive_B",dims = "PE-A",gating_method ="mindensity",gating_args = "min=0,max=1"
  )
)

template = rbind(
  template, add_pop(
    gs,pop="+",alias = "Naive_IgM+",parent = "Naive_B",dims = "PE-A",gating_method ="mindensity",gating_args = "min=0,max=1"
  )
)


ggcyto(gs[100:130],mapping = aes(x = "PE-A",y = "SSC-A"),subset = "Naive_B") +
  geom_hex(bins=100) + geom_gate("IGD_only_naive")+ geom_gate("Naive_Bcells")+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.9,size=2.7)


pdf("Bcells_norm/Naive_IgM+.pdf")
ggcyto(gs[1:40],mapping = aes(x = "PE-A",y = "SSC-A"),subset = "Naive_B") +
  geom_hex(bins=100)+geom_gate("Naive_IgM+")
dev.off()

library(ggridges)
pdf("Bcells_norm/Naive_IgM+_BN.pdf")
ggcyto(gs_pop_get_data(gs[1:40], "Naive_IgM+"), aes(x="PE-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()

gs <- normalize(gs, populations=c("Naive_IgM+"), bwFac = 3,
                dims=c("PE-A"), minCountThreshold = 50)


pdf("Bcells_norm/Naive_IgM+_AN.pdf")
ggcyto(gs_pop_get_data(gs[1:40], "Naive_B"), aes(x="PE-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()
dev.off()


md <- data.frame(sampleNames(gs))
names(md)[1]="file_name"
md$sample_id<-gsub("Bcells_*",'',md$file_name)
md$sample_id=gsub("_001.fcs",".fcs",md$sample_id)
md$sample_id=gsub("_2.fcs",".fcs",md$sample_id)
md$sample_id<-gsub("*.fcs",'',md$sample_id)

metadata=read.delim("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/metadata_samples.txt",sep="\t",header = T)
md=merge(md,metadata,by="sample_id",all.x=T)
table(duplicated(md$sample_id))


save_gs(gs,path="Bcells_gs")


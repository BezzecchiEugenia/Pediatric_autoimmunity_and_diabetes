
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


path<-"/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/myFiles/"
files <- list.files(pattern = "NKcells",path, full =F,recursive = F)
files=files[1:15]

fs<-read.flowSet(files,path,transformation = F,truncate_max_range = F)
gs <- GatingSet(fs)
rm(fs)

Comp=read.table("../spillover_autospill/autospill-spillover-NK2.csv",sep=",",header=T,row.names = 1)
colnames(Comp)=row.names(Comp)


for (i in 1:length(gs)){
  gs[[i]]<-cyto_compensate(gs[[i]],Comp)
}

a=gs_get_compensations(gs)
gs<-cyto_transform(gs,type = "logicle")


polygonGate=function(fr=NA, pp_res=NULL, channels=NA, ...){
  flowCore:::polygonGate(fr, stains=channels, ...)}
registerPlugins(fun=polygonGate,methodName="polygonGate",dep="flowCore","gating")

template <- add_pop(gs,
                    parent = "root",
                    alias = "Singlets",
                    dims="FSC-A,FSC-H",
                    gating_method ="polygonGate",gating_args ="boundaries=structure(c(0,260000,260000,20000,0,150000,200000,50000),.Dim = c(4L, 2L), .Dimnames = list(NULL, c('FSC-A', 'FSC-H')))")




template = rbind(
  template, add_pop(gs,
                    parent = "Singlets",
                    alias = "white cells",
                    dims="BV510-A",
                    gating_method ="mindensity",gating_args ="min=1.7,max=2"))


template = rbind(
  template, add_pop(
    gs,alias = "Lymphocytes",parent = "white cells",dims = "FSC-A,SSC-A",gating_method =
      "flowClust.2d",preprocessing_method = "prior_flowClust",preprocessing_args="K=3",collapseDataForGating = FALSE,gating_args =
      "K=3, target=c(1e5,2.5e4), quantile=0.98,min=c(50000,0),max=c(1.5e5,3e4)"
    #"K=5, quantile=0.9"
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
    gs,alias = "Tcells",parent = "Lymphocytes",dims = "PerCP-A",gating_method = "mindensity"
  )
)

template = rbind(
  template, add_pop(
    gs,pop= "-",alias ="CD3-",parent = "Lymphocytes",dims = "PerCP-A",gating_method = "mindensity"))


template <- rbind(template, add_pop(gs,
                                    parent = "Tcells",
                                    alias = "CD4_Tcells",
                                    dims="PE-Cy7-A,Pacific Blue-A",
                                    gating_method ="boundary",gating_args ="min=c(-0.7,1.8),max=c(2.2,3.2),quantile=0.95"))


template <- rbind(template, add_pop(gs,
                                    parent = "Tcells",
                                    alias = "CD8_Tcells",
                                    dims="PE-Cy7-A,Pacific Blue-A",
                                    gating_method ="boundary",gating_args ="min=c(2.2,-0.5),max=c(4.4,1.8),quantile=0.95"))


gs <- normalize(gs, populations=c("CD8_Tcells"),bwFac=1,
                dims=c("PE-Cy7-A","Pacific Blue-A"), minCountThreshold = 50)

gs <- normalize(gs, populations=c("CD4_Tcells"),bwFac=1,
                dims=c("Pacific Blue-A","PE-Cy7-A"), minCountThreshold = 50)


template = rbind(
  template, add_pop(
    gs,alias = "Tcells_CD27+",parent = "Tcells",dims = "APC-Cy7-A",gating_method ="mindensity",gating_args = "min=1.8,max=2.1"
   
  )
)

gs <- normalize(gs, populations=c("Tcells_CD27+"),bwFac = 3,
                dims=c("APC-Cy7-A"), minCountThreshold = 0)

template = rbind(
  template, add_pop(
    gs,alias = "Tcells_CD28+",parent = "Tcells",dims = "FITC-A",gating_method ="mindensity",gating_args = "min=1.4,max=1.6"

  )
)


gs <- normalize(gs, populations=c("Tcells_CD28+"),bwFac = 3,
                dims=c("FITC-A"), minCountThreshold = 0)



template = rbind(
  template, add_pop(
    gs,alias = "CD3-CD56+",pop='+',parent = "CD3-",dims = "PE-A",gating_method ="mindensity",gating_args = "min=1,max=1.2"
  )
)

gs <- normalize(gs, populations=c("CD3-CD56+"),bwFac =1,
                dims=c("PE-A"), minCountThreshold = 0)


save_gs(gs,path="NK_gs")


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

path<-"/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/myFiles/"
files <- list.files(pattern = "Tcells",path, full =F,recursive = F)
files=files[1:15]

fs<-read.flowSet(files,path,transformation = F,truncate_max_range = F)

gs <- GatingSet(fs)

Comp=read.table("../spillover_autospill/autospill-spillover-T2.csv",sep=",",header=T,row.names = 1)
colnames(Comp)=row.names(Comp)


for (i in 1:length(gs)){
  gs[[i]]<-cyto_compensate(gs[[i]],Comp)
}

a=gs_get_compensations(gs)
gs<-cyto_transform(gs,type = "logicle")


samples=as.character(c("T cells_CD014.fcs","T cells_CD019.fcs","T cells_CD029.fcs","T cells_preT1D309-5F.fcs","T cells_T1D493.fcs","T cells_T1D550.fcs"))
q=match(samples,sampleNames(gs))
gs=gs[-q]



polygonGate=function(fr=NA, pp_res=NULL, channels=NA, ...){
  flowCore:::polygonGate(fr, stains=channels, ...)}
registerPlugins(fun=polygonGate,methodName="polygonGate",dep="flowCore","gating")

template=data.frame()

template <- add_pop(gs,
                    parent = "root",
                    alias = "Singlets",
                    dims="FSC-A,FSC-H",
                    gating_method ="polygonGate",gating_args ="boundaries=structure(c(0,260000,260000,20000,0,150000,200000,50000),.Dim = c(4L, 2L), .Dimnames = list(NULL, c('FSC-A', 'FSC-H')))")

template <- rbind(
  template, add_pop(gs,
                    parent = "Singlets",
                    alias = "white cells",
                    dims="BV510-A",
                    gating_method ="mindensity",gating_args ="min=1.7,max=2"))


template = rbind(
  template, add_pop(
    gs,alias = "Lymphocytes",parent = "white cells",dims = "FSC-A,SSC-A",gating_method =
      "flowClust.2d",preprocessing_method = "prior_flowClust",preprocessing_args="K=3",collapseDataForGating = FALSE,gating_args =
      "K=3, target=c(1e5,2.5e4), quantile=0.97,min=c(50000,0),max=c(1.5e5,3e4)"
  )
)


nm=sampleNames(gs)

mat<- matrix(c(78346624,24069920,9969920,25426784), nrow = 2, dimnames = list(c("FSC-A","SSC-A"), c("FSC-A","SSC-A")))

for(n in nm){
  mygate = getGate(gs[[n]],"Lymphocytes") #extract a gates
  setGate(gs[[n]],"Lymphocytes",mygate) #set the gate for the cell population
  recompute(gs[[n]]) # recompute the statistics.
}


template = rbind(
  template, add_pop(
    gs,alias = "Tcells",parent = "Lymphocytes",dims = "PerCP-A",gating_method = "mindensity",gating_args ="max=3"
  )
)

save_gs(gs,path="Tcells_gs")


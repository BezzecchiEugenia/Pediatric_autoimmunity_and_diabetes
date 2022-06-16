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


path<-"/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/myFiles/"
files <- list.files(pattern = "TREG",path, full =F,recursive = F)
files=files[1:15]

fs<-read.flowSet(files,path,transformation = F,truncate_max_range = F)


gs <- GatingSet(fs)
rm(fs)

Comp=read.table("../spillover_autospill/autospill-spillover-TREG2.csv",sep=",",header=T,row.names = 1)
colnames(Comp)=row.names(Comp)


for (i in 1:length(gs)){
  gs[[i]]<-cyto_compensate(gs[[i]],Comp)
}

#####STORE COMPENSATION
a=gs_get_compensations(gs)
gs<-cyto_transform(gs,type = "logicle")


polygonGate=function(fr=NA, pp_res=NULL, channels=NA, ...){
  flowCore:::polygonGate(fr, stains=channels, ...)}
registerPlugins(fun=polygonGate,methodName="polygonGate",dep="flowCore","gating")

template <- add_pop(gs,
                    parent = "root",
                    alias = "Singlets",
                    dims="FSC-A,FSC-H",
                    gating_method ="polygonGate",gating_args ="boundaries=structure(c(0,260000,260000,20000,0,165000,220000,40000),.Dim = c(4L, 2L), .Dimnames = list(NULL, c('FSC-A', 'FSC-H')))")


template = rbind(
  template, add_pop(gs,
                    parent = "Singlets",
                    alias = "white cells",
                    dims="BV510-A",
                    gating_method ="mindensity",gating_args ="min=1.7,max=2"))


template = rbind(
  template, add_pop(
    gs,alias = "Lymphocytes",parent = "white cells",dims = "FSC-A,SSC-A",gating_method =
      "flowClust.2d",preprocessing_method = "prior_flowClust",preprocessing_args="K=2",collapseDataForGating = FALSE,gating_args =
      "K=2, target=c(1e5,2.5e4), quantile=0.95,max=c(1.8e5,1.5e5)"
  )
)


template = rbind(
  template, add_pop(
    gs,alias ="CD4_Tcells",parent = "Lymphocytes",dims = "PerCP-A,Pacific Blue-A",gating_method = "boundary",gating_args ="min=c(2,2),max=c(3.5,3.2)"))


template = rbind(
  template, add_pop(
    gs,alias = "CD25+",parent = "CD4_Tcells",dims = "APC-A",gating_method ="mindensity",gating_args = "max=1.5"

  )
)


gs <- normalize(gs, populations=c("CD25+"),bwFac = 3,
                dims=c("APC-A"), minCountThreshold = 0)




template = rbind(
  template, add_pop(
    gs,pop="+",alias = "FoxP3+",parent = "CD4_Tcells",dims = "FITC-A",gating_method ="mindensity",gating_args = "min=0.5,max=1.5"

  )
)


gs <- normalize(gs, populations=c("FoxP3+"), bwFac = 3,
                dims=c("FITC-A"), minCountThreshold = 0)



template = rbind(
  template, add_pop(
    gs,alias = "CD127+",parent = "CD4_Tcells",dims = "PE-Cy7-A",gating_method ="mindensity",gating_args = "max=1"

  )
)


gs <- normalize(gs, populations=c("CD127+"), bwFac = 3,
                dims=c("PE-Cy7-A"), minCountThreshold = 0)


template = rbind(
  template, add_pop(
    gs,alias = "CD45RA+",parent = "CD4_Tcells",dims = "PE-A",gating_method ="mindensity"

  )
)

gs <- normalize(gs, populations=c("CD45RA+"), bwFac = 3, 
                dims=c("PE-A"), minCountThreshold = 0)

ggcyto(gs_pop_get_data(gs, "CD45RA+"), aes(x="PE-A")) +
  geom_density_ridges(mapping = aes(y = name), alpha = 0.4)+ facet_null()

ggcyto(gs,mapping = aes(x = "APC-A",y = "SSC-A"),subset = "CD4_Tcells") +
  geom_hex(bins=200)+geom_gate("CD25+")

save_gs(gs,path="TREG_gs")



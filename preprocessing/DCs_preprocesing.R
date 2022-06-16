
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
files <- list.files(pattern = "DCs",path, full =F,recursive = F)
files=files[1:15]


fs<-read.flowSet(files,path,transformation = F,truncate_max_range = F)
gs <- GatingSet(fs)


Comp=read.table("../spillover_autospill/autospill-spillover-DC.csv",sep=",",header=T,row.names = 1)
colnames(Comp)=row.names(Comp)


for (i in 1:length(gs)){
  gs[[i]]<-cyto_compensate(gs[[i]],Comp)
}

a=gs_get_compensations(gs)
gs<-cyto_transform(gs,type = "logicle")

olygonGate=function(fr=NA, pp_res=NULL, channels=NA, ...){
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
    gs,alias = "Mono_Lymphocytes",parent = "white cells",dims = "FSC-A,SSC-A",
    gating_method ="polygonGate",gating_args ="boundaries=structure(c(40000,260000,260000,40000,0,0,110000,40000),.Dim = c(4L, 2L), .Dimnames = list(NULL, c('FSC-A', 'SSC-A')))"))


template = rbind(
  template, add_pop(
    gs,alias ="HLADR+LIN-",parent = "Mono_Lymphocytes",dims = "APC-A,APC-Cy7-A",gating_method = "boundary",gating_args ="min=c(-0.5,1.7),max=c(1.9,4)"))

gs=normalize(gs
             , populations = "HLADR+LIN-"
             , dims = c("APC-A","APC-Cy7-A")
             , minCountThreshold = 100
             , nPeaks = list('HLADR+LIN-' = 1)
             , chunksize = 10
             , bwFac = 1
)

ggcyto(gs[255:285],mapping = aes(x = "APC-A",y = "APC-Cy7-A"),subset = "Mono_Lymphocytes") +
geom_hex(bins=200)+ geom_gate("HLADR+LIN-")+ggcyto_par_set(limits = list(x = c(0, 4),y=c(0,4.5)))+axis_x_inverse_trans()+axis_y_inverse_trans()+ geom_stats(type = c("percent"),digits = 2,label.padding = unit(0.01, "lines"),adjust=0.5,size=3)

template = rbind(
  template, add_pop(
    gs,pop="-",alias ="CD14-",parent = "HLADR+LIN-",dims = "FITC-A",gating_method = "mindensity",gating_args ="min=2.2,max=2.5"))



template = rbind(
  template, add_pop(
    gs,alias ="CD14-_PE-A+",parent = "CD14-",dims = "PE-A",gating_method = "mindensity",gating_args ="min=1.8,max=2.1"))



gs <- normalize(gs, populations=c("CD14-_PE-A+"),bwFac = 2,
                dims=c("PE-A"), minCountThreshold = 0)


template <- rbind(template, add_pop(gs,
                                    parent = "CD14-",
                                    alias = "pDCs",
                                    dims="PE-A,PE-Cy7-A",
                                    gating_method ="boundary",gating_args ="min=c(0,2.7),max=c(1.8,4)"))


template <- rbind(template, add_pop(gs,
                                    parent = "CD14-",
                                    alias = "mDCs",
                                    dims="PE-A,PE-Cy7-A",
                                    gating_method ="boundary",gating_args ="min=c(1.8,1),max=c(3.5,3.2)"))

template = rbind(
  template, add_pop(
    gs,alias ="HLADR+LIN-CD14+",parent = "HLADR+LIN-",dims = "FITC-A",gating_method = "mindensity",gating_args ="min=1.8,max=2.1"))


gs <- normalize(gs, populations=c("HLADR+LIN-CD14+"),bwFac = 2,
                dims=c("FITC-A"), minCountThreshold = 0)



md <- data.frame(sampleNames(gs))
names(md)[1]="file_name"
md$sample_id<-gsub("DCs_*",'',md$file_name)
md$sample_id=gsub("_001.fcs",".fcs",md$sample_id)
md$sample_id=gsub("_2.fcs",".fcs",md$sample_id)
md$sample_id<-gsub("*.fcs",'',md$sample_id)

metadata=read.delim("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/metadata_samples.txt",sep="\t",header = T)
md=merge(md,metadata,by="sample_id",all.x=T)
table(duplicated(md$sample_id))

save_gs(gs,path="DCs_gs")

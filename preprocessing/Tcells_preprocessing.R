
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

# here I normalized some channels to provide as input to the unsupervised analysis a normalized gating set object

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


save_gs(gs,path="Tcells_gs")



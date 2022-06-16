
setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data/")

dir.create("Bcells_unsup")

setwd("Bcells_unsup")


library(openCyto)
library(CytoML)
library(flowWorkspace)
library(BiocManager)
library(devtools) #load it
library(reshape2)
library(stats)
library(cytolib)
library(flowCore)
library(CytoExploreR)
library(na.tools)
library(CytoML)
library(ggcyto)
library(data.table)
library(ggplot2)
library(ggridges)
library(dplyr)
library(CytoML)
library(CytoExploreR)
library("openCyto")
library(flowWorkspace)
library(flowCore)
library("CATALYST")
library(tidyverse)
library(data.table)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(data.table)
library(stats)
library(moments)
library(rstatix)
library(ggpubr)
library("bestNormalize")
library(car)
library("multcomp")
library("gtools")

gs=load_gs("../Bcells_gs/")
plot(gs)


md <- data.frame(sampleNames(gs))
names(md)[1]="file_name"
md$sample_id<-gsub("Bcells_*",'',md$file_name)
md$sample_id=gsub("_001.fcs",".fcs",md$sample_id)
md$sample_id=gsub("_2.fcs",".fcs",md$sample_id)
md$sample_id<-gsub("*.fcs",'',md$sample_id)

metadata=read.delim("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/metadata_samples.txt",sep="\t",header = T)
md=merge(md,metadata,by="sample_id",all.x=T)

write.table(md,"md.txt",sep="\t",quote=F,col.names = NA)

md$sample=gsub("(.*)-.*","\\1",md$sample_id)


Bcells <- getData(gs,"Bcells")
trans <- cyto_transformer_extract(gs)


inv <- transformList(names(trans), lapply(trans, `[[`, "inverse"))
Bcells=cytoset_to_flowSet(Bcells)
Bcells<- transform(Bcells, inv)
fs=Bcells

saveRDS(fs, file = "fs.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

fs=readRDS("fs.rds")
md=read.delim("md.txt")

metaDATA=as.data.frame(cbind("file_name"=md$file_name,"sample_id"=md$sample_id,"condition"=md$condition))
row.names(metaDATA)=metaDATA$file_name

panel<-guessPanel(fs[[1]])
all(panel$name %in% colnames(fs))


panel$marker_class <-ifelse(grepl("(FALSE)",panel$use_channel),"state",ifelse(grepl("(CD45RA)",panel$antigen),"type",ifelse(grepl("(CD45)",panel$antigen),"state",ifelse(grepl("(HLADR)",panel$antigen),"state",ifelse(grepl("(LIN)",panel$antigen),"state",ifelse(grepl("(CD19)",panel$antigen),"state","type"))))))
panel$antigen <- ifelse(is.na(panel$antigen), panel$fcs_colname,panel$antigen)

panel<-panel[,c(1,3,5)]


sce<-prepData(
  fs,
  panel = panel,
  md = metaDATA,
  features = panel$fcs_colname,
  transform = T,
  cofactor =150,
  panel_cols = list(channel = "fcs_colname", antigen = "antigen",marker_class="marker_class"),
  md_col = list(factors = "condition"),
  FACS = T
)


pdf("plot_counts.pdf",10,8)
plotCounts(sce, group_by = "sample_id", color_by = "condition")
dev.off()

fsApply(fs, identifier)
keyword(fs, "FILENAME")
colnames(fs)



pdf("expr_heatmap.pdf",25,12)
plotExprHeatmap(sce,features="type",by="sample_id",scale = "last",
                hm_pal = rev(hcl.colors(10, "YlGnBu")))
dev.off()

xdim =12
ydim =12
clusters ="meta20"

set.seed(12354)
sce <- cluster(sce, features = "type",
               xdim = xdim, ydim = ydim, maxK = 40, seed = 5678)


pdf("area_under_curve2.pdf")
delta_area(sce)
dev.off()


pdf("ExprHeatmap_clusters_last.pdf",9,7)
plotExprHeatmap(sce, features = "type", scale="last",col_clust = F,
                by = "cluster_id", k = clusters, bin_anno = T,
                bars = TRUE, perc = TRUE)
dev.off()

pdf("ExprHeatmap_clusters_never.pdf",9,7)
plotExprHeatmap(sce, features = "type", scale="never",col_clust = F,
                by = "cluster_id", k = clusters, bin_anno = T,
                bars = TRUE, perc = TRUE)
dev.off()

pdf("ExprHeatmap_clusters_first.pdf",9,7)
plotExprHeatmap(sce, features = "type", scale="first",col_clust = F,
                by = "cluster_id", k = clusters, bin_anno = T,
                bars = TRUE, perc = TRUE)
dev.off()

pdf("ClusterExprs.pdf",9,7)
plotClusterExprs(sce, k = clusters, features ="type")
dev.off()

pdf("plotAbundances.pdf",15,7)
plotAbundances(sce, k = clusters, by = "sample_id", group_by = "condition")
dev.off()


pdf("plotAbundances_by_cluster2.pdf",9,7)
plotAbundances(sce, k = clusters, by = "cluster_id", 
               group_by = "condition", shape_by = "sample_id")
dev.off()



pdf("Multimap.pdf",15,12)
plotMultiHeatmap(sce, 
                 hm1 = "type", hm2 = "abundances", 
                 k = clusters, row_anno = T, bars = TRUE, perc = TRUE,
                 col_dend = c(F, TRUE))
dev.off()

p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 6                   
p



set.seed(1234)
sce <- runDR(sce, "TSNE", cells = 500, features = "type")

set.seed(1234)
sce <- runDR(sce, "UMAP", cells = 500, features = "type")


pdf("tsne.pdf",10.8)
plotDR(sce, "TSNE", color_by =clusters)
dev.off()

pdf("UMAP.pdf",10.8)
plotDR(sce, "UMAP", color_by =clusters)
dev.off()


pdf("tsne-byCond.pdf",10,8)
plotDR(sce, "TSNE", color_by = clusters, facet_by = "condition")
dev.off()


pdf("UMAP-byCond2.pdf",10,8)
plotDR(sce, "UMAP", color_by = clusters, facet_by = "condition")
dev.off()


cells_by_sample <- split(seq_len(ncol(sce)),sce@colData$sample_id)

cells_by_cond <- split(seq_len(ncol(sce)), sce@colData$condition)

len=t(as.data.frame(rapply(cells_by_cond, length, how="list")))
min=min(len)
N_max <- min

cells_to_use <- sapply(cells_by_cond, function(cs) {
  Ni <- length(cs)
  sample(cs, N_max)
})


pdf("tsne-byCond_eq_size.pdf",10,8)
plotDR(sce[,cells_to_use], "TSNE", color_by =clusters, facet_by = "condition",scale=T)
dev.off()


pdf("UMAP-byCond_eq_size.pdf",10,8)
plotDR(sce[,cells_to_use], "UMAP", color_by = clusters, facet_by = "condition")
dev.off()


markers<-(rownames(sce))[-c(1:4,13)]

pdf("tSNE_markers.pdf",15,10)
plotDR(sce, color_by = markers, ncol = 5,scale = TRUE)#,facet_by = "condition")
dev.off()


pdf("UMAP_markers.pdf",15,10)
plotDR(sce, color_by = markers, ncol = 4,scale = TRUE)#,facet_by = "condition")
dev.off()



pdf("plot_Abundance2.pdf",10,8)
plotAbundances(sce, k = clusters, by = "cluster_id")
dev.off()


merging_table1=read.table("merging_Bcells.txt",sep="\t",header = T)

clusters ="meta20"

sce <- mergeClusters(sce, k = clusters, 
                     table = merging_table1, id = "merging1",overwrite = T)



library(RColorBrewer)
pdf("ExprHeatmap_merged_last_merging.pdf",15,10)
plotExprHeatmap(sce, features = "type", scale="last",col_clust = F,
                by = "cluster_id", k ="merging1", bin_anno = T,
                bars = TRUE, perc = TRUE)
dev.off()


pdf("UMAP_merged.pdf",10.8)
plotDR(sce, "UMAP", color_by ="merging1")
dev.off()

pdf("merging_strategy.pdf",6,9)
plotExprHeatmap(sce, features = "type", scale="last",col_clust = F,
                by = "cluster_id", k =clusters,m= "merging1", bin_anno = F,
                bars = TRUE, perc = F)
dev.off()

pdf("Bcells_Heatmap.pdf",10,3.4)
  plotExprHeatmap(sce, features = "type", scale="last",col_clust = F,
                  by = "cluster_id", k ="merging1", bin_anno = F,
                  bars = TRUE, perc = TRUE)
  dev.off()

 
  tiff("UMAP_Bcells.tiff",width = 8, height = 4,units = 'in',res=200)
  plotDR(sce, "UMAP", color_by ="merging1")
  dev.off()
  
pdf("ClusterExprs_merged.pdf",9,7)
plotClusterExprs(sce, k = "merging1", features ="type")
dev.off()


saveRDS(sce, file = "sce_Bcells.rds", ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL)


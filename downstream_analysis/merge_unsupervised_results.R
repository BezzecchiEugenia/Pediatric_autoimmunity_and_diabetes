library(data.table)
library(dplyr)
library(plyr)

setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data")
path="/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/data"
files=list.files(path,"*unsup_stats.txt",recursive = T)

files

list=list()

for (i in 1:5) {
 # i=1
  list[[i]]=read.delim(files[i])
  list[[i]]= list[[i]][!grepl("sex", list[[i]]$variable),]
  #list[[i]]= list[[i]][!grepl("_001.fcs", list[[i]]$file_name),]
  #list[[i]]= list[[i]][!grepl("_2.fcs", list[[i]]$file_name),]
  list[[i]]= list[[i]][!duplicated(list[[i]]),]
  list[[i]]=unique(setDT(list[[i]]), by = c("sample_id", "variable"))
  #list[[i]]= list[[i]][!duplicated(list[[i]]),]
}


samp1=list[[1]]$sample_id
samp1=samp1[!duplicated(samp1)]

samp2=list[[2]]$sample_id
samp2=samp2[!duplicated(samp2)]

samp3=list[[3]]$sample_id
samp3=samp3[!duplicated(samp3)]

samp4=list[[4]]$sample_id
samp4=samp4[!duplicated(samp4)]

samp5=list[[5]]$sample_id
samp5=samp5[!duplicated(samp5)]

intersection=Reduce(intersect, list(samp1,samp2,samp3,samp4,samp5))

list[[1]]=list[[1]][list[[1]]$sample_id %in% intersection,]
list[[2]]=list[[2]][list[[2]]$sample_id %in% intersection,]
list[[3]]=list[[3]][list[[3]]$sample_id %in% intersection,]
list[[4]]=list[[4]][list[[4]]$sample_id %in% intersection,]
list[[5]]=list[[5]][list[[5]]$sample_id %in% intersection,]

for (i in 1:5) {
  list[[i]]=data.frame(list[[i]])
}


Bcells=list()

samples=list[[1]]$file_name
samples=samples[!duplicated(samples)]
n=length(samples)

for (i in 1:n) {
#i=1
j=samples[i]
Bcells[[i]]=list[[1]][grepl(j,list[[1]]$file_name),]
Bcells[[i]]= Bcells[[i]][!duplicated(Bcells[[i]]),]
q=as.character(Bcells[[i]][1,1])
Bcells[[i]]=Bcells[[i]][,c(2,3)]
colnames(Bcells[[i]])[2]=q
Bcells[[i]]=data.frame(Bcells[[i]])
vector_names=Bcells[[i]]$variable
row.names(Bcells[[i]])=vector_names
Bcells[[i]]$variable=NULL
Bcells[[i]]$pop=rownames(Bcells[[i]])
Bcells[[i]]=Bcells[[i]][order(Bcells[[i]][,'pop']), ]
Bcells[[i]]$pop=NULL
}

Bcells=bind_cols(Bcells)


DCs=list()

samples=list[[2]]$file_name
samples=samples[!duplicated(samples)]
n=length(samples)

for (i in 1:n) {
#i=1
j=samples[i]
DCs[[i]]=list[[2]][grepl(j,list[[2]]$file_name),]
DCs[[i]]= DCs[[i]][!duplicated(DCs[[i]]),]
q=as.character(DCs[[i]][1,1])
DCs[[i]]=DCs[[i]][,c(2,3)]
colnames(DCs[[i]])[2]=q
DCs[[i]]=data.frame(DCs[[i]])
vector_names=DCs[[i]]$variable
row.names(DCs[[i]])=vector_names
DCs[[i]]$variable=NULL
DCs[[i]]$pop=rownames(DCs[[i]])
DCs[[i]]=DCs[[i]][order(DCs[[i]][,'pop']), ]
DCs[[i]]$pop=NULL
}

DCs=bind_cols(DCs)

NK=list()

samples=list[[3]]$file_name
samples=samples[!duplicated(samples)]
n=length(samples)

for (i in 1:n) {
#i=1
j=samples[i]
NK[[i]]=list[[3]][grepl(j,list[[3]]$file_name),]
NK[[i]]= NK[[i]][!duplicated(NK[[i]]),]
q=as.character(NK[[i]][1,1])
NK[[i]]=NK[[i]][,c(2,3)]
colnames(NK[[i]])[2]=q
NK[[i]]=data.frame(NK[[i]])
vector_names=NK[[i]]$variable
row.names(NK[[i]])=vector_names
NK[[i]]$variable=NULL
NK[[i]]$pop=rownames(NK[[i]])
NK[[i]]=NK[[i]][order(NK[[i]][,'pop']), ]
NK[[i]]$pop=NULL
}

NK=bind_cols(NK)

Tcells=list()

samples=list[[4]]$file_name
samples=samples[!duplicated(samples)]
n=length(samples)

for (i in 1:n) {
#i=1
j=samples[i]
Tcells[[i]]=list[[4]][grepl(j,list[[4]]$file_name),]
Tcells[[i]]= Tcells[[i]][!duplicated(Tcells[[i]]),]
q=as.character(Tcells[[i]][1,1])
Tcells[[i]]=Tcells[[i]][,c(2,3)]
colnames(Tcells[[i]])[2]=q
Tcells[[i]]=data.frame(Tcells[[i]])
vector_names=Tcells[[i]]$variable
row.names(Tcells[[i]])=vector_names
Tcells[[i]]$variable=NULL
#order=sort(rownames(Tcells[[i]]))
Tcells[[i]]$pop=rownames(Tcells[[i]])
Tcells[[i]]=Tcells[[i]][order(Tcells[[i]][,'pop']), ]
Tcells[[i]]$pop=NULL
}

Tcells=bind_cols(Tcells)

TREG=list()

samples=list[[5]]$file_name
samples=samples[!duplicated(samples)]
n=length(samples)

for (i in 1:n) {
#i=1
j=samples[i]
TREG[[i]]=list[[5]][grepl(j,list[[5]]$file_name),]
TREG[[i]]= TREG[[i]][!duplicated(TREG[[i]]),]
q=as.character(TREG[[i]][1,1])
TREG[[i]]=TREG[[i]][,c(2,3)]
colnames(TREG[[i]])[2]=q
TREG[[i]]=data.frame(TREG[[i]])
vector_names=TREG[[i]]$variable
row.names(TREG[[i]])=vector_names
TREG[[i]]$variable=NULL
TREG[[i]]$pop=rownames(TREG[[i]])
TREG[[i]]=TREG[[i]][order(TREG[[i]][,'pop']), ]
TREG[[i]]$pop=NULL
}

TREG=bind_cols(TREG)

merged=rbind(Tcells,TREG,NK,DCs,Bcells)


write.table(merged,"merged_unsup_intersect.txt",sep="\t",quote = F,row.names = T)

library(pROC)
library(ROCR)
library(data.table)
#BiocManager::install("PAGWAS")
#BiocManager::install("verification")
library("verification")
library("PRROC")#per le curve

library("olsrr")


setwd("/beegfs/scratch/ric.petrelli/ric.petrelli/flowCytometry/OPENCYTO_final/ULTIMATE_analysis/data")

#######ROC curve with my data

##this line of code uses merged unsup results (saved in IDAA1C/AUC_single_pop)
data=t(read.delim("merged_unsup_intersect.txt"))

data=as.data.frame(data)
data$sample_id=row.names(data)

META=read.delim("../meta_to_merge_final_V2.txt")
Idaa1C=META[,c(1,32)]
Idaa1C=na.omit(Idaa1C)

BMI=read.delim("../BMI.txt")
BMI$BMI=as.numeric(BMI$BMI)
META=merge(META,BMI,by="sample_id")

data=merge(data,META,by="sample_id")
DF=data
DF=DF[!grepl("OSR",DF$sample_id),]
DF=DF[!grepl("preT1D",DF$sample_id),]
DF=DF[grepl("T1D",DF$sample_id),]

DF=DF[!is.na(DF$IDAA1C),]
DF$labels=ifelse(DF$IDAA1C>9,1,0)

DF$IDAA1C=NULL
row.names(DF)=DF$sample_id
DF$sample_id=NULL

ind = sample(2, nrow(DF), replace=TRUE, prob=c(0.7, 0.3))
train = DF[ind==1,]
test = DF[ind==2,]

DF=DF

r=data.frame(table(colnames(DF)))
r$Var1=as.factor(r$Var1)
dim=dim(r)[1]

DF[,dim]=DF$labels
AUC=list()

     
 i=58 ## extract naive CD4
 var=as.character(r$Var1[i])
 a=DF[,var]
 IDAA1C=DF$labels
 data=cbind(a,IDAA1C)
 data=data.frame(data)
 data=na.omit(data)
 mod1=glm(IDAA1C~a,family = binomial(link = "logit"),data=data)

summary(mod1)

pred <- prediction(predict(mod1), data$IDAA1C)
perf <- performance(pred,"tpr","fpr")

 
pdf(paste0(var,"_AUC.pdf"),5,4)
PRROC_obj <- roc.curve(scores.class0 = predict(mod1), weights.class0=data$IDAA1C,
                       curve=TRUE)
plot(PRROC_obj,main=var)
dev.off()



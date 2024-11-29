## exploratory data analysis and ICA
rm(list=ls())
##================================================ 
## Libraries
setwd("Z:/users/petr.nazarov/scripts")
source("Lib/plotPCA.r")
source("Lib/plotDataPDF.r")
source("Lib/drawTable.r")
source("Lib/summarizeRows.r")
source("Lib/violinplot.r")
source("Lib/LibDEA.r")
#source("Lib/LibICA.r")
#source("Lib/enrichGO.r")
library(pheatmap)
library(consICA)
library(SummarizedExperiment)
##================================================
## Constants
path = "Z:/users/petr.nazarov/results/multiomicsMAE_2024-03-22"

setwd(path)

##================================================
## Analysis mRNA 
load("data_mRNA.RData")

## PCA
pdf("PCA_mRNA.pdf",width=8.3, height=11.7,onefile=TRUE)
PC = prcomp(t(mRNA$X),scale. = FALSE)
par(mfcol=c(3,2))
col = c(FEMALE="#FF000088",MALE="#0000FF88")[as.character(mRNA$Var$Sex)]
cex = 1+0.1*as.integer(mRNA$Var$LibSize)
plotPCA(mRNA$X,pc=PC,col=col,cex=cex,main="PCA (1,2)",cex.main=1,cex.names=0.1)
plotPCA(mRNA$X,pc=PC,ipc=c(3:4),col=col,cex=cex,main="PCA (3,4)",cex.main=1,cex.names=0.1)
plotPCA(mRNA$X,pc=PC,ipc=c(5:6),col=col,cex=cex,main="PCA (5,6)",cex.main=1,cex.names=0.1)
plotPCA(mRNA$X,pc=PC,ipc=c(7:8),col=col,cex=cex,main="PCA (7,8)",cex.main=1,cex.names=0.1)
plotPCA(mRNA$X,pc=PC,ipc=c(9:10),col=col,cex=cex,main="PCA (9,10)",cex.main=1,cex.names=0.1)
explvar = PC$sdev^2 / sum(PC$sdev^2)
names(explvar) = colnames(PC$x)
barplot(explvar[1:10],las=2,cex.names =0.6,col="skyblue",main="Explained variability",cex.main=1)
dev.off()

## here 1 - rerun ICA, 0 - load the data
if(1){
  se <- SummarizedExperiment(assays = List(X=mRNA$X), colData = mRNA$Var)
  set.seed(2024)
  cica <- consICA(se, ncomp = 20, ntry = 100, ncores = 4)
  save(cica,file=sprintf("IC%d(100)_mRNA.RData",ncol(cica$S)))
  saveReport(cica,Var=mRNA$Var, genenames = mRNA$genes, file=sprintf("report_mRNA_IC%d.pdf",ncol(cica$S)))
  cica <- getGO(cica, alpha = 0.01, genenames=mRNA$genes, ncores = 4)
  saveReport(cica,Var=mRNA$Var, genenames = mRNA$genes, file=sprintf("report_mRNA_IC%d.pdf",ncol(cica$S)))
  save(cica,file = sprintf("IC%d(100)_mRNA.RData",ncol(cica$S)))
}else{
  load("IC20(100)_mRNA.RData")
  saveReport(cica,Var=mRNA$Var,genenames=mRNA$genes,file=sprintf("report_mRNA_IC%d.pdf",ncol(cica$S)))
}


##================================================
## Analysis miR 
load("data_miR.RData")

## PCA
pdf("PCA_miR.pdf",width=8.3, height=11.7,onefile=TRUE)
PC = prcomp(t(miR$X),scale. = FALSE)
par(mfcol=c(3,2))
col = c(FEMALE="#FF000088",MALE="#0000FF88")[as.character(miR$Var$Sex)] 
cex = 1+0.1*as.integer(miR$Var$LibSize)
plotPCA(miR$X,pc=PC,col=col,cex=cex,main="PCA (1,2)",cex.main=1,cex.names=0.1)
plotPCA(miR$X,pc=PC,ipc=c(3:4),col=col,cex=cex,main="PCA (3,4)",cex.main=1,cex.names=0.1)
plotPCA(miR$X,pc=PC,ipc=c(5:6),col=col,cex=cex,main="PCA (5,6)",cex.main=1,cex.names=0.1)
plotPCA(miR$X,pc=PC,ipc=c(7:8),col=col,cex=cex,main="PCA (7,8)",cex.main=1,cex.names=0.1)
plotPCA(miR$X,pc=PC,ipc=c(9:10),col=col,cex=cex,main="PCA (9,10)",cex.main=1,cex.names=0.1)
explvar = PC$sdev^2 / sum(PC$sdev^2)
names(explvar) = colnames(PC$x)
barplot(explvar[1:10],las=2,cex.names =0.6,col="skyblue",main="Explained variability",cex.main=1)
dev.off()

## consICA 
## here 1 - rerun ICA, 0 - load the data
if(1){
  se <- SummarizedExperiment(assays = List(X=miR$X), colData = miR$Var)
  set.seed(2024)
  cica <- consICA(se, ncomp = 20, ntry = 100, ncores = 4)
  save(cica,file=sprintf("IC%d(100)_miR.RData",ncol(cica$S)))
  saveReport(cica,Var=miR$Var, file=sprintf("report_miR_IC%d.pdf",ncol(cica$S)))
}else{
  load("IC20(100)_miR.RData")
  saveReport(cica,Var=miR$Var, file=sprintf("report_miR_IC%d.pdf",ncol(cica$S)))
}

##================================================
## Analysis Protein data 
load("data_Prot.RData")

## PCA
pdf("PCA_Prot.pdf",width=8.3, height=11.7,onefile=TRUE)
PC = prcomp(t(Prot$X),scale. = FALSE)
par(mfcol=c(3,2))
col = c(FEMALE="#FF000088",MALE="#0000FF88")[as.character(Prot$Var$Sex)]
cex = apply(Prot$X,2,sum)
cex = 1 + (cex-min(cex))/(max(cex)-min(cex))

plotPCA(Prot$X,pc=PC,col=col,cex=cex,main="PCA (1,2)",cex.main=1,cex.names=0.1)
plotPCA(Prot$X,pc=PC,ipc=c(3:4),col=col,cex=cex,main="PCA (3,4)",cex.main=1,cex.names=0.1)
plotPCA(Prot$X,pc=PC,ipc=c(5:6),col=col,cex=cex,main="PCA (5,6)",cex.main=1,cex.names=0.1)
plotPCA(Prot$X,pc=PC,ipc=c(7:8),col=col,cex=cex,main="PCA (7,8)",cex.main=1,cex.names=0.1)
plotPCA(Prot$X,pc=PC,ipc=c(9:10),col=col,cex=cex,main="PCA (9,10)",cex.main=1,cex.names=0.1)
explvar = PC$sdev^2 / sum(PC$sdev^2)
names(explvar) = colnames(PC$x)
barplot(explvar[1:10],las=2,cex.names =0.6,col="skyblue",main="Explained variability",cex.main=1)
dev.off()


##consICA
## here 1 - rerun ICA, 0 - load the data
if(1){
  se <- SummarizedExperiment(assays = List(X=Prot$X), colData = Prot$Var)
  set.seed(2024)
  cica <- consICA(se, ncomp = 20, ntry = 100, ncores = 4)
  save(cica,file=sprintf("IC%d(100)_Prot.RData",ncol(cica$S)))
  saveReport(cica,Var=Prot$Var, file=sprintf("report_Prot_IC%d.pdf",ncol(cica$S)))
}else{
  load("IC20(100)_Prot.RData")
  saveReport(cica,Var=Prot$Var, file=sprintf("report_Prot_IC%d.pdf",ncol(cica$S)))
}

##================================================
## Analysis Methylation data 
load("data_Methyl.RData")

## PCA
pdf("PCA_Methyl.pdf",width=8.3, height=11.7,onefile=TRUE)
PC = prcomp(t(Methyl$X),scale. = FALSE)
par(mfcol=c(3,2))
col = c(FEMALE="#FF000088",MALE="#0000FF88")[as.character(Methyl$Var$Sex)]
cex = apply(Methyl$X,2,sum)
cex = 1 + (cex-min(cex))/(max(cex)-min(cex))

plotPCA(Methyl$X,pc=PC,col=col,cex=cex,main="PCA (1,2)",cex.main=1,cex.names=0.1)
plotPCA(Methyl$X,pc=PC,ipc=c(3:4),col=col,cex=cex,main="PCA (3,4)",cex.main=1,cex.names=0.1)
plotPCA(Methyl$X,pc=PC,ipc=c(5:6),col=col,cex=cex,main="PCA (5,6)",cex.main=1,cex.names=0.1)
plotPCA(Methyl$X,pc=PC,ipc=c(7:8),col=col,cex=cex,main="PCA (7,8)",cex.main=1,cex.names=0.1)
plotPCA(Methyl$X,pc=PC,ipc=c(9:10),col=col,cex=cex,main="PCA (9,10)",cex.main=1,cex.names=0.1)
explvar = PC$sdev^2 / sum(PC$sdev^2)
names(explvar) = colnames(PC$x)
barplot(explvar[1:10],las=2,cex.names =0.6,col="skyblue",main="Explained variability",cex.main=1)
dev.off()

rm(PC); gc()

##consICA
## here 1 - rerun ICA, 0 - load the data
if(1){
  se <- SummarizedExperiment(assays = List(X=Methyl$X), colData = Methyl$Var)
  set.seed(2024)
  cica <- consICA(se, ncomp = 10, ntry = 100, ncores = 4)
  save(cica,file=sprintf("IC%d(100)_Methyl.RData",ncol(cica$S)))
  saveReport(cica,Var=Methyl$Var, file=sprintf("report_Methyl_IC%d.pdf",ncol(cica$S)))
  rm(cica);gc()
  
  cica <- consICA(se, ncomp = 20, ntry = 100, ncores = 4)
  save(cica,file=sprintf("IC%d(100)_Methyl.RData",ncol(cica$S)))
  saveReport(cica,Var=Methyl$Var, file=sprintf("report_Methyl_IC%d.pdf",ncol(cica$S)))
  rm(cica);gc()
  
  cica <- consICA(se, ncomp = 30, ntry = 100, ncores = 4)
  save(cica,file=sprintf("IC%d(100)_Methyl.RData",ncol(cica$S)))
  saveReport(cica,Var=Methyl$Var, file=sprintf("report_Methyl_IC%d.pdf",ncol(cica$S)))
  rm(cica);gc()
  
}else{
  load("IC20(100)_Methyl.RData")
  saveReport(cica,Var=Methyl$Var, file=sprintf("report_Methyl_IC%d.pdf",ncol(cica$S)))
}

##================================================
## Analysis Metabolomics data 
load("data_Metab.RData")

## PCA
pdf("PCA_Metab.pdf",width=8.3, height=11.7,onefile=TRUE)
PC = prcomp(t(Metab$X),scale. = FALSE)
par(mfcol=c(3,2))
col = c(FEMALE="#FF000088",MALE="#0000FF88")[as.character(Metab$Var$Sex)]
cex = apply(Metab$X,2,sum)
cex = 1 + (cex-min(cex))/(max(cex)-min(cex))

plotPCA(Metab$X,pc=PC,col=col,cex=cex,main="PCA (1,2)",cex.main=1,cex.names=0.1)
plotPCA(Metab$X,pc=PC,ipc=c(3:4),col=col,cex=cex,main="PCA (3,4)",cex.main=1,cex.names=0.1)
plotPCA(Metab$X,pc=PC,ipc=c(5:6),col=col,cex=cex,main="PCA (5,6)",cex.main=1,cex.names=0.1)
plotPCA(Metab$X,pc=PC,ipc=c(7:8),col=col,cex=cex,main="PCA (7,8)",cex.main=1,cex.names=0.1)
plotPCA(Metab$X,pc=PC,ipc=c(9:10),col=col,cex=cex,main="PCA (9,10)",cex.main=1,cex.names=0.1)
explvar = PC$sdev^2 / sum(PC$sdev^2)
names(explvar) = colnames(PC$x)
barplot(explvar[1:10],las=2,cex.names =0.6,col="skyblue",main="Explained variability",cex.main=1)
dev.off()


##consICA
if(1){
  se <- SummarizedExperiment(assays = List(X=Metab$X), colData = Metab$Var)
  set.seed(2024)
  cica <- consICA(se, ncomp = 20, ntry = 100, ncores = 4)
  save(cica,file=sprintf("IC%d(100)_Metab.RData",ncol(cica$S)))
  saveReport(cica,Var=Metab$Var, file=sprintf("report_Metab_IC%d.pdf",ncol(cica$S)))
}else{
  load("IC20(100)_Metab.RData")
  saveReport(cica,Var=Metab$Var, file=sprintf("report_Metab_IC%d.pdf",ncol(cica$S)))
}

##================================================
## Analysis CGH data 
load("data_aCGH.RData")

## PCA
pdf("PCA_aCGH.pdf",width=8.3, height=11.7,onefile=TRUE)
PC = prcomp(t(aCGH$X),scale. = FALSE)
par(mfcol=c(3,2))
col = c(FEMALE="#FF000088",MALE="#0000FF88")[as.character(aCGH$Var$Sex)]
cex = apply(aCGH$X,2,sum)
cex = 1 + (cex-min(cex))/(max(cex)-min(cex))

plotPCA(aCGH$X,pc=PC,col=col,cex=cex,main="PCA (1,2)",cex.main=1,cex.names=0.1)
plotPCA(aCGH$X,pc=PC,ipc=c(3:4),col=col,cex=cex,main="PCA (3,4)",cex.main=1,cex.names=0.1)
plotPCA(aCGH$X,pc=PC,ipc=c(5:6),col=col,cex=cex,main="PCA (5,6)",cex.main=1,cex.names=0.1)
plotPCA(aCGH$X,pc=PC,ipc=c(7:8),col=col,cex=cex,main="PCA (7,8)",cex.main=1,cex.names=0.1)
plotPCA(aCGH$X,pc=PC,ipc=c(9:10),col=col,cex=cex,main="PCA (9,10)",cex.main=1,cex.names=0.1)
explvar = PC$sdev^2 / sum(PC$sdev^2)
names(explvar) = colnames(PC$x)
barplot(explvar[1:10],las=2,cex.names =0.6,col="skyblue",main="Explained variability",cex.main=1)
dev.off()

rm(PC); gc()

## aCGH excluded from deconvolution
## consICA
# if(1){
#   se <- SummarizedExperiment(assays = List(X=aCGH$X), colData = aCGH$Var)
#   set.seed(2024)
#   cica <- consICA(se, ncomp = 20, ntry = 100, ncores = 4)
#   save(cica,file=sprintf("IC%d(100)_aCGH.RData",ncol(cica$S)))
#   saveReport(cica,Var=aCGH$Var, file=sprintf("report_aCGH_IC%d.pdf",ncol(cica$S)))
# }else{
#   load("IC20(100)_aCGH.RData")
#   saveReport(cica,Var=aCGH$Var, file=sprintf("report_aCGH_IC%d.pdf",ncol(cica$S)))
# }

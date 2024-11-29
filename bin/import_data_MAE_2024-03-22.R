## exploratory data analysis, import from MultiAssayExperiment
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
source("Lib/Venn.r")
library(pheatmap)
library(MultiAssayExperiment)
library(S4Vectors)
##================================================
## Constants
path = "Z:/users/petr.nazarov/results/multiomicsMAE_2024-03-22"
dir.create(path)
setwd(path)


##================================================
## Load MAE

mae = readRDS("Z:/omics/MultiAssayExperiments/h5MAE/mae_mae.rds")
names(mae)

MAE <- loadHDF5MultiAssayExperiment(dir = "Z:/omics/MultiAssayExperiments/h5MAE",prefix = NULL)

experiments(MAE)
names(experiments(MAE))

##------------------------------------------------
## mRNA (choose one) -> g
mae_mrna = "mRNA-seq-2 | batch-adjusted | vst | corrected"
file_mrna_anno = "Z:/users/petr.nazarov/libs/mart_export_hs_221111.txt"
thr.mrna = -2

mae_mirna = "miRNA-seq, mature | batch-adjusted | corrected"
thr.mir = -2

## Protein (choose one) -> p
mae_prot = "Proteomics | imputed missing values"
file_prot_anno = "Z:/omics/Proteomics_IMTM/Upload_2023-02-02/proteins_annotation.txt"
 
## Methylation -> e
mae_methyl = "EM-seq | cell-type adjusted | 100,000 most variable CpG sites"
#mae_methyl = "EM-seq | 100,000 random CpG sites"


## Metabolomics -> m
# [1] "Acylcarnitines | with missing values"                                    
# [2] "Acylcarnitines | unscaled"                                               
# [3] "Amino acids | with missing values"                                       
# [4] "Amino acids | unscaled"                                                  
# [5] "Very long chain fatty acids"                                             
# [6] "Very long chain fatty acids | unscaled"                                  
# [7] "Acylcarnitines | batch-adjusted | with missing values"                   
# "Amino acids | batch-adjusted | with missing values"                      

mae_metab = c(
    acylcarnitines = "Acylcarnitines | batch-adjusted | imputed missing values",
    aminoacids = "Amino acids | batch-adjusted | imputed missing values",
    fattyacids = "Very long chain fatty acids | batch-adjusted",
    lipidspos = "Lipidomics, positive | transformed",                                      
    lipidsneg = "Lipidomics, negative | transformed")

## CGH -> c
mae_acgh = "aCGH"


##================================================
## Phenotype Data Import 

Pheno = as.data.frame(colData(MAE))
rownames(Pheno) = Pheno[[1]]
save(Pheno,file="data_Pheno.RData")

##================================================
## mRNA Data Import 

mRNA=list()
mRNA$name = mae_mrna
mRNA$Anno = NA ## feature annotation
mRNA$genes = NA ## vector of feature names
mRNA$Pheno = NA ## sample annotation
mRNA$Var = NA ## variables of sample
mRNA$X0 = NA ## original data
mRNA$X  = NA ## normalized data

##---------------------------------
## read counts
mRNA$X0 = as.matrix(assay(MAE,mae_mrna,withDimnames=TRUE))
range(mRNA$X0)
plot(density(mRNA$X0),main="Raw data",lwd=2)
plot(mRNA$X0["ENSG00000229807",])
#plot(mRNA$X0["ENSG00000258992",])
#plot(mRNA$X0["ENSG00000258992",])
rownames(mRNA$X0)

boxplot(mRNA$X0)

plot(apply(mRNA$X0,1,mean))
plot(apply(mRNA$X0,1,sd))

#X = scale(mRNA$X0)
#min(X)
#max(X)

#mRNA$X0[min(mRNA$X0)<0] = 0
sum(is.na(mRNA$X0))

#plot(density(log2(1+ mRNA$X0)),main="Raw data",lwd=2)

png("mRNA.png")
par(mfcol=c(1,2))
plotDataPDF(mRNA$X0,xlim=c(-2,2),ylim=c(0,3))
boxplot(mRNA$X0,las=2,outline=FALSE,col=1:ncol(mRNA$X0))
dev.off()

## normalize, log-transform, filter
#mRNA$X = norm.Counts(mRNA$X0,method="DESeq",doLog=T)
mRNA$X = scale(mRNA$X0)
pdf("import_mRNA.pdf", onefile = TRUE)
#plot(density(log2(1+mRNA$X0)),main="Raw data",lwd=2)
plot(density(mRNA$X0),main="Raw(no!) data",lwd=2)
abline(v=thr.mrna,col=2,lty=2)
#boxplot(log2(1+mRNA$X0),main="Raw data",outline=FALSE,las=2)
boxplot(mRNA$X0,main="Raw(no!) data",outline=FALSE,las=2)
ikeep = apply(mRNA$X,1,max)>thr.mrna
table(ikeep)
plot(density(mRNA$X[ikeep,]),main="Normalized / filtered data",lwd=2,col=4)
abline(v=thr.mrna,col=2,lty=2)
boxplot(mRNA$X[ikeep,],main="Normalized / filtered data",outline=FALSE,las=2,col=4)
mRNA$X0 = mRNA$X0[ikeep,]
mRNA$X = mRNA$X[ikeep,]

##----------------------------------- 
## gene annotation
Tab = read.table(file_mrna_anno,header=TRUE,sep="\t",quote="\"")
Tab = Tab[Tab[[1]] %in% rownames(mRNA$X0),]
MADS = summarizeRows(
  data = Tab[,c("Gene.start..bp.","Gene.end..bp.","Strand","Gene...GC.content","Transcript.count")],
  anno = Tab[,c("Gene.stable.ID","Gene.name","Gene.description","Gene.type","Strand","Gene...GC.content","Transcript.count","Chromosome.scaffold.name")],
  sum.by.col = 1)
Tab = data.frame(MADS$anno,MADS$data)
mRNA$Anno = Tab[rownames(mRNA$X),]
mRNA$genes = mRNA$Anno$Gene.name
names(mRNA$genes) = rownames(mRNA$Anno)

##----------------------------------- 
## sample annotation
mRNA$Pheno = Pheno[colnames(mRNA$X),-1]
mRNA$Var = df2factor(mRNA$Pheno, nlev = 8)
mRNA$Var$LibSize = num2fact(apply(mRNA$X0,2,sum),8)
str(mRNA$Var)

par(mfcol=c(2,2))
for (i in 1:ncol(mRNA$Var))
  barplot(summary(mRNA$Var[[i]]),main=names(mRNA$Var)[i],las=2)

dev.off()

save(mRNA,file="data_mRNA.RData")
#load(file="data_mRNA.RData")

##================================================
## miRNA Data Import 
miR=list()
miR$Anno = NA ## feature annotation
miR$Pheno = NA ## sample annotation
miR$Var = NA ## variables of sample
miR$X0 = NA ## original data
miR$X  = NA ## normalized data

##---------------------------------
## read counts (no! 2024-02-13 it looks like Lorenz now)
miR$X0 = as.matrix(assay(MAE,mae_mirna,withDimnames=TRUE))
range(miR$X0)
plot(density(miR$X0))
boxplot(miR$X0)
#miR$X0[miR$X0<0]=0
plot(density(miR$X0))
boxplot(miR$X0)
tmpfeature=rownames(miR$X0)
rownames(miR$X0) = gsub("[.]","-",gsub("^hsa.|_.+","",rownames(miR$X0)))

## normalize, log-transform, filter
#miR$X = norm.Counts(miR$X0,method="DESeq",doLog=F)
miR$X = scale(miR$X0)
pdf("import_miR.pdf", onefile = TRUE)
plot(density((miR$X0)),main="Raw data (no!)",lwd=2)
abline(v=thr.mir,col=2,lty=2)
boxplot((miR$X0),main="Raw data (no!)",outline=FALSE,las=2)
ikeep = apply(miR$X,1,max)>thr.mir
table(ikeep)
plot(density(miR$X[ikeep,]),main="Normalized / filtered data",lwd=2,col=4)
abline(v=thr.mir,col=2,lty=2)
boxplot(miR$X[ikeep,],main="Normalized / filtered data",outline=FALSE,las=2,col=4)
miR$X0 = miR$X0[ikeep,]
miR$X = miR$X[ikeep,]
tmpfeature = tmpfeature[ikeep]
##----------------------------------- 
## miR annotation
miR$Anno = t(as.data.frame(strsplit(sub("Homo_sapiens_","",tmpfeature),"_")))
rownames(miR$Anno) = rownames(miR$X)

str(miR)
##----------------------------------- 
## sample annotation
miR$Pheno = Pheno[colnames(miR$X),-1]
miR$Var = df2factor(miR$Pheno, nlev = 8)
miR$Var$LibSize = num2fact(apply(miR$X0,2,sum),8)
str(miR$Var)

par(mfcol=c(2,2))
for (i in 1:ncol(miR$Var))
  barplot(summary(miR$Var[[i]]),main=names(miR$Var)[i],las=2)

dev.off()

save(miR,file="data_miR.RData")


##================================================
## Protein Data Import 
Prot=list()
Prot$name = mae_prot
Prot$Anno = NA ## feature annotation
Prot$Pheno = NA ## sample annotation
Prot$Var = NA ## variables of sample
Prot$X0 = NA ## original data
Prot$X  = NA ## normalized data

##---------------------------------
## reading prot
Prot$X0 = as.matrix(assay(MAE,mae_prot,withDimnames=TRUE))
c(min(Prot$X0),max(Prot$X0))
Prot$X0[Prot$X0<0]=0
## normalize, log-transform, filter ?
plot(density(Prot$X0))
boxplot(Prot$X0)
Prot$X = scale(Prot$X0)

pdf("import_Prot.pdf", onefile = TRUE)
plot(density(Prot$X),main="Externaly normalized data",lwd=2)
boxplot(Prot$X,main="Externaly normalized data",outline=FALSE,las=2,col=4)
##----------------------------------- 
## Prot annotation
Prot$Anno = read.table(file_prot_anno,sep="\t",header=TRUE)
rownames(Prot$Anno) = Prot$Anno$Accession 
Prot$Anno = Prot$Anno[rownames(Prot$X),c(2,2,3,1,4,5)]
names(Prot$Anno)[2]="Gene"
Prot$Anno$Gene[grep("GN=",Prot$Anno$Description)] = sub(" .+","",sub(".+GN=","",Prot$Anno$Description[grep("GN=",Prot$Anno$Description)]))
str(Prot)
##----------------------------------- 
## sample annotation
Prot$Pheno = Pheno[colnames(Prot$X),-1]
Prot$Var = df2factor(Prot$Pheno, nlev = 8)
Prot$Var$LibSize = num2fact(apply(Prot$X,2,sum),8)
str(Prot$Var)

par(mfcol=c(2,2))
for (i in 1:ncol(Prot$Var))
  barplot(summary(Prot$Var[[i]]),main=names(Prot$Var)[i],las=2)

dev.off()

save(Prot,file="data_Prot.RData")

##================================================
## Mehylation Data Import 
Methyl=list()
Methyl$name = mae_methyl
Methyl$Anno = NA ## feature annotation
Methyl$Pheno = NA ## sample annotation
Methyl$Var = NA ## variables of sample
Methyl$X  = NA ## normalized data

mae_methyl
Methyl$X = as.matrix(assay(MAE,mae_methyl,withDimnames=TRUE))
str(Methyl$X)
sum(is.na(Methyl$X))
##----------------------------------- 
## feature annotation
Methyl$Anno = data.frame(chr = sub("_.+","",rownames(Methyl$X)),coord = sub(".+_","",rownames(Methyl$X))  )
str(Methyl)
pdf("import_Methyl.pdf", onefile = TRUE)
par(mfcol=c(1,1))
plot(density(Methyl$X,width=0.1),main="Externaly normalized data",lwd=2)
boxplot(Methyl$X,main="Externaly normalized data",outline=FALSE,las=2,col=4)

##----------------------------------- 
## sample annotation
Methyl$Pheno = Pheno[colnames(Methyl$X),-1]
Methyl$Var = df2factor(Methyl$Pheno, nlev = 8)
Methyl$Var$LibSize = num2fact(apply(Methyl$X,2,sum),8)
str(Methyl$Var)

par(mfcol=c(2,2))
for (i in 1:ncol(Methyl$Var))
  barplot(summary(Methyl$Var[[i]]),main=names(Methyl$Var)[i],las=2)

dev.off()

save(Methyl,file="data_Methyl.RData")


##================================================
## Metabolomic (targeted) Data Import 
Metab=list()
Metab$Anno = NA ## feature annotation
Metab$Pheno = NA ## sample annotation
Metab$Var = NA ## variables of sample
Metab$X0 = NA ## original data
Metab$X  = NA ## combined data

##---------------------------------
## read data
Metab$X0 = list()
Metab$X  = NA 
ifl=1
pat= NULL

for (ifl in 1: length(mae_metab)){
  #Metab$X0[[ifl]] = as.matrix(read.table(file_metab[ifl],header=TRUE,sep=",",row.names=1,check.names = FALSE))  
  print(mae_metab[ifl])
  Metab$X0[[ifl]] = as.matrix(assay(MAE,mae_metab[ifl],withDimnames=TRUE))
  names(Metab$X0)[ifl] = names(mae_metab)[ifl]
  pat=c(pat,colnames(Metab$X0[[ifl]]))
}
str(Metab)
par(mfcol=c(3,3))
for (ifl in 1: length(mae_metab)){
  plot(density(Metab$X0[[ifl]]),main=names(Metab$X0)[ifl])
}
pat=sort(unique(pat))

censcale = function(x){
  x = (x - mean(x)) / sd(x)
}

ifl=1
Metab$X  = NA
Metab$Anno  = NA
for (ifl in 1:length(Metab$X0)){
  if (is.na(Metab$X)[1]){
    Metab$X = censcale(Metab$X0[[ifl]][,pat])
    Metab$Anno = data.frame(metabolite = rownames(Metab$X0[[ifl]]),
                            class=rep(names(Metab$X0)[ifl],nrow((Metab$X0[[ifl]])))
                            ) 
  }else{
    Metab$X = rbind(Metab$X, censcale(Metab$X0[[ifl]][,pat]))
    Metab$Anno = rbind(Metab$Anno,
                       data.frame(metabolite = rownames(Metab$X0[[ifl]]),
                                  class=rep(names(Metab$X0)[ifl],nrow((Metab$X0[[ifl]]))))
                      )
  }
}
str(Metab$X)
str(Metab$Anno)

## any NA?
sum(is.na(Metab$X))

pdf("import_Metab.pdf", onefile = TRUE)
plot(density(Metab$X,na.rm=TRUE),main="Externaly normalized data",lwd=2)
boxplot(Metab$X,main="Externaly normalized data",outline=FALSE,las=2,col=4)

##----------------------------------- 
## sample annotation
Metab$Pheno = Pheno[colnames(Metab$X),-1]
Metab$Var = df2factor(Metab$Pheno, nlev = 8)
Metab$Var$LibSize = num2fact(apply(Metab$X,2,sum),8)
str(Metab$Var)

par(mfcol=c(2,2))
for (i in 1:ncol(Metab$Var))
  barplot(summary(Metab$Var[[i]]),main=names(Metab$Var)[i],las=2)

dev.off()

save(Metab,file="data_Metab.RData")


##================================================
## CGH Data Import 
aCGH=list()
aCGH$name = mae_acgh
aCGH$Anno = NA ## feature annotation
aCGH$Pheno = NA ## sample annotation
aCGH$Var = NA ## variables of sample
aCGH$X  = NA ## normalized data

mae_acgh
aCGH$X = as.matrix(assay(MAE,mae_acgh,withDimnames=TRUE))

str(aCGH$X)
sum(is.na(aCGH$X))
##----------------------------------- 
## feature annotation
aCGH$Anno = data.frame(id = rownames(aCGH$X), chr = NA  )
str(aCGH)
pdf("import_aCGH.pdf", onefile = TRUE)
par(mfcol=c(1,1))
plot(density(aCGH$X,width=0.1),main="Externaly normalized data",lwd=2)
boxplot(aCGH$X,main="Externaly normalized data",outline=FALSE,las=2,col=4)

##----------------------------------- 
## sample annotation
aCGH$Pheno = Pheno[colnames(aCGH$X),-1]
aCGH$Var = df2factor(aCGH$Pheno, nlev = 8)
aCGH$Var$LibSize = num2fact(apply(aCGH$X,2,sum),8)
str(aCGH$Var)

par(mfcol=c(2,2))
for (i in 1:ncol(aCGH$Var))
  barplot(summary(aCGH$Var[[i]]),main=names(aCGH$Var)[i],las=2)

dev.off()

save(aCGH,file="data_aCGH.RData")

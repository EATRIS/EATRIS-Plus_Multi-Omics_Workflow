## Integration based on ICA 
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
source("Lib/LibML.r")
source("Lib/Venn.r")
library(pheatmap)
library(consICA)
library(SummarizedExperiment)
library(reshape2)
library(igraph)
##================================================
## Constants
path = "Z:/users/petr.nazarov/results/multiomicsMAE_2024-03-22"

##================================================
## Data Import 
setwd(path)

load("data_Pheno.RData")
## combine data
Data = list()
## mRNA - code: g
load("data_mRNA.RData")
Data$mRNA = mRNA; rm(mRNA)
## microRNA - code: u
load("data_miR.RData")
Data$miR = miR; rm(miR)
## Protein - code: p
load("data_Prot.RData")
Data$Prot = Prot; rm(Prot)
## Metabol - code: m
load("data_Metab.RData")
Data$Metab = Metab; rm(Metab)
## Methyl - code: e
load("data_Methyl.RData")
Data$Methyl = Methyl; rm(Methyl)
# ## aCGH - code: c
# load("data_aCGH.RData")
# Data$aCGH = aCGH; rm(aCGH)

## combine ICA
ICA  = list()
load("IC20(100)_mRNA.RData")
ICA$mRNA = cica
ICA$mRNA$X=NULL
ICA$mRNA$X_num=NULL
rownames(ICA$mRNA$M) = sub("ic[.]","g",rownames(ICA$mRNA$M))

load("IC20(100)_miR.RData")
ICA$miR = cica; rm(cica)
ICA$miR$X=NULL
ICA$miR$X_num=NULL
rownames(ICA$miR$M) = sub("ic[.]","u",rownames(ICA$miR$M))

load("IC20(100)_Prot.RData")
ICA$Prot = cica; rm(cica)
ICA$Prot$X=NULL
ICA$Prot$X_num=NULL
rownames(ICA$Prot$M) = sub("ic[.]","p",rownames(ICA$Prot$M))

load("IC20(100)_Metab.RData")
ICA$Metab = cica; rm(cica)
ICA$Metab$X=NULL
ICA$Metab$X_num=NULL
rownames(ICA$Metab$M) = sub("ic[.]","m",rownames(ICA$Metab$M))

load("IC20(100)_Methyl.RData")
ICA$Methyl = cica; rm(cica)
ICA$Methyl$X=NULL
ICA$Methyl$X_num=NULL
rownames(ICA$Methyl$M) = sub("ic[.]","e",rownames(ICA$Methyl$M))

# load("IC20(100)_aCGH.RData")
# ICA$aCGH = IC; rm(IC)
# rownames(ICA$aCGH$M) = sub("ic[.]","c",rownames(ICA$aCGH$M))

R2=list()
i1=1;i2=4
for (i1 in 1:(length(ICA)-1)){
  for (i2 in i1:length(ICA)){
    nm = c(names(ICA)[i1],names(ICA)[i2])
    nm = sub("aCGH","c",sub("Methyl","e",sub("Metab","m",sub("Prot","p",sub("miR","u",sub("mRNA","g",nm))))))
    nm = sort(nm)
    nm = paste(nm,collapse="")
    print(nm)
    pat = intersect(colnames(ICA[[i1]]$M),colnames(ICA[[i2]]$M))
    if (length(pat)==0) next
    R2[[nm]] = cor(t(ICA[[i1]]$M[,pat]),t(ICA[[i2]]$M[,pat]))^2
  }
}
str(R2)

if(0){
  save(Data,ICA,R2,file="integrate_5mod20comp.RData")
}else{
  load("integrate_5mod20comp.RData")
}

pdf("R2heatmaps_5mod.pdf",onefile = TRUE)
for (i in 1:length(R2)){
  pheatmap(R2[[i]],main=names(R2)[i])
}
dev.off()

plot(density(unlist(R2)))
r2 = unlist(R2)
r2[r2>0.99999]=NA
sum(r2>0.2, na.rm=TRUE)
sum(r2>0.15, na.rm=TRUE)
sum(r2>0.1, na.rm=TRUE)
plot(density(r2,na.rm=TRUE))

which.max(r2)

thr = 0.1

par(mfcol=c(3,2))
for (i in 1:length(R2)){
  plot(density(R2[[i]]),lwd=2,col=4,main=names(R2)[i],xlim=c(0,0.5))
  points(as.vector(R2[[i]]),as.vector(R2[[i]])*0,col=1)
  abline(v=thr,col=2,lty=2)
}

i=1
Links = NULL
for (i in 1:length(R2)){
  link = melt(R2[[i]])
  link[[1]] = as.character(link[[1]])
  link[[2]] = as.character(link[[2]])
  link = link[abs(link[,3])>thr,]
  link = link[link[[1]]!=link[[2]],]
  if (nrow(link) == 0) next
  if (is.null(Links)){
    Links=link
  }else{
    Links=rbind(Links,link)
  }
}
str(Links)
Links = Links[!duplicated(apply(apply(Links[,1:2],1,sort),2,paste,collapse="-")),]
write.table(Links,file = "similarity.txt",sep="\t",quote=F,row.names=F)

pdf(sprintf("graph_R2=%g_5mod20comp.pdf",thr),width=10,height=10,onefile=TRUE)
G = graph_from_edgelist(as.matrix(Links[,1:2]), directed = FALSE)
V(G)$color = c(g = "#AAAAFF", u = "#AAFFAA", p="#FFAAAA", m="#FFFFAA", e="#AAFFFF", c = "#AAAAAA")[substr(V(G)$name,1,1)]
V(G)$size = 10 #+ nn[V(G)$name]*4
V(G)$label.col = "black"
V(G)$label.cex=0.8
V(G)$label.font=2
E(G)$width = (Links$value^2)*10+1
E(G)$color = "#00000055"#ifelse(rall$value>0,"red","blue")
par(mfcol=c(1,1))
plot(G,layout = layout_nicely)
dev.off()

#####################
pat = intersect(colnames(Data$mRNA$X),colnames(Data$Metab$X))
plot(ICA$mRNA$M[16,pat],ICA$Metab$M[8,pat])

#####################
v=Venn(list(mRNA  = colnames(Data$mRNA$X),
            miR   = colnames(Data$miR$X),
            Prot  = colnames(Data$Prot$X),
            Metab = colnames(Data$Metab$X)))

rx = matrix(NA,nrow=20*20,ncol=92)
rx[,] = rnorm(ncol(rx)*nrow(rx))
r2 = cor(t(rx))^2
r2[r2==1]=NA
plot(density(r2,na.rm=TRUE),main="Random (n=92)",lwd=2)
x=sort(as.vector(r2),decreasing=TRUE)[1:1000]
points(x,x*0,main="Random (n=92)",lwd=2)

## Future steps:
## 1) ICA on each "omics" data
## 2) Data integration via ICA
## 3) Age prediction
## 4) Prediction of one data set with another



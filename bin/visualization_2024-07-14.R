## Data and ICA are read in integrate_5mod_2024-03-22.r
rm(list=ls())
library(ggplot2)
library(pheatmap)
library(ggcorrplot)

path = "Z:/users/petr.nazarov/results/figures_2024-07"
dir.create(path)
setwd(path)

if(0){
  save(list=ls(),file="workspace.RData")
}else{
  #load("workspace.RData")
  load("../multiomicsMAE_2024-03-22/integrate_5mod20comp.RData")
  #load("../multiomics_2023-05-13/data_Pheno.RData")
}

####################################################################
## mRNA
str(Data$mRNA)
df = data.frame(scale(t(ICA$mRNA$M)), Data$mRNA$Var,Data$mRNA$Pheno)
str(df)

pdf("mRNA.pdf", onefile = TRUE,  width = 4, height=4)

mod = t.test(g1~Sex,data=df)
tit=sprintf("mRNA ICA linked to Sex, t=%.1f (pval=%.1e)",
            mod$statistic,mod$p.value)
ggplot(df, aes(x=Sex, y=g1, fill=Sex))+
      geom_violin(trim = TRUE) +
      stat_summary(fun.data="mean_sdl", fun.args = list(mult = 1),width=0.05,geom="pointrange",col="black")+
      scale_fill_manual(values=c('#FFAAAA','#AAAAFF'))+
      labs(title=tit,
           x = "Sex",
           y = "Standardized weight (component g1)")


mod = lm(Age.1~g7,data=df)
tit=sprintf("mRNA ICA correlated with Age, R2=%.1f (pval=%.1e)",
                  summary(mod)$r.squared,summary(mod)$coefficients[2,4])
ggplot(df, aes(y=Age.1, x=g7)) +
  geom_point(shape=18, col="black") +
  geom_smooth(method=lm, color="blue", fill="blue") +
  labs(title=tit,
       y="Age (years)",
       x="Standardized weight (component g7)")


mod = lm(Leukocytes.1 ~ g14,data=df)
tit=sprintf("mRNA ICA correlated with Leukocytes, R2=%.2f (pval=%.1e)",
            summary(mod)$r.squared,summary(mod)$coefficients[2,4])
ggplot(df, aes(y=Leukocytes.1, x=g14)) +
  geom_point(shape=18, col="black") +
  geom_smooth(method=lm, color="blue", fill="blue") +
  labs(title=tit,
       y="Leukocytes abundance",
       x="Standardized weight (component g14)")

dev.off()

## get GO for components
ls()
names(ICA$mRNA$GO$GOBP)
str(ICA$mRNA$GO$GOBP[[1]])
ICA$mRNA$GO

extractGO = function(GO,alpha=0.01){
  TabGO = GO[[1]][[1]]$pos[,c(1,2,7)]
  TabGO = TabGO[-(1:nrow(TabGO)),]
  TabGO = data.frame(Component=integer(0),
                     Gene.Subset=character(0),
                     Domain=character(0),
                     TabGO[-(1:nrow(TabGO)),])
  TabGO
  ncomp = length(GO[[1]])
  icomp=1
  for (icomp in 1:ncomp){
    goname="GOBP";d="pos";
    for (goname in c(names(GO))){
      for (d in c("pos","neg")){
        igo = which(GO[[goname]][[icomp]][[d]]$FDR<alpha)
        if (length(igo)>0){
          TabGO = rbind(TabGO,
                  data.frame(Component=icomp,
                           Gene.Subset=d,
                           Domain=sub("GO","",goname),
                           GO[[goname]][[icomp]][[d]][igo,c("GO.ID","Term","FDR")])
                  )
        } ## if some signif
      } ## for pos and neg
    } ## for domains
  } ## for components
  return(TabGO)
}

write.table(extractGO(ICA$mRNA$GO),
            file="mRNA_ICA-GO.tsv",
            sep="\t",row.names=FALSE,quote=FALSE)


###############################################################
## Methyl
str(Data$Methyl)
df = data.frame(scale(t(ICA$Methyl$M)), Data$Methyl$Var,Data$Methyl$Pheno)
str(df)

pdf("Methyl.pdf", onefile = TRUE,  width = 4, height=4)

mod = t.test(e3~Sex,data=df)
tit=sprintf("Methyl ICA linked to Sex, t=%.1f (pval=%.1e)",
            mod$statistic,mod$p.value)
ggplot(df, aes(x=Sex, y=e3, fill=Sex))+
  geom_violin(trim = TRUE) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill="#FFFFFF")+
  #stat_summary(fun.data="mean_sdl",mult=1,width=0.05,geom="crossbar",width=0.2)
  stat_summary(fun.data="mean_sdl", fun.args = list(mult = 1),width=0.05,geom="pointrange",col="black")+
  scale_fill_manual(values=c('#FFAAAA','#AAAAFF'))+
  labs(title=tit,
       x = "Sex",
       y = "Standardized weight (component e3)")


mod = lm(e19~Lymphocytes.1,data=df)
tit=sprintf("(t)Methyl ICA correlated with Lymphocytes, R2=%.2f (pval=%.1e)",
            summary(mod)$r.squared,summary(mod)$coefficients[2,4])
ggplot(df, aes(x=Lymphocytes.1, y=e19)) +
  geom_point(shape=18, col="black") +
  geom_smooth(method=lm, color="blue", fill="blue") +
  labs(title=tit,
       x="Lymphocyte abundance",
       y="Standardized weight (component e19)")


mod = lm(Lymphocytes.1~e19,data=df)
tit=sprintf("Methyl ICA correlated with Lymphocytes, R2=%.2f (pval=%.1e)",
            summary(mod)$r.squared,summary(mod)$coefficients[2,4])
ggplot(df, aes(y=Lymphocytes.1, x=e19)) +
  geom_point(shape=18, col="black") +
  geom_smooth(method=lm, color="blue", fill="blue") +
  labs(title=tit,
       y="Lymphocyte abundance",
       x="Standardized weight (component e19)")


mod = lm(Neutrophils.1~e19,data=df)
tit=sprintf("Methyl ICA correlated with Neutrophils, R2=%.2f (pval=%.1e)",
            summary(mod)$r.squared,summary(mod)$coefficients[2,4])
ggplot(df, aes(y=Neutrophils.1, x=e19)) +
  geom_point(shape=18, col="black") +
  geom_smooth(method=lm, color="blue", fill="blue") +
  labs(title=tit,
       y="Neutrophils abundance",
       x="Standardized weight (component e19)")


mod = lm(Age.1~e17,data=df)
tit=sprintf("Methyl ICA correlated with Age, R2=%.1f (pval=%.1e)",
            summary(mod)$r.squared,summary(mod)$coefficients[2,4])
ggplot(df, aes(y=Age.1, x=e17)) +
  geom_point(shape=18, col="black") +
  geom_smooth(method=lm, color="blue", fill="blue") +
  labs(title=tit,
       y="Age (years)",
       x="Standardized weight (component e17)")


mod = lm(Age.1~e17+e4,data=df)
x = predict(mod,df)
tit=sprintf("Methyl ICA correlated with Age, R2=%.1f (pval=%.1e, %.1e)",
            summary(mod)$r.squared,summary(mod)$coefficients[2,4],summary(mod)$coefficients[3,4])
ggplot(df, aes(y=Age.1, x=x)) +
  geom_point(shape=18, col="black") +
  geom_smooth(method=lm, color="blue", fill="blue") +
  labs(title=tit,
       y="Age (years)",
       x="Standardized weight (component e17 + e4)")


dev.off()


##########################################################################

Ph = list()
M = list()

for (idata in 1:length(Data)){
  M[[names(Data)[idata]]] = scale(t(ICA[[idata]]$M))
  Ph[[names(Data)[idata]]] = Data[[idata]]$Pheno
  
  for (i in 1:ncol(Ph[[idata]])){
    if (class(Ph[[idata]][[i]])[1]%in%c("character","ordered","factor")){
      Ph[[idata]][[i]] = as.integer(factor(Ph[[idata]][[i]]))
    }
  }
  Ph[[idata]]=Ph[[idata]][,-c(3:5,9,ncol(Ph[[idata]])-1,ncol(Ph[[idata]]))]
}


str(Data$mRNA$Pheno)

pdf("heatmaps_pearson.pdf", onefile = TRUE)

for (idata in 1:length(Data)){
  pheatmap(cor(M[[idata]], Ph[[idata]], method="pearson", use="pairwise.complete.obs")^2, 
         main=sprintf("R2 (Pearson) %s", names(Data)[idata]), 
         cluster_cols = FALSE)
}

dev.off()


pdf("heatmaps_spearman.pdf", onefile = TRUE)

for (idata in 1:length(Data)){
  pheatmap(cor(M[[idata]], Ph[[idata]], method="spearman", use="pairwise.complete.obs")^2, 
           main=sprintf("R2 (Spearman) %s", names(Data)[idata]), 
           cluster_cols = FALSE)
}

dev.off()

###############################################################################
## regression for phenotype
## keep only common patients 

pat = summary(factor(unlist(lapply(Ph,rownames))),maxsum = Inf)
pat = names(pat)[pat==5]

pdf("bestcomp1.pdf", onefile = TRUE,  width = 10, height=10)
par(mfcol=c(4,4))
iph = 1
names(Ph[[1]])
for (iph in 1:ncol(Ph[[1]])){
  # select the best predictor for each omics
  Y = Ph[[1]][pat,iph]
  nmY=names(Ph[[1]])[iph]
  YX = data.frame(Y=Y)
  idata=1
  for (idata in 1:length(Data)){
    imax = which.max(cor(Y,M[[idata]][pat,],use="pairwise.complete.obs")^2)
    YX[[names(Data)[idata]]] = M[[idata]][pat,imax]
  }

  str(YX)
  res = c(
    mRNA = summary(lm(Y~mRNA,data=YX))$r.squared,
    miRNA = summary(lm(Y~miR,data=YX))$r.squared,
    Protein= summary(lm(Y~Prot,data=YX))$r.squared,
    Metabol. = summary(lm(Y~Metab,data=YX))$r.squared,
    Methyl. = summary(lm(Y~Methyl,data=YX))$r.squared,
   #aCGH = summary(lm(Y~aCGH,data=YX))$r.squared,
    Multiomics = summary(lm(Y~.,data=YX))$r.squared)
  
  barplot(res, las=2, 
        col=c("#AAAAFF", "#AAFFAA","#FFAAAA","#FFFFAA","#AAFFFF", "#AAAAAA"),
        main=sprintf("R2: Best vs %s",nmY), 
        cex.main=1, cex.axis =1, cex.names=1)

}
dev.off()


pdf("bestcomp2.pdf", onefile = TRUE,  width = 7, height=4)
par(mfrow=c(2,4))
iph = 1
names(Ph[[1]])

for (iph in c(1,2,5,9,13,14,16,8) ){
  # select the best predictor for each omics
  #Y = Ph[[1]][pat,iph]
  Y = Ph$Metab[pat,iph]
  ## if factor (1,2) -> (0,1) to be consistant with point biserial correlation
  if (all(Y%in%c(1,2))) Y = Y-1 
  nmY=names(Ph[[1]])[iph]
  YX = data.frame(Y=Y)
  idata=1
  for (idata in 1:length(Data)){
    imax = which.max(cor(Y,M[[idata]][pat,],use="pairwise.complete.obs")^2)
    YX[[names(Data)[idata]]] = M[[idata]][pat,imax]
  }
  
  str(YX)
  res = c(
    "mRNA-seq" = summary(lm(Y~mRNA,data=YX))$r.squared,
    "miRNA-seq" = summary(lm(Y~miR,data=YX))$r.squared,
    "Proteomics" = summary(lm(Y~Prot,data=YX))$r.squared,
    "Metabolomics"  = summary(lm(Y~Metab,data=YX))$r.squared,
    "EM-seq" = summary(lm(Y~Methyl,data=YX))$r.squared,
    "Multiomics" = summary(lm(Y~.,data=YX))$r.squared)
  r2adj = summary(lm(Y~.,data=YX))$adj.r.squared
  
  barplot(res, las=2, 
          col=c("#AAAAFF", "#AAFFAA","#FFAAAA","#FFFFAA","#AAFFFF", "#AAAAAA"),
          main=sprintf("%s",nmY), 
          ylab="R2 (best comp.)",
          cex.main=1.2, cex.axis =1, cex.names=1)
  #points(6*1.2-0.5,r2adj,pch="+")
  lines(6*1.2-0.5+c(-0.5,0.5),c(r2adj,r2adj), pch="+",col=2,lwd=2)

}
dev.off()

####################################################################
## reproduce MOFA plot with best ICA
#install.packages("ggcorrplot")
library(ggcorrplot)

pdf("corrplots.pdf", onefile = TRUE,  width = 12, height=8)

par(mfrow=c(1,1))
iph = 1
idx =  c(1:16)
names(Ph$Metab)[idx]
## final correlation and PV tables (most correlated) 
R = matrix(nrow=5,ncol=length(idx))
rownames(R) = names(Data)
colnames(R) = names(Ph[[1]])[idx]
PV = R
IC = R

## all correlations and PV tables in a list (for each omics) 
lst_r = list()
r = matrix(nrow=20,ncol=length(idx)) # temp variable to be put to list
rownames(r) = 1:20
colnames(r) = names(Ph[[1]])[idx]
lst_pv=lst_r

for (idata in 1:length(Data)){
  rownames(r) = colnames(M[[names(Data)[idata]]])
  lst_r[[names(Data)[idata]]] = r
  lst_pv[[names(Data)[idata]]] = r
}

## loop for pheno variables 
for (iph in idx){
  ## get variable (metab is the most complete)
  Y = Ph$Metab[pat,iph]
  ## if factor (1,2) -> (0,1) to be consistent with point biserial correlation
  if (all(Y%in%c(1,2))) Y = Y-1 
  nmY=names(Ph[[1]])[iph]
  YX = data.frame(Y=Y)
  idata=1
  # cal correlations and select the best predictor for each omics 
  for (idata in 1:length(Data)){
    imax = which.max(cor(Y,M[[idata]][pat,],use="pairwise.complete.obs")^2)
    YX[[names(Data)[idata]]] = M[[idata]][pat,imax]
    IC[idata,nmY] = paste(names(Data)[idata],imax)
    lst_r[[names(Data)[idata]]][,iph] = cor(Y,M[[idata]][pat,],use="pairwise.complete.obs")
    lst_pv[[names(Data)[idata]]][,iph] = cor_pmat(cbind(Y,M[[idata]][pat,]))[1,-1]
  }
  R[,nmY] = cor(YX[,1],YX[,-1],use="pairwise.complete.obs")
  PV[,nmY] = cor_pmat(YX)[1,-1] 
}
IC = sub("mRNA ","g",IC)
IC = sub("miR ","u",IC)
IC = sub("Metab ","m",IC)
IC = sub("Methyl ","e",IC)
IC = sub("Prot ","p",IC)
R = t(R)
PV = t(PV)
IC = t(IC)
for (idata in 1:length(Data)){
  lst_r[[idata]] = t(lst_r[[idata]])  
  lst_pv[[idata]] = t(lst_pv[[idata]])  
}  

#ggcorrplot(R, p.mat=PV, title="IC components in each omics, most correlated with variables")

for (idata in 1:length(Data)){
  print(
    ggcorrplot(lst_r[[idata]], p.mat=lst_pv[[idata]], lab=F, 
            title=sprintf("Correlation of component weights in %s",names(Data)[idata]))
  )  
}
dev.off()

pdf("corrplot_bestcomp.pdf", onefile = TRUE,  width = 8, height=5)

## rename omics and reorder
colnames(R) = sub("mRNA","mRNA-seq",colnames(R))
colnames(R) = sub("miR","miRNA-seq",colnames(R))
colnames(R) = sub("Prot","Proteomics",colnames(R))
colnames(R) = sub("Metab","Metabolomics",colnames(R))
colnames(R) = sub("Methyl","EM-seq",colnames(R))
colnames(PV) = colnames(R)
iord = c(4,3,2,1,5)
R = R[,iord]
PV = PV[,iord]
IC = IC[,iord]

## Do correction within vizualized 
FDR = PV
FDR[,] = p.adjust(PV,"fdr")

cor_plot = ggcorrplot(R, p.mat=PV, lab=F, title="Components in each omics dataset, most correlated with phenotype")

print(
cor_plot+
  geom_text(aes(x=Var1,
                y=Var2),
          label=as.vector(IC)
          )
)

cor_plot = ggcorrplot(R, p.mat=FDR, lab=F, title="Components in each omics dataset, most correlated with phenotype")

print(
  cor_plot+
    geom_text(aes(x=Var1,
                  y=Var2),
              label=as.vector(IC)
    )
  )

dev.off()

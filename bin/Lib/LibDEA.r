###############################################################################
## Collection of various warp-ups for DEA methods: RNAseq and microarrays
## (c) GNU GPL Petr Nazarov, Luxembourg Institute of Health, petr.nazarov[at]lih.lu
## last modification 2017-02-14
###############################################################################

## install packages required
{
print("----------------------------------------------------------------")
print("libDEA: Differential expression analysis")
packages = c("limma","edgeR","DESeq2","caTools","compiler","sva")
print("Needed:")
print(packages)
print("Absent:")
print(packages[!packages %in% rownames(installed.packages())])

if (!requireNamespace("BiocManager", quietly = TRUE)) 
	install.packages("BiocManager")

for (pk in packages){
	if (! pk %in% rownames(installed.packages()))
		BiocManager::install(pk)
}
rm(packages,pk)
}

## ToDO:
## add combat correction 
## = ComBat(Data$LXN, batch = as.integer(Data$meta$Replicate%in%c("M2","M3")))


library(compiler)
## I can use dir.create() !!!
requireFolder = function(folder=NULL,relative=TRUE){
	if (is.null(folder)) return(FALSE)
	if (folder=="") return(FALSE)
	if (!file.exists(folder)) {
		cat(sprintf("Folder '%s' is not found. Creating...\n",folder))
		if (Sys.info()[1] == "Windows"){
			if (relative){
				#try(shell(paste("mkdir",gsub("/","\\\\",file.path(getwd(),folder)))))
				try(shell(paste("mkdir \"",gsub("/","\\\\",file.path(getwd(),folder)),"\"",sep="")))
				## ToDo: try here parameter tranaslate=TRUE for slash->backslash
			}else{
				try(shell(paste("mkdir",gsub("/","\\\\",folder))))
			}
		}
		if (Sys.info()[1] == "Linux")
			try(system(  paste("mkdir",file.path(getwd(),folder)) ))    
		if (!file.exists(folder))
				stop("Cannot create folder to store dowloaded data!\n")
	}
	#return(TRUE)
}

getTopIdx=function(x,n){
	return(order(x,na.last=TRUE,decreasing=TRUE)[1:n])
}
getConcord = function(list1,list2){
	return( (length(intersect(list1,list2))/length(list1) + length(intersect(list1,list2))/length(list2))/2)
}
getJaccard = function(list1,list2){
	return( length(intersect(list1,list2))/(length(union(list1,list2))))
}
getMeanCI = function(x,alpha=0.05,do.print=FALSE){
	n=length(x[!is.na(x)])
	out = list()
	out$m = mean(x,na.rm=TRUE)
	out$me = -qt(alpha/2,n-1)*sd(x,na.rm=TRUE)/sqrt(n)
	
	if (do.print) print(sprintf("Mean(x[%d]) = %g +/- %g (CI %g%%)",n,out$m,out$me,100*(1-alpha)))
	return(out)
}

##=============================================================================
## calc.AUC - calculate Area Under roc-Curve 
##	x - data vector
##	group - classes vector
##	key0 - id for class 0 in group
##	key1 - id for class 1 in group
##=============================================================================
calc.AUC=function(x,group,key0=NULL,key1=NULL,fun=NULL){
	require("caTools")
	nkey = nlevels(factor(as.character(group)))
	lvs = levels(factor(as.character(group)))
	if (nkey == 2 & (is.null(key0) | is.null(key1))) {
		key0 = lvs[1]
		key1 = lvs[2]
	}
	if (is.null(key0) | is.null(key1)){
		## compare all combinations
		if (class(x)[1] == "numeric"){
			auc = double(nrow = nkey*(nkey-1)/2)
		}else{
			auc = matrix(nrow = nkey*(nkey-1)/2, ncol = nrow(x))
			rownames(auc) =1:nrow(auc)
			colnames(auc) =rownames(x)
			k = 1
			for (i in 1:(nkey-1)){
				for (j in (i+1):nkey) {
					rownames(auc)[k] = paste(lvs[i],"vs.",lvs[j])
					auc[k,] = calc.AUC(x,group,key0 = lvs[i],key1=lvs[j])
					k=k+1
				}
			}
		}
	} else {
		## simple case - comparing 2 groups
		idx0 = which(group == key0)
		idx1 = which(group == key1)
		group = group[c(idx0,idx1)]
		if (class(x)[1] =="numeric" || class(x)[1] =="integer"){
			auc = colAUC(t(x[c(idx0,idx1)]), factor(as.character(group)), plotROC=FALSE)[1,1]
		}else{
			auc = colAUC(t(x[,c(idx0,idx1)]), factor(as.character(group)), plotROC=FALSE)
		}
	}
	if (!is.null(fun)) auc = apply(auc,2,fun)
	return(auc)
}
##=============================================================================
## Normalization
## method = "DESeq","edgeR","voom"
##=============================================================================
if (FALSE){
	method="DESeq"
	doLog=TRUE
	group=NULL
}
norm.Counts=function(count, group=NULL, method="DESeq", doLog = FALSE){
	## ToDo: add voom!
	
	if (!method %in% c("DESeq","edgeR","voom")) stop("in norm.Counts {LibDEA}: unknown method",method,"\n")

	if (is.null(group)) group = rep("one",ncol(count))
	group = factor(group)
	
	if (method =="voom"){
	  library(limma)
		if(!doLog) print("Voom automatically does log2 transformation!")
		x = voom(counts = count)$E
	}
	if (method =="edgeR"){
	  require(edgeR)
		dge = DGEList(counts=round(count),group=group)
		dge = calcNormFactors(dge,method = "TMM")
		x = cpm(dge, normalized.lib.sizes=TRUE)
		if (doLog) x = log2(0.5 + x)
	}
	if (method =="DESeq"){
	  require(DESeq2)
		colData = data.frame(condition=group)
		if (nlevels(group)==1) {
			des = DESeqDataSetFromMatrix(round(count),colData, formula(~1))
		}else{
			des = DESeqDataSetFromMatrix(round(count),colData,formula(~ condition))
		}
		#des = estimateSizeFactors(des,type="iterate")
		des = estimateSizeFactors(des)
		x = counts(des,normalize=T)
		if (doLog) x = log2(1 + x)
		## ToDo: try VST (check `getNormalizedData`)
	}
	return(x)
}

Count2FPKM = function(X, len){
	FPKM = X * 0
	for (i in 1:ncol(X)) FPKM[,i] = X[,i] / len / sum(X[,i]) * 1e9
	return(FPKM)
}

Count2TPM= function(X, len){
	TPM = X * 0
	for (i in 1:ncol(X)){
		TPM[,i] = X[,i] / len
		TPM[,i] = TPM[,i] /sum(TPM[,i]) * 1e6
	}
	return(TPM)
}
##=============================================================================
## DESeq
##=============================================================================
if (FALSE){
	count = NGS$counts
	group = paste(meta$Type,meta$State,sep=".")
	key0 = paste(type,"Normal",sep=".")
	key1 = paste(type,"Tumor",sep=".")
	name="zzz"
	folder="DEA"
	pair=NULL
}
DEA.DESeq = function(count,group,key0,key1,pair=NULL,name=NULL,folder="",return.auc=FALSE,adjust = "BH"){
	require(DESeq2)
	requireFolder(folder)
	if (ncol(count)!=length(group)) 
		stop("Length of 'group' should corresponds to the number of data columns")
	if (!is.null(pair))
		if (ncol(count)!=length(pair)) 
			stop("Length of 'pair' should corresponds to the number of data columns")
	idx0 = which(key0==group)
	idx1 = which(key1==group)
	group = factor(group[c(idx0,idx1)])
	group = relevel(group,key0)
	if (!is.null(pair)) {
		pat = factor(pair[c(idx0,idx1)])
		colData = data.frame(condition=group,pair=pat)
		des = formula(~ condition + pair)
	}else{
		colData = data.frame(condition=group)
		des = formula(~ condition)
	}
	des = DESeqDataSetFromMatrix(count[,c(idx0,idx1)],colData,design = des)
	des = DESeq(des)
	resultsNames(des)
	#v. <1.3 - prev.version!!!!!!! check versions!
	#res = results(des,sprintf("condition_%s_vs_%s",levels(group)[2],levels(group)[1]),pAdjustMethod = "BH")
	#v. >=1.3
	res = results(des,contrast=c("condition",levels(group)[2],levels(group)[1]),pAdjustMethod = adjust)
	#str(res)
	Tab = data.frame(id=res@rownames,res@listData,stringsAsFactors=F)
	names(Tab)=sub("pvalue","PValue",sub("log2FoldChange","logFC",sub("padj","FDR",names(Tab))))
	rownames(Tab)=Tab[,1]
	Tab$FDR[is.na(Tab$FDR)]=1
	res = Tab
	cat(sprintf("\nDESeq,%s%s: %d DEG (FDR<0.05), %d DEG (FDR<0.01)\n\n",
					ifelse(is.null(pair),"unpaired","paired"),
					ifelse(is.null(name),"",paste(" on",name)),
					sum(res$FDR<0.05),sum(res$FDR<0.01)))
	if (return.auc){
		if (require(caTools)){
			res$AUC = calc.AUC(x = count[,c(idx0,idx1)],group = group, key0,key1)[1,]
			## count[,..] not all count as gtoup is reassigned
		}
	}
	if (!is.null(name)){
		if (!is.null(pair)){
			file.name = file.path(folder,sprintf("DEA(DESeq,paired)=%s_%s-%s.txt",name,key1,key0))
		}else{
			file.name = file.path(folder,sprintf("DEA(DESeq)=%s_%s-%s.txt",name,key1,key0))
		}
		write.table(res,file=file.name,sep="\t",quote=F,row.names=F)
	}
	
	res$logFC[is.na(res$logFC)] = 0
	return(res)
}


##=============================================================================
## edgeR
##=============================================================================
if (FALSE){

#	ResP = DEA.edgeR(X,group=group,key0=key0,key1=key1,pair=xmeta$Batch,name=NULL,folder="DEA")
	count = X
	group = group
	key0 = key0
	key1 = key1
	pair=xmeta$Batch
	prior.df=10
	name=NULL
	folder=""
	return.auc=FALSE
	return.model=TRUE
	prefer.GLM=TRUE
}
##--------------------------------
DEA.edgeR=function(count,group,key0,key1,pair=NULL,prior.df=10,name=NULL,folder="",return.auc=FALSE,return.model=FALSE, prefer.GLM=FALSE,adjust="BH"){
	require(edgeR)
	requireFolder(folder)
	idx0 = which(group == key0)
	idx1 = which(group == key1)
	group = factor(group[c(idx0,idx1)])
	group = relevel(group,key0)
	dge = DGEList(counts=count[,c(idx0,idx1)],group=group)
	if (!is.null(pair) | prefer.GLM){ ## use GLM
		if (!is.null(pair)){
			pair = factor(pair[c(idx0,idx1)])
			design = model.matrix(~pair+group)  ## important factor order for glmLRT - (1)pat (2)group
		}else{
			design = model.matrix(~group)
		}
		dge = calcNormFactors(dge)
		dge = estimateGLMCommonDisp(dge,design)
		dge = estimateGLMTrendedDisp(dge,design)
		dge = estimateGLMTagwiseDisp(dge,design)
		fit = glmFit(dge, design)
		lrt = glmLRT(fit)
		res = topTags(lrt,n=nrow(count),adjust=adjust,sort.by="none")[[1]]
		res = data.frame(id=rownames(res),res,stringsAsFactors=F)
	}else{ ## use exactTest
		dge = calcNormFactors(dge)
		dge = estimateCommonDisp(dge)
		dge = estimateTagwiseDisp(dge,prior.df=prior.df)
		res= exactTest(dge)[[1]]
		res$FDR=p.adjust(res$PValue,"fdr")
		res = data.frame(id=rownames(res),res,stringsAsFactors=F)
	}
	cat(sprintf("\nedgeR,%s%s: %d DEG (FDR<0.05), %d DEG (FDR<0.01)\n\n",
					ifelse(is.null(pair),"unpaired","paired"),
					ifelse(is.null(name),"",paste(" on",name)),
					sum(res$FDR<0.05),sum(res$FDR<0.01)))	
	
	if (return.auc){
		if (require(caTools)){
			res$AUC = calc.AUC(x = count[,c(idx0,idx1)],group = group, key0,key1)[1,]
			## count[,..] not all count as group is reassigned
		}
	}
	
	if (!is.null(name)){
		if (!is.null(pair)){
			file.name = file.path(folder,sprintf("DEA(edgeR,paired)=%s_%s-%s.txt",name,key1,key0))
		}else{
			file.name = file.path(folder,sprintf("DEA(edgeR)=%s_%s-%s.txt",name,key1,key0))
		}
		write.table(res,file=file.name,sep="\t",quote=F,row.names=F)
	}
	
	if (return.model){
		TMP = res
		res=list()
		res$name = ifelse(is.null(pair),"unpaired","paired")
		if (exists("fit")){
			res$dge = dge
			res$lmFit = fit
			res$lmLRT = lrt
		}else{
			res$model = "exactTest was used. No model."
		}
		res$topTable = TMP
	}
	
	return(res)
}

##=============================================================================
## LIMMA: simple 2- or n-group comparison
## Batch correction is possible
##=============================================================================
DEA.limma=function(	data=NULL,		## matrix with genes in rows
					group,			## vector of conditions
					key0=NULL,		## for contrast: control condition
					key1=NULL,		## for contrast: experimental condition
					#ToDo:use=c("specific","all")[1],## is model based on all data or only specific key0, key1,
					pair=NULL,		## vector of pairs for paired analysis
					name=NULL,		## name for saving results
					folder="", 		## folder for saving results
					counted=FALSE,	## does data contain counts? 
					norm="none",	## normalizations: none, scale, quantile
					return.model=FALSE, ## should linear model be returned?
					return.auc=FALSE, ## should AUC be calculated?
					adjust="BH",	## adjustment method. Default - Benjamini-Hochberg
					silent = FALSE,	## avoid reporting messages
					anno=NULL		## if additional annotation should be given
					){
	require(limma)
	requireFolder(folder)
	if (is.null(rownames(data))) rownames(data) = sprintf("row%05d",1:nrow(data))
	
	keep.treatment = names(summary(factor(group))[summary(factor(group))>1])
	data = data[,group %in% keep.treatment]
	group = group[group %in% keep.treatment]
	
	
	## 2 groups (t-test) or n groups (F-test) ?
	if (is.null(key0)){ 
		## F-test
		group=factor(as.character(group))
		design=model.matrix(~ 0 + group)
		colnames(design)=sub("group","",colnames(design))
		comp.list=NULL;k=0
		for (i in 2: nlevels(group)){
			for (j in 1:(i-1) ){
				k=k+1
				comp.list[k] = sprintf("%s-%s",levels(group)[i],levels(group)[j])
			}
		}
		print(comp.list)
		fit = lmFit(data,design=design)
		contr = makeContrasts(contrasts = comp.list,levels=levels(group))
		c.fit = contrasts.fit(fit, contr)
		EB = eBayes(c.fit)
		res = data.frame(id=rownames(data),topTable(EB, number=nrow(data),adjust=adjust,sort.by="none"),stringsAsFactors=F) 
	}else{
		idx0 = which(group == key0)
		idx1 = which(group == key1)
		group = factor(group[c(idx0,idx1)])
		group = relevel(group,key0)
		## Paired (with batch) or unpaired?
		if (!is.null(pair)){ 
			## paired
			pair = factor(pair[c(idx0,idx1)])
			design = model.matrix(~pair+group)  ## important factor order for glmLRT - (1)pat (2)group
			x = data[,c(idx0,idx1)]
			if (counted) {
				x = voom(x,design=design,plot=FALSE,normalize.method=norm)
			}else{
				if (norm=="scale") x = scale(x)
				if (norm=="quantile") x = normalizeQuantiles(x, ties=TRUE)
			}
			fit = lmFit(x,design=design)
			EB = eBayes(fit)
			res = data.frame(id=rownames(data),topTable(EB,coef=ncol(design),number=nrow(data),adjust=adjust,sort.by="none"),stringsAsFactors=F)
		}else{ 
			## unpaired
			design=model.matrix(~ group)
			x = data[,c(idx0,idx1)]
			if (counted) {
				x = voom(x,design=design,plot=FALSE,normalize.method=norm)
			}else{
				if (norm=="scale") x = scale(x)
				if (norm=="quantile") x = normalizeQuantiles(x, ties=TRUE)
			}
			fit = lmFit(x,design=design)
			EB = eBayes(fit)
			res = data.frame(id=rownames(data),topTable(EB,coef=2,number=nrow(data),adjust=adjust,sort.by="none"),stringsAsFactors=FALSE)
		} ## if (!is.null(pair))
	} ## if (is.null(key0))
	
	names(res)=sub("P.Value","PValue",sub("adj.P.Val","FDR",names(res)))
	res$FDR[is.na(res$FDR)] = 1
	
	if (!silent)
	cat(sprintf("\nLimma,%s%s: %d,%d,%d,%d DEG (FDR<0.05, 1e-2, 1e-3, 1e-4)\n\n",
					ifelse(is.null(pair),"unpaired","paired"),
					ifelse(is.null(name),"",paste(" on",name)),
					sum(res$FDR<0.05),sum(res$FDR<0.01),sum(res$FDR<0.001),sum(res$FDR<0.0001)))
	
	if (return.auc){
		if (require(caTools)){
			if (is.null(key0)){
				res$AUC = calc.AUC(x = data,group = group, fun="mean")
			}else{
				res$AUC = calc.AUC(x = data[,c(idx0,idx1)],group = group, key0,key1,fun="mean")
			}
		}
	}
	if (!is.null(name)){
		if (!is.null(pair)) file.name = file.path(folder,sprintf("DEA(limma,paired)=%s_%s-%s.txt",name,key1,key0))
		if (is.null(pair) & !is.null(key0)) file.name =file.path(folder,sprintf("DEA(limma)=%s_%s-%s.txt",name,key1,key0))
		if (is.null(pair) & is.null(key0)) file.name =file.path(folder,sprintf("DEA(limma,Fstat)=%s.txt",name))
		write.table(res,file=file.name,sep="\t",quote=F,row.names=F)
	}
	
	if (!is.null(anno)){
		res= data.frame(id = res$id,anno,res[,-1],stringsAsFactors=FALSE)
	}
	
	if (return.model){
		TMP = res
		res=list()
		res$name = ifelse(is.null(pair),"unpaired","paired")
		res$lmFit = fit
		res$eBayes = EB
		res$topTable = TMP
	}
		
	return(res)
}

getPHeatmapFont = function(n) return(ifelse (n < 40/2,20,400/n))

plotVolcano=function(	res,
						thr.fdr=0.05,
						thr.lfc=1,
						pch=19,
						col=c(all="#AAAAAA",pos="#FF0000",neg="#0000FF"),
						cex=1,
						main="Volcano",
						usePValue=FALSE,
						add.names=c("none","singificant","all")[1],
						cex.names=0.8,
						col.names="black"){
	if (usePValue) {
		fdr = c(res$PValue,res$P.Value)
	}else{
		fdr = res$FDR
	}
	lfc = res$logFC
	ipos = (fdr<thr.fdr)&(lfc > thr.lfc)
	ineg = (fdr<thr.fdr)&(lfc < -thr.lfc)
	inon = !(ipos|ineg)
	plot(lfc[inon],-log10(fdr[inon]),pch=pch,col=col[1],cex=cex,
		main=main,xlim=range(lfc),ylim=c(0,max(-log10(fdr))),
		xlab="log2FC",ylab=ifelse(usePValue,"-log10(p-value)","-log10(FDR)"))
	abline(v=c(-thr.lfc,thr.lfc),lty=2,col="grey")
	abline(h=-log10(thr.fdr),lty=2,col="grey")
	if (sum(ipos)>0)	points(lfc[ipos],-log10(fdr[ipos]),pch=pch,col=col[2],cex=cex)
	if (sum(ineg)>0)	points(lfc[ineg],-log10(fdr[ineg]),pch=pch,col=col[3],cex=cex)
	
	if (add.names == "significant"){
		if (sum(ipos)>0)  text(lfc[ipos],-log10(fdr[ipos]),res[[1]][ipos],col=col.names,cex=cex.names,adj=c(0.5,0))
		if (sum(ineg)>0) text(lfc[ineg],-log10(fdr[ineg]),res[[1]][ineg],col=col.names,cex=cex.names,adj=c(0.5,0))
	}
	if (add.names == "all"){
		text(lfc,-log10(fdr),res[[1]],col=col.names,cex=cex.names,adj=c(0.5,0))
	}

}
#plotVolcano(Res[[2]],usePValue=TRUE,thr.lfc=0.5)

plotDotBox = function(x,group,col=NULL, border="#888888", pt.col=NULL,pt.cex=1.5,main="Dot-Box plot",ylab="value",xlab="groups",cex.axis=1.5,cex.lab=1.5){
	group = factor(group)
	gr = levels(group)[1]
	if (is.null(col)) 		col = rainbow(nlevels(group),alpha=0.3)
	if (length(col)==1) 	col=rep(col,nlevels(group))
	if (is.null(names(col)))names(col)=levels(group)
	
	if (is.null(pt.col)) 	pt.col = col
	if (length(pt.col)==1)	pt.col=rep(pt.col,nlevels(group))
	if (is.null(names(pt.col)))names(pt.col)=levels(group)
	
	boxplot(x~group,col=col,border=border,las=1,cex.axis=cex.axis,cex.lab=cex.lab,xlab=xlab,ylab=ylab,main=main,outline = FALSE)
	for (gr in levels(group)){
		idx = group == gr
		v = x[idx]
		h = rep(which(levels(group)==gr),sum(idx)) + runif(sum(idx),-0.25,0.25)#(1:length(v)-length(v)/2-0.5)/80
		points(h,v,col=pt.col[gr],pch=19,cex=pt.cex)
		points(h,v,col=1,pch=1,cex=pt.cex)
	}
}

##########################################################################
## GSEA

#members = sig$ccGO
#score = apply(Data$LX[sig$all,Var$Subgroup=="GBM_G34"],1,mean) - apply(Data$LX[sig$all,Var$Subgroup=="GBM_K27"],1,mean)
#show= TRUE
runGSEA = function(score, genes=NULL, members, show=FALSE,nperm=10000){
	if (is.null(genes)) genes = names(score)
	members = members[members %in% genes]
	if (! length(genes)>10 & length(genes)==length(score) & length(members)>1) stop("Wrong lenghts of vectors")
	idx = order(score,decreasing=TRUE)
	score=score[idx]
	genes=genes[idx]
	ng = length(genes)
	
	if (show){
		par(fig=c(0,1,0.75,1))
		plot(NA,NA,xlim=c(1,ng),ylim=c(0,1),main="Genes ordered by score",xlab="",ylab="",yaxt="n",xaxt="n",cex.main=1.2)
		abline(v=which(genes %in% members),col="#00000088")
		text(1,0,"up-regulated",col=2,font=2,adj=c(0,-1))
		text(ng,0,"down-regulated",col=4,font=2,adj=c(1,-1))
	}
	S = double(ng) + ifelse(genes %in% members, 1/length(members),-1/(ng-length(members)))
	ES = max(abs(cumsum(S)))
	
	## permutations
	ES0 = double(nperm)
	for (i in 1:length(ES0)){
		S0 = sample(S)
		ES0[i] = max(abs(cumsum(S0)))
	}
	pv = sum(ES<ES0)/nperm
	
	if (show){
		par(fig=c(0,1,0.4,0.8),new=TRUE)
		plot(cumsum(S),type="l",col=3,lwd=2,xlab="Position in the ranked gene list",ylab="S(g)",
			main = sprintf("\nES = %.2f vs %d permutations: p-value = %.2e",ES,nperm,pv,cex.main=1.2))
		abline(h=0,lty=2)
		
		par(fig=c(0,1,0.0,0.4),new=TRUE)
		plot(density(ES0),lwd=2,col=4,xlim=range(c(ES0,ES)),main="Observed and Permutation ES")
		abline(v=ES,lwd=2,col=2)
		legend("bottomleft",c("observed","permutation"),col=c(2,4),lwd=2,bg="white")
	}
	list(ES=ES, meanES0=mean(ES0), pvalue=pv, nperm=nperm)
}






###########################################################################
## TESTING
###########################################################################

if (FALSE){

key0="hNSC100.ctrl"
key1="hNSC100.TMZ"
count=NGS$counts
logx=HTA$expr
group=paste(meta$Cells,meta$Treatment,sep=".")
pair=NULL
name=NULL
rm(list=c("key0","key1","count","logx","group","pair","name"))

res = DEA.DESeq(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"hNSC100.ctrl","hNSC100.TMZ",name="hNSC100")
 res = DEA.DESeq(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"hNSC100.ctrl","hNSC100.TMZ",name="hNSC100",pair=meta$Replicate)
res = DEA.DESeq(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"NCH644.ctrl","NCH644.TMZ",name="NCH644")
 res = DEA.DESeq(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"NCH644.ctrl","NCH644.TMZ",name="NCH644",pair=meta$Replicate)
res = DEA.DESeq(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"NCH421k.ctrl","NCH421k.TMZ",name="NCH421k")
 res = DEA.DESeq(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"NCH421k.ctrl","NCH421k.TMZ",name="NCH421k",pair=meta$Replicate)
res = DEA.DESeq(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"NCH601.ctrl","NCH601.TMZ",name="NCH601")
 res = DEA.DESeq(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"NCH601.ctrl","NCH601.TMZ",name="NCH601",pair=meta$Replicate)

res = DEA.edgeR(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"hNSC100.ctrl","hNSC100.TMZ",name="hNSC100")
 res = DEA.edgeR(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"hNSC100.ctrl","hNSC100.TMZ",name="hNSC100",pair=meta$Replicate)
res = DEA.edgeR(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"NCH644.ctrl","NCH644.TMZ",name="NCH644")
 res = DEA.edgeR(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"NCH644.ctrl","NCH644.TMZ",name="NCH644",pair=meta$Replicate)
res = DEA.edgeR(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"NCH421k.ctrl","NCH421k.TMZ",name="NCH421k")
 res = DEA.edgeR(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"NCH421k.ctrl","NCH421k.TMZ",name="NCH421k",pair=meta$Replicate)
res = DEA.edgeR(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"NCH601.ctrl","NCH601.TMZ",name="NCH601")
 res = DEA.edgeR(NGS$counts,paste(meta$Cells,meta$Treatment,sep="."),"NCH601.ctrl","NCH601.TMZ",name="NCH601",pair=meta$Replicate)

res = DEA.limma(HTA$expr,paste(meta$Cells,meta$Treatment,sep="."),"hNSC100.ctrl","hNSC100.TMZ",name="hNSC100")
 res = DEA.limma(HTA$expr,paste(meta$Cells,meta$Treatment,sep="."),"hNSC100.ctrl","hNSC100.TMZ",name="hNSC100",pair=meta$Replicate)
res = DEA.limma(HTA$expr,paste(meta$Cells,meta$Treatment,sep="."),"NCH644.ctrl","NCH644.TMZ",name="NCH644")
 res = DEA.limma(HTA$expr,paste(meta$Cells,meta$Treatment,sep="."),"NCH644.ctrl","NCH644.TMZ",name="NCH644",pair=meta$Replicate)
res = DEA.limma(HTA$expr,paste(meta$Cells,meta$Treatment,sep="."),"NCH421k.ctrl","NCH421k.TMZ",name="NCH421k")
 res = DEA.limma(HTA$expr,paste(meta$Cells,meta$Treatment,sep="."),"NCH421k.ctrl","NCH421k.TMZ",name="NCH421k",pair=meta$Replicate)
res = DEA.limma(HTA$expr,paste(meta$Cells,meta$Treatment,sep="."),"NCH601.ctrl","NCH601.TMZ",name="NCH601")
 res = DEA.limma(HTA$expr,paste(meta$Cells,meta$Treatment,sep="."),"NCH601.ctrl","NCH601.TMZ",name="NCH601",pair=meta$Replicate)


	setwd("b:/data/")
	download.file("http://data.sablab.net/astro.RData",mode="wb", destfile="astro.RData")
	load("astro.RData")
	str(X)
	X = X[,grep("hta",colnames(X))]
	group = sub("[.].+","",colnames(X))
	data = X
	key0=NULL;key1=NULL;pair=NULL;name=NULL;folder="";counted=FALSE;norm="none";return.model=FALSE;return.auc=FALSE;adjust="BH";silent = FALSE
	str(DEA.limma(X,group))

} 
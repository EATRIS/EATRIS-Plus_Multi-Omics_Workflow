## A warp-up function to plot PCA 
##     x - data frame or matrix. Objects(samples) in columns, properties(genes) in rows
##     scale - do we scale the data before PCA (makes all features of the same importance)?
##     ipc - vector of 2 elements: which PC should be plotted (x,y)
##     col - dot color (scalar or vector)
##     pch - dot shape (scalar or vector)
##     cex - dot size (scalar or vector)
##     cex.names - if not NULL, - the size of the text to be written on top of dots
##     col.names - color of the names

plotPCA = function(	x=NULL,
					pc=NULL,
					scale=TRUE,
					imputed = NULL,
					ipc=1:2,
					randomize=FALSE,
					col="blue",
					pch=19,
					cex=1,
					main="PCA",
					cex.names=NULL,
					col.names="black",
					cex.main=1.2,
					xlim=NULL,
					ylim=NULL,
					return.object=FALSE){
	if (is.null(x)){
		if (is.null(pc)){
			stop("Either data `x` or prcomp-object `pc` should be defined!")
		}
		x = t(pc$x)
	}
	x = as.matrix(x)
	if (length(col)==1) col=rep(col,ncol(x))
	if (length(cex)==1) cex=rep(cex,ncol(x))
	if (length(pch)==1) pch=rep(pch,ncol(x))
	if (length(cex.names)==1) cex.names=rep(cex.names,ncol(x))
	if (length(col.names)==1) col.names=rep(col.names,ncol(x))
	
	if (is.null(pc)){
		## if there are NA - impute
		if (sum(is.na(x))>0){
			if (is.null(imputed)) imputed = median(x,na.rm=TRUE)
			x[is.na(x)]=imputed
		}
		## if variance is 0 - exclude object
		if (sum(apply(x,1,var)==0)>0) x = x[apply(x,1,var)>0,]
		##calculate PC
		pc=prcomp(t(x),scale = scale)
	}
	## plot
	i1 = ipc[1]
	i2 = ipc[2]
	if (randomize){
		idx = sample(1:ncol(x))
	}else{
		idx = 1:ncol(x)
	}
	plot(pc$x[idx,i1],pc$x[idx,i2],col=col[idx],pch=pch[idx],cex=cex[idx],cex.main=cex.main,
		xlim=xlim,ylim=ylim,
		main=sprintf("%s (%d%% variability)",main,round((pc$sdev[i1]^2+pc$sdev[i2]^2)/sum(pc$sdev^2)*100)),
		xlab=sprintf("PC%d,  %d%% variability",i1,round(pc$sdev[i1]^2 /sum(pc$sdev^2)*100)),
		ylab=sprintf("PC%d,  %d%% variability",i2,round(pc$sdev[i2]^2 /sum(pc$sdev^2)*100)))
	
	if(!is.null(cex.names)){
		text(pc$x[idx,i1],pc$x[idx,i2],colnames(x)[idx],cex=cex.names[idx],col=col.names[idx])
	}
	
	if (return.object) return(pc)
}


plotTSNE = function(x=NULL,
					tsne=NULL,
					perplexity=30,
					niter=1000,
					imputed = NULL,
					col="blue",
					pch=19,
					cex=1,
					main="tSNE",
					cex.names=NULL,
					col.names="black",
					cex.main=1.2,
					xlim=NULL,
					ylim=NULL,
					return.object=FALSE){
	if (is.null(x)){
		if (is.null(tsne)){
			stop("Either data `x` or tsne-coordinates `tsne` should be defined!")
		}
		x = t(tsne)
	}
	x = as.matrix(x)
	if (length(col)==1) col=rep(col,ncol(x))
	if (length(cex)==1) cex=rep(cex,ncol(x))
	if (length(pch)==1) pch=rep(pch,ncol(x))
	if (length(cex.names)==1) cex.names=rep(cex.names,ncol(x))
	if (length(col.names)==1) col.names=rep(col.names,ncol(x))
	
	if (is.null(tsne)){
		## if there are NA - impute
		if (sum(is.na(x))>0){
			if (is.null(imputed)) imputed = median(x,na.rm=TRUE)
			x[is.na(x)]=imputed
		}
		## if variance is 0 - exclude object
		if (sum(apply(x,1,var)==0)>0) x = x[apply(x,1,var)>0,]
		##calculate PC
		require(Rtsne)
		TS = Rtsne(t(x),perplexity=perplexity,max_iter=niter,check_duplicates =FALSE)$Y
	}else{
		TS=tsne
	}
	
	## plot
	if (0){
		idx = sample(1:ncol(x))
	}else{
		idx = 1:ncol(x)
	}
	plot(TS[idx,1],TS[idx,2],col=col[idx],pch=pch[idx],cex=cex[idx],cex.main=cex.main,
		xlim=xlim,ylim=ylim,
		main=sprintf("%s",main),
		xlab=sprintf("ts1"),
		ylab=sprintf("ts2"))
	
	if(!is.null(cex.names)){
		text(TS[idx,1],TS[idx,2],colnames(x)[idx],cex=cex.names[idx],col=col.names[idx])
	}
	
	if (return.object) return(TS)
}


## From Tony Kaoma - knee plot
### pca = prcomp(X)
# eiv = pca$sdev^2
# eiv.per = eiv/sum(eiv)
# eiv.per.cumsum = cumsum(eiv.per)
# dF = cbind(pc=1:length(eiv.per.cumsum),eiv=eiv.per.cumsum)
# plot(dF)

# ### 2 parameters : step & delta
# step = 50 # 2
# delta  = eiv.per.cumsum[-1] - eiv.per.cumsum[-length(eiv.per.cumsum)]
# ave.diff = sapply(
    # 1:(length(eiv.per.cumsum)-1),function(x){
    # enD = min(length(eiv.per.cumsum),x+step)
    # mean(abs(eiv.per.cumsum[(x+1):enD] - eiv.per.cumsum[x]))
    # }
# )
# delta_ = 25 # length(delta)
# cutoff = tail(which(ave.diff > mean(delta[1:delta_])),1)

# plot(
    # dF,
    # xlab="nbr of PC(s)",
    # ylab="% of var",
    # type="l"
# )
# abline(v=cutoff,h=eiv.per.cumsum[cutoff],lty=2,col="grey")

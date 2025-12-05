################################################################################
## Plots P.D.F. of dataframe (each column as line)
##------------------------------------------------------------------------------
## data - data frame
## xlim, ylim - starting coordinates
## main - title
## col - colors
##------------------------------------------------------------------------------
## (c) GNU GPL P.Nazarov, 2011.  petr.nazarov[at]crp-sante.lu
################################################################################
# data = Data
# col=NULL#rainbow(ncol(data))
# ylim=NULL
# xlim=NULL
# main="P.D.F."
 
plotDataPDF = function (data,main="P.D.F.", xlab="Data value", ylim=NULL, xlim=NULL, col=NULL,cex=1,cex.main=1.2, add.legend=F,x.legend = "topright",lwd=1,add.total=FALSE,new=F) {
	## get columns for showing
	if (class(data)[1] == "data.frame") {
		cl = character(ncol(data))
		for (ima in 1:ncol(data))
			cl[ima] = class(data[,ima])
		idx.col=grep("numeric",cl)
		if (length(idx.col)==0) 
			idx.col=grep("integer",cl)
	} else{
		if (class(data)[1] == "matrix") rownames(data) = NULL
		data = data.frame(data)
		idx.col = 1:ncol(data)
	}
	if (length(idx.col)==0) return("Unable to plot P.D.F - no numeric columns")
	dens = density(as.matrix(data[,idx.col]),na.rm=T)
	
	if (length(ylim) == 0) ylim = c(0,max(dens$y)*2)
	if (length(xlim) == 0) xlim = c(min(dens$x),max(dens$x))
	if (length(col) == 0) col = (1:ncol(data[,idx.col]))+1
	if (length(col) < ncol(data[,idx.col])) col = rep(col,1+ncol(data[,idx.col]) / length(col) )[1:ncol(data[,idx.col])]
	
	if(!new) plot(0,0,col=0,ylim=ylim,xlim=xlim,main=main,cex.main=cex.main,xlab=xlab,ylab="Probability density")
	for (ima in 1:ncol(data[,idx.col])){
		lines(density(data[,idx.col][,ima],na.rm=T),col=col[ima],lwd=lwd)
	}
	if (add.total) lines(dens,lwd=lwd+1)
	
	if (add.legend) 
		if (add.total){
			legend(x=x.legend, legend=c(names(data)[idx.col],"Total"),pch=19,col=c(col,1),cex=cex)
		}else{
			legend(x=x.legend, legend=names(data)[idx.col],pch=19,col=c(col,1),cex=cex)
		}
}
################################################################################
## Function for summarization of the data  (by rows)
## -----------------------------------------------------------------------------
## data - matrix of doubles
## anno - vector/matrix of row annotation
## sum.by.col - column number in anno, where IDs are stored
## sum.func - function for summarization. Mean is the fastest
## get.sd - whether stadard deviations for the replicate should be given
## keep.factor.order - summarized data are sorted alphabeticaly (FALSE), or by appearance in anno (TRUE)
## anno.collapse - which character to use collapsing annotations. Set NA to speed-up avoiding anno-summarization
## -----------------------------------------------------------------------------
## (c)GNU GPL P.Nazarov, 2009.  petr.nazarov[at]crp-sante.lu
## (c)GNU GPL P.Nazarov, 2015.  petr.nazarov[at]lih.lu
##  2012-02-16: correction: bug with NA and "mean" function
##  2015-10-30: added summarization of annotation. When sevral variants appear they are colleased to "v1,v2,v3"
################################################################################

### value assignment for testing reason
 # data = Data
 # anno = Anno
 # sum.by.col="Gene.Symbol"
 # sum.func="mean"
 # get.sd=F
 # get.n=F
 # scaled.sd=F
 # anno.collapse=","
 # keep.factor.order=FALSE
###

##============================================================================
## pasteUnique helping function to deal with summarized annotation strings
pasteUnique= function(x,collapse=","){
	ux = unique(x)
	if (length(ux) == 1){return(as.character(ux))}else{return(paste(ux,collapse=collapse))}
}
##============================================================================
summarizeRows=function(data,anno,
                       sum.by.col=1,sum.func="mean",
                       get.sd=FALSE,get.n=FALSE,scaled.sd=FALSE,
					   anno.collapse=",",
                       keep.factor.order=FALSE) 
{
	METH.FLAG = F
	if (sum.func=="mean") METH.FLAG = T
	if (class(data)[1]!="matrix") data=as.matrix(data)
	## if ID is a vector - transform to matrix
	if (sum.by.col<1) sum.by.col=1
	if (is.null(ncol(anno))) {
		anno=matrix(anno,ncol=1,nrow=length(anno))
		sum.by.col=1
	}
	if (class(anno)[1]!="matrix") anno=as.matrix(anno)
	ng = nrow(data)
	
	## make a id-factor vector for summarization
	items = unique(anno[,sum.by.col])
	if (!keep.factor.order) items = sort(items)
	id.idx = as.integer(factor(anno[,sum.by.col],levels = items))
	
	## sort by id
	idx.sort = sort(id.idx,index.return=T)$ix
	anno[] = anno[idx.sort,]
	id.idx = id.idx[idx.sort]
	data[] = data[idx.sort,]

	## create data storage
	DataSum=list()
	DataSum$ng=max(id.idx)
	DataSum$nma=ncol(data)
	DataSum$anno=matrix(nr=DataSum$ng,nc=ncol(anno))
	DataSum$data=matrix(nr=DataSum$ng,nc=ncol(data))
	if (get.sd) DataSum$sd=matrix(nr=DataSum$ng,nc=ncol(data))
	if (get.n) DataSum$nsum=matrix(nr=DataSum$ng,nc=ncol(data))
	rownames(DataSum$anno) = items
	rownames(DataSum$data) = items
	
	print(sprintf("Performing summarization: %d rows -> %d rows by <%s>",ng,DataSum$ng,sum.func))
	flush.console()
	pb = txtProgressBar(min = 0, max = DataSum$ng, style = 3)
	t0=Sys.time()
	j=1 ## initial index in original set is 1
	
	## main loop - select summarized genes
	for (i in 1:DataSum$ng){
		DataSum$anno[i,]=anno[j,]
		k=0;
		if (j+k+1<=ng){
			while (id.idx[j+k+1]==i & j+k+1<=ng) k=k+1
		}
		if (k>0) {
			if (METH.FLAG){DataSum$data[i,]=colMeans(data[j:(j+k),],na.rm = T)
			}else {DataSum$data[i,]=apply(data[j:(j+k),],2,sum.func,na.rm=T)}
			if (get.sd & scaled.sd) {     
				tmp = t(scale(t(data[j:(j+k),])))
				DataSum$sd[i,]=apply(tmp,2,sd,na.rm=T)
			} else {
				if (get.sd) DataSum$sd[i,]=apply(data[j:(j+k),],2,sd,na.rm=T)
			}  
			if (get.n) DataSum$nsum[i,]=apply(!is.na(data[j:(j+k),]),2,sum)
			if (!is.na(anno.collapse)) DataSum$anno[i,] = apply(anno[j:(j+k),],2,pasteUnique,anno.collapse)
		} else {
			DataSum$data[i,]=data[j,]
			if (get.sd) DataSum$sd[i,]=NA
			if (get.n) DataSum$nsum[i,]=1
		}
		j=j+k+1;
		setTxtProgressBar(pb, i)
	}
	t1=Sys.time()
	print("done in")
	print(difftime(t1,t0))
	
	if (!is.null(colnames(data))) colnames(DataSum$data)=colnames(data)
	if (!is.null(colnames(anno))) colnames(DataSum$anno)=colnames(anno)

	return(DataSum)
}

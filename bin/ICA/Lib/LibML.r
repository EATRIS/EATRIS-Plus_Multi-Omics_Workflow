###############################################################################
## Collection of various warp-ups for Machine Learning:
## (c) GNU GPL Petr Nazarov, Luxembourg Institute of Health, petr.nazarov[at]lih.lu
## last modification 2020-09-02
###############################################################################
## ToDo:
## 1) generalize logistic regression (multinomial case)
## 2) allow svm / rf to work on a single feature?

## install packages required
{
print("----------------------------------------------------------------")
print("libML: Machine learning library")
packages = c("randomForest","e1071")
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

## Calculate confusion matrix
getConfusionMatrix = function(gr, gr.pred) {
	if (class(gr) != "factor") {gr = factor(gr); gr.pred = factor(gr.pred)}
	nm = unique(levels(gr),levels(gr.pred))
	Tab = matrix(nc = length(nm), nr = length(nm))
	rownames(Tab) = paste("pred", nm, sep = ".")
	colnames(Tab) = nm
	for (i in 1:length(nm) )
		for (j in 1:length(nm))
			Tab[i,j] = sum((gr.pred == nm[i]) & (gr== nm[j]))
	return(Tab)
}

# Computes the misclassification error from a confusion matrix.
getMCError = function(CM, balanced=FALSE) {
	if (!balanced){
		return(1-sum(diag(CM))/sum(CM))
	}else{
		err = double(ncol(CM))
		for (i in 1:ncol(CM)) err[i] = 1 - CM[i,i] / sum(CM[,i])
		return(mean(err))
	}
}

## Computes accuracy
getAccuracy = function(CM, balanced=FALSE) {
	if (!balanced){
		sum(diag(CM))/sum(CM)
	}else{
		acc = double(ncol(CM))
		for (i in 1:ncol(CM)) acc[i] = CM[i,i] / sum(CM[,i])
		return(mean(acc))
	}
}


## Calculate sensitivity (for each class)
getSensitivity = function(CM, do.mean = TRUE) {
	sens = double(nrow(CM))
	for (ic in 1:nrow(CM)){
		sens[ic] = CM[ic,ic] / sum(CM[,ic])
	}
	if (do.mean) return(mean(sens))
	if (!do.mean) return(sens)
}

## Calculate sensitivity (for each class)
getSpecificity = function(CM, do.mean = TRUE) {
	spec = double(nrow(CM))
	for (ic in 1:nrow(CM)){
		spec[ic] = (sum(CM) - sum(CM[,ic]) - sum(CM[ic,]) + CM[ic,ic]) / ( sum(CM) - sum(CM[,ic]))
	}
	if (do.mean) return(mean(spec))
	if (!do.mean) return(spec)
}

## Computes F1: TP - diagonal. FP+FN - nondiagonal
getF1 = function(CM) {
	2*sum(diag(CM)) / (sum(diag(CM)) + sum(CM) )
}



## Uniform random numer - integer
randInt = function(min,max){floor(runif(1,min=min,max=max+1))}
## Runs n-fold cross-validation and returns misclassification error (1-accuracy)
## X - matrix of features (features in columns, samples in row)
## Y - vector of labels
## use nfold=0 for LOOCV

## method=c("rf","svm","glm")[3];nfold=5;ntry=10;echo=TRUE;do.mean=TRUE;save.pred=TRUE

runCrossValidation=function(X, Y, method=c("rf","svm","glm")[1], nfold=5, ntry=10, echo=FALSE, do.mean=TRUE,save.pred=FALSE){ #, up.sample=FALSE){
	## check class of Y: char/fact -> classification, numeric -> regression
	if(class(Y) == "character") Y = factor(Y)
	## if 1 variable - put vector x to a matrix X=[x,1] and use glm / lm
	if(is.null(ncol(X))){
		if (echo) cat("Single feature - use Logistic regression (method='glm')\n")
		X = cbind(X,1)
		if (class(Y)=="factor") {
			method = "glm"
		}else{
			method = "lm"
		}
	}
	## check correspondance of X and Y
	if(nrow(X)!=length(Y)) {
		print("ERROR: Number of X-rows should be the same as the length of Y! Stopped.")
		return(FALSE)
	}
	## remove NA labels
	ikeep = !is.na(Y)
	X = X[ikeep,]
	Y = Y[ikeep]
	## annotate X if needed
	if (is.null(rownames(X))) rownames(X) = paste0("sample.",1:nrow(X))
	if (is.null(colnames(X))) colnames(X) = paste0("feature.",1:ncol(X))
	## report
	if (echo) {
		cat(sprintf("Classification cross-validation by `%s` (v.2021-02-28):",method),"\n")
		cat(sprintf("%d features in %d samples (%d samples removed because of NA label)",ncol(X),sum(ikeep),sum(!ikeep)),"\n")
	}
	if (method == "svm") require(e1071)
	if (method == "rf") require(randomForest)
	
	## define groups for n-fold CV. if nfold==0 - use leave-one-out-CV
	if (nfold>0){
		folds = cut(seq(1,length(Y)),breaks=nfold,labels=FALSE)
	}else{ ## LOOCV
		folds = cut(seq(1,length(Y)),breaks=length(Y),labels=FALSE)
	}
	
	Res = list()
	Res$method = method
	Res$nfold = nfold
	Res$ntry = ntry
	## in case of classificarion
	if (class(Y)=="factor"){
		Res$error = double(ntry)
		Res$accuracy = double(ntry)
		Res$b_error = double(ntry)
		Res$b_accuracy = double(ntry)
		Res$specificity = double(ntry)
		Res$sensitivity = double(ntry)
		Res$CM = 0
	}else{ ## regression
		Res$mse = double(ntry)
		Res$r2 = double(ntry)
		Res$rp = double(ntry)
		Res$rs = double(ntry)
	}
	if (save.pred) {
		Res$P = matrix(nrow = nrow(X), ncol = ntry)
		rownames(Res$P) = rownames(X)
		colnames(Res$P) = paste0("run.",1:ntry)
	}
	itry=1
	for (itry in 1:ntry){
		# if(up.sample & class(Y) == "factor"){
			# imax = which.max(summary(Y))
			# cl = "1"
			# for (cl in levels(Y)){
				# if (cl == names(imax)) next
				# toadd = max(summary(Y)) - summary(Y)[cl]
			# }
		# }
		
		idx = sample(length(Y))
		X = X[idx,]
		Y = Y[idx]
		P = Y; P[] = NA
		fold = 1
		k=0
		if (echo) {
			if (nfold>0) cat(sprintf("try %d. %d-fold CV:",itry,nfold))
			if (nfold<=0) cat(sprintf("try %d. LOOCV:\n",itry))
		}
		for (fold in unique(folds)){ ## cross-validation loop
			idx.test = which(fold == folds)
			idx.train = (1:nrow(X))[-idx.test]
			## training
			if (method == "rf") model = randomForest(x=X[idx.train,],y=Y[idx.train])
			if (method == "svm") model = svm(x=X[idx.train,],y=Y[idx.train])
			if (method == "glm") {
				if (class(Y)=="factor"){
					if (nlevels(Y)>2) {
						## ToDo!
						stop("Logistic regression for >2 classes is not implemented. ToDo!")
						## one hot encoding
						Yoh=matrix(0,ncol=nlevels(Y),nrow=length(idx.train))
						Yoh[cbind(1:length(idx.train),as.integer(Y[idx.train]))]=1
						colnames(Yoh) = levels(Y)
						train = data.frame(x = X[idx.train,],y = Yoh)
						names(train)[1:ncol(X)] = paste0("x",1:ncol(X))
						model = list()
						for (i in 1:nlevels(Y))
							model[[levels(Y)[i]]] = glm(y ~ .,data=train,family="binomial")
							#.....bla..... problem with y! exclued from "."!
					}
					if (nlevels(Y)==2) {
						train = data.frame(x = X[idx.train,],y = as.integer(Y[idx.train])-1)
						model = glm(y ~ .,data=train,family="binomial")
					}
				}else{
					stop("GLM for regression is not implemented here yet")
				}
			}
			if (method == "lm") {
				if (class(Y)=="factor"){
					stop("LM for classification is not implemented here yet")
				}else{
					train = data.frame(x = X[idx.train,],y = Y[idx.train])
					names(train)[1:ncol(X)] = paste0("x",1:ncol(X))
					model = lm(y ~ .,data=train)
				}
			}
			
			## prediction
			if (method %in% c("svm","rf")) P[idx.test] = predict(model,X[idx.test,])
			if (method == "glm"){
				test = data.frame(x = X[idx.test,],y=NA)
				#names(test)[1:ncol(X)] = paste0("x",1:ncol(X))
				test$y = predict(model,test,type="response")
				P[idx.test] = factor(levels(Y)[1+as.integer(test$y >0.5)],levels = levels(Y))
			}
			if (method == "lm" & class(Y) != "factor"){
				test = data.frame(x = X[idx.test,],y=NA)
				names(test)[1:ncol(X)] = paste0("x",1:ncol(X))
				P[idx.test] = predict(model,test,type="response")
			}
			
			## output
			if (echo) {
				cat(".")
				k=k+1
				if (k%%50==0) {
					cat(sprintf("[%d of %d]\n",fold,length(unique(folds))))
				}
				flush.console()
			}
		}## for fold
		
		## if classification
		if (class(Y) == "factor"){
			CM=getConfusionMatrix(Y,P)
			Res$error[itry] = getMCError(CM)
			Res$accuracy[itry] = getAccuracy(CM)
			Res$b_error[itry] = getMCError(CM,balanced=TRUE)
			Res$b_accuracy[itry] = getAccuracy(CM,balanced=TRUE)
			Res$specificity[itry] = getSpecificity(CM)
			Res$sensitivity[itry] = getSensitivity(CM)
			Res$CM = Res$CM + CM
			if (echo) 
				cat(sprintf(" accuracy = %.3f (balanced = %.3f), error=%.3f (%.3f)\n",Res$accuracy[itry],Res$b_accuracy[itry],Res$error[itry],Res$b_error[itry]))
		} else { ## regression
			Res$mse[itry] = sum((Y-P)^2) / length(Y) ## mse
			Res$rp[itry] = cor(Y,P,use="pairwise.complete.obs")
			Res$r2[itry] = Res$rp[itry]^2
			Res$rs[itry] = cor(Y,P,method="spearman",use="pairwise.complete.obs")
			if (echo) 
				cat(sprintf(" R2 = %.3f, MSE=%.3f\n",Res$r2[itry],Res$mse[itry]))
		}
		if(save.pred) Res$P[rownames(X),itry] = P
	} ## for itry
	if (echo) {
		if (class(Y) =="factor"){
			cat(sprintf("Error =%.4f +/- %.4f\n",mean(Res$error),-qt(0.025,ntry-1)*sd(Res$error)/sqrt(ntry) ))
		} else {
			cat(sprintf("MSE =%.4f +/- %.4f, r2 =%.4f +/- %.4f\n",mean(Res$mse),-qt(0.025,ntry-1)*sd(Res$mse)/sqrt(ntry), mean(Res$r2),-qt(0.025,ntry-1)*sd(Res$r2)/sqrt(ntry) ))
		}
	}
	if (class(Y) =="factor"){ Res$CM = Res$CM / ntry }
	if (do.mean & class(Y) =="factor"){
		Res$sd.error = sd(Res$error)
		Res$error = mean(Res$error)
		Res$b_error = mean(Res$b_error)
		Res$sd.accuracy = sd(Res$accuracy)
		Res$accuracy = mean(Res$accuracy)
		Res$b_accuracy = mean(Res$b_accuracy)
		Res$sd.specificity = sd(Res$specificity)
		Res$specificity = mean(Res$specificity)
		Res$sd.sensitivity = sd(Res$sensitivity)
		Res$sensitivity = mean(Res$sensitivity)
	}
	if (do.mean & class(Y) !="factor"){
		Res$sd.mse = sd(Res$mse)
		Res$mse = mean(Res$mse)
		Res$sd.r2 = sd(Res$r2)
		Res$r2= mean(Res$r2)
		Res$sd.rp = sd(Res$rp)
		Res$rp= mean(Res$rp)
		Res$sd.rs = sd(Res$rs)
		Res$rs= mean(Res$rs)
	}
	return(Res)
}

if (FALSE){
	## test regression
	X = iris[,1:3]
	Y = iris[,4]
	method="lm"
	echo=TRUE
	ntry=3
	nfold=5
	save.pred=TRUE
	res = runCrossValidation(X,Y,method="rf", echo="TRUE",save.pred="TRUE")
}


###############################################
## survival cross-validation

## Features - a vector or matrix of predictive features (features in columns, samples in row)
## SurvData - a data frame with 2 columns: time and event
# ikeep = RVS$Var$dataset == tcga
# Features = t(RIC$M[,ikeep])
# SurvData = RVS$Sur[ikeep,]
# weights = apply(RIC$stab,2,mean)

runSurvivalCrossValidation = function(Features, SurvData, weights=1, nfold = 5, ntry=10, echo=TRUE, thr.fdr=0.05, scaled=TRUE){
	require(survival)
	if (class(Features)[1]!="matrix") Features = as.matrix(Features)
	if (class(SurvData)!="data.frame") SurvData = as.data.frame(SurvData)
	if (!(all(names(SurvData)%in%c("time","event")))) {
		print("ERROR:Survival data should have 2 columns named 'time' and 'event'.")
		return(FALSE)
	}
	if(nrow(Features)!=nrow(SurvData)) {
		print("ERROR: The number of Feature rows should be the same as number of rows of SurvData.")
		return(FALSE)
	}
	if(length(weights)<ncol(Features)) weights = rep(weights,ncol(Features))[1:ncol(Features)]
	
	#ikeep = !is.na(SurvData$time) & SurvData$time > 0 & SurvData$time
	#Features = Features[ikeep,]
	#SurvData = SurvData[ikeep,]
	
	if (is.null(rownames(Features))) rownames(Features) = paste0("sample.",1:nrow(Features))
	if (is.null(colnames(Features))) colnames(Features) = paste0("feature.",1:ncol(Features))
	
	if (echo) {
		cat(sprintf("Survival prediction cross-validation (v.2020-09-14):"),"\n")
		cat(sprintf("%d features in %d samples (%d samples have NA or 0 survival time)",ncol(Features),nrow(Features),sum(is.na(SurvData$time)|SurvData$time<=0)),"\n")
	}
	
	RES = list()
	RES$score = matrix(0,nrow=nrow(Features),ncol=ntry)
	rownames(RES$score) = rownames(Features)
	colnames(RES$score) = paste0("try",1:ntry)
	RES$pv = double(ntry)+NA
	RES$lhr = double(ntry)+NA
	folds = cut(seq(1,nrow(Features)),breaks=nfold,labels=FALSE)
	itry = 1
	for (itry in 1:ntry){
		if (echo) cat(sprintf("Try #%d",itry))
		irand = sample(nrow(Features))  ## random shuffle
		X = Features[irand,]
		Y = SurvData[irand,]
		## scale the features
		if (scaled) X[]=t(scale(t(X)))
		## risk score vector
		
		## ToDo: glitchi in returned score - randomly permutated. Need to change!
		
		fold = 1 
		for (fold in unique(folds)){ ## cross-validation loop
			idx.test = which(fold == folds)
			idx.train = (1:nrow(X))[-idx.test]
			## remove NA in training set
			idx.train = idx.train[!is.na(Y$event[idx.train]*Y$time[idx.train])]
			## do Cox for each features on training set
			SurvTrain = Surv(time = Y$time[idx.train], event =Y$event[idx.train])
			SurvTest = Surv(time = Y$time[idx.test], event =Y$event[idx.test])
			Res = data.frame(matrix(nrow=ncol(X),ncol=4))
			names(Res) = c("W","LHR","PV","FDR")
			rownames(Res) = colnames(X)
			ic=1
			if (sum(c(0,SurvTrain[,2]),na.rm=TRUE) > 5){
				for (ic in 1:ncol(X)){
					cox = coxph(SurvTrain ~ X[idx.train,ic])
					Res$W[ic] = weights[ic]
					Res$PV[ic] = summary(cox)$logtest["pvalue"]
					Res$LHR[ic] = log((summary(cox)$conf.int)[1])
					if (!is.finite(Res$LHR[ic])){
						Res$PV[ic] = 1
						Res$LHR[ic] = 0
					}
				}
				Res$FDR = p.adjust(Res$PV, "fdr")
				## put values into risk score
				for (ic in which(Res$FDR<thr.fdr)) 
					RES$score[idx.test,itry] = RES$score[idx.test,itry] + Res$LHR[ic] * Res$W[ic] * X[idx.test,ic]
				if (length(which(Res$FDR<thr.fdr))==0) {
					if (echo) cat("p")
					ic = which.min(Res$PV)
					RES$score[idx.test,itry] = RES$score[idx.test,itry] + Res$LHR[ic] * Res$W[ic] * X[idx.test,ic]
				}
			}else{
				Res$W = weights
				Res$PV[] = 1
				Res$LHR[] = 0
				Res$FDR[] = 1
				RES$score[idx.test,itry] = NA
				if (echo) cat("na")
			}
			if (echo) cat(".")
			flush.console()
		}
		## folds
		## get p-value of "score" association to survival
		if (sum(c(0,Y$event),na.rm=TRUE) > 5){
			cox = coxph(Surv(time = Y$time, event=Y$event) ~ RES$score[,itry])
			RES$pv[itry] = summary(cox)$logtest["pvalue"]
			RES$lhr[itry] = log((summary(cox)$conf.int)[1])
			cat(sprintf(", LHR = %.3f (pv = %.2e)",RES$lhr[itry],RES$pv[itry]),"\n")
		}else{
			RES$pv[itry] = 1
			RES$lhr[itry] = NA
			cat(sprintf(", no data\n"))
		}
	}
	cat(sprintf("Result for %d permutations: LHR = %.4f +/- %.4f (geometric mean p-value = %.2e)\n",
			ntry,
			mean(RES$lhr,na.rm=TRUE),
			sd(RES$lhr,na.rm=TRUE)/sqrt(sum(!is.na(RES$lhr))),
			exp(mean(log(RES$pv),na.rm=TRUE))
			))
	return(RES)
}
##############################################################################################

## plot univariable predictore
## x - matrix with 3 columns: mean, q0.025, q0.975
plotPredictors = function(x=NULL,is.log=FALSE,xlab="Coefficient", main="Predictive Features"){
	if (is.null(nrow(x))) stop("in plotPredictors(): x should be a matrix ot data.frame")
	if (is.null(rownames(x))) rownames(x)=paste0("feature",1:nrow(x))
	xlim = round(range(x[,2:3]),1)+c(-0.1,0.1)
	xlim[xlim> 20] = 20
	xlim[xlim< -20]= -20
	barplot(x[,1],
			names.arg=rownames(x),
			horiz=TRUE,las=1,font=2,font.lab=2,
			xlab=xlab,
			main=main,
			col=0,border=0,
			xlim=xlim)
	abline(v=ifelse(is.log,0,1),col=2,lty=2)
	for (i in 1:nrow(x)){
		ix = c(x[i,2],x[i,3])
		iy = rep(i*1.2-0.5,2)
		lines(ix,iy,lwd=2)
		lines(ix[c(1,1)],iy + c(0.05,-0.05)*5/nrow(x),lwd=2)
		lines(ix[c(2,2)],iy + c(0.05,-0.05)*5/nrow(x),lwd=2)
	}
	points(x[,1],1:nrow(x)*1.2-0.5,pch=19,col="grey",cex=1.2)
	points(x[,1],1:nrow(x)*1.2-0.5,pch=1,col=1,cex=1.5)
}

##=============================================================================
## num2fact - transforms numeric data to factors
##-----------------------------------------------------------------------------
num2fact = function(x, nlev = 4,digits=1){
	if (!class(x) %in% c("numeric","double","integer")) return(factor(x))
	f = rep(ifelse(digits==1,sprintf("q%d",1),sprintf("q%02d",1)),length(x))
	thr = quantile(x,seq(1/nlev,1,by = 1/nlev),na.rm=TRUE)
	for (ilev in 2:nlev)
		f[x > thr[ilev-1]] = ifelse(digits==1,sprintf("q%d",ilev),sprintf("q%02d",ilev))
	f[is.na(x)]=NA
	f=factor(f)
	return(f)
}
##=============================================================================
## df2fact - transforms data.frame to factors. Excludes factors with too small or large number of treatments
##-----------------------------------------------------------------------------
df2factor = function(Var,skip=TRUE,minlev =2, maxlev=Inf, NAs = c("","NA"), nlev=4, max.na.prop=0.9){
	if (class(Var)=="matrix") Var=as.data.frame(Var,stringsAsFactors=FALSE)
	if (!class(Var)=="data.frame") return(NULL)
	ikeep = logical(ncol(Var))|TRUE
	for (i in 1:ncol(Var)){
		if (sum(is.na(Var[[i]]))/nrow(Var) > max.na.prop) ikeep[i]= FALSE
		if (class(Var[[i]])%in%c("integer","numeric")) {
			Var[[i]] = num2fact(Var[[i]],nlev=nlev)
		}else{
			Var[[i]][Var[[i]] %in% NAs ] =NA
			Var[[i]] = factor(Var[[i]])
		}
		if (nlevels(Var[[i]])<minlev) ikeep[i]= FALSE
		if (nlevels(Var[[i]])>maxlev) ikeep[i]= FALSE
	}
	if (skip) Var=Var[,ikeep] 
	return(Var)
}

## DEBUGING:
if(FALSE){
	CM = cbind(c(90,5),c(200,800))
	getAccuracy(CM)
	getSensitivity(CM,T)
	getSpecificity(CM,T)
	
	library(e1071)
	model = svm(Species~Sepal.Length+Sepal.Width+Petal.Length+Petal.Width,data=iris[sample(1:150,50),])
	CM = getConfusionMatrix(iris[,5], predict(model,iris[,-5]))
	getAccuracy(CM)
	getSensitivity(CM,T)
	getSpecificity(CM,T)
}


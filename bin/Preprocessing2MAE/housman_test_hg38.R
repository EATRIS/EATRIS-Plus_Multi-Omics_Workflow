library(nlme)

### script from the publication Houseman et al 2012 BMC Bioinformatics

#### arguments needed

# reference dataset in methrix format
# included regions
# pivot table for the reference
# dataset of interest
# included cell types -> it has to overlap with the reference dataset and the pivot table

# add a small script that opens the pivot table and gives the user a set of cell types that can be included

list_cell_types <- function(genome_assembly="hg19", pivot=NULL){

  
  if (genome_assembly=="hg19" & is.null(pivot)){
  pivot <- readRDS("N:/data/Reference_methylomes/analysis/pivot_hg19.RDS")
  } else if (genome_assembly=="hg38" & is.null(pivot)) {
    pivot <- readRDS("N:/data/Reference_methylomes/analysis/pivot_hg38.RDS")
  }
  

  return(colnames(pivot))
}





housman_cell_type <- function(m=NULL, reference_meth=NULL, included_regions=NULL, pivot=NULL,
                               included_cell_types=NULL, genome_assembly="hg19"){

  if (!(genome_assembly %in% c("hg19", "hg38"))){
    stop("This genome assembly is not yet supported. ")
  }
  #later modify this to use the data from methrix. Or download it? Depends on the size.
  if (is.null(reference_meth) & genome_assembly=="hg19"){
    reference_meth <- readRDS("N:/data/Reference_methylomes/analysis/processed_methrix_hg19.RDS")
  } else if (is.null(reference_meth) & genome_assembly=="hg38"){
    reference_meth <- readRDS("N:/data/Reference_methylomes/analysis/processed_methrix.RDS")
  }

  if (is.null(included_regions)){
    included_regions <- readRDS( "N:/data/Reference_methylomes/support_data/processed_cell_spec_hg19.RDS")
  } else if (is.null(included_regions) & genome_assembly=="hg38"){
    included_regions <- readRDS("N:/data/Reference_methylomes/support_data/processed_cell_spec.RDS")
  }

  if (is.null(pivot)){
    pivot <- readRDS("N:/data/Reference_methylomes/analysis/pivot_hg19.RDS")
  } else if (is.null(pivot) & genome_assembly=="hg38"){
    pivot <- readRDS("N:/data/Reference_methylomes/analysis/pivot_hg38.RDS")
  }

  ## maybe generalize it for the future, a so a beta table would be enough.
  if (!is(reference_meth, "methrix")){
    stop("A valid methrix object needs to be supplied as a reference.")
  }
  if (is.null(included_cell_types)){
  message("As there are no selected cell types, all of them will be included. \n")
  included_cell_types <- colnames(pivot)
  } else {
    if (!all(included_cell_types %in% colnames(pivot)))
    stop("All included cell types have to appear in the pivot table as separate columns.", call. = FALSE)
  }

  if (!all(rownames(pivot)==colnames(reference_meth))){
    stop("All reference samples have to appear in the pivot table as rows in the same order.", call. = FALSE)
  }
  if (!all(rownames(pivot)==colnames(reference_meth))){
    stop("All reference samples have to appear in the pivot table as rows in the same order.", call. = FALSE)
  }

  if (!is(m, "methrix")){
    stop("A valid methrix object needs to be supplied.")
  }


  theModel = as.formula(paste0("y~", paste(included_cell_types, collapse = "+"))) # Linear model for our samples
  sizeModel = length(included_cell_types)+1 #Number of coefficients in "theModel"

  reference_matrix <- as.data.frame(get_region_summary(reference_meth[,rownames(pivot)[rowSums(pivot[,included_cell_types])>0]], regions=included_regions))
  rownames(reference_matrix) <- paste0(reference_matrix[,1], "_", reference_matrix[,2])
  reference_matrix <- reference_matrix[,-(1:5)]

  pivot <- pivot[rownames(pivot)[rowSums(pivot[,included_cell_types])>0],]
  pivot <- pivot[,included_cell_types]

  test_meth <- as.data.frame(get_region_summary(m, regions=included_regions))
  rownames(test_meth) <- paste0(test_meth[,1], "_", test_meth[,2])

  test_meth <- test_meth[,-(1:5)]
  CpGSelection <- intersect(rownames(test_meth[complete.cases(test_meth),]),
                            rownames(reference_matrix[complete.cases(reference_matrix),]))
  reference_matrix <- reference_matrix[CpGSelection,]
  test_meth <- test_meth[CpGSelection,]
  # Linear transformation of coefficient vector
  #  representing contrast to test F statistic
  L.forFstat = diag(sizeModel)[-1,]  #All non-intercept coefficients
  # Initialize various containers
  sigmaResid = sigmaIcept = nObserved = nClusters = Fstat = rep(NA, nrow(reference_matrix))
  coefEsts = matrix(NA, nrow(reference_matrix), sizeModel)
  coefVcovs =list()

  

  for(j in 1:nrow(reference_matrix)){ # For each CpG

    #Remove missing methylation values
    ii = !is.na(reference_matrix[j,])
    nObserved[j] = sum(ii)
    pivot$y = unlist(reference_matrix[j,])

    try({ # Try to fit a mixed model to adjust for plate

      fit = lm(theModel, data=pivot[ii,])
      fitCoef = fit$coef
      sigmaResid[j] = summary(fit)$sigma
      sigmaIcept[j] = 0
      nClusters[j] = 0
      coefEsts[j,] = fitCoef
      coefVcovs[[j]] = vcov(fit)
      useCoef = L.forFstat %*% fitCoef
      useV = L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
      Fstat[j] = (t(useCoef) %*% solve(useV, useCoef))/sizeModel
    })
  }

  # Name the rows so that they can be easily matched to the target data set
  rownames(coefEsts) = rownames(reference_matrix)
  colnames(coefEsts) = names(fitCoef)

  ####### Projections for mixture experiment data


  Lwbc = diag(sizeModel)[-1,]
  Lwbc[,1] = 1; #Lwbc[1:2,2] = 1 ; # First column is "interception"
  rownames(Lwbc) = colnames(coefEsts)[-1]
  colnames(Lwbc) = colnames(coefEsts)

  ObsMix = projectWBC(
    test_meth,
    coefEsts[ rownames(coefEsts),-ncol(coefEsts)], # due to NA value of linear model
    Lwbc[, -ncol(Lwbc)])

print(round(100*ObsMix,1))


return(ObsMix)



}

############################################
# WBC INFERENCE
############################################

# Permission to use, copy, modify, and distribute this software and its
# documentation for any purpose other than its incorporation into a
# commercial product is hereby granted without fee, provided that the
# above copyright notice appear in all copies and that both that
# copyright notice and this permission notice appear in supporting
# documentation, and that the name of Brown University not be used in
# advertising or publicity pertaining to distribution of the software
# without specific, written prior permission.

# BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
# PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
# ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

projectWBC = function(Y, coefWBC, contrastWBC=NULL, nonnegative=TRUE){

  if(is.null(contrastWBC)) Xmat = coefWBC
  else Xmat = coefWBC %*% t(contrastWBC)

  nCol = dim(Xmat)[2]
  nSubj = dim(Y)[2]

  mixCoef = matrix(0, nSubj, nCol)
  rownames(mixCoef) = colnames(Y)
  colnames(mixCoef) = colnames(Xmat)

  if(nonnegative){
    library(quadprog)

    Amat = diag(nCol)
    b0vec = rep(0,nCol)

    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i]))
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec)$sol
    }
  }
  else{
    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i]))
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
    }
  }

  return(mixCoef)
}

############################################

inferWBCbyLme = function(Y, pheno, modelFix, modelRand, coefWBC,
                         contrastWBC=NULL, detailLevel=1, silentErrors=TRUE){

  M = dim(coefWBC)[1]
  n = dim(Y)[2]
  if(dim(Y)[1] != M) stop("Y does not match coefWBC in dimension\n")

  if(detailLevel == -1){
    lmeVcomp = list(
      mu = matrix(NA,M,n),
      e = matrix(NA,M,n),
      a = list()
    )
  }

  sigmaResid = sigmaIcept = nObserved = nClusters = rep(NA, M)
  coefEsts = list()
  for(j in 1:M){
    ii = !is.na(Y[j,])
    nObserved[j] = sum(ii)
    pheno$y = Y[j,]

    try({
      fit = try(lme(modelFix, random=modelRand, data=pheno[ii,]),
                silent=silentErrors)

      if(inherits(fit,"try-error")){
        fit = lm(modelFix, data=pheno[ii,])
        fitCoef = fit$coef
        sigmaResid[j] = summary(fit)$sigma
        sigmaIcept[j] = 0
        nClusters[j] = 0
        if(detailLevel == -1){
          lmeVcomp$mu[j,ii] = predict(fit)
          lmeVcomp$e[j,ii] = residuals(fit)
          lmeVcomp$a[[j]] = list()
        }
      }
      else{
        fitCoef = fit$coef$fixed
        sigmaResid[j] = fit$sigma
        sigmaIcept[j] = sqrt(getVarCov(fit)[1])
        nClusters[j] = length(fit$coef$random[[1]])
        if(detailLevel == -1){
          lmeVcomp$mu[j,ii] = predict(fit,level=0)
          lmeVcomp$e[j,ii] = residuals(fit,level=1)
          lmeVcomp$a[[j]] = fit$coef$random
        }
      }
      coefEsts[[j]] = fitCoef
    })
  }
  coefEsts = t(matrix(unlist(coefEsts), length(fitCoef), M))
  rownames(coefEsts) = rownames(coefWBC)
  colnames(coefEsts) = names(fitCoef)

  if(is.null(contrastWBC)) Z = cbind(1,coefWBC)
  else Z = cbind(1, coefWBC %*% t(contrastWBC))
  colnames(Z)[1] = "<Intercept>"

  ZtZ = t(Z) %*% Z
  ZtZqr = qr(ZtZ)
  G = solve(ZtZqr, t(Z) %*% coefEsts)

  out = list(method="lme", GammaMatrix=G, Beta1Matrix=coefEsts,
             sigma=cbind(resid=sigmaResid, intercept=sigmaIcept),
             N=cbind(obs=nObserved,clusters=nClusters))

  if(detailLevel>=1){
    out$var.unscaled=solve(ZtZqr)

    d = dim(Z)[2]
    out$predicted = Z %*% G
    residual2 = coefEsts - out$predicted
    out$Sigma = (1/(M-d)) * t(residual2)%*%residual2

    X = model.matrix(modelFix, pheno)
    Xbar = apply(X,2,mean)
    Xctr = (X - matrix(1,n,1) %*% Xbar)
    residual2z = Xctr %*% t(residual2)
    out$ssu = apply(residual2z*residual2z, 2, sum)

    residual1 = Y - coefEsts %*% t(X)
    out$sse = apply(residual1*residual1, 1, sum, na.rm=TRUE)

    residual0 = Y - apply(Y,1,mean,na.rm=TRUE) %*% matrix(1,1,n)
    out$ss0 = apply(residual0*residual0, 1, sum, na.rm=TRUE)

    if(detailLevel>=1){
      out$residuals = list(stage1=residual1, stage2=residual2, inner=t(residual2z))
    }
  }
  if(detailLevel==-1){
    trueLme = which(sapply(lmeVcomp$a, length)>0)
    ch = sort(unique(unlist(lapply(lmeVcomp$a[trueLme], function(u)rownames(u[[1]])))))
    nch = length(ch)
    a = matrix(0, M, nch)
    colnames(a) = ch
    for(j in trueLme){
      aa = lmeVcomp$a[[j]][[1]]
      a[j,rownames(aa)] = aa
    }
    lmeVcomp$ch = names(lmeVcomp$a[[trueLme[1]]])
    lmeVcomp$a = a
    out$lmeVcomp = lmeVcomp
  }

  class(out) = "inferWBC"
  out
}

summary.inferWBC = function(x, xboot=NULL){

  out = list(Coef=x$Gamma)
  class(out) = "inferWBCsummary"

  if(is.null(x$var.unscaled)){
    return(out)
  }

  sUnscaled = sqrt(diag(x$var.unscaled))
  sigma = sqrt(diag(x$Sigma))
  se = outer(sUnscaled,sigma,"*")
  rownames(se) = rownames(x$Gamma)
  out$StdErrNaive=se

  class(out) = "inferWBCsummary"

  if(!is.null(x$ss0)){
    out$RsquareTotal = 1-sum(x$ssu+x$sse)/sum(x$ss0)
    out$RsquareStage1 = 1-sum(x$ssu)/sum(x$ss0-x$sse)
  }

  if(is.null(xboot)) return(out)

  sboot = summary(xboot)

  out$Bias.Single = x$Gamma - sboot$Mean.Single
  out$Bias.Double = x$Gamma - sboot$Mean.Double
  out$SD.Single = sboot$SD.Single
  out$SD.Double = sboot$SD.Double
  out$numBoots = sboot$numBoots

  out
}

print.inferWBC = function(x, digits=4){
  cat("WBC inference object (method = ",x$method,")\n\n",sep="")
  print(summary(x), digits=digits)
}

print.inferWBCsummary = function(x, digits=4, displayFactor=100){
  vars = setdiff(colnames(x$Coef),"(Intercept)")
  if(length(x)==1) {
    cat("(no inference data)\n")
    print(x$Coef[,vars])
    return()
  }

  tag = flag = "?"
  tab = NULL
  for(v in vars){
    cat(v,"\n")
    tab = cbind(Est=x$Coef[,v], StdErr0=x$StdErrNaive[,v])
    if(is.null(x$SD.Single)){
      flag = "StdErr0"
      tag = "naive"
    }
    else{
      tab = cbind(tab, StdErr1=x$SD.Single[,v])
      if(is.null(x$SD.Double)){
        flag = "StdErr1"
        tag = "single bootstrap"
      }
      else{
        tab = cbind(tab, StdErr2=x$SD.Double[,v])
        flag = "StdErr2"
        tag = "double bootstrap"
      }
    }
    tab = tab*displayFactor
    tab = cbind(tab, Zscore = as.vector(x$Coef[,v]/tab[,flag])*displayFactor)
    tab = cbind(tab, Pvalue = 2*pnorm(-abs(tab[,"Zscore"])))
    print(tab, digits=digits)
    cat("\n")
  }

  cat("Inference based on ", tag, " standard errors (", flag, ")\n",sep="")
  if(dim(tab)[2]>4) cat("\tfrom",x$numBoots,"bootstrap iterations\n\n")
  else cat("\n")

  cat("Proportion of total variation explained by WBC:",
      format(x$RsquareTotal,digits=digits),"\n")
  cat("Proportion of stage 1 model explained by WBC:",
      format(x$RsquareStage1,digits=digits),"\n")

}

bootInferWBCbyLme = function(Y, pheno, modelFix, modelRand,
                             coefWBC, contrastWBC=NULL, strata=NULL, R=10,
                             vcovWBC=NULL, # Supply this as a list to compute double-bootstrap
                             degFree=NULL,  # Supply this as a vector to use t-distribution for 2-boot
                             saveBetaCoefficients = FALSE,
                             silentErrors = TRUE
){

  M = dim(coefWBC)[1]
  d = dim(coefWBC)[2]

  if(is.null(contrastWBC)) {
    dBeta0 = d
    contrastWBC = diag(d)
    rownames(contrastWBC) = colnames(coefWBC)
  }
  else {
    dBeta0 = dim(contrastWBC)[1]
  }

  n = dim(pheno)[1]
  if(is.null(strata)) stratList = list(1:n)
  else{
    stratList = split(1:n, strata)
  }
  nStrata = length(stratList)
  stratListN = sapply(stratList, length)

  doDoubleBoot = FALSE
  if(!is.null(vcovWBC)){
    doDoubleBoot = TRUE
    oneR = matrix(1,R,1)
    Beta0 = array(NA, dim=c(M, dBeta0, R))
    for(j in 1:M){
      vchol = chol(vcovWBC[[j]])
      if(is.null(degFree)) noise = rnorm(R*d)
      else noise = rt(R*d, degFree[j])

      bootWBC = t(oneR %*% coefWBC[j,] + matrix(noise, R, d) %*% vchol)
      Beta0[j,,] = contrastWBC %*% bootWBC
    }
  }

  fit0 = inferWBCbyLme(Y, pheno, modelFix, modelRand,
                       coefWBC, contrastWBC, detailLevel=-1, silentErrors=silentErrors)
  nchip = dim(fit0$lmeVcomp$a)[2]
  chips = as.character(pheno[[fit0$lmeVcomp$ch]])

  GammaList = Beta1List = list()
  for(r in 1:R){
    bootsL = list()
    for(s in 1:nStrata){
      bootsL[[s]] = sample(stratList[[s]], stratListN[s], replace=TRUE)
    }
    boots = unlist(bootsL)

    Ystar = fit0$lmeVcomp$mu + fit0$lmeVcomp$e[,boots]
    astar = fit0$lmeVcomp$a[,sample(1:nchip, nchip, replace=TRUE)]

    colnames(astar) = colnames(fit0$lmeVcomp$a)

    Ystar = Ystar + astar[,chips]

    GammaObject = inferWBCbyLme(Ystar, pheno, modelFix, modelRand,
                                coefWBC, contrastWBC, detailLevel=0, silentErrors=silentErrors)

    if(r %% 10==0) cat(r,"\n")
    GammaList[[r]] = GammaObject$GammaMatrix
    Beta1List[[r]] = GammaObject$Beta1Matrix
  }

  Z1 = cbind(1, coefWBC %*% t(contrastWBC))

  out = list()
  #out = list(numBoots=R, tDistribution=!is.null(degFree))

  out$Gamma1 = array(NA, dim=c(dBeta0+1, dim(Beta1List[[1]])[2], R))
  for(r in 1:R) {
    out$Gamma1[,,r] = solve(t(Z1)%*%Z1, t(Z1) %*% Beta1List[[r]])
  }

  dimnames(out$Gamma1) = c(dimnames(GammaList[[1]]), list(1:R))

  if(doDoubleBoot){
    out$Gamma2 = array(NA, dim=c(dBeta0+1, dim(Beta1List[[1]])[2], R))
    for(r in 1:R) {
      Z2 = cbind(1, Beta0[,,r])
      out$Gamma2[,,r] = solve(t(Z2)%*%Z2, t(Z2) %*% Beta1List[[r]])
    }
    dimnames(out$Gamma2) = c(dimnames(GammaList[[1]]), list(1:R))
  }

  if(saveBetaCoefficients) {
    if(doDoubleBoot) out$Beta0 = Beta0
    out$Beta1 = array(NA, dim=c(dim(Beta1List[[1]]),R))
    for(r in 1:R) {
      out$Beta1[,,r] = Beta1List[[r]]
    }
    dimnames(out$Beta1) = c(dimnames(Beta1List[[1]]), list(1:R))
  }

  class(out) = "inferWBCBoot"
  out
}

summary.inferWBCBoot = function(x){
  out = list(numBoots = dim(x$Gamma1)[3])
  out$Mean.Single = apply(x$Gamma1,1:2,mean)
  out$SD.Single = apply(x$Gamma1,1:2,sd)
  if(!is.null(x$Gamma2)){
    out$Mean.Double = apply(x$Gamma2,1:2,mean)
    out$SD.Double = apply(x$Gamma2,1:2,sd)
  }
  out
}

print.inferWBCBoot = function(x, digits=NULL){
  print(summary(x), digits=digits)
}

inferWBCbyLm = function(Y, pheno, modelFix, coefWBC,
                        contrastWBC=NULL, detailLevel=1, silentErrors=TRUE){

  M = dim(coefWBC)[1]
  n = dim(Y)[2]
  if(dim(Y)[1] != M) stop("Y does not match coefWBC in dimension\n")

  if(detailLevel == -1){
    lmVcomp = list(
      mu = matrix(NA,M,n),
      e = matrix(NA,M,n)
    )
  }

  sigmaResid = sigmaIcept = nObserved = nClusters = rep(NA, M)
  coefEsts = list()
  for(j in 1:M){
    ii = !is.na(Y[j,])
    nObserved[j] = sum(ii)
    pheno$y = Y[j,]

    fit = lm(modelFix, data=pheno[ii,])
    fitCoef = fit$coef
    sigmaResid[j] = summary(fit)$sigma
    sigmaIcept[j] = 0
    nClusters[j] = 0
    if(detailLevel == -1){
      lmVcomp$mu[j,ii] = predict(fit)
      lmVcomp$e[j,ii] = residuals(fit)
    }
    coefEsts[[j]] = fitCoef
  }
  coefEsts = t(matrix(unlist(coefEsts), length(fitCoef), M))
  rownames(coefEsts) = rownames(coefWBC)
  colnames(coefEsts) = names(fitCoef)

  if(is.null(contrastWBC)) Z = cbind(1,coefWBC)
  else Z = cbind(1, coefWBC %*% t(contrastWBC))
  colnames(Z)[1] = "<Intercept>"

  ZtZ = t(Z) %*% Z
  ZtZqr = qr(ZtZ)
  G = solve(ZtZqr, t(Z) %*% coefEsts)

  out = list(method="lm", GammaMatrix=G, Beta1Matrix=coefEsts,
             sigma=sigmaResid,N=nObserved)

  if(detailLevel>=1){
    out$var.unscaled=solve(ZtZqr)

    d = dim(Z)[2]
    out$predicted = Z %*% G
    residual2 = coefEsts - out$predicted
    out$Sigma = (1/(M-d)) * t(residual2)%*%residual2

    X = model.matrix(modelFix, pheno)
    Xbar = apply(X,2,mean)
    Xctr = (X - matrix(1,n,1) %*% Xbar)
    residual2z = Xctr %*% t(residual2)
    out$ssu = apply(residual2z*residual2z, 2, sum)

    residual1 = Y - coefEsts %*% t(X)
    out$sse = apply(residual1*residual1, 1, sum, na.rm=TRUE)

    residual0 = Y - apply(Y,1,mean,na.rm=TRUE) %*% matrix(1,1,n)
    out$ss0 = apply(residual0*residual0, 1, sum, na.rm=TRUE)

    if(detailLevel>=1){
      out$residuals = list(stage1=residual1, stage2=residual2, inner=t(residual2z))
    }
  }
  if(detailLevel==-1){
    out$lmVcomp = lmVcomp
  }

  class(out) = "inferWBC"
  out
}

bootInferWBCbyLm = function(Y, pheno, modelFix,
                            coefWBC, contrastWBC=NULL, strata=NULL, R=10,
                            vcovWBC=NULL, # Supply this as a list to compute double-bootstrap
                            degFree=NULL,  # Supply this as a vector to use t-distribution for 2-boot
                            saveBetaCoefficients = FALSE,
                            silentErrors = TRUE
){

  M = dim(coefWBC)[1]
  d = dim(coefWBC)[2]

  if(is.null(contrastWBC)) {
    dBeta0 = d
    contrastWBC = diag(d)
    rownames(contrastWBC) = colnames(coefWBC)
  }
  else {
    dBeta0 = dim(contrastWBC)[1]
  }

  n = dim(pheno)[1]
  if(is.null(strata)) stratList = list(1:n)
  else{
    stratList = split(1:n, strata)
  }
  nStrata = length(stratList)
  stratListN = sapply(stratList, length)

  doDoubleBoot = FALSE
  if(!is.null(vcovWBC)){
    doDoubleBoot = TRUE
    oneR = matrix(1,R,1)
    Beta0 = array(NA, dim=c(M, dBeta0, R))
    for(j in 1:M){
      vchol = chol(vcovWBC[[j]])
      if(is.null(degFree)) noise = rnorm(R*d)
      else noise = rt(R*d, degFree[j])

      bootWBC = t(oneR %*% coefWBC[j,] + matrix(noise, R, d) %*% vchol)
      Beta0[j,,] = contrastWBC %*% bootWBC
    }
  }

  fit0 = inferWBCbyLm(Y, pheno, modelFix,
                      coefWBC, contrastWBC, detailLevel=-1, silentErrors=silentErrors)

  GammaList = Beta1List = list()
  for(r in 1:R){
    bootsL = list()
    for(s in 1:nStrata){
      bootsL[[s]] = sample(stratList[[s]], stratListN[s], replace=TRUE)
    }
    boots = unlist(bootsL)
    Ystar = fit0$lmVcomp$mu + fit0$lmVcomp$e[,boots]

    GammaObject = inferWBCbyLm(Ystar, pheno, modelFix,
                               coefWBC, contrastWBC, detailLevel=0, silentErrors=silentErrors)

    if(r %% 10==0) cat(r,"\n")
    GammaList[[r]] = GammaObject$GammaMatrix
    Beta1List[[r]] = GammaObject$Beta1Matrix
  }

  Z1 = cbind(1, coefWBC %*% t(contrastWBC))

  out = list()

  out$Gamma1 = array(NA, dim=c(dBeta0+1, dim(Beta1List[[1]])[2], R))
  for(r in 1:R) {
    out$Gamma1[,,r] = solve(t(Z1)%*%Z1, t(Z1) %*% Beta1List[[r]])
  }

  dimnames(out$Gamma1) = c(dimnames(GammaList[[1]]), list(1:R))

  if(doDoubleBoot){
    out$Gamma2 = array(NA, dim=c(dBeta0+1, dim(Beta1List[[1]])[2], R))
    for(r in 1:R) {
      Z2 = cbind(1, Beta0[,,r])
      out$Gamma2[,,r] = solve(t(Z2)%*%Z2, t(Z2) %*% Beta1List[[r]])
    }
    dimnames(out$Gamma2) = c(dimnames(GammaList[[1]]), list(1:R))
  }

  if(saveBetaCoefficients) {
    if(doDoubleBoot) out$Beta0 = Beta0
    out$Beta1 = array(NA, dim=c(dim(Beta1List[[1]]),R))
    for(r in 1:R) {
      out$Beta1[,,r] = Beta1List[[r]]
    }
    dimnames(out$Beta1) = c(dimnames(Beta1List[[1]]), list(1:R))
  }

  class(out) = "inferWBCBoot"
  out
}
############################################





residual_beta <- function(meth, pheno, adjustment_columns){
  meth <- as.data.frame(meth)
  residual_values <- apply(meth, 1, function(x) {
  
    pheno$x <- x
    form <- as.formula(paste0("x ~ ", paste(adjustment_columns, collapse = "+")))
    mod <- lm(form, data=pheno)
    pval <- pf(summary(mod)$fstatistic[1], summary(mod)$fstatistic[2], summary(mod)$fstatistic[3], lower.tail = FALSE)
    if (!is.na(pval) && pval < 0.01){
      residuals.lm(mod)
    } else {
      scale(x, scale=F)
    }
    
  })

  return(t(residual_values))
}

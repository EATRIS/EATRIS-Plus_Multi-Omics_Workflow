library(methrix)
source("Z:/omics/EM-seq_UU/scripts/calculate_residuals.R")
meth <- load_HDF5_methrix(dir="Z:/omics/EM-seq_UU/h5_QC_SNP_filtered/")

meth <- coverage_filter(meth, cov_thr = 1, prop_samples = 1)
meth <- methrix::subset_methrix(meth, contigs = 1:22)

beta <- get_matrix(meth, add_loci = F)
rownames(beta) <- paste0(meth@elementMetadata$chr, "_", meth@elementMetadata$start)


#Filtering steps
#1. Random filtering
set.seed(3451)
meth_filt_rand <- meth[sample(1:nrow(meth), 100000, replace = F), ]
meth_filt_rand <-convert_HDF5_methrix(meth_filt_rand)
saveRDS(meth_filt_rand, "Z:/omics/EM-seq_UU/filtered_datasets/meth_filt_rand.RDS")

#2. most variable sites


sd_filter <- order(DelayedMatrixStats::rowSds(beta, na.rm=T), decreasing = T)[1:100000]

meth_filt_sd <- meth[sd_filter, ]
meth_filt_sd <-convert_HDF5_methrix(meth_filt_sd)
saveRDS(meth_filt_sd, "Z:/omics/EM-seq_UU/filtered_datasets/meth_filt_sd.RDS")

pheno <- as.data.frame(readRDS("Z:/omics/EM-seq_UU/scripts/cell_type_composition.RDS"))
all(colnames(meth)==rownames(pheno))

#3. random filtering, cell type correction
set.seed(16633)

beta_filt_rand <- beta[sample(1:nrow(beta), 100000, replace = F), ]
res_beta_filt_rand <- residual_beta(beta_filt_rand, pheno, colnames(pheno))  
rownames(res_beta_filt_rand) <- rownames(beta_filt_rand)
saveRDS(res_beta_filt_rand, "Z:/omics/EM-seq_UU/filtered_datasets/adjusted_filtered_rand.RDS")


#4. SD filtering, after cell type correction

#residuals <- residuals(beta, pheno, colnames(pheno))
library(doParallel)
library(foreach)

doParallel::registerDoParallel(4)

m <- 1000000
res_list <- list()

foreach (i=1:ceiling(nrow(beta)/1000000)) %dopar% {
  
  beta1 <- as.matrix(beta[(((i-1)*m)+1):min(nrow(beta), (((i-1)*m)+m)),])
  res <- residual_beta(beta1, pheno, colnames(pheno))     
  #saveRDS(residuals, paste0("Z:/omics/EM-seq_UU/residuals/res_", i, ".rds"))
  res_list[[i]] <- res[order(rowSds(res, na.rm=T), decreasing = T)[1:100000],] 
  colnames(res_list[[i]]) <- colnames(beta1)
  res <- NULL
  gc()
  
}

saveRDS(res_list, "Z:/omics/EM-seq_UU/residuals/residuals_list_29_06.RDS")

res <- do.call(rbind.data.frame,res_list)
##exclude sex chromosomes, because they were not excluded in this run. 
#The exclusion was done later, see the script above
res <- res[-grep("X|Y_", rownames(res)),]
res <- res[order(rowSds(as.matrix(res), na.rm=T), decreasing = T)[1:100000],]
saveRDS(res, "Z:/omics/EM-seq_UU/filtered_datasets/adjusted_filtered_sd.RDS")



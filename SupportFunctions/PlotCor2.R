## Iterate over data frame with significant linear associations
##
## @param input_df Data frame with significant associations from the RLM models
## @return Correlation ggplot
plot_cor <- function(input_df) {
  
  ## Determine number of rows (number of significant associations)
  rows <- (1 : nrow(input_df))

  ## Iterate over rows and plot correlations
  for (i in rows){
    
    ## Feature and Assay names
    assay <- input_df[i, 'assay']
    assay_name <- strsplit(assay, ' | ')[[1]][[1]]
    feature <- input_df[i, 'feature']
    feature_name <- feature
    

    # Add protein name or gene symbol if applicable
    if (startsWith(assay, 'Proteomics')) {
      feature2 <- protein_mappings[protein_mappings$Feature == feature, 2]
      if (nchar(feature2) > 40) {
        feature2 <- substr(feature2, 1, 30)
      }

      feature_name <- paste0(feature_name, '\n', feature2)
    }
    
    if (startsWith(assay, 'mRNA')) {
      feature2 <- mrna_mappings[mrna_mappings$ensembl_IDs == feature, 2]
      feature_name <- paste0(feature_name, '\n', feature2)
    }
    
    if (startsWith(assay, 'Lipidomics')) {
      feature2 <- lipid_mappings[lipid_mappings$IMTM_Feature_ID == feature, 2]
      feature_name <- paste0(feature_name, '\n', feature2)
    }
    
    
    # Retrieve -omics data
    omics <- assay_list[[assay]]
    
    # Subset pheno df on row names in omics data
    pheno_df <- pheno_df[rownames(omics),]
    
    # Add Age and Outcome columns to omics data frame
    ## Outcome
    omics$Age <- pheno_df[, 'Age']
    omics$BLUP_outcome <- pheno_df[, 'Age_acc_BLUP']
    omics$Levine_outcome <- pheno_df[, 'Age_acc_Levine']
    omics$y_feature <- omics[, feature]
    
    plt1 <- ggplot(data = omics, mapping = aes(y = y_feature, x = BLUP_outcome))  +  
      set_theme + geom_point(shape = 21,  color = color_palette_discrete[[2]], size = 3) + 
      sm_statCorr() + ylab(paste0(assay_name, '\n', feature_name)) + xlab('Age acceleration (BLUP)') +
      scale_fill_discrete(color_palette_discrete) 
    
    plt2 <- ggplot(data = omics, mapping = aes(y = y_feature, x = Levine_outcome))  +  
      set_theme + geom_point(shape = 21, fill = color_palette_discrete[[3]], color = color_palette_discrete[[4]], size = 3) + 
      sm_statCorr() + ylab(paste0(assay_name, '\n', feature_name)) + xlab('Age acceleration (Levine)')
    
    plt3 <- ggplot(data = omics, mapping = aes(y = y_feature, x = Age))  +  
      set_theme + geom_point(shape = 21, fill = color_palette_discrete[[5]], color = color_palette_discrete[[6]], size = 3) + 
      sm_statCorr() + ylab(paste0(assay_name, '\n', feature_name)) 
    
    # Save and print
    print(plt1 + plt2 + plt3)
    ggsave(plt1 + plt2 + plt3, filename = paste0(feature, "_", "_age_acc_cor.png"))
  }
}
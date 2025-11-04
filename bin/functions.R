################################################################################
## Script name:    functions.R
## Script purpose: helper functions to read EATRIS-Plus multi-omics data and 
##                 metadata
## Date created:   Sep 28, 2022
## Author(s):      Anna Niehues
## Email(S):       Anna.Niehues@radboudumc.nl
################################################################################
## Notes:
##
################################################################################


#' Read phenotype and hematology data from phenopackets file
#' 
#' @description
#' `read_phenopackets` returns a data frame with phenotype information
#' 
#' @details
#' The information is read from a Phenopackets json file using term labels from
#' controlled vocabularies
#' 
#' @param phenopackets_file_path Path to Phenopackets json file
#' 
#' @returns Data frames 
read_phenopackets <- function(phenopackets_file_path) {
  library(jsonlite)
  pheno_json <- fromJSON(phenopackets_file, flatten = TRUE)
  pheno_df <- do.call(rbind, apply(pheno_json[["members"]], 1, function(x) {
    df1 <- data.frame(
      `subject.id` = x$subject.id,
      `Sex` = x$subject.sex,
      `Age` = sapply( # convert ISO time stamp to numeric age in years
        x$subject.timeAtLastEncounter.age.iso8601duration, function(iso) {
          as.numeric(gsub("P", "", gsub("Y", "", iso)))}),
      `Smoking Status` = x$measurements[
        x$measurements$description == "Smoking Status", 
        "value.ontologyClass.label"],
      `ABO Blood Group` = x$measurements[
        x$measurements$description == "ABO Blood Group", 
        "value.ontologyClass.label"],
      `Rh Blood Group` = x$measurements[
        x$measurements$description == "Rh Blood Group", 
        "value.ontologyClass.label"],
      `Height` = x$measurements[
        x$measurements$description == "Height", 
        "value.quantity.value"],
      `Weight` = x$measurements[
        x$measurements$description == "Weight", 
        "value.quantity.value"],
      `BMI` = x$measurements[
        x$measurements$description == "BMI", 
        "value.quantity.value"]
    )
    if (is.null(x$phenotypicFeatures)) {
      df1[["PhenotypicFeature"]] = NA
    } else {
      df1[["PhenotypicFeature"]] = x$phenotypicFeatures$description
    }
    # hemtalogy
    df2 <- x[["biosamples"]][[
      "measurements"]][[1]][, c("description", "value.quantity.value")]
    df2 <- setNames(data.frame(t(df2[,-1])), df2[,1])
    cbind(df1, df2)
  }))
  # convert categorical variables to factor
  pheno_df$Sex <- as.factor(pheno_df$Sex)
  pheno_df$Smoking.Status <- as.factor(pheno_df$Smoking.Status)
  pheno_df$ABO.Blood.Group <- as.factor(pheno_df$ABO.Blood.Group)
  pheno_df$Rh.Blood.Group <- as.factor(pheno_df$Rh.Blood.Group)
  return(pheno_df)
}


#' Get targeted metabolomics data from directory
#' 
#' @description
#' `get_metabolomics_data` returns a list of data frames containing measured 
#' metabolite levels.
#' 
#' @details
#' The data is extracted from Metabolite Annotation Files (MAF) using sample 
#' names (columns) retrieved from ISA Assay files
#' 
#' @param data_dir Path to directory containing ISA and MAF files
#' 
#' @returns List of data frames 
get_metabolomics_data <- function(data_dir) {
  # # get file names via investigation files
  # isa_investigation_file <- file.path(data_dir, "i_investigation.txt")
  # isa_investigation <- read.delim(
  #   isa_investigation_file, header = FALSE, 
  #   fill = TRUE, stringsAsFactors = FALSE)
  # isa_assay_file_names <- isa_investigation[
  #   isa_investigation[1] == "Study Assay File Name"]
  # isa_assay_file_names <- isa_assay_file_names[
  #   isa_assay_file_names != "" & isa_assay_file_names != "Study Assay File Name"]
  
  # get file names from directory
  isa_assay_file_names <- list.files(
    path = data_dir, pattern = "a_metabolomics")
  maf_file_names <- list.files(
    path = data_dir, pattern = "m_metabolomics")
  
  # assay names
  assay_names <- sapply(isa_assay_file_names, function(x) {
    gsub(".txt", "", gsub("a_", "", x))})
  stopifnot(assay_names == sapply(maf_file_names, function(x) {
    gsub("_v2_maf.tsv", "", gsub("m_", "", x))}))
  names(assay_names) <- assay_names
  names(isa_assay_file_names) <- assay_names
  names(maf_file_names) <- assay_names
  
  # read ISA assay metadata 
  isa_assay <- lapply(isa_assay_file_names, function(x) {
    read.delim(file.path(data_dir, x))})
  
  # read metabolites annotation files
  maf <- lapply(maf_file_names, function(x) {
    read.delim(file.path(data_dir, x),
               row.names = "metabolite_identification",)
  })
  
  # get sample names from assay files
  sample_names <- lapply(isa_assay, function(a) {
    sapply(a$Sample.Name, function(x) {make.names(x)})
  })
  
  # get measured metabolite levels from MAF files using samples names
  assay_data <- lapply(assay_names, function(x) {
    df <- maf[[x]][, sample_names[[x]]]
    if (substr(names(df)[[1]], 1, 1) == "X") {
      colnames(df) <- gsub("X", "CZC", colnames(df))}
    df
  })
  
  assay_data
}


#' Get targeted metabolomics feature metadata from directory
#' 
#' @description
#' `get_metabolomics_feature_metadata` returns a list of data frames containing 
#' metabolite annotations.
#' 
#' @details
#' The metabolite annotations areextracted from Metabolite Annotation Files 
#' (MAF) 
#' 
#' @param data_dir Path to directory containing ISA and MAF files
#' 
#' @returns List of data frames 
get_metabolomics_feature_metadata <- function(data_dir) {
  # get file names from directory
  isa_assay_file_names <- list.files(
    path = data_dir, pattern = "a_metabolomics")
  maf_file_names <- list.files(
    path = data_dir, pattern = "m_metabolomics")
  
  # assay names
  assay_names <- sapply(isa_assay_file_names, function(x) {
    gsub(".txt", "", gsub("a_", "", x))})
  stopifnot(assay_names == sapply(maf_file_names, function(x) {
    gsub("_v2_maf.tsv", "", gsub("m_", "", x))}))
  names(assay_names) <- assay_names
  names(isa_assay_file_names) <- assay_names
  names(maf_file_names) <- assay_names
  
  # read ISA assay metadata 
  #isa_assay <- lapply(isa_assay_file_names, function(x) {
  #  read.delim(file.path(data_dir, x))})
  
  # read metabolites annotation files
  maf <- lapply(maf_file_names, function(x) {
    read.delim(file.path(data_dir, x),
               row.names = "metabolite_identification",)
  })
  
  # get sample names from assay files
  #sample_names <- lapply(isa_assay, function(a) {
  #  sapply(a$Sample.Name, function(x) {make.names(x)})
  #})
  
  # get measured metabolite levels from MAF files using samples names
  feature_metadata <- lapply(assay_names, function(x) {
    df <- maf[[x]][, c(
      "database_identifier", "chemical_formula", "smiles", "inchi",
      "database", "database_version", "uri")]
    df <- cbind(data.frame(metabolite_identification = rownames(df)), df)
  })
  
  feature_metadata
}


#' Get targeted metabolomics experimental metadata from directory
#' 
#' @description
#' `get_metabolomics_metadata` returns a list of data frames containing ISA 
#' Assay metadata.
#' 
#' @param data_dir Path to directory containing ISA and MAF files
#' 
#' @returns List of data frames 
get_metabolomics_assay_metadata <- function(data_dir) {
  # # get file names via investigation files
  # isa_investigation_file <- file.path(data_dir, "i_investigation.txt")
  # isa_investigation <- read.delim(
  #   isa_investigation_file, header = FALSE, 
  #   fill = TRUE, stringsAsFactors = FALSE)
  # isa_assay_file_names <- isa_investigation[
  #   isa_investigation[1] == "Study Assay File Name"]
  # isa_assay_file_names <- isa_assay_file_names[
  #   isa_assay_file_names != "" & isa_assay_file_names != "Study Assay File Name"]
  
  # get file names from directory
  isa_assay_file_names <- list.files(
    path = data_dir, pattern = "a_metabolomics")
  
  # assay names
  assay_names <- sapply(isa_assay_file_names, function(x) {
    gsub(".txt", "", gsub("a_", "", x))})
  names(assay_names) <- assay_names
  names(isa_assay_file_names) <- assay_names
  
  # read ISA assay metadata 
  isa_assay <- lapply(isa_assay_file_names, function(x) {
    read.delim(file.path(data_dir, x))})
  
  isa_assay
}
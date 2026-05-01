library(methrix)
source("Z:/omics/EM-seq_UU/scripts/housman_test_hg38.R")

pivot <- readRDS("Z:/omics/EM-seq_UU/scripts/pivot_hg38.RDS")
reference_meth <- readRDS("Z:/omics/EM-seq_UU/scripts/processed_methrix.RDS")
included_regions <- readRDS("Z:/omics/EM-seq_UU/scripts/processed_cell_spec.RDS")

cell_types_blood <- list_cell_types(pivot=pivot)[grep("Blood|Bone_marrow", list_cell_types(pivot=pivot))]

meth <- load_HDF5_methrix(dir="Z:/omics/EM-seq_UU/h5_QC_filtered/")


meth@elementMetadata$chr <- paste0("chr", meth@elementMetadata$chr)


res <- housman_cell_type(m=meth, pivot = pivot, reference_meth = reference_meth, included_cell_types = cell_types_blood, included_regions = included_regions, genome_assembly = "hg38")

saveRDS(res, "Z:/omics/EM-seq_UU/scripts/cell_type_composition.RDS")

#pheatmap::pheatmap(get_matrix(reference_meth))



rr <- meth_f@elementMetadata
rr$end <- rr$start+2
rr <- makeGRangesFromDataFrame(rr)
reference_meth <- subset_methrix(reference_meth, regions = rr)

meth_mat <- cbind(get_matrix(meth_f), get_matrix(reference_meth))
meth_mat <- meth_mat[complete.cases(meth_mat),]
pheatmap::pheatmap(meth_mat)

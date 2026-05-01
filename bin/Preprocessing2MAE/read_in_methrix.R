library(methrix)

bdg_files <- list.files("Z:/omics/EM-seq_UU/coverage/", full.names = T)
hg38_cpgs <- extract_CPGs(ref_genome = "BSgenome.Hsapiens.NCBI.GRCh38")
meth <- read_bedgraphs(files=bdg_files, ref_cpgs = hg38_cpgs, pipeline = "Bismark_cov", 
                       zero_based = FALSE, stranded = T, collapse_strands = T, h5 = T, h5temp = "h5_temp/", 
                       h5_dir = "h5_object/")

save_HDF5_methrix(meth, dir="h5_final/", replace = T)



library(methylclock)
library(methrix)
library(tidyverse)

setwd("Z:/omics/EM-seq_UU/analysis/scripts")
#loadhistory("Z:/omics/EM-seq_UU/analysis/scripts/methylClocks.RData")
meth <- load_HDF5_methrix(dir="Z:/omics/EM-seq_UU/h5_QC_filtered")


HM450 <- data.table::fread("HM450.hg38.manifest.tsv.gz")
HM450 <- as.data.frame(HM450)[!is.na(HM450$CpG_chrm),]
HM450$CpG_chrm <- gsub("chr", "", HM450$CpG_chrm)

HM450 <- makeGRangesFromDataFrame(as.data.frame(HM450), seqnames.field="CpG_chrm", start.field="CpG_beg", 
                                 end.field="CpG_end", na.rm=T, keep.extra.columns = T)

#beta_cg <- methrix::subset_methrix(meth, regions = HM450, overlap_type = "any")
beta_cg <- methrix::get_region_summary(meth, regions=HM450, overlap_type = "any")
beta_cg$cgid <- HM450$Probe_ID[beta_cg$rid]

#beta_cg$cgid <- gsub("_.*", "", beta_cg$cgid)
beta_cg <- beta_cg[!duplicated(beta_cg$cgid),]

methylclock_mat <- as.data.frame(beta_cg)[,-(1:5)]
methylclock_mat <- methylclock_mat %>% 
  column_to_rownames("cgid")

load("methylClocks.RData")

missing_cpgs <- checkClocks(as.matrix(methylclock_mat))

methylclock_mat2 <- methylclock_mat[rownames(methylclock_mat) %in% cpgs.all,]



age <- methylclock_mat2 %>% 
  rownames_to_column("cgid") %>%
  DNAmAge_own(age = meth@colData$Age, cell.count = T, cell.count.data = cell_types, fastImp = T)
saveRDS(age, "../Age_prediction/age.RDS")
age_long <- age %>%
  rename_with(., ~ ifelse(grepl("Acc|id", .x), .x, paste0("pred.", .x, recycle0 = TRUE))) %>%
  pivot_longer(cols = !id) %>%
  separate(name, into = c("type", "clock"), sep="\\." ) %>%
  pivot_wider(names_from = "clock", values_from = "value")

age_very_long <-  age %>%
  rename_with(., ~ ifelse(grepl("Acc|id", .x), .x, paste0("pred.", .x, recycle0 = TRUE))) %>%
  pivot_longer(cols = !id) %>%
  separate(name, into = c("type", "clock"), sep="\\." ) %>%
  inner_join(., as.data.frame(meth@colData), by=c("id", "subject.id"))


age_very_long %>%
ggplot(aes(clock, value, fill=id))+
#geom_violin()+
  geom_point(shape=21)+
  facet_wrap(~type)+theme_bw()+theme(legend.position = "none") 

age_very_long %>%
  filter(type=="ageAcc3") %>%
  ggplot(aes(value, Age, fill=clock))+
  #geom_smooth(method = "lm", se = FALSE, 
  #            color = "black") + 
  #ggpmisc::stat_poly_eq( aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), parse = TRUE)+
  #geom_violin()+
  geom_point(shape=21, color="black")+
  facet_wrap(~clock)+theme_bw()+theme(legend.position = "none") 

age_very_long %>%
  filter(type=="ageAcc3") %>%
  ggplot(aes(Sex,value,  fill=clock))+
  #geom_smooth(method = "lm", se = FALSE, 
  #            color = "black") + 
  #ggpmisc::stat_poly_eq( aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), parse = TRUE)+
  geom_violin()+
  geom_point(shape=21, color="black")+
  facet_wrap(~clock)+theme_bw()+theme(legend.position = "none") 


age_very_long %>%
  filter(type=="ageAcc3") %>%
  ggplot(aes(Smoking.Status,value,  fill=clock))+
  #geom_smooth(method = "lm", se = FALSE, 
  #            color = "black") + 
  #ggpmisc::stat_poly_eq( aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), parse = TRUE)+
  geom_violin()+
  geom_jitter(shape=21, color="black", width = 0.1)+
  facet_wrap(~clock)+theme_bw()+theme(legend.position = "none") 

age_very_long %>%
  filter(type=="ageAcc3") %>%
  ggplot(aes(BMI,value,  fill=clock))+
  geom_smooth(method = "lm", se = FALSE, 
              color = "black") + 
  ggpmisc::stat_poly_eq( aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), parse = TRUE)+
  #geom_violin()+
  geom_jitter(shape=21, color="black", width = 0.1)+
  facet_wrap(~clock)+theme_bw()+theme(legend.position = "none") 

age_very_long %>%
  filter(type=="ageAcc3") %>%
  ggplot(aes(Rh.Blood.Group,value,  fill=clock))+
  #geom_smooth(method = "lm", se = FALSE, 
  #            color = "black") + 
  #ggpmisc::stat_poly_eq( aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), parse = TRUE)+
  #geom_violin()+
  geom_jitter(shape=21, color="black", width = 0.1)+
  facet_wrap(~clock)+theme_bw()+theme(legend.position = "none") 

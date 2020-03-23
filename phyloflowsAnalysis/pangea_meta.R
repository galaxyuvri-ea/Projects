rm(list=ls())
load("MRCPopSample_phsc_stage2_output_newali_250_HKC_analysis/MRCUVRI_phscnetworks_190827.rda")

pairs_conf<-dconfpairs[,c("H1","H2", "LINK_12")]
#View(pairs_conf)
df<-read.csv("metadata_age_gender.csv", row.names = 1)
View(df)

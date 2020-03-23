rm(list=ls())
load("MRCPopSample_phsc_stage2_output_newali_250_HKC_analysis/MRCUVRI_phscnetworks_190827.rda")

dnet<-data.frame(dnet)
dnet<-dnet[c("H1","H2","CLU_SIZE")]
dnet<-dnet[!duplicated(dnet),]
dnet$pair<-paste(dnet$H1,dnet$H2, sep="_")

patristic<-dw[dw$ADJACENT=="TRUE",]
patristic<-dw[c("H1","H2","PATRISTIC_DISTANCE")]
pat<-aggregate(PATRISTIC_DISTANCE ~ ., data=patristic, FUN = mean)
pat$pair<-paste(pat$H1,pat$H2, sep="_")

dnet_pat<-merge(dnet, pat, by="pair")
dnet_pat<-dnet_pat[c("CLU_SIZE","PATRISTIC_DISTANCE")]

clu_size=unique(dnet_pat$CLU_SIZE)

metrics<-data.frame()
for (size in clu_size) {
  clustr<-dnet_pat[dnet_pat$CLU_SIZE==size, ]
  pat_mean<-aggregate(PATRISTIC_DISTANCE ~ ., data = clustr, FUN = mean)
  pat_intqt<-aggregate(PATRISTIC_DISTANCE ~ ., data = clustr, FUN = quantile)
  tmp<-data.frame(cluster_size=size ,mean=pat_mean$PATRISTIC_DISTANCE, interQuartile=pat_intqt$PATRISTIC_DISTANCE)
  metrics<-rbind(metrics, tmp)
  }
metrics


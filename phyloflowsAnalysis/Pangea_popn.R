rm(list=ls())
load("MRCPopSample_phsc_stage2_output_newali_250_HKC_analysis/MRCUVRI_phscnetworks_190827.rda")
#Cluster size distribution
clu_df<-dnet[c("IDCLU","CLU_SIZE")]
clu_df<-clu_df[!duplicated(clu_df),]
clu_df<-data.frame(table(clu_df$CLU_SIZE))

require(ggplot2)
p <- ggplot(data=clu_df, aes(x=Var1, y=Freq, group=1)) + geom_line(color="blue")+geom_point() + theme_bw()
p<- p+ylab("Number of clusters with size n, N(n)")+xlab("Cluster size (n)") 
p<-p+theme(axis.text= element_text(size = 16, face="bold"), axis.title = element_text(size = 18, face = "bold"))
p

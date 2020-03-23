rm(list=ls())
load("MRCPopSample_phsc_stage2_output_newali_250_HKC_analysis/MRCUVRI_phscnetworks_190827.rda")
metadata<- read.csv("metadata_age_gender.csv")
metadata$pair<-paste(metadata$H1, metadata$H2, sep = "_")

dnet<-data.frame(dnet)
dnet<-dnet[c("H1","H2","CLU_SIZE")]
dnet<-dnet[!duplicated(dnet),]
dnet$pair<-paste(dnet$H1, dnet$H2, sep = "_")

meta_net<-merge(metadata ,dnet, by="pair")

metadata<-meta_net

metadata$H1_age<-cut(metadata$H1_age, breaks = c(17, 24, 59), labels = c("18-24", "25-59"))
metadata$H2_age<-cut(metadata$H2_age, breaks = c(17, 24, 59), labels = c("18-24", "25-59"))

patristic<-dw[c("H1","H2","PATRISTIC_DISTANCE")]
pat<-aggregate(PATRISTIC_DISTANCE ~ ., data=patristic, FUN = mean)
pat$pair<-paste(pat$H1,pat$H2, sep="_")

pat_net<-merge(metadata, pat, by="pair")

# Here, we start extracting info per cluster
net2nods<-pat_net[pat_net$CLU_SIZE==2,]
net3_5nods<-pat_net[pat_net$CLU_SIZE%in%c(3,4,5),]
net3_5nods$CLU_SIZE<-"3-5"
net11_15nods<-pat_net[pat_net$CLU_SIZE%in%c(11,15),]
net11_15nods$CLU_SIZE<-"11-15"
#for pairs
net<-net11_15nods
sex<-c(net["H1_sex"], net["H2_sex"])
popn<-c(net["H1_popn"], net["H2_popn"])
age<-c(net["H1_age"], net["H2_age"])
nodes<-c(net["H1"], net["H2"])

df<-data.frame(Nodes=unlist(nodes), Popn=unlist(popn), 
               Age=unlist(age), Sex=unlist(sex))
df$Popn<-as.numeric(factor(df$Popn,levels=c("FF","GP","WHR"), labels=c(1,2,3)))
df$Age<-as.numeric(factor(df$Age,levels=c("18-24","25-59"), labels=c(1,2)))

df<-df[!duplicated(df),]
nodes<-df

edgeList<-net[2:3]

require(igraph)
G<-graph_from_data_frame(d=edgeList, vertices = nodes, directed = F)

assortativity_nominal(G, V(G)$Sex)
assortativity_nominal(G, V(G)$Age)
assortativity_nominal(G, V(G)$Popn)

pat_net<-net[c("CLU_SIZE","PATRISTIC_DISTANCE")]

pat_mean<-aggregate(PATRISTIC_DISTANCE ~ ., data = pat_net, FUN = median)
pat_median<-aggregate(PATRISTIC_DISTANCE ~ ., data = pat_net, FUN = mean)

interqt<-function(x){
  quantile(x, c(0.25, 0.75))
}

pat_intqt<-aggregate(PATRISTIC_DISTANCE ~ ., data = pat_net, FUN = interqt)

metrics<-data.frame(cluster_size=pat_mean$CLU_SIZE ,mean=pat_mean$PATRISTIC_DISTANCE, median=pat_median$PATRISTIC_DISTANCE, interQuartile=pat_intqt$PATRISTIC_DISTANCE)
metrics


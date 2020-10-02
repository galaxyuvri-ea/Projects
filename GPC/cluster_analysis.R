rm(list = ls())

#import the data
gpc<-read.csv("All_GPC.csv")

#clean ART
gpc$ART[gpc$ART%in%c("N","No")]<-"ART NAIVE"
gpc$ART[gpc$ART=="Yes"]<-"On ART"
gpc$ART<-as.character(gpc$ART)
table(gpc$ART)

#clean clustered
gpc$clustered[is.na(gpc$clustered)]<-0
table(gpc$clustered)

#clean source/location
gpc$sourceLocation[gpc$sourceLocation=="non-GPC"]<-"Non-GPC"
gpc$sourceLocation[gpc$sourceLocation==" Central"]<-"Central"
gpc$sourceLocation<-as.character(gpc$sourceLocation)
table(gpc$sourceLocation)

#clean subtype
gpc$Subtype <- sub("^([unassigned | 01_]).*", "ISR", gpc$Subtype)
table(gpc$Subtype)

gpc$sex[gpc$sex%in%c("f","Female")]<-"F"
gpc$sex[gpc$sex%in%c("m","Male")]<-"M"
gpc$sex<-as.character(gpc$sex)
table(gpc$sex)

#Define age groups
gpc$AgeGroups<-cut(gpc$Age, breaks = c(14, 24, 34, 44, 54, 84), labels = c("15-24","25-34","35-44", "45-54", ">54"))
table(gpc$AgeGroups)

#number of clusters in this dataset
gpc$Cluster<-paste0("C",gpc$ClusterNumber)
length(unique(gpc$Cluster))
clu_df<-gpc[c("Cluster","ClusterSize")]
clu_df<-clu_df[!duplicated(clu_df),]
cluster_distrn<-data.frame(table(clu_df$ClusterSize))
cluster_distrn$per<-cluster_distrn$Freq/sum(cluster_distrn$Freq)

require(ggplot2)
p<-ggplot(cluster_distrn, aes(x=Var1, y=Freq, group=1)) + geom_line(color="blue") + geom_point() + theme_classic()
p<-p+xlab("Cluster size") +ylab("Number of clusters")
p

#Here we use cluster information to generate and edgelist
edgeList<-data.frame()
#get all clusters
clus<-as.character(na.omit(unique(gpc$ClusterNumber)))
#iterate over different clusters
for (clu in clus) {
  this_clu<-as.character(na.omit(gpc$SequenceID[gpc$ClusterNumber==clu]))
  #iterate over sequences in each cluster and create a pairwise combinations
  for (seq in this_clu) {
    pair<-combn(this_clu,2)
    #print(length(pairs))
    for (j in 1:dim(pair)[2]) {
      this_pair<-pair[,j]
      seq1<-this_pair[1]
      seq2<-this_pair[2]
      tmp<-data.frame(seq1=seq1,seq2=seq2, clu=clu)
      edgeList<-rbind(edgeList, tmp)
      }
  }
}
#remove duplicate edges
edgeList<-edgeList[!duplicated(edgeList),]
#pick on sequences and leave out the cluster number
edgeList<-edgeList[c(1,2)]

require(igraph)
#generate a graph object from the edge list dataframe
g<-graph_from_data_frame(edgeList,vertices = gpc, directed = F)
#as_data_frame(g, what="vertices")

V(g)$sex<-as.numeric(as.factor(V(g)$sex))
V(g)$sourceLocation<-as.numeric(as.factor(V(g)$sourceLocation))
V(g)$Age<-as.numeric(as.factor(V(g)$AgeGroups))

##Assortativity for entire network for location, age and gender
assortativity_nominal(g, V(g)$sourceLocation)
assortativity_nominal(g, V(g)$Age)
assortativity_nominal(g, V(g)$sex)

#get sequence IDs for different cluster sizes
clus_2<-as.character(
  na.omit(gpc[gpc$ClusterSize==2,]$SequenceID))
clus_3<-as.character(
  na.omit(gpc[gpc$ClusterSize==3,]$SequenceID))
clus_4<-as.character(
  na.omit(gpc[gpc$ClusterSize==4,]$SequenceID))
clus_5<-as.character(
  na.omit(gpc[gpc$ClusterSize==5,]$SequenceID))
clus_6<-as.character(
  na.omit(gpc[gpc$ClusterSize==6,]$SequenceID))
clus_7<-as.character(
  na.omit(gpc[gpc$ClusterSize==7,]$SequenceID))

#get sub-graphs for each cluster size
g2<-induced_subgraph(g, v=clus_2)
g3<-induced_subgraph(g, v=clus_3)
g4<-induced_subgraph(g, v=clus_4)
g5<-induced_subgraph(g, v=clus_5)
g6<-induced_subgraph(g, v=clus_6)
g7<-induced_subgraph(g, v=clus_7)

##################
g<-g7 #change for different sizes; g2, g3, g4, g5, g6
assortativity_nominal(g, V(g)$sourceLocation)
assortativity_nominal(g, V(g)$Age)
assortativity_nominal(g, V(g)$sex)


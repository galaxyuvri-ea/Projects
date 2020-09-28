setwd("C:/Users/ProGadgets.co.ug/Desktop/GPC_Analysis")

#clear environment
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

#SECTION 1. Obtain univariate stats btn clustering and other varibles
#here we use Chi-square tests to test association between
#cluster membership and demographics/clinical variables 

#location Vs clustering
m<-table(gpc$sourceLocation, gpc$clustered)
m
chisq.test(m)

#subtype Vs clustering
m<-table(gpc$Subtype, gpc$clustered)
m 
#We realise that our contigency table contains zeros.
#As such we use simulate.p.value in the Chi-square test
#We obtain a result equivalent to that obtained using
#Fisher's exact test.
chisq.test(m,simulate.p.value = T)

#Age
m<-table(gpc$AgeGroups, gpc$clustered)
m
chisq.test(m)

#Sex
m<-table(gpc$sex, gpc$clustered)
m
chisq.test(m)

#ART
m<-table(gpc$ART, gpc$clustered)
m
#
m<-matrix(data = c(178, 791, 2419, 406), nrow = 2, ncol = 2, byrow = T)
chisq.test(m)

#Multivariable analysis of cluster membership Vs 
#demographics
# logitM<-glm(clustered~Subtype, data = gpc, family = binomial(link = "logit"))
# summary(logitM)
# logitM<-glm(clustered~sourceLocation, data = gpc, family = binomial(link = "logit"))
# summary(logitM)
# logitM<-glm(clustered~ART, data = gpc, family = binomial(link = "logit"))
# summary(logitM)
gpc<-gpc[gpc$ART%in%c("ART NAIVE","On ART"), ]
dim(gpc)

#
gpc$Subtype<-as.factor(gpc$Subtype)
gpc$Subtype<-relevel(gpc$Subtype, ref="A1")
gpc$sex<-as.factor(gpc$sex)
gpc$sex<-relevel(gpc$sex, ref="M")
gpc$ART<-as.factor(gpc$ART)
gpc$ART<-relevel(gpc$ART, ref="ART NAIVE")
gpc$AgeGroups<-as.factor(gpc$AgeGroups)
gpc$AgeGroups<-relevel(gpc$AgeGroups, ref="15-24")

logitM<-glm(clustered~sourceLocation+sex+AgeGroups+ART+Subtype, data = gpc, family = binomial(link = "logit"))
a<-summary(logitM)
a
b<-exp(cbind(OR = coef(logitM), confint(logitM)))
ab<-cbind(b, pvalue=coef(a)[,'Pr(>|z|)'])
View(ab)

#Get overall effect of the factors using wald.test
require(aod)
#Location
wald.test(b = coef(logitM), Sigma = vcov(logitM), Terms = 2:4)
#Sex
wald.test(b = coef(logitM), Sigma = vcov(logitM), Terms = 5:5)
#Age
wald.test(b = coef(logitM), Sigma = vcov(logitM), Terms = 6:9)
#ART
wald.test(b = coef(logitM), Sigma = vcov(logitM), Terms = 10:10)
#Subtype
wald.test(b = coef(logitM), Sigma = vcov(logitM), Terms = 11:16)

#Obtain overall contribution of variables to the pred
#location, subtype vs clustering
#m<-table(gpc$sourceLocation, gpc$clustered, gpc$Subtype)
#m

#pick clusters of size 5 and above
big_clusters <- gpc[gpc$ClusterSize%in%c(5:7),]
clusters<-unique(big_clusters$Cluster)
for (cluster in clusters) {
  tmp<-big_clusters[big_clusters$Cluster==cluster,]
  print(unique(tmp$Cluster))
  print(unique(tmp$ClusterSize))
  print(quantile(tmp$Age))
  print(table(tmp$sex))
  print(table(tmp$ART))
  print(table(tmp$sourceLocation))
  print(table(tmp$Subtype))
}

gpc<-gpc[gpc$Subtype%in%c("A1","D","ISR", "C"),]
gpc<-gpc[gpc$clustered=="1",]
gpc$Subtype<-as.character(gpc$Subtype)
table(gpc$Subtype, gpc$sourceLocation)

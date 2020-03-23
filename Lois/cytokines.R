rm(list = ls())
setwd("~/Desktop/Recent/lois/lois_merged")

meta_data<-read.csv("clustered_samples_meta.csv", row.names = 1)

check_it<-function(x){
  if(x=="C1" | x=="C2"){
    "C1-C2"
  }else if(x=="C4" | x=="C3"){
    "C3-C4"
  }
}

meta_data$Cluster<-sapply(meta_data$Cluster, check_it)

cytokines2<-read.csv("cytokines.csv", row.names = 1)
cytokines1<-read.csv("cytokines_plate1.csv", row.names = 1)

quick_fn<-function(x){
  as.numeric(as.character(x))
}
cytokines1[1:6]<-sapply(cytokines1[1:6], quick_fn)

str(cytokines1)

cytokines2<-cytokines2[rownames(cytokines2)%in%rownames(meta_data),]
cytokines2_meta_data<-meta_data[rownames(meta_data)%in%rownames(cytokines2),]
cytokines2_meta_data<-cytokines2_meta_data[c("labhivtst","wksameno2","Cluster")]
cytokine2_new<- cbind(cytokines2, cytokines2_meta_data)

cytokines1<-cytokines1[rownames(cytokines1)%in%rownames(meta_data),]
cytokines1_meta_data<-meta_data[rownames(meta_data)%in%rownames(cytokines1),]
cytokines1_meta_data<-cytokines1_meta_data[c("labhivtst","wksameno2","Cluster")]
cytokine1_new<- cbind(cytokines1, cytokines1_meta_data)


cytokine_new<- rbind(cytokine1_new, cytokine2_new)

#write.csv(cytokine_new, "cytokines_with_cluster_info.csv")

cytokine_hiv_positive <- cytokine_new[cytokine_new$labhivtst=="Positive",]
cytokine_hiv_negative <- cytokine_new[cytokine_new$labhivtst=="Negative",]

p <- c(25, 50, 75)/100
quantile(cytokine_hiv_negative$TNF.ALPHA,  probs = p, na.rm=T)
quantile(cytokine_hiv_negative$IL.6,  probs = p,na.rm=T)
quantile(cytokine_hiv_negative$IL.8,  probs = p,na.rm=T)
quantile(cytokine_hiv_negative$IL.10,  probs = p,na.rm=T)
quantile(cytokine_hiv_negative$IL.1BETA,  probs = p,na.rm=T)
quantile(cytokine_hiv_negative$IL.1RA,  probs = p,na.rm=T)

quantile(cytokine_hiv_positive$TNF.ALPHA,  probs = p,na.rm=T)
quantile(cytokine_hiv_positive$IL.6,  probs = p,na.rm=T)
quantile(cytokine_hiv_positive$IL.8,  probs = p,na.rm=T)
quantile(cytokine_hiv_positive$IL.10,  probs = p,na.rm=T)
quantile(cytokine_hiv_positive$IL.1BETA,  probs = p,na.rm=T)
quantile(cytokine_hiv_positive$IL.1RA,  probs = p,na.rm=T)

################ compare cytokine levels with cluster and trimester #######

require(ggpubr)
ylim1 = boxplot.stats(cytokine_new$IL.10)$stats[c(1, 5)]
ylim1
groups<-list(c("C1","C2"),c("C2","C3"), c("C3","C4"))#,c("C4","C1"), c("C3","C1"), c("C4","C2")
p<-ggboxplot(cytokine_new, x="Cluster2", y="IL.10", color="Cluster", add="jitter", 
             ylim=c(1.26, 3.4),xlab = "")
p + stat_compare_means(comparisons = groups, 
                       label.y = c(2.2, 2.6, 3.1),
                       tip.length = c(0.002, 0.002,0.002),
                       label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 3.4)+ylab("IL-10 (pg/ml)") # Add global p-value

###################
ylim1 = boxplot.stats(cytokine_new$IL.6)$stats[c(1, 5)]
ylim1
groups<-list(c("C1","C3"),c("C1","C4"),c("C2","C3"))#, c("C3","C1"), c("C4","C2")
p<-ggboxplot(cytokine_new, x="Cluster", y="IL.6", color="Cluster", add="jitter", 
             ylim=c(1.5, 125),xlab = "", outlier.size=-1)
p + stat_compare_means(comparisons = groups, 
                       label.y = c(75, 90, 105),
                       tip.length = c(0.001, 0.001, 0.001),
                       label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 120, show.legend = FALSE)+ylab("IL-6 (pg/ml)") # Add global p-value
compare_means(IL.6~Cluster, comparisons = groups, data = cytokine_new)

###################
ylim1 = boxplot.stats(cytokine_new$IL.8)$stats[c(1, 5)]
ylim1
p<-ggboxplot(cytokine_new, x="Cluster", y="IL.8", color="Cluster", add="jitter", 
             ylim=c(8, 8000),xlab = "")
compare_means(IL.8~Cluster, comparisons = groups, data = cytokine_new)
groups<-list(c("C1","C3"),c("C1","C4"),c("C2","C3"), c("C2","C4"))#, c("C4","C2")

p<-p + stat_compare_means(comparisons = groups, 
                       label.y = c(6000, 6500, 7000, 7500),
                       label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 8000, show.legend = FALSE)+ylab("IL-8 (pg/ml)") # Add global p-value
p

###################
ylim1 = boxplot.stats(cytokine_new$IL.1BETA)$stats[c(1, 5)]
ylim1
p<-ggboxplot(cytokine_new, x="Cluster", y="IL.1BETA", color="Cluster", add="jitter", 
             ylim=c(8, 8000),xlab = "", outlier.size=-1)
compare_means(IL.1BETA~Cluster, comparisons = groups, data = cytokine_new)
groups<-list(c("C1","C3"),c("C1","C4"),c("C2","C3"), c("C2","C4"))#, c("C4","C2")

p<-p + stat_compare_means(comparisons = groups, 
                          label.y = c(6000, 6400, 6800, 7200),
                          label = "p.signif"
)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 8000, show.legend = FALSE) +ylab("IL-1BETA (pg/ml)")# Add global p-value
p


###################
ylim1 = boxplot.stats(cytokine_new$IL.1RA)$stats[c(1, 5)]
ylim1
p<-ggboxplot(cytokine_new, x="Cluster", y="IL.1RA", color="Cluster", add="jitter", 
             ylim=c(6300, 50000),xlab = "", outlier.size=-1)
p
compare_means(IL.1RA~Cluster, comparisons = groups, data = cytokine_new)
groups<-list(c("C1","C3"),c("C1","C4"),c("C2","C3"), c("C2","C4"))#, c("C4","C2")

p<-p + stat_compare_means(comparisons = groups, 
                          label.y = c(42000, 44000, 46000, 48000),
                          label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50000, show.legend = FALSE)+ylab("IL-1RA (pg/ml)") # Add global p-value
p

###################
ylim1 = boxplot.stats(cytokine_new$TNF.ALPHA)$stats[c(1, 5)]
ylim1
cytokine_new_alpha<-cytokine_new[cytokine_new$TNF.ALPHA<10,]
compare_means(TNF.ALPHA~Cluster, comparisons = groups, data = cytokine_new_alpha)
groups<-list(c("C1","C4"),c("C2","C4"))#,c("C4","C1"), c("C3","C1"), c("C4","C2")
p<-ggboxplot(na.omit(cytokine_new_alpha), x="Cluster", y="TNF.ALPHA", color="Cluster", add="jitter", 
             ylim=c(3, 11),xlab = "")
p
p + stat_compare_means(comparisons = groups, 
                       label.y = c(9, 10),label = "p.signif")+
  stat_compare_means(label.y = 11, show.legend = FALSE) +ylab("TNF-ALPHA (pg/ml)")# Add global p-value
dev.off()

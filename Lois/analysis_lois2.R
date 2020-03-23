rm(list = ls())
setwd("~/Desktop/Recent/lois/lois_merged")
physeq<-readRDS("merged.physeq2")

require(phyloseq)
require(ggplot2)
#Standardising by sample depth
total<-median(sample_sums(physeq))
standf<-function(x, t=total) round(t * (x / sum(x)))
M.std<-transform_sample_counts(physeq, standf)
ntaxa(M.std)
#M.f <- filter_taxa(M.std,function(x) sum(x > 10) > (0.01*length(x)) | sum(x) > 0.001*total, TRUE)
#ntaxa(M.f)
physeq<-M.std

################# Alpha diversity ##################
p<-plot_richness(physeq,x = "labhivtst",color = "labhivtst",measures=c("Simpson","Shannon"))
p<-p+geom_boxplot()+geom_point()+theme_bw()+theme(axis.text=element_text(size=14),
                                                  axis.title=element_text(size=14),
                                                  legend.text=element_text(size=13),
                                                  legend.title=element_text(size = 13))+geom_point(size=2)
p$labels$colour<-"HIV status"
p<-p+xlab("")+geom_jitter(width = 0.25)
p
p_hiv<-p$data
p_hiv<-p_hiv[c("samples","variable","value", "labhivtst")]
write.csv(p_hiv, "Alpha_div_hiv_status.csv")
#ggplotly(p)

est <- estimate_richness(physeq, split = TRUE, measures = c("Shannon", "Simpson"))
hiv_div <- cbind(est,sample_data(physeq)[,"labhivtst"])
aggregate(Simpson~labhivtst, FUN=median, data=hiv_div)
aggregate(Shannon~labhivtst, FUN=median, data=hiv_div)
t <- wilcox.test(hiv_div$Shannon~hiv_div$labhivtst)
t

t <- wilcox.test(hiv_div$Simpson~hiv_div$labhivtst)
t

######################## Collapse datset according to variable and use this downstream ############################
physeq_wks<-subset_samples(physeq, wksameno2%in%c("Second Trimester (12-27)","Third Trimester (28-40)")) #"First Trimester (<12)", 
physeq_wks
div_measures<-estimate_richness(physeq_wks, measures=c("Observed", "simpson", "shannon"))
div_measures$samples<-rownames(div_measures)

est <- estimate_richness(physeq_wks, split = TRUE, measures = c("Shannon", "Simpson"))
hiv_div <- cbind(est,sample_data(physeq_wks)[,"wksameno2"])
hiv_div$wksameno2<-as.factor(hiv_div$wksameno2)
aggregate(Simpson~wksameno2, FUN=median, data=hiv_div)
aggregate(Shannon~wksameno2, FUN=median, data=hiv_div)
t <- wilcox.test(hiv_div$Shannon~hiv_div$wksameno2)
t

t <- wilcox.test(hiv_div$Simpson~hiv_div$wksameno2)
t

#require(dunn.test)
#dunn.test::dunn.test(hiv_div$Simpson, hiv_div$wksameno)

p<-plot_richness(physeq_wks,x = "wksameno2",color = "wksameno2",measures=c("Simpson","Shannon"))
p<-p+geom_boxplot()+geom_point()+theme_bw()+theme(axis.text=element_text(size=14),
                                                  axis.title=element_text(size=14),
                                                  legend.text=element_text(size=13),
                                                  legend.title=element_text(size = 13))+geom_point(size=2)
p$labels$colour<-"Weeks of Amenorrhea"
p<-p+xlab("")+theme(axis.text.x = element_text(angle = 90))
p<-p+geom_jitter(width=0.25)
p

p_trimester<-p$data
p_trimester<-p_trimester[c("samples","variable","value", "wksameno")]
#write.csv(p_trimester, "Alpha_div_wksameno.csv")

############ last born year ##################
sample_data(physeq)$lstbornyr<-as.character(sample_data(physeq)$lstbornyr)
physeq_lstborn<-subset_samples(physeq, lstbornyr%in%c("1","2","3")) #"First Trimester (<12)", 
physeq_lstborn
div_measures<-estimate_richness(physeq_lstborn, measures=c("Observed", "simpson", "shannon"))
div_measures$samples<-rownames(div_measures)

est <- estimate_richness(physeq_lstborn, split = TRUE, measures = c("Shannon", "Simpson"))
hiv_div <- cbind(est,sample_data(physeq_lstborn)[,"lstbornyr"])
hiv_div$lstbornyr<-as.factor(hiv_div$lstbornyr)
aggregate(Simpson~lstbornyr, FUN=median, data=hiv_div)
aggregate(Shannon~lstbornyr, FUN=median, data=hiv_div)
t <- kruskal.test(hiv_div$Shannon~hiv_div$lstbornyr)
t

t <- kruskal.test(hiv_div$Simpson~hiv_div$lstbornyr)
t

require(dunn.test)
dunn.test::dunn.test(hiv_div$Simpson, hiv_div$lstbornyr)

change_lstborn <- function(x){
  if(x=="1"){
    return("1-2")
  }
  else if(x=="2"){
    return("3-5")
  }
  else{
    return("Above 5")
  }
}
sample_data(physeq_lstborn)$lstbornyr<-sapply(sample_data(physeq_lstborn)$lstbornyr, change_lstborn)
p<-plot_richness(physeq_lstborn,x = "lstbornyr",color = "lstbornyr",measures=c("Simpson","Shannon"))
p<-p+geom_boxplot()+geom_point()+theme_bw()+theme(axis.text=element_text(size=14),
                                                  axis.title=element_text(size=14),
                                                  legend.text=element_text(size=13),
                                                  legend.title=element_text(size = 13))+geom_point(size=2)
p$labels$colour<-"Age of youngest child"
p<-p+xlab("")+geom_jitter(width = 0.25)
p

p_lstborn<-p$data
p_lstborn<-p_lstborn[c("samples","variable","value", "lstbornyr")]
write.csv(p_lstborn, "Alpha_div_lstbornyr.csv")

#ggplotly(p)
#Reduce the dataset to GS
require(microbiomeSeq)
phy_Gs<-taxa_level(physeq, "GS")
df1<-data.frame(otu_table(phy_Gs))
df_samples<-data.frame()
for (i in 1:dim(df1)[1]) {
  cor_col<-colnames(df1)[df1[i,]==max(df1[i,])][1]
  sample_name<-rownames(df1)[i]
  tmp<-data.frame(sample=sample_name, taxa=cor_col)
  #print(tmp)
  df_samples<-rbind(df_samples, tmp)
}

TopNOTUs<-names(sort(taxa_sums(phy_Gs), TRUE)[1:10])
#TopNOTUs<-TopNOTUs[TopNOTUs!="__Unknowns__"]
topTaxa<-prune_taxa(TopNOTUs, phy_Gs)
div_measures<-estimate_richness(topTaxa, measures=c("Observed", "simpson", "shannon"))
topTaxa<-normalise_data(topTaxa, "relative")

df1<-data.frame(otu_table(topTaxa))
#df1<-df1[which(rownames(df1)!=" "),]

C1<-rownames(df1[df1$Lactobacillus.other.than.L.iners>=0.75,])
C2<-rownames(df1[df1$Lactobacillus.iners>=0.75,])
C3<-rownames(df1[df1$Gardnerella>=0.60,])
#C32<-rownames(df1[df1$Atopobium.Atopobium.vaginae>=0.5,])
#C3<-c(C31, C32)
grps<-c(C1, C2, C3)
C4<-rownames(df1)[!rownames(df1)%in%grps]

metadata<-data.frame(sample_data(topTaxa))
metadata$Cluster[rownames(metadata)%in%C1]<-"C1 (11, 6%)"
metadata$Cluster[rownames(metadata)%in%C2]<-"C2 (24, 13%)"
metadata$Cluster[rownames(metadata)%in%C3]<-"C3 (87, 49%)"
metadata$Cluster[rownames(metadata)%in%C4]<-"C4 (57, 32%)"
sample_data(topTaxa)<-metadata
sample_data(M.std)<-metadata

samples_cluster<-metadata[c("Cluster")]
#write.csv(samples_cluster, "sample_clusters.csv")

p<-plot_bar(topTaxa, fill="OTU")+facet_grid(.~Cluster, scales = "free_x", space = "free_x")+theme_bw()
p$labels$fill<-"Taxa"
p<-p+scale_fill_manual(values= c("green", "navyblue", "#E67F80", "orange", "cyan", "#F333FF", "pink","yellow","red", "#1C5846", "blue", "red", "white", "#717D7E","#7FB3D5","#784212","#5B2C6F",""))
p<-p+ylab("Abundance proportion")+xlab("Participants")+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p<-p+theme(axis.text=element_text(size=14),
           axis.title=element_text(size=14),legend.text=element_text(size=13),
           legend.title=element_text(size = 13))
p<-p+theme(panel.spacing.x = unit(0.8, "lines"))
p

################ Beta diversity #############
GP.ord <- ordinate(M.std, "PCoA", "bray")
p1<-plot_ordination(M.std, GP.ord, color="Cluster", title="")+theme_bw()
p1<-p1+theme_classic()
p1

require(vegan)
bray_dist<-phyloseq::distance(M.std, method = "bray")
tmp<-data.frame(sample_data(M.std))
adonis(bray_dist~Cluster, data=tmp)
adonis(bray_dist~labhivtst, data=tmp)
#adonis(bray_dist~wksameno, data=tmp)

###########################
phy_Gs<-taxa_level(M.std, "GS")
TopNOTUs<-names(sort(taxa_sums(phy_Gs), TRUE)[1:51])
TopNOTUs<-TopNOTUs[TopNOTUs!="__Unknowns__"]
topTaxa<-prune_taxa(TopNOTUs, phy_Gs)
#################################################
require(pheatmap)
df1<-data.frame(t(otu_table(topTaxa)))
df<-data.frame(sample_data(topTaxa))[,c("Cluster","wksameno2","labhivtst")]
colnames(df)<-c("Clusters", "Trimesters", "HIV status")
p<-pheatmap(log2(df1+1), cluster_rows=FALSE, 
            show_rownames=TRUE,show_colnames=F, 
            cluster_cols=TRUE, annotation_col=df,
            annotation_row=NA, cluster_row = TRUE)
pdf(width = 15, height = 10)
p
dev.off()

############# Differential Abundance labhivst ################
require(DESeq2)
otu_table(phy_Gs)<-otu_table(phy_Gs)+1
dds<-phyloseq_to_deseq2(phy_Gs,~labhivtst)
dds<-DESeq(dds, fitType="local")
res<-results(dds)

alpha<-0.05
sigtab<-as.data.frame(res[which(res$padj < alpha & abs(res$log2FoldChange)>0), ])
sigtab$Abundant_Group <- levels(df$`HIV status`)[as.numeric(sigtab$log2FoldChange>0) + 1]
#sigtab=cbind(as(sigtab, "data.frame"), as(tax_table(phy_Gs)[rownames(sigtab), ], "matrix"))
dim(sigtab)
table_sig<-sigtab[c("log2FoldChange","padj")]
write.csv(table_sig, "Diff_abund_HIV_status.csv", row.names = F)
#select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:51]

#create a waterfall plot for differentially abundant taxa
df_wf<-sigtab[c("log2FoldChange","padj", "Abundant_Group")]
df_wf$taxa<-rownames(df_wf)
df_wf<-df_wf[base::order(-df_wf$log2FoldChange),]
wf<-ggplot(data = df_wf, aes(x=taxa, y=log2FoldChange, fill=Abundant_Group))
wf<-wf+geom_bar(stat="identity", width=0.7, position = position_dodge(width=0.4))
wf<-wf+theme_classic()+theme(axis.text.x = element_text(angle=90,hjust = 1, vjust=0.5),
                        axis.text=element_text(size=14),
                        axis.title=element_text(size=14),
                        legend.text=element_text(size=13),
                        legend.title=element_text(size = 13))+xlab("")
wf$labels$fill<-"HIV status"
wf$data$taxa<-factor(wf$data$taxa, levels=df_wf$taxa)
wf

############# Differential Abundance wksameno ################
require(DESeq2)
physeq_wks<-subset_samples(M.std, wksameno2%in%c("Second Trimester (12-27)","Third Trimester (28-40)")) #"First Trimester (<12)", 
physeq_wks

phy_Gs<-taxa_level(physeq_wks, "GS")
otu_table(phy_Gs)<-otu_table(phy_Gs)+1
dds<-phyloseq_to_deseq2(phy_Gs,~wksameno)
dds<-DESeq(dds, fitType="local")
res<-results(dds)

alpha<-0.05
sigtab<-as.data.frame(res[which(res$padj < alpha & abs(res$log2FoldChange)>0), ])
df<-data.frame(sample_data(phy_Gs))[,c("Cluster","wksameno","labhivtst")]
sigtab$Abundant_Group <- levels(df$wksameno)[as.numeric(sigtab$log2FoldChange>0) + 1]
dim(sigtab)
table_sig<-sigtab[c("log2FoldChange","padj")]
write.csv(table_sig, "Diff_abund_wksameno.csv", row.names = F)
#select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:51]

#create a waterfall plot for differentially abundant taxa
df_wf<-sigtab[c("log2FoldChange","padj", "Abundant_Group")]
df_wf$taxa<-rownames(df_wf)
df_wf<-df_wf[base::order(-df_wf$log2FoldChange),]
wf<-ggplot(data = df_wf, aes(x=taxa, y=log2FoldChange, fill=Abundant_Group))
wf<-wf+geom_bar(stat="identity", width=0.7, position = position_dodge(width=0.4))
wf<-wf+theme_classic()+theme(axis.text.x = element_text(angle=90,hjust = 1, vjust=0.5),
                             axis.text=element_text(size=14),
                             axis.title=element_text(size=14),
                             legend.text=element_text(size=13),
                             legend.title=element_text(size = 13))+xlab("")
wf$labels$fill<-"Weeks of Amenorrhea"
wf$data$taxa<-factor(wf$data$taxa, levels=df_wf$taxa)
wf

############# Differential Abundance, last born year ################
physeq_lstborn<-subset_samples(physeq, lstbornyr%in%c("1","4")) 
physeq_lstborn
phy_Gs<-taxa_level(physeq_lstborn, "GS")
otu_table(phy_Gs)<-otu_table(phy_Gs)+1
dds<-phyloseq_to_deseq2(phy_Gs,~lstbornyr)
dds<-DESeq(dds, fitType="local")
res<-results(dds)

alpha<-0.05
sigtab<-as.data.frame(res[which(res$padj < alpha & abs(res$log2FoldChange)>0), ])
df<-data.frame(sample_data(phy_Gs))[,c("Cluster","wksameno","lstbornyr")]
sigtab$Abundant_Group <- levels(as.factor(df$lstbornyr))[as.numeric(sigtab$log2FoldChange>0) + 1]
dim(sigtab)
table_sig<-sigtab[c("log2FoldChange","padj")]
write.csv(table_sig, "Diff_abund_lst_born_yr.csv", row.names = F)
#create a waterfall plot for differentially abundant taxa
df_wf<-sigtab[c("log2FoldChange","padj", "Abundant_Group")]
df_wf$taxa<-rownames(df_wf)
df_wf<-df_wf[base::order(-df_wf$log2FoldChange),]
wf<-ggplot(data = df_wf, aes(x=taxa, y=log2FoldChange, fill=Abundant_Group))
wf<-wf+geom_bar(stat="identity", width=0.7, position = position_dodge(width=0.4))
wf<-wf+theme_classic()+theme(axis.text.x = element_text(angle=90,hjust = 1, vjust=0.5),
                             axis.text=element_text(size=14),
                             axis.title=element_text(size=14),
                             legend.text=element_text(size=13),
                             legend.title=element_text(size = 13))+xlab("")
wf$labels$fill<-"Last born year"
wf$data$taxa<-factor(wf$data$taxa, levels=df_wf$taxa)
wf

#################### examine associations btn microbiome clusters and STIs ############
require(rcompanion)

meta_sti<-data.frame(sample_data(M.std))
fisher.test(meta_sti$Cluster, meta_sti$labhivtst)
tmp<-table(meta_sti$Cluster, meta_sti$labhivtst)# meta_sti[c("Cluster","labhivtst")]
PT<-pairwiseNominalIndependence(tmp,fisher = TRUE,gtest = F, chisq=F, digits=3)
PT

fisher.test(meta_sti$Cluster, meta_sti$labhsv2tst)
tmp<-table(meta_sti$Cluster, meta_sti$labhsv2tst)# meta_sti[c("Cluster","labhivtst")]
PT<-pairwiseNominalIndependence(tmp,fisher = TRUE,gtest = F, chisq=F, digits=3)
PT

fisher.test(meta_sti$Cluster, meta_sti$percln)
fisher.test(meta_sti$Cluster, meta_sti$douch)
fisher.test(meta_sti$Cluster, meta_sti$lstsexp)

write.csv(meta_sti, "clustered_samples_meta.csv")
################################################
# obj=phyloseq_to_metagenomeSeq(phy_Gs)
# require(metagenomeSeq)
# trials = pData(obj)$lstsexp
# heatmapColColors = brewer.pal(12, "Set3")[as.integer(factor(trials))] 
# heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
# plotMRheatmap(obj = obj, n = 30, cexRow = 0.4, cexCol = 0.4, trace = "none", col = heatmapCols, ColSideColors = heatmapColColors)
# 
#M.f <- filter_taxa(meta_phy.kv,function(x) sum(x > 10) > (0.01*length(x)) | sum(x) > 0.001*total, TRUE)
#ntaxa(M.f)

#a<-super.fitZig.kv(physeq = p,factor = "labhivtst",outDir = ".",FileName =c("hiv_status"),
#                   heatmap.descriptor=c("tax_annot"), main=c("Dog G vs. B, taxa merged"), subt=c("subt = FDR < 0.05,|coeff| >= 1.25, >20%+ in either group"), 
#                  ordered=TRUE, p=0.05, FC = 1.25, perc=0.1, extra.cols = c("wksameno"))
#rownames(meta_sti)[!rownames(meta_sti)%in%rownames(pcr_res)]

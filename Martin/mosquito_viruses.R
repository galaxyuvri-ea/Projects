rm(list = ls())


require(reshape2)
require(phyloseq)
require(ggplot2)

blast<-read.csv("~/Downloads/mosquito_viruses.csv")

strsplit_extended <- function(x){
  strsplit(x,"_")[[1]][1]
}

blast$Sample2<-sapply(blast$Sample, strsplit_extended)
blast$Virus<-trimws(blast$Virus)
blast$Species<-gsub("_"," ",blast$Species)
blast$Species<-trimws(blast$Species)
blast$Species[blast$Species=="Culex dutton"]<-"Culex duttoni"
blast$Number<-1

#the sample data
metadata<-blast[c("Sample2","Species", "Site")]
metadata$Site2<-sapply(metadata$Site, strsplit_extended)

arua<-c("Oniba","Yedu","Ambala","Barize","Arua","AruaS4","AruaS5")
kasese<-c("Kyondo","Kidodo","Barrier","Kirembe")  

metadata$District<-metadata$Site2
metadata$District[metadata$Site2%in%arua]<-"Arua"
metadata$District[metadata$Site2%in%kasese]<-"Kasese"
metadata<-metadata[metadata$Distric%in%c("Arua","Kasese"),]
metadata$Site<-NULL
metadata$Site2<-NULL
metadata<-metadata[!duplicated(metadata),]
rownames(metadata)<-metadata$Sample2

############ the abundance table #######################
blast<-blast[!blast$Virus%in%c("Rhabdoviridae environmental sample","Bunyaviridae environmental sample"),]
tmp <- blast[c("Sample2","Virus", "Number")]
tmp<-aggregate(Number~Sample2+Virus, data=tmp, FUN=sum)
abund<-dcast(tmp, Sample2 ~ Virus, value.var = "Number")
abund[is.na(abund)] <- 0
rownames(abund)<-abund$Sample2
abund<-abund[rownames(abund)%in%rownames(metadata),]
abund$Sample2<-NULL
rowSums(abund)
colSums(abund)

########### The taxonomy information ###################
taxon<-blast[c("Virus", "Virusfamily", "Segmentation")]
#do some cleaning here, assign correct family for Yonago Culex iflavirus Unclassified Bunyavirales
taxon$Virusfamily[taxon$Virus=="Yonago Culex iflavirus"]<-"Iflaviviridae"
taxon$Virusfamily[taxon$Virus=="Wuhan Mosquito Virus 2"]<-"Orthomyxoviridae"
taxon$Segmentation[taxon$Virusfamily=="Iflaviviridae"]<-"Non-segmented viruses"
taxon<-taxon[!duplicated(taxon),]
rownames(taxon)<-taxon$Virus

#Check that we have common dimensions of the sample data, abundance data and taxonomy 
dim(metadata)
dim(abund)
dim(taxon)

OTU<-otu_table(abund, taxa_are_rows = F)
TAX<-tax_table(as.matrix(taxon))
SAM<-sample_data(metadata)
physeq<-phyloseq(OTU, SAM, TAX)
physeq

#physeq<-taxa_level(physeq, "Virusfamily")
phy<-physeq
total <- median(sample_sums(phy))
standf <- function(x, t = total) round(t * (x/sum(x)))
M.std <- transform_sample_counts(phy, standf)
M.std <- filter_taxa(M.std, function(x) sum(x > 10) > (0.01 * length(x)) | sum(x) > 0.001 * total, TRUE)
ntaxa(M.std)

#Alpha div
p <- plot_richness(M.std, x = "District", color = "District", measures = c("Observed", "Simpson"))
p <- p + geom_boxplot(outlier.shape = NA) + theme_bw() + 
  theme(axis.text = element_text(size = 10,face = "bold"), axis.title = element_text(size = 12, face = "bold")) + geom_point(size = 2)
p<-p+geom_jitter()+xlab("")
p$labels$colour<-"District"
print(p)

est <- estimate_richness(M.std, split = TRUE, measures = c("Simpson","Observed"))
div_df <- cbind(est,sample_data(M.std)[,"District"])
shapiro.test(div_df$Simpson)
shapiro.test(div_df$Observed)

aggregate(.~District, FUN=median, data=div_df)
kruskal.test(Observed~District, data = div_df)
kruskal.test(Simpson~District, data = div_df)

df <- melt(div_df, id=c("District"), variable_name = "variable")
p1<-ggplot(data=df, aes(x=District,y=value, color=District))+
  geom_boxplot(na.rm = T) + geom_point()+  geom_jitter(height = 0.5)+
  scale_color_manual("District",values=c("#728FCE","#EBCC2A")) +
  xlab('')+ ylab('Alpha diversity measure')+#ggtitle("Lineages and mutations")+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill="white",colour = "black"),
        strip.text = element_text(color = "black", size=15))+
  xlab("Sampling locations")+ 
  facet_wrap(~variable, scales = "free")
p1
#ggsave("~/Desktop/MISC/Alpha_diversity_plot.pdf", plot=p, device="pdf",height=6, width=10)

require(pheatmap)

abund1<-abund[rowSums(abund)>=5,]
abund1<-abund1[,colSums(abund1)>5]
#make sure we are picking the same samples after filtering
metadata1<-metadata[rownames(metadata)%in%rownames(abund1),]
#Now order the metadata with respect to District and mosquito species
metadata1<-metadata1[order(metadata1$District, metadata1$Species),]
abund2<- abund1[rownames(metadata1), ]

p2<-pheatmap(as.matrix(t(log(abund2+1))), 
         cluster_cols = F,
         cluster_rows = T,
         annotation_row = taxon[c(2,3)],
         annotation_col = metadata[c(2,3)],
         treeheight_row = 0)
p2

#Beta diversity, PCOA plots and PERMANOVA
PCoA.ord.bray <- ordinate(M.std, "PCoA", "bray")
#title = c("PCoA of 16S microbiome, Bray-Curtis distance")
#cancer status
PCoA.ord <- plot_ordination(M.std, PCoA.ord.bray, color = "District")
p3<-ggplot(data=PCoA.ord$data, aes(x=Axis.1, y=Axis.2,color=District)) +
  geom_point(size=4) +
  xlab(paste("Axis.1:(",round(100*PCoA.ord.bray$values$Relative_eig[1],1),"%)", sep="")) +
  ylab(paste("Axis.2:(", round(100*PCoA.ord.bray$values$Relative_eig[2],1), "%)", sep="")) +
  scale_color_manual("District",values=c("#728FCE","#EBCC2A")) +
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.title = element_text(color="black", size=15),
        legend.text = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA))+
  stat_ellipse()+ggtitle("PCoA based on Bray-Curtis distance")
p3
bottom<-cowplot::plot_grid(p1, p3, labels = c("B","C"), nrow = 1)
top<-cowplot::plot_grid(p2[[4]], labels = "A", nrow = 1)
all<-cowplot::plot_grid(top, bottom,nrow=2, rel_heights = c(4,3))
ggsave("~/Desktop/MISC/Mosquito_viruses_diversity.pdf", plot = all, device = "pdf", height = 20, width = 17, dpi = 300)
dev.off()
#plot_heatmap(M.std,sample.label = "Species", taxa.label = "Virus")
#tmp$Family<-blast$Virusfamily[match(df2$Virus, blast$Virus)]
# df2$Segmentation<-dt$Segmentation[match(df2$Virus, dt$Virus)]
# 
# df3<-data.frame(table(dt$Virus))
# df3<-df3[df3$Freq>5,]
# 
# df2<-df2[df2$Virus%in%df3$Var1,]
# 
# require(ggplot2)
# segment=unique(df2$Segmentation)[1]
# df2<-df2[df2$Segmentation==segment,]
# 
# df2$VirusID<-as.numeric(as.factor(df2$Virus))
# df2<-df2[df2$Family!="Phasmaviridae", ]
# #df2$number<-log(df2$number+1)
# 
# p<-ggplot(data = df2, aes(x=Species, y=Virus, fill=number))
# p<-p+geom_tile(color = "black") +
#    facet_wrap(~Family,scales = 'free')
# p<-p+theme(axis.text.x = element_text(color="black", size=16,angle = 90, vjust = 0., hjust = 1))+ #vertical alignment of dates and changing size
#   theme(axis.title.x = element_text(color="black", size=20, face="bold"))+
#   theme(axis.text.x = element_text(color="black", size=15))+
#   theme(axis.title.y = element_text(color="black", size=20, face="bold"))+
#   theme(axis.text.y = element_text(color="black", size=15))+
#   theme(legend.text = element_text(size=15))+
#   theme(legend.title = element_text(size=20, face="bold"),
#         strip.text = element_text(size = 25, face = "bold", colour = "black"))+
#   #scale_x_date(date_labels = "%d-%b-%y", date_breaks = "5 days")+ #The date labels will have 4 weeks apart, this reduces cluter on the axis
#         # panel.border = element_rect(fill=NA),
#   xlab('Mosquito species')+
#   ylab('Virus species')
# p$labels$fill<-"Abundance"
# p
# 
# #ggsave(filename = "~/Desktop/non-segmented.png", plot=p, width=20, height=12, device='png', dpi=300)
# ggsave(filename = "~/Desktop/segmented.png", plot=p, width=20, height=12, device='png', dpi=300)

# #Beta diversity, PCOA plots and PERMANOVA
# PCoA.ord.bray <- ordinate(M.std, "PCoA", "bray")
# #title = c("PCoA of 16S microbiome, Bray-Curtis distance")
# #cancer status
# PCoA.ord <- plot_ordination(M.std, PCoA.ord.bray, color = "Species") 
# #title = title)
# PCoA.ord <- PCoA.ord + theme(axis.text = element_text(size = 16, face = "bold"), 
#                              axis.title = element_text(size = 18, face = "bold"), legend.title = element_text(size = 14)) + 
#   theme_bw() + labs(color = "District") + geom_point(size = 5)
# PCoA.ord
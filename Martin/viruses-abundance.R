rm(list = ls())

require(reshape2)
require(phyloseq)
require(ggplot2)
require(grid)
require(gridExtra)

setwd("~/Desktop/Martin")
metadata<-read.csv("metadata.csv")
viruses<-read.csv("mosquito_viruses.csv")

viruses$count<-1
viruses$district<-metadata$District[match(viruses$Pool, metadata$Sample2)]
rownames(metadata)<-metadata$Sample2
viruses$site<-metadata$Site2[match(viruses$Pool, metadata$Sample2)]
viruses<-na.omit(viruses)

tmp<-aggregate(count~Pool+genus, data=viruses, FUN=sum)
abund<-dcast(tmp, Pool ~ genus, value.var = "count")
abund[is.na(abund)] <- 0
rownames(abund)<-abund$Pool
abund$Pool<-NULL
rowSums(abund)
colSums(abund)

taxon<-viruses[c("genus")]
taxon<-taxon[!duplicated(taxon),]
taxon<-data.frame(genus=taxon)
rownames(taxon)<-taxon$genus

abund<-abund[rownames(abund)%in%rownames(metadata),]
metadata<-metadata[rownames(metadata)%in%rownames(abund),]

OTU<-otu_table(abund, taxa_are_rows = F)
TAX<-tax_table(as.matrix(taxon))
SAM<-sample_data(metadata)
physeq<-phyloseq(OTU, SAM, TAX)
physeq

M.std<-physeq
est <- estimate_richness(M.std, split = TRUE, measures = c("Simpson","Observed"))
div_df <- cbind(est,sample_data(M.std)[,c("District","Site2")])
shapiro.test(div_df$Simpson)
shapiro.test(div_df$Richness)

aggregate(.~District, FUN=median, data=div_df)
kruskal.test(Richness~District, data = div_df)
kruskal.test(Simpson~District, data = div_df)

p1<-ggplot(data=div_df, aes(x=District,y=Simpson, color=District))+
  geom_boxplot(na.rm = T) + geom_point()+  geom_jitter(height = 0.5)+
  scale_color_manual("District",values=c("brown","navyblue")) +
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill="white",colour = "black"),
        strip.text = element_text(color = "black", size=15))+
  xlab("District")+ ylab("Simpson's diversity index")
p1

p2<-ggplot(data=div_df, aes(x=District,y=Observed, color=District))+
  geom_boxplot(na.rm = T) + geom_point()+  geom_jitter(height = 0.5)+
  scale_color_manual("District",values=c("brown","navyblue")) +
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill="white",colour = "black"),
        strip.text = element_text(color = "black", size=15))+
  xlab("District")+ ylab("Richness")
p2

##################### Mosquito species diversity analysis ###########################
est <- estimate_richness(M.std, split = TRUE, measures = c("Simpson","Observed"))
div_df <- cbind(est,sample_data(M.std)[,"Species"])
shapiro.test(div_df$Simpson)
shapiro.test(div_df$Observed)

aggregate(.~Species, FUN=median, data=div_df)
kruskal.test(Observed~Species, data = div_df)
kruskal.test(Simpson~Species, data = div_df)

p3<-ggplot(data=div_df, aes(x=Species,y=Simpson))+
  geom_boxplot(na.rm = T) + geom_point()+ geom_jitter(height = 0.5)+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill="white",colour = "black"),
        strip.text = element_text(color = "black", size=15))+
  coord_flip()+xlab('Mosquito species')+ ylab("Simpson's diversity index")
p3

p4<-ggplot(data=div_df, aes(x=Species,y=Observed))+
  geom_boxplot(na.rm = T) + geom_point()+ geom_jitter(height = 0.5)+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill="white",colour = "black"),
        strip.text = element_text(color = "black", size=15))+
  coord_flip()+xlab('Mosquito species')+ ylab("Richness")
p4

grid.arrange(p1,p2,nrow=1,top = textGrob("Genus level alpha diversity across sampling sites",gp=gpar(fontsize=15, fontface="bold")))
grid.arrange(p3,p4,nrow=1,top = textGrob("Genus level alpha diversity across mosquito species",gp=gpar(fontsize=15, fontface="bold")))

##################################
require(pheatmap)
abund$unclassified<-NULL
abund$Unassigned<-NULL
abund$viruses<-NULL
p5<-pheatmap(as.matrix(t(log(abund+1))), 
             cluster_cols = F,
             cluster_rows = T,
             annotation_col = metadata[4],
             treeheight_row = 0,
             main="Virus genera detected amongst mosquito pools")
p5

# top<-cowplot::plot_grid(p1, p2, nrow = 1)
# bottom<-cowplot::plot_grid(p3, p4, nrow = 1)
# all<-cowplot::plot_grid(top, bottom,nrow=2, rel_heights = c(4,3))
# all
#ggsave("~/Desktop/MISC/Alpha_diversity_plot.pdf", plot=p, device="pdf",height=6, width=10)

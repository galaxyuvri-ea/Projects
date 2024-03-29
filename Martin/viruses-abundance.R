# R script to generate alpha diversity plots for the mosquito-viruses dataset
# Author: Alfred Ssekagiri
# Link: https://github.com/galaxyuvri-ea/Projects/blob/master/Martin/viruses-abundance.R
# Created: 07th-August-2022
# Modified: 26th-August-2022

# clean environment
rm(list = ls())

# load required packages
require(reshape2)
require(phyloseq)
require(ggplot2)
require(grid)
require(gridExtra)
require(pheatmap)

# set working directory
setwd("~/Desktop/Martin/")
# import data into R environment
metadata<-read.csv("metadata.csv")
metadata<-metadata[metadata$Mosquito.genera!="Eretmopodites",]

viruses<-read.csv("mosquito_viruses.csv")

viruses$count<-1 # add a dummy counter to each row
# add district info to the virus data frame
viruses$district<-metadata$District[match(viruses$Pool, metadata$Sample2)]
rownames(metadata)<-metadata$Sample2
viruses$site<-metadata$Site2[match(viruses$Pool, metadata$Sample2)]
# get rid of incomplete records, records with NA for some columns
viruses<-na.omit(viruses)
# following Marc's comment, we now put all unclassified genera in one group - unclassified
viruses$genus[viruses$genus%in%c("unclassified","Unassigned","viruses",
                                 "Riboviria","Bunyaviridae","Rhabdoviridae",
                                 "Picornavirales","Densovirinae")] <- "Unclassified"

# obtain total number of hits per genus in each sample/pool
tmp<-aggregate(count~Pool+genus, data=viruses, FUN=sum)
# make a wide version of tmp, here we shall have a counts table with pools/samples in rows and genera in columns
# this will be the abundance table for the phyloseq object we shall create later
abund<-dcast(tmp, Pool ~ genus, value.var = "count")
abund[is.na(abund)] <- 0
rownames(abund)<-abund$Pool
abund$Pool<-NULL
rowSums(abund)
colSums(abund)

# create a taxonomy dataframe to be used for the phyloseq object, we are interested in the genus level at the moment
taxon<-viruses[c("genus")]
taxon<-taxon[!duplicated(taxon),]
taxon<-data.frame(genus=taxon)
rownames(taxon)<-taxon$genus


# create a taxonomy dataframe at family to be used for the phyloseq object, we are interested in the genus level at the moment
taxon_fam<-viruses[c("genus","family")]
taxon_fam<-taxon_fam[!duplicated(taxon_fam),]
taxon_fam<-taxon_fam[taxon_fam$genus!="Unclassified",]

# make sure we have matching samples/pools in abundance and metadta dataframes
abund<-abund[rownames(abund)%in%rownames(metadata),]
metadata<-metadata[rownames(metadata)%in%rownames(abund),]

# create the components of the phyloseq object
OTU<-otu_table(abund, taxa_are_rows = F)
TAX<-tax_table(as.matrix(taxon))
SAM<-sample_data(metadata)
physeq<-phyloseq(OTU, SAM, TAX)
physeq

#================== TASK 1: Generate a heatmap showing the counts of hits for different virus genera detected across districts =================#

# first we sort the metadata by district
metadata <- metadata[order(metadata$District, metadata$Mosquito.genera), ]
# then re-order abund table to be in order created above
abund <- abund[rownames(metadata), ]
# similarly, we sort want to sort the genera by family
taxon_fam <- taxon_fam[order(taxon_fam$family), ]
# then re-order abund table to be in order created above
abund <- abund[,taxon_fam$genus]
# generate the heatmap

#library("RColorBrewer")
#brewer.pal(5, "Set2")
#"#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854"

annotation_colors = list(
  District = c(Arua="brown", Kasese="navyblue"),
  `Mosquito genera` = c(Aedes="#66C2A5", Anopheles="#FC8D62", Coquillettidia="#8DA0CB",
            Culex="#E78AC3", Mansonia="#A6D854"))

metadata_pheat<-metadata
colnames(metadata_pheat)<-sub("\\.", " ", colnames(metadata_pheat))
pheatmap(as.matrix(t(log(abund+1))), 
             cluster_cols = F,
             cluster_rows = F,
         annotation_col = metadata_pheat[c("Mosquito genera","District")],
         treeheight_row = 0,
         annotation_colors = annotation_colors,
            # main="Virus genera detected among mosquito pools from Arua and Kasese districts"
         )

#================== TASK 2: Alpha diversity (Richness and Simpson's diversity) across districts =================#

# compute genera richness and Simpson's diversity index across districts
est <- estimate_richness(physeq, split = TRUE, measures = c("Simpson","Observed"))
div_df <- cbind(est,sample_data(physeq)[,c("District")])

# test for normality of the distribution of the alpha diversity measures
shapiro.test(div_df$Simpson)
shapiro.test(div_df$Observed)

# compute median of alpha diversity measures across districts
aggregate(.~District, FUN=median, data=div_df)

# test for significance of differences of alpha diversity measures across districts
kruskal.test(Observed~District, data = div_df)
kruskal.test(Simpson~District, data = div_df)

# generate a boxplot showing the distribution of genera Simpson's diversity across districts
p1<-ggplot(data=div_df, aes(x=District,y=Simpson, color=District))+
  geom_boxplot(na.rm = T) + geom_point()+  geom_jitter(height = 0.25)+
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
  xlab("District")+ ylab("Simpson's diversity index")+
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))

# generate a boxplot showing the distribution of genera richness across districts
p2<-ggplot(data=div_df, aes(x=District,y=Observed, color=District))+
  geom_boxplot(na.rm = T) + geom_point()+ geom_jitter(height = 0.5)+
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
  xlab("District")+ ylab("Richness")+
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10))

# create a side by side plot of the richness and Simpson's diversity plots across districts
#grid.arrange(p1,p2,nrow=1,top = textGrob("Genus level alpha diversity across districts",gp=gpar(fontsize=15, fontface="bold")))
grid.arrange(p1,p2,nrow=1)
#================== TASK 3: Alpha diversity (Richness and Simpson's diversity) across sampling sites =================#

# compute genera richness and Simpson's diversity index across districts
est <- estimate_richness(physeq, split = TRUE, measures = c("Simpson","Observed"))
div_df <- cbind(est,sample_data(physeq)[,c("Site2","District")])
div_df<-div_df[!div_df$Site2%in%c("Arua-mixed","Kasese-mixed"),]
# test for normality of the distribution of the alpha diversity measures
shapiro.test(div_df$Simpson)
shapiro.test(div_df$Observed)

# compute median of alpha diversity measures across districts
aggregate(Observed~Site2, FUN=median, data=div_df)
aggregate(Simpson~Site2, FUN=median, data=div_df)

# test for significance of differences of alpha diversity measures across districts
kruskal.test(Observed~Site2, data = div_df)
kruskal.test(Simpson~Site2, data = div_df)

# generate a boxplot showing the distribution of genera Simpson's diversity across districts
p3<-ggplot(data=div_df, aes(x=Site2,y=Simpson, color=District))+
  geom_boxplot(na.rm = T) + geom_point()+ geom_jitter(height = 0.5)+
  scale_color_manual("Site2",values=c("brown","navyblue")) +
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+ facet_wrap(~District, scales = "free")+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill="white",colour = "black"),
        strip.text = element_text(color = "black", size=15))+
  xlab("Sampling sites")+ ylab("Simpson's diversity index")+
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))

# generate a boxplot showing the distribution of genera richness across districts
p4<-ggplot(data=div_df, aes(x=Site2,y=Observed, color=District))+
  geom_boxplot(na.rm = T) + geom_point()+ geom_jitter(height = 0.5)+
  scale_color_manual("Site2",values=c("brown","navyblue")) +
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+facet_wrap(~District, scales = "free")+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill="white",colour = "black"),
        strip.text = element_text(color = "black", size=15))+
  xlab("Sampling sites")+ ylab("Richness")+
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10))

# create a side by side plot of the richness and Simpson's diversity plots across districts
# grid.arrange(p3,p4,nrow=2,top = textGrob("Genus level alpha diversity across sampling sites",gp=gpar(fontsize=15, fontface="bold")))
grid.arrange(p3,p4,nrow=2)
#================== TASK 3: Genus level alpha diversity (Richness and Simpson's diversity) across mosquito genera =================#

# compute genera richness and Simpson's diversity index
est <- estimate_richness(physeq, split = TRUE, measures = c("Simpson","Observed"))
div_df <- cbind(est,sample_data(physeq)[,"Mosquito.genera"])
div_df<-div_df[div_df$Mosquito.genera!="Eretmopodites",]
# test for normality of the distribution of the alpha diversity measures
shapiro.test(div_df$Simpson)
shapiro.test(div_df$Observed)

# compute median of alpha diversity measures among mosquito species
aggregate(.~Genus, FUN=median, data=div_df)

# test for significance of differences of alpha diversity measures among mosquito species
kruskal.test(Observed~Genus, data = div_df)
kruskal.test(Simpson~Genus, data = div_df)
dunn.test::dunn.test(div_df$Simpson, div_df$Mosquito.genera)
dunn.test::dunn.test(div_df$Observed, div_df$Mosquito.genera)

# generate a boxplot showing the distribution of genera Simpson's diversity among mosquito species
p5<-ggplot(data=div_df, aes(x=Mosquito.genera,y=Simpson,color=Mosquito.genera))+
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
        xlab('Mosquito species')+ ylab("Simpson's diversity index")+
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))+
  scale_color_brewer(palette = "Set2")

# generate a boxplot showing the distribution of genera richness among mosquito species
p6<-ggplot(data=div_df, aes(x=Mosquito.genera,y=Observed, color=Mosquito.genera))+
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
  xlab('Mosquito species')+ ylab("Richness")+
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10))+
  scale_color_brewer(palette = "Set2")

# create a side by side plot of the richness and Simpson's diversity plots among mosquito species
#grid.arrange(p5,p6,nrow=1,top = textGrob("Genus level alpha diversity among mosquito genera",gp=gpar(fontsize=15, fontface="bold")))
grid.arrange(p5,p6,nrow=1)

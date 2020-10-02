setwd('~/Desktop/marin-prossy')
martin<-read.csv('blastn_summary.csv')
martin$abundance<-1
head(martin)

#Heat map 
df_heatmap<-aggregate(abundance~virus+mosquito+district, 
                      data=martin,sum)

require(ggplot2)
df_heatmap$Abundance<-df_heatmap$abundance
p<-ggplot(data = df_heatmap,
          aes(x=mosquito,y=virus, fill=Abundance))
p<-p+facet_grid(~district,scales = "free_x", space = "free_x")
p<-p+geom_tile()
p<-p+ylab("")+xlab("")+scale_fill_gradient('Abundance', limits=c(0, 4), breaks = c(0, 1, 2, 3, 4),  low = "lightblue", high = "darkblue",na.value = 'black')
p<-p+theme_classic()
p<-p+ theme(text = element_text(size=25),axis.text.x = element_text(angle = 90, hjust = 1))
p
# require(reshape2)
# m<-martin[c('pool', 'virus', 'abundance')]
# m_wide <- dcast(m, pool ~ virus, sum, value.var="abundance")
# rownames(m_wide)<-m_wide$pool
# m_wide$pool<-NULL
# 
# ###metadata
# metadata<-martin[c("pool","mosquito","district")]
# metadata<-metadata[!duplicated(metadata),]
# rownames(metadata)<-metadata$pool
# metadata$pool<-NULL
# 
# require(phyloseq)
# OTU<-otu_table(m_wide, taxa_are_rows = F)
# SAM<-sample_data(metadata)
# 
# phy<-merge_phyloseq(OTU, SAM)
# 
# total = median(sample_sums(phy))
# standf = function(x, t=total) round(t * (x / sum(x)))
# phy = transform_sample_counts(phy, standf)
# 
# div_measures<-estimate_richness(phy)
# df<-div_measures[c("Shannon","Simpson", "Observed")]
# df
# 
# p<-plot_richness(phy, x = "district", measures="Shannon"
#                  ,color = "district")
# p<-p+geom_boxplot()+theme_classic()
# p<-p + theme(legend.position = "none") + xlab('District')
# p
# 
# ###############
# 
# 
# pd<-data.frame(p$data)
# #write.csv(pd, "diversity.csv")
# p<-ggplot(data=pd, aes(x='district', y='value'))
# p<-p+geom_boxplot()
# p

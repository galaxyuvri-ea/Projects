#pcr results
#ignore results with -ve internal control
#summarise the metadata
#age, marital status, income, childno, famplan
setwd("~/Desktop/Recent/lois/lois_merged")

rm(list = ls())
checkt<-function(x){
  return(strtrim(x,6))
}

feature.table <- read.delim("feature-table.tsv", row.names=1)
#View(feature.table)
#colnames(feature.table)<-unlist(lapply(colnames(feature.table), checkt))
taxonomy <- read.delim("taxonomy.tsv", row.names=1)
#View(taxonomy)
metadata <- read.csv("~/Desktop/Recent/lois/metadata.csv", row.names = 1)
metadata$lstsexp<-as.character(metadata$lstsexp)
metadata$percln<-as.character(metadata$percln)
metadata$percln2<-as.character(metadata$percln2)
metadata$douch<-as.character(metadata$douch)
metadata$labhsv2tst<-as.character(metadata$labhsv2tst)
metadata$hivkno<-as.character(metadata$hivkno)
metadata$labhivtst<-as.character(metadata$labhivtst)
metadata$lastmensp<-as.Date(metadata$lastmensp, "%d/%m/%Y")
metadata$coldate<-as.Date(metadata$coldate, "%d/%m/%Y")
metadata[is.na(metadata$wksameno),]$wksameno<-round((metadata[is.na(metadata$wksameno),]$coldate-metadata[is.na(metadata$wksameno),]$lastmensp)/30, 0)
#rownames(metadata)<-metadata$stid_mod

label_status<-function(x){
  if(x=="1"){
    return("Negative")
  }
  else if(x=="2"){
    return("Positive")
  }
}
metadata$labhivtst<-as.factor(sapply(metadata$labhivtst, label_status))

levels_wksameno<-function(x){
  if(is.na(x)){return(NA)}
  else if(x<12){return("First Trimester (<12)")}
  else if(x>11 & x<28){return("Second Trimester (12-27)")}
  else if(x>27 & x<41){return("Third Trimester (28-40)")}
}
metadata$wksameno2<-as.factor(as.character(sapply(metadata$wksameno, levels_wksameno)))
sum(table(metadata$wksameno))
taxonomy[]<-lapply(taxonomy, function(x) gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__|D_7__|D_8__|D_9__|D_10__|D_11__|D_12__|D_13__|D_14__", "", x))
#subset abundance table
#abund_table<-feature.table[colnames(feature.table)%in%rownames(split_log) ,]
#subset taxonomy table
#taxonomy<-taxonomy[rownames(taxonomy)%in%rownames(abund_table), ]
require(splitstackshape)
tmp<-as.data.frame(cSplit(taxonomy, "Taxon", ";"))
rownames(tmp)<-rownames(taxonomy)
tmp<-subset(tmp, Confidence>=0.90)
tmp<-as.data.frame(tmp[, 1:8])
names(tmp)<-c("Confidence","Kingdom","Phylum","Class","Order","Family","Genus","Species")

tmp$Species[is.na(tmp$Species)]<-NA
tmp$Genus[is.na(tmp$Genus)]<-NA
tmp$GS<-""
t<-as.matrix(tmp)
t[,9]<-ifelse(is.na(t[,8]), t[,7], t[,8])
t[,9]<-ifelse((t[,9]=="Lactobacillus"), "Lactobacillus other than L.iners", t[,9])
t[,9]<-ifelse((t[,9]=="Lactobacillus iners AB-1"), "Lactobacillus iners", t[,9])
t[,9]<-ifelse((t[,9]=="Atopobium"), "Atopobium other than A.vaginae", t[,9])
t[,9]<-ifelse((t[,9]=="Aerococcus"), "Aerococcus other than A.christensenii", t[,9])
t[is.na(t)]<-""
tax<-data.frame(t, row.names = rownames(tmp))

feature.table<-feature.table[rownames(feature.table)%in%rownames(tax),]

require(vegan)
require(phyloseq)
require(ggplot2)

OTU<-otu_table(feature.table, taxa_are_rows = TRUE)
TAX<-tax_table(as.matrix(tax))
physeq<-merge_phyloseq(OTU,TAX)
physeq<-prune_samples(sample_sums(physeq)>5000, physeq)
reads <- sample_sums(physeq)

extend_split<-function(x){
  return(strsplit(x, split=".rcb")[[1]][1])
}
merged_samples<-unname(sapply(names(reads), extend_split))
colnames(otu_table(physeq))<-merged_samples

metadata<-metadata[which(rownames(metadata)%in%colnames(otu_table(physeq))),]
sample_data(physeq)<-metadata

physeq<-prune_taxa(taxa_sums(physeq)>0, physeq)
saveRDS(physeq, "merged.physeq2")

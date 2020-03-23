rm(list = ls())
meta_data<-read.csv("clustered_samples_meta.csv", row.names = 1)
#Get number of participants in certain 
#groups stratified by HIV status

#first obtain HIV +ve individuals
hiv_positive <- meta_data[meta_data$labhivtst=="Positive", ]
hiv_negative <- meta_data[meta_data$labhivtst=="Negative", ]

#get median and range of maternal age for HIV + and HIV -ve individuals
range(hiv_positive$age)
median(hiv_positive$age)

range(hiv_negative$age)
median(hiv_negative$age)

#get median and range of gestational age for HIV + and HIV -ve individuals
range(na.omit(hiv_positive$wksameno))
median(na.omit(hiv_positive$wksameno))

range(na.omit(hiv_negative$wksameno))
median(na.omit(hiv_negative$wksameno))

#Obtain parity (number of children) for HIV + and HIV -ve individuals
table(hiv_negative$childno) #HIV -ve
table(hiv_positive$childno) # HIV +ve

table(is.na(hiv_positive$childno)) # HIV positive individuals with no child no record
table(is.na(hiv_negative$childno)) #HIV negative individuals with no child no record

#Obtain marital status (number of children) for HIV + and HIV -ve individuals
table(hiv_negative$marstat) #HIV -ve
table(hiv_positive$marstat) # HIV +ve

#Obtain marital status (number of children) for HIV + and HIV -ve individuals
table(hiv_negative$incosat) #HIV -ve
table(hiv_positive$incosat) # HIV +ve

#Lab diagnoses with respect to HIV status
#HSV2
table(hiv_positive$labhsv2tst)
table(hiv_positive$labhsv2tst)/179

table(hiv_negative$labhsv2tst)
table(hiv_negative$labhsv2tst)/179

###################
table(hiv_positive$Cluster)
table(hiv_negative$Cluster)

pcr_res<-read.csv("PCR_results.csv", row.names = 1)
#pcr_res<-pcr_res[pcr_res$internal_control=="Positive ", ]

#get participants with more than one STI
allthree_stis <- subset(pcr_res, N_gonorrhea=="Positive" & Chlamydia_trachomatis=="Positive" & Mycoplasma_genitalium=="Positive")
gon_chla <- subset(pcr_res, N_gonorrhea=="Positive" & Chlamydia_trachomatis=="Positive")
gon_myco <- subset(pcr_res, N_gonorrhea=="Positive" & Mycoplasma_genitalium=="Positive")
chla_myco <- subset(pcr_res ,Chlamydia_trachomatis=="Positive" & Mycoplasma_genitalium=="Positive")

pcr_res<-pcr_res[rownames(pcr_res)%in%rownames(meta_data),]
meta_data<-meta_data[rownames(meta_data)%in%rownames(pcr_res),]
meta_sti_new<- cbind(pcr_res, meta_data)
meta_sti_new<-meta_sti_new[c("N_gonorrhea","Chlamydia_trachomatis","Mycoplasma_genitalium",
                             "Cluster","labhivtst", "wksameno2","labhsv2tst")]

#categorising clusters by STI status
table(meta_sti_new$Cluster, meta_sti_new$N_gonorrhea)
table(meta_sti_new$Cluster, meta_sti_new$Chlamydia_trachomatis)
table(meta_sti_new$Cluster, meta_sti_new$Mycoplasma_genitalium)

#categorising STI status by HIV status
table(meta_sti_new$labhivtst, meta_sti_new$N_gonorrhea)
table(meta_sti_new$labhivtst, meta_sti_new$Chlamydia_trachomatis)
table(meta_sti_new$labhivtst, meta_sti_new$Mycoplasma_genitalium)

#compare number of individuals in each cluster to PCR STI results
fisher.test(meta_sti_new$Cluster, meta_sti_new$N_gonorrhea)
tmp<-table(meta_sti_new$Cluster, meta_sti_new$N_gonorrhea)# meta_sti[c("Cluster","labhivtst")]
PT<-pairwiseNominalIndependence(tmp,fisher = TRUE,gtest = F, chisq=F, digits=3)
PT

fisher.test(meta_sti_new$Cluster, meta_sti_new$Chlamydia_trachomatis)
fisher.test(meta_sti_new$Cluster, meta_sti_new$Mycoplasma_genitalium)
combi_trans <- combined
combi_trans$rs3025039 <- as.factor(ifelse(combi_trans$Polymorphism == 1 & combi_trans$Genotype == "homo", "homo",
                                       ifelse(combi_trans$Polymorphism == 1 & combi_trans$Genotype == "het", "het",
                                              ifelse(combi_trans$Polymorphism == 4, "het", 0))))

combi_trans$Polymorphism <- as.factor(ifelse(combi_trans$Polymorphism == 1 | combi_trans$Polymorphism == 4, 1, 0))


## get the counts
rs30_dat <- combi_trans %>%  
  group_by(Preeclampsia_Status, Polymorphism, rs3025039) %>% 
  summarise(Total = n())


## homozygote OR
rs30_hom <- matrix(c(0, 110, 1, 115), ncol = 2, byrow = F)


## pvalues for homozygote
rs30_odd_h <- oddsratio(rs30_hom)
rs30_odd_hp <- rs30_odd_h$p.value
rs30_odd_hconf <- rs30_odd_h$conf.int
rs30_odd_hest <- rs30_odd_h$estimate



## pvalue for heterozygote OR genotype
rs30_het <- matrix(c(15, 110, 9, 115), ncol = 2, byrow = F)
rs30_odd_ht <- oddsratio(rs30_het)
rs30_odd_htp <- rs30_odd_ht$p.value
rs30_odd_htconf <- rs30_odd_ht$conf.int
rs30_odd_htest <- rs30_odd_ht$estimate


## pvalue for the allele
rs30_all <- matrix(c(15, 235, 11, 239), ncol = 2, byrow =F)
rs30_odd_hta <- oddsratio(rs30_all)
rs30_odd_htpa <- rs30_odd_hta$p.value
rs30_odd_htconfa <- rs30_odd_hta$conf.int
rs30_odd_htesta <- rs30_odd_hta$estimate

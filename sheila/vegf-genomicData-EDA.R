# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script is for downstream analysis of VEGF genomic data                             #
# with associated sample/participant data.                                                #
# Author: Alfred Ssekagiri (assekagiri@gmail.com)                                         #
# Date created: 20th April 2022                                                           #
# Date last modified: 15th May 2022                                                       #
# Link: https://github.com/galaxyuvri-ea/Projects/blob/master/sheila/vegf-genomicData.R   #
# # # # # # # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # # # # # # # # # #  #

# clear environment
rm(list = ls())

# load required packages 
require(reshape)
require(ggplot2)
require(cowplot)

# load data into R
GenomicData <- read.csv("~/Desktop/MISC/sheila_updated.csv")
# create a dataframe with varibles of interest
tmp_vars<-c("No","Gestational_age","Maternal_age", "VEGF","Status", "Genotype")
tmp_df<-GenomicData[tmp_vars]
# obtain the mean gestation and maternal ages for cases and controls
aggregate(Gestational_age~Status, tmp_df, FUN=mean)
aggregate(Maternal_age~Status, tmp_df, FUN=mean)
#create seperate dataframes for cases and controls 
cases<-GenomicData[GenomicData$Status=="Cases",]
controls<-GenomicData[GenomicData$Status=="Controls",]
# perform a normality test on the distribution of maternal and gestation age for cases
shapiro.test(cases$Gestational_age)
shapiro.test(cases$Maternal_age)
# perform a normality test on the distribution of maternal and gestation age for controls
shapiro.test(controls$Gestational_age)
shapiro.test(controls$Maternal_age)
# peform wilcoxon test for difference in medians of maternal and gestational ages between cases and controls 
wilcox.test(Gestational_age~Status,data = tmp_df)
wilcox.test(Maternal_age~Status,data = tmp_df)
# obtain the quantiles of gestation and maternal ages for cases and controls
aggregate(Gestational_age~Status, tmp_df, FUN=quantile)
aggregate(Maternal_age~Status, tmp_df, FUN=quantile)
# create a long representation of tmp_df in preparation for visualing using ggplot2
df <- melt(tmp_df, id=c("Status","No","Genotype"), variable_name = "variable")

# Firstly, we visualize gestational age in cases and controls
df1<-df[df$variable=="Gestational_age",]
panelA<-ggplot(data=df1, aes(x=Status,y=value,color=Status))+
  geom_boxplot(na.rm = T)+geom_point()+
  geom_jitter(height = 0.5)+
  scale_color_manual("Status",values=c("#728FCE","#EBCC2A")) +
  ylab('Gestation age (Weeks)')+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA))+xlab("Preeclampsia status")+
  annotate(geom = "text", x=1.4, y = 45,label="p-value=0.5294",size=5)
panelA

# Secondly, we visualize maternal age in cases and controls
df2<-df[df$variable=="Maternal_age",]
panelB<-ggplot(data=df2, aes(x=Status,y=value,color=Status))+
  geom_boxplot(na.rm = T)+geom_point()+
  geom_jitter(height = 0.5)+
  scale_color_manual("Status",values=c("#728FCE","#EBCC2A")) +
  ylab('Maternal age (Years)')+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA))+xlab("Preeclampsia status")+
  annotate(geom = "text", x=1.4, y = 40,label="p-value=0.9923",size=5)
panelB

# obtain VEGF quantiles and compare the VEGF levels between cases and controls
df3<-df[df$variable=="VEGF",]
df3<-na.omit(df3)
aggregate(value~Status, data=df3, FUN = quantile)
kruskal.test(value~Status, data=df3)
# for visualization, lets do a log-2 transform of the data
df3$value<-log2(df3$value)
panelC<-ggplot(data=df3, aes(x=Status,y=value,color=Status))+
  geom_boxplot(na.rm = T)+geom_point()+#facet_wrap(~variable, scales = "free")
  geom_jitter(height = 0.1)+
  scale_color_manual("Status",values=c("#728FCE","#EBCC2A")) +
  ylab('Log2-transformed  plasma concentration')+ggtitle("All genotypes")+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA))+xlab("Preeclampsia status")+
  annotate(geom = "text", x=1.4, y = 4,label="p-value=0.0616",size=5)
panelC

# get VEGF levels for CC genotype in controls and cases
df3<-df[df$variable=="VEGF",]
df3<-df3[df3$Genotype=="CC",]
df3<-na.omit(df3)
wilcox.test(value~Status, data=df3)
df3$value<-log2(df3$value)
panelD<-ggplot(data=df3, aes(x=Status,y=value,color=Status))+
  geom_boxplot(na.rm = T)+geom_point()+#facet_wrap(~variable, scales = "free")
  geom_jitter(height = 0.1)+
  scale_color_manual("Status",values=c("#728FCE","#EBCC2A")) +
  ylab('Log2-transformed plasma concentration')+ggtitle("Genotype CC")+#ggtitle("Lineages and mutations")+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA))+xlab("Preeclampsia status")+
  annotate(geom = "text", x=1.4, y = 4,label="p-value=0.0437",size=5)

panelD

# get VEGF levels for non CC genotypes in controls and cases
df3<-df[df$variable=="VEGF",]
df3<-df3[df3$Genotype!="CC",]
df3<-na.omit(df3)
kruskal.test(value~Status, data=df3)
df3$value<-log2(df3$value)
#df3<-df3[df3$value<40,]
panelE<-ggplot(data=df3, aes(x=Status,y=value,color=Status))+
  geom_boxplot(na.rm = T)+geom_point()+#facet_wrap(~variable, scales = "free")
  geom_jitter(height = 0.1)+
  scale_color_manual("Status",values=c("#728FCE","#EBCC2A")) +
  ylab('Log-transformed plasma concentration')+ggtitle("Genotypes CT and TT")+#ggtitle("Lineages and mutations")+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA))+xlab("Preeclampsia status")+
  annotate(geom = "text", x=1.4, y = 2,label="p-value=0.7903",size=5)
panelE

# combine panels 
fig1<-cowplot::plot_grid(panelC, panelD, panelE, labels = c("A","B","C"), nrow = 1)
fig1
fig2<-cowplot::plot_grid(panelC, panelD, labels = c("A","B"), nrow = 1)
#ggsave("Figure2.pdf", plot=fig2, device="pdf", height=7, width = 12)
#all<-cowplot::plot_grid(panelA, panelB, panelC, panelD, labels = c("A","B","C", "D"), nrow = 2)
#ggsave("Figure1.pdf", plot=all, device="pdf", height=10, width = 12)

# Detailed analysis of VEGF levels between cases and controls
vegf<-GenomicData[c("VEGF","Genotype","Status")]
shapiro.test(vegf$VEGF)
vegf<-na.omit(vegf)
dim(vegf)
# quantiles of VEGF across genotypes
aggregate(VEGF~Genotype, data = vegf, FUN = quantile)
kruskal.test(VEGF~Genotype, data = vegf)
# VEGF levels for cases with CT genotype
cases_vegf<-vegf[vegf$Genotype=="CT",]
shapiro.test(cases$VEGF)
cases_vegf<-na.omit(cases_vegf)
dim(cases_vegf)
table(cases_vegf$Status)
aggregate(VEGF~Status, data = cases_vegf, FUN = quantile)
wilcox.test(VEGF~Status, data = cases_vegf)
kruskal.test(VEGF~Status, data = cases_vegf)
# maternal age categoris
age_category<-function(x){
  if(x>12 & x< 18){return("13-17")}
  else if(x>17 & x< 23){return("18-22")}
  else if(x>22 & x< 28){return("23-27")}
  else if(x>27 & x< 33){return("28-32")}
  else if(x>32 & x< 38){return("33-37")}
  else if(x>37 & x< 43){return("38-42")}
  else if(x>42 & x< 48){return("43-47")}
}
GenomicData$Maternal_age_category<-sapply(GenomicData$Maternal_age, age_category)
# gestational age categories
gestational_age_category<-function(x){
  if(x>19 & x< 29){return("20-28")}
  else if(x>28 & x< 34){return("29-33")}
  else if(x>33 & x< 38){return("34-37")}
  else if(x>37 & x< 44){return("38-43")}
}
GenomicData$gestational_age_category<-sapply(GenomicData$Gestational_age, gestational_age_category)
GenomicData2<-GenomicData[!colnames(GenomicData)%in%c("No","Gestational_age","Maternal_age","VEGF", "STUDY_ID","tribe.mother", "tribe.father")]
GenomicData2<-data.frame(sapply(GenomicData2, as.character))
# obtain counts for each categorical variables and compare them by PE status
for (i in 1:dim(GenomicData2)[2]) {
  print(paste("===================",names(GenomicData2)[i], "======================"))
  m<-table(GenomicData2[,i], GenomicData2$Status)
  print(m)
  print(chisq.test(m))
}
GenomicData2$No<-GenomicData$No
# change PE levels to 0 and 1 for controls and 0 respectively, this is required for clogit later 
GenomicData2$case[GenomicData2$Status=="Cases"]<-1
GenomicData2$case[GenomicData2$Status=="Controls"]<-0
# set ref levels for hiv status and rs3025039
GenomicData2$rs3025039<-relevel(factor(GenomicData2$rs3025039, levels = c(1,2), label=c("Present","Absent")), ref='Absent')
GenomicData2$hiv_status<-relevel(factor(GenomicData2$hiv_status), ref='2')
# function to run clogit over a number of variables, one at a time
results_CLR <- function(data, var) {
  event<-data[,"case"]
  id<-data[,"No"]
  var_<-data[,var]
  mod<-clogit(event~var_ + strata(id))
  print(var)
  print(summary(mod))
}

# get all variable names
all_vars<-names(GenomicData2)
# get only vars on interests from above list
all_vars<-all_vars[!all_vars%in%c("No","case","STUDY_ID","Status","Maternal_age_category","gestational_age_category")]
all_vars
# run clogit over all variables of interest, one at a time
for(var_ in all_vars){
  results_CLR(GenomicData2, var_)
}
# run clogit on multiple variables at a ago
mod<-clogit(case~rs3025039+
              Family.history.of.preeclampsia+
              Family.history.of.hypertension+
              Family.history.of.diabetes.Mellitus+
              first.pregnancy+
              #diagnosis.with.hypertension.in.previous.pregnancy+
              type.of.pregnancy+
              hiv_status+
              strata(No), data = GenomicData2)
summary(mod)

# generate plots for variant frequency by PE status
vars2<-c("Status","rs3025039", "Genotype")
plot_df<-GenomicData[vars2]
plot_df$Number<-1
plot_df2<-aggregate(Number~., data=plot_df, FUN=sum)
plot_df2$Genotype2[plot_df2$Genotype%in%c("CC","TT")]<-"Homozygous"
plot_df2$Genotype2[plot_df2$Genotype%in%c("CT")]<-"Heterozygous"
plot_df2$variant<-factor(plot_df2$rs3025039,
                           levels = c(1,2),
                           labels = c("rs3025039","Wild"))
p<-ggplot(data=plot_df2, aes(x=Status,y=Number,fill=variant))
p<-p+geom_bar(stat = "identity",position=position_dodge())
p<-p+geom_text(aes(label=Number),position=position_dodge(width=0.9), vjust=-0.25,size=5)
p<-p+facet_wrap(~Genotype2)
p$labels$fill<-"Polymorphism"
p<-p+theme(axis.title.x = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=15))+
  theme(strip.background = element_rect(fill="white",colour = "black"),
        strip.text = element_text(size = 15, face = "bold", colour = "black"))+
  theme(legend.position="right", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        legend.text = element_text(size=12, colour = "black"),
        legend.title = element_text(size=15, colour = "black", face = "bold"),
        panel.border = element_rect(fill=NA))+xlab("Preeclampsia status")+
  ylab("Number of participants")
p
#ggsave("~/Desktop/MISC/Figure1.tiff", plot=p, device="tiff", height=08, width = 10, dpi = 300)
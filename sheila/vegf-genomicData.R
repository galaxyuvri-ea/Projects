rm(list = ls())
GenomicData <- read.csv("~/Desktop/MISC/sheila_updated.csv")

#################
numeric_vars<-c("No","Gestational_age","Maternal_age", "VEGF","Status", "Genotype")
age_df<-GenomicData[numeric_vars]
aggregate(Gestational_age~Status, age_df, FUN=mean)
aggregate(Maternal_age~Status, age_df, FUN=mean)

cases<-GenomicData[GenomicData$Status=="Cases",]
controls<-GenomicData[GenomicData$Status=="Controls",]

shapiro.test(cases$Gestational_age)
shapiro.test(cases$Maternal_age)

shapiro.test(controls$Gestational_age)
shapiro.test(controls$Maternal_age)

wilcox.test(Gestational_age~Status,data = age_df)
wilcox.test(Maternal_age~Status,data = age_df)

aggregate(Gestational_age~Status, age_df, FUN=quantile)
aggregate(Maternal_age~Status, age_df, FUN=quantile)


require(reshape)
#age_df$VEGF<-log10(age_df$VEGF)
df <- melt(age_df, id=c("Status","No","Genotype"), variable_name = "variable")
require(ggplot2)
#vegf<-vegf[vegf$VEGF<10,]
df1<-df[df$variable=="Gestational_age",]
panelA<-ggplot(data=df1, aes(x=Status,y=value,color=Status))+
  geom_boxplot(na.rm = T)+geom_point()+#facet_wrap(~variable, scales = "free")
  geom_jitter(height = 0.5)+
  scale_color_manual("Status",values=c("#728FCE","#EBCC2A")) +
  ylab('Gestation age (Weeks)')+#ggtitle("Lineages and mutations")+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA))+xlab("Preeclampsia status")+
  annotate(geom = "text", x=1.4, y = 45,label="p-value=0.5294",size=5)


panelA

df2<-df[df$variable=="Maternal_age",]
panelB<-ggplot(data=df2, aes(x=Status,y=value,color=Status))+
  geom_boxplot(na.rm = T)+geom_point()+#facet_wrap(~variable, scales = "free")
  geom_jitter(height = 0.5)+
  scale_color_manual("Status",values=c("#728FCE","#EBCC2A")) +
  ylab('Maternal age (Years)')+#ggtitle("Genotypes CC and CT")+
  theme(axis.title.x = element_text(color="black", size=15))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(legend.position="none", axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        panel.border = element_rect(fill=NA))+xlab("Preeclampsia status")+
  annotate(geom = "text", x=1.4, y = 40,label="p-value=0.9923",size=5)


panelB

df3<-df[df$variable=="VEGF",]
df3<-na.omit(df3)
aggregate(value~Status, data=df3, FUN = quantile)
kruskal.test(value~Status, data=df3)
df3$value<-log2(df3$value)
df3<-df3[df3$value<40,]
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

df3<-df[df$variable=="VEGF",]
df3<-df3[df$Genotype=="CC",]
df3<-na.omit(df3)
wilcox.test(value~Status, data=df3)
df3$value<-log2(df3$value)
#df3<-df3[df3$value<40,]
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

df3<-df[df$variable=="VEGF",]
df3<-df3[df$Genotype!="CC",]
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

require(cowplot)
fig1<-cowplot::plot_grid(panelC, panelD, panelE, labels = c("A","B","C"), nrow = 1)
fig1
fig2<-cowplot::plot_grid(panelC, panelD, labels = c("A","B"), nrow = 1)
ggsave("~/Desktop/MISC/SheilaFigure2-updated.pdf", plot=fig1, device="pdf", height=7, width = 12)
#ggsave("~/Desktop/MISC/SheilaFigure2.pdf", plot=fig2, device="pdf", height=7, width = 12)

#all<-cowplot::plot_grid(panelA, panelB, panelC, panelD, labels = c("A","B","C", "D"), nrow = 2)
#all
#ggsave("~/Desktop/MISC/SheilaFigure1_updated.pdf", plot=all, device="pdf", height=10, width = 12)
################### VEGF analysis #############
vegf<-GenomicData[c("VEGF","Genotype","Status")]
shapiro.test(vegf$VEGF)
vegf<-na.omit(vegf)
dim(vegf)
aggregate(VEGF~Genotype, data = vegf, FUN = quantile)
kruskal.test(VEGF~Genotype, data = vegf)

cases_vegf<-vegf[vegf$Genotype=="CT",]
shapiro.test(cases$VEGF)
cases_vegf<-na.omit(cases_vegf)
dim(cases_vegf)
table(cases_vegf$Status)
aggregate(VEGF~Status, data = cases_vegf, FUN = quantile)
wilcox.test(VEGF~Status, data = cases_vegf)
kruskal.test(VEGF~Status, data = cases_vegf)

#################
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

gestational_age_category<-function(x){
  if(x>19 & x< 29){return("20-28")}
  else if(x>28 & x< 34){return("29-33")}
  else if(x>33 & x< 38){return("34-37")}
  else if(x>37 & x< 44){return("38-43")}
}
GenomicData$gestational_age_category<-sapply(GenomicData$Gestational_age, gestational_age_category)

GenomicData2<-GenomicData[!colnames(GenomicData)%in%c("No","Gestational_age","Maternal_age","VEGF", "STUDY_ID","tribe.mother", "tribe.father")]
GenomicData2<-data.frame(sapply(GenomicData2, as.character))


for (i in 1:dim(GenomicData2)[2]) {
  print(paste("===================",names(GenomicData2)[i], "======================"))
  m<-table(GenomicData2[,i], GenomicData2$Status)
  print(m)
  print(chisq.test(m))
}

GenomicData2$No<-GenomicData$No

GenomicData2$case[GenomicData2$Status=="Cases"]<-1
GenomicData2$case[GenomicData2$Status=="Controls"]<-0

GenomicData2$rs3025039<-relevel(factor(GenomicData2$rs3025039, levels = c(1,2), label=c("Present","Absent")), ref='Absent')
GenomicData2$hiv_status<-relevel(factor(GenomicData2$hiv_status), ref='2')

require(survival)

results_CLR <- function(data, var) {
  event<-data[,"case"]
  id<-data[,"No"]
  var_<-data[,var]
  mod<-clogit(event~var_ + strata(id))
  print(var)
  print(summary(mod))
}

ab<-names(GenomicData2)
ab<-ab[!ab%in%c("No","case","STUDY_ID","Status","Maternal_age_category","gestational_age_category")]
ab

for(var_ in ab){
  results_CLR(GenomicData2, var_)
}

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


genotype_mat <- table(GenomicData$Genotype, GenomicData$Status)
genotype_mat
fisher.test(genotype_mat)
allelic_mat <- matrix(data=c(235, 239, 15, 11), nrow = 2, ncol = 2, byrow = T, dimnames = list(c("C","T")))
allelic_mat
fisher.test(allelic_mat)

gen_mat <- matrix(data=c(110, 115, 15, 9), nrow = 2, ncol = 2, byrow = T, dimnames = list(c("CC","CT")))
gen_mat
fisher.test(gen_mat)

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
p<-p+geom_text(aes(label=Number),
               position=position_dodge(width=0.9), vjust=-0.25,
               size=5)
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

ggsave("~/Desktop/MISC/Figure1.tiff", plot=p, device="tiff", height=08, width = 10, dpi = 300)

rm(list = ls())
stella <- read.csv("stella.csv")

require(ggpubr)

names(stella)

stella$meanod1<-stella$meanod
stella[7:33]<-sapply(stella[7:33], as.character)

stella$B_temp[which(stella$B_temp<36)]<-36

stella$Temp_ranges <- cut(stella$B_temp, breaks = c(35, 37.5, 38, 38.5, 39), labels = c("<37.6","37.6-38.0","38.1-38.5","38.6-39.0"))

stella$Age_groups <- cut(stella$Age, breaks = c(7, 14, 24, 34, 44, 54, 84), labels = c("<14", "15-24", "25-34", "35-44", "45-54", ">54"))

###### descriptive stats #######
stella$Interpret<-as.character(factor(stella$Interpret,levels=c("1","2"), labels=c("Positive","Negative")))
stella$Sex<-as.character(factor(stella$Sex,levels=c("1","2"), labels=c("Male","Female")))
table(stella$status, stella$Sex)
table(stella$status, stella$Age_groups)
#table(stella$Subcounty)

#meanod by gender
ylim1 = boxplot.stats(stella$meanod1)$stats[c(1, 5)]
ylim1
groups<-list(c("Male","Female"))#,c("C2","C3"), c("C3","C4"))#,c("C4","C1"), c("C3","C1"), c("C4","C2")
p<-ggboxplot(stella, x="Sex", y="meanod1", color="Sex", palette = c("blue","red"), add="jitter", 
             ylim=c(0,0.4), xlab = "")
p + stat_compare_means(label.y = 0.4)+ylab("MeanOD (nm)") # Add global p-value

#meanod by status
groups<-list(c("Controls","Farmers"))#,c("C2","C3"), c("C3","C4"))#,c("C4","C1"), c("C3","C1"), c("C4","C2")
stella$Status<-stella$status
p<-ggboxplot(stella, x="Status", y="meanod1", color="Status", palette = c("blue","red"), add="jitter", 
             ylim=c(0,0.4), xlab = "")
p + stat_compare_means(label.y = 0.4)+ylab("MeanOD (nm)") # Add global p-value
p$labels$color<-"Status"

#meanod by age groups
groups<-list(c("<14","15-24","25-34","35-44","45-54",">54"))#,c("C2","C3"), c("C3","C4"))#,c("C4","C1"), c("C3","C1"), c("C4","C2")
stella$Age<-stella$Age_groups
p<-ggboxplot(stella, x="Age", y="meanod1", color="Age", add="jitter", 
             ylim=c(0,0.5), xlab = "")
p + stat_compare_means(label.y = 0.5)+ylab("MeanOD (nm)") # Add global p-value
#sero prevalence by status
m<-table(stella$status,stella$Interpret)
fisher.test(m)

m<-table(stella$Sex, stella$Interpret)
fisher.test(m)

m<-table(stella$Age_groups, stella$Interpret)
fisher.test(m, simulate.p.value=TRUE)

#sero prevalence and risk of participants
m<-table(stella$Rearinglivestock, stella$Interpret)
fisher.test(m)
m<-table(stella$Tickcrushing, stella$Interpret)
fisher.test(m)
m<-table(stella$Wlife_exp, stella$Interpret)
fisher.test(m)
m<-table(stella$Temp_ranges, stella$Interpret)
fisher.test(m)
m<-table(stella$Petcaring, stella$Interpret)
fisher.test(m)
m<-table(stella$Bitingflies, stella$Petcaring)
m
fisher.test(m)
m<-table(stella$Wlife_bites, stella$Petcaring)
m
fisher.test(m)

m<-table(stella$Ticks, stella$Petcaring)
m
fisher.test(m)

#Graphics
m<-table(stella$Interpret, stella$status)
m
fisher.test(m)
m<-as.data.frame(m)
require(ggplot2)
p<-ggplot(m, aes(x = Var2, y = Freq))+
  geom_bar(
    aes(fill = Var2), stat = "identity", color = "white",
    position = position_dodge(0.9)
  )+
  facet_wrap(~Var1)+theme_bw()+ylab("Number of Participants")+xlab("")
p$labels$fill<-""
p


m<-table(stella$Interpret, stella$Sex)
m
fisher.test(m)
m<-as.data.frame(m)
require(ggplot2)
p<-ggplot(m, aes(x = Var2, y = Freq))+
  geom_bar(
    aes(fill = Var2), stat = "identity", color = "white",
    position = position_dodge(0.9)
  )+
  facet_wrap(~Var1)+theme_bw()+ylab("Number of Participants")+xlab("")
p$labels$fill<-"Gender"
p

m<-table(stella$Interpret, stella$Age_groups)
m
fisher.test(m, simulate.p.value = TRUE)
m<-as.data.frame(m)
p<-ggplot(m, aes(x = Var2, y = Freq))+
  geom_bar(
    aes(fill = Var2), stat = "identity", color = "white",
    position = position_dodge(0.9)
  )+
  facet_wrap(~Var1)+theme_bw()+ylab("Number of Participants")+xlab("")
p$labels$fill<-"Age groups"
p

stella$Med_history<-as.character(factor(stella$Med_history,levels=c("1","2"), labels=c("Yes","No")))
m<-table(stella$Interpret, stella$Med_history)
m
fisher.test(m)
m<-as.data.frame(m)

p<-ggplot(m, aes(x = Var2, y = Freq))+
  geom_bar(
    aes(fill = Var2), stat = "identity", color = "white",
    position = position_dodge(0.9)
  )+
  facet_wrap(~Var1)+theme_bw()+ylab("Number of Participants")+xlab("")
p$labels$fill<-"Medical history"
p

#Risks
stella$Butch_slaughter<-as.character(factor(stella$Butch_slaughter,levels=c("1","2"), labels=c("Yes","No")))
m<-table(stella$Interpret, stella$Butch_slaughter)
m
fisher.test(m)
m<-as.data.frame(m)
p<-ggplot(m, aes(x = Var2, y = Freq))+
  geom_bar(
    aes(fill = Var2), stat = "identity", color = "white",
    position = position_dodge(0.9)
  )+
  facet_wrap(~Var1)+theme_bw()+ylab("Number of Participants")+xlab("")
p$labels$fill<-"Buther/Slaughter"
p

stella$Rearinglivestock<-as.character(factor(stella$Rearinglivestock,levels=c("1","2"), labels=c("Yes","No")))
m<-table(stella$Interpret, stella$Rearinglivestock)
m
fisher.test(m)
m<-as.data.frame(m)
p<-ggplot(m, aes(x = Var2, y = Freq))+
  geom_bar(
    aes(fill = Var2), stat = "identity", color = "white",
    position = position_dodge(0.9)
  )+
  facet_wrap(~Var1)+theme_bw()+ylab("Number of Participants")+xlab("")
p$labels$fill<-"Livestock rearing"
p

stella$Tickcrushing<-as.character(factor(stella$Tickcrushing,levels=c("1","2"), labels=c("Yes","No")))
m<-table(stella$Interpret, stella$Tickcrushing)
m
fisher.test(m)
m<-as.data.frame(m)
p<-ggplot(m, aes(x = Var2, y = Freq))+
  geom_bar(
    aes(fill = Var2), stat = "identity", color = "white",
    position = position_dodge(0.9)
  )+
  facet_wrap(~Var1)+theme_bw()+ylab("Number of Participants")+xlab("")
p$labels$fill<-"Tick crushing"
p

stella$Petcaring<-as.character(factor(stella$Petcaring,levels=c("1","2"), labels=c("Yes","No")))
m<-table(stella$Interpret, stella$Petcaring)
m
fisher.test(m)
m<-as.data.frame(m)
p<-ggplot(m, aes(x = Var2, y = Freq))+
  geom_bar(
    aes(fill = Var2), stat = "identity", color = "white",
    position = position_dodge(0.9)
  )+
  facet_wrap(~Var1)+theme_bw()+ylab("Number of Participants")+xlab("")
p$labels$fill<-"Pet caring"
p

stella$Wlife_exp<-as.character(factor(stella$Wlife_exp,levels=c("1","2"), labels=c("Yes","No")))
m<-table(stella$Interpret, stella$Wlife_exp)
m
fisher.test(m)
m<-as.data.frame(m)
p<-ggplot(m, aes(x = Var2, y = Freq))+
  geom_bar(
    aes(fill = Var2), stat = "identity", color = "white",
    position = position_dodge(0.9)
  )+
  facet_wrap(~Var1)+theme_bw()+ylab("Number of Participants")+xlab("")
p$labels$fill<-"Wild life exposure"
p
ylim1 = boxplot.stats(stella$B_temp)$stats[c(1, 5)]
ylim1
stella$Seropositivity<-stella$Interpret
groups<-list(c("Negative","Positive"))#,c("C2","C3"), c("C3","C4"))#,c("C4","C1"), c("C3","C1"), c("C4","C2")
p<-ggboxplot(stella, x="Seropositivity", y="B_temp", color="Seropositivity", add="jitter", 
             ylim=c(35.5, 40),xlab = "")
p + stat_compare_means(label.y=40)+ylab("Body temperature (Degrees centigrade)") # Add global p-value
# m<-table(stella$Subcounty)
# m<-as.data.frame(m)
# p<-ggplot(m, aes(x = Var1, y = Freq))+
#   geom_bar(
#     aes(fill = Var1), stat = "identity", color = "white",
#     position = position_dodge(0.9)
#   )+theme_bw()+ylab("Number of Participants")+xlab("")
# p$labels$fill<-"Sub-counties"
# p
ylim1 = boxplot.stats(stella$meanod)$stats[c(1, 5)]
ylim1
stella$Seropositivity<-stella$Interpret
groups<-list(c("Controls","Farmers"))#,c("C2","C3"), c("C3","C4"))#,c("C4","C1"), c("C3","C1"), c("C4","C2")
p<-ggboxplot(stella, x="status", y="meanod", color="status", add="jitter", 
             xlab = "")
p + stat_compare_means(label.y=40)+ylab("MeanOD (nm)") # Add global p-value


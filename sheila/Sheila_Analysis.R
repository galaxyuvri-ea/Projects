rm(list = ls())
dev.off()
sheila<-read.csv("clean_data/PE_labdata.csv")
sheila$sFlt1_P1GF<-sheila$sFlt1/sheila$PlGF
#View(sheila)
#shapiro.test(sheila$sFlt1[sheila$Status=='Controls'])
controls<-sheila[sheila$Status=='Controls',]
cases<-sheila[sheila$Status=='Cases',]

wilcox.test(controls$sFlt1,cases$sFlt1)
wilcox.test(controls$PlGF,cases$PlGF)
wilcox.test(controls$VEGF,cases$VEGF)
wilcox.test(controls$sFlt1_P1GF,cases$sFlt1_P1GF)

#get quantiles of the analytes in ctls and cases
quantile(controls$PlGF)
quantile(controls$sFlt1)
quantile(controls$VEGF)
quantile(controls$sFlt1_P1GF)

quantile(cases$PlGF)
quantile(cases$sFlt1)
quantile(cases$VEGF)
quantile(cases$sFlt1_P1GF)

sheilaOrig<-sheila
sheila$sFlt1<-log(sheila$sFlt1)
sheila$PlGF<-log(sheila$PlGF)
sheila$VEGF<-log(sheila$VEGF)
sheila$sFlt1_P1GF<-log(sheila$sFlt1_P1GF)
######### create facet plots #########
df<-sheila[c(2,4:8)]
names(df)<-c("studyID","sFlT1","PlGF","VEGF","Status","sFlT1/PlGF")
require(reshape2)
df1<-melt(df, id.vars=c("studyID","Status"))
p <- ggplot(df1, aes(Status, value, colour = factor(Status))) + geom_boxplot()+geom_jitter()
p<-p + facet_grid(cols = vars(variable))+theme_classic()
p<-p+xlab('')+ylab('Plasma levels')
p$labels$colour<-"Status"
p

############ conditional logistic regression ########
clinic<-read.csv("clean_data/PE_demodata.csv")
screenGA<-read.csv('clean_data/PE_screenGA.csv')
screen_ids<-read.csv('clean_data/screen_recruit.csv')
#add study id to GA by matching ids
screenGA$studyID<-screen_ids$studyid[match(screenGA$screening_no, screen_ids$screeningn)]
#Let us start by defining age groups and wksameno
#13-17, 18-22, 23-27, 28-32, 33-37, 38-42, 43-47, 48-52 years and 
#20-28, 29-33, 34-37, 38-42 weeks for gestational age
clinic$agegroup<-cut(clinic$age, breaks = c(12, 17, 22, 27, 32, 37, 42, 47, 52), 
                labels = c("13-17", "18-22", "23-27", "28-32", "33-37", "38-42", "43-47", "48-52"))
clinic$GA<-screenGA$Gestational_age2[match(clinic$studyid, screenGA$studyID)]

clinic$wksameno<-cut(clinic$GA, 
                     breaks = c(19, 28, 33, 37, 43), 
                     labels = c("20-28", "29-33", "34-37", "38-43"))
#make study id similar to the one in the labdata
#Just add a prefix, PE for Pre-Eclampsia
clinic$studyID<-paste0("PE",clinic$studyid)
clinic$pregnancy7months[is.na(clinic$pregnancy7months)]<-0
clinic$pregnancy7months<-cut(clinic$pregnancy7months, 
                             breaks = c(-1,0,2,8),
                             labels=c('None','1-2','Above 2'))

hypert<-read.csv('clean_data/dxhypertension.csv')
hypert$studyid<-paste0('PE',hypert$studyid)
hypert$status<-NULL
names(hypert)<-c("studyID","dxhypertension")
clinic<-merge(clinic, hypert, by='studyID')
#make all characters and then factors 
library(dplyr)
clinic <- clinic %>%
  mutate_all(as.character)
#str(clinic)
clinic <- clinic %>%
  mutate_all(as.factor)
str(clinic)

#Merge clinical and lab data
merged<-merge(sheilaOrig, clinic, by="studyID",all.x = T)

require(survival)
merged$case[merged$Status.x=="Cases"]<-1
merged$case[merged$Status.x=="Controls"]<-0
rownames(merged)<-merged$studyID

#compare clinical variables btn cases and controls
m<-table(merged$dxhypertension, merged$case)
m
fisher.test(m)

m<-table(merged$hivstatus, merged$case)
m
fisher.test(m)

m<-table(merged$historypreeclampsia, merged$case)
m
fisher.test(m)

m<-table(merged$historyhypertension, merged$case)
m
fisher.test(m)

m<-table(merged$firstpregnancy, merged$case)
m
fisher.test(m)

m<-table(merged$condomuse, merged$case)
m
fisher.test(m)

m<-table(merged$historystroke, merged$case)
m
fisher.test(m)

m<-table(merged$diabetesmellitus, merged$case)
m
fisher.test(m)

m<-table(merged$heartattack, merged$case)
m
fisher.test(m)

m<-table(merged$differentpartner, merged$case)
m
fisher.test(m)

m<-table(merged$pregnancy7months,merged$case)
m
fisher.test(m)

############# compare demographics ###########
m<-table(merged$agegroup, merged$case)
m
fisher.test(m)

m<-table(merged$wksameno, merged$case)
m
fisher.test(m)

m<-table(merged$education, merged$case)
m
fisher.test(m)

m<-table(merged$religion, merged$case)
m
fisher.test(m)

m<-table(merged$maritalstatus, merged$case)
m
fisher.test(m)

m<-table(merged$alcohol, merged$case)
m
fisher.test(m) 

m<-table(merged$smoking, merged$case)
m
fisher.test(m) 

#univariate conditional logistic regression analysis
modData<-merged[ , -which(names(merged) %in% c("Status.x","Status.y","age","Sample",
                                               "studyID","studyid.x","studyid.y",
                                               "gestationagewks"))]

# modData$sFlt1<-log2(modData$sFlt1)
# modData$PlGF<-log2(modData$PlGF)
# modData$VEGF<-log2(modData$VEGF)
# modData$sFlt1_P1GF<-log2(modData$sFlt1_P1GF)

results_CLR <- function(data, var) {
  event<-data[,"case"]
  id<-data[,"order"]
  var_<-data[,var]
  mod<-clogit(event~var_ + strata(id))
  print(var)
  print(summary(mod))
}

ab<-names(modData)
ab<-ab[!ab%in%c("order","case")]

for(var_ in ab){
  results_CLR(modData, var_)
}

########## Multivariate conditional logistic regression ##########
modall<-clogit(case~VEGF+sFlt1+PlGF+dxhypertension + strata(order), data = modData)
summary(modall)
modData$pregnancy7months<-relevel(modData$pregnancy7months, ref='None')
allsig<-clogit(case~VEGF+sFlt1+PlGF+pregnancy7months+
                 historypreeclampsia+historyhypertension+
                 hivstatus+condomuse+strata(order), data=modData)
summary(allsig)

############## make plots of plasma levels
#stratified bu GA and PE status ###############
rm(sheila)
sheila<-merged
sheila$sFlt1<-log(sheila$sFlt1)
sheila$PlGF<-log(sheila$PlGF)
sheila$VEGF<-log(sheila$VEGF)
sheila$sFlt1_P1GF<-log(sheila$sFlt1_P1GF)
sheila<-sheila[!is.na(sheila$wksameno),]
sheila$Status[sheila$case==1]<-"Cases"
sheila$Status[sheila$case==0]<-"Controls"

######### get other plots #####
df2<-sheila[c("studyID", "Status", "sFlt1", "PlGF", "VEGF", "sFlt1_P1GF", "wksameno")]
names(df2)<-c("studyID","Status", "sFlT1", "PlGF","VEGF", "sFlT1/PlGF","wksameno")
df3<-melt(df2, id.vars=c("studyID","Status","wksameno"))
p <- ggplot(df3, aes(wksameno, value, colour = factor(Status))) + geom_boxplot()+  geom_point(position = position_jitterdodge())
p<-p + facet_grid(cols = vars(variable))+theme_classic()
p<-p+xlab('')+ylab('Plasma levels')
p$labels$colour<-"Status"
p

########### Heatmap ###########
ann_col<-df2[c("Status","wksameno")]
plasma<-df2[3:6]

#require(pheatmap)
#pheatmap(plasma, annotation_col = ann_col)
ph<-ggplot(data = df3,
           aes(x=studyID,y=variable,fill=value))
ph<-ph+geom_tile()+facet_grid(~wksameno+Status,scales = "free_x", space = "free_x")+scale_fill_gradient(name = "Plasma levels",
                                                                                                        low = "#FFFFFF",
                                                                                                        high = "#012345")
ph<-ph+ylab('')+xlab('Study participants')
ph<-ph+theme_classic()+theme(axis.text.x=element_blank(),
             axis.ticks.x=element_blank())
ph

ph<-ggplot(data = df3,
           aes(x=studyID,y=variable,fill=value))
ph<-ph+geom_tile()+facet_grid(~Status,scales = "free_x", space = "free_x")+scale_fill_gradient(name = "Plasma levels",
                                                                                                        low = "#FFFFFF",
                                                                                                        high = "#012345")
ph<-ph+ylab('')+xlab('Study participants')
ph<-ph+theme_classic()+theme(axis.text.x=element_blank(),
                             axis.ticks.x=element_blank())
ph
############split into wksameno grps############
rm(sheila)
sheila<-merged
sheila<-sheila[!is.na(sheila$wksameno),]

sheila$Status[sheila$case==1]<-"Cases"
sheila$Status[sheila$case==0]<-"Controls"

first_wk<-sheila[sheila$wksameno=="20-28",]
second_wk<-sheila[sheila$wksameno=="29-33",]
third_wk<-sheila[sheila$wksameno=="34-37",]
fourth_wk<-sheila[sheila$wksameno=="38-43",]
####compare plasma concn btn cases and ctrls by gestational age

#start by comparing PlGF
wilcox.test(PlGF~Status, data = first_wk)
wilcox.test(PlGF~Status, data = second_wk)
wilcox.test(PlGF~Status, data = third_wk)
wilcox.test(PlGF~Status, data = fourth_wk)
aggregate(PlGF~Status, data = first_wk, quantile)
aggregate(PlGF~Status, data = second_wk, quantile)
aggregate(PlGF~Status, data = third_wk, quantile)
aggregate(PlGF~Status, data = fourth_wk, quantile)
#VEGF
wilcox.test(VEGF~Status, data = first_wk)
wilcox.test(VEGF~Status, data = second_wk)
wilcox.test(VEGF~Status, data = third_wk)
wilcox.test(VEGF~Status, data = fourth_wk)
aggregate(VEGF~Status, data = first_wk, quantile)
aggregate(VEGF~Status, data = second_wk, quantile)
aggregate(VEGF~Status, data = third_wk, quantile)
aggregate(VEGF~Status, data = fourth_wk, quantile)
#sFlt1
wilcox.test(sFlt1~Status, data = first_wk)
wilcox.test(sFlt1~Status, data = second_wk)
wilcox.test(sFlt1~Status, data = third_wk)
wilcox.test(sFlt1~Status, data = fourth_wk)
aggregate(sFlt1~Status, data = first_wk, quantile)
aggregate(sFlt1~Status, data = second_wk, quantile)
aggregate(sFlt1~Status, data = third_wk, quantile)
aggregate(sFlt1~Status, data = fourth_wk, quantile)
#ratio
wilcox.test(sFlt1_P1GF~Status, data = first_wk)
wilcox.test(sFlt1_P1GF~Status, data = second_wk)
wilcox.test(sFlt1_P1GF~Status, data = third_wk)
wilcox.test(sFlt1_P1GF~Status, data = fourth_wk)
aggregate(sFlt1_P1GF~Status, data = first_wk, quantile)
aggregate(sFlt1_P1GF~Status, data = second_wk, quantile)
aggregate(sFlt1_P1GF~Status, data = third_wk, quantile)
aggregate(sFlt1_P1GF~Status, data = fourth_wk, quantile)

######### lets us do normal logistic regression ####
#
library(pROC)
mydata <- sheila# read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
mylogit <- glm(case ~ VEGF, data = mydata, family = "binomial")
summary(mylogit)
prob=predict(mylogit,type=c("response"))
mydata$prob=prob
gvegf <-roc(case ~ prob, data = mydata)
plot(gvegf)
# library(caret)
# # Use your model to make predictions, in this example newdata = training set, but replace with your test set    
# pdata <- predict(logitMod, newdata = train, type = "response")
# 
# # use caret and compute a confusion matrix
# confusionMatrix(data = as.numeric(pdata>0.5), reference = train$LoanStatus_B)

gvegf_roc<-data.frame(ses=gvegf$sensitivities, spc=gvegf$specificities)
gvegf_roc$Analyte<-"VEGF"

mylogit <- glm(case ~ sFlt1, data = mydata, family = "binomial")
summary(mylogit)
prob=predict(mylogit,type=c("response"))
mydata$prob=prob
gsflt1 <- roc(case ~ prob, data = mydata)
gsflt1_roc<-data.frame(ses=gsflt1$sensitivities, spc=gsflt1$specificities)
gsflt1_roc$Analyte<-"sFlT1"

mylogit <- glm(case ~ PlGF, data = mydata, family = "binomial")
summary(mylogit)
prob=predict(mylogit,type=c("response"))
mydata$prob=prob
gplgf <- roc(case ~ prob, data = mydata)
gplgf_roc<-data.frame(ses=gplgf$sensitivities, spc=gplgf$specificities)
gplgf_roc$Analyte<-"PlGF"

mylogit <- glm(case ~ sFlt1_P1GF, data = mydata, family = "binomial")
summary(mylogit)
prob=predict(mylogit,type=c("response"))
mydata$prob=prob
gratio <- roc(case ~ prob, data = mydata)
gratio_roc<-data.frame(ses=gratio$sensitivities, spc=gratio$specificities)
gratio_roc$Analyte<-"sFlT1/PlGF"

roc_combn<-rbind(gratio_roc, gplgf_roc,
                 gvegf_roc, gsflt1_roc)
roc_combn$spc<-1-roc_combn$spc
p<-ggplot(data = roc_combn, aes(x=spc, y=ses, color=Analyte))
p<-p+geom_line()+xlab('1-Specificity')+ylab('Sensitivity')
p<-p+theme_classic()
p

# allsig<-glm(case~VEGF+sFlt1+PlGF+
#                  historypreeclampsia+firstpregnancy+historyhypertension+
#                  hivstatus+condomuse+wksameno, data=sheila, family = "binomial")
# summary(allsig)

###### Overall sensitivity ########
set.seed(101) # Set Seed so that same sample can be reproduced in future also
# Now Selecting 75% of data as sample from total 'n' rows of the data  
sample <- sample.int(n = nrow(sheila), size = floor(.75*nrow(sheila)), replace = F)

train <- sheila[sample, ]
test  <- sheila[-sample, ]

dim(train)
dim(test)

library(caret)
train$case<-factor(train$case)
train$case<-relevel(train$case, ref = '0')

test$case<-factor(test$case)
test$case<-relevel(test$case, ref = '0')

mylogit <- glm(case ~ sFlt1_P1GF, data = train, family = binomial(link="logit"))
pdata <- predict(mylogit, newdata = test, type = "response")
confusionMatrix(data = as.factor(as.numeric(pdata>0.5)), reference = test$case, positive = '1')

mylogit <- glm(case ~ sFlt1, data = train, family = binomial(link="logit"))
pdata <- predict(mylogit, newdata = test, type = "response")
confusionMatrix(data = as.factor(as.numeric(pdata>0.5)), reference = test$case, positive = '1')

mylogit <- glm(case ~ PlGF, data = train, family = binomial(link="logit"))
pdata <- predict(mylogit, newdata = test, type = "response")
confusionMatrix(data = as.factor(as.numeric(pdata>0.5)), reference = test$case, positive = '1')

mylogit <- glm(case ~ VEGF, data = train, family = binomial(link="logit"))
pdata <- predict(mylogit, newdata = test, type = "response")
confusionMatrix(data = as.factor(as.numeric(pdata>0.5)), reference = test$case, positive = '1')


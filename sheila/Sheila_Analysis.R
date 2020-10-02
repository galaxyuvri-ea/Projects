rm(list = ls())
dev.off()
sheila<-read.csv("PE_labdata.csv")
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
require(reshape2)
df<-sheila[c(2,4:8)]
names(df)<-c("studyID","sFlt1","PlGF","VEGF","Status","sFlt1/PlGF")
require(reshape2)
require(tidyr)
df1<-melt(df, id.vars=c("studyID","Status"))

df1<-df1 %>% 
  dplyr::mutate(variable = factor(variable,levels = c('VEGF','PlGF','sFlt1','sFlt1/PlGF'))) 

require(ggplot2)
#analyte='sFlt1'
#df1s<-df1[df1$variable==analyte,]
p <- ggplot(df1, aes(Status, value, colour = factor(Status))) + geom_boxplot()+geom_jitter()
p<-p+theme_classic(base_family = "Arial Black")
p<-p + facet_wrap(~variable, scales = 'free_y')+theme_classic()
p<-p+xlab('')+ylab('Log-transformed plasma concentration')
p<-p+theme(text = element_text(size=20, face='bold'),
           axis.text.x = element_text(size=20, face='bold'),
           legend.position = 'none'
           )
#p$labels$colour<-"Status"
p


############ conditional logistic regression ########
clinic<-read.csv('clinic.csv')
clinic$wksameno<-cut(clinic$GA, 
                     breaks = c(19, 28, 33, 37, 43), 
                     labels = c("20-28", "29-33", "34-37", "38-43"))

hypert<-read.csv('dxhypertension.csv')
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

#summarise ctrls and cases by gestational age group
table(merged$Status.x,merged$wksameno)
#trends 
trend<-merged[c('studyID','Status.x',"sFlt1", "PlGF", "VEGF", "sFlt1_P1GF", "wksameno")]
names(trend)<-c('studyID','Status.x',"sFlt1", "PlGF", "VEGF", "sFlt1/PlGF", "wksameno")
trend<-melt(trend, id.vars=c("studyID","Status.x","wksameno"))
#trend$value<-log(trend$value)
trend_mean<-aggregate(value~wksameno+Status.x+variable, data = trend, median)
trend_mean$value<-log(trend_mean$value)
trend_sd<-aggregate(value~wksameno+Status.x+variable, data = trend, sd)
trend_mean$sd<-trend_sd$value

trend_mean<-trend_mean %>% 
  dplyr::mutate(variable = factor(variable,levels = c('VEGF','PlGF','sFlt1','sFlt1/PlGF'))) 

p <- ggplot(trend_mean, aes(x=wksameno, y=value, group=Status.x, color = factor(Status.x))) + 
  geom_point()+geom_line()+theme_classic()
#p<-p+geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
 #                  position=position_dodge(0.05))
p<-p + facet_grid(~variable, scales = 'free')#+theme_classic()
p<-p+xlab('Gestational age (weeks)')+ylab('Log-transformed plasma concentration')
  #+theme(text = element_text(size=30))
p<-p+theme(text = element_text(size=20, face = 'bold'),
           axis.text.x = element_text(angle=90,hjust = 1, 
                                      size=20,vjust=0.5, face = 'bold'))

p$labels$colour<-"Status"
p


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

modData$sFlt1<-log2(modData$sFlt1)
modData$PlGF<-log2(modData$PlGF)
modData$VEGF<-log2(modData$VEGF)
modData$sFlt1_P1GF<-log2(modData$sFlt1_P1GF)

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
modall<-clogit(case~VEGF+sFlt1+PlGF + strata(order), data = modData)
summary(modall)

modData$pregnancy7months<-relevel(modData$pregnancy7months, ref='None')
allsig<-clogit(case~VEGF+sFlt1+PlGF+historyhypertension+
                 historypreeclampsia
                 +strata(order), data=modData)
summary(allsig)

############## make plots of plasma levels stratified bu GA and PE status ###############
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
names(df2)<-c("studyID","Status", "sFlt1", "PlGF","VEGF", "sFlt1/PlGF","wksameno")
df3<-melt(df2, id.vars=c("studyID","Status","wksameno"))

############ Draw Heatmap ###########
#change order of angiogenic factors A-VEGF, B-PlGF, C-SFlt1
require(dplyr)
df3<-df3 %>% 
  dplyr::mutate(variable = factor(variable,levels = rev(c('VEGF','PlGF','sFlt1','sFlt1/PlGF')))) 
  
ph<-ggplot(data = df3,
          aes(x=studyID,y=variable,fill=value))
ph<-ph+geom_raster()+facet_grid(~wksameno+Status, scales = "free")+
 #scale_fill_gradient(name = "Plasma levels",low = "#FFFFFF",high = "#012345")
 scale_fill_distiller(palette = 'RdYlBu')
ph<-ph+ylab('')+xlab('')
ph<-ph+theme_classic(base_family = "Arial Black")+theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                            text = element_text(size=20, face='bold'))
ph$labels$fill<-""
ph

############ split into data into gestational age groups ############
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

###### Overall sensitivity ########hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhnj
library(caret)
library(pROC)

sheila<-modData
sheila$case<-ifelse(sheila$case=='0','Control', 'Case')
sheila$case<-factor(sheila$case)
sheila$case<-relevel(sheila$case, ref = 'Control')

set.seed(3456) 
train.control <- trainControl(method = "cv", 
                              number = 10,
                              savePredictions = T,
                              #repeats = 5,
                              classProbs = T,
                              summaryFunction = twoClassSummary)
#Train the sFlt1/P1GF model
model_ratio <- train(case ~ sFlt1_P1GF, data=sheila, 
               trControl=train.control, method="glm",
               family=binomial(),
               metric='ROC')
#print(model_ratio)
pred.ratio<-predict(model_ratio)
confusionMatrix(pred.ratio, sheila$case, positive = 'Case')
#Draw roc curve
prob=predict(model_ratio$finalModel,type=c("response"))
sheila$ratio_prob=prob
ratio_roc <- roc(case ~ ratio_prob, data = sheila)
auc(ratio_roc)
ci.auc(ratio_roc)
ratio_roc_df<-data.frame(ses=ratio_roc$sensitivities, spc=ratio_roc$specificities)
ratio_roc_df$Analyte<-"sFlt1/PlGF"
coords(ratio_roc, "best", ret="threshold", best.method="youden", transpose = T)


#Train the sFlt1 model
model_sflt1 <- train(case ~ sFlt1, data=sheila, 
                     trControl=train.control, method="glm",
                     family=binomial(),
                     metric='ROC')
#print(model_ratio)
pred.sflt1<-predict(model_sflt1)
confusionMatrix(pred.sflt1, sheila$case, positive = 'Case')

#Draw roc curve
prob=predict(model_ratio$finalModel,type=c("response"))
sheila$sflt1_prob=prob
sflt1_roc <- roc(case ~ sflt1_prob, data = sheila)
auc(sflt1_roc)
ci.auc(sflt1_roc)
sflt1_roc_df<-data.frame(ses=sflt1_roc$sensitivities, spc=sflt1_roc$specificities)
sflt1_roc_df$Analyte<-"sFlt1"
coords(sflt1_roc, "best", ret="threshold", best.method="youden", transpose = T)

#Train the PlGF model
model_gplf <- train(case ~ PlGF, data=sheila, 
                     trControl=train.control, method="glm",
                     family=binomial(),
                     metric='ROC')
#print(model_ratio)
pred.gplf<-predict(model_gplf)
confusionMatrix(pred.gplf, sheila$case, positive = 'Case')

#Draw roc curve
prob=predict(model_gplf$finalModel,type=c("response"))
sheila$plgf_prob=prob
plgf_roc <- roc(case ~ plgf_prob, data = sheila)
auc(plgf_roc)
ci.auc(plgf_roc)
plgf_roc_df<-data.frame(ses=plgf_roc$sensitivities, spc=plgf_roc$specificities)
plgf_roc_df$Analyte<-"PlGF"
coords(plgf_roc, "best", ret="threshold", best.method="youden", transpose = T)

#Train the VEGF model
model_vegf <- train(case ~ VEGF, data=sheila, 
                    trControl=train.control, method="glm",
                    family=binomial(),
                    metric='ROC')
#print(model_ratio)
pred.vegf<-predict(model_vegf)
confusionMatrix(pred.vegf, sheila$case, positive = 'Case')

#Draw roc curve
prob=predict(model_vegf$finalModel,type=c("response"))
sheila$vegf_prob=prob
vegf_roc <- roc(case ~ vegf_prob, data = sheila)
auc(vegf_roc)
ci.auc(vegf_roc)
vegf_roc_df<-data.frame(ses=vegf_roc$sensitivities, spc=vegf_roc$specificities)
vegf_roc_df$Analyte<-"VEGF"
#ci.coords(vegf_roc, "best", ret="threshold", best.method="youden", 
 #         transpose = T, best.policy = 'random')
coords(vegf_roc, "best", ret="threshold", best.method="youden", 
          transpose = T)

################## Combine and plot ROC ###########
df4<-rbind(sflt1_roc_df, ratio_roc_df, plgf_roc_df, vegf_roc_df)
df4$spc<-1-df4$spc
df4<-df4 %>% 
  dplyr::mutate(Analyte = factor(Analyte,levels = c('VEGF','PlGF','sFlt1','sFlt1/PlGF'))) 
p<-ggplot(data = df4, aes(x=spc, y=ses, color=Analyte))
p<-p+geom_line(size=1.2)+xlab('1-Specificity')+ylab('Sensitivity')
p<-p+theme_classic()+theme(text = element_text(size=20,face='bold'))
p$labels$colour<-"Angiogenic factors"
p

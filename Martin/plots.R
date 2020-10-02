rm(list = ls())

#tanoti results summary
tanoti<-read.csv('tanoti_results.csv')
tanoti<-tanoti[tanoti$Virus!="",]
unique(tanoti$Virus)

p<-ggplot(data = tanoti,aes(x=PoolID,y=AverageDepth))
p<-p+geom_bar(stat="identity", fill='blue')
p<-p+facet_wrap(~Virus, scales = 'free')
p<-p+theme_classic()
p<-p+ylab("Average Depth")+xlab("Mosquito pools")#+scale_fill_gradient('Abundance', limits=c(0, 4), breaks = c(0, 1, 2, 3, 4),  low = "lightblue", high = "darkblue",na.value = 'black')
p<-p+ theme(text = element_text(size=25),axis.text.x = element_text(angle = 90, hjust = 1))
p


for (virus in unique(tanoti$Virus)) {
  tmp<-tanoti[tanoti$Virus==virus,]
}

###################============================#################################
viruses<-read.csv("~/Desktop/marin-prossy/viruses_research.csv")
viruses$year<-factor(viruses$year, levels = viruses$year[order(viruses$year)])

require(ggplot2)
p<-ggplot(data = viruses,aes(x=year,y=count))+ylim(0, 14)
p<-p+geom_bar(stat="identity", fill='blue')
p<-p+theme_classic()
p<-p+theme(strip.text = element_text(size = 10, colour = "black", angle = 90,vjust=0)) #vjust=0 aligns the labels to the bottom in this case
p<-p+ylab("Number of Viruses discovered")+xlab("Time (Years)")#+scale_fill_gradient('Abundance', limits=c(0, 4), breaks = c(0, 1, 2, 3, 4),  low = "lightblue", high = "darkblue",na.value = 'black')
#p<-p+theme_classic()
p<-p+ theme(text = element_text(size=25),axis.text.x = element_text(angle = 90, hjust = 1))
p

#clinical manifestation
manifest <- read.csv("~/Desktop/marin-prossy/clinic_manifestation.csv", row.names = 1)
df<-NULL
rn<-rownames(manifest)
cn<-colnames(manifest)
for (i in 1:dim(manifest)[1]) {
  for (j in 1:dim(manifest)[2]) {
    virus<-rn[i]
    Symptoms<-cn[j]
    value<-manifest[i,j]
    tmp<-data.frame(virus, Symptoms, value)
    df<-rbind(df, tmp)
  }
}

df<-df[!df$value%in%c('x','X'),]
df<-df[c(1,2)]

df<-data.frame(table(df$Symptoms))
df$Var1<-gsub('.', ' ', df$Var1, fixed = T)
df$Var1<-gsub('CNS signs encephalitis', 'CNS signs-encephalitis', df$Var1, fixed = T)
p<-ggplot(data = df,aes(x=Freq,y=reorder(Var1, Freq)))
p<-p+geom_bar(stat="identity", fill='blue')
p<-p+theme_classic()
p<-p+ylab("")+xlab("Frequency of symptom in different viruses")#+scale_fill_gradient('Abundance', limits=c(0, 4), breaks = c(0, 1, 2, 3, 4),  low = "lightblue", high = "darkblue",na.value = 'black')
#p<-p+theme_classic()
p<-p+ theme(text = element_text(size=25))+geom_text(aes(label = Freq), nudge_x = 0.5)
p

######## Animal Host ###########
virus_host <- read.csv("~/Desktop/marin-prossy/virus_animal_host.csv", row.names = 1)
virus_host$Family<-NULL
df<-NULL
rn<-rownames(virus_host)
cn<-colnames(virus_host)
for (i in 1:dim(virus_host)[1]) {
  for (j in 1:dim(virus_host)[2]) {
    virus<-rn[i]
    Host<-cn[j]
    value<-virus_host[i,j]
    tmp<-data.frame(virus, Host, value)
    df<-rbind(df, tmp)
  }
}

df$value[df$value=='x']<-'Not isolated'
df$value[df$value=='v']<-'Isolated'
df$value[df$value=='n']<-'Unknown'

#df$value<-as.numeric(df$value)

ph<-ggplot(data = df,aes(x=Host,y=virus,fill=value))
ph<-ph+geom_tile()+scale_colour_manual(values = c("yellow", "blue", "red"))
ph<-ph+ylab('')+xlab('')
ph<-ph+theme_classic()+theme(text = element_text(size=25, face='bold'),
                             axis.text.x = element_text(angle = 90, 
                                                        size=25, face='bold',hjust = 1))
ph$labels$fill<-""
ph

#df<-df[c(1,2)]

# df <- df %>%
#   group_by(virus, Symptoms) %>%
#   summarise(count = n()) %>%
#   mutate(cut.count = sum(count),
#          prop = count/sum(count)) %>%
#   ungroup()
# 
# p<-ggplot(df,
#        aes(x = virus, y = prop, width = cut.count, fill = Symptoms)) +
#   geom_bar(stat = "identity",colour = "black") +theme_classic()+
# 
#   # geom_text(aes(label = scales::percent(prop)), position = position_stack(vjust = 0.5)) + # if labels are desired
#   facet_grid(~virus, scales = "free_x") #+
#   #scale_fill_brewer(palette = "RdYlGn") +
#   #theme(panel.spacing.x = unit(0, "npc")) #+ # if no spacing preferred between bars
#   #theme_void() 
# p<-p+theme(strip.text = element_text(size = 25, colour = "black", angle = 90,vjust=0)) #vjust=0 aligns the labels to the bottom in this case
# p<-p+ylab("")+xlab("Viruses")#+scale_fill_gradient('Abundance', limits=c(0, 4), breaks = c(0, 1, 2, 3, 4),  low = "lightblue", high = "darkblue",na.value = 'black')
# p<-p+ theme(text = element_text(size=25),axis.text.x = element_blank(), 
#             axis.text.y = element_blank())
# p

#### Generates figures and relevant statistics

#Load these libraries 
library(ggplot2)
library(ez)
library(nlme)
library(gridExtra)
library(plyr)
library(lme4)
library(lmerTest)
library(psych)

#### Set working directory
workingdir='/Users/Matilde/Dropbox/Vortex_bet/GitHub/Data/RPlotData/'
setwd (workingdir)

outputdir='/Users/Matilde/Dropbox/Vortex_bet/GitHub/Plots/'

##############################
# ZSCORES ALL TRIALS 
##############################
data<-read.table('VortexBet_LR_CONF_ZSCORES.txt', header=FALSE,sep=",")
colnames(data)<-c("zValues", "TrialType", "Sub","Group")

data$Group[data$Group==1]<- "CTL"
data$Group[data$Group==0]<- "OCD"
data$TrialType[data$TrialType==1]<- "zLR"
data$TrialType[data$TrialType==2]<- "zCONF"


options(scipen=999)
ezANOVA(data,
        dv = .(zValues),
        wid = .(Sub),  # subject
        within = .(TrialType),  # confidence, lr 
        between=.(Group),
        type=2, 
        # detailed=TRUE,
        return_aov=TRUE)

library(cowplot)
options( scipen = 10 )

ggplot(data, aes(x=TrialType, y=zValues, fill=factor(Group))) + 
  stat_summary(fun.y=mean, geom="bar",position=position_dodge(1)) + 
  scale_color_discrete("Group")+
  scale_colour_manual(values=c("grey28", "darkred"))+
  scale_fill_manual(values=c("gainsboro", "red3"))+
  labs(y = "zValues")+
  theme(axis.text=element_text(size=12 ),
        axis.title=element_text(size=12), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        axis.title.x=element_blank())

#ZSCORES DIVIDING IN BINS OF 30 TRIALS [10 BINS, 1-30, 31-60 ETC...]

datatime<-read.table('VortexBet_LR_CONF_ZSCORE_TIME.txt', header=FALSE,sep=",")
colnames(datatime)<-c("zValues","Time", "TrialType", "Sub","Group")

datatime$Group[datatime$Group==1]<- "CTL"
datatime$Group[datatime$Group==0]<- "OCD"
datatime$TrialType[datatime$TrialType==1]<- "zLR"
datatime$TrialType[datatime$TrialType==2]<- "zCONF"

datatime <- within(datatime, {
  timeFactor <- as.factor(Time) #Force ProgCont to be discrete
})

options(scipen=999)
ezANOVA(datatime,
        dv = .(zValues),
        wid = .(Sub),  # subject
        within = .(TrialType, timeFactor),  # confidence, lr 
        between=.(Group),
        type=2, 
        # detailed=TRUE,
        return_aov=TRUE)

ggplot(datatime, aes(x=(Time), y=zValues, fill=factor(Group))) + 
  facet_grid(.~TrialType)+
  stat_summary(fun.y=mean, geom="bar", position=position_dodge(1)) + 
  scale_color_discrete("Group")+
  scale_colour_manual(values=c("grey28", "darkred"))+
  scale_fill_manual(values=c("gainsboro", "red3"))+
  labs(y = "zValues")+
  theme(axis.text=element_text(size=12 ),
        axis.title=element_text(size=12), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        axis.title.x=element_blank())+
  scale_x_continuous(breaks=seq(0, 10.5, 1), limits = c(0, 10.5))

##############################
#DATA FROM CHANGE POINT 
###############################


################################
#Figure as original paper but with zscores also for learning rate to compare 
################################
##############################
# Change points  - Figure 2
##############################
changePoint<-read.table('VortexBet_ZVALUES_ChangePoint.txt', header=TRUE,sep="")
changePoint$Cpoint<- c("Cpoint1", "Cpoint2", "Cpoint3", "Cpoint4", "Cpoint5", "Cpoint6", "Cpoint7", "Cpoint8", "Cpoint9")


library(cowplot)
humanLR <- ggplot(data = changePoint) + 
  geom_errorbar(aes(x=Cpoint, ymin=CTL_LRydat-CTL_LRedat, ymax=CTL_LRydat+CTL_LRedat), width=.15,lwd=0.5,colour="black")+
  geom_point(aes(x = Cpoint, y = CTL_LRydat),pch=21, fill="gainsboro",colour='grey28', size=3.5) +
  
  geom_errorbar(aes(x=Cpoint, ymin=OCD_LRydat-OCD_LRedat, ymax=OCD_LRydat+OCD_LRedat), width=.15,lwd=0.5,colour="red3")+
  geom_point(aes(x = Cpoint, y = OCD_LRydat),pch=21, fill="red3",colour='darkred',size=3.5) +
  
  geom_vline(xintercept = 5.5, colour="black", size=0.5, linetype="dashed")+
  scale_x_discrete(breaks = c("Cpoint1","Cpoint5", "Cpoint9"),labels = c("-4","0","4"))+
  labs(x = "Trials after change point")+ 
  labs(y = "Human learning rate, a.u.")+
  theme(axis.text=element_text(size=12 ))

humanCF<- ggplot(data = changePoint) + 
  geom_errorbar(aes(x=Cpoint, ymin=CTL_CFydat-CTL_CFedat, ymax=CTL_CFydat+CTL_CFedat), width=.15,lwd=0.5,colour="black")+
  geom_point(aes(x = Cpoint, y = CTL_CFydat),pch=21, fill="gainsboro",colour='grey28', size=3.5) +
  
  geom_errorbar(aes(x=Cpoint, ymin=OCD_CFydat-OCD_CFedat, ymax=OCD_CFydat+OCD_CFedat), width=.15,lwd=0.5,colour="red3")+
  geom_point(aes(x = Cpoint, y = OCD_CFydat),pch=21, fill="red3",colour='darkred',size=3.5) +
  
  geom_vline(xintercept = 5.5, colour="black", size=0.5, linetype="dashed")+
  scale_x_discrete(breaks = c("Cpoint1","Cpoint5", "Cpoint9"),labels = c("-4","0","4"))+
  labs(x = "Trials after change point")+ 
  labs(y = "Human confidence, a.u.")+
  theme(axis.text=element_text(size=12 ))



################################
#Stats analyses 
################################

dataCPP<-read.table('VortexBet_ZVALUES_LRCF_CPP.txt', header=FALSE,sep=",")
colnames(dataCPP)<-c("C1","C2", "C3", "C4", "C5", "C6", "C7", "C8","C9", "Sub", "TrialType")

dataCPP$TrialType[dataCPP$TrialType==1]<- "zLR"
dataCPP$TrialType[dataCPP$TrialType==2]<- "zCONF"

dataCPP$Group[substr(dataCPP$Sub, start=1, stop=1)=="7"]<-"CTL"
dataCPP$Group[substr(dataCPP$Sub, start=1, stop=1)=="8"]<-"OCD"

library(reshape)
dataCPP.m = melt(dataCPP, id=c("TrialType", "Sub", "Group"))

names(dataCPP.m )[names(dataCPP.m ) == 'variable'] <- 'SamplesCPP'
names(dataCPP.m )[names(dataCPP.m ) == 'value'] <- 'zValues'

#ANOVA considering all data points before and after 
options(scipen=999)
ezANOVA(dataCPP.m,
        dv = .(zValues),
        wid = .(Sub),  # subject
        within = .(TrialType, SamplesCPP),  # confidence, lr 
        between=.(Group),
        type=2, 
        # detailed=TRUE,
        return_aov=TRUE)

ggplot(dataCPP.m, aes(x=(SamplesCPP), y=zValues, fill=factor(Group))) + 
  facet_grid(.~TrialType)+
  stat_summary(fun.y=mean, geom="bar", position=position_dodge(1)) + 
  scale_color_discrete("Group")+
  scale_colour_manual(values=c("grey28", "darkred"))+
  scale_fill_manual(values=c("gainsboro", "red3"))+
  labs(y = "zValues")+
  theme(axis.text=element_text(size=12 ),
        axis.title=element_text(size=12), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        axis.title.x=element_blank())+
  scale_x_continuous(breaks=seq(0, 10.5, 1), limits = c(0, 10.5))

#ANOVA considering only BEFORE CHANGE POINT
dataCPPbefore.m<-subset(dataCPP.m, (SamplesCPP=="C1"|SamplesCPP=="C2"|SamplesCPP=="C3" | SamplesCPP=="C4"|SamplesCPP=="C5"))
options(scipen=999)
ezANOVA(dataCPPbefore.m,
        dv = .(zValues),
        wid = .(Sub),  # subject
        within = .(TrialType, SamplesCPP),  # confidence, lr 
        between=.(Group),
        type=2, 
        # detailed=TRUE,
        return_aov=TRUE)

dataCPPafter.m<-subset(dataCPP.m, (SamplesCPP=="C7"|SamplesCPP=="C8"|SamplesCPP=="C9"))
options(scipen=999)
ezANOVA(dataCPPafter.m,
        dv = .(zValues),
        wid = .(Sub),  # subject
        within = .(TrialType, SamplesCPP),  # confidence, lr 
        between=.(Group),
        type=2, 
        # detailed=TRUE,
        return_aov=TRUE)

dataCPPafter.m<-subset(dataCPP.m, (SamplesCPP=="C8"|SamplesCPP=="C9"))
options(scipen=999)
ezANOVA(dataCPPafter.m,
        dv = .(zValues),
        wid = .(Sub),  # subject
        within = .(TrialType, SamplesCPP),  # confidence, lr 
        between=.(Group),
        type=2, 
        # detailed=TRUE,
        return_aov=TRUE)








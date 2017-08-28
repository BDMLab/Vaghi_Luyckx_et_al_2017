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
# Learning Rate - Figure 1 
##############################
data<-read.table('VortexBet_LR.txt', header=TRUE,sep=",")
colnames(data)<-c("LR", "Sub", "Group")

data$Group[data$Group==1]<- "CTL"
data$Group[data$Group==0]<- "OCD"


library(car)
leveneTest(data$LR, data$Group)
t.test(data$LR~data$Group, var.equal = FALSE)


dsdata<-ddply(data, .(Group), summarise, 
              LRN = sum(!is.na(LR)),
              LRmean =mean(LR, na.rm=TRUE), 
              LRsd=sd(LR, na.rm=TRUE), 
              LRse =LRsd/sqrt(LRN))

library(cowplot)
Plot_LearningRate<- ggplot(dsdata,aes(x=Group, y=LRmean, fill=Group))+ 
  geom_point(data=data,aes(x=Group, y=LR, fill=Group, colour=Group),
             position=position_jitterdodge(jitter.width = 0.6, jitter.height=0, dodge.width=0.5), 
             shape=21, size=2.5, alpha=0.8)+
  scale_colour_manual(values=c("grey28", "darkred"))+
  scale_fill_manual(values=c("gainsboro", "red3"))+
  geom_errorbar(aes(ymin=LRmean-LRse, ymax=LRmean+LRse, color=Group), width=.1,lwd=0.5,
                position=position_dodge(width = 0.5), colour="blue1")+
  stat_summary(data=data,aes(x=Group, y=LR, fill=Group, colour=Group),
               fun.y="mean", size=2.5, geom="point", 
               position=position_dodge(width = 0.5), colour='blue1')+
  labs(y = expression(paste ("Estimated learning rate, " , hat(alpha) )))+
  theme(axis.text=element_text(size=12 ),
        axis.title=element_text(size=12), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        axis.title.x=element_blank(),
        legend.position='none')+
  scale_y_continuous(breaks=seq(0, 1.05, 0.25), limits = c(0, 1.05))

figname= paste(outputdir, 'Plot_LearningRate.tiff')
tiff(file = figname, width = 2, height = 3, units = "in", res = 300)
Plot_LearningRate
dev.off()



##############################
# Error Magnitude - Figure 1 
##############################

ErrorMagn<-read.table('VortexBet_ErrorMag.txt', header=FALSE,sep=",")
colnames(ErrorMagn)<-c("Bin", "XVALUES", "MEAN", "SEM", "Group","CAT")
ErrorMagn$Group[ErrorMagn$Group==1]<- "CTL"
ErrorMagn$Group[ErrorMagn$Group==0]<- "OCD"


library(cowplot)
Plot_ErrMagnitude<-ggplot(ErrorMagn,aes(x=XVALUES, y=MEAN,fill=Group))+ 
  annotate("rect", xmin=0, xmax=9, ymin=0, ymax=Inf, alpha=0.1, fill="gray66") +
  annotate("rect", xmin=9.1, xmax=20, ymin=0, ymax=Inf, alpha=0.2, fill="gray66") +
  annotate("rect", xmin=20.1, xmax=Inf, ymin=0, ymax=Inf, alpha=0.3, fill="gray66") +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=MEAN-SEM, ymax=MEAN+SEM, color=Group), width=1.25,lwd=0.5)+
  scale_colour_manual(values=c("grey28", "darkred"))+
  geom_point(aes(fill=Group, colour=Group),  pch= 21, size= 2.5)+
  # geom_line(aes(group = Group),  lty = 3, lwd = 1.3, color='black')+
  scale_fill_manual(name="", values=c("gainsboro" , "red3"))+
  labs(y = expression(paste ("Estimated learning rate, " , hat(alpha) )))+
  labs(x = "Error magnitude (i.e., spatial prediction error)")+ 
  theme(axis.text=element_text(size=12 ),
        axis.title=element_text(size=12), 
        legend.position='none',
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))+
  #scale_x_continuous(breaks=seq(0, 60, 10), limits = c(0, 60))+
  # scale_y_continuous(breaks=seq(0, 1.05, 0.25), limits = c(0, 1.05))+
  geom_segment(aes(x=1.5, y=1, xend=51.5, yend=1), size=0.05) +
  geom_segment(aes(x=1.5, y=0.65, xend=1.5, yend=1), size=0.05) +
  geom_segment(aes(x=51.5, y=0.97, xend=51.5, yend=1), size=0.05) +
  geom_segment(aes(x=54, y=0.97, xend=56, yend=0.97), size=0.05) +
  geom_segment(aes(x=56, y=0.97, xend=56, yend=0.05), size=0.05)+
  geom_segment(aes(x=54, y=0.05, xend=56, yend=0.05), size=0.05)+
  annotate("text", x = 4, y = 1.04, label = "low", size=3)+
  annotate("text", x = 14, y = 1.04, label = "medium", size=3)+
  annotate("text", x = 30, y = 1.04, label = "high", size=3)+
  annotate("text", x = 20, y = 1.01, label = "***", size=4)+
  annotate("text", x = 57.5, y = 0.5, label = "***", size=4,  angle = 90)

figname= paste(outputdir, 'Plot_ErrorMagnitude.tiff')
tiff(file = figname, width = 4.5, height = 3, units = "in", res = 300)
Plot_ErrMagnitude
dev.off()


## ANOVA based on 3 bins #####

ErrorMag3<-read.table('ErrorMagn_3bins.txt', header=FALSE,sep=",")
colnames(ErrorMag3)<-c("Low", "Medium", "High", "Sub", "Group")
ErrorMag3$Group[ErrorMag3$Group==1]<- "CTL"
ErrorMag3$Group[ErrorMag3$Group==0]<- "OCD"


leveneTest(ErrorMag3$Low, ErrorMag3$Group)
t.test(ErrorMag3$Low~ErrorMag3$Group, var.equal = FALSE)

leveneTest(ErrorMag3$Medium, ErrorMag3$Group)
t.test(ErrorMag3$Medium~ErrorMag3$Group, var.equal = FALSE)

leveneTest(ErrorMag3$High, ErrorMag3$Group)
t.test(ErrorMag3$High~ErrorMag3$Group, var.equal = TRUE)


library(reshape)
data.m = melt(ErrorMag3, id=c('Group', 'Sub'))

options(scipen=999)
ezANOVA(data.m,
        dv = .(value),
        wid = .(Sub),  # subject
        within = .(variable),  # high, medium, low
        between=.(Group),
        type=2, 
        # detailed=TRUE,
        return_aov=TRUE)

##############################
### STATS FOR 20 BINS 
##############################
ErrorMag20<-read.table('ErrorMagn_20bins.txt', header=FALSE,sep=",")
colnames(ErrorMag20)<-c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10", "Q11", "Q12", 
                       "Q13", "Q14", "Q15", "Q16", "Q17", "Q18", "Q19", "Q20", "Sub", "Group")
ErrorMag20$Group[ErrorMag20$Group==1]<- "CTL"
ErrorMag20$Group[ErrorMag20$Group==0]<- "OCD"

library(reshape)
data.m = melt(ErrorMag20, id=c('Group', 'Sub'))

options(scipen=999)
ezANOVA(data.m,
        dv = .(value),
        wid = .(Sub),  # subject
        within = .(variable),  # high, medium, low
        between=.(Group),
        type=2, 
        # detailed=TRUE,
        return_aov=TRUE)


##############################
### STATS FOR 3 bins but excluding values of first Q1 BASED ON20 BINS ?? 
##############################

ErrorMag3<-read.table('ErrorMagn_3bins_add.txt', header=FALSE,sep=",")
colnames(ErrorMag3)<-c("Low", "Medium", "High", "Sub", "Group")
ErrorMag3$Group[ErrorMag3$Group==1]<- "CTL"
ErrorMag3$Group[ErrorMag3$Group==0]<- "OCD"

leveneTest(ErrorMag3$Low, ErrorMag3$Group)
t.test(ErrorMag3$Low~ErrorMag3$Group, var.equal = FALSE)

library(reshape)
data.m = melt(ErrorMag3, id=c('Group', 'Sub'))

options(scipen=999)
ezANOVA(data.m,
        dv = .(value),
        wid = .(Sub),  # subject
        within = .(variable),  # high, medium, low
        between=.(Group),
        type=2, 
        # detailed=TRUE,
        return_aov=TRUE)






##############################
# Change points  - Figure 2
##############################
changePoint<-read.table('VortexBet_ChangePoint.txt', header=TRUE,sep="")
changePoint$Cpoint<- c("Cpoint1", "Cpoint2", "Cpoint3", "Cpoint4", "Cpoint5", "Cpoint6", "Cpoint7", "Cpoint8", "Cpoint9")

changePointModel<-read.table('VortexBet_ChangePointModel.txt', header=TRUE,sep="")
changePointModel$Cpoint<- c("Cpoint1", "Cpoint2", "Cpoint3", "Cpoint4", "Cpoint5", "Cpoint6", "Cpoint7", "Cpoint8", "Cpoint9")


library(cowplot)
modelLR<- ggplot(data=changePoint, aes( x=Cpoint, y=Model_LRydat, group=1))+ 
  geom_line(colour="black")+
  geom_point(size=3.5, colour="black")+
  geom_errorbar(aes(ymin=Model_LRydat-Model_LRedat, ymax=Model_LRydat+Model_LRedat, color=Group), width=.15,lwd=0.5,colour="black")+
  geom_vline(xintercept = 5.5, colour="black", size=0.5, linetype="dashed")+
  scale_x_discrete(breaks = c("Cpoint1","Cpoint5", "Cpoint9"),labels = c("-4","0","4"))+
  labs(x = "Trials after change point")+ 
  labs(y = "Model learning rate, a.u.")+
  theme(axis.text=element_text(size=12 ))+
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits = c(0, 1))

modelCF<- ggplot(data=changePoint, aes( x=Cpoint, y=Model_CFydat, group=1))+ 
  geom_line(colour="black")+
  geom_point(size=3.5, colour="black")+
  geom_errorbar(aes(ymin=Model_CFydat-Model_CFedat, ymax=Model_CFydat+Model_CFedat, color=Group), width=.15,lwd=0.5,colour="black")+
  geom_vline(xintercept = 5.5, colour="black", size=0.5, linetype="dashed")+
  scale_x_discrete(breaks = c("Cpoint1","Cpoint5", "Cpoint9"),labels = c("-4","0","4"))+
  labs(x = "Trials after change point")+ 
  labs(y = "Model confidence, a.u.")+
  theme(axis.text=element_text(size=12 ))+
  scale_y_continuous(breaks=seq(0.4, 0.9, 0.1), limits = c(0.4, 0.9))

humanLR <- ggplot(data = changePoint) + 
  
  geom_errorbar(aes(x=Cpoint, ymin=CTL_LRydat-CTL_LRedat, ymax=CTL_LRydat+CTL_LRedat), width=.15,lwd=0.5,colour="black")+
  geom_point(aes(x = Cpoint, y = CTL_LRydat),pch=21, fill="gainsboro",colour='grey28', size=3.5) +
  
  geom_errorbar(aes(x=Cpoint, ymin=OCD_LRydat-OCD_LRedat, ymax=OCD_LRydat+OCD_LRedat), width=.15,lwd=0.5,colour="red3")+
  geom_point(aes(x = Cpoint, y = OCD_LRydat),pch=21, fill="red3",colour='darkred',size=3.5) +
  
  geom_vline(xintercept = 5.5, colour="black", size=0.5, linetype="dashed")+
  scale_x_discrete(breaks = c("Cpoint1","Cpoint5", "Cpoint9"),labels = c("-4","0","4"))+
  labs(x = "Trials after change point")+ 
  labs(y = "Human learning rate, a.u.")+
  theme(axis.text=element_text(size=12 ))+
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits = c(0, 1))
  

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

  figname= paste(outputdir, 'Plot_ChangePoint.tiff')
  tiff(file = figname, width = 7, height = 5.5, units = "in", res = 300)
  grid.arrange(modelLR, humanLR, modelCF, humanCF, nrow=2)
  dev.off()

#Figure with only first change point for model 
changePointModel<-read.table('VortexBet_ChangePointModel.txt', header=TRUE,sep="")
changePointModel$Cpoint<- c("Cpoint1", "Cpoint2", "Cpoint3", "Cpoint4", "Cpoint5", "Cpoint6", "Cpoint7", "Cpoint8", "Cpoint9")
  
  
  library(cowplot)
  modelLR<- ggplot(data=changePointModel, aes( x=Cpoint, y=Model_LRydatFirstCP, group=1))+ 
    geom_line(colour="black")+
    geom_point(size=3.5, colour="black")+
    geom_vline(xintercept = 5.5, colour="black", size=0.5, linetype="dashed")+
    scale_x_discrete(breaks = c("Cpoint1","Cpoint5", "Cpoint9"),labels = c("-4","0","4"))+
    labs(x = "Trials after change point")+ 
    labs(y = "Model learning rate, a.u.")+
    theme(axis.text=element_text(size=12 ))+
    scale_y_continuous(breaks=seq(0, 1, 0.1), limits = c(0, 1))
  
  modelCF<- ggplot(data=changePointModel, aes( x=Cpoint, y=Model_CFydatFirstCP, group=1))+ 
    geom_line(colour="black")+
    geom_point(size=3.5, colour="black")+
    geom_vline(xintercept = 5.5, colour="black", size=0.5, linetype="dashed")+
    scale_x_discrete(breaks = c("Cpoint1","Cpoint5", "Cpoint9"),labels = c("-4","0","4"))+
    labs(x = "Trials after change point")+ 
    labs(y = "Model confidence, a.u.")+
    theme(axis.text=element_text(size=12 ))+
    scale_y_continuous(breaks=seq(0.4, 0.9, 0.1), limits = c(0.4, 0.9))
  
  
  figname= paste(outputdir, 'ModelLR_SINGLECP.tiff')
  tiff(file = figname, width =4, height = 3, units = "in", res = 300)
  modelLR
  dev.off()
  
  figname= paste(outputdir, 'ModelCF_SINGLECP.tiff')
  tiff(file = figname, width =4, height = 3, units = "in", res = 300)
  modelCF
  dev.off()
  

  
##############################
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

##############################
# Regressions  - Figure 2
##############################

##############################
library(cowplot)
#### Confidence 
RegConf<-read.table('VortexBet_RegConfidence.txt', header=TRUE,sep="")
colnames(RegConf)<-c("Int", "|PE|", "CPP", "(1-CPP)*(1-MC)", "Hit","Group")

RegConf$Group[RegConf$Group==1]<- "CTL"
RegConf$Group[RegConf$Group==0]<- "OCD"
RegConf$Int<-NULL

library(reshape)
conf.m<-melt(RegConf, id=c('Group'))

desconf<-summarySE(conf.m, measurevar="value",  groupvars=c( "variable", "Group"), na.rm=FALSE, conf.interval=.95)

regconf<- ggplot(data=desconf, aes(x=variable, y=value, fill=Group))+  
          geom_hline(yintercept = 0.0, colour="black", size=0.5, linetype="dashed")+
          geom_errorbar(aes(ymin=value-se, ymax=value+se, colour=Group), width=.15,lwd=0.5, position=position_dodge(width = 0.5))+
          geom_point(data=desconf,aes(x=variable, y=value, fill=Group, colour=Group), size=3.5, position=position_dodge(width = 0.5), shape=21)+
          scale_fill_manual(values=c("gainsboro", "red3"), name="")+
          scale_colour_manual(values=c("grey28", "darkred"),  name="")+
          scale_y_continuous(breaks=seq(-0.5, 0.5, 0.1), limits = c(-0.5, 0.5))+
          labs(x = "Predictors")+ 
          labs(y = "Regression coefficients \n beta weights")+
  theme(axis.text=element_text(size=12 ),
        axis.title=element_text(size=12, face="bold"), 
        legend.text=element_text(size=12),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        legend.position=c(.90,.95) )


figname= paste(outputdir, 'Plot_RegConfidence.tiff')
tiff(file = figname, width = 4.5, height = 3, units = "in", res = 300)
regconf
dev.off()

#### Learning Rate 
RegLearn<-read.table('VortexBet_RegLearningRate.txt', header=TRUE,sep="")
colnames(RegLearn)<-c("Int", "|PE|", "CPP", "(1-CPP)*(1-MC)", "Hit","Group")

RegLearn$Group[RegLearn$Group==1]<- "CTL"
RegLearn$Group[RegLearn$Group==0]<- "OCD"
RegLearn$Int<-NULL


lr.m<-melt(RegLearn, id=c('Group'))

deslr<-summarySE( lr.m, measurevar="value",  groupvars=c( "variable", "Group"), na.rm=FALSE, conf.interval=.95)

reglr<- ggplot(data=deslr, aes(x=variable, y=value, fill=Group))+  
        geom_hline(yintercept = 0.0, colour="black", size=0.5, linetype="dashed")+
        geom_errorbar(aes(ymin=value-se, ymax=value+se, colour=Group), width=.15,lwd=0.5, position=position_dodge(width = 0.5))+
        geom_point(data=deslr,aes(x=variable, y=value, fill=Group, colour=Group), size=3.5, position=position_dodge(width = 0.5), shape=21)+
        scale_fill_manual(values=c("gainsboro", "red3"), name="")+
       scale_colour_manual(values=c("grey28", "darkred"),  name="")+
       scale_y_continuous(breaks=seq(-1, 1.5, 0.5), limits = c(-1, 1.5))+
       labs(x = "Predictors")+ 
       labs(y = "Regression coefficients \n beta weights")+
       theme(axis.text=element_text(size=12 ),
        axis.title=element_text(size=12, face="bold"), 
        legend.text=element_text(size=12),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        legend.position=c(.90,.95) )+

  geom_segment(aes(x=0.9, y=0.9, xend=1.1, yend=0.9), size=0.4) +
  geom_segment(aes(x=0.9, y=0.85, xend=0.9, yend=0.9), size=0.4)+ 
  geom_segment(aes(x=1.1, y=0.85, xend=1.1, yend=0.9), size=0.4)+ 
  geom_segment(aes(x=1.9, y=0.9, xend=2.1, yend=0.9), size=0.4) +
  geom_segment(aes(x=1.9, y=0.85, xend=1.9, yend=0.9), size=0.4)+
  geom_segment(aes(x=2.1, y=0.85, xend=2.1, yend=0.9), size=0.4)+
  annotate("text", x = 1, y = 1, label = "*", size=7)+
  annotate("text", x = 2, y = 1, label = "**", size=7)

figname= paste(outputdir, 'Plot_RegLearning.tiff')
tiff(file = figname, width = 4.5, height = 3, units = "in", res = 300)
reglr
dev.off()

##############################
# Regressions  - Figure 3
##############################

### Discrepancy 
RegDisc<-read.table('VortexBet_RegDiscrepancy.txt', header=TRUE,sep="")
colnames(RegDisc)<-c("Int", "DIFFCONF", "Group", "OCI", "GroupCheck")

RegDisc$Group[RegDisc$Group==1]<- "CTL"
RegDisc$Group[RegDisc$Group==0]<- "OCD"
RegDisc$Int<-NULL
RegDisc$GroupCheck<-NULL


disc.m<-melt(RegDisc, id=c('Group', 'OCI'))

desdisc<-summarySE( disc.m, measurevar="value",  groupvars=c( "variable", "Group"), na.rm=FALSE, conf.interval=.95)

regdisc<- ggplot(data=desdisc, aes(x=Group, y=value, fill=Group))+  
          geom_hline(yintercept = 0.0, colour="black", size=0.5, linetype="dashed")+
          geom_point(data=RegDisc, aes(x=Group, y= DIFFCONF, fill=Group, colour=Group), 
             position=position_jitterdodge(jitter.width=0.6, jitter.height=0, dodge.width = 0.5), 
             shape=21, size=2, alpha=0.8)+
          scale_colour_manual(values=c("grey28", "darkred"))+
          scale_fill_manual(values=c("gainsboro", "red3"))+
          geom_errorbar(aes(ymin=value-se, ymax=value+se, color=Group), width=.15,lwd=0.5,
                position=position_dodge(width = 0.5), colour='blue1')+
          stat_summary(data=RegDisc, aes(x=Group, y=DIFFCONF, fill=Group, colour=Group), 
               fun.y="mean", size=2, geom='point', position=position_dodge(width=0.5), colour="blue1")+
          scale_y_continuous(breaks=seq(-0.4, 0.52, 0.1), limits = c(-0.4, 0.52))+
     labs(y = "Regression coefficients \n beta weights")+
     theme(axis.text=element_text(size=12 ),
        axis.title=element_text(size=12, face="bold"), 
        legend.text=element_text(size=12),
        axis.title.x=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"), 
        legend.position="none")+
  geom_segment(aes(x=1, y=0.5, xend=2, yend=0.5), size=0.4) +
  geom_segment(aes(x=1, y=0.47, xend=1, yend=0.5), size=0.4)+ 
  geom_segment(aes(x=2, y=0.47, xend=2, yend=0.5), size=0.4)+
  annotate("text", x = 1.5, y = 0.51, label = "*", size=7)

figname= paste(outputdir, 'Plot_RegDiscrepancy.tiff')
tiff(file = figname, width = 2, height = 3.5, units = "in", res = 300)
regdisc
dev.off()


##############################
# Correlation  - Figure 3
##############################

valuesOCD <- RegDisc [which (RegDisc$Group =='OCD'),]
cor.test(valuesOCD$DIFFCONF, valuesOCD$OCI)

Correlation<- ggplot (valuesOCD, aes(x = DIFFCONF , y= OCI))+
              geom_smooth ( method = lm,  colour ="red3")+
              geom_point(size = 2.5, pch =21, colour="darkred", fill="red3")+
              scale_y_continuous(limits=c(0, 63), breaks=seq(0,63,10))+
              scale_x_continuous(limits=c(-0.1, 0.2))+
        labs(y = "OCD symptom severity (OCI-R)")+
        labs(x = "Regression coefficients beta weights \n coupling confidence-action")+
        theme(axis.text=element_text(size=12 ),
        axis.title=element_text(size=12, face="bold"), 
        legend.text=element_text(size=12),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"), 
        legend.position="none")

figname= paste(outputdir, 'Plot_Correlation.tiff')
tiff(file = figname, width = 4, height = 3.5, units = "in", res = 300)
Correlation
dev.off()
#'Lasius neglectus workers- checking survival to oxidative stress using paraquat
#'survival checked everyday for 14 days.
#'3 colonies Lane23,Lane24,Lane26, split for 4 weeks. 
#'Each queenright colony had random number of queens


library(readxl)
library(e1071)
library(survival)
library(survminer)
library(coxme)

Lane_surv2<- read_excel("SupplementaryTable5_OxStressParaquat_Data_Results.xlsx", 
                        sheet = "RawData_Survival", col_types = c("text", 
                                                           "text", "text", "text", "numeric", 
                                                           "numeric"))
Lane_surv2<-as.data.frame(Lane_surv2)
#View(Lane_surv)
table(Lane_surv2$Colony,Lane_surv2$Queen)
table(Lane_surv2$Location,Lane_surv2$Treatment,Lane_surv2$Queen)
# REMOVE NAs

Lane_surv2<-na.exclude(Lane_surv2)###removed na
summary(Lane_surv2)

#Full model:
Lane2_Survobject<-Surv(Lane_surv2$Day,Lane_surv2$Status)
Lane2me1<-coxme(Lane2_Survobject~
                  Location*Queen+Location*Treatment+Queen*Treatment+(1|Colony),data=Lane_surv2)
summary(Lane2me1)
car::Anova(Lane2me1)
ranef(Lane2me1)

#Just Water:
LaneW<-Lane_surv2[Lane_surv2$Treatment=="Control",]
Lane_Survobject_W<-Surv(LaneW$Day,LaneW$Status)
COXLane_W<-coxme(Lane_Survobject_W~Queen*Location+(1|Colony),data=LaneW)
summary(COXLane_W)
car::Anova(COXLane_W)
#Does a TYPE II Anova (where order doesn't matter and main effect factors are tested against each other.
ranef(COXLane_W)

#Paraquat:
LaneP<-Lane_surv2[Lane_surv2$Treatment!="Control",]
Lane_Survobject_P<-Surv(LaneP$Day,LaneP$Status)
COXLane_P<-coxme(Lane_Survobject_P~Queen*Location+(1|Colony),data=LaneP)
summary(COXLane_P)
car::Anova(COXLane_P)
ranef(COXLane_P)

####PLOTS###
#scale_color_manual(values = c("Young" = "mediumseagreen", "Old" = "mediumpurple4"))
fitLaneL<-survfit(Surv(Day,Status)~Location, data=Lane_surv2)
fitLaneL

survcurvLnegW<-ggsurvplot_facet(fitLaneL,data=Lane_surv2,#surv.median.line="hv",
                                  palette=c("mediumseagreen","mediumpurple4"),
                                  ylab="Proportion of workers alive", #legend="none",
                                  size=3,
                                  xlab="Day",size=2,ggtheme = theme_bw()+
                                    theme(plot.title = element_text(hjust = 0.5, face = "bold")),
                                  break.time.by=1.5,font.title=25,
                                  facet.by =c("Treatment"),
                                  panel.labs=list(Location=c("Control","Paraquat")),
                                  panel.labs.font.x= list(face = "bold", color = NULL, size = 14,angle = NULL),
                                  panel.labs.font.y= list(face = "bold", color = NULL, size = 14,angle = NULL),
                                  legend.labs=c("Inside workers","Outside workers"),
                                  #legend = c(0.5,0.5),
                                  legend="top",
                                  font.tickslab=c(14,"bold","black"),font.legend=20,
                                  xlim=c(0,15),
                                  ylim=c(0,1),font.y=c(20,"bold","black"),
                                  panel.labs.background = list(colour="black",fill="white"),
                                  legend.title="Types of workers:",
                                  font.x=c(20,"bold","black"))#,censor.shape="X", censor.size = 8)

survcurvLnegW


fitLaneQ<-survfit(Surv(Day,Status)~Queen, data=Lane_surv2)
fitLaneQ
#Blue- "#6699CC", Red "#D73027"
survcurvLnegQ<-ggsurvplot_facet(fitLaneQ,data=Lane_surv2,#surv.median.line="hv",
                                palette=c("firebrick","royalblue"),
                                ylab="Proportion of workers alive", #legend="none",
                                size=3,
                                xlab="Day",size=2,ggtheme = theme_bw()+
                                  theme(plot.title = element_text(hjust = 0.5, face = "bold")),
                                break.time.by=1.5,font.title=25,
                                facet.by =c("Treatment"),
                                panel.labs=list(Location=c("Control","Paraquat")),
                                panel.labs.font.x= list(face = "bold", color = NULL, size = 14,angle = NULL),
                                panel.labs.font.y= list(face = "bold", color = NULL, size = 14,angle = NULL),
                                legend.labs=c("Queen presence","Queen absence"),
                                #legend = c(0.5,0.5),
                                legend="top",
                                font.tickslab=c(14,"bold","black"),font.legend=20,
                                xlim=c(0,15),
                                ylim=c(0,1),font.y=c(20,"bold","black"),
                                panel.labs.background = list(colour="black",fill="white"),
                                legend.title="Types of workers:",
                                font.x=c(20,"bold","black"))#,censor.shape="X", censor.size = 8)

survcurvLnegQ

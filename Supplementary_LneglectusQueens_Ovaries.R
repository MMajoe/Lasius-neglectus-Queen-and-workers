#ovaries-Queen Lasius neglectus
#' firstly most of the queens from our metadata are not photographed 
#' and many don't have good enough photos or even scales so 

library(readxl)
library(car)
library(lme4)
library(ggplot2)
library(glmm)
#library(emmmeans)
library(moments)
library(fitdistrplus)
library(betareg)

Lneg_Q_Ov_mm <- read_excel("Results_OvariesAvil.xlsx", sheet="ToUse",
                            col_types = c("text", "text", 
                                          "numeric", "text","numeric","text"))
Lneg_Q_Ov_mm<-as.data.frame(Lneg_Q_Ov_mm)
mynamestheme<-theme(plot.title=element_text(face="bold",size=(28)),
                    legend.title=element_text(colour="black",face="bold"),
                    legend.text=element_text(colour="black",face="italic"),
axis.title=element_text(colour="black",face="bold",size = 24),
axis.text=element_text(colour="black",face="italic",size=24))



hist(Lneg_Q_Ov_mm$MeanOvLength)
hist(log10(Lneg_Q_Ov_mm$MeanOvLength))
skewness(Lneg_Q_Ov_mm$MeanOvLength)
skewness(log10(Lneg_Q_Ov_mm$MeanOvLength))
descdist(log10(Lneg_Q_Ov_mm$MeanOvLength),discrete=FALSE)
#cor(Lneg_Q_Ov_mm$Age_Group,log10(Lneg_Q_Ov_mm$MeanOvLength),use='everything',method='spearman')

m2.1<-lm(log10(MeanOvLength)~Age_Group,data=Lneg_Q_Ov_mm)
Anova(m2.1) #
#Not significantly associated with age

m2plot_1<-ggplot(data=Lneg_Q_Ov_mm,aes(x=reorder(Age_Group,-MeanOvLength),y=log10(MeanOvLength),
                                       color=Age_Group))+
  #x=reorder(Age_Group,-MeanOvLength) reorders the boxplot so young is before old
  geom_boxplot(size=2)+scale_color_manual(#breaks = Lneg_Q_Ov_mm$Age_Group,
                                   values = c(  "grey70","#56B4E9"))+ #,"#56B4E9","#0072B2","#E69F00"
  geom_point(size=2,shape=19,colour='black')+theme_bw()+
  labs(x="Queen age",y='log10(Average ovariole length-mm)')
#, 
  #     title='Linear correlation between queen ovariole length and fungal volume')+
  #theme(plot.title = element_text(hjust=0.5,face='bold',size=16))
m2plot_1+mynamestheme


###FOR THESIS####
#Decide to randomly keep the first 5 available YOUNG (WITH SCALE) from excel and all 5 old (just to check)
Lneg_Subset<-Lneg_Q_Ov_mm[Lneg_Q_Ov_mm$Consider=="KEEP",]
hist(Lneg_Subset$MeanOvLength)
hist(log10(Lneg_Subset$MeanOvLength))
skewness(Lneg_Subset$MeanOvLength)
skewness(log10(Lneg_Subset$MeanOvLength))
descdist(log10(Lneg_Subset$MeanOvLength),discrete=FALSE)
#cor(Lneg_Subset$Age_Group,log10(Lneg_Subset$MeanOvLength),use='everything',method='spearman')

m2.2<-lm(log10(MeanOvLength)~Age_Group,data=Lneg_Subset)
summary(m2.2)
Anova(m2.1) #
#Not significantly associated with age
m_B<-betareg(log10(MeanOvLength)~Age_Group,data=Lneg_Subset)
summary(m_B)
Anova(m_B)

m22plot_1<-ggplot(data=Lneg_Subset,aes(x=reorder(Age_Group,-MeanOvLength),y=log10(MeanOvLength),
                                       color=Age_Group))+
  #x=reorder(Age_Group,-MeanOvLength) reorders the boxplot so young is before old
  geom_boxplot(size=3)+scale_color_manual(#breaks = Lneg_Subset$Age_Group,
    values = c(  "grey70","firebrick"))+ #,"#56B4E9","#0072B2","#E69F00"
  geom_point(size=4,shape=23,colour='black')+theme_classic(base_line_size=2)+
  labs(x="Queen age",y='log10(Average ovariole length-mm)')+
  theme(legend.position = "none")
  #theme(legend.text = element_text(NULL))
#, 
#     title='Linear correlation between queen ovariole length and fungal volume')+
#theme(plot.title = element_text(hjust=0.5,face='bold',size=16))
m22plot_1+mynamestheme

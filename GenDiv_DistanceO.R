rm(list=ls())

library(gridExtra)
library(dplyr)

setwd("C:/Users/thais/Documents/R/Data/TenchQC_Ch3/GenDiv/")

#Import data
IndGenDiv <- read.csv("IndGenDiv.csv",sep=",")
IndGenDiv <- IndGenDiv[1:195,]
LCWdist <- read.csv("../Geospatial/LCW_distance.csv",sep=" ")
ODist <- LCWdist[-1,1] #to keep pairwise dist between each sample and "origin"
IBD_GD <- readRDS("../Popstructure/sGD/IBD_GD_rad220.rds")
mycoord <- read.csv("../Final_Data/Final_Coordinates.csv")

#Create a column specifying whether it was North or South (Champlain) of putative origin
mycoord2 <- mycoord %>%
  mutate(InvDir = ifelse(Loc %in% c("LKChamplain","Richelieu") &
                           Latitude<mycoord$Latitude[1],
                         "South","North"))
InvDir <- mycoord2$InvDir[2:196]
GenDiv_ODist <- cbind(IndGenDiv,ODist,InvDir,IBD_GD)

#Produce some summary stats
library(dplyr)
GenDiv_ODist2[,c(14,16,17)] %>% summarise_if(is.numeric, mean)
GenDiv_ODist2[,c(14,16,17)] %>% summarise_if(is.numeric, min)
GenDiv_ODist2[,c(14,16,17)] %>% summarise_if(is.numeric, max)

#remove the influentional point (127)
GenDiv_ODist2 <- GenDiv_ODist[-c(1,126),]

###################
#Ind-level metrics#
###################

#Model fitting and selection for IR
lm_Inrel2 <- lm(GenDiv_ODist2$Inrel~GenDiv_ODist2$InvDir*GenDiv_ODist2$ODist)
lm_Inrel2b <- lm(GenDiv_ODist2$Inrel~as.factor(GenDiv_ODist2$InvDir)+GenDiv_ODist2$ODist)
anova(lm_Inrel2,lm_Inrel2b)# the interaction is significant
summary(lm_Inrel2)
plot(lm_Inrel2)
confint_Inrel2 <- confint(lm_Inrel2) 

#Model fitting and selection for MLH
lm_MLH2 <- lm(GenDiv_ODist2$MLH~GenDiv_ODist2$ODist*GenDiv_ODist2$InvDir)
lm_MLH2b <- lm(GenDiv_ODist2$MLH~GenDiv_ODist2$ODist+GenDiv_ODist2$InvDir) 
anova(lm_MLH2,lm_MLH2b) # no sig difference
lm_MLH2b1 <- lm(GenDiv_ODist2$MLH~GenDiv_ODist2$ODist) 
lm_MLH2b2 <- lm(GenDiv_ODist2$MLH~GenDiv_ODist2$InvDir) 
lm_MLH2b3 <- lm(GenDiv_ODist2$MLH~1)
anova(lm_MLH2b,lm_MLH2b1)
anova(lm_MLH2b,lm_MLH2b2)
anova(lm_MLH2b1,lm_MLH2b3)
plot(lm_MLH2b1)
summary(lm_MLH2b1)

#############################
#Neighbourhood-level metrics#
#############################

Ar <- lm(GenDiv_ODist2$Ar~GenDiv_ODist2$InvDir*GenDiv_ODist2$ODist)
Arb <- lm(GenDiv_ODist2$Ar~as.factor(GenDiv_ODist2$InvDir)+GenDiv_ODist2$ODist)
anova(Ar,Arb)# the interaction is significant
summary(Ar)
plot(Ar)
confint_Ar <- confint(Ar) 

Hs <- lm(GenDiv_ODist2$Hs~GenDiv_ODist2$InvDir*GenDiv_ODist2$ODist)
Hsb <- lm(GenDiv_ODist2$Hs~as.factor(GenDiv_ODist2$InvDir)+GenDiv_ODist2$ODist)
anova(Hs,Hsb)# the interaction is significant
summary(Hs)
plot(Hs)
confint_Hs <- confint(Hs) 

Ho <- lm(GenDiv_ODist2$Ho~GenDiv_ODist2$InvDir*GenDiv_ODist2$ODist)
Hob <- lm(GenDiv_ODist2$Ho~as.factor(GenDiv_ODist2$InvDir)+GenDiv_ODist2$ODist)
anova(Ho,Hob)# the interaction is significant
summary(Ho)
plot(Ho)
confint_Ho <- confint(Ho)

GenDiv_ODist2$FIS <- 1-(GenDiv_ODist2$Ho/GenDiv_ODist2$Hs)
FIS <- lm(GenDiv_ODist2$FIS~GenDiv_ODist2$InvDir*GenDiv_ODist2$ODist)
FISb <- lm(GenDiv_ODist2$FIS~as.factor(GenDiv_ODist2$InvDir)+GenDiv_ODist2$ODist)
anova(FIS,FISb)# the interaction is significant
summary(FIS)
plot(FIS)
confint_FIS <- confint(FIS)

###########
#Make plot#
###########
library(RgoogleMaps)
library(ggmap)

Lat <- c(min(GenDiv_ODist2$Y)-0.23,max(GenDiv_ODist2$Y)+0.25)
Long <- c(min(GenDiv_ODist2$X)-0.23,max(GenDiv_ODist2$X)+0.25)
center <- c(mean(Long),mean(Lat))
#for the following, need to register for an API key, and make sure that billing + some of the APIs are enables (geolocation, openmaps,and a couple of others).
register_google(key="AIzaSyB4L1wVSxspit6FPCSRKOhHmvh6B86Ky0o") 
back <- get_googlemap(center=center,zoom=7,maptype="satellite")

#IR
IRmap <- ggmap(back)+
  scale_x_continuous(limits=c(-75,-72))+
  scale_y_continuous(limits=c(43.5,46.5))+
  geom_point(data=GenDiv_ODist2,aes(x=X,y=Y,fill=Inrel),size=3,shape=21)+
  scale_fill_gradient2("IR",low="blue",mid="white",high="red",midpoint=mean(GenDiv_ODist2$Inrel))+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  geom_point(data=mycoord[1,],aes(x=Longitude,y=Latitude),shape=18,size=4)

IRmap+
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8))
  guides(color = guide_legend(override.aes = list(size = 0.5)))

#MLH
MLHmap <- ggmap(back)+
  scale_x_continuous(limits=c(-75,-72))+
  scale_y_continuous(limits=c(43.5,46.5))+
  geom_point(data=GenDiv_ODist2,aes(x=X,y=Y,fill=MLH),size=3,shape=21)+
  scale_fill_gradient2(low="red",mid="white",high="blue",midpoint=mean(GenDiv_ODist2$MLH))+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  geom_point(data=mycoord[1,],aes(x=Longitude,y=Latitude),shape=18,size=4)

pdf("Map_indGenDIv.pdf")
grid.arrange(IRmap,MLHmap,Armap,Homap,ncol=2,nrow=2)
dev.off()

#allelic richness
Armap <- ggmap(back)+
  scale_x_continuous(limits=c(-75,-72))+
  scale_y_continuous(limits=c(43.5,46.5))+
  geom_point(data=GenDiv_ODist2,aes(x=X,y=Y,fill=Ar),size=3,shape=21)+
  scale_fill_gradient2(low="red",mid="white",high="blue",midpoint=mean(GenDiv_ODist2$Ar))+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  geom_point(data=mycoord[1,],aes(x=Longitude,y=Latitude),shape=18,size=4)

#observed heterozygosity
Homap <- ggmap(back)+
  scale_x_continuous(limits=c(-75,-72))+
  scale_y_continuous(limits=c(43.5,46.5))+
  geom_point(data=GenDiv_ODist2,aes(x=X,y=Y,fill=Ho),size=3,shape=21)+
  scale_fill_gradient2(low="red",mid="white",high="blue",midpoint=mean(GenDiv_ODist2$Ho))+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  geom_point(data=mycoord[1,],aes(x=Longitude,y=Latitude),shape=18,size=4)

#FIS
FISmap <- ggmap(back)+
  scale_x_continuous(limits=c(-75,-72))+
  scale_y_continuous(limits=c(43.5,46.5))+
  geom_point(data=GenDiv_ODist2,aes(x=X,y=Y,fill=FIS),size=3,shape=21)+
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=mean(GenDiv_ODist2$FIS))+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  geom_point(data=mycoord[1,],aes(x=Longitude,y=Latitude),shape=18,size=4)

pdf("Map_indGenDIv.pdf")
grid.arrange(IRmap,MLHmap,Armap,Homap,ncol=2,nrow=2)
dev.off()

pdf("Map_indGenDIv2.pdf")
grid.arrange(FISmap,nrow=2,ncol=2)
dev.off()




P1 <- ggplot(data=GenDiv_ODist2,aes(x=ODist,y=Inrel))+
  ylab("Internal Relatedness")+ 
  xlab("")+
  #Add points
  geom_point(aes(fill=InvDir),size=2,shape=21)+
  scale_fill_manual(values=c("North"="#4B0055","South"="orange"))+
  #North 
  geom_abline(slope = coef(lm_Inrel2)[[3]],intercept = coef(lm_Inrel2)[[1]],colour="#4B0055",size=2)+ 
  #South
  geom_abline(slope = coef(lm_Inrel2)[[4]]+coef(lm_Inrel2)[[3]], intercept = coef(lm_Inrel2)[[2]]+coef(lm_Inrel2)[[1]],colour="orange",size=2)+ #South
  #Esthetic stuff
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14,face="bold"),
        axis.text= element_text(size=10,colour="black"),
        legend.key = element_blank(),
        legend.text=element_text(size=10))+
  guides(fill=guide_legend(override.aes=list(size=5)))

P2 <- ggplot(data=GenDiv_ODist2,aes(x=ODist,y=MLH))+
  geom_point(aes(fill=InvDir),size=2,shape=21)+
  scale_fill_manual(values=c("North"="#4B0055","South"="orange"))+
  ylab("Multilocus Heterozygosity (%)")+
  xlab("Distance from the origin (Km)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14,face="bold"),
        axis.text= element_text(size=10,colour="black"),
        legend.key = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=10))+
  guides(fill=guide_legend(override.aes=list(size=5)))

P3 <- ggplot(data=GenDiv_ODist2,aes(x=ODist,y=Ar))+
  geom_point(aes(fill=InvDir),size=2,shape=21)+
  scale_fill_manual(values=c("North"="#4B0055","South"="orange"))+
  ylab("Allelic richness")+
  xlab("Distance from the origin (Km)")+
  geom_abline(slope = coef(Ar)[[3]],intercept = coef(Ar)[[1]],colour="#4B0055",size=2)+ 
  geom_abline(slope = coef(Ar)[[4]]+coef(Ar)[[3]], intercept = coef(Ar)[[2]]+coef(Ar)[[1]],colour="orange",size=2)+ #South
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14,face="bold"),
        axis.text= element_text(size=10,colour="black"),
        legend.key = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=10))+
  guides(fill=guide_legend(override.aes=list(size=5)))

P4 <- ggplot(data=GenDiv_ODist2,aes(x=ODist,y=Ho))+
  geom_point(aes(fill=InvDir),size=2,shape=21)+
  scale_fill_manual(values=c("North"="#4B0055","South"="orange"))+
  ylab("Observed heterozygosity")+
  xlab("Distance from the origin (Km)")+
  geom_abline(slope = coef(Ho)[[3]],intercept = coef(Ho)[[1]],colour="#4B0055",size=2)+ 
  geom_abline(slope = coef(Ho)[[4]]+coef(Ho)[[3]], intercept = coef(Ho)[[2]]+coef(Ho)[[1]],colour="orange",size=2)+ #South
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14,face="bold"),
        axis.text= element_text(size=10,colour="black"),
        legend.key = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=10))+
  guides(fill=guide_legend(override.aes=list(size=5)))

FISp <- ggplot(data=GenDiv_ODist2,aes(x=ODist,y=FIS))+
  geom_point(aes(fill=InvDir),size=2,shape=21)+
  scale_fill_manual(values=c("North"="#4B0055","South"="orange"))+
  ylab("FIS")+
  xlab("Distance from the origin (Km)")+
  geom_abline(slope = coef(FIS)[[3]],intercept = coef(FIS)[[1]],colour="#4B0055",size=2)+ 
  geom_abline(slope = coef(FIS)[[4]]+coef(FIS)[[3]], intercept = coef(FIS)[[2]]+coef(FIS)[[1]],colour="orange",size=2)+ #South
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14,face="bold"),
        axis.text= element_text(size=10,colour="black"),
        legend.key = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=10))+
  guides(fill=guide_legend(override.aes=list(size=5)))

pdf("GenDiv_cline.pdf")
grid_arrange_shared_legend(P1,P2,nrow=2,ncol=1,position="top")
dev.off()









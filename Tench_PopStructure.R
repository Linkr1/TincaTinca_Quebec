rm(list=ls())

library(vcfR)
library(adegenet)
library(vegan)
library(memgene)

library(dplyr)


#for plotting
library(grDevices)
library(ggplot2)
library(ggmap)
library(gridExtra)

#for gendiv
library(hierfstat)

setwd("C:/Users/thais/Documents/R/Data/TenchQC_Ch3/Popstructure/")

####import and prepare data####
####import and prepare data####
#import the vcf and convert to genind for Adegenet
vcfRobj <- read.vcfR("../Final_Data/Final_Filtered.recode.vcf",convertNA = TRUE)
genindobj <- vcfR2genind(vcfRobj,return.alleles = TRUE)

#Download coordinates
mycoord <- read.csv("../Final_Data/Final_Coordinates.csv")
mycoord <- mycoord[-1,]
indNames(genindobj)

#Check that the IDs are in the same order as in the coord files. 
all.equal(mycoord$Sample,indNames(genindobj))
#if not; 
which(mycoord$Sample!=indNames(genindobj))

#Attach coordinates to the genind object.

genindobj@other$xy <- mycoord[,c(2,3)]
genindobj@pop <- as.factor(mycoord[,5])
genindobj #now @other has the xy coordinates and the sampling locations!

#retrieve allelic freq and replace missing value (will need it for PCA and DAPC)
x.dat <- tab(genindobj,freq=TRUE,NA.method="mean")


######################################
# Principal Component Analysis (PCA) #
######################################

# Run PCA  --------------------------
pca.dat <- dudi.pca(x.dat,center=TRUE,scale=FALSE)
##decided to go with 3; 
barplot(pca.dat$eig,main="PCA eigenvalues",col=spectral(length(pca.dat$eig)))

pca.dat <- dudi.pca(df = x.dat, center = TRUE, scale = FALSE, scannf = FALSE, nf = 3)

#look at how many PCs have eigenvalues > 1
barplot(pca.dat$eig,ylab="PC's eigenvalues")
abline(h=1,col="red",lwd=2) # many!

#Look at how many PCs explain more than 1% of the variance
percent <- pca.dat$eig/sum(pca.dat$eig)*100
barplot(percent,ylab="Variance proportion (%) explained by PCs")
abline(h=1,col="red",lwd=2) #many!

#allele contrib
PCA_allele <- loadingplot(pca.dat$c1^2,thres=0.003)
#check if the alleles overlap with DAPC


# Plot PCA  --------------------------------
#create a df with ind coordinates
ind_coords <- as.data.frame(pca.dat$li)
#add col with site ID (reorganize arranging to dispersal path)
ind_coords$Site <- factor(genindobj$pop)
#calculate the average position for each population (centroid)
centroid <- aggregate(cbind(Axis1,Axis2,Axis3)~Site,data=ind_coords,FUN=mean)
#add coordinates to df
ind_coords <- left_join(ind_coords,centroid,by="Site",suffix=c("",".cen"))
#create color palette. Using a diverging one here.
cols <- hcl.colors(nPop(genindobj),palette="viridis")
#create custom axes labels
xlab <- paste("Axis 1 (",format(round(percent[1],1),nsmall=1), "%)",sep="")
ylab <- paste("Axis 2 (",format(round(percent[2],1),nsmall=1), "%)",sep="")

PCAplot <- ggplot(data=ind_coords,aes(x=Axis1,y=Axis2))+
  #adds the spiders
  geom_segment(aes(xend=Axis1.cen,yend=Axis2.cen,colour=Site),show.legend=FALSE)+
  #adds axes - careful with order!
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  #pch 21 is the one with both fill and borders
  geom_point(aes(fill=Site),shape=21,size=3,show.legend=FALSE)+
  geom_label(data=centroid,aes(label=Site,fill=Site),size=3,show.legend=FALSE)+
  #specify colours for both lines and points
  scale_fill_manual(values=cols)+
  scale_colour_manual(values=cols)+
  labs(x=xlab,y=ylab)+
  ggtitle("PCA")+
  theme(axis.text.y=element_text(colour="black",size=16),
        axis.text.x=element_text(colour="black",size=16),
        axis.title=element_text(colour="black",size=20),
        panel.border=element_rect(colour="black",fill=NA,size=1),
        panel.background=element_blank(),
        plot.title=element_text(colour="black",size=24,face="bold",hjust=0.5))

pdf("PCA_Tench.pdf")
PCAplot
dev.off() 


########
# DAPC #
########

# 1. DAPC with a-priori groups---------------------------------------------
#do cross validation to select the ideal number of PCs
xval <- xvalDapc(x.dat,genindobj$pop,result="groupMean",xval.plot=TRUE)
xval[5] #PCs with best stats; lower score=better
xval[4]
xval[6]
numPCs <- as.numeric(xval[6])
#run DAPC, using site ID as prior
dapc_prior <- dapc(x.dat,genindobj$pop,n.pca=numPCs,n.da=3)
#calculate euclidian distance between groups
dist(dapc_prior$grp.coord,method="euclidian")
#analyze proportion explained by axes
percent2 <- dapc_prior$eig/sum(dapc_prior$eig)*100
barplot(percent2,ylab= " % Genetic variance explained by eigenvectors ")


# 2. DAPC, no prior -------------------------------------------------
grp_man <- find.clusters(x.dat,n.pca=200,max.n.clust=10)
#look at the plot. BIC increases almost linearly from 1-10!
grp_min <- find.clusters(x.dat,criterion=c("min"),n.pca=200,max.n.clust=10,choose.n.clust=FALSE)
grp_min
#according to this criteria, best K=1
grp <- find.clusters(x.dat,criterion=c("diffNgroup"),n.pca=200,max.n.clust=10,choose.n.clust=FALSE)
grp
#according to this criteria, best K=2. BUT it is not designed to detect K=1!

xval2 <- xvalDapc(x.dat,grp$grp,result="groupMean",xval.plot=TRUE)
xval2[5] #PCs with best stats; lower score=better
xval2[4]
xval2[6]
numPCs2 <- as.numeric(xval2[6])

dapc_noprior <- dapc(x.dat,grp$grp,n.pca=numPCs2)
scatter(dapc_noprior)

#euclidian distance between clusters
dist(dapc_noprior$grp.coord,method="euclidian")

#look at the individuals assigned to the clusters; 
dapc_noprior$posterior

# 3. plot DAPC ------------------------
ind_coords2 <- as.data.frame(dapc_noprior$ind.coord)
#rename columns
colnames(ind_coords2) <- c("Axis1","Axis2","Axis3")
#add col with individuals
ind_coords2$Ind <- indNames(genindobj)
#add col with site ID (reorganize arranging to dispersal path)
ind_coords2$Site <- factor(genindobj$pop)
#add col with group ID 
ind_coords2$grp <- factor(grp$grp)
#calculate the average position for each grp (centroid)
centroid2 <- aggregate(cbind(Axis1,Axis2,Axis3)~grp,data=ind_coords2,FUN=mean)
#add coordinates to df
ind_coords2 <- left_join(ind_coords2,centroid2,by="grp",suffix=c("",".cen"))
#create color palette. Using a diverging one here (same as above).
cols <- hcl.colors(nPop(genindobj),palette="viridis")
#create custom axes labels
percent2 <- dapc_noprior$eig/sum(dapc_noprior$eig)*100
xlab2 <- paste("Axis 1 (",format(round(percent2[1],1),nsmall=1), "%)",sep="")
ylab2 <- paste("Axis 2 (",format(round(percent2[2],1),nsmall=1), "%)",sep="")

pdf("DAPC_TenchQC.pdf")
ggplot(data=ind_coords2,aes(x=Axis1,y=Axis2))+
  #adds the spiders
  geom_segment(aes(xend=Axis1.cen,yend=Axis2.cen,colour=Site),show.legend=FALSE)+
  #adds axes - careful with order!
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  #pch 21 is the one with both fill and borders
  geom_point(aes(fill=Site),shape=21,size=3,show.legend=TRUE)+
  geom_label(data=centroid2,aes(label=grp),size=3,show.legend=FALSE)+
  #specify colours for both lines and points
  scale_fill_manual(values=cols)+
  scale_colour_manual(values=cols)+
  labs(x=xlab2,y=ylab2)+
  ggtitle("DAPC")+
  theme(axis.text.y=element_text(colour="black",size=12),
        axis.text.x=element_text(colour="black",size=12),
        axis.title=element_text(colour="black",size=12),
        plot.title=element_text(colour="black",size=16,face="bold",hjust=0.5),
        panel.border=element_rect(colour="black",fill=NA,size=1),
        panel.background=element_blank())
dev.off()


########
# sPCA #
########

# FOr this analysis, we only want the contemporary samples. So, remove the last 8 rows of the genind object.
removeInd <- c("TTo121_210915","TTo122_210915", "TTo123_210915" ,"TTo124_210915" ,"TTo127_210915", "TTo128_210915" ,"TTo129_210915" ,"TTo130_210915")
# 1. Transform geo coordinates into Cartesian coordinates ----------------------
genindobj2 <- genindobj[indNames(genindobj)[1:195]]
#load the geographic distance between samples(in geospatial script)
geodist <- read.table("../Geospatial/LCW_distance.csv", sep=" ")
geodist <- geodist[-1,-1] #remove the "Origin" row

#convert distances to cartesian coordinates
cartcoord <- metaMDS(geodist,k=2,trymax=1000)

#check the convergence
cartcoord #there is low stress (<0.01)= low disagreement between 2D conf and predicted values)
stressplot(cartcoord,pch=1) #relationship is 1-1

#check the cartesian coordinates of the samples 
ordiplot(cartcoord)

#Convert those coordinates to df
cartxy <- as.data.frame(cartcoord$points)

#save, so you dont have to rerun it all the time. 
write.table(cartxy,file="cartxy.csv",quote=FALSE,col.names=TRUE)

#introduce some jitter to avoid duplicates
cartxy$MDS2 <- jitter(cartxy$MDS2)
cartxy$MDS1 <- jitter(cartxy$MDS1)

#attach to genindobj
genindobj2_cartxy <- genindobj2
genindobj2_cartxy@other$xy <- cartxy


# 2. Run sPCA using the Neighbourhood-by-distance method -------------------------------
mySpca <- spca(genindobj2_cartxy,type=5,d1=0,d2="dmin") 
#This takes a little while (~10mn)

#I chose 2 positive and three negative axes.
plot(mySpca)
summary(mySpca)
str(mySpca)

#test local and global structure
tab.x <- tab(genindobj2_cartxy,freq=T,NA.method="mean")#replace missing value
myGtest <- global.rtest(tab.x,mySpca$lw,nperm=999)
myGtest 
plot(myGtest)#not significant (which suggests the existence of no global spatial structure)

myLtest <- local.rtest(tab.x,mySpca$lw,nperm=999)
myLtest 
plot(myLtest) #not significant

#run PCA with the correct number of axes.
mySpca2 <- spca(genindobj2_cartxy,type=5,d1=0,d2="dmin") 
#keep 3 global and 3 local

#check results
plot(mySpca2)
summary(mySpca2)
head(mySpca2)


# 3. Plot sPCA results --------------------

#plot eigenvalues
barplot(mySpca2$eig,main="sPCA eigenvalues",col=spectral(length(mySpca2$eig)))
legend("bottom",fill=spectral(2),leg=c("Global structures","Local structures"),cex=0.75,bty="n",x.intersp=0.5,horiz=T)
abline(h=0,col="black")
screeplot(mySpca2)

#this can be used to plot on cartesian coordinates. 
plot(genindobj2_cartxy@other$xy)
s.value(genindobj2_cartxy@other$xy,mySpca2$ls[,1],add.p=T,cleg=0) 

#now I want to plot the loadings on geographic coordinates (easier to interpret)
#make data frame
loadings <- as.data.frame(mySpca2$ls)
#rename col 
colnames(loadings) <- c("Axis1","Axis2","Axis3","Axis188","Axis189","Axis190")
#add col with individuals
loadings$ID <- row.names(mySpca2$ls) 
#add col with coordinates
loadings$Lat <- (genindobj2$other$xy)[,1]
loadings$Long <- (genindobj2$other$xy)[,2]

#create background map
library(RgoogleMaps)
QCrange <- c(left=min(loadings$Long)-0.23,right=max(loadings$Long)+0.23,bottom=min(loadings$Lat)-0.23,top=max(loadings$Lat)+0.23)
Lat <- c(min(loadings$Lat)-0.23,max(loadings$Lat)+0.25)
Long <- c(min(loadings$Long)-0.23,max(loadings$Long)+0.25)
center <- c(mean(Long),mean(Lat))
#for the following, need to register for an API key, and make sure that billing + some of the APIs are enables (geolocation, openmaps,and a couple of others).
register_google(key="AIzaSyB4L1wVSxspit6FPCSRKOhHmvh6B86Ky0o") 
back <- get_googlemap(center=center,zoom=7,maptype="satellite")

#plot and export as pdf
mypal <- c("white","black")
mypal2 <- c("black","white")
p1 <- ggmap(back)+
  scale_x_continuous(limits=c(min(loadings$Long)-0.23,max(loadings$Long)+0.25))+
  scale_y_continuous(limits=c(min(loadings$Lat)-0.23,max(loadings$Lat)+0.25))+
  geom_jitter(data=loadings,aes(x=Long,y=Lat,size=abs(Axis1),fill=Axis1>0,colour=Axis1>0),shape=22,width=0.10,height=0.1,show.legend=FALSE)+
  scale_size_continuous(range=c(1,13))+
  scale_fill_manual(values=mypal)+
  scale_colour_manual(values=mypal2)+
  ggtitle("Axis 1 (global structure)")+
  labs(x="Longitude",y="Latitude")+
  theme(axis.text.y=element_text(colour="black",size=12),
        axis.text.x=element_text(colour="black",size=12),
        axis.title=element_text(colour="black",size=12),
        plot.title=element_text(colour="black",size=14,face="bold",hjust=0.5),
        panel.border=element_rect(colour="black",fill=NA,size=1),
        panel.background=element_blank())

p2 <- ggmap(back)+
  scale_x_continuous(limits=c(min(loadings$Long)-0.23,max(loadings$Long)+0.25))+
  scale_y_continuous(limits=c(min(loadings$Lat)-0.23,max(loadings$Lat)+0.25))+
  geom_jitter(data=loadings,aes(x=Long,y=Lat,size=abs(Axis2),fill=Axis2>0,colour=Axis2>0),shape=22,width=0.10,height=0.1,show.legend=FALSE)+
  scale_size_continuous(range=c(1,13))+
  scale_fill_manual(values=mypal)+
  scale_colour_manual(values=mypal2)+
  ggtitle("Axis 2 (global structure)")+
  labs(x="Longitude",y="Latitude")+
  theme(axis.text.y=element_text(colour="black",size=12),
        axis.text.x=element_text(colour="black",size=12),
        axis.title=element_text(colour="black",size=12),
        plot.title=element_text(colour="black",size=14,face="bold",hjust=0.5),
        panel.border=element_rect(colour="black",fill=NA,size=1),
        panel.background=element_blank())

p3 <- ggmap(back)+
  scale_x_continuous(limits=c(min(loadings$Long)-0.23,max(loadings$Long)+0.25))+
  scale_y_continuous(limits=c(min(loadings$Lat)-0.23,max(loadings$Lat)+0.25))+
  geom_jitter(data=loadings,aes(x=Long,y=Lat,size=abs(Axis3),fill=Axis3>0,colour=Axis3>0),shape=22,width=0.10,height=0.1,show.legend=FALSE)+
  scale_size_continuous(range=c(1,13))+
  scale_fill_manual(values=mypal)+
  scale_colour_manual(values=mypal2)+
  ggtitle("Axis 1 (global structure)")+
  labs(x="Longitude",y="Latitude")+
  theme(axis.text.y=element_text(colour="black",size=12),
        axis.text.x=element_text(colour="black",size=12),
        axis.title=element_text(colour="black",size=12),
        plot.title=element_text(colour="black",size=14,face="bold",hjust=0.5),
        panel.border=element_rect(colour="black",fill=NA,size=1),
        panel.background=element_blank())

p4 <- ggmap(back)+
  scale_x_continuous(limits=c(min(loadings$Long)-0.23,max(loadings$Long)+0.25))+
  scale_y_continuous(limits=c(min(loadings$Lat)-0.23,max(loadings$Lat)+0.25))+
  geom_jitter(data=loadings,aes(x=Long,y=Lat,size=abs(Axis188),fill=Axis188>0,colour=Axis188>0),shape=22,width=0.10,height=0.1,show.legend=FALSE)+
  scale_size_continuous(range=c(1,13))+
  scale_fill_manual(values=mypal)+
  scale_colour_manual(values=mypal2)+
  ggtitle("Axis -3 (local structure)")+
  labs(x="Longitude",y="Latitude")+
  theme(axis.text.y=element_text(colour="black",size=12),
        axis.text.x=element_text(colour="black",size=12),
        axis.title=element_text(colour="black",size=12),
        plot.title=element_text(colour="black",size=14,face="bold",hjust=0.5),
        panel.border=element_rect(colour="black",fill=NA,size=1),
        panel.background=element_blank())

p5 <- ggmap(back)+
  scale_x_continuous(limits=c(min(loadings$Long)-0.23,max(loadings$Long)+0.25))+
  scale_y_continuous(limits=c(min(loadings$Lat)-0.23,max(loadings$Lat)+0.25))+
  geom_jitter(data=loadings,aes(x=Long,y=Lat,size=abs(Axis189),fill=Axis189>0,colour=Axis189>0),shape=22,width=0.10,height=0.1,show.legend=FALSE)+
  scale_size_continuous(range=c(1,13))+
  scale_fill_manual(values=mypal)+
  scale_colour_manual(values=mypal2)+
  ggtitle("Axis -2 (local structure)")+
  labs(x="Longitude",y="Latitude")+
  theme(axis.text.y=element_text(colour="black",size=12),
        axis.text.x=element_text(colour="black",size=12),
        axis.title=element_text(colour="black",size=12),
        plot.title=element_text(colour="black",size=14,face="bold",hjust=0.5),
        panel.border=element_rect(colour="black",fill=NA,size=1),
        panel.background=element_blank())

p6 <- ggmap(back)+
  scale_x_continuous(limits=c(min(loadings$Long)-0.23,max(loadings$Long)+0.25))+
  scale_y_continuous(limits=c(min(loadings$Lat)-0.23,max(loadings$Lat)+0.25))+
  geom_jitter(data=loadings,aes(x=Long,y=Lat,size=abs(Axis190),fill=Axis190>0,colour=Axis190>0),shape=22,width=0.10,height=0.1,show.legend=FALSE)+
  scale_size_continuous(range=c(1,13))+
  scale_fill_manual(values=mypal)+
  scale_colour_manual(values=mypal2)+
  ggtitle("Axis -1 (local structure)")+
  labs(x="Longitude",y="Latitude")+
  theme(axis.text.y=element_text(colour="black",size=12),
        axis.text.x=element_text(colour="black",size=12),
        axis.title=element_text(colour="black",size=12),
        plot.title=element_text(colour="black",size=14,face="bold",hjust=0.5),
        panel.border=element_rect(colour="black",fill=NA,size=1),
        panel.background=element_blank())

#Export the first two (significant, global structure)
pdf("sPCA_TenchQC.pdf")
grid.arrange(p1,p2,ncol=2,nrow=2)
dev.off()
  
#####MEMGENE####

#will use genind2, but needs some tweaking first
gen <- as.data.frame(genindobj2@tab)

#rename the loci so the software does not get confused
nloc<- ncol(gen)/2
namdf <- rep(1:nloc,each=2)
suf <- rep(c("a","b"), length.out=nloc)
colnames(gen) <- paste0(namdf,suf)
colnames(gen)[1:10] # lookin' good <3

#produce a proportion of shared alleles genetic distance matrix
radialDM <- codomToPropShared(gen)

#get the spatial data in!
radialXY <- as.data.frame(cbind(genindobj2@other$xy[2],genindobj2@other$xy[1]))


###Export data to run on cluster. 
write.table(gen,file="gen_memgene.csv",quote=FALSE,col.names=TRUE)
write.table(radialXY,file="radialXY.csv",quote=FALSE,col.names=TRUE)



#This is a resistance raster (0=land, 1=water)
resmap <- raster("Reformat_Finalmap.asc")
resmap2 <- resmap
resmap2[resmap2==0] <- 10000
#Extract MEMGENE variables
radialAnalysis <- mgQuick(radialDM,radialXY)

mgMap(radialXY,radialAnalysis$memgene[,1:2])
radialAnalysis$RsqAdj #very small, which reflects that only a small proportion of the patterns are attributed



# to spatial patterns. 
TenchMEMGENEprop <- radialAnalysis$sdev/sum(radialAnalysis$sdev)
format(signif(TenchMEMGENEprop,3 )[1:3],scientific=FALSE)

#Now, text whether the resistance layer explains structure better. 
comparetwo <- mgLandscape(resmap2,radialDM,radialXY,euclid=TRUE,forwardPerm = 500, finalPerm=1000)

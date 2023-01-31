# ptdByAge.R

setwd('C:\\SSRS\\YBCU\\LCR\\2018_YBCU_LCR')

#excel file has more data but something wrong with it, thinks all the numbers are factors....
#library(xlsx)
#det<-read.xlsx("tdByAge.xlsx", sheetName="tdByAge",header=TRUE) 

#original dataset used in presentation
det <- read.table("tphByAge2.txt",header=TRUE,as.is=TRUE)
det <- subset(det, (ptd>-1) & (is.na(note) | (note != "presv")))
#det <- subset(det, !(Site %in% c('CNCRN','CVCA03','PVER01','PVER08'))) # remove the stupid sites or ones with little data
det$big<-ifelse(det$cha2>99,1,0) # make big depend on overall area of habitat the site is in 

#Territories by age
ttl <- "Cuckoo Density by Stand Age"
ylab<-"Terr./20 ha"
ylab<-"    Territories/20 ha"
xlab="Patch Age (Years)"
det$y <- sqrt(det$ptd) # sum(PRB, COB) - best?

#Detections by age
ttl <- "Cuckoo Detections by Stand Age"
ylab<-"Dets/ha"
xlab="Patch Age (Years)"
det$y <- sqrt(det$dph) # detections per ha

#Territories by year
ttl <- "Cuckoo Density by Year"
xlab="Year"
ylab<-"    Territories/Year"
det$y <- sqrt(det$ptd) 

#Powerpoint
cex.PPT=1.5
cex.lab.PPT=1.5

pcex=cex.PPT
pcex.lab=cex.lab.PPT
pcex.lab.title=cex.lab.PPT

#normal
pcex=1
pcex.lab=1
pcex.lab.title=1

# Variables:
# Area (factor)
# Site (factor)
# Year (factor)
# Age (of the patch in years since planted, integer)
# dph (detections per hectare, numeric)
# contig (is the patch contiguous with another patch?, integer)
# big (is the patch within a large planted area (>100 ha)?, integer)
# svdets (number of survey detections of cuckoos, integer)
# nests (number of nests found, integer)
# ha (size of the patch of that age in hectares, numeric)
# cha (size of the contiguous patch, numeric)
# cha2 (size of the contiguous patch - treating PVER as 1 contiguous patch, numeric)
# inc (include in analysis? 1=yes, factor)
# youngest (is this patch the youngest patch around?, factor)
# nAges (number of different ages of patches within the contiguous patch, integer)
# nAges2 (nAges - but treat ages > 5 years as the same age, integer)
# tp20h (estimated territories per 20 hectares)
# PO PR CO (possible, probable and confirmed breeding territories, integer)
# ptd (probable territory density: (PR + CO only)/ha * 20

str(det)
View(det)
summary(det)

wsize<-5
windows(wsize,wsize)		# Figure 1 (Plot dph vs Age)
par(las=1,oma=c(0,0,0,0),mar=c(2.5,2.5,2.5,1),mgp=c(1.1,0.1,0),tcl=0.2)
xlm <- c(0,max(det$Age))
plot(det$Age,det$dph,xlim=xlm,xlab="Age of Patch (years)",ylab="")
ttl <- "Shannon's Cuckoos"
mtext(side=3,line=0.1,adj=-0.08,"Detections/ha")
mtext(side=3,line=0.1,adj=1,ttl)
n <- nrow(det)
legend("topright",leg=paste("n:",n),bty="n")

windows(wsize,wsize)		# Figure 1.5 (Plot ptd vs Age)
par(las=1,oma=c(0,0,0,0),mar=c(2.5,2.5,2.5,1),mgp=c(1.1,0.1,0),tcl=0.2)
par(cex.axis=pcex, cex.lab=pcex, cex.main=pcex)
xlm <- c(0,max(det$Age))
ylm <- c(0,max(det$ptd))
plot(det$Age,det$ptd,xlim=xlm,xlab="Patch Age (Years)",ylab=ylab,cex=pcex,cex.axis=pcex.lab,cex.lab=pcex.lab)
mtext(side=3,line=0.1,adj=-0.08,ylab,cex=pcex.lab)
mtext(side=3,line=0.1,adj=1,ttl,cex=pcex.lab)
n <- nrow(det)
legend("topright",leg=paste("n:",n),bty="n",cex=pcex.lab)



#Repeat above with y...
xlm <- c(0,max(det$Age))
ylm <- c(0,max(det$ptd))
windows(wsize,wsize)		# Figure 1.5 (Plot tph vs Age)
par(las=1,oma=c(0,0,0,0),mar=c(2.5,2.5,2.5,1),mgp=c(1.1,0.1,0),tcl=0.2)
par(cex.axis=pcex, cex.lab=pcex, cex.main=pcex)
plot(det$Age,det$y,xlim=xlm,xlab="Patch Age (Years)",yaxt="n",ylab=ylab,cex=pcex,cex.axis=pcex.lab,cex.lab=pcex.lab)

#or just pver...
pver<-subset(det,Area=="PVER")
plot(pver$Age,pver$y,xlim=xlm,xlab="Patch Age (Years)",yaxt="n",ylab=ylab,cex=pcex,cex.axis=pcex.lab,cex.lab=pcex.lab)

laby <- c(0,5,10,20,40,60,80)/10
axis(side=2,at=sqrt(laby),lab=laby,cex=pcex)
mtext(side=3,line=0.1,adj=1,cex=pcex,ttl)
n <- nrow(det)
legend("topright",leg=paste("n:",n),bty="n",cex=pcex.lab)
mod <- lm(data=det,y~Age)
summary(mod) # P-value = 0.00211, Adj. r-squared = 0.08152
abline(mod$coef,col=2)

#-------------------------------------------------Pause 1
#
# Repeat with big
#
windows(wsize,wsize)		# Figure 2 (Plot sqrt(dph) vs Age)
par(las=1,oma=c(0,0,0,0),mar=c(2.5,3.0,2.5,0.5),mgp=c(1.1,0.1,0),tcl=0.2)
par(cex.axis=pcex, cex.lab=pcex, cex.main=pcex)
x <- det$Age+0.2*det$big-0.1
xlm <- c(0,max(det$Age))
pt <- ifelse(det$big==0,20,21)
plot(x,det$y,type="n",xlim=xlm,xlab="Patch Age (Years)",ylab=ylab,yaxt="n",cex=pcex,cex.axis=pcex.lab,cex.lab=pcex.lab)
points(x,det$y,pch=pt,cex=pcex)

laby <- c(0,5,10,20,40,60,80)/10
axis(side=2,at=sqrt(laby),lab=laby,cex=pcex)
#mtext(side=3,line=0.1,adj=-0.1,cex=pcex,ylab)
mtext(side=3,line=0.1,adj=1,cex=pcex,ttl)

det0b <- subset(det,big==0); n0 <- nrow(det0b)
det1b <- subset(det,big==1); n1 <- nrow(det1b)
lg1 <- paste("Small (n=",n0,")",sep="")
lg1 <- c(lg1,paste("Large (n=",n1,")",sep=""))
legend("topright",leg=lg1,pch=c(20,1),inset=c(0.01,0.01),bg="ivory",cex=pcex)

mod1b <- lm(data=det1b,y~Age)		# fit models to subgroups
summary(mod1b)
mod0b <- lm(data=det0b,y~Age)
summary(mod0b)

#abline(mod$coef,col=1)
abline(mod1b$coef,col=2)
abline(mod0b$coef,col=2,lty=2)

lg2 <- round(c(mod1b$coef[2],mod0b$coef[2]),4)
legend("bottomright",leg=lg2,inset=c(0.01,0.08),lty=c(1,2),col=2,
       bg="ivory",title="Slope",x.intersp=0.2,y.intersp=0.9,cex=1.1)

#-------------------------------------------------Pause 2
#
# 
c <- det$big

det$age0 <- ifelse(c==0,det$Age,0)
det$age1 <- ifelse(c==1,det$Age,0)

mod <- lm(data=det,y~big+age0+age1)	# fit model to all data
summary(mod)
rsq <- round(100*summary(mod)$adj.r.squared,1)

windows(wsize,wsize)		# Figure 3 (Q-Q plot of studentized residuals)
par(las=1,oma=c(0,0,0,0),mar=c(2.5,2.5,2.5,1),mgp=c(1.1,0.1,0),tcl=0.2)

library(MASS)
qqnorm(studres(mod),ylab="",pch=pt,main="",cex=pcex)
qqline(studres(mod),col=2)
mtext(side=3,line=0.1,adj=-0.1,"Studentized Residuals")
mtext(side=3,line=0.1,adj=1,ttl)
legend("topleft",leg=paste("r-sq: ",rsq,"%",sep=""),bty="n")

#-------------------------------------------------Pause 3

det$z <- studres(mod)
det2 <- subset(det,z<(-0.9)) # =outliers

windows(wsize,wsize)		# Figure 4 (Plot sqrt(dph) vs Age for three groups)
par(las=1,oma=c(0,0,0,0),mar=c(2.5,2.5,2.5,1),mgp=c(1.1,0.1,0),tcl=0.2)
x <- det$Age+0.2*det$big-0.1

pt <- ifelse(det$big==0,20,21)
plot(x,det$y,type="n",xlim=xlm,xlab="Patch Age (Years)",ylab=ylab,yaxt="n",cex=pcex,cex.axis=pcex.lab,cex.lab=pcex.lab)
points(x,det$y,pch=pt,cex=pcex)
x2 <- det2$Age+0.2*det2$big-0.1
points(x2,det2$y,pch=3,cex=pcex)
laby <- c(0,5,10,20,40,60,80)/10
axis(side=2,at=sqrt(laby),cex=pcex,cex.axis=pcex,lab=laby)
mtext(side=3,line=0.1,adj=1,cex=pcex,ttl)
det0 <- subset(det,big==0 & z>=(-0.9)); n0 <- nrow(det0)
det1 <- subset(det,big==1 & z>=(-0.9)); n1 <- nrow(det1)
n2 <- nrow(det2)
lg1 <- paste("Small (n=",n0,")",sep="")
lg1 <- c(lg1,paste("Large (n=",n1,")",sep=""))
lg1 <- c(lg1,paste("Outliers (n=",n2,")",sep=""))
legend("topright",leg=lg1,pch=c(20,1,3),inset=c(0.01,0.01),bg="ivory",cex=1.2)

mod0 <- lm(data=det0,y~Age)
mod1 <- lm(data=det1,y~Age)		# fit models to subgroups
mod2 <- lm(data=det2,y~Age)
abline(mod1$coef,col=2)
abline(mod0$coef,col=2,lty=2)
lg2 <- round(c(mod1$coef[2],mod0$coef[2]),4)
legend("bottomright",leg=lg2,inset=c(0.01,0.07),lty=c(1,2),col=c(2,2),
       bg="ivory",title="Slope",x.intersp=0.2,y.intersp=0.8,cex=1.2)

det3 <- subset(det,z>(-0.9))
mod <- lm(data=det3,y~big+age0+age1)	# fit model to all reduced data table
summary(mod)
rsq <- round(100*summary(mod)$adj.r.squared,1)

#-------------------------------------------------Pause 4

windows(wsize,wsize)		# Figure 5 (corresponding Q-Q plot of studentized residuals)
par(las=1,oma=c(0,0,0,0),mar=c(2.5,2.5,2.5,1),mgp=c(1.1,0.1,0),tcl=0.2)
qqnorm(studres(mod),ylab="",pch=pt,main="")
qqline(studres(mod),col=2)
mtext(side=3,line=0.1,adj=-0.1,"Studentized Residuals")
mtext(side=3,line=0.1,adj=1,ttl)
legend("topleft",leg=paste("r-sq: ",rsq,"%",sep=""),bty="n")

#
# Make new plot without the outliers! det3 = reduced
#
windows(wsize,wsize)		# Figure 4 (Plot sqrt(dph) vs Age for three groups)
par(las=1,oma=c(0,0,0,0),mar=c(2.5,2.5,2.5,1),mgp=c(1.1,0.1,0),tcl=0.2)
x <- det3$Age+0.2*det3$big-0.1
pt <- ifelse(det3$big==0,20,21)
plot(x,det3$y,type="n",xlim=xlm,xlab="Patch Age (Years)",ylab=ylab,yaxt="n",cex=pcex,cex.axis=pcex.lab,cex.lab=pcex.lab)
points(x,det3$y,pch=pt,cex=pcex)
laby <- c(0,5,10,20,40,60,80)/10
axis(side=2,at=sqrt(laby),lab=laby,cex=pcex,cex.axis=pcex)
mtext(side=3,line=0.1,adj=1,ttl,cex=pcex)
det0 <- subset(det3,big==0 & z>=(-0.9)); n0 <- nrow(det0)
det1 <- subset(det3,big==1 & z>=(-0.9)); n1 <- nrow(det1)
lg1 <- paste("Small (n=",n0,")",sep="")
lg1 <- c(lg1,paste("Large (n=",n1,")",sep=""))
legend("topright",leg=lg1,pch=c(20,1,3),inset=c(0.01,0.01),bg="ivory",cex=pcex)
mod0 <- lm(data=det0,y~Age)
mod1 <- lm(data=det1,y~Age)		# fit models to subgroups
abline(mod1$coef,col=2)
abline(mod0$coef,col=2,lty=2)
lg2 <- round(c(mod1$coef[2],mod0$coef[2]),4)
legend("bottomleft",leg=lg2,inset=c(0.01,0.07),lty=c(1,2),col=c(2,2),
       bg="ivory",title="Slope",x.intersp=0.2,y.intersp=0.8,cex=1.1)

mod <- lm(data=det3,y~big+age0+age1)	# fit model to all reduced data table
summary(mod)
rsq <- round(100*summary(mod)$adj.r.squared,1)
rlb <- paste("R Sq.: ", rsq, "%")
text(x=10, y=sqrt(0.6), lab = rlb)
#
#------------------------------------------------end of program

#Repeat above with density/year...
xlm <- c(min(det$Year),max(det$Year))
ylm <- c(0,max(det$ptd))
windows(wsize,wsize)		# Figure 1.5 (Plot tph vs Age)
par(las=1,oma=c(0,0,0,0),mar=c(2.5,2.5,2.5,1),mgp=c(1.1,0.1,0),tcl=0.2)
par(cex.axis=pcex, cex.lab=pcex, cex.main=pcex)
plot(det$Year,det$y,xlim=xlm,xlab=xlab,yaxt="n",ylab=ylab,cex=pcex,cex.axis=pcex.lab,cex.lab=pcex.lab)
laby <- c(0,5,10,20,40,60,80)/10
axis(side=2,at=sqrt(laby),lab=laby,cex=pcex)
mtext(side=3,line=0.1,adj=1,cex=pcex,ttl)
n <- nrow(det)
legend("topright",leg=paste("n:",n),bty="n",cex=pcex.lab)
mod <- lm(data=det,y~Year)
summary(mod) # P-value = 0.5784, Adj. r-squared = -0.006871 - no linear pattern 
abline(mod$coef,col=2)
#
#------------------------------------------------
# Repeat with big
#
windows(wsize,wsize)		# Figure 2 (Plot sqrt(dph) vs Age)
par(las=1,oma=c(0,0,0,0),mar=c(2.5,3.0,2.5,0.5),mgp=c(1.1,0.1,0),tcl=0.2)
par(cex.axis=pcex, cex.lab=pcex, cex.main=pcex)
x <- det$Year+0.2*det$big-0.1
xlm <- c(min(det$Year),max(det$Year))
pt <- ifelse(det$big==0,20,21)
plot(x,det$y,type="n",xlim=xlm,xlab=xlab,ylab=ylab,yaxt="n",cex=pcex,cex.axis=pcex.lab,cex.lab=pcex.lab)
points(x,det$y,pch=pt,cex=pcex)

laby <- c(0,5,10,20,40,60,80)/10
axis(side=2,at=sqrt(laby),lab=laby,cex=pcex)
mtext(side=3,line=0.1,adj=1,cex=pcex,ttl)

det0b <- subset(det,big==0); n0 <- nrow(det0b)
det1b <- subset(det,big==1); n1 <- nrow(det1b)
lg1 <- paste("Small (n=",n0,")",sep="")
lg1 <- c(lg1,paste("Large (n=",n1,")",sep=""))
legend("topright",leg=lg1,pch=c(20,1),inset=c(0.01,0.01),bg="ivory",cex=pcex)

mod1b <- lm(data=det1b,y~Year)		# fit models to subgroups
summary(mod1b)
mod0b <- lm(data=det0b,y~Year)
summary(mod0b)

#abline(mod$coef,col=1)
abline(mod1b$coef,col=2)
abline(mod0b$coef,col=2,lty=2)

lg2 <- round(c(mod1b$coef[2],mod0b$coef[2]),4)
legend("bottomright",leg=lg2,inset=c(0.01,0.08),lty=c(1,2),col=2,
       bg="ivory",title="Slope",x.intersp=0.2,y.intersp=0.9,cex=1.1)

#------------------------------------------------







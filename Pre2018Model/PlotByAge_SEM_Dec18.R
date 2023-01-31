# ptdByAge_/awn.R
# Oct 15, 2018
# Authors - John Stanek, Shannon McNeil Southern Sierra Research Station

setwd('C:\\SSRS\\YBCU\\LCR\\2018_YBCU_LCR')

library(lme4)       # for glmer
library(lmerTest)   # for p value
library(dplyr)      # for cleaner written code


#  library(AED)  #This package is no longer available by Highstat Statistics, 
#     It is a set of functions used in data visualization (such as corvif() below) and exploration
# I copied the AED functions into a new file and use it to call the functions
if(!exists("AED_functions.R", mode="function")) source("AED_functions.R")


### to plot figures
library(jtools)    
library(ggplot2)
library(sjstats)
library(sjPlot)

# #SEM updated data set
 det <- data.frame(read.table("detsByAge.txt",header=TRUE,as.is=TRUE))
 det <- subset(det, (svdets>-1)) # Remove any missing data
 
 #These 3 sites have too few observations to use as random effects, removed from the analysis
 det <- filter(det, det$Site != 'CNHFS' &
                 det$Site != 'CVCA07' &
                 det$Site != 'CVCA08')
 
 
# #original dataset supplied by Shannon McNeil, SSRS
# det <- data.frame(read.table("tphByAge2.txt",header=TRUE,as.is=TRUE))
# det <- subset(det, (ptd>-1) & (is.na(note) | (note != "presv")))
# 

# Data Description:
# Only Area and Age are used in this model.
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
# 
# det$youngest <- factor(det$youngest)
# det$COBPRB <- det$CO + det$PR

# #######    Data Exploration     ##########
str(det)
summary(det)
table(det$Site)

dotchart(det$ha, main = "area (ha)")
dotchart(det$Age, main = "age")
dotchart(det$svdets, main = "survey detections")

hist(det$ha, labels = TRUE, main = "area (ha)")
hist(det$Age, labels = TRUE,main = "age")

boxplot(det$svdets)
boxplot(det$ha)

table(det$Site)
# 
# head(det)

#data visualization plot - AED functions used here
Z <- cbind(det$svdets, # potential response variable
           det$ha, # potential predictor variables
           det$Age) # potential predictor variables

colnames(Z)<-c("SvDets", "ha", "Age")  # no spaces in col names
pairs(Z, lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist);    #pairwise scatterplots and correlation coefficients to examine for colinearity

#check for correlation and VIF - VIF less than 3 is good.
corvif(Z[,c(-1,-8)])   # from AED script file, gets correlation coefficients between pairs and VIF for each


# #################  Territory count model with poisson distribution   ###########

# R4  <- lmer(svdets ~ ha*Age+ (1|Site)+(1|Year), data = det, family = poisson)
# summary(R4)

#tansformed the response variable, no longer need a poisson distribution for count data
det$y<-sqrt(det$svdets)
R4  <- lmer(y ~ ha*Age+ (1|Site)+(1|Year) , data = det)
summary(R4)


############   Model Diagnositcs    ############

#Residuals look pretty good
qqnorm(residuals(R4))

#Section 7.4.4 from Zuur 2013 "A Beginner's Guide to GLM and GLMM in R"
#overdispersion - variation is greater than that predicted by model

#Overdispersion result from model,  0.5740415 less than 1, not overdispersed
E1 <- residuals(R4)
p1 <- length(fixef(R4)) + 1
Overdisp1 <- sum(E1^2) / (nrow(det) - p1)
Overdisp1


#Checking for Deviance Residual GOF, a check for overdispersion
phi <- sum(resid(R4, type = "pearson")^2) / df.residual(R4)
phi
#  0.5862551 -  data is not overdispersed


#Residual diagnostic plots
op <- par(mfrow=c(2,3),mar=c(5,4,1,2))
E<-resid(R4)
fit <-fitted(R4)
hist(E,xlab="Residuals",main="")   #histogram to check for normality
qqnorm(E)                           #QQ plot to check for normality
plot(fit,E, main="Residuals vs. Fitted Values", ylab="Residual", xlab="Fitted values")    #fitted vs residuals to check homogenaity/heteroskedacity - look for patterns
plot(det$ha,E,xlab="ha",ylab="Residuals")                  #residuals vs each explanatory variable to check for independence
plot(det$Age,E,xlab="age",ylab="Residuals")                #residuals vs each explanatory variable to check for independence
plot(det$nAge,E,xlab="n age",ylab="Residuals")             #residuals vs each explanatory variable to check for independence
par(op)



# Looking at the various intercept values for random effect variables
# coef(R4)
# coefs <-coef(R4)$Site[,"(Intercept)"]
# plot(coefs)
# boxplot(coefs)
# 
# coef(R4)
# coefs <-coef(R4)$Year[,"(Intercept)"]
# plot(coefs)
# boxplot(coefs)

######    Plot of the model    ############

InteractionPlot = plot_model(R4, type = "pred",  colors = "bw", terms = c("Age", "ha [20,50,80]"), title = "", axis.title = "Square Root of Detections" )
InteractionPlot  + theme_bw(base_size = 20) + scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12)) 



# ptdByAge_JRS.R
# Oct 15, 2018
# Author - John Stanek, Southern Sierra Research Station

setwd('C:\\SSRS\\YBCU\\LCR\\2018_YBCU_LCR')

library(lme4)       # for glmer
library(car)
library(MASS)
library(MuMIn)   # for AICc function, R2 (r.squaredGLMM(R4))

#  library(AED)  #This package is no longer available by Highstat Statistics, 
#     It is a set of functions used in data visualization (such as corvif() below) and exploration
# I copied the AED functions into a new file and use it to call the functions
if(!exists("AED_functions.R", mode="function")) source("AED_functions.R")

#For Model Averaging - not used
# library(AICcmodavg)

### to plot figures
library(jtools)
library(ggplot2)
library(sjstats)
library(sjPlot)

#original dataset supplied by Shannon McNeil, SSRS
det <- data.frame(read.table("tphByAge2.txt",header=TRUE,as.is=TRUE))
det <- subset(det, (ptd>-1) & (is.na(note) | (note != "presv")))

# Data Description:
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

det$youngest <- factor(det$youngest)
det$COBPRB <- det$CO + det$PR

#######    Data Exploration     ##########
str(det)
summary(det)
table(det$Site)

dotchart(det$ha, main = "area (ha)")
dotchart(det$Age, main = "age")
dotchart(det$svdets, main = "survey detections")
dotchart(det$cha2, main = "contiguous area size")
dotchart(det$COBPRB, main = "counts of COB+PRB")


hist(det$ha, labels = TRUE, main = "area (ha)")
hist(det$Age, labels = TRUE,main = "age")
hist(det$cha2, labels = TRUE,main = "Continuous Ha")
hist(det$nAges, labels = TRUE,main = "nAges")
hist(det$COBPRB, labels = TRUE,main = "COBPRB")
boxplot(det$COBPRB)
table(det$Site)

head(det)

#data visualization plot - AED functions used here
Z <- cbind(det$COBPRB, # potential response variable
           det$Year, # potential predictor variables
           det$ha, # potential predictor variables
           det$Age, # potential predictor variables
           det$cha2, # potential predictor variables
           det$nAges,# potential predictor variables
           det$youngest)# potential predictor variables
colnames(Z)<-c("COBPRB", "Year", "ha", "Age", "cha2", "nAges", "youngest")  # no spaces in col names
pairs(Z, lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist);    #pairwise scatterplots and correlation coefficients to examine for colinearity

#check for correlation and VIF - A VIF less than 3 is good. 
corvif(Z[,c(-1,-8)])   # from AED script file, gets correlation coefficients between pairs and VIF for each 



#################  Territory count models with poisson distribution   ###########
## All combinations of global model (R1) run to give equal credence to each covariate

## Rerescaling of cha2 variable leads to exact same results and removes warnings,
#     but rescaling adds difficulty in interpreting results, analyzed as w/o rescaling as the results are the same either way
det$cha2 <- scale(det$cha2)
# det$cha2<- (cha2 - mean(cha2)) / sd(cha2)  ## does the same as scale() function above


R1  <- glmer(COBPRB ~ ha + Age + cha2 + nAges + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R2  <- glmer(COBPRB ~ ha + Age + cha2 + nAges +(1|Site)+(1|Year), data = det, family = poisson)
R3  <- glmer(COBPRB ~ ha + Age + cha2 +(1|Site)+(1|Year), data = det, family = poisson)
R4  <- glmer(COBPRB ~ ha + Age +(1|Site)+(1|Year), data = det, family = poisson)
R5  <- glmer(COBPRB ~ ha + cha2 + nAges + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R6  <- glmer(COBPRB ~ ha + cha2 + nAges +(1|Site)+(1|Year), data = det, family = poisson)
R7  <- glmer(COBPRB ~ ha + cha2 + youngest+(1|Site)+(1|Year), data = det, family = poisson)
R8  <- glmer(COBPRB ~ ha + cha2 +(1|Site)+(1|Year), data = det, family = poisson)
R9  <- glmer(COBPRB ~ ha + Age + nAges + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R10  <- glmer(COBPRB ~ ha + Age + nAges +(1|Site) +(1|Year), data = det, family = poisson)
R11  <- glmer(COBPRB ~ ha + Age + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R12  <- glmer(COBPRB ~ ha + nAges + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R13  <- glmer(COBPRB ~ ha + nAges +(1|Site)+(1|Year), data = det, family = poisson)
R14  <- glmer(COBPRB ~ ha + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R15  <- glmer(COBPRB ~ ha +(1|Site)+(1|Year), data = det, family = poisson)

R16  <- glmer(COBPRB ~ Age + cha2 + nAges + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R17  <- glmer(COBPRB ~ Age + cha2 + nAges +(1|Site)+(1|Year), data = det, family = poisson)

R18  <- glmer(COBPRB ~ ha + Age + cha2 + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R19  <- glmer(COBPRB ~ Age + cha2 + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R20  <- glmer(COBPRB ~ Age + nAges + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R21  <- glmer(COBPRB ~ Age + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R22  <- glmer(COBPRB ~ Age + cha2 +(1|Site)+(1|Year), data = det, family = poisson)
R23  <- glmer(COBPRB ~ Age + nAges +(1|Site)+(1|Year), data = det, family = poisson)
R24  <- glmer(COBPRB ~ Age +(1|Site)+(1|Year), data = det, family = poisson)

R25 <- glmer(COBPRB ~ cha2 + nAges + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R26 <- glmer(COBPRB ~ cha2 + nAges +(1|Site)+(1|Year), data = det, family = poisson)
R27 <- glmer(COBPRB ~ cha2 + youngest +(1|Site)+(1|Year), data = det, family = poisson)
R28 <- glmer(COBPRB ~ cha2 +(1|Site)+(1|Year), data = det, family = poisson)

R29 <- glmer(COBPRB ~ nAges + youngest +(1|Site)+(1|Year), data = det, family = poisson)

R30 <- glmer(COBPRB ~ nAges +(1|Site)+(1|Year), data = det, family = poisson)
R31 <- glmer(COBPRB ~ youngest +(1|Site)+(1|Year), data = det, family = poisson)


#############
AICs <- AICc(R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29, R30, R31)
models <- c(R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29, R30, R31)
logs <- lapply(models, function(i) logLik(i)[1])
like <- data.frame(Reduce(rbind, logs))

# Display a table of Model#, AIC, delta AIC, & Akaike wt
MyDf <- AICs[,1]
AICsNum <- AICs[,2]
minAW <- min(AICsNum)
Delta <- AICsNum-minAW
RL <- exp(-0.5 * Delta)
wi <- RL / sum(RL)
Z <- data.frame(MyDf, AICsNum, Delta, wi, like)
Z <- round(Z, digits = 3)
colnames(Z)<- c("Df", "AICc", " Delta_AICc", "Akaike_weights", "log_likelihood")
Z <- Z[order(Z[,3]),]
Z


#Top models <2 delta AIC
summary(R10)
summary(R2)
# Delta AICc close to 2
summary(R1)
summary(R9)


############   Model Diagnositcs    ############

# icc(R10)
# icc_Site = 0.2576
# icc_Year = 0.1034
# icc_sum = icc_Site + icc_Year
# AveSamplesPerSite = 6.375
# ApproxDesignEffect = 1 +(AveSamplesPerSite-1) * icc_sum
# approxESS = 102/ApproxDesignEffect
# approxESS

#Residuals look pretty good
qqnorm(residuals(R10))

#Section 7.4.4 from Zuur 2013 "A Beginner's Guide to GLM and GLMM in R"
#overdispersion - variation is greater than that predicted by model

#Overdispersion result from model 1 (Global Model): 0.8586421 less than 1, so no overdispersion
E1 <- residuals(R1)
p1 <- length(fixef(R1)) + 1
Overdisp1 <- sum(E1^2) / (nrow(det) - p1)
Overdisp1


#Overdispersion result from model 10: 0.8679573 less than 1, so no overdispersion
E1 <- residuals(R10)
p1 <- length(fixef(R10)) + 1
Overdisp1 <- sum(E1^2) / (nrow(det) - p1)
Overdisp1


#Checking for Deviance Residual GOF, a check for overdispersion
phi <- sum(resid(R1, type = "pearson")^2) / df.residual(R1)
phi
# 0.7196588 -  data is not overdispersed


#Residual diagnostic plots
op <- par(mfrow=c(2,3),mar=c(5,4,1,2))
E<-resid(R17)
fit <-fitted(R17)
hist(E,xlab="Residuals",main="")   #histogram to check for normality
qqnorm(E)                           #QQ plot to check for normality
plot(fit,E, main="Residuals vs. Fitted Values", ylab="Residual", xlab="Fitted values")    #fitted vs residuals to check homogenaity/heteroskedacity - look for patterns
plot(det$ha,E,xlab="ha",ylab="Residuals")                  #residuals vs each explanatory variable to check for independence
plot(det$Age,E,xlab="age",ylab="Residuals")                #residuals vs each explanatory variable to check for independence
plot(det$nAge,E,xlab="n age",ylab="Residuals")             #residuals vs each explanatory variable to check for independence
par(op)


# Looking at the various intercept values for random effect variables
coef(R10)
coefs <-coef(R10)$Site[,"(Intercept)"]
plot(coefs)
boxplot(coefs)

coef(R10)
coefs <-coef(R10)$Year[,"(Intercept)"]
plot(coefs)
boxplot(coefs)


###########  Model Averaging      ##### 
# not used- we selected the most parsimonious top model
# # model.set <- dredge(R1)
# # top.models <- get.models(model.set, subset = delta< 2)
# # # top.models <- get.models(model.set,cumsum(weight)< 0.95)#
# # aveModel <- model.avg(top.models, methods = "NA")
# # summary(aveModel)
# # confint(aveModel)
# 
# # model averaging
# M1=model.avg(R10, R2, R9) # get averaged coefficients from delta AIC less than 2
# summary(M1)
# confint.default(M1)

######## Plot the results #############
# Model 10: glmer(COBPRB ~ ha + Age + nAges +(1|Site), data = det, family = poisson)
summ(R10)

#jTools to plot the predicted values, and ggplot2 themes to change the non-data formating in plot
SizePlot = effect_plot(R10, pred = ha, interval = TRUE, x.label = "Site Size (ha)", y.label = "Territory Abundance")
SizePlot + theme_bw(base_size = 22)


AgePlot = effect_plot(R10, pred = Age, interval = TRUE, x.label = "Age (Years)", y.label = "Territory Abundance")
# AgePlot + theme_bw(base_size = 22)
AgePlot  + theme_bw(base_size = 22) + scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12))

nAgesPlot = effect_plot(R10, pred = nAges, interval = TRUE, x.label = "Adjacent Site Age Variation", y.label = "Territory Abundance")
nAgesPlot + theme_bw(base_size = 22)

#manually set the plot window to 500 x 500 pixels, pasted plot into clipboard, then saved as pic for report.


#JRS new plot Dec 17, 2018
plot_model(R4, type = "pred", terms = c("Age", "ha"))
rm(list = ls())
library(data.table)
library(stringr)
library(limma)
library(ggbiplot)
library(dplyr)
library(sm)
# ResponderSubtype
clinical.data <- read.table(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\ClinicalDataFinal.txt', sep ="\t", header=TRUE)

clinical.data2 <- clinical.data[complete.cases(clinical.data$DeltaImmune),]
sm.density.compare(clinical.data2$DeltaImmune, clinical.data2$ResponderType, model="equal", nboot= 500, ngrid= 100, xlab = 'Delta Immune')
title(main="Distribution by sample subtype")
legend("topright", levels(clinical.data$ResponderType), fill=2+(0:nlevels(clinical.data$ResponderType)))

summary(clinical.data$DeltaImmune[clinical.data$ResponderType == 'U'])
summary(clinical.data$DeltaImmune[clinical.data$ResponderType == 'D'])

aov1 = aov(clinical.data$NES ~ clinical.data$RecurrentSubtype)
summary(aov1)

clinical.data2 <- clinical.data[complete.cases(clinical.data$NES),]
sm.density.compare(clinical.data2$NES, clinical.data2$RecurrentSubtype, model="equal", nboot= 500, ngrid= 100, xlab = 'NES')
title(main="Distribution by recurrent subtype")
legend("topleft", levels(clinical.data$RecurrentSubtype), fill=2+(0:nlevels(clinical.data$RecurrentSubtype)))





sm.density.compare(clinical.data$PrimaryPurity, clinical.data$ResponderType, model="equal", nboot= 500, ngrid= 100, xlab = 'Recurrent Purity')
title(main="Distribution by sample subtype")
legend("topright", levels(clinical.data$ResponderType), fill=2+(0:nlevels(clinical.data$ResponderType)))

#t-tests instead of ANOVA
ttest1 <- t.test(clinical.data$NES ~ clinical.data$ResponderType)
ttest1
ttest2 <- t.test(clinical.data$MGMT_ExpP ~ clinical.data$ResponderType)
ttest2
ttest3 <- t.test(clinical.data$Age ~ clinical.data$ResponderType)
ttest3
ttest4 <- t.test(clinical.data$PFS_m ~ clinical.data$ResponderType)
ttest4
ttest5 <- t.test(clinical.data$OS_m ~ clinical.data$ResponderType)
ttest5
ttest6 <- t.test(clinical.data$PrimaryPurity ~ clinical.data$ResponderType)
ttest6
ttest7 <- t.test(clinical.data$RecurrentPurity ~ clinical.data$ResponderType)
ttest7
ttest8 <- t.test(clinical.data$DeltaPurity ~ clinical.data$ResponderType)
ttest8
ttest9 <- t.test(clinical.data$PrimaryImmuneScore ~ clinical.data$ResponderType)
ttest9
ttest10 <- t.test(clinical.data$RecurrentImmuneScore ~ clinical.data$ResponderType)
ttest10
ttest11 <- t.test(clinical.data$DeltaImmune ~ clinical.data$ResponderType)
ttest11

#NES
aov1 = aov(clinical.data$NES ~ clinical.data$ResponderType)
summary(aov1)
#MGMT-Exp
aov2 = aov(clinical.data$MGMT_ExpP ~ clinical.data$ResponderType)
summary(aov2)
#Gender
clinical.data2 <- clinical.data[complete.cases(clinical.data$Gender),]
chsq1 <- chisq.test(clinical.data2$Gender, clinical.data2$ResponderType)
chsq1$p.value
#Age
aov3 = aov(clinical.data$Age ~ clinical.data$ResponderType)
summary(aov3)
#LocationPrimary
clinical.data2 <- clinical.data[complete.cases(clinical.data$LocationPrimary),]
chsq2 <- chisq.test(clinical.data2$LocationPrimary, clinical.data2$ResponderType)
chsq2$p.value
#Location Recurrent
clinical.data2 <- clinical.data[complete.cases(clinical.data$LocationRecurrence),]
chsq3 <- chisq.test(clinical.data2$LocationRecurrence, clinical.data2$ResponderType)
chsq3$p.value
#PFS
aov4 = aov(clinical.data$PFS_m ~ clinical.data$ResponderType)
summary(aov4)
#OS
aov5 = aov(clinical.data$OS_m ~ clinical.data$ResponderType)
summary(aov5)
#Primary Purity
aov6 = aov(clinical.data$PrimaryPurity ~ clinical.data$ResponderType)
summary(aov6)
#Recurrent Purity
aov7 = aov(clinical.data$RecurrentPurity ~ clinical.data$ResponderType)
summary(aov7)
#Delta Purity
aov8 = aov(clinical.data$DeltaPurity ~ clinical.data$ResponderType)
summary(aov8)
#Primary Immune Score
aov9 = aov(clinical.data$PrimaryImmuneScore ~ clinical.data$ResponderType)
summary(aov9)
#Recurrent Immune Score
aov10 = aov(clinical.data$RecurrentImmuneScore ~ clinical.data$ResponderType)
summary(aov10)
#Delta Immune Score
aov11 = aov(clinical.data$DeltaImmune ~ clinical.data$ResponderType)
summary(aov11)

#####################
###Other Variables###
#####################

#MGMT-PS
clinical.data2 <- clinical.data[complete.cases(clinical.data$MGMT_PS),]
chsq4 <- chisq.test(clinical.data2$MGMT_PS, clinical.data2$ResponderType)
chsq4$p.value
#Cohort
clinical.data2 <- clinical.data[complete.cases(clinical.data$Cohort),]
chsq5 <- chisq.test(clinical.data2$Cohort, clinical.data2$ResponderType)
chsq5$p.value
#IDH
clinical.data2 <- clinical.data[complete.cases(clinical.data$IDH),]
clinical.data2$IDH <- as.factor(clinical.data2$IDH)
chsq6 <- chisq.test(clinical.data2$IDH, clinical.data2$ResponderType)
chsq6$p.value
#Primary Subtype
clinical.data2 <- clinical.data[complete.cases(clinical.data$PrimarySubtype),]
chsq7 <- chisq.test(clinical.data2$PrimarySubtype, clinical.data2$ResponderType)
chsq7$p.value
#Primary Unanimous
clinical.data2 <- clinical.data[complete.cases(clinical.data$PrimaryUnanimous),]
chsq8 <- chisq.test(clinical.data2$PrimaryUnanimous, clinical.data2$ResponderType)
chsq8$p.value
#Recurrent Subtype
clinical.data2 <- clinical.data[complete.cases(clinical.data$RecurrentSubtype),]
chsq9 <- chisq.test(clinical.data2$RecurrentSubtype, clinical.data2$ResponderType)
chsq9$p.value
#Recurrent Unanimous
clinical.data2 <- clinical.data[complete.cases(clinical.data$RecurrentUnanimous),]
chsq10 <- chisq.test(clinical.data2$RecurrentUnanimous, clinical.data2$ResponderType)
chsq10$p.value
#Subtype Stable
clinical.data2 <- clinical.data[complete.cases(clinical.data$SubtypeStable),]
chsq11 <- chisq.test(clinical.data2$SubtypeStable, clinical.data2$ResponderType)
chsq11$p.value

###############################
###Progression-free Survival###
###############################
cor.test(clinical.data$PFS_m, clinical.data$MGMT_ExpP)
cor.test(clinical.data$PFS_m, clinical.data$Age)
cor.test(clinical.data$PFS_m, clinical.data$NES)
cor.test(clinical.data$PFS_m, clinical.data$PrimaryPurity)
cor.test(clinical.data$PFS_m, clinical.data$RecurrentPurity)
cor.test(clinical.data$PFS_m, clinical.data$DeltaPurity)
cor.test(clinical.data$PFS_m, clinical.data$PrimaryImmuneScore)
cor.test(clinical.data$PFS_m, clinical.data$RecurrentImmuneScore)
cor.test(clinical.data$PFS_m, clinical.data$DeltaImmune)

clinical.data$IDH <- as.factor(clinical.data$IDH)
lm <- lm(formula = PFS_m ~ IDH, data = clinical.data)
summary(lm)

lm <- lm(formula = PFS_m ~ Gender, data = clinical.data)
summary(lm)

lm <- lm(formula = PFS_m ~ LocationPrimary, data = clinical.data)
summary(lm)

lm <- lm(formula = PFS_m ~ LocationRecurrence, data = clinical.data)
summary(lm)

lm <- lm(formula = PFS_m ~ PrimarySubtype, data = clinical.data)
summary(lm)

lm <- lm(formula = PFS_m ~ RecurrentSubtype, data = clinical.data)
summary(lm)

lm <- lm(formula = PFS_m ~ Cohort, data = clinical.data)
summary(lm)

lm <- lm(formula = PFS_m ~ PrimaryUnanimous, data = clinical.data)
summary(lm)

lm <- lm(formula = PFS_m ~ RecurrentUnanimous, data = clinical.data)
summary(lm)

lm <- lm(formula = PFS_m ~ SubtypeStable, data = clinical.data)
summary(lm)

###############################
###Overall Survival###
###############################
cor.test(clinical.data$OS_m, clinical.data$MGMT_ExpP)
cor.test(clinical.data$OS_m, clinical.data$Age)
cor.test(clinical.data$OS_m, clinical.data$NES)
cor.test(clinical.data$OS_m, clinical.data$PrimaryPurity)
cor.test(clinical.data$OS_m, clinical.data$RecurrentPurity)
cor.test(clinical.data$OS_m, clinical.data$DeltaPurity)
cor.test(clinical.data$OS_m, clinical.data$PrimaryImmuneScore)
cor.test(clinical.data$OS_m, clinical.data$RecurrentImmuneScore)
cor.test(clinical.data$OS_m, clinical.data$DeltaImmune)

clinical.data$IDH <- as.factor(clinical.data$IDH)
lm <- lm(formula = OS_m ~ IDH, data = clinical.data)
summary(lm)

lm <- lm(formula = OS_m ~ Gender, data = clinical.data)
summary(lm)

lm <- lm(formula = OS_m ~ LocationPrimary, data = clinical.data)
summary(lm)

lm <- lm(formula = OS_m ~ LocationRecurrence, data = clinical.data)
summary(lm)

lm <- lm(formula = OS_m ~ PrimarySubtype, data = clinical.data)
summary(lm)

lm <- lm(formula = OS_m ~ RecurrentSubtype, data = clinical.data)
summary(lm)

lm <- lm(formula = OS_m ~ Cohort, data = clinical.data)
summary(lm)

lm <- lm(formula = OS_m ~ PrimaryUnanimous, data = clinical.data)
summary(lm)

lm <- lm(formula = OS_m ~ RecurrentUnanimous, data = clinical.data)
summary(lm)

lm <- lm(formula = OS_m ~ SubtypeStable, data = clinical.data)
summary(lm)

#NES#
cor.test(clinical.data$NES, clinical.data$MGMT_ExpP)
cor.test(clinical.data$NES, clinical.data$Age)
cor.test(clinical.data$NES, clinical.data$OS_m)
cor.test(clinical.data$NES, clinical.data$PFS_m)
cor.test(clinical.data$NES, clinical.data$PrimaryPurity)
cor.test(clinical.data$NES, clinical.data$RecurrentPurity)
cor.test(clinical.data$NES, clinical.data$DeltaPurity)
cor.test(clinical.data$NES, clinical.data$PrimaryImmuneScore)
cor.test(clinical.data$NES, clinical.data$RecurrentImmuneScore)
cor.test(clinical.data$NES, clinical.data$DeltaImmune)

clinical.data$IDH <- as.factor(clinical.data$IDH)
lm <- lm(formula = NES ~ IDH, data = clinical.data)
summary(lm)

lm <- lm(formula = NES ~ Gender, data = clinical.data)
summary(lm)

lm <- lm(formula = NES ~ LocationPrimary, data = clinical.data)
summary(lm)

lm <- lm(formula = NES ~ LocationRecurrence, data = clinical.data)
summary(lm)

lm <- lm(formula = NES ~ PrimarySubtype, data = clinical.data)
summary(lm)

lm <- lm(formula = NES ~ RecurrentSubtype, data = clinical.data)
summary(lm)

lm <- lm(formula = NES ~ PrimaryUnanimous, data = clinical.data)
summary(lm)

lm <- lm(formula = NES ~ RecurrentUnanimous, data = clinical.data)
summary(lm)

lm <- lm(formula = NES ~ SubtypeStable, data = clinical.data)
summary(lm)

lm <- lm(formula = NES ~ Cohort, data = clinical.data)
summary(lm)


plot(NES ~ RecurrentSubtype, data = clinical.data)
abline(lm(formula = RecurrentPurity ~ NES, data = clinical.data))


clinical.dataLowNES <- clinical.data[clinical.data$NES < 0,]
clinical.dataLowNES <- clinical.dataLowNES[complete.cases(clinical.dataLowNES$NES),]

clinical.dataHighNES <- clinical.data[clinical.data$NES > 0,]
clinical.dataHighNES <- clinical.dataHighNES[complete.cases(clinical.dataHighNES$NES),]

sm.density.compare(clinical.dataLowNES$PrimaryPurity, clinical.dataLowNES$RecurrentPurity, model="equal", nboot= 500, ngrid= 100, xlab = 'Delta Immune')
title(main="Distribution by sample subtype")
legend("topright", levels(clinical.data$ResponderType), fill=2+(0:nlevels(clinical.data$ResponderType)))

plot(density(clinical.dataLowNES$PrimaryPurity), xlab = 'Purity', main = 'Purity in primary/recurrent samples with Low NES')
lines(density(clinical.dataLowNES$RecurrentPurity), col = 'red')
legend("topleft", legend=c('PrimaryPurity', 'RecurrentPurity'), fill=1:2)

plot(density(clinical.dataHighNES$PrimaryPurity), xlab = 'Purity', main = 'Purity in primary/recurrent samples with High NES', ylim = range(0:3))
lines(density(clinical.dataHighNES$RecurrentPurity), col = 'red')
legend("topleft", legend=c('PrimaryPurity', 'RecurrentPurity'), fill=1:2)


summary(clinical.dataLowNES$PrimaryPurity)
summary(clinical.dataLowNES$RecurrentPurity)

summary(clinical.dataHighNES$PrimaryPurity)
summary(clinical.dataHighNES$RecurrentPurity)


clinical.data$ResponderType[clinical.data$RecurrentSubtype == 'Proneural']
clinical.data$ResponderType[clinical.data$RecurrentSubtype == 'Classical']
clinical.data$ResponderType[clinical.data$RecurrentSubtype == 'Mesenchymal']




lm <- lm(formula = NES ~ PrimarySubtype, data = clinical.data)
summary(lm)
plot(NES ~ PrimarySubtype, data = clinical.data)
abline(lm)


summary(clinical.data$DeltaPurity[clinical.data$ResponderType == 'U'])
summary(clinical.data$DeltaPurity[clinical.data$ResponderType == 'D'])

t.test(clinical.data$DeltaPurity[clinical.data$ResponderType == 'U'], clinical.data$DeltaPurity[clinical.data$ResponderType == 'D'])

library(vars)
causalityvars <- subset(clinical.data, select=c(14, 22, 23, 26, 18))


causalityvars <- subset(clinical.data, select=c(22, 23, 26))
causality(varcause, cause = 'RecurrentPurity', boot=TRUE, boot.runs=1000)

summary(causalityvars$ResponderType)
causalityvars$RecurrentSubtype


lm <- lm(formula = RecurrentPurity ~ ResponderType + DeltaPurity + DeltaImmune, data = clinical.data)
summary(lm)
plot(RecurrentPurity ~ ResponderType + DeltaPurity + DeltaImmune, data = clinical.data)

par(mfrow=c(2,2))
plot(clinical.data$PrimarySubtype[clinical.data$ResponderType == 'U'],col="darkgreen", main = 'Primary U samples') 
plot(clinical.data$RecurrentSubtype[clinical.data$ResponderType == 'U'],col="red", main = 'Recurrent U samples')
plot(clinical.data$PrimarySubtype[clinical.data$ResponderType == 'D'],col="darkgreen", main = 'Primary D samples') 
plot(clinical.data$RecurrentSubtype[clinical.data$ResponderType == 'D'],col="red", main = 'Recurrent D samples')


clinical.dataU <- clinical.data[clinical.data$ResponderType == 'U',]
clinical.dataD <- clinical.data[clinical.data$ResponderType == 'D',]
plot1 <- ggplot(clinical.dataU, aes(x="", fill=clinical.dataU$PrimarySubtype)) + ggtitle('Primary U samples') + geom_bar(width = 1) + coord_polar("y") + theme_void() + guides(fill=guide_legend(title="Legend"))
plot2 <- ggplot(clinical.dataU, aes(x="", fill=clinical.dataU$RecurrentSubtype)) + ggtitle('Recurrent U samples') + geom_bar(width = 1) + coord_polar("y") + theme_void() + guides(fill=guide_legend(title="Legend"))
plot3 <- ggplot(clinical.dataD, aes(x="", fill=clinical.dataD$PrimarySubtype)) + ggtitle('Primary D samples') + geom_bar(width = 1) + coord_polar("y") + theme_void() + guides(fill=guide_legend(title="Legend"))
plot4 <- ggplot(clinical.dataD, aes(x="", fill=clinical.dataD$RecurrentSubtype)) + ggtitle('Recurrent D samples') + geom_bar(width = 1) + coord_polar("y") + theme_void() + guides(fill=guide_legend(title="Legend"))
require(gridExtra)
grid.arrange(plot1, plot2, plot3, plot4, ncol=2)


#NES
aov = aov(clinical.data$RecurrentPurity ~ clinical.data$RecurrentSubtype)
summary(aov)

par(mfrow=c(1,2))
plot(clinical.data$PrimaryPurity ~ clinical.data$PrimarySubtype, data = clinical.data, ylab = 'Primary Purity', xlab = 'Primary Subtype')
plot(clinical.data$RecurrentPurity ~ clinical.data$RecurrentSubtype, data = clinical.data, ylab = 'Recurrent Purity', xlab = 'Recurrent Subtype')


plot(clinical.data$PrimaryPurity ~ clinical.data$RecurrentSubtype, data = clinical.data)
plot(clinical.data$RecurrentPurity ~ clinical.data$PrimarySubtype, data = clinical.data)


t.test(clinical.data$PrimaryPurity[clinical.data$ResponderType == 'D'], clinical.data$RecurrentPurity[clinical.data$ResponderType == 'D'], paired = T)

par(mfrow=c(1,1))
plot(clinical.data$PrimaryPurity[clinical.data$ResponderType == 'U'], clinical.data$RecurrentPurity[clinical.data$ResponderType == 'U'])
cor.test(clinical.data$PrimaryPurity[clinical.data$ResponderType == 'U'], clinical.data$RecurrentPurity[clinical.data$ResponderType == 'U'])







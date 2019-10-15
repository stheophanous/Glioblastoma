rm(list = ls())
library(data.table)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)
# Clinical data
clinical.data <- read.table(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\ClinicalDataFinal.txt', sep ="\t", header=TRUE)
summary(clinical.data)
clinical.data$IDH <- as.factor(clinical.data$IDH)
clinical.data$event <- NA
for (i in 1:nrow(clinical.data)){
  if (clinical.data$Status[i] == 'DECEASED'){
    clinical.data$event[i] <- 1
  } else {
    clinical.data$event[i] <- 0
  }
}
clinical.data$ResponderType <- factor(clinical.data$ResponderType, levels = c("U", "D"))

# Fit survival data using the Kaplan-Meier method for OS
surv_objectOS <- Surv(time = clinical.data$OS_m, event = clinical.data$event)
surv_objectOS
fit1 <- survfit(surv_objectOS ~ ResponderType, data = clinical.data)
summary(fit1)
surv_median(fit1)
ggsurvplot(fit1, data = clinical.data, pval = TRUE)
autoplot(fit1, main = 'Overall Survival')
ggsurvplot(
  fit1,                     # survfit object with calculated statistics.
  pval = T,             # show p-value of log-rank test.
  xlab = "Time in Months",   # customize X axis label
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = F,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = T,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  palette = c("#E7B800", "#2E9FDF") # custom color palettes.
)

surv_diff <- survdiff(Surv(OS_m) ~ ResponderType, data = clinical.data)
surv_diff
surv_diff <- survdiff(Surv(PFS_m) ~ ResponderType, data = clinical.data)
surv_diff

# Fit survival data using the Kaplan-Meier method for PFS
surv_objectPFS <- Surv(time = clinical.data$PFS_m, event = clinical.data$event)
surv_objectPFS
fit2 <- survfit(surv_objectPFS ~ ResponderType, data = clinical.data)
summary(fit2)
ggsurvplot(fit2, data = clinical.data, pval = TRUE)

#Multivariate analysis for OS - linear regression
lmfit <- lm(OS_m ~ MGMT_ExpP + Age + IDH, data = clinical.data)
summary(lmfit)

# Univariate Cox regression
# Fit a Cox proportional hazards model for PFS
fit.coxph <- coxph(surv_objectPFS ~ MGMT_ExpP, data = clinical.data)
summary(fit.coxph)
ggforest(fit.coxph, data = clinical.data)

# Fit a Cox proportional hazards model for OS
fit.coxph <- coxph(surv_objectOS ~ Cohort, data = clinical.data)
summary(fit.coxph)
ggforest(fit.coxph, data = clinical.data)


# Multivariate Cox regression
# Fit a Cox proportional hazards model for PFS
fit.coxph <- coxph(surv_objectPFS ~ PrimaryImmuneScore + IDH + MGMT_ExpP + LocationPrimary + LocationRecurrence + Cohort, data = clinical.data)
ggforest(fit.coxph, data = clinical.data)
survcox <- survfit(fit.coxph)
autoplot(survcox)

# Fit a Cox proportional hazards model for OS
fit.coxph <- coxph(surv_objectOS ~ MGMT_ExpP + Age + IDH + LocationPrimary + LocationRecurrence, data = clinical.data)
ggforest(fit.coxph, data = clinical.data)
survcox <- survfit(fit.coxph)
autoplot(survcox)


surv_objectPFS <- Surv(time = clinical.data$PFS_m, event = clinical.data$event)
surv_objectPFS
fit1 <- survfit(surv_objectPFS ~ ResponderType, data = clinical.data)
summary(fit1)
surv_median(fit1)
ggsurvplot(fit1, data = clinical.data, pval = TRUE)
autoplot(fit1, main = 'Progression-Free Survival')
ggsurvplot(
  fit1,                     # survfit object with calculated statistics.
  pval = T,             # show p-value of log-rank test.
  xlab = "Time in Months",   # customize X axis label
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = F,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = T,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  palette = c("#E7B800", "#2E9FDF") # custom color palettes.
)




rm(list = ls())
library(data.table)
library(stringr)
library(limma)
library(ggbiplot)
library(dplyr)
FPKM.merged <- as.data.frame(fread(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\mergedFPKM_ALLsamples.txt', stringsAsFactors = F))
FPKM.merged <- FPKM.merged[,-ncol(FPKM.merged)]
FPKM.merged <- as.data.frame(t(FPKM.merged))
colnames(FPKM.merged) <- as.character(unlist(FPKM.merged[1,]))
FPKM.merged <- FPKM.merged[-1,]
FPKM.merged[] <- lapply(FPKM.merged, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
sapply(FPKM.merged, class)

# ResponderSubtype
clinical.data <- read.table(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\ClinicalDataFinal.txt', sep ="\t", header=TRUE)
responder.subtype <- as.data.frame(clinical.data[,c(1,14)])
colnames(responder.subtype) <- c('Sample', 'Group')
responder.subtype$Group <- as.factor(responder.subtype$Group)

#Remove batch effect
batch <- c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1)

FPKM.merged.fixed <- removeBatchEffect(FPKM.merged, batch)

par(mfrow=c(1,2))
boxplot(as.data.frame(FPKM.merged),main="Original")
boxplot(as.data.frame(FPKM.merged.fixed),main="Batch corrected")


# Dataset
FPKM.merged.fixed <- as.data.frame(t(FPKM.merged.fixed))
FPKM.merged.fixed$ResponderType <- NA
for (i in 1:nrow(FPKM.merged.fixed)){
  FPKM.merged.fixed$ResponderType[i] <- as.character(responder.subtype$Group[rownames(FPKM.merged.fixed)[i] == responder.subtype$Sample])
}
FPKM.merged.fixed$ResponderType <- as.factor(FPKM.merged.fixed$ResponderType)
#write.table(FPKM.merged.fixed, file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\FPKMallsamplesNoBatchEffects')

FPKM.pca <- prcomp(FPKM.merged.fixed[,-ncol(FPKM.merged.fixed)], center = TRUE, scale. = TRUE)
summary(FPKM.pca)
ggbiplot(FPKM.pca, labels=rownames(FPKM.merged.fixed), var.axes = FALSE, groups = FPKM.merged.fixed$ResponderType, ellipse = TRUE, choices = c(1,2))

#=======================================================================================================================================================
DEgenes <- read.table(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\DEgenes66samples.txt')
#DEgenes$x <- gsub("\\..*","",DEgenes$x)
FPKM.merged.DE <- FPKM.merged.fixed[colnames(FPKM.merged.fixed) %in% DEgenes$x]

# Dataset
FPKM.merged.DE$ResponderType <- NA
for (i in 1:nrow(FPKM.merged.DE)){
  FPKM.merged.DE$ResponderType[i] <- as.character(responder.subtype$Group[rownames(FPKM.merged.DE)[i] == responder.subtype$Sample])
}
FPKM.merged.DE$ResponderType <- as.factor(FPKM.merged.DE$ResponderType)

FPKM.pca <- prcomp(FPKM.merged.DE[,-ncol(FPKM.merged.DE)], center = TRUE, scale. = TRUE)
summary(FPKM.pca)
ggbiplot(FPKM.pca, labels=rownames(FPKM.merged.DE), var.axes = FALSE, groups = FPKM.merged.DE$ResponderType, ellipse = TRUE, choices = c(1,2))

#=======================================================================================================================================================
#Separate samples based on batch
FPKM.merged2 <- as.data.frame(t(FPKM.merged))
FPKM.merged2$batch <- batch
FPKM.merged.batch1 <- FPKM.merged2[FPKM.merged2$batch == 1,]
FPKM.merged.batch2 <- FPKM.merged2[FPKM.merged2$batch == 2,]
FPKM.merged.batch1 <- FPKM.merged.batch1[,-ncol(FPKM.merged.batch1)]
FPKM.merged.batch2 <- FPKM.merged.batch2[,-ncol(FPKM.merged.batch2)]

FPKM.merged.batch1$ResponderType <- NA
for (i in 1:nrow(FPKM.merged.batch1)){
  FPKM.merged.batch1$ResponderType[i] <- as.character(responder.subtype$Group[rownames(FPKM.merged.batch1)[i] == responder.subtype$Sample])
}
FPKM.merged.batch1$ResponderType <- as.factor(FPKM.merged.batch1$ResponderType)

FPKM.merged.batch2$ResponderType <- NA
for (i in 1:nrow(FPKM.merged.batch2)){
  FPKM.merged.batch2$ResponderType[i] <- as.character(responder.subtype$Group[rownames(FPKM.merged.batch2)[i] == responder.subtype$Sample])
}
FPKM.merged.batch2$ResponderType <- as.factor(FPKM.merged.batch2$ResponderType)

FPKM.pca.batch1 <- prcomp(FPKM.merged.batch1[,-ncol(FPKM.merged.batch1)], center = TRUE, scale. = F)
summary(FPKM.pca.batch1)
ggbiplot(FPKM.pca.batch1, labels=rownames(FPKM.merged.batch1), var.axes = FALSE, groups = FPKM.merged.batch1$ResponderType, ellipse = TRUE, choices = c(1,2))

FPKM.pca.batch2 <- prcomp(FPKM.merged.batch2[,-ncol(FPKM.merged.batch2)], center = TRUE, scale. = F)
summary(FPKM.pca.batch2)
ggbiplot(FPKM.pca.batch2, labels=rownames(FPKM.merged.batch2), var.axes = FALSE, groups = FPKM.merged.batch2$ResponderType, ellipse = TRUE, choices = c(2,5))

#=======================================================================================================================================================
DEgenes <- read.table(file = 'C:\\Users\\medsthe\\Desktop\\edgeR.txt')
#DEgenes$x <- gsub("\\..*","",DEgenes$x)
FPKM.merged <- as.data.frame(t(FPKM.merged))
FPKM.merged.DE <- FPKM.merged[colnames(FPKM.merged) %in% DEgenes$x]

# Dataset
FPKM.merged.DE$ResponderType <- NA
for (i in 1:nrow(FPKM.merged.DE)){
  FPKM.merged.DE$ResponderType[i] <- as.character(responder.subtype$Group[rownames(FPKM.merged.DE)[i] == responder.subtype$Sample])
}
FPKM.merged.DE$ResponderType <- as.factor(FPKM.merged.DE$ResponderType)

FPKM.pca <- prcomp(FPKM.merged.DE[,-ncol(FPKM.merged.DE)], center = TRUE, scale. = TRUE)
summary(FPKM.pca)
ggbiplot(FPKM.pca, labels=rownames(FPKM.merged.DE), var.axes = FALSE, groups = FPKM.merged.DE$ResponderType, ellipse = TRUE, choices = c(1,2))



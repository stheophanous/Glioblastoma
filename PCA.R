rm(list = ls())
library(ggbiplot)
library(data.table)
# ResponderSubtype
clinical.data <- read.table(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\ClinicalDataFinal.txt', sep ="\t", header=TRUE)
responder.subtype <- as.data.frame(clinical.data[,c(1,14)])
colnames(responder.subtype) <- c('Sample', 'Group')
responder.subtype$Group <- as.factor(responder.subtype$Group)
# Dataset
FPKM.merged <- as.data.frame(fread(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\mergedFPKM_ALLsamples.txt'))
rownames(FPKM.merged) <- FPKM.merged$V1
FPKM.merged <- FPKM.merged[,-1]
FPKM.merged$ResponderType <- NA
for (i in 1:nrow(FPKM.merged)){
  FPKM.merged$ResponderType[i] <- as.character(responder.subtype$Group[rownames(FPKM.merged)[i] == responder.subtype$Sample])
}
FPKM.merged$ResponderType <- as.factor(FPKM.merged$ResponderType)

FPKM.pca <- prcomp(FPKM.merged[,-ncol(FPKM.merged)], center = TRUE, scale. = TRUE)
summary(FPKM.pca)
ggbiplot(FPKM.pca, labels=rownames(FPKM.merged), var.axes = FALSE, groups = FPKM.merged$ResponderType, ellipse = TRUE, choices = c(1,2))

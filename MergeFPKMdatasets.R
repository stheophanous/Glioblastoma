rm(list = ls())

#Loading the required libraries
library('stringr')
library('data.table')
library('dplyr')
library('tidyverse')

#Load datasets
# 3. ResponderSubtype
clinical.data <- read.table(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\ClinicalDataFinal.txt', sep ="\t", header=TRUE)
responder.subtype <- as.data.frame(clinical.data[,c(1,14)])
colnames(responder.subtype) <- c('Sample', 'Group')
responder.subtype$Group <- as.factor(responder.subtype$Group)
#Old dataset
primary.FPKM.tr <- as.data.frame(fread(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\ForStelios\\PrimaryGBM_FPKMsClean.csv', sep = ','))
colnames(primary.FPKM.tr)[1] <- 'Sample'
colnames(primary.FPKM.tr) <- sub('\\..*', "", colnames(primary.FPKM.tr))
#New dataset
primary.FPKM2 <- as.data.frame(fread(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\ForStelios\\GSEA\\NewK_STDLOCAL_FPKMs.txt', sep = '\t'))
rownames(primary.FPKM2) <- primary.FPKM2$EnsID
colnames(primary.FPKM2) <- sub("_P_FPKM", "", colnames(primary.FPKM2))
primary.FPKM2 <- primary.FPKM2[,-c(1,2)]
primary.FPKM2 <- primary.FPKM2[,c(1:8)]
primary.FPKM.tr2 <- as.data.frame(t(primary.FPKM2))
sample.names <- data.frame(Sample = rownames(primary.FPKM.tr2))
primary.FPKM.tr2 <- cbind(sample.names, primary.FPKM.tr2)
#Non standard treatment dataset
primary.FPKM3 <- as.data.frame(fread(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\ForStelios\\GSEA\\NewK_LOCALnSTD_FPKMs.txt', sep = '\t'))
rownames(primary.FPKM3) <- primary.FPKM3$EnsID
colnames(primary.FPKM3) <- sub("_P_FPKM", "", colnames(primary.FPKM3))
primary.FPKM3 <- primary.FPKM3[,-c(1,2)]
primary.FPKM3 <- primary.FPKM3[,c(1:10)]
primary.FPKM.tr3 <- as.data.frame(t(primary.FPKM3))
sample.names2 <- data.frame(Sample = rownames(primary.FPKM.tr3))
primary.FPKM.tr3 <- cbind(sample.names2, primary.FPKM.tr3)
#Merge datasets
merged.dataset <- dplyr::bind_rows(primary.FPKM.tr, primary.FPKM.tr2)
merged.dataset <- dplyr::bind_rows(merged.dataset, primary.FPKM.tr3)
merged.dataset[is.na(merged.dataset)] <- 0
merged.dataset <- merged.dataset[ order(merged.dataset$Sample), ]
rownames(merged.dataset) <- merged.dataset$Sample
merged.dataset <- merged.dataset[,-1]

#Preprocessing
#preProcValues <- preProcess(primary.FPKM.tr, method = c("knnImpute","center","scale"))
#FPKM.processed <- predict(preProcValues, primary.FPKM.tr)
#x <- scale(merged.dataset[,-nearZeroVar(merged.dataset)])
FPKM.processed <- as.data.frame(merged.dataset)

#Add responder type
FPKM.processed$ResponderType <- NA
for (i in 1:nrow(FPKM.processed)){
  FPKM.processed$ResponderType[i] <- responder.subtype$Group[rownames(FPKM.processed)[i] == responder.subtype$Sample]
}
FPKM.processed$ResponderType <- as.factor(FPKM.processed$ResponderType)

savetable <- as.data.frame(FPKM.processed)
write.table(savetable, file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\mergedFPKM_ALLsamples.txt', col.names = T, row.names = T, quote = F)



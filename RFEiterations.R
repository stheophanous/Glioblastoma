rm(list = ls())

#Loading the required libraries
library('caret')
library('caretEnsemble')
library('dplyr')
library('tidyr')

#Responder subtype
clinical.data <- read.table(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\ClinicalDataFinal.txt', sep ="\t", header=TRUE)
responder.subtype <- as.data.frame(clinical.data[,c(1,14)])
colnames(responder.subtype) <- c('Sample', 'Group')
responder.subtype$Group <- as.factor(responder.subtype$Group)
#All samples
FPKM.merged <- as.data.frame(fread(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\mergedFPKM_ALLsamples.txt'))
rownames(FPKM.merged) <- FPKM.merged$V1
FPKM.merged <- FPKM.merged[,-1]

#Remove highly correlated genes
#corrm <- cor(FPKM.merged)
#remove <- findCorrelation(corrm, cutoff = 0.9, exact = FALSE)
#FPKM.merged <- FPKM.merged[,-c(remove)]

#Add responder subtype
FPKM.merged$ResponderType <- NA
for (i in 1:nrow(FPKM.merged)){
  FPKM.merged$ResponderType[i] <- as.character(responder.subtype$Group[rownames(FPKM.merged)[i] == responder.subtype$Sample])
}
FPKM.merged$ResponderType <- as.factor(FPKM.merged$ResponderType)

results <- data.frame(matrix(NA, nrow = 10, ncol = 1))
i <- 1
for (i in 1:100){
  set.seed(i)
  #Split to training and test sets
  index <- createDataPartition(FPKM.merged$ResponderType, p=0.83, list=FALSE)
  trainSet <- FPKM.merged[ index,]
  testSet <- FPKM.merged[-index,]
  #Feature selection using rfe
  control <- rfeControl(functions = rfFuncs, method = "repeatedcv", repeats = 5, verbose = TRUE)
  outcomeName<-'ResponderType'
  predictors<-names(trainSet)[!names(trainSet) %in% outcomeName]
  subset.sizes <- c(5,10,15,20)
  ResponderType_Pred_Profile <- rfe(trainSet[,predictors], trainSet[,outcomeName], rfeControl = control, sizes = subset.sizes)
  features <- head(ResponderType_Pred_Profile$optVariables, 10)
  iter <- data.frame(features)
  results <- cbind(results, iter)
  cat('Iteration', i, 'DONE')
  i <- i+1
}

results.clean <- results[,-1]
results.clean <- data.frame(lapply(results.clean, as.character), stringsAsFactors=FALSE)
pos1 <- as.list(results.clean[1,])
pos1 <- rep(pos1, 10)
pos2 <- as.list(results.clean[2,])
pos2 <- rep(pos2, 9)
pos3 <- as.list(results.clean[3,])
pos3 <- rep(pos3, 8)
pos4 <- as.list(results.clean[4,])
pos4 <- rep(pos4, 7)
pos5 <- as.list(results.clean[5,])
pos5 <- rep(pos5, 6)
pos6 <- as.list(results.clean[6,])
pos6 <- rep(pos6, 5)
pos7 <- as.list(results.clean[7,])
pos7 <- rep(pos7, 4)
pos8 <- as.list(results.clean[8,])
pos8 <- rep(pos8, 3)
pos9 <- as.list(results.clean[9,])
pos9 <- rep(pos9, 2)
pos10 <- as.list(results.clean[10,])
scored.list <- as.character(c(pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8, pos9, pos10))
most.frequent <- sort(table(scored.list),decreasing=TRUE)
most.frequent <- as.data.frame(most.frequent)
most.frequent$score <- NA
most.frequent$score <- most.frequent$Freq/10
write.table(most.frequent, file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\mostfrequentRFEgenes100rep2.txt')


results.list <- c(results.clean[1,], results.clean[2,], results.clean[3,], results.clean[4,], results.clean[5,], results.clean[6,], results.clean[7,], results.clean[8,], results.clean[9,], results.clean[10,])
results.list <- as.character(results.list)
sort(table(results.list), decreasing = T)








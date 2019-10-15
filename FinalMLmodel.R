rm(list = ls())

#Loading the required libraries
library('caret')
library('caretEnsemble')
library('dplyr')
library('data.table')
set.seed(1)

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
corrm <- cor(FPKM.merged)
remove <- findCorrelation(corrm, cutoff = 0.8, exact = FALSE)
FPKM.merged <- FPKM.merged[,-c(remove)]

#Add responder subtype
FPKM.merged$ResponderType <- NA
for (i in 1:nrow(FPKM.merged)){
  FPKM.merged$ResponderType[i] <- as.character(responder.subtype$Group[rownames(FPKM.merged)[i] == responder.subtype$Sample])
}
FPKM.merged$ResponderType <- as.factor(FPKM.merged$ResponderType)

#Split to training and test sets
index <- createDataPartition(FPKM.merged$ResponderType, p=0.83, list=FALSE)
trainSet <- FPKM.merged[ index,]
testSet <- FPKM.merged[-index,]

#Feature selection using rfe
control <- rfeControl(functions = rfFuncs, method = "repeatedcv", repeats = 10, verbose = TRUE)
outcomeName<-'ResponderType'
predictors<-names(trainSet)[!names(trainSet) %in% outcomeName]
subset.sizes <- c(4, 8, 12, 16, 20, 25, 50, 100, 200)
ResponderType_Pred_Profile <- rfe(trainSet[,predictors], trainSet[,outcomeName], rfeControl = control, sizes = subset.sizes)
ResponderType_Pred_Profile
#**Remember to change number of features**
features <- head(ResponderType_Pred_Profile$optVariables, 100) #Change number of features depending on the prediction profile
trainSet <- trainSet %>% select(features)


trainSet$ResponderType <- NA
for (i in 1:nrow(trainSet)){
  trainSet$ResponderType[i] <- as.character(clinical.data$ResponderType[rownames(trainSet)[i] == clinical.data$Patient.ID])
}
trainSet$ResponderType <- factor(trainSet$ResponderType)

#Classification using Neural Networks
#Accuracy and Kappa
control <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = 'all', classProbs = T)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="Accuracy", trControl=control, tuneLength = 5)
#LogLoss
control2 <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = 'all', classProbs = T, summaryFunction = mnLogLoss)
nn.fit2 <- train(ResponderType~., data=trainSet, method="nnet", metric="logLoss", trControl=control2, tuneLength = 5)
#ROC and prediction on test set
control3 <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = 'all', classProbs = T, summaryFunction = twoClassSummary)
nn.fit3 <- train(ResponderType~., data=trainSet, method="nnet", metric="ROC", trControl=control3, tuneLength = 5)

nn.fit
nn.fit2
nn.fit3

#Prediction on validation/test set
predictors <- colnames(trainSet)[-length(trainSet)]
nn.pred <- predict(object = nn.fit,testSet[,predictors])
#Checking the accuracy of classification of test set
confusionMatrix(testSet$ResponderType,nn.pred)

#========================================================
#Effect of changing seed on the prediction accuracy (on test set)
results <- data.frame(seed=numeric(), Accuracy=numeric(), Sensitivity=numeric(), Specificity=numeric(), BalancedAccuracy=numeric())
j <- 1
#ROC over 1000 seeds
for (j in 1:10){
  set.seed(j)
  #Spliting training set into two parts based on outcome: 75% and 25%
  index <- createDataPartition(FPKM.merged$ResponderType, p=0.75, list=FALSE)
  trainSet <- FPKM.merged[ index,]
  testSet <- FPKM.merged[-index,]
  trainSet <- trainSet %>% select(features)
  trainSet$ResponderType <- NA
  for (i in 1:nrow(trainSet)){
    trainSet$ResponderType[i] <- as.character(FPKM.merged$ResponderType[rownames(trainSet)[i] == rownames(FPKM.merged)])
  }
  trainSet$ResponderType <- factor(trainSet$ResponderType)
  #ROC
  control <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = 'all', classProbs = T, summaryFunction = twoClassSummary)
  nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="ROC", trControl=control, tuneLength = 5)
  #Prediction on test set
  predictors <- colnames(trainSet)[-length(trainSet)]
  nn.pred <- predict(object = nn.fit,testSet[,predictors])
  #Checking the accuracy of classification of test set
  confmatrix <- confusionMatrix(testSet$ResponderType,nn.pred)
  seed <- j
  accuracy <- confmatrix$overall[[1]]
  sensitivity <- confmatrix$byClass[[1]]
  specificity <- confmatrix$byClass[[2]]
  balaccuracy <- confmatrix$byClass[[11]]
  iter <- data.frame(seed, accuracy, sensitivity, specificity, balaccuracy)
  results <- rbind(results, iter)
  cat('Iteration', j, 'DONE')
  j <- j+1
}
summary(results)
hist <- hist(results$balaccuracy)
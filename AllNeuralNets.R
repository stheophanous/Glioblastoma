rm(list = ls())

#Loading the required libraries
library('caret')
library('caretEnsemble')
library('stringr')
library('data.table')
library('dplyr')
library('tidyverse')
#Setting the seed
set.seed(1)

#Import dataset (old and new samples together)
full.dataset <- fread(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\mergedFPKM.txt')
full.dataset <- full.dataset %>% remove_rownames %>% column_to_rownames(var="V1")
full.dataset$ResponderType <- as.factor(full.dataset$ResponderType)

#Spliting training set into two parts based on outcome: 75% and 25%
index <- createDataPartition(full.dataset$ResponderType, p=0.75, list=FALSE)
trainSet <- full.dataset[ index,]
testSet <- full.dataset[-index,]

#Select features and add responder type again
#commonnew <- read.table(file='C:\\Users\\medsthe\\Desktop\\Glioblastoma\\DifferentialGeneExpressionAnalysis\\OnlineTool\\comparison_resultsMERGED\\edger_deseq2_overlap_only.txt', sep= "\t")
#commonnew <- as.character(commonnew$V1)
#commonnew <- as.list(commonnew)[-1]

commonnew <- read.table(file='C:\\Users\\medsthe\\Desktop\\Glioblastoma\\DifferentialGeneExpressionAnalysis\\FinalResults\\commonmerged.txt', sep= " ")
commonnew <- as.list(as.character(commonnew$x))

featurelist <- unique(commonnew)
featurelist <- sub('\\..*', "", featurelist)

trainSet <- trainSet %>% select(featurelist)

trainSet$ResponderType <- NA
for (i in 1:nrow(trainSet)){
  trainSet$ResponderType[i] <- as.character(full.dataset$ResponderType[rownames(trainSet)[i] == rownames(full.dataset)])
}
trainSet$ResponderType <- factor(trainSet$ResponderType)

#========================================================================================================================================
#Accuracy and Kappa
control <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = 'all', classProbs = T)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="Accuracy", trControl=control, tuneLength = 5)
avnn.fit <- train(ResponderType~., data=trainSet, method="avNNet", metric="Accuracy", trControl=control, tuneLength = 5)
mlpnn.fit <- train(ResponderType~., data=trainSet, method="monmlp", metric="Accuracy", trControl=control, tuneLength = 5)
pcann.fit <- train(ResponderType~., data=trainSet, method="pcaNNet", metric="Accuracy", trControl=control, tuneLength = 5)

nn.fit
avnn.fit
mlpnn.fit
pcann.fit

#======================================================================================================================================
#LogLoss
control <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = 'all', classProbs = T, summaryFunction = mnLogLoss)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="logLoss", trControl=control, tuneLength = 5)
avnn.fit <- train(ResponderType~., data=trainSet, method="avNNet", metric="logLoss", trControl=control, tuneLength = 5)
mlpnn.fit <- train(ResponderType~., data=trainSet, method="monmlp", metric="logLoss", trControl=control, tuneLength = 5)
pcann.fit <- train(ResponderType~., data=trainSet, method="pcaNNet", metric="logLoss", trControl=control, tuneLength = 5)

nn.fit
avnn.fit
mlpnn.fit
pcann.fit

#======================================================================================================================================
#ROC and prediction on test set
control <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = 'all', classProbs = T, summaryFunction = twoClassSummary)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="ROC", trControl=control, tuneLength = 5)
avnn.fit <- train(ResponderType~., data=trainSet, method="avNNet", metric="ROC", trControl=control, tuneLength = 5)
mlpnn.fit <- train(ResponderType~., data=trainSet, method="monmlp", metric="ROC", trControl=control, tuneLength = 5)
pcann.fit <- train(ResponderType~., data=trainSet, method="pcaNNet", metric="ROC", trControl=control, tuneLength = 5)

nn.fit
avnn.fit
mlpnn.fit
pcann.fit

#======================================================================================================================================
#Prediction on validation/test set
predictors <- colnames(trainSet)[-length(trainSet)]
nn.pred <- predict(object = nn.fit,testSet[,predictors])
avnn.pred <- predict(object = avnn.fit,testSet[,predictors])
mlpnn.pred <- predict(object = mlpnn.fit,testSet[,predictors])
pcann.pred <- predict(object = pcann.fit,testSet[,predictors])

#Checking the accuracy of classification of test set
confusionMatrix(testSet$ResponderType,nn.pred)
confusionMatrix(testSet$ResponderType,avnn.pred)
confusionMatrix(testSet$ResponderType,mlpnn.pred)
confusionMatrix(testSet$ResponderType,pcann.pred)

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
svm.fit <- train(ResponderType~., data=trainSet, method="svmRadial", metric="Accuracy", trControl=control, tuneLength = 10)
svm2.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="Accuracy", trControl=control, tuneLength = 10)
gbm.fit <- train(ResponderType~., data=trainSet, method="gbm", metric="Accuracy", trControl=control, tuneLength = 10, verbose = F, bag.fraction = 0.65)
rf.fit <- train(ResponderType~., data=trainSet, method="rf", metric="Accuracy", trControl=control, tuneLength = 10)
bab.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="Accuracy", trControl=control, tuneLength = 10)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="Accuracy", trControl=control, tuneLength = 10, verbose = F)
cart.fit <- train(ResponderType~., data=trainSet, method="rpart", metric="Accuracy", trControl=control, tuneLength = 10)
c5.fit <- train(ResponderType~., data=trainSet, method="C5.0", metric="Accuracy", trControl=control, tuneLength = 10)
adaboost.fit <- train(ResponderType~., data=trainSet, method="adaboost", metric="Accuracy", trControl=control, tuneLength = 3)
glm.fit <- train(ResponderType~., data=trainSet, method="glm", metric="Accuracy", trControl=control, tuneLength = 10)
glm2.fit <- train(ResponderType~., data=trainSet, method="glmStepAIC", metric="Accuracy", trControl=control, tuneLength = 10, verbose = F)
mars.fit <- train(ResponderType~., data=trainSet, method="earth", metric="Accuracy", trControl=control, tuneLength = 10)
knn.fit <- train(ResponderType~., data=trainSet, method="kknn", metric="Accuracy", trControl=control, tuneLength = 10)
nb.fit <- train(ResponderType~., data=trainSet, method="naive_bayes", metric="Accuracy", trControl=control, tuneLength = 10)

svm.fit
svm2.fit
gbm.fit
rf.fit
bab.fit
nn.fit
cart.fit
c5.fit
adaboost.fit
glm.fit
glm2.fit
mars.fit
knn.fit
nb.fit

#======================================================================================================================================
#LogLoss
control <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = 'all', classProbs = T, summaryFunction = mnLogLoss)
svm.fit <- train(ResponderType~., data=trainSet, method="svmRadial", metric="logLoss", trControl=control, tuneLength = 10)
svm2.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="logLoss", trControl=control, tuneLength = 10)
gbm.fit <- train(ResponderType~., data=trainSet, method="gbm", metric="logLoss", trControl=control, tuneLength = 10, verbose = F, bag.fraction = 0.65)
rf.fit <- train(ResponderType~., data=trainSet, method="rf", metric="logLoss", trControl=control, tuneLength = 10)
bab.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="logLoss", trControl=control, tuneLength = 10)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="logLoss", trControl=control, tuneLength = 10, verbose = F)
cart.fit <- train(ResponderType~., data=trainSet, method="rpart", metric="logLoss", trControl=control, tuneLength = 10)
c5.fit <- train(ResponderType~., data=trainSet, method="C5.0", metric="logLoss", trControl=control, tuneLength = 10)
adaboost.fit <- train(ResponderType~., data=trainSet, method="adaboost", metric="logLoss", trControl=control, tuneLength = 3)
glm.fit <- train(ResponderType~., data=trainSet, method="glm", metric="logLoss", trControl=control, tuneLength = 10)
glm2.fit <- train(ResponderType~., data=trainSet, method="glmStepAIC", metric="logLoss", trControl=control, tuneLength = 10, verbose = F)
mars.fit <- train(ResponderType~., data=trainSet, method="earth", metric="logLoss", trControl=control, tuneLength = 10)
knn.fit <- train(ResponderType~., data=trainSet, method="kknn", metric="logLoss", trControl=control, tuneLength = 10)
nb.fit <- train(ResponderType~., data=trainSet, method="naive_bayes", metric="logLoss", trControl=control, tuneLength = 10)

svm.fit
svm2.fit
gbm.fit
rf.fit
bab.fit
nn.fit
cart.fit
c5.fit
adaboost.fit
glm.fit
glm2.fit
mars.fit
knn.fit
nb.fit

#======================================================================================================================================
#ROC and prediction on test set
control <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = 'all', classProbs = T, summaryFunction = twoClassSummary)
svm.fit <- train(ResponderType~., data=trainSet, method="svmRadial", metric="ROC", trControl=control, tuneLength = 10)
svm2.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="ROC", trControl=control, tuneLength = 10)
gbm.fit <- train(ResponderType~., data=trainSet, method="gbm", metric="ROC", trControl=control, tuneLength = 10, verbose = F, bag.fraction = 0.65)
rf.fit <- train(ResponderType~., data=trainSet, method="rf", metric="ROC", trControl=control, tuneLength = 10)
bab.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="ROC", trControl=control, tuneLength = 10)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="ROC", trControl=control, tuneLength = 10, verbose = F)
cart.fit <- train(ResponderType~., data=trainSet, method="rpart", metric="ROC", trControl=control, tuneLength = 10)
c5.fit <- train(ResponderType~., data=trainSet, method="C5.0", metric="ROC", trControl=control, tuneLength = 10)
adaboost.fit <- train(ResponderType~., data=trainSet, method="adaboost", metric="ROC", trControl=control, tuneLength = 3)
glm.fit <- train(ResponderType~., data=trainSet, method="glm", metric="ROC", trControl=control, tuneLength = 10)
glm2.fit <- train(ResponderType~., data=trainSet, method="glmStepAIC", metric="ROC", trControl=control, tuneLength = 10, verbose = F)
mars.fit <- train(ResponderType~., data=trainSet, method="earth", metric="ROC", trControl=control, tuneLength = 10)
knn.fit <- train(ResponderType~., data=trainSet, method="kknn", metric="ROC", trControl=control, tuneLength = 10)
nb.fit <- train(ResponderType~., data=trainSet, method="naive_bayes", metric="ROC", trControl=control, tuneLength = 10)

svm.fit
svm2.fit
gbm.fit
rf.fit
bab.fit
nn.fit
cart.fit
c5.fit
adaboost.fit
glm.fit
glm2.fit
mars.fit
knn.fit
nb.fit

#Prediction on validation/test set
predictors <- colnames(trainSet)[-length(trainSet)]
svm.pred <- predict(object = svm.fit,testSet[,predictors])
svm2.pred <- predict(object = svm2.fit,testSet[,predictors])
gbm.pred <- predict(object = gbm.fit,testSet[,predictors])
rf.pred <- predict(object = rf.fit,testSet[,predictors])
bab.pred <- predict(object = bab.fit,testSet[,predictors])
nn.pred <- predict(object = nn.fit,testSet[,predictors])
cart.pred <- predict(object = cart.fit,testSet[,predictors])
c5.pred <- predict(object = c5.fit,testSet[,predictors])
adaboost.pred <- predict(object = adaboost.fit,testSet[,predictors])
glm.pred <- predict(object = glm.fit,testSet[,predictors])
glm2.pred <- predict(object = glm2.fit,testSet[,predictors])
mars.pred <- predict(object = mars.fit,testSet[,predictors])
knn.pred <- predict(object = knn.fit,testSet[,predictors])
nb.pred <- predict(object = nb.fit,testSet[,predictors])

#Checking the accuracy of classification of test set
confusionMatrix(testSet$ResponderType,svm.pred)
confusionMatrix(testSet$ResponderType,svm2.pred)
confusionMatrix(testSet$ResponderType,gbm.pred)
confusionMatrix(testSet$ResponderType,rf.pred)
confusionMatrix(testSet$ResponderType,bab.pred)
confusionMatrix(testSet$ResponderType,nn.pred)
confusionMatrix(testSet$ResponderType,cart.pred)
confusionMatrix(testSet$ResponderType,c5.pred)
confusionMatrix(testSet$ResponderType,adaboost.pred)
confusionMatrix(testSet$ResponderType,glm.pred)
confusionMatrix(testSet$ResponderType,glm2.pred)
confusionMatrix(testSet$ResponderType,mars.pred)
confusionMatrix(testSet$ResponderType,knn.pred)
confusionMatrix(testSet$ResponderType,nb.pred)



#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#Lucy's Method
#Split dataset into training and testing
trainSet <- full.dataset[unlist(tapply(1:nrow(full.dataset),full.dataset$ResponderType,function(x) sample(x,14))),]
testSet <- full.dataset[setdiff(rownames(full.dataset),rownames(trainSet)),]

trainSet <- trainSet %>% select(featurelist)

trainSet$ResponderType <- NA
for (i in 1:nrow(trainSet)){
  trainSet$ResponderType[i] <- as.character(full.dataset$ResponderType[rownames(trainSet)[i] == rownames(full.dataset)])
}
trainSet$ResponderType <- factor(trainSet$ResponderType)

#Accuracy and Kappa
control <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = 'all', classProbs = T)
svm.fit <- train(ResponderType~., data=trainSet, method="svmRadial", metric="Accuracy", trControl=control, tuneLength = 10)
svm2.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="Accuracy", trControl=control, tuneLength = 10)
gbm.fit <- train(ResponderType~., data=trainSet, method="gbm", metric="Accuracy", trControl=control, tuneLength = 10, verbose = F, bag.fraction = 0.9)
rf.fit <- train(ResponderType~., data=trainSet, method="rf", metric="Accuracy", trControl=control, tuneLength = 10)
bab.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="Accuracy", trControl=control, tuneLength = 10)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="Accuracy", trControl=control, tuneLength = 10, verbose = F)
cart.fit <- train(ResponderType~., data=trainSet, method="rpart", metric="Accuracy", trControl=control, tuneLength = 10)
c5.fit <- train(ResponderType~., data=trainSet, method="C5.0", metric="Accuracy", trControl=control, tuneLength = 10)
adaboost.fit <- train(ResponderType~., data=trainSet, method="adaboost", metric="Accuracy", trControl=control, tuneLength = 3)
glm.fit <- train(ResponderType~., data=trainSet, method="glm", metric="Accuracy", trControl=control, tuneLength = 10)
glm2.fit <- train(ResponderType~., data=trainSet, method="glmStepAIC", metric="Accuracy", trControl=control, tuneLength = 10, verbose = F)
mars.fit <- train(ResponderType~., data=trainSet, method="earth", metric="Accuracy", trControl=control, tuneLength = 10)
knn.fit <- train(ResponderType~., data=trainSet, method="kknn", metric="Accuracy", trControl=control, tuneLength = 10)
nb.fit <- train(ResponderType~., data=trainSet, method="naive_bayes", metric="Accuracy", trControl=control, tuneLength = 10)

svm.fit
svm2.fit
gbm.fit
rf.fit
bab.fit
nn.fit
cart.fit
c5.fit
adaboost.fit
glm.fit
glm2.fit
mars.fit
knn.fit
nb.fit

#======================================================================================================================================
#LogLoss
control <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = 'all', classProbs = T, summaryFunction = mnLogLoss)
svm.fit <- train(ResponderType~., data=trainSet, method="svmRadial", metric="logLoss", trControl=control, tuneLength = 10)
svm2.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="logLoss", trControl=control, tuneLength = 10)
gbm.fit <- train(ResponderType~., data=trainSet, method="gbm", metric="logLoss", trControl=control, tuneLength = 10, verbose = F, bag.fraction = 0.65)
rf.fit <- train(ResponderType~., data=trainSet, method="rf", metric="logLoss", trControl=control, tuneLength = 10)
bab.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="logLoss", trControl=control, tuneLength = 10)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="logLoss", trControl=control, tuneLength = 10, verbose = F)
cart.fit <- train(ResponderType~., data=trainSet, method="rpart", metric="logLoss", trControl=control, tuneLength = 10)
c5.fit <- train(ResponderType~., data=trainSet, method="C5.0", metric="logLoss", trControl=control, tuneLength = 10)
adaboost.fit <- train(ResponderType~., data=trainSet, method="adaboost", metric="logLoss", trControl=control, tuneLength = 3)
glm.fit <- train(ResponderType~., data=trainSet, method="glm", metric="logLoss", trControl=control, tuneLength = 10)
glm2.fit <- train(ResponderType~., data=trainSet, method="glmStepAIC", metric="logLoss", trControl=control, tuneLength = 10, verbose = F)
mars.fit <- train(ResponderType~., data=trainSet, method="earth", metric="logLoss", trControl=control, tuneLength = 10)
knn.fit <- train(ResponderType~., data=trainSet, method="kknn", metric="logLoss", trControl=control, tuneLength = 10)
nb.fit <- train(ResponderType~., data=trainSet, method="naive_bayes", metric="logLoss", trControl=control, tuneLength = 10)

svm.fit
svm2.fit
gbm.fit
rf.fit
bab.fit
nn.fit
cart.fit
c5.fit
adaboost.fit
glm.fit
glm2.fit
mars.fit
knn.fit
nb.fit

#======================================================================================================================================
#ROC and prediction on test set
control <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = 'all', classProbs = T, summaryFunction = twoClassSummary)
svm.fit <- train(ResponderType~., data=trainSet, method="svmRadial", metric="ROC", trControl=control, tuneLength = 10)
svm2.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="ROC", trControl=control, tuneLength = 10)
gbm.fit <- train(ResponderType~., data=trainSet, method="gbm", metric="ROC", trControl=control, tuneLength = 10, verbose = F, bag.fraction = 0.65)
rf.fit <- train(ResponderType~., data=trainSet, method="rf", metric="ROC", trControl=control, tuneLength = 10)
bab.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="ROC", trControl=control, tuneLength = 10)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="ROC", trControl=control, tuneLength = 10, verbose = F)
cart.fit <- train(ResponderType~., data=trainSet, method="rpart", metric="ROC", trControl=control, tuneLength = 10)
c5.fit <- train(ResponderType~., data=trainSet, method="C5.0", metric="ROC", trControl=control, tuneLength = 10)
adaboost.fit <- train(ResponderType~., data=trainSet, method="adaboost", metric="ROC", trControl=control, tuneLength = 3)
glm.fit <- train(ResponderType~., data=trainSet, method="glm", metric="ROC", trControl=control, tuneLength = 10)
glm2.fit <- train(ResponderType~., data=trainSet, method="glmStepAIC", metric="ROC", trControl=control, tuneLength = 10, verbose = F)
mars.fit <- train(ResponderType~., data=trainSet, method="earth", metric="ROC", trControl=control, tuneLength = 10)
knn.fit <- train(ResponderType~., data=trainSet, method="kknn", metric="ROC", trControl=control, tuneLength = 10)
nb.fit <- train(ResponderType~., data=trainSet, method="naive_bayes", metric="ROC", trControl=control, tuneLength = 10)

svm.fit
svm2.fit
gbm.fit
rf.fit
bab.fit
nn.fit
cart.fit
c5.fit
adaboost.fit
glm.fit
glm2.fit
mars.fit
knn.fit
nb.fit

#Prediction on validation/test set
predictors <- colnames(trainSet)[-length(trainSet)]
svm.pred <- predict(object = svm.fit,testSet[,predictors])
svm2.pred <- predict(object = svm2.fit,testSet[,predictors])
gbm.pred <- predict(object = gbm.fit,testSet[,predictors])
rf.pred <- predict(object = rf.fit,testSet[,predictors])
bab.pred <- predict(object = bab.fit,testSet[,predictors])
nn.pred <- predict(object = nn.fit,testSet[,predictors])
cart.pred <- predict(object = cart.fit,testSet[,predictors])
c5.pred <- predict(object = c5.fit,testSet[,predictors])
adaboost.pred <- predict(object = adaboost.fit,testSet[,predictors])
glm.pred <- predict(object = glm.fit,testSet[,predictors])
glm2.pred <- predict(object = glm2.fit,testSet[,predictors])
mars.pred <- predict(object = mars.fit,testSet[,predictors])
knn.pred <- predict(object = knn.fit,testSet[,predictors])
nb.pred <- predict(object = nb.fit,testSet[,predictors])

#Checking the accuracy of classification of test set
confusionMatrix(testSet$ResponderType,svm.pred)
confusionMatrix(testSet$ResponderType,svm2.pred)
confusionMatrix(testSet$ResponderType,gbm.pred)
confusionMatrix(testSet$ResponderType,rf.pred)
confusionMatrix(testSet$ResponderType,bab.pred)
confusionMatrix(testSet$ResponderType,nn.pred)
confusionMatrix(testSet$ResponderType,cart.pred)
confusionMatrix(testSet$ResponderType,c5.pred)
confusionMatrix(testSet$ResponderType,adaboost.pred)
confusionMatrix(testSet$ResponderType,glm.pred)
confusionMatrix(testSet$ResponderType,glm2.pred)
confusionMatrix(testSet$ResponderType,mars.pred)
confusionMatrix(testSet$ResponderType,knn.pred)
confusionMatrix(testSet$ResponderType,nb.pred)

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#Dave's Method - leave-one-out CV
trainSet <- full.dataset %>% select(featurelist)
trainSet$ResponderType <- NA
for (i in 1:nrow(trainSet)){
  trainSet$ResponderType[i] <- as.character(full.dataset$ResponderType[rownames(trainSet)[i] == rownames(full.dataset)])
}
trainSet$ResponderType <- factor(trainSet$ResponderType)

control <- trainControl(method="LOOCV", number=10, savePredictions = 'all', classProbs = T)
svm.fit <- train(ResponderType~., data=trainSet, method="svmRadial", metric="Accuracy", trControl=control, tuneLength = 10)
svm2.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="Accuracy", trControl=control, tuneLength = 10)
gbm.fit <- train(ResponderType~., data=trainSet, method="gbm", metric="Accuracy", trControl=control, tuneLength = 10, verbose = F, bag.fraction = 0.65)
rf.fit <- train(ResponderType~., data=trainSet, method="rf", metric="Accuracy", trControl=control, tuneLength = 10)
bab.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="Accuracy", trControl=control, tuneLength = 10)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="Accuracy", trControl=control, tuneLength = 10, verbose = F)
cart.fit <- train(ResponderType~., data=trainSet, method="rpart", metric="Accuracy", trControl=control, tuneLength = 10)
c5.fit <- train(ResponderType~., data=trainSet, method="C5.0", metric="Accuracy", trControl=control, tuneLength = 10)
adaboost.fit <- train(ResponderType~., data=trainSet, method="adaboost", metric="Accuracy", trControl=control, tuneLength = 3)
glm.fit <- train(ResponderType~., data=trainSet, method="glm", metric="Accuracy", trControl=control, tuneLength = 10)
glm2.fit <- train(ResponderType~., data=trainSet, method="glmStepAIC", metric="Accuracy", trControl=control, tuneLength = 10, verbose = F)
mars.fit <- train(ResponderType~., data=trainSet, method="earth", metric="Accuracy", trControl=control, tuneLength = 10)
knn.fit <- train(ResponderType~., data=trainSet, method="kknn", metric="Accuracy", trControl=control, tuneLength = 10)
nb.fit <- train(ResponderType~., data=trainSet, method="naive_bayes", metric="Accuracy", trControl=control, tuneLength = 10)

svm.fit
svm2.fit
gbm.fit
rf.fit
bab.fit
nn.fit
cart.fit
c5.fit
adaboost.fit
glm.fit
glm2.fit
mars.fit
knn.fit
nb.fit

#======================================================================================================================================
#LogLoss
control <- trainControl(method="LOOCV", number=10, savePredictions = 'all', classProbs = T, summaryFunction = mnLogLoss)
svm.fit <- train(ResponderType~., data=trainSet, method="svmRadial", metric="logLoss", trControl=control, tuneLength = 10)
svm2.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="logLoss", trControl=control, tuneLength = 10)
gbm.fit <- train(ResponderType~., data=trainSet, method="gbm", metric="logLoss", trControl=control, tuneLength = 10, verbose = F, bag.fraction = 0.65)
rf.fit <- train(ResponderType~., data=trainSet, method="rf", metric="logLoss", trControl=control, tuneLength = 10)
bab.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="logLoss", trControl=control, tuneLength = 10)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="logLoss", trControl=control, tuneLength = 10, verbose = F)
cart.fit <- train(ResponderType~., data=trainSet, method="rpart", metric="logLoss", trControl=control, tuneLength = 10)
c5.fit <- train(ResponderType~., data=trainSet, method="C5.0", metric="logLoss", trControl=control, tuneLength = 10)
adaboost.fit <- train(ResponderType~., data=trainSet, method="adaboost", metric="logLoss", trControl=control, tuneLength = 3)
glm.fit <- train(ResponderType~., data=trainSet, method="glm", metric="logLoss", trControl=control, tuneLength = 10)
glm2.fit <- train(ResponderType~., data=trainSet, method="glmStepAIC", metric="logLoss", trControl=control, tuneLength = 10, verbose = F)
mars.fit <- train(ResponderType~., data=trainSet, method="earth", metric="logLoss", trControl=control, tuneLength = 10)
knn.fit <- train(ResponderType~., data=trainSet, method="kknn", metric="logLoss", trControl=control, tuneLength = 10)
nb.fit <- train(ResponderType~., data=trainSet, method="naive_bayes", metric="logLoss", trControl=control, tuneLength = 10)

svm.fit
svm2.fit
gbm.fit
rf.fit
bab.fit
nn.fit
cart.fit
c5.fit
adaboost.fit
glm.fit
glm2.fit
mars.fit
knn.fit
nb.fit

#======================================================================================================================================
#ROC and prediction on test set
control <- trainControl(method="LOOCV", number=10, savePredictions = 'all', classProbs = T, summaryFunction = twoClassSummary)
svm.fit <- train(ResponderType~., data=trainSet, method="svmRadial", metric="ROC", trControl=control, tuneLength = 10)
svm2.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="ROC", trControl=control, tuneLength = 10)
gbm.fit <- train(ResponderType~., data=trainSet, method="gbm", metric="ROC", trControl=control, tuneLength = 10, verbose = F, bag.fraction = 0.65)
rf.fit <- train(ResponderType~., data=trainSet, method="rf", metric="ROC", trControl=control, tuneLength = 10)
bab.fit <- train(ResponderType~., data=trainSet, method="svmLinear", metric="ROC", trControl=control, tuneLength = 10)
nn.fit <- train(ResponderType~., data=trainSet, method="nnet", metric="ROC", trControl=control, tuneLength = 10, verbose = F)
cart.fit <- train(ResponderType~., data=trainSet, method="rpart", metric="ROC", trControl=control, tuneLength = 10)
c5.fit <- train(ResponderType~., data=trainSet, method="C5.0", metric="ROC", trControl=control, tuneLength = 10)
adaboost.fit <- train(ResponderType~., data=trainSet, method="adaboost", metric="ROC", trControl=control, tuneLength = 3)
glm.fit <- train(ResponderType~., data=trainSet, method="glm", metric="ROC", trControl=control, tuneLength = 10)
glm2.fit <- train(ResponderType~., data=trainSet, method="glmStepAIC", metric="ROC", trControl=control, tuneLength = 10, verbose = F)
mars.fit <- train(ResponderType~., data=trainSet, method="earth", metric="ROC", trControl=control, tuneLength = 10)
knn.fit <- train(ResponderType~., data=trainSet, method="kknn", metric="ROC", trControl=control, tuneLength = 10)
nb.fit <- train(ResponderType~., data=trainSet, method="naive_bayes", metric="ROC", trControl=control, tuneLength = 10)

svm.fit
svm2.fit
gbm.fit
rf.fit
bab.fit
nn.fit
cart.fit
c5.fit
adaboost.fit
glm.fit
glm2.fit
mars.fit
knn.fit
nb.fit















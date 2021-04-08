###############################################
###### Analysis of high-dimensional data ######
######            Project                ######
######  Bosser, Dokuz, Thardak, Yunus    ######
###############################################

##Useful libraries 
library("factoextra")
library(glmnet)
library("ROCR")
library(DMwR)
library('PRROC')
library('pls')
library(hydroGOF)
library("rlist")
library(caret)
library(tidyverse)

##Loading datasets, observation and data cleaning

#DATA SET-I

load("X_GSE21374.rda")
GeneExpression<-t(X_GSE21374) 
dim(GeneExpression)

head(GeneExpression[,1:10])

colnames(GeneExpression)
row.names(GeneExpression)

class(GeneExpression)

#DATA SET-II

load("RejectionStatus.rda") 
dim(RejectionStatus)
head(RejectionStatus)

row.names(RejectionStatus) <- RejectionStatus$Patient_ID
RejectionStatus <- RejectionStatus[2]

class(RejectionStatus)

##Merging the two datasets

MergedData <- merge(RejectionStatus, GeneExpression,  by="row.names")
dim(MergedData)

head(MergedData[,1:10])

rownames(MergedData) <- MergedData[,"Row.names"]
MergedData <- MergedData[,2:ncol(MergedData)]
head(MergedData[,1:5])

#**********************************
#RESEARCH Q1
#**********************************

##All data PCA

PCA <- prcomp(MergedData[,2:ncol(MergedData)], center = TRUE, scale. = TRUE)

##Converting the response variable from continuous to class variable

MergedData[,1] <- as.character(MergedData[,1])


##Visualisation for class variable.

plot(PCA$x[,1], PCA$x[,2], xlab = "PC1 (17%)", ylab = "PC2 (7.5%)", main = "PC1 / PC2 - Plot", col = c("grey","red"))
plot(PCA$x[,1], PCA$x[,3], xlab = "PC1 (17%)", ylab = "PC3 (5%)", main = "PC1 / PC3 - Plot", col = c("grey","red"))
plot(PCA$x[,2], PCA$x[,3], xlab = "PC2 (7.5%)", ylab = "PC3 (5%)", main = "PC2 / PC3 - Plot", col = c("grey","red"))
plot(PCA$x[,3], PCA$x[,4], xlab = "PC3 (5%)", ylab = "PC4 (4%)", main = "PC3 / PC4 - Plot", col = c("grey","red"))


##Percentage of variances explained by each principal component

fviz_screeplot(PCA, addlabels = TRUE, ylim = c(0, 50))


fviz_pca_ind(PCA, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = MergedData[,1], 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Rejection Status") +
  ggtitle("2D PCA-plot") +
  theme(plot.title = element_text(hjust = 5))

#From plots, it was concluded that it is hard to say that the variablitiy is associated with kidney rejection.

#**********************************
  #RESEARCH Q2
#**********************************

##Change binary response to integer type (method doesn't accept characters)

MergedData[,1] <- as.integer(MergedData[,1])

##FDR method

library('fuzzySim')
fdr <- FDR(MergedData, 
    sp.cols = 1, 
    var.cols = 2:54676, 
    family = 'binomial', 
    q = 0.1)

#When controlling the FDR at the 10% level, 17805 variables are selected,
#and 36870 are discarded. When choosing a gene with a given q-value (let's say 5%),
#we can conclude that 5% of the genes that present a p-value more extreme (than the one
#of the gene selected) are false positives. Controlling at the 10% level makes sure that
#if a gene Y has a q-value of 0.1, then 10% of the genes with a more extreme p-value than the gene Y
#are false positive 

#**********************************
#RESEARCH Q3
#***************************************************************************************************************************

##Converting the response variable back from class to continuous variable to use in the model

MergedData[,1] <- as.numeric(MergedData[,1])

##Scaling the x-variables 

Scaled <- scale(MergedData[,2:ncol(MergedData)], center = TRUE, scale = TRUE)
MergedData <- data.frame(cbind(MergedData$Reject_Status,Scaled))
colnames(MergedData)[1] <- 'Reject_Status'

##Splitting data as train and test
#Select 70% of data as train, 30% as test data

set.seed(1)
n<- dim(MergedData)[1]
nTrain <- round(0.7*n)
trainID <- sample(1:n, nTrain)


##Training data

train <-  MergedData[trainID,]

train$Reject_Status <- as.factor((train$Reject_Status))

##Random Oversampling minority class 
# --> Best results obtained with ration 1:1 (around 5% increase of AUCPRC)

Positive <- which(train$Reject_Status == 1) 
Negative <- which(train$Reject_Status == 0) 
numPositive <- length(which(train$Reject_Status == 1))
numNegative <- length(which(train$Reject_Status == 0))
nInstances <- round(1*(numNegative - numPositive))

NewSamples_idx <- sample(Positive, nInstances, replace=TRUE)

NewSamples <- train[NewSamples_idx, ]

NewData <- rbind(train, NewSamples)

##SMOTE oversampling
# --> Worse perf. than rd.oversampling (even than no resampling)
# --> SMOTE is known not to perform good on very high-dimensional data

#Neg <- train[Negative, ]

#NewData <- SMOTE(form = Reject_Status ~ ., data = train, 
# perc.over = 200, 
# perc.under = 200)

#NewData <- rbind(NewData, Neg)


x <- as.matrix(MergedData[2:ncol(MergedData)])
y <- as.matrix(MergedData[1])


##Original train data - No oversampling methods used

trainX <- x[trainID,]
trainY <- y[trainID,]

##Balanced train data - Random oversampling

trainX_OS <- as.matrix(NewData[2:ncol(NewData)])
trainY_OS <- as.matrix(NewData[1])

##Test data

testX <- x[-trainID,] 
testY <- y[-trainID,]

#%27 are one (not once oversampling is performed)

table(y)
table(trainY)
table(trainY_OS)

#Among observations in test data 59 are zero, 26 are 1 (~30%). 

table(testY)

#*********************#
# LASSO REGRESSION    #
#*********************#

##Selection of Gamma 

m.cv <- cv.glmnet(x=trainX_OS, y=trainY_OS, alpha=1, type.measure = "auc" , family = "binomial", nfolds = 10)
plot(m.cv, xlab = "log(gamma)")
m.cv

#lambdamin:0.0046
#lambda1se:0.0392

##Refit the data with selected penalty parameter= min

m <- glmnet(x=trainX_OS, y=trainY_OS, alpha=1,lambda = m.cv$lambda.min ,family = "binomial")
plot(summary(coef(m))[,1], summary(coef(m))[,3], cex=1, pch=16,
     xlab="variables", ylab="beta-hat", cex.lab=1.5, cex.axis=1.5)
abline(h=0, col=2, lwd=2)

#non-zero beta estimates in the model

summary(coef(m))

#Probabilities and predicted values

pred.m <- prediction(predict(m, newx = testX, s= m.cv$lambda.min, type = "response"), testY)
probs.m <- predict(m, newx = testX, s= m.cv$lambda.min, type = "response")
pred.val.m <- predict(m, newx = testX, s= m.cv$lambda.min, type = "class")

#Area under the ROC curve

performance(pred.m, 'auc')

#Confusion Matrix

table(testY, pred.val.m)

#ROC and PR curves

fg.m <- probs.m[testY == 1]
bg.m <- probs.m[testY == 0]

pr.m <- pr.curve(scores.class0 = fg.m,
               scores.class1 = bg.m,
               curve=TRUE)
roc.m <- roc.curve(scores.class0 = fg.m,
                 scores.class1 = bg.m,
                 curve=TRUE)

plot(pr.m) #0.605
plot(roc.m) #0.780

# No oversampling : AUROC : 0.77, AUPRC : 0.55 
# Oversampling : AUROC : 0.78, AUPCR : 0.605

#######

##Refit the data with selected penalty parameter= 1se

m2 <- glmnet(x=trainX_OS, y=trainY_OS, alpha=1,lambda = m.cv$lambda.1se ,family = "binomial")
plot(summary(coef(m2))[,1], summary(coef(m2))[,3], cex=1, pch=16,
     xlab="variables", ylab="beta-hat", cex.lab=1.5, cex.axis=1.5)
abline(h=0, col=2, lwd=2)

#non-zero beta estimates in the model. m2 is a simpler model comparing to m1

summary(coef(m2))

#Probabilities and predicted values

pred.m2 <- prediction(predict(m2, newx = testX, s= m.cv$lambda.1se, type = "response"), testY)
probs.m2 <- predict(m2, newx = testX, s= m.cv$lambda.1se, type = "response")
pred.val.m2 <- predict(m2, newx = testX, s= m.cv$lambda.1se, type = "class")

#Aurea under the ROC curve

performance(pred.m2, 'auc')

#Confusion Matrix 

table(testY, pred.val.m2)

#ROC and PR curves

fg.m2 <- probs.m2[testY == 1]
bg.m2 <- probs.m2[testY == 0]

pr.m2 <- pr.curve(scores.class0 = fg.m2,
                 scores.class1 = bg.m2,
                 curve=TRUE)
roc.m2 <- roc.curve(scores.class0 = fg.m2,
                   scores.class1 = bg.m2,
                   curve=TRUE)

plot(pr.m2) #0.59
plot(roc.m2) #0782

# No oversampling : AUROC : 0.769, AUPCR : 0.54
# Oversampling : AUROC : 0.783 , AUPCR : 0.597


#***************************************************************************************************************************


#*******************#
#RIDGE REGRESSION   #
#*******************#

##Selection of Gamma 

r.cv <- cv.glmnet(x=trainX, y=trainY, alpha=0, type.measure = "auc" , family = "binomial")
plot(r.cv, xlab = "log(gamma)")
r.cv

#lambdamin: 41.6
#lambda1se: 243.6

##Refit the data with selected penalty parameter= min

r <- glmnet(x=trainX, y=trainY, alpha=0,lambda = r.cv$lambda.min ,family = "binomial")
plot(summary(coef(r))[,1], summary(coef(r))[,3], cex=1, pch=16,
     xlab="variables", ylab="beta-hat", cex.lab=1.5, cex.axis=1.5)
abline(h=0, col=2, lwd=2)

#Many non-zero beta estimates in the model, however mostly they are very close to zero. One of them is not.
summary(coef(r))

#Probabilities and predicted values

pred.r <- prediction(predict(r, newx = testX, s= r.cv$lambda.min, type = "response"), testY)
probs.r <- predict(r, newx = testX, s= r.cv$lambda.min, type = "response")
pred.val.r <- predict(r, newx = testX, s= r.cv$lambda.min, type = "class")

#Area under the ROC curve 

performance(pred.r, 'auc')

#Confusion Matrix 

table(testY, pred.val.r)

#ROC and PR curves

fg.r <- probs.r[testY == 1]
bg.r <- probs.r[testY == 0]

pr.r <- pr.curve(scores.class0 = fg.r,
                  scores.class1 = bg.r,
                  curve=TRUE)
roc.r <- roc.curve(scores.class0 = fg.r,
                    scores.class1 = bg.r,
                    curve=TRUE)

plot(pr.r) #0.595
plot(roc.r) #0.762

#No oversampling : AUROC : 0.76, AUPCR : 0.593
#Oversampling : AUROC : 0.755, AUPCR : 0.584

#Refit the data with selected penalty parameter= 1se

r2 <- glmnet(x=trainX, y=trainY, alpha=0,lambda = r.cv$lambda.1se ,family = "binomial")
plot(summary(coef(r2))[,1], summary(coef(r2))[,3], cex=1, pch=16,
     xlab="variables", ylab="beta-hat", cex.lab=1.5, cex.axis=1.5)
abline(h=0, col=2, lwd=2)

#Mmany non-zero beta estimates in the model, however mostly they are very close to zero. One of them is not.

summary(coef(r2))

#Probabilities and Predicted values

pred.r2 <- prediction(predict(r2, newx = testX, s= r.cv$lambda.1se, type = "response"), testY)
probs.r2 <- predict(r2, newx = testX, s= r.cv$lambda.1se, type = "response")
pred.val.r2 <- predict(r2, newx = testX, s= r.cv$lambda.1se, type = "class")

#Area under the ROC curve 

performance(pred.r2, 'auc')

#Confusion Matrix

table(testY, pred.val.r2)

#ROC and PR curves

fg.r2 <- probs.r2[testY == 1]
bg.r2 <- probs.r2[testY == 0]

pr.r2 <- pr.curve(scores.class0 = fg.r2,
               scores.class1 = bg.r2,
               curve=TRUE)
roc.r2 <- roc.curve(scores.class0 = fg.r2,
                 scores.class1 = bg.r2,
                 curve=TRUE)

plot(pr.r2) #0.57
plot(roc.r2) #0.755

# No over-sampling : AUROC : 0.75, AUPCR : 0.576
# Oversampling : AUROC: 0.75, AUPCR : 0.6

#*********************************#
# PRINCIPAL COMPONENT REGRESSION  #
#*********************************#

scaled.data <- scale(MergedData[,2:ncol(MergedData)], center = TRUE, scale = TRUE)
Xsvd <- svd(scaled.data)

##Plots of eigen values by principal components

plot(Xsvd$d^2, type = "b", ylab="eigenvalue", xlab = "PC", cex=2, cex.axis=1.5, cex.lab=1.5)
barplot(Xsvd$d^2, names.arg = 1:282,xlab = "PC", ylab="eigenvalue", cex.lab=1.5)
barplot(cumsum(Xsvd$d^2) / sum(Xsvd$d^2), names.arg = 1:282,xlab = "PC", ylab="cumulative eigen value", cex.lab=1.5)
abline(h=0.8,col=2,lty=5)

#Based on eigenvalues, k=120 dimensions are to be selected as 80% of variance is retained.

set.seed(1)
k <- 120
Uk <- Xsvd$u[,1:k]
Dk <- diag(Xsvd$d[1:k])
Zk <- Uk%*%Dk
Vk <- Xsvd$v[,1:k]

##Training data

trainZ <- Zk[trainID,]
trainY <- y[trainID,]

##Test data

testZ <- (Zk[-trainID,]) 
testY <- y[-trainID, ]

suppressMessages(library(boot))
suppressMessages(library(pROC))

##Selection of number of PCs

max.n.comps <- 120
cv.pcr.errors <- rep(NA, max.n.comps)

##MSE function definition

MSE <- function(obs, pred){ 
  mse <- mean((obs - pred)^2)
  return(mse) }

##Compute MSE of model given number of PC ranging from 1 to 120

for(i in 1:max.n.comps){
  pcr.train <- pcr(trainY ~ trainZ, ncomp = i) # Calculate the predictions
  ypred <- predict(pcr.train, newdata = testZ) # Calculate the extra-sample error 
  cv.pcr.errors[i] <- MSE(testY,ypred)
}

##Visualization of PCs vs MSE

nPC.at.min <- which.min(cv.pcr.errors) 
plot(1:max.n.comps, cv.pcr.errors, type = "l",
     ylab = "MSE", xlab = "number PCs") 
abline(v = nPC.at.min, col = "red")

##nPC.at.min which is 22 giving minimum MSE is selected

k <- 22
Uk <- Xsvd$u[,1:k]
Dk <- diag(Xsvd$d[1:k])
Zk <- Uk%*%Dk
Vk <- Xsvd$v[,1:k]

##Training data

trainZ <- Zk[trainID,]
trainY <- y[trainID,]

##Test data

testZ <- (Zk[-trainID,]) 
testY <- y[-trainID, ]

##Recombine to dataframes for the glm method 

test.pcr <- data.frame(cbind(testY, testZ))
train.pcr <- data.frame(cbind(trainY, trainZ))

colnames(test.pcr) <- paste("X", 1:dim(test.pcr)[2], sep = "")
colnames(train.pcr) <- paste("X", 1:dim(train.pcr)[2], sep = "")

##Oversampling 
#No evidence of increase in performance, wathever the ratio (ratio <- 0) for no oversampling

Positive <- which(train.pcr$X1 == 1)
Negative <- which(train$X1 == 0)
numPositive <- length(which(train$X1 == 1))
numNegative <- length(which(train$X1 == 0))
ratio <- 0
nInstances <- round(ratio*(numNegative - numPositive))

NewSamples_idx <- sample(Positive, nInstances, replace=TRUE)

NewSamples <- train.pcr[NewSamples_idx, ]

train.pcr.over <- data.frame(rbind(train.pcr, NewSamples))

##Fit the model

pcr_mod <- glm(X1 ~ ., data=train.pcr.over, 
               family = "binomial", 
               control = list(maxit = 50))

summary(pcr_mod)

#Probabilities and predictions

pred.pcr <- prediction(predict(pcr_mod, newdata  = test.pcr, type = "response"), testY)
probs.pcr <- predict(pcr_mod, newdata = test.pcr, type = "response")

#Area under the ROC curve 

performance(pred.pcr, 'auc')

##Convert probabilities to predicted values (given threshold)

pred.val.pcr <- probs.pcr 
threshold <- 0.5
id_val <- which(pred.val.pcr > threshold)
pred.val.pcr[id_val] <- 1
pred.val.pcr[-id_val] <- 0

#Confusion matrix

table(testY, pred.val.pcr)

#ROC and PRC plots 

fg.pcr <- probs.pcr[testY == 1]
bg.pcr <- probs.pcr[testY == 0]

pr.pcr <- pr.curve(scores.class0 = fg.pcr,
                 scores.class1 = bg.pcr,
                 curve=TRUE)
roc.pcr <- roc.curve(scores.class0 = fg.pcr,
                   scores.class1 = bg.pcr,
                   curve=TRUE)

plot(pr.pcr) #0.529
plot(roc.pcr) #0.714
#****************************************************************************************************************************



#***********************************************#
# Threshold selection - Selected model is Lasso #
#***********************************************#

##Refit Lasso Regression using CV

#If you run this part indenpendently, select the set.seed, otherwise you'll have random results
set.seed(1)

##Function to generate Cross-validation folds

fold_generator <- function(dataframe, n_fold){
  n<- dim(dataframe)[1]
  ratio <- round(n/n_fold) 
  folds <- list()
  Idx <- 1:n
  Idx2 <- Idx
  for (k in 1:(n_fold-1)){
    foldID <- sample(Idx2, ratio)
    Idx2 <- Idx2[!Idx2 %in% foldID]
    folds <- list.append(folds, foldID)
  }
  folds <- list.append(folds, Idx2)
  
  return(folds)
}

##Function that returns the sensitivity and specificity given probabilities

Sens_Spec <- function(proba, test_val){
  c.vec <- seq(0,1,0.01)
  sens.vec <- rep(NA, length(c.vec))
  spec.vec <- rep(NA, length(c.vec))
  for (i in seq(length(c.vec))){
    predicted <- ifelse(proba > c.vec[i], 1, 0)
    actual <- test_val
    confusion_matrix <- table(actual, predicted)
    if(length(which(predicted == 1)) == 0){
      sens.vec[i] <- 0
      spec.vec[i] <- 1
    }
    else if(length(which(predicted == 0)) == 0){
      sens.vec[i] <- 1
      spec.vec[i] <- 0 
    }
    else {
      sensitivity <- confusion_matrix[2,2] / (confusion_matrix[2,2] + confusion_matrix[2,1])
      specificity <- confusion_matrix[1,1] / (confusion_matrix[1,1] + confusion_matrix[1,2])
      sens.vec[i] <- sensitivity
      spec.vec[i] <- specificity
    }
  }
  SS <- list(sens.vec, spec.vec)
  return(SS) 
}

##Generate 10-folds for CV

folds <- fold_generator(train, 10)

c.vals = rep(NA, length(folds)) #Empty vector for best threshold within each fold
c.vec <- seq(0,1,0.01) #Thresholds

##Cross-validation loop with random oversampling

for(i in 1:length(folds)){
  
  #Train/validation split (Oversampling)
  
  train_CV <- train[-folds[[i]], ]
  
  Positive <- which(train_CV$Reject_Status == 1)
  Negative <- which(train_CV$Reject_Status == 0)
  numPositive <- length(which(train_CV$Reject_Status == 1))
  numNegative <- length(which(train_CV$Reject_Status == 0))
  nInstances <- round(1*(numNegative - numPositive))
  NewSamples_idx <- sample(Positive, nInstances, replace=TRUE)
  NewSamples <- train_CV[NewSamples_idx, ]
  NewData <- rbind(train_CV, NewSamples)
  
  trainX_OS_CV <- as.matrix(NewData[2:ncol(NewData)])
  trainY_OS_CV <- as.matrix(NewData[1])
  
  testX_CV <- trainX[folds[[i]], ]
  testY_CV <- trainY[folds[[i]]]
  
  #Fit the models
  
  m2.CV <- glmnet(x=trainX_OS_CV, y=trainY_OS_CV, alpha=1,lambda = m.cv$lambda.1se ,family = "binomial")
  
  #Compute pobabilities
  
  proba.CV <- predict(m2.CV, newx = testX_CV, s= m.cv$lambda.1se, type = "response")
  
  #Compute sensi/speci for each thresholds
  
  SS <- Sens_Spec(proba.CV, testY_CV)
  
  #Sum sensi/speci for each thresholds
  
  SS.sum <- SS[[1]] + SS[[2]]
  
  #Compute Younden index and find corresponding optimal threshold
  
  Younden.index <- max(SS.sum)
  
  YI.i <- which(SS.sum == Younden.index)
  
  opt.threshold <- c.vec[YI.i]
  
  c.vals[i] <- opt.threshold
  
}

##Compute mean threshold

c.mean <- mean(c.vals) #Mean c = 0.38

##Function for misclassification error

MISERR <- function(obs, pred, cutoff = 0.5){
  ypred <- as.numeric(pred > cutoff) # translates TRUE/FALSE to 1/0 
  tab <- table(obs, ypred)
  miserr <- 1 - sum(diag(tab))/sum(tab)
  return(miserr)
}

##Refit the model using all the training set

m2.final <- glmnet(x=trainX_OS, y=trainY_OS, alpha=1,lambda = m.cv$lambda.1se ,family = "binomial")


##Performance on training set

c <- c.mean 
predicted_train <- ifelse(predict(m2.final, newx = trainX, s= m.cv$lambda.1se, type = "response") > c, 1, 0)
actual_train <- trainY
confusion_matrix <- table(actual_train, predicted_train)
confusion_matrix

sensitivity_train <- confusion_matrix[2,2] / (confusion_matrix[2,2] + confusion_matrix[2,1])
specificity_train  <- confusion_matrix[1,1] / (confusion_matrix[1,1] + confusion_matrix[1,2])
accuracy_train <- (confusion_matrix[1,1] + confusion_matrix[2,2]) / (confusion_matrix[1,1]+confusion_matrix[1,2]+confusion_matrix[2,1]+confusion_matrix[2,2])
#sens:1
#spec: 0.94
#acc:0.95

MISERR(actual_train, predicted_train, cutoff = c)
#misserr: 0.05

##Performance on test set

predicted_test <- ifelse(predict(m2.final, newx = testX, s= m.cv$lambda.1se, type = "response") > c, 1, 0)
actual_test <- testY
confusion_matrix <- table(actual_test, predicted_test)
confusion_matrix

sensitivity_test <- confusion_matrix[2,2] / (confusion_matrix[2,2] + confusion_matrix[2,1])
specificity_test <- confusion_matrix[1,1] / (confusion_matrix[1,1] + confusion_matrix[1,2])
accuracy_test <- (confusion_matrix[1,1] + confusion_matrix[2,2]) / (confusion_matrix[1,1]+confusion_matrix[1,2]+confusion_matrix[2,1]+confusion_matrix[2,2])
#sens:0.73
#spec: 0.64
#acc:0.67

MISERR(actual_test, predicted_test, cutoff = c)
#misserr: 0.33

##Importance of features in the model

ftrs <- varImp(m2, lambda = m.cv$lambda.1se)
m.f <- data.frame(rownames(ftrs)[which(ftrs != 0)], ftrs[which(ftrs != 0),])
m.f


m.f <- m.f %>% 
  rename(
    'Gene' = 'rownames.ftrs..which.ftrs....0..',
    'beta' = 'ftrs.which.ftrs....0....'
  )

m.f <- m.f[order(m.f$beta, decreasing = TRUE),]

##Most important genes in the model

barplot(m.f$beta[1:20], names.arg = m.f$Gene[1:20], ylim = c(0,0.5), las = 2, xlab = 'Genes', ylab = 'beta coef.')


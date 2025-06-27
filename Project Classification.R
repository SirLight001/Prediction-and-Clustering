#######################################################################################################
# PROJECT LECTURER: PROF. ADAM ZAGDANSKI                                                              #
# COURSE TITLE: DATA MINING                                                                           #
# PROJECT TOPIC: BREAST CANCER PREDICTION USING ROC BASED FEATURE FILTERING:                          #
#                A COMPARATIVE STUDY OF CLASSIFICATION METHODS                                        #
# STUDENTS AND ID: SEGUN LIGHT JEGEDE (257389) and ISAAC AKOJI PAUL (257388)                          #
#######################################################################################################
#-------------------------------
## LOADING THE DATA INTO R
#----------------------------------------------------------------------------------------------------------------------------
#Reading the data into r and renaming the variable names in a readable manner.
col.names=c("id_number","clump_thickness","uniformity_cell_size","uniformity_cell_shape",
            "marginal_adhesion","single_epithelial_cell_size","bare_nuclei",
            "bland_chromatin","normal_nucleoli","mitoses","class")
bcw <- read.csv("C:/Users/jeged/Downloads/breast-cancer-wisconsin.data", header=FALSE, col.names=col.names)
#View(bcw)
attach(bcw)
#-------------------------------
## DATA PREPARATION AND CLEANING
#----------------------------------------------------------------------------------------------------------------------------
library(DataExplorer)

#Checking the data type of each columns
str(bcw)
bcw$class = as.factor(bcw$class) #convert the class to factor with "2" as benign and "4" as malignant
levels(bcw$class)[levels(bcw$class)=="2"] <- "benign"
levels(bcw$class)[levels(bcw$class)=="4"] <- "malignant"
bcw$bare_nuclei = suppressWarnings(as.numeric(bcw$bare_nuclei)) #convert bare_nuclei to numeric because it is a measurement like others
bcw$id_number = as.character(bcw$id_number) #id_number is nothing but a string of cells identification number
str(bcw) #Every attributes is now in their respective perfect form.

#Handling missing data by deleting the corresponding rows if the missing observations are not too much
t(introduce(bcw))
sum(is.na(bcw)) #check for missing observations
plot_intro(bcw)
plot_missing(bcw) 
bcw<-na.omit(bcw) #2.29% of the bare_nuclei measurement are missing variables, thus, we decided to remove any form of missing observation.
sum(is.na(bcw))
plot_missing(bcw) 
nrow(bcw) #The data reduced from 699 to 683, we suppose we did not lose too much information, just about 3%
#View(bcw)

#-------------------------------
##EXPLORATORY DATA ANALYSIS
#----------------------------------------------------------------------------------------------------------------------------
#Describing the Grouping Variable
library(tidyverse)
library(scales)
bcwnew <- bcw %>% group_by(class) %>% 
  summarize(count = n()) %>%  # count records by species
  mutate(percentage = count/sum(count))  # find percent of total

ggplot(bcwnew, aes(class, percentage, fill = class)) + 
  geom_bar(stat='identity') + 
  geom_text(aes(label=scales::percent(percentage)), position = position_stack(vjust = .5))+
  scale_y_continuous(labels = scales::percent)

#Describing the features
my.summary <- function(df)
{
  results <- matrix(, nrow = 9, ncol = ncol(df))
  for (i in 1:ncol(df)){
    X=df[,i]
    results[,i] <- rbind(min(X),quantile(X,0.25), median(X), mean(X), quantile(X,0.75), max(X), var(X), sd(X), IQR(X))
  }
    rownames(results) <- c("min", "Q1", "median", "mean", "Q3", "max", "var", "sd", "IQR")
    colnames(results) <-names(df)
    return(results)
}
ms<-my.summary(bcw[,2:10]) 
ms
write.table(ms, file = "summary statistics.txt", sep = ",", quote = FALSE, row.names = F)

#construct the plots three by three
#construct the histogram plots
library(ggpubr)
ha<-fg<-ggplot(bcw, aes(x = clump_thickness, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
hb<-ggplot(bcw, aes(x = uniformity_cell_size, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
hc<-ggplot(bcw, aes(x = uniformity_cell_shape, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
hd<-fg<-ggplot(bcw, aes(x = marginal_adhesion, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
ggarrange(ha,hb,hc,hd,labels = c("A", "B", "C","D"),ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

hd<-fg<-ggplot(bcw, aes(x = marginal_adhesion, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
he<-ggplot(bcw, aes(x = single_epithelial_cell_size, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
hf<-ggplot(bcw, aes(x = bare_nuclei, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
hg<-fg<-ggplot(bcw, aes(x = bland_chromatin, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
hh<-ggplot(bcw, aes(x = normal_nucleoli, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
ggarrange(he,hf,hg,hh,labels = c("E", "F","G", "H"),ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

hi<-ggplot(bcw, aes(x = mitoses, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
ggarrange(hi,labels = c("I"),ncol = 2, nrow = 2) #, common.legend = TRUE, legend="bottom")

'plot_histogram(bcw)'

#construct the density plots
plot_density(bcw[2:4])
plot_density(bcw[5:7])
plot_density(bcw[8:10])

#construct the normal qq plot
plot_qq(bcw[2:4])
plot_qq(bcw[5:7])
plot_qq(bcw[8:10])

#construct the barplots
bpa <- ggplot(bcw, aes(x = class, y = clump_thickness))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
bpb <- ggplot(bcw, aes(x = class, y = uniformity_cell_size))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
bpc <- ggplot(bcw, aes(x = class, y = uniformity_cell_shape))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
bpd <- ggplot(bcw, aes(x = class, y = marginal_adhesion))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
ggarrange(bpa,bpb,bpc,bpd,labels = c("A", "B", "C","D"),ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

bpe <- ggplot(bcw, aes(x = class, y = single_epithelial_cell_size))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
bpf <- ggplot(bcw, aes(x = class, y = bare_nuclei))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
bpg <- ggplot(bcw, aes(x = class, y = bland_chromatin))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
bph <- ggplot(bcw, aes(x = class, y = normal_nucleoli))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
ggarrange(bpe,bpf,bpg,bph,labels = c("E", "F", "G","H"),ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

bpi <- ggplot(bcw, aes(x = class, y = mitoses))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
ggarrange(bpi,labels = c("I"),ncol = 2, nrow = 2) #, common.legend = TRUE, legend="bottom")

'plot_boxplot(bcw, by="class")'

#Construct Boxplot without grouping
library(reshape)
bcwData <- melt(bcw)
par(mar=c(10,7,1,1))
boxplot(data=bcwData, value~variable, las=2)
pairs(bcw[2:10], pch = 21, bg = c("#d95f02", "#7570b3")[unclass(bcw$class)])

plot_correlation(bcw, type = "continuous") #correlation plot
#we can also do sum data preparation here; it is an error if the measurement does not fall within 1-10 i.e Minimum and Maximum

#-------------------------------
## NORMALIZATION
#----------------------------------------------------------------------------------------------------------------------------
a<-(bcw$clump_thickness-min(bcw$clump_thickness))/(max(bcw$clump_thickness)-min(bcw$clump_thickness))
b<-(bcw$uniformity_cell_size-min(bcw$uniformity_cell_size))/(max(bcw$uniformity_cell_size)-min(bcw$uniformity_cell_size))
c<-(bcw$uniformity_cell_shape-min(bcw$uniformity_cell_shape))/(max(bcw$uniformity_cell_shape)-min(bcw$uniformity_cell_shape))
d<-(bcw$marginal_adhesion-min(bcw$marginal_adhesion))/(max(bcw$marginal_adhesion)-min(bcw$marginal_adhesion))
e<-(bcw$single_epithelial_cell_size-min(bcw$single_epithelial_cell_size))/(max(bcw$single_epithelial_cell_size)-min(bcw$single_epithelial_cell_size))
f<-(bcw$bare_nuclei-min(bcw$bare_nuclei))/(max(bcw$bare_nuclei)-min(bcw$bare_nuclei))
g<-(bcw$bland_chromatin-min(bcw$bland_chromatin))/(max(bcw$bland_chromatin)-min(bcw$bland_chromatin))
h<-(bcw$normal_nucleoli-min(bcw$normal_nucleoli))/(max(bcw$normal_nucleoli)-min(bcw$normal_nucleoli))
i<-(bcw$mitoses-min(bcw$mitoses))/(max(bcw$mitoses)-min(bcw$mitoses))
bcw1<-data.frame(a,b,c,d,e,f,g,h,i,bcw[,11])
colnames(bcw1)<-c("clump_thickness","uniformity_cell_size","uniformity_cell_shape","marginal_adhesion",
                      "single_epithelial_cell_size","bare_nuclei","bland_chromatin","normal_nucleoli","mitoses","class")
str(bcw1)
View(bcw1)


#-------------------------------
## FEATURE SELECTION
library(mlbench)
library(ROCR)
#----------------------------------------------------------------------------------------------------------------------------
#***************ROC-based assessment of the importance of variables/features****************************
features <- c("clump_thickness","uniformity_cell_size","uniformity_cell_shape",
              "marginal_adhesion","single_epithelial_cell_size","bare_nuclei",
              "bland_chromatin","normal_nucleoli","mitoses")
f.iter <- 1 #iterator
# auxiliary function
fun <- function(f)
{
  X <- bcw1[,f]
  pred <- prediction(X, labels = bcw1$class)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf, col=f.iter, add=(f.iter>1))
  f.iter <- f.iter + 1
  auc <- performance(pred, "auc")@y.values[[1]]
}
features.auc <- sapply(features, fun)
features.auc
legend("bottomright", legend=paste(features,": AUC=",round(features.auc,2)), lwd=1, col=1:4, bg="azure")
grid()
lines(c(0,1), c(0,1), lty=2, lwd=2) #random classifier
title("ROC-based assessment of the variables importance")

#-------------------------------
## CLASSIFICATION
## ALL FEATURES
#----------------------------------------------------------------------------------------------------------------------------
class.labels <-  bcw1$class # class
n <- length(class.labels) # number of objects
K <- length(levels(class.labels)) # number of classes

library(class)
learning.set.index <- sample(1:n,2/3*n) # random split of data into training and test sets (ratio 1:2)
learning.set <- bcw1[learning.set.index,] # we create learning sets
test.set     <- bcw1[-learning.set.index,] # we create test sets
real.labels <- test.set$class # real class labels


## ACCURACY MEASURES
cm.stats <- function(cm) 
{
  # function to construct different accuracy statistics
  tp<-cm[1,1]
  fn<-cm[1,2]
  fp<-cm[2,1]
  tn<-cm[2,2]
  err<-(fp*fn)/(tp+tn+fn+fp)
  tpr<-tp/(tp+fn)
  tnr<-tn/(tn+fp)
  ppv<-tp/(tp+fp)
  fpr<-fp/(tn+fp)
  result<-cbind(err,tpr,tnr,ppv,fpr)
  colnames(result)<-cbind("error rate","sensitivity","specificity","precision","false positive rate")
  return(result)
}

draw_confusion_matrix <- function(cm) {
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, 'Benign', cex=1.2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, 'Malignant', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, 'Benign', cex=1.2, srt=90)
  text(140, 335, 'Malignant', cex=1.2, srt=90)
  
  # confusion matrix results 
  tp<-cm[1,1]
  fn<-cm[1,2]
  fp<-cm[2,1]
  tn<-cm[2,2]
  text(195, 400, tp, cex=1.6, font=2, col='white')
  text(195, 335, fp, cex=1.6, font=2, col='white')
  text(295, 400, fn, cex=1.6, font=2, col='white')
  text(295, 335, tn, cex=1.6, font=2, col='white')
  
  # specifics and accuracy information
  n<-tp+tn+fn+fp
  err<-(fp*fn)/n
  acc<-(tp+tn)/n
  tpr<-tp/(tp+fn)
  tnr<-tn/(tn+fp)
  ppv<-tp/(tp+fp)
  fpr<-fp/(tn+fp)
  miserr <- (n-(tn+tp))/n
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, "error rate", cex=1.2, font=2)
  text(10, 70, round((err), 4), cex=1.2)
  text(30, 85, "sensitivity", cex=1.2, font=2)
  text(30, 70, round(tpr, 4), cex=1.2)
  text(50, 85, "specificity", cex=1.2, font=2)
  text(50, 70, round(tnr, 4), cex=1.2)
  text(70, 85, "precision", cex=1.2, font=2)
  text(70, 70, round(ppv, 4), cex=1.2)
  text(90, 85, "fpr", cex=1.2, font=2)
  text(90, 70, round(fpr, 4), cex=1.2)
  
  text(30, 35, "accuracy", cex=1.2, font=2)
  text(30, 20, round(acc, 4), cex=1.2)
  text(70, 35, "misclassification error", cex=1.2, font=2)
  text(70, 20, round(miserr, 4), cex=1.2)
}  

# ***************************************** 1.0 Linear Regression ***********************************************
# Fitting the vector regression model
X <- cbind(rep(1,n), bcw1[,1:9])
X <- as.matrix(X)
head(X)
Y <- matrix(0, nrow=n, ncol=K) # initialization (fill Y matrix with zeros)
labels.num <- as.numeric(class.labels) #convert labels to numeric
for (k in 1:K)  
  Y[labels.num==k, k] <- 1
B <- solve(t(X)%*%X) %*% t(X) %*% Y
Y.hat <- X%*%B # Predicted values (i.e. predicted probabilities of belonging to particular classes)
rowSums(Y.hat) # Do the predicted probabilities add up to 1?
matplot(Y.hat, main="Predictions (Y.hat)",xlab="id", ylim=c(-.5,2)) #Y.hat plot
abline(h=(0.5), lty=2, col="blue")
legend(x="topright", legend=paste(1:2,levels(bcw1$class)), col=1:2, text.col=1:2, bg="azure2")
classes <- levels(bcw1$class) # Conversion of predicted Y.hat values to class labels
maks.ind <- apply(Y.hat, 1, FUN=function(x) which.max(x)) # For each row we check which element of Y.hat (1 or 2) is the maximum
predicted.labels <- classes[maks.ind] # Conversion to class labels
real.labels <- class.labels
confusion.matrix <- table(real.labels, predicted.labels) #confusion matrix and accuracy
draw_confusion_matrix(confusion.matrix)

pred.ROCR.reg <- prediction(Y.hat[,2], real.labels)
perf.ROCR.reg <- performance(pred.ROCR.reg, "tpr", "fpr")
plot(perf.ROCR.reg)
plot(perf.ROCR.reg,print.cutoffs.at=seq(0.1, 1, 0.1), colorize=TRUE, lwd=2)
lines(c(0,1), c(0,1), lwd=3, lty=2) # add ROC curve for random classifier
grid()
legend("bottomright",lty=c(1,2), lwd=2, col=c("green","black"),
       legend=c("Linear Regression", "random classifier"), bg="azure2", cex=0.7)

AUC.reg  <- performance(pred.ROCR.reg, "auc")@y.values[[1]]
AUC.reg

# ************************ 2.0 LR (Logistic Regression) *********************************** 
model.logit <-  glm(class~., data=bcw1, family=binomial(link="logit"))
summary(model.logit) #model information: coefficients, test values and diagnostics
# Diagnostics include: AIC criterion, deviance (=-2log(likelihood)) and deviance for the reference model (without independent variables)

plot(residuals(model.logit)) # residuals plot
hist(residuals(model.logit)) # residuals histogram
pred.prob <- predict(model.logit, bcw1 , type = "response") # we predict the posterior probability i.e.  Pr('positive'|x)
pp <- ncol(bcw1) # plot (color = class)
colors <- character(n)
colors[bcw1$class=="malignant"] <- "red"
colors[bcw1$class=="benign"] <- "blue"
plot(pred.prob, col=colors)
legend("topright",legend=c("malignant","benign"), col=c("red","blue"), pch="o", bg="azure2")
hist(pred.prob)   # histogram of predicted probabilities
pred.logodds <- predict(model.logit, bcw1, type = "link") # we predict log(odds)
plot(pred.logodds, col=colors) # plot (color = class)
hist(pred.logodds)   # histogram of predicted  log(odds)
prob2labels <- function(probs, cutoff) 
{
  # auxiliary function to convert probabilities to class labels for a given cutoff
  classes <- rep("benign",length(probs))
  classes[probs>cutoff] <- "malignant"
  return(as.factor(classes))
}
lr.pred.labels <- prob2labels(probs=pred.prob, cutoff=0.5) # real and predicted labels
real.labels <- bcw1$class
conf.matrix <- table(lr.pred.labels, real.labels) # confusion matrix
draw_confusion_matrix(conf.matrix)

pred.ROCR.logit <- prediction(pred.prob, real.labels)
perf.ROCR.logit <- performance(pred.ROCR.logit, "tpr", "fpr")
plot(perf.ROCR.logit)
plot(perf.ROCR.logit, print.cutoffs.at=seq(0.1, 1, 0.1), colorize=TRUE, lwd=2) #ROC curve
lines(c(0,1), c(0,1), lwd=3, lty=2) # add ROC curve for random classifier
legend("bottomright",lty=c(1,2), lwd=2, col=c("green","black"),
       legend=c("logistic regression","random classifier"), bg="azure2", cex=0.7)
grid()
title("Logistic ROC curve")
prop.table(table(real.labels)) # compare optimal cut-off with the prior probability
AUC.logit <- performance(pred.ROCR.logit, "auc")@y.values[[1]]
AUC.logit

# ************************ 3.0 k-NN implementation available in the R-package 'ipred' ***********************************
library(caret)
trainSet<-learning.set[,1:9]
trainClass<-learning.set[,10]
cvControl <- trainControl(method="cv", number=10, returnResamp="final", classProbs=TRUE)
'knnFit <- train(trainSet, trainClass, method="knn", tuneLength=10, trControl=cvControl)
kres<-knnFit$results
kres
write.table(kres, file = "Optimal k for kNN.txt", sep = ",", quote = FALSE, row.names = F)
plot(knnFit)'
k.grid <- data.frame(k=1:15)
knnFitGrid <- train(trainSet, trainClass, method="knn", trControl=cvControl, tuneGrid=k.grid)
#write.table(knnFitGrid, file = "Optimal k for kNN.txt", sep = ",", quote = FALSE, row.names = F)
plot(knnFitGrid)
#run the model for 6
library(ipred)
model.knn.6 <- ipredknn(class ~ ., data=learning.set, k=6)
predicted.label6 <- predict(model.knn.6,test.set, type="class") # we predict class labels for test set on the basis of learning set
real.labels <- test.set$class # real class labels
confusion.matrix6 <- table(predicted.label6, real.labels) # confusion matrix and misclassification error
draw_confusion_matrix(confusion.matrix6)
predicted.label.prob6 <- predict(model.knn.6, test.set, type="prob") #we also obtain the proportions
predicted.label.prob6[predicted.label6=="benign"] <- 1 - predicted.label.prob6[predicted.label6=="benign"]
pred.ROCR.knn6 <- prediction(predicted.label.prob6, real.labels)
perf.ROCR.knn6 <- performance(pred.ROCR.knn6, "tpr", "fpr")
plot(perf.ROCR.knn6)
plot(perf.ROCR.knn6,print.cutoffs.at=seq(0.1, 1, 0.1), colorize=TRUE, lwd=2)
lines(c(0,1), c(0,1), lwd=3, lty=2) # add ROC curve for random classifier
grid()
legend("bottomright",lty=c(1,2), lwd=2, col=c("green","black"),
       legend=c("6-NN", "random classifier"), bg="azure2", cex=0.7)

AUC.knn6  <- performance(pred.ROCR.knn6, "auc")@y.values[[1]]
AUC.knn6

# ************************ 4.0 LDA (Linear Discriminant Analysis) ***********************************
#we keep our random split of data: training and test sets
library(MASS)
bcw1.lda  <- lda(class~., data=bcw1, subset=learning.set.index)
print(bcw1.lda)
plot(bcw1.lda)
prediction.lda  <-  predict(bcw1.lda, test.set) # prediction for the test set
pred.prob.lda <- prediction.lda$posterior # predicted posterior probabilities
rowSums(pred.prob.lda)
pred.labels.lda <- prediction.lda$class # predicted class labels
#real.labels <- bcw1$class[-learning.set.index] # real labels for objects from the test set
conf.mat.lda <- table(pred.labels.lda, real.labels) # confusion matrix and misclassification error
draw_confusion_matrix(conf.mat.lda)
real.labels<-test.set$class
pred.ROCR.lda <- prediction(prediction.lda$posterior[,2], real.labels)
perf.ROCR.lda <- performance(pred.ROCR.lda, "tpr", "fpr")
plot(perf.ROCR.lda)
plot(perf.ROCR.lda,print.cutoffs.at=seq(0.1, 1, 0.1), colorize=TRUE, lwd=2)
lines(c(0,1), c(0,1), lwd=3, lty=2) # add ROC curve for random classifier
grid()
legend("bottomright",lty=c(1,2), lwd=2, col=c("green","black"),
       legend=c("LDA", "random classifier"), bg="azure2", cex=0.7)
AUC.lda  <- performance(pred.ROCR.lda, "auc")@y.values[[1]]
AUC.lda

# ************************ 5.0 QDA (Quadratic Discriminant Analysis) ***********************************
bcw.qda  <- qda(class~., data=bcw1, subset=learning.set.index) # construction of the classification rule for all variables
prediction.qda  <-  predict(bcw.qda, test.set) # prediction for the test set
str(prediction.qda)
pred.labels.qda <- prediction.qda$class # predicted class labels
conf.mat.qda <- table(pred.labels.qda, real.labels) #confusion matrix and misclassification error
draw_confusion_matrix(conf.mat.qda)
pred.prob.qda <- prediction.qda$posterior # predicted posterior probabilities
rowSums(pred.prob.qda)
real.labels<-test.set$class
pred.ROCR.qda <- prediction(prediction.qda$posterior[,2], real.labels)
perf.ROCR.qda <- performance(pred.ROCR.qda, "tpr", "fpr")
plot(perf.ROCR.qda)
plot(perf.ROCR.qda,print.cutoffs.at=seq(0.1, 1, 0.1), colorize=TRUE, lwd=2)
lines(c(0,1), c(0,1), lwd=3, lty=2) # add ROC curve for random classifier
grid()
legend("bottomright",lty=c(1,2), lwd=2, col=c("green","black"),
       legend=c("QDA", "random classifier"), bg="azure2", cex=0.7)

AUC.qda  <- performance(pred.ROCR.qda, "auc")@y.values[[1]]
AUC.qda

# ************************ 6.0 Feature selection - stepwise method ***********************************
library(klaR)
lda.forward.selection <- stepclass(class~., data=learning.set, method="lda", direction="forward", improvement=0.01)
qda.forward.selection <- stepclass(class~., data=learning.set, method="qda", direction="forward", improvement=0.01)


#-------------------------------
##FIVE SELECTED FEATURES (IN ORDER OF IMPORTANCE)
#----------------------------------------------------------------------------------------------------------------------------
#make a new dataframe from bcw1
bcw2 <- data.frame(bcw1[,2:3],bcw1[,5:7],bcw1[,10])
colnames(bcw2)<-c("uniformity_cell_size","uniformity_cell_shape","single_epithelial_cell_size",
                  "bare_nuclei","bland_chromatin","class")
View(bcw2)
str(bcw2)

class.labels <-  bcw2$class # class
n <- length(class.labels) # number of objects
K <- length(levels(class.labels)) # number of classes

library(class)
learning.set.index <- sample(1:n,2/3*n) # random split of data into training and test sets (ratio 1:2)
learning.set <- bcw2[learning.set.index,] # we create learning sets
test.set     <- bcw2[-learning.set.index,] # we create test sets
real.labels <- test.set$class # real class labels

# ***************************************** 1.0 Linear Regression ***********************************************
# Fitting the vector regression model
X <- cbind(rep(1,n), bcw2[,1:5])
X <- as.matrix(X)
head(X)
Y <- matrix(0, nrow=n, ncol=K) # initialization (fill Y matrix with zeros)
labels.num <- as.numeric(class.labels) #convert labels to numeric
for (k in 1:K)  
  Y[labels.num==k, k] <- 1
B <- solve(t(X)%*%X) %*% t(X) %*% Y
Y.hat <- X%*%B # Predicted values (i.e. predicted probabilities of belonging to particular classes)
rowSums(Y.hat) # Do the predicted probabilities add up to 1?
matplot(Y.hat, main="Predictions (Y.hat)",xlab="id", ylim=c(-.5,2)) #Y.hat plot
abline(h=(0.5), lty=2, col="blue")
legend(x="topright", legend=paste(1:2,levels(bcw2$class)), col=1:2, text.col=1:2, bg="azure2")
classes <- levels(bcw2$class) # Conversion of predicted Y.hat values to class labels
maks.ind <- apply(Y.hat, 1, FUN=function(x) which.max(x)) # For each row we check which element of Y.hat (1 or 2) is the maximum
predicted.labels <- classes[maks.ind] # Conversion to class labels
real.labels <- class.labels
confusion.matrix <- table(real.labels, predicted.labels) #confusion matrix and accuracy
draw_confusion_matrix(confusion.matrix)

pred.ROCR.reg <- prediction(Y.hat[,2], real.labels)
perf.ROCR.reg <- performance(pred.ROCR.reg, "tpr", "fpr")
plot(perf.ROCR.reg)
plot(perf.ROCR.reg,print.cutoffs.at=seq(0.1, 1, 0.1), colorize=TRUE, lwd=2)
lines(c(0,1), c(0,1), lwd=3, lty=2) # add ROC curve for random classifier
grid()
legend("bottomright",lty=c(1,2), lwd=2, col=c("green","black"),
       legend=c("Linear Regression", "random classifier"), bg="azure2", cex=0.7)

AUC.reg  <- performance(pred.ROCR.reg, "auc")@y.values[[1]]
AUC.reg

# ************************ 2.0 LR (Logistic Regression) *********************************** 
library(mlbench)
library(ROCR)
model.logit <-  glm(class~., data=bcw2, family=binomial(link="logit"))
summary(model.logit) #model information: coefficients, test values and diagnostics
# Diagnostics include: AIC criterion, deviance (=-2log(likelihood)) and deviance for the reference model (without independent variables)

plot(residuals(model.logit)) # residuals plot
hist(residuals(model.logit)) # residuals histogram
pred.prob <- predict(model.logit, bcw2 , type = "response") # we predict the posterior probability i.e.  Pr('positive'|x)
pp <- ncol(bcw2) # plot (color = class)
colors <- character(n)
colors[bcw2$class=="malignant"] <- "red"
colors[bcw2$class=="benign"] <- "blue"
plot(pred.prob, col=colors)
legend("topright",legend=c("malignant","benign"), col=c("red","blue"), pch="o", bg="azure2")
hist(pred.prob)   # histogram of predicted probabilities
pred.logodds <- predict(model.logit, bcw2, type = "link") # we predict log(odds)
plot(pred.logodds, col=colors) # plot (color = class)
hist(pred.logodds)   # histogram of predicted  log(odds)
prob2labels <- function(probs, cutoff) 
{
  # auxiliary function to convert probabilities to class labels for a given cutoff
  classes <- rep("benign",length(probs))
  classes[probs>cutoff] <- "malignant"
  return(as.factor(classes))
}
lr.pred.labels <- prob2labels(probs=pred.prob, cutoff=0.5) # real and predicted labels
real.labels <- bcw2$class
conf.matrix <- table(lr.pred.labels, real.labels) # confusion matrix
draw_confusion_matrix(conf.matrix)

pred.ROCR.logit <- prediction(pred.prob, real.labels)
perf.ROCR.logit <- performance(pred.ROCR.logit, "tpr", "fpr")
plot(perf.ROCR.logit)
plot(perf.ROCR.logit, print.cutoffs.at=seq(0.1, 1, 0.1), colorize=TRUE, lwd=2) #ROC curve
lines(c(0,1), c(0,1), lwd=3, lty=2) # add ROC curve for random classifier
legend("bottomright",lty=c(1,2), lwd=2, col=c("green","black"),
       legend=c("logistic regression","random classifier"), bg="azure2", cex=0.7)
grid()
title("Logistic ROC curve")
prop.table(table(real.labels)) # compare optimal cut-off with the prior probability
AUC.logit <- performance(pred.ROCR.logit, "auc")@y.values[[1]]
AUC.logit

# ************************ 3.0 k-NN implementation available in the R-package 'ipred' ***********************************
library(caret)
trainSet<-learning.set[,1:5]
trainClass<-learning.set[,6]
#run the model for k=6 which we used before
library(ipred)
model.knn.6 <- ipredknn(class ~ ., data=learning.set, k=6)
predicted.label6 <- predict(model.knn.6,test.set, type="class") # we predict class labels for test set on the basis of learning set
real.labels <- test.set$class # real class labels
confusion.matrix6 <- table(predicted.label6, real.labels) # confusion matrix and misclassification error
draw_confusion_matrix(confusion.matrix6)
predicted.label.prob6 <- predict(model.knn.6, test.set, type="prob") #we also obtain the proportions
predicted.label.prob6[predicted.label6=="benign"] <- 1 - predicted.label.prob6[predicted.label6=="benign"]
pred.ROCR.knn6 <- prediction(predicted.label.prob6, real.labels)
perf.ROCR.knn6 <- performance(pred.ROCR.knn6, "tpr", "fpr")
plot(perf.ROCR.knn6)
plot(perf.ROCR.knn6,print.cutoffs.at=seq(0.1, 1, 0.1), colorize=TRUE, lwd=2)
lines(c(0,1), c(0,1), lwd=3, lty=2) # add ROC curve for random classifier
grid()
legend("bottomright",lty=c(1,2), lwd=2, col=c("green","black"),
       legend=c("6-NN", "random classifier"), bg="azure2", cex=0.7)

AUC.knn6  <- performance(pred.ROCR.knn6, "auc")@y.values[[1]]
AUC.knn6

# ************************ 4.0 LDA (Linear Discriminant Analysis) ***********************************
#we keep our random split of data: training and test sets
bcw2.lda  <- lda(class~., data=bcw2, subset=learning.set.index)
print(bcw2.lda)
plot(bcw2.lda)
prediction.lda  <-  predict(bcw2.lda, test.set) # prediction for the test set
pred.prob.lda <- prediction.lda$posterior # predicted posterior probabilities
rowSums(pred.prob.lda)
pred.labels.lda <- prediction.lda$class # predicted class labels
#real.labels <- bcw1$class[-learning.set.index] # real labels for objects from the test set
conf.mat.lda <- table(pred.labels.lda, real.labels) # confusion matrix and misclassification error
draw_confusion_matrix(conf.mat.lda)
real.labels<-test.set$class
pred.ROCR.lda <- prediction(prediction.lda$posterior[,2], real.labels)
perf.ROCR.lda <- performance(pred.ROCR.lda, "tpr", "fpr")
plot(perf.ROCR.lda)
plot(perf.ROCR.lda,print.cutoffs.at=seq(0.1, 1, 0.1), colorize=TRUE, lwd=2)
lines(c(0,1), c(0,1), lwd=3, lty=2) # add ROC curve for random classifier
grid()
legend("bottomright",lty=c(1,2), lwd=2, col=c("green","black"),
       legend=c("LDA", "random classifier"), bg="azure2", cex=0.7)
AUC.lda  <- performance(pred.ROCR.lda, "auc")@y.values[[1]]
AUC.lda

# ************************ 5.0 QDA (Quadratic Discriminant Analysis) ***********************************
bcw.qda  <- qda(class~., data=bcw2, subset=learning.set.index) # construction of the classification rule for all variables
prediction.qda  <-  predict(bcw.qda, test.set) # prediction for the test set
str(prediction.qda)
pred.labels.qda <- prediction.qda$class # predicted class labels
conf.mat.qda <- table(pred.labels.qda, real.labels) #confusion matrix and misclassification error
draw_confusion_matrix(conf.mat.qda)
pred.prob.qda <- prediction.qda$posterior # predicted posterior probabilities
rowSums(pred.prob.qda)
real.labels<-test.set$class
pred.ROCR.qda <- prediction(prediction.qda$posterior[,2], real.labels)
perf.ROCR.qda <- performance(pred.ROCR.qda, "tpr", "fpr")
plot(perf.ROCR.qda)
plot(perf.ROCR.qda,print.cutoffs.at=seq(0.1, 1, 0.1), colorize=TRUE, lwd=2)
lines(c(0,1), c(0,1), lwd=3, lty=2) # add ROC curve for random classifier
grid()
legend("bottomright",lty=c(1,2), lwd=2, col=c("green","black"),
       legend=c("QDA", "random classifier"), bg="azure2", cex=0.7)

AUC.qda  <- performance(pred.ROCR.qda, "auc")@y.values[[1]]
AUC.qda

# ************************ 6.0 Feature selection - stepwise method ***********************************
library(klaR)
lda.forward.selection <- stepclass(class~., data=learning.set, method="lda", direction="forward", improvement=0.01)
qda.forward.selection <- stepclass(class~., data=learning.set, method="qda", direction="forward", improvement=0.01)








rname <- list()
for (i in 1:length(bcw.real.class.labels)){
  b <- paste0("obs_",toString(i))
  rname <- append(rname, b)
}
rownames(bcw) <- bcw.names
selectedFeatures <- setdiff(colnames(bcw), c("id_number","class"))
bcw.selected  <- bcw[ ,selectedFeatures]

bcw.DissimilarityMatrix <- daisy(bcw.features)
# Conversion to matrix
bcw.DissimilarityMatrix.mat <- as.matrix(bcw.DissimilarityMatrix)

fviz_dist(bcw.DissimilarityMatrix.mat, order = FALSE) # without ordering
fviz_dist(bcw.DissimilarityMatrix.mat, order = TRUE) # after ordering







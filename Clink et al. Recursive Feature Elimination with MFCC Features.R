##### R code for Clink et al.

##### Application of a semi-automated vocal fingerprinting approach to monitor Bornean
##### gibbon females in an experimentally fragmented landscape in Sabah, Malaysia

#####Part 2 Recursive Feature Elimination ####

####Set source code (downloaded from https://github.com/johncolby/SVM-RFE)
source("/Users/denasmacbook/Downloads/SVM-RFE-master/msvmRFE.R")

####Read in MFCC data with standardized number of time window
mfcc.standard.windows <-
  read.csv(
    "/Users/denasmacbook/Downloads/mfcc.with.delta.8.windows.without.1.updated.Nov.12.csv"
  )

### For SVM RFE he first column needs to include class labels; remove other id columns
for.svm.rfe <-
  subset(mfcc.standard.windows, select = -c(X, filehandles, call.id))

###We want k=10 for the k-fold cross validation as the “multiple” part of mSVM-RFE.
svm.rfe.output <- svmRFE(for.svm.rfe, k = 10, halve.above = 100)
str(svm.rfe.output)

###Reorder the data so highest ranked feature is first
new.svm.rfe <-
  for.svm.rfe[, 2:ncol(for.svm.rfe)][, dput(svm.rfe.output)]
str(new.svm.rfe)

###Create a list to store cross-validation accuracies
accuracy.list <- list()

n.females <- length(unique(mfcc.standard.windows$female.id))

####Loop to add features one by one and calculate accuracy using leave-one-out cross-validation
for (j in 2:length(svm.rfe.output)) {
  svm.rfe.for.lda <- new.svm.rfe[1:j]
  
  fit.svm.rfe <- lda(
    svm.rfe.for.lda,
    center = TRUE,
    prior = rep(1 / n.females, n.females),
    
    scale. = TRUE,
    grouping = mfcc.standard.windows$female.id,
    CV = T
  )
  
  ####Assess how well the leave one out cross validation did
  ct <- table(grouping = for.svm.rfe$female.id, fit.svm.rfe$class)
  
  # total percent correct
  percent <- sum(diag(prop.table(ct)))
  print(percent)
  accuracy.list[[j]] <- percent
}

###Create a figure
plot(unlist(accuracy.list), xlab = "Number of Features", ylab = "Classification Accuracy")

###Find which number of features provides the maximum classification accuracy
max.feature <- which.max(unlist(accuracy.list)) + 1

###Subset the highest ranked variables which yield the highest accuracy
svm.rfe.for.lda <- new.svm.rfe[1:max.feature]

### Combine class labels with new subset of features into a data frame for analysis
svm.rfe.for.classification <-
  cbind.data.frame(mfcc.standard.windows$female.id, svm.rfe.for.lda)
colnames(svm.rfe.for.classification)[1] <- "female.id"
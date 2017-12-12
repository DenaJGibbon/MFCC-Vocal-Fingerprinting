##### R code for Clink et al.

##### Application of a semi-automated vocal fingerprinting approach to monitor Bornean
##### gibbon females in an experimentally fragmented landscape in Sabah, Malaysia

#####Part 4 Random iterations of SVM ####
## NOTE: To save time the code only runs over 10 iterations

##### Load required libraries
library(MASS)
library(e1071)

##### Load Data

mfcc.data.all <-
  read.csv("/Users/denasmacbook/Desktop/Clink_et_at_MFCC/MFCC.data.csv")
str(mfcc.data.all)

### Create empty list
list.dens.linear = list()

### Loop to randomly choose calls as test and training set for linear kernel
for (j in 1:10) {
  mfcc.data.all$ran.num <- runif(376, 0, 1)
  mfcc.train <- mfcc.data.all[mfcc.data.all$ran.num < .8, ]
  mfcc.test <- mfcc.data.all[mfcc.data.all$ran.num >= .8, ]
  
  linear.tune <-
    tune(
      svm,
      mfcc.train[, 4:25],
      mfcc.train$female.id,
      kernel = "linear",
      tunecontrol = tune.control(cross = 5),
      ranges = list(
        cost = c(0.001, 0.01, 0.1, 1, 2, 10, 10, 100),
        gamma = c(0.01, 0.1, 0.5, 1.0, 2.0)
      )
    )
  
  svm.model <-
    svm(
      mfcc.train[, 4:25],
      mfcc.train$female.id,
      kernel = "linear",
      cost = linear.tune$best.parameters$cost,
      cross = 5
    )
  
  mod.pred <- predict(svm.model, mfcc.test[, 4:25])
  confusion.mat <-
    table(pred = mod.pred, true = mfcc.test$female.id)
  print(sum(diag(prop.table(confusion.mat))))
  new.percent <- sum(diag(prop.table(confusion.mat)))
  list.dens.linear[[j]] <- new.percent
}



### Create empty list
list.dens.poly = list()

### Loop to randomly choose calls as test and training set for polynomial kernel
for (j in 1:10) {
  mfcc.data.all$ran.num <- runif(376, 0, 1)
  mfcc.train <- mfcc.data.all[mfcc.data.all$ran.num < .8, ]
  mfcc.test <- mfcc.data.all[mfcc.data.all$ran.num >= .8, ]
  
  tune.poly <-
    tune(
      svm,
      mfcc.train[, 4:25],
      mfcc.train$female.id,
      kernel = "polynomial",
      tunecontrol = tune.control(cross = 5),
      ranges = list(
        cost = c(0.001, 0.01, 0.1, 1, 2, 10, 10, 100),
        gamma = c(0.01, 0.1, 0.5, 1.0, 2.0)
      )
    )
  
  svm.model <-
    svm(
      mfcc.train[, 4:25],
      mfcc.train$female.id,
      kernel = "polynomial",
      gamma = tune.poly$best.parameters$gamma,
      cost = tune.poly$best.parameters$cost,
      cross = 5
    )
  
  mod.pred <- predict(svm.model, mfcc.test[, 4:25])
  confusion.mat <-
    table(pred = mod.pred, true = mfcc.test$female.id)
  print(sum(diag(prop.table(confusion.mat))))
  new.percent <- sum(diag(prop.table(confusion.mat)))
  list.dens.poly[[j]] <- new.percent
}


### Create empty list
list.dens.sig = list()

### Loop to randomly choose calls as test and training set for sigomidal kernel
for (j in 1:10) {
  mfcc.data.all$ran.num <- runif(376, 0, 1)
  mfcc.train <- mfcc.data.all[mfcc.data.all$ran.num < .8, ]
  mfcc.test <- mfcc.data.all[mfcc.data.all$ran.num >= .8, ]
  
  tune.sig <-
    tune(
      svm,
      mfcc.train[, 4:25],
      mfcc.train$female.id,
      kernel = "sigmoid",
      tunecontrol = tune.control(cross = 5),
      ranges = list(
        cost = c(0.001, 0.01, 0.1, 1, 2, 10, 10, 100),
        gamma = c(0.01, 0.1, 0.5, 1.0, 2.0)
      )
    )
  
  svm.model <-
    svm(
      mfcc.train[, 4:25],
      mfcc.train$female.id,
      kernel = "sigmoid",
      gamma = tune.sig$best.parameters$gamma,
      cost = tune.sig$best.parameters$cost,
      cross = 5
    )
  
  mod.pred <- predict(svm.model, mfcc.test[, 4:25])
  confusion.mat <-
    table(pred = mod.pred, true = mfcc.test$female.id)
  
  print(sum(diag(prop.table(confusion.mat))))
  
  new.percent <- sum(diag(prop.table(confusion.mat)))
  list.dens.sig[[j]] <- new.percent
}



list.dens.rad <- list()

### Loop to randomly choose calls as test and training set for sigomidal kernel
for (j in 1:10) {
  mfcc.data.all$ran.num <- runif(376, 0, 1)
  mfcc.train <- mfcc.data.all[mfcc.data.all$ran.num < .8, ]
  mfcc.test <- mfcc.data.all[mfcc.data.all$ran.num >= .8, ]
  
  tune.rad <-
    tune(
      svm,
      mfcc.train[, 4:25],
      mfcc.train$female.id,
      kernel = "radial",
      tunecontrol = tune.control(cross = 5),
      ranges = list(
        cost = c(0.001, 0.01, 0.1, 1, 2, 10, 10, 100),
        gamma = c(0.01, 0.1, 0.5, 1.0, 2.0)
      )
    )
  
  svm.model <-
    svm(
      mfcc.train[, 4:25],
      mfcc.train$female.id,
      kernel = "radial",
      gamma = tune.rad$best.parameters$gamma,
      cost = tune.rad$best.parameters$cost,
      cross = 5
    )
  
  mod.pred <- predict(svm.model, mfcc.test[, 4:25])
  confusion.mat <-
    table(pred = mod.pred, true = mfcc.test$female.id)
  print(sum(diag(prop.table(confusion.mat))))
  new.percent <- sum(diag(prop.table(confusion.mat)))
  list.dens.rad[[j]] <- new.percent
}

plot(density(unlist(list.dens.rad)))

### Prep data to go to dataframe
mfcc.den.lin <-
  cbind.data.frame(unlist(list.dens.linear), rep("Linear", 10))
mfcc.den.rad <-
  cbind.data.frame(unlist(list.dens.rad), rep("Radial", 10))
mfcc.den.poly <-
  cbind.data.frame(unlist(list.dens.poly), rep("Polynomial", 10))
mfcc.den.sig <-
  cbind.data.frame(unlist(list.dens.sig), rep("Sigmoidal", 10))

### Add column names
colnames(mfcc.den.lin) <- c("Accuracy", "Kernel")
colnames(mfcc.den.rad) <- c("Accuracy", "Kernel")
colnames(mfcc.den.poly) <- c("Accuracy", "Kernel")
colnames(mfcc.den.sig) <- c("Accuracy", "Kernel")

### Create new dataframe for density plot
mfcc.den.df <-
  rbind.data.frame(mfcc.den.lin, mfcc.den.rad, mfcc.den.poly, mfcc.den.sig)
mfcc.den.df$Accuracy <- mfcc.den.df$Accuracy * 10

#### R code for Clink et al. 
#### Application of a semi-automated vocal fingerprinting approach to monitor Bornean
#### gibbon females in an experimentally fragmented landscape in Sabah, Malaysia 

#### Part 1 Two methods of MFCC feature extraction: averaging across all time windows  
#### and creating a standardized number of time windows for each call

#### Part 2 Recursive Feature Elimination

#### Part 3 LDA vs. SVM

#### Part 4 Random iterations of SVM


# Load required libraries
library(tuneR)
library(seewave)
library(MASS)
library(e1071)
library(stringr)
library(ggplot2)
library(viridis)



#### Part 1 MFCC Estimation ####

### Save all .wav files and set working directory to access the files 
### There are nine .wav files, three from three different females (SAFA, SAFBA and VJRN_01)
setwd("/Users/denasmacbook/Downloads/MFCC-Vocal-Fingerprinting-master")

### Set the input directory to loop over the wave files
input.dir <-"/Users/denasmacbook/Downloads/MFCC-Vocal-Fingerprinting-master"

### List all .wav files in directory
L = list.files(input.dir, pattern="*.wav", full.names=FALSE)
L

### Extract file names 
filehandles <- str_split_fixed(L, pattern=".wav", n=2)[,1]

### MFCC Feature Extraction Method 1: Averaging over time windows ###

### Create empty list to hold MFCC values
mfcc.vector.list = list()

### Loop to calculate MFCC for each .wav file in the directory
for (j in 1:length(filehandles)) { 
  filehandle <-  L[j]
  filename <- paste(input.dir, "/", filehandle, sep="")
  print(paste("processing", filehandle))
  
  # Read in wav file 
  w<- readWave(filename)
  
  # Calculate 12 MFCCs for each 0.25 ms window
  melfcc.output <- melfcc(w, minfreq=400, maxfreq=2000, wintime = 0.25,fbtype="mel",numcep=12)
  
  #Add MFCC values to list
  mfcc.vector.list[[j]] <- melfcc.output
}


### Check structure of MFCC list; for each call there is matrix with 12 MFCC for each 0.25 s frame index
### Calls of different length will have different number of frames 
str(mfcc.vector.list)


### Create an empty list to store the mean and sd for each mel-frequency bin
mean.sd.list = list()

### Loop to calculate MFCC mean and sd 
for(j in 1:length(mfcc.vector.list)) {
  tmp.list <- mfcc.vector.list[[j]]
  list.elements<- lapply(1:ncol(tmp.list), function(x) {
    vec.temp <- tmp.list[,x]
    vec.mean <- mean(vec.temp)
    vec.sd <- sd(vec.temp)
    data.vec <- c(vec.mean,vec.sd)
    data.vec
  })
  mean.sd.list[[j]] <- list.elements
}

### Convert the list to a vector 
vec <- unlist(mean.sd.list)

### Convert the vector into a matrix
data.matrix.all <- matrix(vec,nrow = length(filehandles), ncol = 24, byrow=T)
data.matrix.all <- as.data.frame(data.matrix.all)

### Split file names into female id and call id 
loc <- str_split_fixed(filehandles, pattern="_", n=2)[,1]
group <- str_split_fixed(filehandles, pattern="_", n=3)[,2]
female.id <- paste(loc,group, sep="_")
call.id <- str_split_fixed(filehandles, pattern="[.]", n=3)[,1]

### Create unique column names for the mean and sd of each Mel-frequency bin
colnames <- rep(c("mean","sd"),12)
nums <- rep(seq(1,12),2)
nums <- sort(nums)
colnames <- paste(colnames,nums)
colnames(data.matrix.all) <- colnames

### Combine into a dataframe for analysis
mfcc.data.all <- cbind.data.frame(data.matrix.all,female.id,call.id,filehandles)

### Check structure of resulting dataframe
str(mfcc.data.all)


### MFCC Feature Extraction Method 2: Standardized number of time windows ###

mfcc.vector.list = list()

for (j in 1:length(filehandles)) { tryCatch({
  filehandle <-  L[j]
  filename <- paste(input.dir, "/", filehandle, sep="")
  print(paste("processing", filehandle))
  
  # Read in wav file and filter
  w<- readWave(filename)
  
  # Find duration of .wav file and divide into 9 windows
  wav.dur <- duration(w)
  win.time <- wav.dur/9
  
  # Calculate MFCCs
  melfcc.output <- melfcc(w, minfreq=400, hoptime=win.time, 
                          maxfreq=2000,numcep=12, wintime = win.time)
  
  # Calculate delta cepstral coefficients
  deltas.output <- deltas(melfcc.output)
  
  # Ensure only 8 time windows are used for MFCC and delta coefficients
  # Also append .wav duration
  mfcc.vector <- c(as.vector(t(melfcc.output[1:8,2:12])),as.vector(t(deltas.output[1:8,2:12])),wav.dur)
  
  # Add to list
  mfcc.vector.list[[j]] <- mfcc.vector
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

### Convert the list to a vector 
vec <- unlist(mfcc.vector.list)

### Split file names into female id and call id 
loc <- str_split_fixed(filehandles, pattern="_", n=2)[,1]
group <- str_split_fixed(filehandles, pattern="_", n=3)[,2]
female.id <- paste(loc,group, sep="_")
call.id <- str_split_fixed(filehandles, pattern="[.]", n=3)[,1]

### Create a matrix with MFCC and delta coefficient values
### The number of columns is a reflection of the number of time windows
data.matrix.all <- matrix(vec,nrow = length(filehandles), ncol = length(mfcc.vector.list[[1]]), byrow=T)
mfcc.data.frame <- cbind.data.frame(female.id,filehandles,call.id,data.matrix.all)
seq.length <- (length(mfcc.vector.list[[1]])-1)/2
col.names <- c("female.id","filehandles","call.id",paste(rep("mfcc",seq.length), seq(seq.length), sep="_"),paste(rep("delta",seq.length), seq(seq.length), sep="_"),"dur")
colnames(mfcc.data.frame)<- col.names
mfcc.data.frame <- as.data.frame(mfcc.data.frame)

### Check data structure 
str(mfcc.data.frame)

#### Part 2 Recursive Feature Elimination #### 

### Set source code (downloaded from https://github.com/johncolby/SVM-RFE)
source("/Users/denasmacbook/Downloads/SVM-RFE-master/msvmRFE.R")

### Read in MFCC data with standardized number of time window
mfcc.data.frame <- read.csv("/Users/denasmacbook/Downloads/mfcc.with.delta.8.windows.without.1.updated.Nov.12.csv")

## The first column needs to include class labels; remove other id columns
for.svm.rfe <- subset(mfcc.data.frame, select = -c(X, filehandles,call.id))

## We want k=10 for the k-fold cross validation as the “multiple” part of mSVM-RFE.
svm.rfe.output <- svmRFE(for.svm.rfe, k=10, halve.above=100)
str(svm.rfe.output)

## Reorder the data so highest ranked feature is first
new.svm.rfe <- for.svm.rfe[,2:ncol(for.svm.rfe)][,dput(svm.rfe.output)]

## Create a list to store cross-validation accuracies
accuracy.list <- list()

### Loop to iteratively add features and calculate accuracy using leave-one-out cross-validation
for (j in 2:length(svm.rfe.output)) { 
  
  svm.rfe.for.lda <- new.svm.rfe[1:j]
  
  fit.svm.rfe <- lda(svm.rfe.for.lda,
                     center = TRUE,
                     prior= rep(1/n.females,n.females), ## want to report
                     scale. = TRUE, grouping=mfcc.data.frame$female.id, CV=T) 
  
  ### Assess how well the leave one out cross validation did
  ct <- table(grouping=for.svm.rfe$female.id, fit.svm.rfe$class)
  
  # total percent correct
  percent <- sum(diag(prop.table(ct))) 
  print(percent)
  accuracy.list[[j]] <-percent
}

## Create a figure 
plot(unlist(accuracy.list), xlab="Number of Features", ylab="Classification Accuracy")

## Find which number of features provides the maximum classification accuracy
max.feature <- which.max(unlist(accuracy.list))+1


svm.rfe.for.lda <- new.svm.rfe[1:max.feature]

fit.svm.rfe <- lda(svm.rfe.for.lda,
                   center = TRUE,
                   prior= rep(1/n.females,n.females), ## want to report
                   scale. = TRUE, grouping=mfcc.data.frame$female.id, CV=T) 

### Assess how well the leave one out cross validation did
ct <- table(grouping=for.svm.rfe$female.id, fit.svm.rfe$class)

# total percent correct
percent <- sum(diag(prop.table(ct))) 
print(percent)


#### Part 3 LDA vs. SVM  ####

### Read in MFCC features averaged across all time windows
mfcc.data.all <-read.csv("/Users/denasmacbook/Desktop/Clink_et_at_MFCC/MFCC.data.csv")
#mfcc.data.all <- read.csv("/Users/denasmacbook/Desktop/Clink_et_at_MFCC/spectral_features_SAFE.data.csv")

### Check structure
str(mfcc.data.all)

### Run LDA
n.females <- length(unique(mfcc.data.all$female.id))

fit <- lda(mfcc.data.all[4:25],
           center = TRUE,
           prior= rep(1/n.females,n.females), ## want to report
           scale. = TRUE, grouping=mfcc.data.all$female.id, CV=T) 

### Assess how well the leave one out cross validation did
ct <- table(grouping=mfcc.data.all$female.id, fit$class)

### Calculate the total percent correct
sum(diag(prop.table(ct))) # 98.4 % accuracy!


### Tune the parameters and run the SVM with radial kernel

tune.radial <- tune(svm,mfcc.data.all[7:36], mfcc.data.all$female.id,kernel= "sigmoid",
                    ranges=list(cost=c(0.001,0.01,0.1,1,5,10,100),gamma=c(0.001,0.01,0.1,0.5,1.0,2.0)))

cost.radial <- tune.radial$best.parameters$cost
gamma.radial <-  tune.radial$best.parameters$gamma

# When cross = number of observations this indicates leave-one-out cross-validation
svm.radial <- svm(mfcc.data.all[7:36], mfcc.data.all$female.id, kernel= "sigmoid", cost=cost.radial, gamma=gamma.radial, cross=nrow(mfcc.data.all))
svm.radial$tot.accuracy # 99.2 % accuracy


### Tune the parameters and run the SVM with linear kernel

tune.linear <- tune(svm,mfcc.data.all[7:36], mfcc.data.all$female.id,kernel= "linear",
                    ranges=list(cost=c(0.001,0.01,0.1,1,5,10,100)))

cost.linear <- tune.linear$best.parameters$cost

# When cross = number of observations this indicates leave-one-out cross-validation
svm.linear <- svm(mfcc.data.all[7:36], mfcc.data.all$female.id, kernel= "linear", cost=cost.linear,k=nrow(mfcc.data.all),cross=376)
svm.linear$tot.accuracy # 99.5% accuracy



#### Part 3 Random iterations of SVM ####
### NOTE: The code below only does 10 iterations to save computing time

mfcc.data.all <- read.csv("/Users/denasmacbook/Downloads/MFCC-Vocal-Fingerprinting-master/MFCC.data.csv")
mfcc.data.all <- mfcc.data.all
### Create empty list
list.dens.linear =list()

### Loop to randomly choose calls as test and training set for linear kernel
for (j in 1:100) { 
  mfcc.data.all$ran.num <- runif(376, 0,1)
  mfcc.train <- mfcc.data.all[mfcc.data.all$ran.num < .8,]
  mfcc.test <- mfcc.data.all[mfcc.data.all$ran.num >= .8,]
  
  linear.tune<- tune(svm,mfcc.train[1:24], mfcc.train$female.id, kernel= "linear", tunecontrol=tune.control(cross=5),
                              ranges=list(cost=c(0.001, 0.01, 0.1, 1, 2, 10, 100, 1000),gamma=c(0.01,0.1,0.5,1.0,2.0)))
  
  svm.model <- svm(mfcc.train[1:24], mfcc.train$female.id,kernel= "linear", 
                   cost=linear.tune$best.parameters$cost,k=nrow(mfcc.data.all))
  
  mod.pred <- predict(svm.model, mfcc.test[1:24])
  confusion.mat <- table(pred = mod.pred, true = mfcc.test$female.id)
  print(sum(diag(prop.table(confusion.mat))))
  new.percent <- sum(diag(prop.table(confusion.mat)))
  list.dens.linear[[j]] <- new.percent
}


### Create empty list
list.dens.radial =list()

### Loop to randomly choose calls as test and training set for radial kernel
for (j in 1:10) { 
  mfcc.data.all$ran.num <- runif(376, 0,1)
  mfcc.train <- mfcc.data.all[mfcc.data.all$ran.num < .4,]
  mfcc.test <- mfcc.data.all[mfcc.data.all$ran.num >= .4,]
  
  tune.radial <- tune(svm,mfcc.data.all[1:24], mfcc.data.all$female.id, kernel= "radial", 
                      ranges=list(cost=c(0.001, 0.01, 0.1, 1, 2, 10, 100, 1000),gamma=c(0.01,0.1,0.5,1.0,2.0)))
  svm.model <- svm(mfcc.train[1:24], mfcc.train$female.id, kernel= "radial",gamma=tune.radial$best.parameters$gamma,cost=tune.radial$best.parameters$cost,
                   k=length(unique(mfcc.data.all.train$female.id.train)),cross=10)
  
  mod.pred <- predict(svm.model, mfcc.test[1:24])
  confusion.mat <- table(pred = mod.pred, true = mfcc.test$female.id)
  diag(prop.table(confusion.mat, 1))
  new.percent <- sum(diag(prop.table(confusion.mat)))
  new.percent
  list.dens.radial[[j]] <- new.percent
}

### Create empty list
list.dens.poly =list()

### Loop to randomly choose calls as test and training set for polynomial kernel
for (j in 1:10) { 
  mfcc.data.all$ran.num <- runif(376, 0,1)
  mfcc.train <- mfcc.data.all[mfcc.data.all$ran.num < .4,]
  mfcc.test <- mfcc.data.all[mfcc.data.all$ran.num >= .4,]
  
  tune.poly <- tune(svm,mfcc.data.all[1:24], mfcc.data.all$female.id, kernel= "polynomial", 
                    ranges=list(cost=c(0.001, 0.01, 0.1, 1, 2, 10, 100, 1000),gamma=c(0.01,0.1,0.5,1.0,2.0)))
  
  svm.model <- svm(mfcc.train[1:24], mfcc.train$female.id, kernel= "polynomial",gamma=tune.poly$best.parameters$gamma,
                   cost=tune.poly$best.parameters$cost,k=length(unique(mfcc.data.all.train$female.id.train)),cross=10)
  
  mod.pred <- predict(svm.model, mfcc.test[1:24])
  confusion.mat <- table(pred = mod.pred, true = mfcc.test$female.id)
  diag(prop.table(confusion.mat, 1))
  new.percent <- sum(diag(prop.table(confusion.mat)))
  list.dens.poly[[j]] <- new.percent
}

### Create enpty list
list.dens.sig =list()

### Loop to randomly choose calls as test and training set for sigomidal kernel
for (j in 1:10) { 
  mfcc.data.all$ran.num <- runif(376, 0,1)
  mfcc.train <- mfcc.data.all[mfcc.data.all$ran.num < .8,]
  mfcc.test <- mfcc.data.all[mfcc.data.all$ran.num >= .8,]

  tune.sig <- tune(svm,mfcc.data.all[1:24], mfcc.data.all$female.id, kernel= "sigmoid", 
                   ranges=list(cost=c(0.001, 0.01, 0.1, 1, 2, 10, 100, 1000),gamma=c(0.01,0.1,0.5,1.0,2.0)))
  
  svm.model <- svm(mfcc.train[1:24], mfcc.train$female.id, kernel= "sigmoid",gamma=tune.sig$best.parameters$gamma,
                   cost=tune.sig$best.parameters$cost,k=length(unique(mfcc.data.all.train$female.id.train)),cross=10)
  mod.pred <- predict(svm.model, mfcc.test[1:24])
  confusion.mat <- table(pred = mod.pred, true = mfcc.test$female.id)
  diag(prop.table(confusion.mat, 1))
  new.percent <- sum(diag(prop.table(confusion.mat)))
  list.dens.sig[[j]] <- new.percent
}


### Prep data to go to dataframe
mfcc.den.lin <- cbind.data.frame(unlist(list.dens.linear), rep("Linear",10))
mfcc.den.rad <- cbind.data.frame(unlist(list.dens.radial),rep("Radial",10))
mfcc.den.poly <- cbind.data.frame(unlist(list.dens.poly),rep("Polynomial",10))
mfcc.den.sig <- cbind.data.frame(unlist(list.dens.sig),rep("Sigmoidal",10))

### Add column names
colnames(mfcc.den.lin) <- c("Accuracy", "Kernel")
colnames(mfcc.den.rad) <- c("Accuracy", "Kernel")
colnames(mfcc.den.poly) <- c("Accuracy", "Kernel")
colnames(mfcc.den.sig) <- c("Accuracy", "Kernel")

### Create new dataframe for density plot
mfcc.den.df <- rbind.data.frame(mfcc.den.lin,mfcc.den.rad,mfcc.den.poly,mfcc.den.sig)
mfcc.den.df$Accuracy <- mfcc.den.df$Accuracy*100

### Plot results
ggplot(mfcc.den.df, aes(x=Accuracy, fill=Kernel)) + geom_density(alpha=.5,bw=1)+
  scale_y_continuous(breaks=NULL)+
  scale_fill_viridis(option="D", discrete=T)+
  guides(fill=guide_legend(title="Kernel Type"))+xlab("Percent Correct Female Identification")+ylab("Count of 1000 iterations")+
  theme_classic()+
  theme(axis.text.y=element_text(size=20),legend.text=element_text(size=18),legend.title = element_text(size=20),
        axis.text.x  = element_text(size=18),axis.title=element_text(size=20,face="bold"))



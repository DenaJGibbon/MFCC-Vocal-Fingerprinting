#### R code for Clink et al. 
#### Part 1 MFCC Estimation
#### Part 2 LDA vs. SVM
#### Part 3 Random iterations of SVM


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
setwd("/Users/denajaneclink/Desktop/MFCC_examples")

### Set the input directory to loop over the wave files
input.dir <-"/Users/denajaneclink/Desktop/MFCC_examples"

### List all .wav files in directory
L = list.files(input.dir, pattern="*.wav", full.names=FALSE)
L

### Extract file names 
filehandles <- str_split_fixed(L, pattern=".wav", n=2)[,1]

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


### Check structure of MFCC list; for each call there is matrix with 12 MFCC for each 0.25 ms frame index
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


#### Part 2 LDA vs. SVM  ####

### Read in MFCC for all females in the dataset
mfcc.data.all <-read.csv("/Users/denajaneclink/Desktop/MFCC_examples/Data/MFCC.data.csv")

### Check structure
str(mfcc.data.all)

hist(mfcc.data.all$V5)

### Run LDA
n.females <- length(unique(mfcc.data.all$female.id))

fit <- lda(mfcc.data.all[3:24],
           center = TRUE,
           prior= rep(1/n.females,n.females), ## want to report
           scale. = TRUE, grouping=mfcc.data.all$female.id)#, CV=T) 

### Assess how well the leave one out cross validation did
ct <- table(grouping=mfcc.data.all$female.id, fit$class)

### Calculate the total percent correct
sum(diag(prop.table(ct))) # 99.7 % accuracy!

ggord(fit,mfcc.data.all$female.id)+ 
  scale_colour_manual(mfcc.data.all$female.id, values = viridis(n.females))+
   theme(legend.position="none")


### Tune the parameters and run the SVM with radial kernel

tune.radial <- tune(svm,mfcc.data.all[1:24], mfcc.data.all$female.id,kernel= "radial",
                    ranges=list(cost=c(0.001,0.01,0.1,1,5,10,100),gamma=c(0.001,0.01,0.1,0.5,1.0,2.0)))

cost.radial <- tune.radial$best.parameters$cost
gamma.radial <-  tune.radial$best.parameters$gamma

# When cross = number of observations this indicates leave-one-out cross-validation
svm.radial <- svm(mfcc.data.all[1:24], mfcc.data.all$female.id, kernel= "radial", cost=cost.radial, gamma=gamma.radial, cross=376)
svm.radial$tot.accuracy # 99.2 % accuracy


### Tune the parameters and run the SVM with linear kernel

tune.linear <- tune(svm,mfcc.data.all[3:24], mfcc.data.all$female.id,kernel= "linear",
                    ranges=list(cost=c(0.001,0.01,0.1,1,5,10,100)))

cost.linear <- tune.linear$best.parameters$cost

# When cross = number of observations this indicates leave-one-out cross-validation
svm.linear <- svm(mfcc.data.all[3:24], mfcc.data.all$female.id, kernel= "linear", cost=cost.linear,k=nrow(mfcc.data.all),cross=376)
svm.linear$tot.accuracy # 99.5% accuracy


#### Part 3 Random iterations of SVM ####
### NOTE: The code below only does 10 iterations to save computing time


### Create empty list
list.dens.linear =list()

### Loop to randomly choose calls as test and training set for linear kernel
for (j in 1:10) { 
  mfcc.data.all$ran.num <- runif(376, 0,1)
  mfcc.train <- mfcc.data.all[mfcc.data.all$ran.num < .4,]
  mfcc.test <- mfcc.data.all[mfcc.data.all$ran.num >= .4,]
  
  linear.tune <- tune(svm,mfcc.data.all[1:24], mfcc.data.all$female.id,k=nrow(mfcc.data.all), 
                      kernel="linear",ranges=list(cost=c(0.001, 0.01, 0.1, 1, 2, 10, 100, 1000)))
  
  svm.model <- svm(mfcc.train[1:24], mfcc.train$female.id,kernel= "linear", 
                   cost=linear.tune$best.parameters$cost,k=nrow(mfcc.data.all),cross=10)
  
  mod.pred <- predict(svm.model, mfcc.test[1:24])
  confusion.mat <- table(pred = mod.pred, true = mfcc.test$female.id)
  diag(prop.table(confusion.mat, 1))
  new.percent <- sum(diag(prop.table(confusion.mat)))
  new.percent
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
  mfcc.train <- mfcc.data.all[mfcc.data.all$ran.num < .4,]
  
  mfcc.test <- mfcc.data.all[mfcc.data.all$ran.num >= .4,]
  nrow(mfcc.train)
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



##### R code for Clink et al.

##### Application of a semi-automated vocal fingerprinting approach to monitor Bornean
##### gibbon females in an experimentally fragmented landscape in Sabah, Malaysia

##### Part 3 LDA

##### Load required libraries
library(MASS)

####Read in MFCC features chosen via recursive feature elimination
svm.rfe.for.classification <-
  read.csv("/Users/denasmacbook/Desktop/Clink_et_at_MFCC/svm.rfe.for.classification.csv")

## Check structure
str(svm.rfe.for.classification)

## Set number of females for prior
n.females <- length(unique(svm.rfe.for.classification$female.id))

###Run LDA on the data subset using RFE
fit.standard.number.windows.svm.rfe <- lda(
  svm.rfe.for.classification[2:ncol(svm.rfe.for.classification)],
  center = TRUE,
  prior = rep(1 / n.females, n.females),
  scale. = TRUE,
  grouping = svm.rfe.for.classification$female.id,
  CV = T
)

####Assess how well the leave one out cross validation did
ct <-
  table(grouping = svm.rfe.for.classification$female.id, fit.standard.number.windows.svm.rfe$class)

####Calculate total percent correct
percent <- sum(diag(prop.table(ct)))
print(percent) # 98.9 % accuracy!

####Read in MFCC features averaged across all time windows
mfcc.average.over.windows <-
  read.csv("/Users/denasmacbook/Desktop/Clink_et_at_MFCC/MFCC.data.csv")

####Check structure
str(mfcc.average.over.windows)

## Set number of females for prior
n.females <- length(unique(mfcc.average.over.windows$female.id))

## Run LDA
fit.average.over.windows <- lda(
  mfcc.average.over.windows[4:25],
  center = TRUE,
  prior = rep(1 / n.females, n.females),
  
  scale. = TRUE,
  grouping = mfcc.average.over.windows$female.id,
  CV = T
)

####Assess how well the leave one out cross validation did
ct <-
  table(grouping = mfcc.average.over.windows$female.id, fit.average.over.windows$class)

####Calculate the total percent correct
sum(diag(prop.table(ct))) # 98.4 % accuracy!


####Read in spectral features estimated from spectrogram (data from Clink et al. 2017)
spectral.features <-
  read.csv("/Users/denasmacbook/Desktop/Clink_et_at_MFCC/spectral_features_SAFE.data.csv")

####Check structure
str(spectral.features)

## Set number of females for prior
n.females <- length(unique(spectral.features$female.id))

####Run LDA
fit <- lda(
  spectral.features[7:36],
  center = TRUE,
  prior = rep(1 / n.females, n.females),
  
  scale. = TRUE,
  grouping = spectral.features$female.id,
  CV = T
)

####Assess how well the leave one out cross validation did
ct <- table(grouping = spectral.features$female.id, fit$class)

####Calculate the total percent correct
sum(diag(prop.table(ct))) # 95.7 % accuracy!




   
    ##### Application of a semi-automated vocal fingerprinting approach to monitor Bornean
    ##### gibbon females in an experimentally fragmented landscape in Sabah, Malaysia
    ##### R code written by Dena J. Clink
    
    ##### Contents:
    ##### Part 1 Two methods of MFCC feature extraction    
    ##### Part 2 Recursive Feature Elimination
    ##### Part 3 Comparison of LDA with MFCC and spectrogram feature extraction methods
    ##### Part 4 Random iterations of SVM
    ##### Part 5 Compare LDA and SVM and use SVM on females recorded over three seasons
    

    ##### BEFORE YOU BEGIN  #####
    ##### Install required packages if not already installed
    ##### NOTE: You must remove comment in front of code to get it to run

    # install.packages("tuneR")
    # install.packages("seewave")
    # install.packages("MASS")
    # install.packages("e1071")
    # install.packages("stringr")
    # install.packages("ggplot2")
    # install.packages("viridis")
        
    ##### Load required libraries
    library(tuneR)
    library(seewave)
    library(MASS)
    library(e1071)
    library(stringr)
    library(ggplot2)
    library(viridis)
    
    ##### Download all required files into a single location on your desktop
    ##### Use find and replace function to change /Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master 
    ##### to local file location

    
    ##### Part 1 Two methods of MFCC feature extraction: averaging across all time windows
    ##### and creating a standardized number of time windows for each call
    
    
    ##### Save all .wav files and set working directory to access the files
    ##### There are nine .wav files, three from three different females (SAFA, SAFBA and VJRN_01)
    setwd("/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master")
    
    ####Set the input directory to loop over the wave files
    input.dir <-
      "/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master"
    
    ####List all .wav files in directory
    L = list.files(input.dir, pattern = "*.wav", full.names = FALSE)
    L
    
    ####Extract file names
    filehandles <- str_split_fixed(L, pattern = ".wav", n = 2)[, 1]
    
    ####MFCC Feature Extraction Method 1: Averaging over time windows ###
    
    ####Create empty list to hold MFCC values
    mfcc.vector.list = list()
    
    ####Loop to calculate MFCC for each .wav file in the directory
    for (j in 1:length(filehandles)) {
      filehandle <-  L[j]
      filename <- paste(input.dir, "/", filehandle, sep = "")
      print(paste("processing", filehandle))
      
      # Read in wav file
      w <- readWave(filename)
      
      # Calculate 12 MFCCs for each 0.25 ms window
      melfcc.output <-
        melfcc(
          w,
          minfreq = 400,
          maxfreq = 2000,
          wintime = 0.25,
          fbtype = "mel",
          numcep = 12
        )
      
      #Add MFCC values to list
      mfcc.vector.list[[j]] <- melfcc.output
    }
    
    
    ####Check structure of MFCC list; for each call there is matrix with 12 MFCC for each 0.25 s frame index
    ####Calls of different length will have different number of frames
    str(mfcc.vector.list)
    
    
    ####Create an empty list to store the mean and sd for each mel-frequency bin
    mean.sd.list = list()
    
    ####Loop to calculate MFCC mean and sd
    for (j in 1:length(mfcc.vector.list)) {
      tmp.list <- mfcc.vector.list[[j]]
      list.elements <- lapply(1:ncol(tmp.list), function(x) {
        vec.temp <- tmp.list[, x]
        vec.mean <- mean(vec.temp)
        vec.sd <- sd(vec.temp)
        data.vec <- c(vec.mean, vec.sd)
        data.vec
      })
      mean.sd.list[[j]] <- list.elements
    }
    
    ####Convert the list to a vector
    vec <- unlist(mean.sd.list)
    
    ####Convert the vector into a matrix
    data.matrix.all <-
      matrix(vec,
             nrow = length(filehandles),
             ncol = 24,
             byrow = T)
    data.matrix.all <- as.data.frame(data.matrix.all)
    
    ####Split file names into female id and call id
    loc <- str_split_fixed(filehandles, pattern = "_", n = 2)[, 1]
    group <- str_split_fixed(filehandles, pattern = "_", n = 3)[, 2]
    female.id <- paste(loc, group, sep = "_")
    call.id <-
      str_split_fixed(filehandles, pattern = "[.]", n = 3)[, 1]
    
    ####Create unique column names for the mean and sd of each Mel-frequency bin
    colnames <- rep(c("mean", "sd"), 12)
    nums <- rep(seq(1, 12), 2)
    nums <- sort(nums)
    colnames <- paste(colnames, nums)
    colnames(data.matrix.all) <- colnames
    
    ####Combine into a dataframe for analysis
    mfcc.data.all <-
      cbind.data.frame(data.matrix.all, female.id, call.id, filehandles)
    
    ####Check structure of resulting dataframe
    str(mfcc.data.all)
    
    
    ####MFCC Feature Extraction Method 2: Standardized number of time windows ###
    
    mfcc.vector.list = list()
    
    for (j in 1:length(filehandles)) {
      tryCatch({
        filehandle <-  L[j]
        filename <- paste(input.dir, "/", filehandle, sep = "")
        print(paste("processing", filehandle))
        
        # Read in wav file and filter
        w <- readWave(filename)
        
        # Find duration of .wav file and divide into 9 windows
        wav.dur <- duration(w)
        win.time <- wav.dur / 9
        
        # Calculate MFCCs
        melfcc.output <- melfcc(
          w,
          minfreq = 400,
          hoptime = win.time,
          maxfreq = 2000,
          numcep = 12,
          wintime = win.time
        )
        
        # Calculate delta cepstral coefficients
        deltas.output <- deltas(melfcc.output)
        
        # Ensure only 8 time windows are used for MFCC and delta coefficients
        # Also append .wav duration
        mfcc.vector <-
          c(as.vector(t(melfcc.output[1:8, 2:12])), as.vector(t(deltas.output[1:8, 2:12])), wav.dur)
        
        # Add to list
        mfcc.vector.list[[j]] <- mfcc.vector
      }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      })
    }
    
    ####Convert the list to a vector
    vec <- unlist(mfcc.vector.list)
    
    ####Split file names into female id and call id
    loc <- str_split_fixed(filehandles, pattern = "_", n = 2)[, 1]
    group <- str_split_fixed(filehandles, pattern = "_", n = 3)[, 2]
    female.id <- paste(loc, group, sep = "_")
    call.id <-
      str_split_fixed(filehandles, pattern = "[.]", n = 3)[, 1]
    
    ####Create a matrix with MFCC and delta coefficient values
    ####The number of columns is a reflection of the number of time windows
    data.matrix.all <-
      matrix(
        vec,
        nrow = length(filehandles),
        ncol = length(mfcc.vector.list[[1]]),
        byrow = T
      )
    
    mfcc.data.frame <-
      cbind.data.frame(female.id, filehandles, call.id, data.matrix.all)
    
    seq.length <- (length(mfcc.vector.list[[1]]) - 1) / 2
    
    col.names <-
      c(
        "female.id",
        "filehandles",
        "call.id",
        paste(rep("mfcc", seq.length), seq(seq.length), sep = "_"),
        paste(rep("delta", seq.length), seq(seq.length), sep = "_"),
        "dur"
      )
    
    colnames(mfcc.data.frame) <- col.names
    
    mfcc.data.frame <- as.data.frame(mfcc.data.frame)
    
    ####Check data structure
    str(mfcc.data.frame)
    
    ##### Part 2 Recursive Feature Elimination ####
    
    ##### Load required libraries
    library(MASS)
    library(e1071)
    
    ##### Set source code (downloaded from https://github.com/johncolby/SVM-RFE)
    source("/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master/msvmRFE.R")
    
    ##### Read in MFCC data with standardized number of time window
    mfcc.standard.windows <-
      read.csv(
        "/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master/mfcc.data.with.eight.standardized.windows.csv"
      )
    
    ##### Check structure of data
    str(mfcc.standard.windows)
    
    ##### For SVM RFE the first column needs to include class labels; remove other id columns
    for.svm.rfe <-
      subset(mfcc.standard.windows, select = -c(filehandles, call.id))
    
    ##### We want k=10 for the k-fold cross validation as the “multiple” part of mSVM-RFE.
    svm.rfe.output <- svmRFE(for.svm.rfe, k = 10, halve.above = 100)
    str(svm.rfe.output)
    
    ##### Reorder the data so highest ranked feature is first
    new.svm.rfe <-
      for.svm.rfe[, 2:ncol(for.svm.rfe)][, dput(svm.rfe.output)]
    str(new.svm.rfe)
    
    ##### Create a list to store cross-validation accuracies
    accuracy.list <- list()
    
    ##### Set prior for LDA so that class membership is equally likely
    n.females <- length(unique(mfcc.standard.windows$female.id))
    
    ##### Loop to add features one by one and calculate accuracy using leave-one-out cross-validation
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
      
      ##### Assess how well the leave one out cross validation did
      ct <- table(grouping = for.svm.rfe$female.id, fit.svm.rfe$class)
      
      ##### total percent correct
      percent <- sum(diag(prop.table(ct)))
      print(percent)
      accuracy.list[[j]] <- percent
    }
    
    ##### Create a figure
    plot(unlist(accuracy.list), xlab = "Number of Features", ylab = "Classification Accuracy")
    
    ##### Find which number of features provides the maximum classification accuracy
    max.feature <- which.max(unlist(accuracy.list)) + 1
    
    ##### Subset the highest ranked variables which yield the highest accuracy
    svm.rfe.for.lda <- new.svm.rfe[1:max.feature]
    
    ##### Combine class labels with new subset of features into a data frame for analysis
    svm.rfe.for.classification <-
      cbind.data.frame(mfcc.standard.windows$female.id, svm.rfe.for.lda)
    colnames(svm.rfe.for.classification)[1] <- "female.id"
    
    ##### Run LDA on the data subset using RFE
    fit.standard.number.windows.svm.rfe <- lda(
      svm.rfe.for.classification[2:ncol(svm.rfe.for.classification)],
      center = TRUE,
      prior = rep(1 / n.females, n.females),
      scale. = TRUE,
      grouping = svm.rfe.for.classification$female.id,
      CV = T
    )
    
    ##### Assess how well the leave one out cross validation did
    ct <-
      table(grouping = svm.rfe.for.classification$female.id,
            fit.standard.number.windows.svm.rfe$class)
    
    ##### Calculate total percent correct
    percent <- sum(diag(prop.table(ct)))
    print(percent) # 98.9 % accuracy!
    
    
    ##### Part 3 Comparison of LDA with MFCC and spectrogram feature extraction methods
    
    ##### Load required libraries
    library(MASS)
    
    ####Read in MFCC features chosen via recursive feature elimination
    svm.rfe.for.classification <-
      read.csv("/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master/svm.rfe.for.classification.csv")
    
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
      read.csv("/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master/mfcc.data.mean.and.sd.over.time.windows.csv")
    
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
    
    
    ####Read in features estimated from spectrogram (data from Clink et al. 2017)
    spectral.features <-
      read.csv("/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master/spectral_features_SAFE.data.csv")
    
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
    
    ##### Part 4 Random iterations of SVM ####
    ## NOTE: To save time the code only runs over 10 iterations
    
    ##### Load required libraries
    library(MASS)
    library(e1071)
    library(ggplot2)
    library(viridis)
    
    ##### Load Data
    
    mfcc.data.all <-
      read.csv("/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master/mfcc.data.mean.and.sd.over.time.windows.csv")
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
    mfcc.den.df$Accuracy <- mfcc.den.df$Accuracy * 100
    
    
    ### Plot results
    ggplot(mfcc.den.df, aes(x = Accuracy, fill = Kernel)) + geom_density(alpha =
                                                                           .5, bw = 1) +
      scale_y_continuous(breaks = NULL) +
      scale_fill_viridis(option = "D", discrete = T) +
      guides(fill = guide_legend(title = "Kernel Type")) + xlab("Percent Correct Female Identification") +
      ylab("Count of 100 iterations") +
      theme_classic() +
      theme(
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.title = element_text(size = 20, face = "bold")
      )
    
    ### Calculate mean accuracy for each kernel type    
    mean(mfcc.den.lin$Accuracy)
    mean(mfcc.den.poly$Accuracy)
    mean(mfcc.den.rad$Accuracy)
    mean(mfcc.den.sig$Accuracy)
    

    ### Read in data to recreate Figure 6 Density Plot
    mfcc.den.df <- read.csv("/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master/mfcc.density.dataframe.feature.extraction.method.1.csv")
    
    ### Plot results
    ggplot(mfcc.den.df, aes(x = Accuracy, fill = Kernel)) + geom_density(alpha =
                                                                           .5, bw = 1) +
      scale_y_continuous(breaks = NULL) +
      scale_fill_viridis(option = "D", discrete = T) +
      guides(fill = guide_legend(title = "Kernel Type")) + xlab("Percent Correct Female Identification") +
      ylab("Count of 100 iterations") +
      theme_classic() +
      theme(
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        axis.text.x  = element_text(size = 18),
        axis.title = element_text(size = 20, face = "bold")
      )

    ### Calculate accuracy mean for each kernel type
    mean(subset(mfcc.den.df,Kernel=="Linear")$Accuracy)
    mean(subset(mfcc.den.df,Kernel=="Radial")$Accuracy)
    mean(subset(mfcc.den.df,Kernel=="Polynomial")$Accuracy)
    mean(subset(mfcc.den.df,Kernel=="Sigmoidal")$Accuracy)
    
    ### Calculate accuracy standard deviation for each kernel type    
    sd(subset(mfcc.den.df,Kernel=="Linear")$Accuracy)
    sd(subset(mfcc.den.df,Kernel=="Radial")$Accuracy)
    sd(subset(mfcc.den.df,Kernel=="Polynomial")$Accuracy)
    sd(subset(mfcc.den.df,Kernel=="Sigmoidal")$Accuracy)
    
    
    #### Part 5a. Compare LDA and SVM     #### 
    #### Part 5b. Use SVM on females recorded over three seasons    #### 
    
    
    ##### Load required libraries
    library(MASS)
    library(e1071)
    
    ####Read in MFCC features averaged across all time windows
    mfcc.average.over.windows <-
      read.csv(
        "/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master/mfcc.data.mean.and.sd.over.time.windows.csv"
      )
    
    ####Check structure
    str(mfcc.average.over.windows)
    
    ####Tune the parameters and run the SVM with sigmoid kernel with feature extraction method 1
    tune.sig.method.1 <-
      tune(
        svm,
        mfcc.average.over.windows[, 4:25],
        mfcc.average.over.windows$female.id,
        kernel = "sigmoid",
        ranges = list(
          cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100),
          gamma = c(0.001, 0.01, 0.1, 0.5, 1.0, 2.0)
        )
      )
    
    cost.sig <- tune.sig.method.1$best.parameters$cost
    gamma.sig <-  tune.sig.method.1$best.parameters$gamma
    
    
    svm.sig.method.1 <-
      svm(
        mfcc.average.over.windows[, 4:25],
        mfcc.average.over.windows$female.id,
        kernel = "sigmoid",
        cost = cost.sig,
        gamma = gamma.sig,
        cross = nrow(mfcc.average.over.windows) # When cross = number of observations this indicates leave-one-out cross-validation
      )
    
    svm.sig.method.1$tot.accuracy # 99.5 % accuracy with sigmoid; 99.2 % with radial; 99.4 % accuracy with linear
    
    ####Read in MFCC features chosen via recursive feature elimination
    svm.rfe.for.classification <-
      read.csv("/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master/svm.rfe.for.classification.csv")
    
    ## Check structure
    str(svm.rfe.for.classification)
    
    ####Tune the parameters and run the SVM with sigmoidal kernel with feature extraction method 2
    
    tune.sig.method.2 <-
      tune(
        svm,
        svm.rfe.for.classification[, 2:ncol(svm.rfe.for.classification)],
        svm.rfe.for.classification$female.id,
        kernel = "sigmoid",
        ranges = list(
          cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100),
          gamma = c(0.001, 0.01, 0.1, 0.5, 1.0, 2.0)
        )
      )
    
    cost.sig <- tune.sig.method.2$best.parameters$cost
    gamma.sig <-  tune.sig.method.2$best.parameters$gamma
    
    
    svm.sig.method.2 <-
      svm(
        svm.rfe.for.classification[, 2:ncol(svm.rfe.for.classification)],
        svm.rfe.for.classification$female.id,
        kernel = "sigmoid",
        cost = cost.sig,
        gamma = gamma.sig,
        cross = nrow(svm.rfe.for.classification) # When cross = number of observations this indicates leave-one-out cross-validation
      )
    
    svm.sig.method.2$tot.accuracy # 96.5 % accuracy with sigmoid; 96.81 with radial; 96.01 with linear
    
    ####Read in spectral features estimated from spectrogram (data from Clink et al. 2017)
    spectral.features <-
      read.csv("/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master/spectral_features_SAFE.data.csv")
    
    ####Check structure
    str(spectral.features)
    
    ####Tune the parameters and run the SVM with sigmoidal kernel with spectral features
    tune.sig.spectral <-
      tune(
        svm,
        spectral.features[7:36],
        spectral.features$female.id,
        kernel = "sigmoid",
        ranges = list(
          cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100),
          gamma = c(0.001, 0.01, 0.1, 0.5, 1.0, 2.0)
        )
      )
    
    cost.sig <- tune.sig.spectral$best.parameters$cost
    gamma.sig <-  tune.sig.spectral$best.parameters$gamma
    
    
    svm.sig.spectral <-
      svm(
        spectral.features[7:36],
        spectral.features$female.id,
        kernel = "sigmoid",
        cost = cost.sig,
        gamma = gamma.sig,
        cross = nrow(spectral.features) # When cross = number of observations this indicates leave-one-out cross-validation
      )
    
    svm.sig.spectral$tot.accuracy # 92.3 % accuracy
    
    
    #### Part 6 Use SVM across multiple seasons
    
    ### Read in dataframe with multiple recording seasons
    multiple.seasons.df <-
      read.csv("/Users/denasmacbook/Desktop/MFCC-Vocal-Fingerprinting-master/mfcc.data.three.seasons.csv")
    
    ### Convert female id names to informative labels
    levels(multiple.seasons.df$female.id.test)[levels(multiple.seasons.df$female.id.test) ==
                                                 "SAFL_07"] <-"Female 1"
    
    levels(multiple.seasons.df$female.id.test)[levels(multiple.seasons.df$female.id.test) ==
                                                 "SAFL_01"] <-"Female 1"
    
    levels(multiple.seasons.df$female.id.test)[levels(multiple.seasons.df$female.id.test) ==
                                                 "SAFBN_01"] <-"Female 2"
    
    levels(multiple.seasons.df$female.id.test)[levels(multiple.seasons.df$female.id.test) ==
                                                 "SAFBN_01"] <-"Female 2"
    
    levels(multiple.seasons.df$female.id.test)[levels(multiple.seasons.df$female.id.test) ==
                                                 "SAFBB_01"] <-"Female 3"
    
    levels(multiple.seasons.df$female.id.test)[levels(multiple.seasons.df$female.id.test) ==
                                                 "SAFBA_01"] <-"Female 4"
    
    ### Subset data based on recording season
    season.1.selected <-
      c("SAFBA_01_156_01",
        "SAFBB_01_081_01",
        "SAFBN_01_004",
        "SAFBN_01_005",
        "SAFL_01_034"
      )
    mfcc.data.season.1 <-
      droplevels(multiple.seasons.df[multiple.seasons.df$recording %in% season.1.selected,])
    
    season.2.selected <-
      c("SAFBA_01_037",
        "SAFBN_01_002_01",
        "SAFBB_01_087_01",
        "SAFL_01_002")
    mfcc.data.season.2 <-
      droplevels(multiple.seasons.df[multiple.seasons.df$recording %in% season.2.selected,])
    
    season.3.selected <-
      c("SAFBA_01_013_09",
        "SAFBB_01_011_09",
        "SAFBN_01_047",
        "SAFL_07_021")
    
    mfcc.data.season.3 <-
      droplevels(multiple.seasons.df[multiple.seasons.df$recording %in% season.3.selected,])
    
    
    ## Create dataframe with two seasons
    two.seasons.train.1.2 <-
      rbind.data.frame(mfcc.data.season.1, mfcc.data.season.2)
    
    ## Use tune function to find optimal parameters
    tune.vocal.fingerprint.seasons.1.2 <-
      tune(
        svm,
        two.seasons.train.1.2[, 4:25],
        two.seasons.train.1.2$female.id.test,
        kernel = "sigmoid",
        ranges = list(
          cost = c(0.001, 0.01, 0.1, 1, 2, 10, 100, 1000),
          gamma = c(0.01, 0.1, 0.5, 1.0, 2.0)
        ),
        cross = 5
      )
    
    ## Run SVM model
    svm.model.seasons.1.2  <-
      svm(
        two.seasons.train.1.2[, 4:25],
        two.seasons.train.1.2$female.id.test,
        kernel = "sigmoid",
        gamma = tune.vocal.fingerprint.seasons.1.2$best.parameters$gamma,
        cost = tune.vocal.fingerprint.seasons.1.2$best.parameters$cost,
        cross = nrow(two.seasons.train.1.2)
      )
    
    ## Predict class membership of season not used in training
    mod.pred.seasons.1.2  <- predict(svm.model.seasons.1.2, mfcc.data.season.3[, 4:25])
    
    ## Create a table of predicted vs. actual and calculate classification accuracy
    table(mod.pred.seasons.1.2, mfcc.data.season.3$female.id.test)
    sum(diag(prop.table(
      table(mod.pred.seasons.1.2, mfcc.data.season.3$female.id.test) ## 94.8 % accuracy
    )))
    
    ## Create dataframe with two seasons
    two.seasons.train.seasons.1.3 <-
      rbind.data.frame(mfcc.data.season.1, mfcc.data.season.3)
    
    ## Use tune function to find optimal parameters
    tune.vocal.fingerprint.seasons.1.3 <-
      tune(
        svm,
        two.seasons.train.seasons.1.3[, 4:25],
        two.seasons.train.seasons.1.3$female.id.test,
        kernel = "sigmoid",
        ranges = list(
          cost = c(0.001, 0.01, 0.1, 1, 2, 10, 100, 1000),
          gamma = c(0.01, 0.1, 0.5, 1.0, 2.0)
        ),
        cross = 5
      )
    
    ## Run SVM model
    svm.model.seasons.1.3 <-
      svm(
        two.seasons.train.seasons.1.3[, 4:25],
        two.seasons.train.seasons.1.3$female.id.test,
        kernel = "sigmoid",
        gamma = tune.vocal.fingerprint.seasons.1.3$best.parameters$gamma,
        cost = tune.vocal.fingerprint.seasons.1.3$best.parameters$cost,
        cross = 5
      )
    
    ## Predict class membership of season not used in training
    mod.pred.seasons.1.3 <- predict(svm.model.seasons.1.3, mfcc.data.season.2[, 4:25])
    
    ## Create a table of predicted vs. actual and calculate classification accuracy
    table(mod.pred.seasons.1.3, mfcc.data.season.2$female.id.test)
    sum(diag(prop.table(
      table(mod.pred.seasons.1.3, mfcc.data.season.2$female.id.test) ## 61.7% accuracy
    )))
    
    
    ## Create dataframe with two seasons
    two.seasons.train.seasons.2.3 <-
      rbind.data.frame(mfcc.data.season.2, mfcc.data.season.3)
    
    ## Use tune function to find optimal parameters
    tune.vocal.fingerprint.seasons.2.3 <-
      tune(
        svm,
        two.seasons.train.seasons.2.3[, 4:25],
        two.seasons.train.seasons.2.3$female.id.test,
        kernel = "sigmoid",
        ranges = list(
          cost = c(0.001, 0.01, 0.1, 1, 2, 10, 100, 1000),
          gamma = c(0.01, 0.1, 0.5, 1.0, 2.0)
        ),
        cross = 5
      )
    
    ## Run SVM model
    svm.model.seasons.2.3 <-
      svm(
        two.seasons.train.seasons.2.3[, 4:25],
        two.seasons.train.seasons.2.3$female.id.test,
        kernel = "sigmoid",
        gamma = tune.vocal.fingerprint.seasons.2.3$best.parameters$gamma,
        cost = tune.vocal.fingerprint.seasons.2.3$best.parameters$cost,
        cross = 5
      )
    
    ## Predict class membership of season not used in training
    mod.pred.seasons.2.3 <- predict(svm.model.seasons.2.3, mfcc.data.season.1[, 4:25])
    
    ## Create a table of predicted vs. actual and calculate classification accuracy
    table(mod.pred.seasons.2.3, mfcc.data.season.1$female.id.test)
    sum(diag(prop.table(
      table(mod.pred.seasons.2.3, mfcc.data.season.1$female.id.test) ## 78.4% accuracy
    )))
    
    

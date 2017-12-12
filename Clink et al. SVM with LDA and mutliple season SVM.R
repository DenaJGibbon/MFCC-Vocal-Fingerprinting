    ##### R code for Clink et al.
    
    ##### Application of a semi-automated vocal fingerprinting approach to monitor Bornean
    ##### gibbon females in an experimentally fragmented landscape in Sabah, Malaysia
    
    
    ##### Load required libraries
    library(MASS)
    library(e1071)

    
    
    #### Part 5a. Compare LDA and SVM
    #### Part 5b. Use SVM on females recorded over three seasons
    
    ####Read in MFCC features averaged across all time windows
    mfcc.average.over.windows <-
      read.csv("/Users/denasmacbook/Desktop/Clink_et_at_MFCC/MFCC.data.csv")
    
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
      read.csv("/Users/denasmacbook/Desktop/Clink_et_at_MFCC/svm.rfe.for.classification.csv")
    
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
      read.csv("/Users/denasmacbook/Desktop/Clink_et_at_MFCC/spectral_features_SAFE.data.csv")
    
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
      read.csv("/Users/denasmacbook/Downloads/mfcc.data.three.seasons.csv")
    
    ### Convert female id names to informative labels
    levels(multiple.seasons.df$female.id.test)[levels(multiple.seasons.df$female.id.test) ==
                                                 "SAFL_07"] <- "Female 1"
    levels(multiple.seasons.df$female.id.test)[levels(multiple.seasons.df$female.id.test) ==
                                                 "SAFL_01"] <- "Female 1"
    levels(multiple.seasons.df$female.id.test)[levels(multiple.seasons.df$female.id.test) ==
                                                 "SAFBN_01"] <- "Female 2"
    levels(multiple.seasons.df$female.id.test)[levels(multiple.seasons.df$female.id.test) ==
                                                 "SAFBN_01"] <- "Female 2"
    levels(multiple.seasons.df$female.id.test)[levels(multiple.seasons.df$female.id.test) ==
                                                 "SAFBB_01"] <- "Female 3"
    levels(multiple.seasons.df$female.id.test)[levels(multiple.seasons.df$female.id.test) ==
                                                 "SAFBA_01"] <- "Female 4"
    
    ### Subset data based on recording season
    season.1.selected <-
      c("SAFBA_01_156_01", "SAFBB_01_081_01","SAFBN_01_004", "SAFBN_01_005","SAFL_01_034")
    mfcc.data.season.1 <-
      droplevels(multiple.seasons.df[multiple.seasons.df$recording %in% season.1.selected, ])
    
    season.2.selected <-
      c("SAFBA_01_037","SAFBN_01_002_01","SAFBB_01_087_01","SAFL_01_002")
    mfcc.data.season.2 <-
      droplevels(multiple.seasons.df[multiple.seasons.df$recording %in% season.2.selected, ])
    
    season.3.selected <-
      c("SAFBA_01_013_09", "SAFBB_01_011_09","SAFBN_01_047", "SAFL_07_021")
       
    mfcc.data.season.3 <-
      droplevels(multiple.seasons.df[multiple.seasons.df$recording %in% season.3.selected, ])
    
    
    ## Create dataframe with two seasons
    two.seasons.train <-
      rbind.data.frame(mfcc.data.season.1, mfcc.data.season.2)
    
    ## Use tune function to find optimal parameters
    tune.vocal.fingerprint <-
      tune(
        svm,
        two.seasons.train[, 4:25],
        two.seasons.train$female.id.test,
        kernel = "sigmoid",
        ranges = list(
          cost = c(0.001, 0.01, 0.1, 1, 2, 10, 100, 1000),
          gamma = c(0.01, 0.1, 0.5, 1.0, 2.0)
        ),
        cross = 5
      )
    
    ## Run SVM model
    svm.model <-
      svm(
        two.seasons.train[, 4:25],
        two.seasons.train$female.id.test,
        kernel = "sigmoid",
        gamma = tune.vocal.fingerprint$best.parameters$gamma,
        cost = tune.vocal.fingerprint$best.parameters$cost,
        cross = nrow(two.seasons.train)
      )
    
    ## Predict class membership of season not used in training
    mod.pred <- predict(svm.model, mfcc.data.season.3[, 4:25])
    
    ## Create a table of predicted vs. actual and calculate classification accuracy
    table(mod.pred, mfcc.data.season.3$female.id.test)
    sum(diag(prop.table(
      table(mod.pred, mfcc.data.season.3$female.id.test) ##
    )))
    
    ## Create dataframe with two seasons
    two.seasons.train <-
      rbind.data.frame(mfcc.data.season.1, mfcc.data.season.3)
    
    ## Use tune function to find optimal parameters
    tune.vocal.fingerprint <-
      tune(
        svm,
        two.seasons.train[, 4:25],
        two.seasons.train$female.id.test,
        kernel = "sigmoid",
        ranges = list(
          cost = c(0.001, 0.01, 0.1, 1, 2, 10, 100, 1000),
          gamma = c(0.01, 0.1, 0.5, 1.0, 2.0)
        ),
        cross = 5
      )
    
    ## Run SVM model
    svm.model <-
      svm(
        two.seasons.train[, 4:25],
        two.seasons.train$female.id.test,
        kernel = "sigmoid",
        gamma = tune.vocal.fingerprint$best.parameters$gamma,
        cost = tune.vocal.fingerprint$best.parameters$cost,
        cross = 5
      )
    
    ## Predict class membership of season not used in training
    mod.pred <- predict(svm.model, mfcc.data.season.2[, 4:25])

    ## Create a table of predicted vs. actual and calculate classification accuracy
    table(mod.pred, mfcc.data.season.2$female.id.test)
    sum(diag(prop.table(
      table(mod.pred, mfcc.data.season.2$female.id.test) ## 61.7% accuracy 
    )))
    
    
    ## Create dataframe with two seasons
    two.seasons.train <-
      rbind.data.frame(mfcc.data.season.2, mfcc.data.season.3)
    
    ## Use tune function to find optimal parameters
    tune.vocal.fingerprint <-
      tune(
        svm,
        two.seasons.train[, 4:25],
        two.seasons.train$female.id.test,
        kernel = "sigmoid",
        ranges = list(
          cost = c(0.001, 0.01, 0.1, 1, 2, 10, 100, 1000),
          gamma = c(0.01, 0.1, 0.5, 1.0, 2.0)
        ),
        cross = 5
      )
    
    ## Run SVM model
    svm.model <-
      svm(
        two.seasons.train[, 4:25],
        two.seasons.train$female.id.test,
        kernel = "sigmoid",
        gamma = tune.vocal.fingerprint$best.parameters$gamma,
        cost = tune.vocal.fingerprint$best.parameters$cost,
        cross = 5
      )
    
    ## Predict class membership of season not used in training
    mod.pred <- predict(svm.model, mfcc.data.season.1[, 4:25])
    
    ## Create a table of predicted vs. actual and calculate classification accuracy
    table(mod.pred, mfcc.data.season.1$female.id.test)
    sum(diag(prop.table(
      table(mod.pred, mfcc.data.season.1$female.id.test) ## 78.4% accuracy
    )))
    
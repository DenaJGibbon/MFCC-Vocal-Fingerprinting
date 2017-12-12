##### R code for Clink et al.

##### Application of a semi-automated vocal fingerprinting approach to monitor Bornean
##### gibbon females in an experimentally fragmented landscape in Sabah, Malaysia

##### Part 1 Two methods of MFCC feature extraction: averaging across all time windows
##### and creating a standardized number of time windows for each call


##### Load required libraries
library(tuneR)
library(seewave)
library(MASS)
library(e1071)
library(stringr)
library(ggplot2)
library(viridis)



##### Part 1 MFCC Estimation ####

##### Save all .wav files and set working directory to access the files
##### There are nine .wav files, three from three different females (SAFA, SAFBA and VJRN_01)
setwd("/Users/denasmacbook/Downloads/MFCC-Vocal-Fingerprinting-master")

####Set the input directory to loop over the wave files
input.dir <-
  "/Users/denasmacbook/Downloads/MFCC-Vocal-Fingerprinting-master"

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

# packages <- c('tuneR', 'seewave', 'fftw','neuralnet', 'caTools', 'randomForest', 'warbleR', 'mice', 'e1071', 'rpart', 'rpart.plot', 'xgboost', 'e1071')
# if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
#   install.packages(setdiff(packages, rownames(installed.packages())))  
# }

library(tuneR)
library(seewave)
library(caTools)
library(rpart)
library(rpart.plot)
library(randomForest)
library(warbleR)
library(mice)
library(e1071)
library(neuralnet)


specan3 <- function(X, bp = c(0,22), wl = 2048, threshold = 5, parallel = 1){
  # To use parallel processing: library(devtools), install_github('nathanvan/parallelsugar')
  if(class(X) == "data.frame") {if(all(c("sound.files", "selec", 
                                         "start", "end") %in% colnames(X))) 
  {
    start <- as.numeric(unlist(X$start))
    end <- as.numeric(unlist(X$end))
    sound.files <- as.character(unlist(X$sound.files))
    selec <- as.character(unlist(X$selec))
  } else stop(paste(paste(c("sound.files", "selec", "start", "end")[!(c("sound.files", "selec", 
                                                                        "start", "end") %in% colnames(X))], collapse=", "), "column(s) not found in data frame"))
  } else  stop("X is not a data frame")
  
  #if there are NAs in start or end stop
  if(any(is.na(c(end, start)))) stop("NAs found in start and/or end")  
  
  #if end or start are not numeric stop
  if(all(class(end) != "numeric" & class(start) != "numeric")) stop("'end' and 'selec' must be numeric")
  
  #if any start higher than end stop
  if(any(end - start<0)) stop(paste("The start is higher than the end in", length(which(end - start<0)), "case(s)"))  
  
  #if any selections longer than 20 secs stop
  if(any(end - start>20)) stop(paste(length(which(end - start>20)), "selection(s) longer than 20 sec"))  
  options( show.error.messages = TRUE)
  
  #if bp is not vector or length!=2 stop
  if(!is.vector(bp)) stop("'bp' must be a numeric vector of length 2") else{
    if(!length(bp) == 2) stop("'bp' must be a numeric vector of length 2")}
  
  #return warning if not all sound files were found
  fs <- list.files(path = getwd(), pattern = ".wav$", ignore.case = TRUE)
  if(length(unique(sound.files[(sound.files %in% fs)])) != length(unique(sound.files))) 
    cat(paste(length(unique(sound.files))-length(unique(sound.files[(sound.files %in% fs)])), 
              ".wav file(s) not found"))
  
  #count number of sound files in working directory and if 0 stop
  d <- which(sound.files %in% fs) 
  if(length(d) == 0){
    stop("The .wav files are not in the working directory")
  }  else {
    start <- start[d]
    end <- end[d]
    selec <- selec[d]
    sound.files <- sound.files[d]
  }
  
  # If parallel is not numeric
  if(!is.numeric(parallel)) stop("'parallel' must be a numeric vector of length 1") 
  if(any(!(parallel %% 1 == 0),parallel < 1)) stop("'parallel' should be a positive integer")
  
  # If parallel was called
  if(parallel > 1)
  { options(warn = -1)
    if(all(Sys.info()[1] == "Windows",requireNamespace("parallelsugar", quietly = TRUE) == TRUE)) 
      lapp <- function(X, FUN) parallelsugar::mclapply(X, FUN, mc.cores = parallel) else
        if(Sys.info()[1] == "Windows"){ 
          cat("Windows users need to install the 'parallelsugar' package for parallel computing (you are not doing it now!)")
          lapp <- pbapply::pblapply} else lapp <- function(X, FUN) parallel::mclapply(X, FUN, mc.cores = parallel)} else lapp <- pbapply::pblapply
  
  options(warn = 0)
  
  if(parallel == 1) cat("Measuring acoustic parameters:")
  x <- as.data.frame(lapp(1:length(start), function(i) { 
    r <- tuneR::readWave(file.path(getwd(), sound.files[i]), from = start[i], to = end[i], units = "seconds") 
    
    b<- bp #in case bp its higher than can be due to sampling rate
    if(b[2] > ceiling(r@samp.rate/2000) - 1) b[2] <- ceiling(r@samp.rate/2000) - 1 
    
    
    #frequency spectrum analysis
    songspec <- seewave::spec(r, f = r@samp.rate, plot = FALSE)
    analysis <- seewave::specprop(songspec, f = r@samp.rate, flim = c(0, 280/1000), plot = FALSE)
    
    #save parameters
    meanfreq <- analysis$mean/1000
    sd <- analysis$sd/1000
    median <- analysis$median/1000
    Q25 <- analysis$Q25/1000
    Q75 <- analysis$Q75/1000
    IQR <- analysis$IQR/1000
    skew <- analysis$skewness
    kurt <- analysis$kurtosis
    sp.ent <- analysis$sh
    sfm <- analysis$sfm
    mode <- analysis$mode/1000
    centroid <- analysis$cent/1000
    
    #Frequency with amplitude peaks
    peakf <- 0#seewave::fpeaks(songspec, f = r@samp.rate, wl = wl, nmax = 3, plot = FALSE)[1, 1]
    
    #Fundamental frequency parameters
    ff <- seewave::fund(r, f = r@samp.rate, ovlp = 50, threshold = threshold, 
                        fmax = 280, ylim=c(0, 280/1000), plot = FALSE, wl = wl)[, 2]
    meanfun<-mean(ff, na.rm = T)
    minfun<-min(ff, na.rm = T)
    maxfun<-max(ff, na.rm = T)
    
    #Dominant frecuency parameters
    y <- seewave::dfreq(r, f = r@samp.rate, wl = wl, ylim=c(0, 280/1000), ovlp = 0, plot = F, threshold = threshold, bandpass = b * 1000, fftw = TRUE)[, 2]
    meandom <- mean(y, na.rm = TRUE)
    mindom <- min(y, na.rm = TRUE)
    maxdom <- max(y, na.rm = TRUE)
    dfrange <- (maxdom - mindom)
    duration <- (end[i] - start[i])
    
    #modulation index calculation
    changes <- vector()
    for(j in which(!is.na(y))){
      change <- abs(y[j] - y[j + 1])
      changes <- append(changes, change)
    }
    if(mindom==maxdom) modindx<-0 else modindx <- mean(changes, na.rm = T)/dfrange
    
    #save results
    return(c(duration, meanfreq, sd, median, Q25, Q75, IQR, skew, kurt, sp.ent, sfm, mode, 
             centroid, peakf, meanfun, minfun, maxfun, meandom, mindom, maxdom, dfrange, modindx))
  }))
  
  #change result names
  
  rownames(x) <- c("duration", "meanfreq", "sd", "median", "Q25", "Q75", "IQR", "skew", "kurt", "sp.ent", 
                   "sfm","mode", "centroid", "peakf", "meanfun", "minfun", "maxfun", "meandom", "mindom", "maxdom", "dfrange", "modindx")
  x <- data.frame(sound.files, selec, as.data.frame(t(x)))
  colnames(x)[1:2] <- c("sound.files", "selec")
  rownames(x) <- c(1:nrow(x))
  
  return(x)
}

voicedata <- read.csv(file="voice.csv", header=TRUE, sep=",")
data <- as.data.frame(voicedata)
colnames(data) <- c("meanfreq", "sd", "median", "Q25", "Q75", "IQR", "skew", "kurt", "sp.ent", "sfm","mode", "centroid","meanfun", "minfun", "maxfun", "meandom", "mindom", "maxdom", "dfrange", "modindx","label")
dataset <- data

# # Create a train and test set.
 set.seed(777)
 spl <- sample.split(data$label, 0.8)
 train <- subset(data, spl == TRUE)
 test <- subset(data, spl == FALSE)
 print("splitting done")

#normalise the data nn ke liye
 # maxs <- apply(train[,1:20], 2, max)       #2 for columns
 # mins <- apply(train[,1:20], 2, min)
 # scaled.data <- as.data.frame(scale(train[,1:20],center = TRUE, scale = maxs-mins))
 Label <- as.numeric(train$label)-1     #put male 1 and female -1
 Label[Label==0]<- -1
 datann <- cbind(Label,train[,1:20])
 print("normalization for nn done")

# generate formula for NN
 feats <- names(scaled.data)
 f<- paste(feats,collapse='+')
 f<- paste('Label ~',f)
 f<- as.formula(f)
 print("formulas generated")

# Build models.
  print("Training neural network")
  genderNn<- neuralnet(f,data = datann, hidden=4,linear.output=FALSE)
  print("Neural net model trained")
 #  genderLog <- glm(label ~ ., data=train, family='binomial')
 #  print("Log model trained")
 #  genderCART <- rpart(label ~ ., data=train, method='class')
 #  print("CART model trained")
 #  genderForest <- randomForest(label ~ ., data=train)
 #  genderTunedForest <- tuneRF(train[, -21], train[, 21], stepFactor=.5, doBest=TRUE)
 #  print("Random forest trained")
 #  genderSvm <- svm(label ~ ., data=train, gamma =0.21, cost=8)
 #  print("SVM trained")

 # # #stacked ensamble model
 #  results1 <- predict(genderSvm, newdata=test)
 #  results2 <- predict(genderTunedForest, newdata=test)
  
  # result3a  <- compute(genderNn,test[,1:20])
  # results3  <- sapply(result3a$net.result,round,digits=0)
  # mf <- c("male","female")
  # results3[result3a$net.result >= 0.5] <- mf[2]
  # results3[result3a$net.result <= 0.5] <- mf[1]
 
  # combo <- data.frame(results1, results2, results3,y = test$label)
#   print(combo)
#   genderStacked <- tuneRF(combo[,-4], combo[,4], stepFactor=.5, doBest=TRUE)
#   print("stacked model trained")

#  #model plots
#   print("CART summary:")
#   summary(genderCART)
#   plot(genderNn)
#   print("RANDOM FOREST SUMMARY :")
#   print(genderTunedForest)
#   print("SVM SUMMARY")
#   print(genderSvm)



# #save models
#   setwd("/home/cloudera/Desktop/MINOR")
#   save(genderCART,file = "CART.bin")
#   save(genderTunedForest,file = "tunedForest.bin")
#   save(genderNn,file = "NeuralNet.bin")
#   save(genderSvm,file = "svm.bin")
#   save(genderLog,file = "Log.bin")
#   save(genderStacked,file="stacked.bin")



# # Accuracy: 0.72
# predictLog <- predict(genderLog, type='response')
# table(train$label, predictLog >= 0.5)
 

# # Accuracy: 0.78
# predictCART2 <- predict(genderCART, newdata=test)
# predictCART2.prob <- predictCART2[,2]
# table(test$label, predictCART2.prob >= 0.5)
 


# # Accuracy: 0.86
# predictForest <- predict(genderForest, newdata=test)
# table(test$label, predictForest)
 


# # Accuracy: 0.85
# predictSvm <- predict(genderSvm, test)
# table(predictSvm, test$label)


# Accuracy : 0.82
  predicted.nn.values <- compute(genderNn,test[,1:20])
  res <- predicted.nn.values
  res[predicted.nn.values$net.result >= 1] <- 'male'
  res[predicted.nn.values$net.result < 1 ]  <- 'female'
  res<- unlist(res)
# nnresult <- sapply(predicted.nn.values$net.result,round,digits=0)
 print( table(test$label,res))# false false

##  Accuracy : 0.95
# results1 <- predict(genderSvm, newdata=test)
# results2 <- predict(genderTunedForest, newdata=test)
# combo <- data.frame(results1, results2, results3)
# result <- predict(genderStacked, newdata=combo)
# table(test$label,result)


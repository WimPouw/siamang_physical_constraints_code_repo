library(rstudioapi)           # package for relative foldering
library(EMD)                  #Empirical Mode decomposition
library(rPraat)               #reading textgrds
library(signal)               # Signal processing
library(dplR)                 # signal processing
library(stringr)              #string manipulation
library(zoo)                  #approximation
library(kza)                  #smoothing
library(wrassp)               #acoustic processing
library(soundgen)             #acoustic analsyes in R
library(pracma)               #time series functions and signal processing
library(tuneR)                #acoustics
library(phonTools)            #acoustics
library(FactoMineR)           #PCA and other data science functions

##############################INFO
#CONTACT: wim.pouw@donders.ru.nl
#Version: 3
#Manuscript: Cross-modal constraints in multimodal vocalizations in Siamang (Syndactylus symphalangus) 
#Script info: This script concerns the preprocessing and aggregation of acoustic and motion tracking data, which will result in a dataset
# that is used for statistical analyses described in a different script.


#########Folder and files
parentfolder <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))  #what is the current folder?

audiofold  <- paste0(parentfolder, '/original_data/audio/')
mtfold  <- paste0(parentfolder, '/original_data/motion_tracking/')
audiofiles <- list.files(audiofold) #this contains all the audio wherefore we handchecked whether the motion tracking is ok (others are excluded)
motfiles <- list.files(mtfold, pattern= '*filtered.csv')
annotations <- paste0(parentfolder, '/R_scripts/processing_and_analyses/annotations/')
proc_fol <- paste0(parentfolder, '/R_scripts/processing_and_analyses/processed_data/')
timeseriesfol <- paste0(proc_fol, 'time_series/')
#latest trained model with deeplabcut
dlcmodel <- 'DLC_resnet50_SiamangSimpleV2Jan31shuffle1_500000_filtered'

#just check what days we need to process (instead having for each day 4 cameras or less)
mf <- str_replace(motfiles, paste0(dlcmodel, ".csv"), "")
mf <- str_replace(mf, 'snip_gopro_orange_opp_', "")
mf <- str_replace(mf, 'snip_gopro_blue_opp_', "")
mf <- str_replace(mf, 'snip_gopro_cyan_opp_', "")
mf <- str_replace(mf, 'snip_gopro_green_opp_', "")
#unique(mf) will now give the unique days that we need to process 

#####################################################FUNCTIONS OPEN
pca.speed <- function(mat) #take a matrix with x,y position data of multiple cameras
{
  size <- ncol(mat) #determine the number of position traces
  speed <- kz(sqrt(diff(mat[,1])^2+diff(mat[,2])^2), 5, 5) #if we only have 1 camera, just compute the speed from 2D info (we also smooth again)

  if(size>2) #if we have more than 1 camera (i.e., > 2 position traces), we can try to disentangle independent components, for x, y, z, for which we calculate a speed magnitude
  {
  PCA = prcomp(mat)
  f1 <- diff(PCA$x[,'PC1'])^2 #differentiate over PC1
  f2 <- diff(PCA$x[,'PC2'])^2 #differentiate over PC2
  f3 <- diff(PCA$x[,'PC3'])^2 #differentiate over PC3
  speed <- kz(c(0, sqrt(f1+f2+f3)),5,5) #take the L2 norm of delta PCx and smooth with kolmogorov zurbenko filter
  }
  return(speed)
}

#FUNCTION loads an event log (begintime,endtime,event_id) into a time series data
#time_ms_rec is the time vector form the time series
#g_d is the event log
##this function loads into a time series an annotation dataframe with columns (begintime, endtime, annotation)
load.in.event <- function(time_original, anno)
{
  output <- rep(NA, length(time_original))
  for(i in 1:length(anno[,1])) #loop through each annotation event
  {
    #print(anno[i,3]) #you can print the value to check if this is running correctly
    output <- ifelse((time_original >= anno[i,1] & time_original <= anno[i,2]), as.character(anno[i,3]), output)
  }
  return(output)
}

#THis function allows you to make event identifiers, so that each run of observations in a time series
# gets an id.
make.event.identifiers <- function(runsval, originaltime)
{ #MAKE SPEECH EVENTS
  rle_x <- rle(runsval)
  values = rle_x$values
  endtime <- cumsum(rle_x$lengths)
  begintime <- c(1, diff(endtime))
  begintime <- (endtime - begintime)+1
  runs <- data.frame(originaltime[begintime], originaltime[endtime], values)
  runs <- subset(runs, values == 1)
  runs <- subset(runs, select =-c(values))
  runs$s_events <- seq.int(nrow(runs))
  events <<- load.in.event(originaltime, runs)
}

#initialize variables for main dataset
videonames<- duration<- angles<- individual<-locomote_posture<- acc_nearest_env_peak<- acc_nearest_f0_peak<-
maxf0<- synchrony<- maxacc <- maxamp <- vector()

#####################################################MAIN PROCESSING ROUTINE
for(snippet in unique(mf))
{
  #progress information
  print(paste0('working on session:', snippet))
  
  #get the annotation and do some text wrangling so we can match things up
  sub <- str_split(paste0(snippet,'_session'), "\\_session", simplify=T)
  annos <- read.csv(paste0(annotations,'opp_',sub[,1], '_session', substr(sub[,2],1,1), '_siamang_locomotionvocal_KePoV1.csv'))
  annos_loco <- read.csv(paste0(annotations,'opp_',sub[,1], '_session', substr(sub[,2],1,1), '_siamang_locomotionvocal_KePoV1_locomotion_sequence.csv'))
  lastcode <- str_split(sub[,2], "\\_", simplify=T)[,3]
  sequencenumber <- as.numeric(str_split(lastcode , "\\.wav", simplify=T)[,1])
  anno <- annos[annos[,4]==sequencenumber,]
  loco <- annos_loco$Sequence_movement[annos_loco[,4]==sequencenumber]
  if(is.null(loco)){loco <- 'undefined'} #make sure that the locomotion types with out entries are coded as undefnined
  if(loco == ""){loco <- 'undefined'}  #make sure that the locomotion types with out entries are coded as undefnined
  
  ############################################ACOUSTICS
  #get the waveform
  snd <- rPraat::snd.read(paste0(audiofold, 'opp_', snippet, '.wav'))
    #produce a smoothed amplitude envelope
  hilb <- seewave::hilbert(snd$sig, f = snd$fs, fftw =FALSE) #apply hilbert
  envelope  <- as.vector(abs(hilb)) #take the complex modulus
  envelope <- dplR::hanning(x= envelope, n = snd$fs/12) #12 hertz hanning window
  envelope[is.na(envelope)] <- 0 #make sure that undefined edges get a value of 0, instead of NA
  
  ############################################KINEMATICS
  #get the cameras we need to process
  cameras_to_process  <- motfiles[which(mf == snippet)]
  
  percentagesg <- vector()
  for(camera in cameras_to_process) #lets check which cameras have reliable tracking
  {
    #load in DLC motion data
    mt <- read.csv(paste0(mtfold, camera))
    mt <- mt[-c(1, 2),]
    colnames(mt) <- c('frame', 'x_thorax', 'y_thorax', 'likelihood_thorax', 'x_bum', 'y_bum', 'likelihood_bum')
    mt <- as.data.frame(apply(mt,2, function(x) as.numeric(as.character(x))))
    #check how many likilhood under .70
    percentagesg <- c(percentagesg, sum(mt$likelihood_thorax>.80)/nrow(mt)) #what is the percentege of data above .80?
  }
  
  cameras_to_process <- cameras_to_process[percentagesg > .95] #only keep those camera percentage with 95% highly reliable (> .80) data
  number_of_angles <- length(cameras_to_process) #what is the number of angles we have left?
  speed <- speed1 <- speed2 <- speed3 <- speed4 <- 0 #intitialize some variables
  for(camera in cameras_to_process)
  {
    mt <- read.csv(paste0(mtfold, camera))
    mt <- mt[-c(1, 2),]
    colnames(mt) <- c('frame', 'x_thorax', 'y_thorax', 'likelihood_thorax', 'x_bum', 'y_bum', 'likelihood_bum')
    mt <- as.data.frame(apply(mt,2, function(x) as.numeric(as.character(x))))
    
    #any likelihood lower than .80 will be interpolated
    mt$x_thorax[mt$likelihood_thorax < .80] <- NA
    mt$y_thorax[mt$likelihood_thorax < .80] <- NA
    mt$x_thorax <- na.approx(mt$x_thorax, na.rm=FALSE) #interpolate
    mt$y_thorax <- na.approx(mt$y_thorax, na.rm=FALSE) #interpolate 
    mt$x_thorax <- na.locf(mt$x_thorax, fromLast = TRUE) #trailing Na's get nearest value
    mt$y_thorax <- na.locf(mt$y_thorax, fromLast = TRUE) #trailing Na's get nearest value
    mt$x_thorax_sm <- kz(mt$x_thorax, 9, 11) #normalize by body size and smooth
    mt$y_thorax_sm <- kz(mt$y_thorax, 9, 11) #normalize by body size and smooth
    mt$x_bum_sm <- kz(mt$x_bum, 9, 11) #normalize by body size and smooth
    mt$y_bum_sm <- kz(mt$y_bum, 9, 11) #normalize by body size and smooth
    
    #get mean body proportion for when likelihoods are 1.0 so that we can normalize
    mts <- mt[mt$likelihood_thorax == 1,]
    bodysize_ts <- sqrt((mts$x_thorax_sm-mts$x_bum_sm)^2+(mts$y_thorax_sm-mts$y_bum_sm)^2) #distance bum-thorax
    bodysize <- mean(bodysize_ts) #this will give the mean body size that we will use to normalize speed
    
    if(nrow(mts)==0){  number_of_angles = number_of_angles-1} # we actually cant use this angle if confidences are so low
    if(nrow(mts)!=0) #if there is no data with data points at least confident normalization
    {
    #smooth traces
    mt$x_thorax_sm <- mt$x_thorax_sm/bodysize #normalize by body size and smooth
    mt$y_thorax_sm <- mt$y_thorax_sm/bodysize #normalize by body size and smooth
      
    #save the info for each camera
    if(which(camera == cameras_to_process)==1){speed1 = cbind(mt$x_thorax_sm, mt$y_thorax_sm)}
    if(which(camera == cameras_to_process)==2){speed2 = cbind(mt$x_thorax_sm, mt$y_thorax_sm)}
    if(which(camera == cameras_to_process)==3){speed3 = cbind(mt$x_thorax_sm, mt$y_thorax_sm)}
    if(which(camera == cameras_to_process)==4){speed4 = cbind(mt$x_thorax_sm, mt$y_thorax_sm)}
    
    #PCA apply to approximate 3D info (uses the PCA function above)
    if(length(speed1) != 1 & length(speed2) == 1 & length(speed3) == 1 & length(speed4) == 1){speed = pca.speed(cbind(speed1))}
    if(length(speed1) != 1 & length(speed2) != 1 & length(speed3) == 1 & length(speed4) == 1) #if two cameras are available
    {
      len <- min(c(nrow(speed1), nrow(speed2))) #there may be a single frame distance due to different frame boundaries of the cameras, so we make sure we take the same amount of frames for each camera
      speed = pca.speed(cbind(speed1[1:len,1:2], speed2[1:len,1:2]))
      }
    if(length(speed1) != 1 & length(speed2) != 1 & length(speed3) != 1 & length(speed4) == 1) #if three cameras are available
    {
      len <- min(c(nrow(speed1), nrow(speed2), nrow(speed3)))  
      speed = pca.speed(cbind(speed1[1:len,1:2], speed2[1:len,1:2], speed3[1:len,1:2]))
    }
    if(length(speed1) != 1 & length(speed2) != 1 & length(speed3) != 1 & length(speed4) != 1) #if four cameras are available
    {
      len <- min(c(nrow(speed1), nrow(speed2), nrow(speed3), nrow(speed4)))  
      speed = pca.speed(cbind(speed1[1:len,1:2], speed2[1:len,1:2], speed3[1:len,1:2], speed4[1:len,1:2]))
    }
      }
    }
  if(length(speed)>1)
  {
  acc <- c(0, kz(diff(speed), 5, 5)) #to get the acceleration we differentiate over the new 3D PCA speed, and smooth once more
  
  mt <- cbind.data.frame(speed, acc)
  mt$time_ms <- seq(0, (nrow(mt)-1)*(1000/50), by = 1000/50)
    
    #resample settings at desired sampling rate (100Hz)
    f <- approxfun(1:(snd$duration*snd$fs),envelope)
    #resample apply
    
    envelope <- f(seq(from=0,to=snd$duration*snd$fs,by=snd$fs/100))
    time_ms <- seq(0, (length(envelope)-1)*1000/100, by=1000/100)
    envelope <- cbind.data.frame(time_ms, envelope)
    ################################################################### F0 and postprocess
      #Note this variable is not of main interest for the study, but we do calculate it for exploratory purposes
    f0 <- ksvF0(paste0(audiofold, 'opp_',snippet, '.wav'),maxF = 800,  minF = 300, windowShift=1000/100, toFile=FALSE)$F0
    f0[f0==0] <- NA 
      #lets approximate with lower than 100 milliseconds
    f0 <- na.approx(f0, maxgap=10)
    f0diff <- rep(NA, times=length(f0)) #derivative of F0
      #identify runs of f0 and smooth
    x <- as.vector(ifelse(is.na(f0), 0, 1))
    p <- make.event.identifiers(x, 1:length(x))
    for(i in unique(p)) #smooth each run of values
    {
      if(!is.na(i))
      {
        f0[which(p==i)] <- kz(f0[which(p==i)], 3, 3)
        f0diff[which(p==i)[-1]]<-  kz(diff(f0[which(p==i)]), 3, 3)
      }
    }
    f0[is.na(f0)] <- 0 #set t0 0 again when there is an NA
    F0m <- max(f0)
    F0diffm <- max(f0diff,na.rm=TRUE)
    time_ms <- seq(1000/100, (length(f0))*1000/100, by=1000/100)
    f0d <- cbind.data.frame(time_ms, f0)
    
    ################################################merge acoustics 
    acoustics <- merge(envelope, f0d, by.x ='time_ms', by.y = 'time_ms',  all = TRUE)
    
    ######################################local peaks in env
    peaksenv <- findpeaks(acoustics$envelope)
    timepeaksenv <- acoustics$time_ms[peaksenv[,2]]
  
    #####################################peaks in F0
    peaksf0 <- findpeaks(acoustics$f0)
    timepeaksf0 <- acoustics$time_ms[peaksf0[,2]]
    #Merging motion tracking and acoustics
    merged <- merge(acoustics, mt, by.x='time_ms', by.y = 'time_ms', all=TRUE)
    merged$speed <- na.approx(merged$speed, x= merged$time_ms) #since sampling times will not perfectly align, we linearly interpolate
    merged$acc <- na.approx(merged$acc, x= merged$time_ms)     #since sampling times will not perfectly align, we linearly interpolate
    merged$envelope[is.na(merged$envelope)] <-0                #for trailing un-interpolatable values at time series edges, we replace with 0
    merged$f0[is.na(merged$f0)] <-0                            #for trailing uninterpolatable values at time series edges, we replace with 0
     #also save some timeseries
     write.csv(merged, paste0(timeseriesfol, snippet, '_TS.csv'))
   
  #global maximum in |acceleration|
  xacc <- max(abs(merged$acc),na.rm=TRUE)
   #time of max acc
  acctime <- merged$time_ms[which.max(abs(merged$acc))] #time of max acceleration
  #make dataset by filling vectors per event (i.e., in this iteration of the loop)
  acc_nearest_env_peak <- c(acc_nearest_env_peak, peaksenv[,1][which.min(timepeaksenv-acctime)])    #add the nearest peak in envelope next to peak acceleration
  acc_nearest_f0_peak <- c(acc_nearest_f0_peak, peaksf0[,1][which.min(timepeaksf0-acctime)])  #add the nearest peak F0
  duration <- c(duration, max(merged$time_ms))                        #duration of the event
  locomote_posture <- c(locomote_posture, loco)                       #annotated posture
  maxamp <- c(maxamp,max(merged$env,na.rm=TRUE))                      #global max amplitude
  maxacc <-  c(maxacc,xacc)                                           #global maximum acceleration
  synchrony <- c(synchrony, min(abs(timepeaksenv-acctime)))           #time between global max acceleration and synchrony
  maxf0 <- c(maxf0, F0m)                                              #global max F0
  individual <- c(individual, anno[,5])                               #individual (baju or fajar)
  videonames <- c(videonames, snippet)                                #add videonames
  angles <- c(angles, number_of_angles)                               #number of angles available
    }
}
#gather the dataset
data <- cbind.data.frame(videonames, duration, angles, individual,locomote_posture, acc_nearest_env_peak, acc_nearest_f0_peak,
                         maxf0, synchrony, maxacc,maxamp)
#save this dataset for input as statistical analyses
write.csv(data, paste0(proc_fol, 'pdatav3.csv'))
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
library(mpmi)                 # mutual information

##############################INFO
#CONTACT: wim.pouw@donders.ru.nl
#Version: 3
#Manuscript: Cross-modal constraints in multimodal vocalizations in Siamang (Syndactylus symphalangus) 
#Script info: This script concerns the preprocessing and aggregation of acoustic and motion tracking data, which will result in a dataset
# that is used for statistical analyses described in a different script.


#########Folder and files
parentfolder <- dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))  #what is the current folder?

audiofold  <- paste0(parentfolder, '/original_data/audio/')
mtfold  <- paste0(parentfolder, '/snippets_processed_MT_audio/motion_tracking/')
audiofiles <- list.files(audiofold) #this contains all the audio wherefore we handchecked whether the motion tracking is ok (others are excluded)
motfiles <- list.files(mtfold, pattern= '*filtered.csv')
annotations <- paste0(parentfolder, '/scripts_auto/analyses_scripts_R/annotations/')
proc_fol <- paste0(parentfolder, '/scripts_auto/analyses_scripts_R/processed_data/')
timeseriesfol <- paste0(proc_fol, 'time_series/')
#latest trained model with deeplabcut
dlcmodel <- 'DLC_resnet50_SiamangSimpleV2Jan31shuffle1_500000_filtered'

#just check what snippets we need to process
mf <- str_replace(motfiles, paste0(dlcmodel, ".csv"), "")
mf <- str_replace(mf, 'snip_gopro_orange_opp_', "")
mf <- str_replace(mf, 'snip_gopro_blue_opp_', "")
mf <- str_replace(mf, 'snip_gopro_cyan_opp_', "")
mf <- str_replace(mf, 'snip_gopro_green_opp_', "")

# we simulate a time series for random comparison
set.seed(123)
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

# compute a mutual information at a lag
calculate_mutual_info_at_lag <- function(x, y, lag) {
  if (lag >= 0) {
    x_lagged <- c(rep(NA, times=lag), x[1:(length(x)-lag)]) #x leads y
    y_lagged <- y[!is.na(x_lagged)]
    x_lagged <- x_lagged[!is.na(x_lagged)]
  } else {
    x_lagged <-  c(x[abs(lag):length(x)], rep(NA, times=abs(lag))) #x follows y
    y_lagged <- y[!is.na(x_lagged )]
    x_lagged <- x_lagged[!is.na(x_lagged)]
  }
  
  # Calculate mutual information
  mut <- cmi.pw(x_lagged, y_lagged)
  
  return(mut$mi)
}

#variables to be made 
maxamp <- maxacc <- maxdec <- synchrony <- envatmaxacc <- accatmaxF0  <-accatmaxenv <- 
  accatmaxF0<- maxf0  <- videonames <- videonamestracked<- individual <- acc_nearest_env_peak <- angles <- angles <- 
  acc_nearest_f0_peak <-  dec_nearest_env_peak <- duration <- dec_nearest_envpeak <- dec_nearest_f0_peak <- locomote_posture <-
  mutinfomax<-mutinfomaxlag <- mutinfomax_ran <- mutinfomaxlag_ran <- vector()

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
    
    #merge acoustics 
    acoustics <- merge(envelope, f0d, by.x ='time_ms', by.y = 'time_ms',  all = TRUE)
    
    #local peaks in env
    peaksenv <- findpeaks(acoustics$envelope)
    timepeaksenv <- acoustics$time_ms[peaksenv[,2]]
  
    #peaks in F0
    peaksf0 <- findpeaks(acoustics$f0)
    timepeaksf0 <- acoustics$time_ms[peaksf0[,2]]
    
  #Merging motion tracking and acoustics
   merged <- merge(acoustics, mt, by.x='time_ms', by.y = 'time_ms', all=TRUE)
   merged$speed <- na.approx(merged$speed, x= merged$time_ms) #since sampling times will not perfectly align, we linearly interpolate
   merged$acc <- na.approx(merged$acc, x= merged$time_ms)     #since sampling times will not perfectly align, we linearly interpolate
   merged$envelope[is.na(merged$envelope)] <-0                #for trailing un-interpolatable values at time series edges, we replace with 0
   merged$f0[is.na(merged$f0)] <-0                            #for trailing uninterpolatable values at time series edges, we replace with 0
   #plot(merged$envelope)
   #also save some timeseres
   write.csv(merged, paste0(timeseriesfol, snippet, '_TS.csv'))
   
   xacc <- max(abs(merged$acc),na.rm=TRUE)
   xdec <- min(merged$acc, na.rm=TRUE)
   #synchrony
   acctime <- merged$time_ms[which.max(abs(merged$acc))]
   dectime <- merged$time_ms[which.min(merged$acc)]
   maxtime <- merged$time_ms[which.max(merged$env)] #at what time was there a max in the envelope
   
    #add the nearest peak in envelope next to peak acceleration
   acc_nearest_env_peak <- c(acc_nearest_env_peak, peaksenv[,1][which.min(timepeaksenv-acctime)])
   acc_nearest_f0_peak <- c(acc_nearest_f0_peak, peaksf0[,1][which.min(timepeaksf0-acctime)])
   
    #deceleration and nearest peak
   dec_nearest_env_peak <- c(dec_nearest_env_peak, peaksenv[,1][which.min(timepeaksenv-dectime)])
   dec_nearest_f0_peak <- c(dec_nearest_f0_peak, peaksf0[,1][which.min(timepeaksf0-dectime)])

  ############################################ Calculate average mutual information
  x <- scale(merged$acc)[,1]
  y <- scale(merged$envelope)[,1]
  lents <- length(x)
  # lets also generate a random time series
  # Parameters of the original series
  mean_original <- mean(x) # close to 0, z-scaled
  sd_original <- sd(x) # close to 1, z-scaled
  acf_original <- acf(x, plot = FALSE)$acf[2]  # Autocorrelation at lag 1
  
  # Simulate a new time series with the same mean, sd, and autocorrelation
  simulated_series <- arima.sim(model = list(order = c(1, 0, 0), ar = acf_original),
                                n = lents, mean = mean_original, sd = sd_original)
  
  # loop through 1:100 lags
    # so for positive values, it means how much does acceleration predict envelope
    # so for negative values, it means how much does envelope predict acceleration
  mutvalues <- mutvalues_ran <- vector()
  lags <- seq(-20,20,1)
  for(lag in lags)
  {
    mutvalues<- c(mutvalues, calculate_mutual_info_at_lag(x,y,lag))
    mutvalues_ran <- c(mutvalues_ran, calculate_mutual_info_at_lag(simulated_series,y,lag))
  }
  
  mut <- max(mutvalues)
  mutlag <- lags[which.max(mutvalues)]*10 # each index is 10 ms
  mut_ran <- max(mutvalues_ran)
  mutlag_ran <- lags[which.max(mutvalues_ran)]*10 # each index is 10 ms

  #####################################################
  
  #make dataset
  duration <- c(duration, max(merged$time_ms))
  locomote_posture <- c(locomote_posture, loco)
  accatmaxenv <-c(accatmaxenv, merged$acc[which.max(merged$env)])
  accatmaxF0 <-c(accatmaxF0, merged$acc[which.max(merged$f0)])
  maxamp <- c(maxamp,max(merged$env,na.rm=TRUE))
  maxacc <-  c(maxacc,xacc)
  maxdec <- c(maxdec,xdec)
  synchrony <- c(synchrony, min(abs(timepeaksenv-acctime)))
  maxf0 <- c(maxf0, F0m)
  individual <- c(individual, anno[,5])
  videonames <- c(videonames, snippet)
  angles <- c(angles, number_of_angles)
  mutinfomax <- c(mutinfomax, mut)
  mutinfomaxlag <- c(mutinfomaxlag, mutlag)
  mutinfomax_ran <- c(mutinfomax_ran, mut_ran)
  mutinfomaxlag_ran <- c(mutinfomaxlag_ran, mutlag_ran)
  }
}

data <- cbind.data.frame(videonames, duration, angles, individual,locomote_posture, acc_nearest_env_peak, acc_nearest_f0_peak, 
                         dec_nearest_env_peak, dec_nearest_f0_peak,
                         maxf0, synchrony, maxacc,maxdec, maxamp, accatmaxF0, accatmaxenv, mutinfomax,mutinfomax_ran,  mutinfomaxlag, mutinfomaxlag_ran)

write.csv(data, paste0(proc_fol, 'pdatav4.csv'))
import numpy as np #arrays
import cv2 #video processing
from matplotlib import pyplot as plt #plotting
import imutils #processing tools
from scipy.signal import savgol_filter #signal processing
from scipy import signal #signal processing
from scipy import interpolate #signal processing
import statistics #some stat measures
import math #some math operations
import os #standard folder stuff
from os import listdir #standard folder stuff
from os.path import isfile, join #standard folder stuff

#########################################
#Wim Pouw (wim.pouw@donders.ru.nl): This script takes in a video and constructs a static bounding box around, around an area
#containing movement in the frame (in the upper region of the frames)
#########################################

# Define a function to interpolate to fill nan values
def fill_nan(A):
    '''
    interpolate to fill nan values
    '''
    A = np.array(A)
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = interpolate.interp1d(inds[good], A[good],bounds_error=False)
    B = np.where(np.isfinite(A),A,f(inds))
    return B.tolist()


# Set seed for random number generator
np.random.seed(42)
#color set
def fixColor(image):
    return(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))

#Main procedure that produces static bounding boxes
mypath = "../snippets_to_process/" #this is your folder with (all) your video(s)
vfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
print(vfiles)
for vid in vfiles:
    print(vid)
    video_stream = cv2.VideoCapture(mypath+vid) #load in the vide
    frames = []
    while(True):
        ret, frame = video_stream.read()
        if ret == False: #if there are no more frames, break the loop
            break
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        frames.append(gray)
        
    median = np.median(frames, axis=0).astype(dtype=np.uint8) 
    #loop through frames of the video
    writer = cv2.VideoWriter("../snippets_output_raw/testoutput"+vid[0:len(vid)-4]+".mp4", 
                           cv2.VideoWriter_fourcc(*"mp4v"), 50,(1500,1080))
    video_stream = cv2.VideoCapture(mypath+vid)
    
    #bounding boxes to save
    #capture all bounding boxes and also put them on raw frames
    framex = []
    framey = []
    framew = []
    frameh = []
    while(True):
        ret, frame = video_stream.read()
        if ret == False: #if there are no more frames, break the loop
            break
        #if ret == False: #if there are no more frames, break the loop
        #  cv2.destroyAllWindows()
        #to greyscale
        gframe = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        #remove background
        dframe = cv2.absdiff(gframe, median)
        #blurred for better performance
        blurred =cv2.GaussianBlur(dframe,(15,15),0)
        #thresholding and dilation
        ret, tframe = cv2.threshold(blurred, 100, 255, cv2.THRESH_BINARY)
        kernel = np.ones((80, 80), 'uint8')
        tframe = cv2.dilate(tframe, kernel, iterations=1)
        #contour finding
        (cnts, _) = cv2.findContours(tframe.copy(), cv2.RETR_EXTERNAL, 
                                 cv2.CHAIN_APPROX_SIMPLE)       
        #For each contour draw the bounding bos
        X = []
        Y = []
        W = []
        H = []
        for cnt in cnts:
            x,y,w,h = cv2.boundingRect(cnt)
            X.append(x)
            Y.append(y)
            W.append(w)
            H.append(h)
        if(len(Y)>0):
            index = Y.index(min(Y))
            framex.append(X[index])
            framey.append(Y[index])
            framew.append(W[index])
            frameh.append(H[index])              
            #save the wraw bounding boxes
            cv2.rectangle(frame,(X[index],Y[index]),(X[index]+W[index],Y[index]+H[index]),(0,255,0),2)
            writer.write(frame) #write to raw output, so that we have a video of the moving bounding box
        if(len(Y)==0):
            framex.append(np.nan)
            framey.append(np.nan)
            framew.append(np.nan)
            frameh.append(np.nan)  
    video_stream.release() 
    writer.release()
    
    ##############################################smooth x and y
    if sum(np.isnan(framex)) < 25:
        #interpolate intermittent NA's
        framex= fill_nan(framex)
        framey= fill_nan(framey)
        frameh= fill_nan(frameh)
        framew= fill_nan(framew)
        #fill beginning NA with nearest
        ##TO DO: fill any nas with adjacent values      
        framex = np.array(framex)
        framex[np.isnan(framex)]=framex[~np.isnan(framex)][0]
        framey = np.array(framey)
        framey[np.isnan(framey)]=framey[~np.isnan(framey)][0]
        frameh = np.array(frameh)
        frameh[np.isnan(frameh)]=frameh[~np.isnan(frameh)][0]
        framew = np.array(framew)
        framew[np.isnan(framew)]=framew[~np.isnan(framew)][0]

        #smooth
        b, a = signal.ellip(6, 0.01, 120, 0.02) 
        framex=signal.filtfilt(b, a, framex.tolist(), method="gust")
        framey=signal.filtfilt(b, a, framey.tolist(), method="gust")
        frameh=signal.filtfilt(b, a, frameh.tolist(), method="gust")
        framew=signal.filtfilt(b, a, framew.tolist(), method="gust")

        writer = cv2.VideoWriter("../snippets_output_smooth/testoutput"+vid[0:len(vid)-4]+".mp4", 
                               cv2.VideoWriter_fourcc(*"mp4v"), 50,(1500,1080))
        video_stream = cv2.VideoCapture(mypath+vid)
        framecount = 0
        
        ############################################static box
        xstatic = (int(min(framex))) #upper left
        ystatic = (int(min(framey)))#upper left
        hstatic = (int(max(framey)))+int(max(frameh))+100 #bottom right
        wstatic = (int(max(framex)))+int(max(framew))+100 #bottom right
        writer2 = cv2.VideoWriter("../snipped_to_process/snip_"+vid[0:len(vid)-4]+".mp4", 
                               cv2.VideoWriter_fourcc(*"mp4v"), 50,(wstatic-xstatic,hstatic-ystatic) )
        while(True):
            ret, frame = video_stream.read()
            if ret == False: #if there are no frames, break
                break     
            trimframe = frame[ystatic:hstatic, xstatic:wstatic]
            cv2.rectangle(frame,(xstatic,ystatic),(wstatic,hstatic),(0,255,0),2)  #draw bounding box
            writer.write(frame)
            writer2.write(trimframe)
            framecount=framecount+1
        writer2.release()
        writer.release()
    cv2.destroyAllWindows()
    video_stream.release()


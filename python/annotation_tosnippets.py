# Importing necessary libraries and functions
import os  # basic foldering functions
from moviepy.video.io.ffmpeg_tools import ffmpeg_extract_subclip  # video clipping function
import pandas  # data reading (csv)
import numpy as np  # saving data in array format
from pydub import AudioSegment  # for audio processing

# Setting up folder paths
basfolder = "../"  # base folder
processed_data = basfolder + "../../SiamangJaderpark_Processed/OpportunisticSampling_Song"  # folder for processed data
annotationdata = processed_data + '/Annotations/csvs/'  # folder for annotations
typedat = dict({"opp": "OpportunisticSampling_Song", "ran": "OpportunisticSampling_Song"})  # dictionary to map folder names
gopros = ['gopro_blue', 'gopro_cyan', 'gopro_green', 'gopro_orange', 'combined']  # list of video types
audio = ['gopros_mics']  # list of audio types
toprocess = os.listdir(processed_data + '/Annotations/csvs/')  # list of files to process
outputfol = basfolder + './snippets/'  # output folder for video clips
outputfol2 = basfolder + './snippets_processed_MT_audio/audio/'  # output folder for processed audio clips

# Video extraction settings
vidformat = ".mp4"  # video file format
addtime = 0.5  # duration to add to video clips

# Extracting video clips
it = 0  # counter for loop
for annot in toprocess:  # loop through each annotation file
    it = it + 1
    ID = annot.replace('_siamang_locomotionvocal_KePoV1.csv', '')  # extracting ID from filename
    typ = ID[0:3]  # extracting type from ID
    typfol = typedat[typ]  # getting corresponding folder name
    session = ID.split('session', 1)[1][0]  # extracting session number
    date = ID.split('_', 3)[1] + '_' + ID.split('_', 3)[2]  # extracting date
    annotationto = annotationdata + annot  # full path to annotation file
    andata = np.array(pandas.read_csv(os.path.abspath((annotationto))))  # read annotation file and store data in numpy array
    for gopro in gopros:  # loop through each video type
        videof = processed_data + '/' + date + '_processed/Session' + session + '/' + gopro + '/'  # path to video folder
        namcurvid = os.listdir(videof)[0]  # get video filename
        namcurvid = videof + namcurvid  # full path to video
        for ch in range(0, andata.shape[0]):  # loop through each row of annotation data
            begintime = andata[ch, 0]  # get start time
            endtime = andata[ch, 1]  # get end time
            annotlabel = andata[ch, 3]  # get annotation label
            ffmpeg_extract_subclip(namcurvid, begintime / 1000, (endtime / 1000) + addtime,
                                   targetname=outputfol + gopro + "_" + ID + "_" + str(ch) + "_" + str(
                                       annotlabel) + vidformat)  # extract video clip

# Audio extraction settings
it = 0  # counter for loop
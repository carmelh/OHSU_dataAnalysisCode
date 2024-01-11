# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 08:57:33 2023

Voltage Clamp Ephys
For voltage holds

Script to run through all the .abf ephys files for one biological repeat

Outputs:
    Dataframe (.csv) of holding_current_values
    avg current pre treatment
    avg current post treatment
    indicies (frame num.) of time averaged over
    figure of entire recorded current time series. Includes stimulation pulse

@author: Carmel Howe
"""

import os
import numpy as np
import pyabf
import matplotlib.pyplot as plt
import pandas as pd 
import sys
sys.path.insert(1, r'~\GitHub\OHSU_dataAnalysisCode\ephys_functions_dont_change')
import ephys_analysis_figure_funcs_dontChange as pf


''' ######## User inputs ######## '''
# insert csX_cellX folder here. You need to run this for each experimental repeat
#e.g. folder=r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\230404_trans_noIncubation_OCT-CAP_x4\cs4_cell1_trans_noIncubation_OCT-CAP_5uM_UV-irradiation_CAP_1uM'
# !!!!!!the r before the file path is very important!!!!!!!
folder=r'INSERT CSX_CELLX FOLDER HERE'
control = 'exp' # indicate whether the experiment was an experiment/ pos control ('exp') or ('neg') control.

# RUN this for voltage holds. UV holds and drug holds need the same code so just run twice 
# for the recording you are analysing
# Change the experiment type below
experiment = 'UV' # or 'drug';

""" End of user inputs """



# finds folders for the experiment type
if experiment == 'UV':
    holdFolder =folder + '\\holds_UV'   # find the folder with the UV hold recordings
elif experiment == 'drug':
    holdFolder =folder + '\\holds_drug'  # find the folder with the drug hold recordings
else:
    print('specify whether this is a UV or drug addition experiment')

filenames=np.array(os.listdir(holdFolder)) # list the files within that folder

abf = pyabf.ABF(holdFolder + '\\' + filenames[0]) # import the hold file

# assign variables from the abf file import
time = abf.sweepX    
current = abf.sweepY
voltage = abf.sweepC
cond= abf.data[2] # either the UV led pulse or drug addition step
sample_rate=abf.sampleRate # samples per second

# get the baseline current 
steadyStateCurrent=np.mean(current[0:30000]) 

# get the pulse information i.e. where the pulse starts and end
#cond is UV pulse or drug addition step

# where the UV pulse starts
# get the index where the UV pulse goes above 0.1 V i.e. a change in step
step_start = next(x for x, val in enumerate(cond)
                                  if val > 0.1) 

# postive current step amplitude 
# in V
step_value = round(cond[step_start] - cond[0])

# where the UV pulse ends
# get the index where the UV pulse goes below 0.1 V i.e. a change in step
step_end = (next(x for x, val in enumerate(cond[step_start:len(cond)])
                                  if val < cond[step_start]-0.1)) + step_start 


# end current is the average over the last 10 seconds of the pulse (UV or drug) and 20 seconds after
# save these indices
if control == 'neg':
    # this is what I used for my negative control
    startAvgIdx = step_end-(sample_rate*10) # last 10 s of the pulse
    endAvgIdx = step_end+(sample_rate*20) # 20s post end of pulse
    endCurrent = np.mean(current[:]) 
    
elif control == 'exp':
    # find the minimum value from the current trace and average across 30s
    minIdx = np.argmin(current)
    startAvgIdx = minIdx-(sample_rate*10) # 10 s before min value
    endAvgIdx = minIdx+(sample_rate*20) # 20s post min value
    endCurrent = np.mean(current[startAvgIdx:endAvgIdx]) 
    
    fig = plt.gcf()    # check this figure makes sense. if not refer to analysis guide
    plt.plot(current[startAvgIdx:endAvgIdx])
    pf.saveFigurePNG(fig,folder + r'\\figures\\' ,'maxCurrent')
    
    ## uncomment and run these lines if your endCurrent doesn't look like the correct value
    # plt.plot(current)
    # startAvgIdx = 200000  # put your start index here
    # endAvgIdx = 800000    # put your end index here
    # endCurrent = np.mean(current[startAvgIdx:endAvgIdx]) 
   
else:
    # error message if control not specified 
    print('you havent specified whether you are running a pos or neg control')


diff = endCurrent - steadyStateCurrent # difference


# saves .csv file in analysedData "holding_current_values_CONDITION.csv"
# folder
# filename
# steadyStateCurrent
# endCurrent
# diff
# condition indices
data = [[folder, filenames[0], steadyStateCurrent, endCurrent, diff, startAvgIdx, endAvgIdx]]
df = pd.DataFrame(data, columns=['folder', 'filename','steadyStateCurrent','endCurrent','deltaCurrent','condIdx1','condIdx2'])
df.to_csv(folder + r'\\analysedData\\' + 'holding_current_values_{}_pA.csv'.format(experiment))    




####### Figure ######

tPlot, axes = plt.subplots(
        nrows=2, ncols=1, sharex=True, sharey=False, 
        gridspec_kw={'height_ratios':[0.25,2]}
        )
    
axes[0].plot(time,np.transpose(cond)/5,linewidth=2,color='#6ba566')
    
axes[1].plot(time,current,linewidth=2,color='k')
axes[1].plot([0,0],[-500,-1500],linewidth=4.0,color='k')
axes[1].plot([0,20],[-1500,-1500],linewidth=4.0,color='k')
pf.noBorders(axes[0])
tPlot.set_size_inches(8,4)
pf.saveFigure(tPlot,folder+'\\'+'figures','timeSeriesExtended_noTRPV1_CAP-1uM_washOnCap_SB-20seconds_1000pA')


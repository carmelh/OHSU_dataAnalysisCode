# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 14:11:07 2023

Current Clamp Ephys
Script to run through all the .abf ephys files for one biological repeat

Outputs:
    Processed steps
    number of spikes
    access resistance
    resting membrane potential

@author: Carmel Howe
"""

import sys
sys.path.insert(1, r'Z:\Labs\Frank Lab\SHARED\000Frank lab shared\Data Analysis Scripts\Python\ephys_functions_dont_change')
import ephys_analysis_funcs_dontChange as eaf
import ephys_analysis_figure_funcs_dontChange as efig



''' ######## User inputs ######## '''
# insert csX_cellX folder here. You need to run this for each experimental repeat
#e.g. folder=r'~\Carmel\Ephys\Hippo\Data\Senktide\231003_senktide_cClamp\cs2_cell1'
folder=r'INSERT_FOLDER_HERE'

""" End of user inputs """



"""
######## current clamp steps ########
"""
# run processing of current clamp steps
output_posSteps, output_negSteps = eaf.runAllSteps(folder)
# notes: this saves the following data file in the analysedData folder:
# processedABFfiles_posSteps.npy
# and/or processedABFfiles_negSteps.npy
# for example output_posSteps data frame is composed of
# output_posSteps[0] = file names of the step files
# output_posSteps[1] = time in seconds
# output_posSteps[2] = injected current in nA
# output_posSteps[3] = recorded voltage in mV
# negSteps is the same format


eaf.numSpikes(output_posSteps,folder)
# saves the file numSpikesThreshold_py.csv
# includes the file name
# how many spikes at each current injection step
# threhold i.e. the current step where the first AP occurs


eaf.access_resistance(output_posSteps, folder)
# calculates the access resistance from the negative current step in your time trace
# saves the file access_resistance_inMOhms_py.csv
# includes the file name
# access resistance in MOhms for each current step
# average of the access res and standard deviation
# you want the deviation to be very small for each file. If it's large then you had a bad patch

eaf.rmp(output_posSteps,folder)
# calculates the resting membrane potential from the start of each current step
# saves the file rmp_py.csv (in mV)
# includes the file name
# rmp, average, and standard deviation
# again want stdev to be really small for each file



##### Analyse negative steps for input resistance if they exist #####
# saves two file:
    # file 1: negativeStepValues_in-mV_py.csv
        # steady state voltage at each current step
    # file 2: negativeStepResults_py.csv
        # linear regression results
        # slope is input resistance in Ohms     
if output_negSteps:
    eaf.inputResistance(output_negSteps, folder)
else:
    print("No negative steps")
    
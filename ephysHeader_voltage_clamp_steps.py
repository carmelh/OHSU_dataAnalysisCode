# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 14:11:07 2023

Voltage Clamp Ephys
Script to run through all the .abf ephys files for one cell from the 'vSteps' folder

Outputs:
    Processed steps
    number of spikes
    access resistance

@author: Carmel Howe
"""

import sys
sys.path.insert(1, r'~\Documents\GitHub\OHSU_dataAnalysisCode\ephys_functions_dont_change')
import ephys_analysis_funcs_dontChange as eaf

''' ######## User inputs ######## '''
# insert csX_cellX folder here. You need to run this for each experimental repeat
#e.g. folder=r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\230404_trans_noIncubation_OCT-CAP_x4\cs4_cell1_trans_noIncubation_OCT-CAP_5uM_UV-irradiation_CAP_1uM'
# !!!!!!the r before the file path is very important!!!!!!!
folder=r'INSERT CSX_CELLX FOLDER HERE'
expDrug = 'DRUG'

''' End of User inputs '''

"""
######## voltage clamp steps ########
"""

# run processing of .abf files from voltage clamp
output_vSteps = eaf.runAllSteps_voltage(folder)
# notes: this saves the following data file in the analysedData folder:
# processedABFfiles_voltageSteps.npy
# for example output_posSteps data frame is composed of
# output_posSteps[0] = file names of the step files
# output_posSteps[1] = time in seconds
# output_posSteps[2] = injected voltage in mV
# output_posSteps[3] = recorded current in pA


eaf.access_resistance_vClamp(output_vSteps, folder)
# calculates the access resistance from the negative voltage step in your time trace
# saves the file access_resistance_inMOhms_py.csv
# includes the file name
# access resistance in MOhms for each voltage step
# average of the access res and standard deviation
# you want the deviation to be very small for each file. If it's large then you had a bad patch


# process current vs, voltage clamp steps
eaf.voltage_IV(output_vSteps, folder, expDrug)
# saves the file steadyStateCurrent_perVoltageStep_in-pA_py.csv
# includes the file name
# steady state current (average of the last half of the voltage step) for each voltage step
# plots IV curve for all files in the vSteps folder








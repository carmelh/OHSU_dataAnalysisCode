# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 13:29:01 2024

@author: howeca
"""

import os
import numpy as np
import pyabf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
from scipy import stats, signal
import sys
sys.path.insert(1, r'~\GitHub\OHSU_dataAnalysisCode\ephys_functions_dont_change')
import plotting_functions as eaf
import general_functions as gf
import time


timestr = time.strftime("%y%m%d") # for saving figures and filenames with current date




# needs stepValues, linRegress_results & recorded_amplitude     
def negativeSteps_figure(data_all, folder, drug):
    df_values = pd.read_csv(folder + r'\\analysedData' + '\\' + 'negativeStepValues_in-mV_py.csv', index_col=0)
    df_lin = pd.read_csv(folder + r'\\analysedData' + '\\' + 'negativeStepResults_py.csv', index_col=0)

    
    


    
    fig = plt.gcf()      
    ax = plt.subplot(111)  
    #baseline
    plt.plot(stepValues[1], linRegress_results[0].intercept + (linRegress_results[0].slope/1e6)*stepValues[1], '--', linewidth=2, color=getExpColors('baseline'), label='Baseline')
    plt.plot(stepValues[1],recorded_amplitude[0],'o',color=getExpColors('baseline'))

    # drug
    plt.plot(stepValues[1], linRegress_results[1].intercept + (linRegress_results[1].slope/1e6)*stepValues[1],  '--', linewidth=2, color=getExpColors(drug), label=drug)
    plt.plot(stepValues[1],recorded_amplitude[1],'o',color=getExpColors(drug))
    pf.lrBorders(ax)
    fig.set_size_inches(4.5,4.5)
    plt.legend(frameon=False)
    plt.xlabel('Current (pA)')
    plt.ylabel('Voltage (mV)')
    pf.saveFigurePNG(fig,folder + r'\\figures\\','inputResistance_py')
    
    return 
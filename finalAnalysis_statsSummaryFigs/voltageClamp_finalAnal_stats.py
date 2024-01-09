# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 15:10:11 2024

different functions (not really but just organising them as such)
for different figures and statistal analyses for different
experiments I did

@author: Carmel Howe
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import sys
sys.path.insert(1, r'C:\Users\howeca\Documents\GitHub\carmel_functions')
import plotting_functions as pf
import general_functions as gf
import seaborn as sns 
import time


timestr = time.strftime("%y%m%d") # for saving figures and filenames with current date



def finalVoltageHolds_HEK_SummaryStats():
    folder = r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\finalAnalysis'
    loc_dataSummary = folder + '\\deltaCurrent_vehicle_noInc-Inc-Cap_CPZ_UV_muant_noTRPV1.csv'
    
    df = pd.read_csv(loc_dataSummary,index_col=0)
      
    
    fig = plt.figure()
    sns.boxplot(x='what',y='diff',data=df, width=0.75, showfliers = False,
                flierprops={"markersize": 8},
                boxprops={"facecolor": (.0, .0, .0, .0)},
                medianprops={"color": "DarkCyan"},)
    sns.stripplot(x='what', y='diff', data=df, alpha=0.5) 
    my_xticks = ['DMSO', 'Cap', 'Unteth', 'Teth', 'CPZ','no TRPV1','Mutant']
    x=np.array([0,1,2,3,4,5,6])
    plt.xticks(x, my_xticks,rotation=30)       
    fig.set_size_inches(5,5) 
     
    # summary stats
    summaryStats = df.groupby(['what']).agg(['mean','median']) 
    summaryStats.to_csv(folder + '\\summaryStats_{}.csv'.format(timestr))
    
    
    # ttest between dmso control and different conditions
    # cannot find a more sophisticated way to do this
    dmso = df[df['what'] == 'DMSO_irrad']
    cap = df[df['what'] == 'capWashOn']
    unteth_OCTCAP = df[df['what'] == 'washOn-OCTCAP']
    inc_OCTCAP = df[df['what'] == 'trans_incub_OCT-CAP_5uM']
    cpz = df[df['what'] == 'CPZ_preInc']
    snapMutant = df[df['what'] == 'snapMutant']
    noTRPV1 = df[df['what'] == 'noTRPV1']

    ttest_dmsoCap = stats.ttest_ind(dmso['diff'],cap['diff'])
    ttest_dmsoUnteth = stats.ttest_ind(dmso['diff'],unteth_OCTCAP['diff'])
    ttest_dmsoTeth = stats.ttest_ind(dmso['diff'],inc_OCTCAP['diff'])
    ttest_dmsoCPZ = stats.ttest_ind(dmso['diff'],cpz['diff'])
    ttest_dmsoMut = stats.ttest_ind(dmso['diff'],snapMutant['diff'])
    ttest_dmso_noTRPV1 = stats.ttest_ind(dmso['diff'],noTRPV1['diff'])
    ttest_teth_noTeth = stats.ttest_ind(unteth_OCTCAP['diff'],inc_OCTCAP['diff'])
    
    
    # send ttests to dictionry then convert to df and save
    ttestDict = {'dmsoCAP': ttest_dmsoCap, 'dmsoUnteth': ttest_dmsoUnteth, 'dmsoTeth': ttest_dmsoTeth,
            'dmsoCPZ': ttest_dmsoCPZ, 'dmsoMut': ttest_dmsoMut, 'dmso_noTRPV1': ttest_dmso_noTRPV1,
            'teth_unTeth': ttest_teth_noTeth}
    
    ttest_df = pd.DataFrame(ttestDict, index=['statistic', 'p-value'])
    ttest_df.to_csv(folder + '\\ttestResults_{}.csv'.format(timestr))
    
    return

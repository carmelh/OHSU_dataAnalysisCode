# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 14:19:11 2023

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
sys.path.insert(1, r'C:\Users\howeca\Documents\GitHub\carmel_functions')
import plotting_functions as pf
import general_functions as gf
import pylab as P

def voltage_IV(output_vSteps, folder):
    folder=r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\230405_trans_incubation-90mins_OCT-CAP_x4\cs2_cell1_trans_90mins-incubation_OCT-CAP_5uM_UV-irradiation_CAP_1uM'
    stepsFolder =folder + '\\vSteps'
    
    filenames=np.array(os.listdir(stepsFolder))
    filenames.sort()
    
    output_vSteps = runAllSteps_voltage(folder)
    stepValues = np.linspace(-100,100,11)
    
    steadyStateCurrent = np.mean(output_vSteps[3][:,:,15000:26000],axis=2)
    
    df_ssc = pd.DataFrame(data=steadyStateCurrent,columns=stepValues,index=filenames)
    df_ssc.to_csv(folder + r'\\analysedData\\' + 'steady_state_current_IV_in-pA_py.csv')


    #figure
    fig = plt.gcf()      
    ax = plt.subplot(111)  
    #baseline
    plt.plot(stepValues,np.transpose(steadyStateCurrent[0]/1000),'k',linewidth=2,label='Baseline')
    
    # unthether
    # plt.plot(stepValues,np.transpose(steadyStateCurrent[1]/1000),color='#f9a19a',linewidth=2,label='Wash on')
    # plt.plot(stepValues,np.transpose(steadyStateCurrent[2]/1000),color='#c3352b',linewidth=2,label='Irradiation')
    # plt.plot(stepValues,np.transpose(steadyStateCurrent[4]/1000),color='#611a15',linewidth=2,label='Capsaicin')
    
    #  tethered
    plt.plot(stepValues,np.transpose(steadyStateCurrent[1]/1000),color='#c3352b',linewidth=2,label='Irradiation')
    plt.plot(stepValues,np.transpose(steadyStateCurrent[3]/1000),color='#611a15',linewidth=2,label='Capsaicin')
   
    ax.spines['left'].set_position(('data', 0))
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.set_size_inches(4.5,4.5)
    plt.legend(frameon=False)
    # plt.xlabel('Voltage (mV)')
    # plt.ylabel('Current (pA)')
    pf.saveFigure(fig,folder + r'\\figures\\','IV_py_currentIn-nA')
    
    return


def voltage_hold():
    folder=r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\230915_SNAP_noTRPV1_noInc_x6\cs5_cell1_noTRPV1' # folder where your single cell data is
    capFolder = folder + '\\holds_drug'   
    
    filenames=np.array(os.listdir(capFolder))
   
    abf = pyabf.ABF(capFolder + '\\' + filenames[0]) # imports the .abf files
    
    # extracts relevant data from the abf
    time = abf.sweepX    
    current = abf.sweepY
    voltage = abf.sweepC
  #  UV_led= abf.data[3]
    drug= abf.data[2]
    
    tPlot, axes = plt.subplots(
        nrows=2, ncols=1, sharex=True, sharey=False, 
        gridspec_kw={'height_ratios':[0.25,2]}
        )
    
    axes[0].plot(time,np.transpose(drug)/5,linewidth=2,color='#6ba566')
    
    axes[1].plot(time,current,linewidth=2,color='k')
    axes[1].plot([0,0],[-500,-1500],linewidth=4.0,color='k')
    axes[1].plot([0,20],[-1500,-1500],linewidth=4.0,color='k')
    pf.noBorders(axes[0])
    pf.noBorders(axes[1])
    tPlot.set_size_inches(8,4)
    pf.saveFigure(tPlot,folder+'\\'+'figures','timeSeriesExtended_noTRPV1_CAP-1uM_washOnCap_SB-20seconds_1000pA')
    
    return

def voltage_hold_oneFile():
    folder=r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\230915_SNAP_noTRPV1_noInc_x6\cs6_cell1_washOCT-CAP' # folder where your single cell data is
    holdFolder =folder + '\\holds_longProt'   
    endFolder = folder + '\\holds_end'   
    
    filenames=np.array(os.listdir(holdFolder))
   
    abf = pyabf.ABF(holdFolder + '\\' + filenames[0]) # imports the .abf files
    
    # extracts relevant data from the abf
    time = abf.sweepX    
    current = abf.sweepY
    voltage = abf.sweepC
    UV_led= abf.data[3]
    drug= abf.data[2]
    
    
    # repeating the first bit to add on the second protocol. if you just have one then ignore this
    filenames1=np.array(os.listdir(endFolder))
    filenames1.sort()

    abf1 = pyabf.ABF(endFolder + '\\' + filenames1[0])
    
    time1 = abf1.sweepX    
    current1 = abf1.sweepY
    voltage1 = abf1.sweepC
    
    totalTimeSize=len(current)+len(current1)
    time =np.linspace(0,totalTimeSize*0.0001,totalTimeSize)
    
    fullCurrent = np.concatenate((current,current1))
    
    UV_led_pad= np.pad(UV_led,(0,len(fullCurrent)-len(UV_led)))
    drug_pad= np.pad(drug,(0,len(fullCurrent)-len(drug)))
    
    
    tPlot, axes = plt.subplots(
        nrows=2, ncols=1, sharex=True, sharey=False, 
        gridspec_kw={'height_ratios':[0.25,2]}
        )
    
   # axes[0].plot(time,np.transpose(UV_led_pad),linewidth=2,color='#a066a5')
    axes[0].plot(time,np.transpose(drug_pad)/5,linewidth=2,color='#6ba566')
    
    axes[1].plot(time,fullCurrent,linewidth=2,color='k')
    #axes[1].plot([0,0],[-500,-1500],linewidth=4.0,color='k')
    #axes[1].plot([0,20],[-1500,-1500],linewidth=4.0,color='k')
    pf.noBorders(axes[0])
    pf.noBorders(axes[1])
    tPlot.set_size_inches(8,4)
    pf.saveFigure(tPlot,folder+'\\'+'figures','timeSeriesExtended_noInc_OCT-CAP-5uM_washOnCap_SB-20seconds_1000pA')
    
    # find average values 
    
    steadyStateBaseline=np.mean(fullCurrent[0:30000])

    uvStart = next(x for x, val in enumerate(UV_led_pad)
                                  if val > 1)
    uvEnd = next(x for x, val in enumerate(UV_led_pad[uvStart:len(UV_led_pad)])
                                  if val < 1)
    uvEnd=uvEnd+uvStart
    
    uvAverage=np.mean(fullCurrent[uvStart:uvEnd])
    
    steadyStateOCTcap=np.mean(fullCurrent[600000:790000])

    capAverage = np.mean(fullCurrent[2000000:2200000])
    
    # if 3 files
    drugFolder =folder + '\\holds_drug'   
    UVfolder =folder + '\\holds_UV'   
    endFolder = folder + '\\holds_end'   
    
    
    filenames=np.array(os.listdir(drugFolder))
   
    abf = pyabf.ABF(drugFolder + '\\' + filenames[0]) # imports the .abf files
    
    # extracts relevant data from the abf
    time = abf.sweepX    
    current = abf.sweepY
    voltage = abf.sweepC
    drug= abf.data[2]
    
    
    # repeating the first bit to add on the second protocol. if you just have one then ignore this
    filenames1=np.array(os.listdir(UVfolder))
    filenames1.sort()

    abf1 = pyabf.ABF(UVfolder + '\\' + filenames1[0])
    
    time1 = abf1.sweepX    
    current1 = abf1.sweepY
    voltage1 = abf1.sweepC
    UV_led= abf1.data[2]
    
    
    filenames2=np.array(os.listdir(endFolder))
    abf2 = pyabf.ABF(endFolder + '\\' + filenames2[0])
    
    time2 = abf2.sweepX    
    current2 = abf2.sweepY
    voltage2 = abf2.sweepC
    
    totalTimeSize=len(current)+len(current1)+len(current2)
    time =np.linspace(0,totalTimeSize*0.0001,totalTimeSize)
    
    fullCurrent = np.concatenate((current,current1,current2))
    
    UV_led_pad= np.pad(UV_led,(0,len(fullCurrent)-len(UV_led)))
    drug_pad= np.pad(drug,(0,len(fullCurrent)-len(drug)))
    
    
    tPlot, axes = plt.subplots(
        nrows=2, ncols=1, sharex=True, sharey=False, 
        gridspec_kw={'height_ratios':[0.25,2]}
        )
    
    #axes[0].plot(time,np.transpose(UV_led_pad),linewidth=2,color='#a066a5')
    axes[0].plot(time,np.transpose(drug_pad)/5,linewidth=2,color='#6ba566')
    
    axes[1].plot(time,fullCurrent,linewidth=2,color='k')
    axes[1].plot([0,0],[-500,-1500],linewidth=4.0,color='k')
    axes[1].plot([0,40],[-1500,-1500],linewidth=4.0,color='k')
    pf.noBorders(axes[0])
    pf.noBorders(axes[1])
    tPlot.set_size_inches(8,4)
    pf.saveFigure(tPlot,folder+'\\'+'figures','timeSeriesExtended_OCT-CAP_washOn_1uM_UV_SB-40seconds_1000pA')
    
    
    fs=10000
    last10 = len(fullCurrent) - 10*fs  #last 10 seconds

    
    steadyStateBaseline=np.mean(current[0:30000])
    steadyStateOCTcap=np.mean(current[50000:150000])
#    steadyStateCap=np.mean(fullCurrent[2500000:2700000])
    endCurrent = np.mean(fullCurrent[last10:len(fullCurrent)])
    
tPlot, axes = plt.subplots(
    nrows=2, ncols=1, sharex=True, sharey=False, 
    gridspec_kw={'height_ratios':[0.25,2]}
    )

#axes[0].plot(time,np.transpose(UV_led_pad),linewidth=2,color='#a066a5')
axes[0].plot(time,np.transpose(drug_pad)/5,linewidth=2,color='#6ba566')

axes[1].plot(time,fullCurrent,linewidth=2,color='k')
    # quickly done tide later
    
    finalAnalysisFolder=r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\finalAnalysis'
    
    
    # for controls. baseline + control e.g. UV
    df = pd.read_csv(finalAnalysisFolder + '\\holding_current_values_trans_CAP_1uM.csv',index_col=0)
    

    baseline = np.array(df['startCurrent'])
    cap = np.array(df['CapWashOn'])
    
    data = np.column_stack((baseline,cap))
    data = data/1e3
    
    boxprops = dict(linewidth=2)
    medianprops = dict(linewidth=2.5,color='DarkCyan')
    flierprops = dict(markersize=10,linestyle='none')
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for repeat in range(len(data)):
         plt.plot([1,2],[data[repeat,0],data[repeat,1]],linewidth=2,color='#bdccd4')
    
    bp=plt.boxplot(data,boxprops=boxprops,medianprops=medianprops,flierprops=flierprops, widths = 0.4)
    my_xticks = ['Baseline', 'CAP']
    x=np.array([1,2])
    plt.xticks(x, my_xticks)
    
    for whisker in bp['whiskers']:
         whisker.set(linewidth=3)
    for cap in bp['caps']:
         cap.set(linewidth=3)  
         
    for i in range(len(data[1])):
          y = data[:,i]
          x = np.random.normal(1+i, 0.04, size=len(y))
          P.plot(x, y, 'r.', markersize=14, alpha=0.6)
     
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    plt.ylabel('Current (nA)')
    # plt.xlim((0.5,3.1))
    plt.tight_layout()
    fig.set_size_inches(2,4.5) 
    pf.saveFigure(fig,finalAnalysisFolder + r'\\figures\\','holdingCurrent_CAP_1uM_poired_230629')

     # stats
     ttest_results = stats.ttest_rel(data[:,0],data[:,1])
     
     stats_baseline = gf.get_stats(data[:,0])
     stats_drug = gf.get_stats(data[:,1])


     tPlot, axes = plt.subplots(
         nrows=2, ncols=1, sharex=True, sharey=False, 
         gridspec_kw={'height_ratios':[0.25,2]}
         )
     
     axes[0].plot(np.transpose(drug[50000:150000])/5,linewidth=2,color='#6ba566')
     
     axes[1].plot(current,linewidth=2,color='k')
     
     
     df_stats = pd.DataFrame(data=ttest_results)
     df_stats.insert(0, "What", ["test stat","p-value"], True)
     df_stats.to_csv(finalAnalysisFolder + r'\\stats_threshold_trans_capWashOn-1uM_control_230629.csv')
    return



folder = r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\231201_trans_concs_1-25uM_x3\cs1_cell3'
def voltage_hold(folder):
    holdFolder =folder + '\\holds_UV'   
    
    
    filenames=np.array(os.listdir(holdFolder))
   
    abf = pyabf.ABF(holdFolder + '\\' + filenames[0])
    
    time = abf.sweepX    
    current = abf.sweepY
    voltage = abf.sweepC
    UV_led= abf.data[2]

    
    fs=10000
    last10 = len(current) - 10*fs  #last 10 seconds

    
    steadyStateCurrent=np.mean(current[0:30000])
    endCurrent = np.mean(current[479844:779844])
    # specifgy
    #endCurrent = np.mean(current[0.6e6:0.7e6])
    
    diff = endCurrent - steadyStateCurrent
    
    
    holdFolder =folder + '\\holds_drug'  
    
    filenames=np.array(os.listdir(holdFolder))
    filenames.sort()
   
    abf = pyabf.ABF(holdFolder + '\\' + filenames[0])

    abf = pyabf.ABF(holdFolder + '\\' + filenames[3])

    
    
    time = abf.sweepX    
    current = abf.sweepY
    voltage = abf.sweepC
    
    abf = pyabf.ABF(holdFolder + '\\' + filenames[1])
    current1 = abf.sweepY

    time =np.linspace(0,(len(current)+len(current1))*0.0001,len(current)+len(current1))
    
    fullCurrent = np.concatenate((current,current1))
    
    fs=10000
    last10 = len(current) - 10*fs  #last 10 seconds
    
    steadyStateCurrent=np.mean(current[0:30000])
    endCurrent = np.mean(current[last10:len(current)])
    diff = endCurrent - steadyStateCurrent
    
    endCurrent = np.mean(fullCurrent[2500000:3000000])    
    
    # ignore for now
    finalAnalysisFolder = r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\finalAnalysis'
    df = pd.read_csv(finalAnalysisFolder + '\\holding_current_values_in-pA_py.csv',index_col=0)
    
    df.loc[len(df.index)] = [filenames[0],steadyStateCurrent, endCurrent, diff] 
    df.to_csv(finalAnalysisFolder + '\\holding_current_values_in-pA_py.csv')
    
    
    # df['filename'] =filenames[0]
    # df['Start current'] = steadyStateCurrent
    # df['Last 10secs'] = endCurrent
    # df['Difference'] = diff
    
    tPlot, axes = plt.subplots(
        nrows=2, ncols=1, sharex=True, sharey=False, 
        gridspec_kw={'height_ratios':[0.25,2]}
        )
    
    axes[0].plot(time,np.transpose(drug),linewidth=2,color='#a066a5')
    # axes[0].plot([0,0],[0,2],linewidth=4.0,color='k')
    
    axes[1].plot(time,np.transpose(current),linewidth=2,color='k')
    axes[1].plot([0,0],[-500,-1500],linewidth=4.0,color='k')
    axes[1].plot([0,20],[-1500,-1500],linewidth=4.0,color='k')
    pf.noBorders(axes[0])
    pf.noBorders(axes[1])
    # plt.ylim((-5500,0))
    tPlot.set_size_inches(8,4)
    pf.saveFigure(tPlot,folder+'\\'+'figures','timeSeries_capWashOn_SB-20seconds_1000pA')
    
    return

def voltage_hold_finalAnal(folder):
    finalAnalysisFolder=r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\finalAnalysis'
    
    
    # for controls. baseline + control e.g. UV
    df = pd.read_csv(finalAnalysisFolder + '\\holding_current_values_trans_OCT-CAP_incubated_diffConcs.csv',index_col=0)
    
    concArray = np.array(df['deltaCurrent'])
    conc125 = concArray[0:1]
    conc25 = concArray[1:4]
    conc5 = concArray[5:11]
    
    conc125_pad= np.pad(conc125,(0,len(conc5)-len(conc125)), 'constant', constant_values=np.nan)
    conc25_pad =  np.pad(conc25,(0,len(conc5)-len(conc25)), 'constant', constant_values=np.nan)

    data = np.column_stack((conc125_pad,conc25_pad,conc5))
    # data = data/1e3
    
    boxprops = dict(linewidth=2)
    medianprops = dict(linewidth=2.5,color='DarkCyan')
    flierprops = dict(markersize=10,linestyle='none')
 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for repeat in range(len(data)):
         plt.plot([1,2],[data[repeat,0],data[repeat,1]],linewidth=2,color='#bdccd4')

    bp=plt.boxplot(data,boxprops=boxprops,medianprops=medianprops,flierprops=flierprops, widths = 0.4)
    my_xticks = ['Baseline', 'UV Cont']
    x=np.array([1,2])
    plt.xticks(x, my_xticks)
    
    for whisker in bp['whiskers']:
         whisker.set(linewidth=3)
    for cap in bp['caps']:
         cap.set(linewidth=3)  
         
    for i in range(len(data[1])):
          y = data[:,i]
          x = np.random.normal(1+i, 0.04, size=len(y))
          P.plot(x, y, 'r.', markersize=14, alpha=0.6)
     
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    plt.ylabel('Current (pA)')
    # plt.xlim((0.5,3.1))
    plt.tight_layout()
    fig.set_size_inches(4.5,4.5) 
    pf.saveFigure(fig,finalAnalysisFolder + r'\\figures\\','holdingCurrent_UV-cont_poired_230516')
    
    
    # stats
    ttest_results = stats.ttest_rel(data[:,0],data[:,1])
    
    stats_baseline = gf.get_stats(data[:,0])
    stats_UV = gf.get_stats(data[:,1])
    stats_drug = gf.get_stats(data[:,2])

    
    
    df_stats = pd.DataFrame(data=ttest_results)
    df_stats.insert(0, "What", ["test stat","p-value"], True)
    df_stats.to_csv(finalAnalysisFolder + r'\\stats_threshold_trans_OCT-CAP_incubated_5uM_230428.csv')
    
    
    
    
    # for incubated 2 things i.e. UV then cap
    df = pd.read_csv(finalAnalysisFolder + '\\holding_current_values_UV-conts.csv',index_col=0)
    
    baseline_UV = np.array(df['startCurrent_UV'])
    octCap_UV = np.array(df['last10secs_UV'])
    cap = np.array(df['last10secs_CAP'])

    
    data = np.column_stack((baseline_UV,octCap_UV,cap))
    data = data/1e3
    
    boxprops = dict(linewidth=2)
    medianprops = dict(linewidth=2.5,color='DarkCyan')
    flierprops = dict(markersize=10,linestyle='none')
 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for repeat in range(len(data)):
         plt.plot([1,2],[data[repeat,0],data[repeat,1]],linewidth=2,color='#bdccd4')
         plt.plot([2,3],[data[repeat,1],data[repeat,2]],linewidth=2,color='#bdccd4')

    bp=plt.boxplot(data,boxprops=boxprops,medianprops=medianprops,flierprops=flierprops, widths = 0.4)
    my_xticks = ['Baseline', 'OCT-Cap', 'Capcasisin']
    x=np.array([1,2,3])
    plt.xticks(x, my_xticks)
    
    for whisker in bp['whiskers']:
         whisker.set(linewidth=3)
    for cap in bp['caps']:
         cap.set(linewidth=3)  
         
    for i in range(len(data[1])):
          y = data[:,i]
          x = np.random.normal(1+i, 0.04, size=len(y))
          P.plot(x, y, 'r.', markersize=14, alpha=0.6)
     
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    plt.ylabel('Current (nA)')
    # plt.xlim((0.5,3.1))
    plt.tight_layout()
    fig.set_size_inches(4.5,4.5) 
    pf.saveFigure(fig,finalAnalysisFolder + r'\\figures\\','holdingCurrent_OCT-CAP_5uM_90mins_incub_1uM-cap_poired_230428')
    
    
    # stats
    ttest_results = stats.ttest_rel(data[:,0],data[:,1])
    
    stats_baseline = gf.get_stats(data[:,0])
    stats_UV = gf.get_stats(data[:,1])
    stats_drug = gf.get_stats(data[:,2])

    
    
    df_stats = pd.DataFrame(data=ttest_results)
    df_stats.insert(0, "What", ["test stat","p-value"], True)
    df_stats.to_csv(finalAnalysisFolder + r'\\stats_threshold_trans_OCT-CAP_incubated_5uM_230428.csv')
    
    
    # for non incubated 3 things i.e. wash on probe UV then cap
    df = pd.read_csv(finalAnalysisFolder + '\\holding_current_values_OCT-CAP_trans_noIncubation_5uM.csv',index_col=0)

    
    baseline_UV = np.array(df['startCurrent_washProbe'])
    octCap_wash = np.array(df['last10secs_washProbe'])
    octCap_UV = np.array(df['last10secs_UV'])
    cap = np.array(df['last10secs_CAP'])

    data = np.column_stack((baseline_UV,octCap_wash,octCap_UV,cap))
    data = data/1e3
    
    boxprops = dict(linewidth=2)
    medianprops = dict(linewidth=2.5,color='DarkCyan')
    flierprops = dict(markersize=10,linestyle='none')
 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for repeat in range(len(data)):
         plt.plot([1,2],[data[repeat,0],data[repeat,1]],linewidth=2,color='#bdccd4')
         plt.plot([2,3],[data[repeat,1],data[repeat,2]],linewidth=2,color='#bdccd4')
         plt.plot([3,4],[data[repeat,2],data[repeat,3]],linewidth=2,color='#bdccd4')


    bp=plt.boxplot(data,boxprops=boxprops,medianprops=medianprops,flierprops=flierprops, widths = 0.4)
    my_xticks = ['Baseline', 'Wash On', 'UV', 'Capcasisin']
    x=np.array([1,2,3,4])
    plt.xticks(x, my_xticks)
    
    for whisker in bp['whiskers']:
         whisker.set(linewidth=3)
    for cap in bp['caps']:
         cap.set(linewidth=3)  
         
    for i in range(len(data[1])):
          y = data[:,i]
          x = np.random.normal(1+i, 0.04, size=len(y))
          P.plot(x, y, 'r.', markersize=14, alpha=0.6)
     
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    plt.ylabel('Current (nA)')
    # plt.xlim((0.5,3.1))
    plt.ylim((-5.2,0.20))
    plt.tight_layout()
    fig.set_size_inches(4.5,4.5) 
    pf.saveFigure(fig,finalAnalysisFolder + r'\\figures\\','holdingCurrent_OCT-CAP_5uM_90mins_noIncub_1uM-cap_poired_230428')
    
    
    
    ttest_base_wash = stats.ttest_rel(data[:,0],data[:,1])
    ttest_wash_UV = stats.ttest_rel(data[:,1],data[:,2])
    ttest_UV_drug = stats.ttest_rel(data[:,2],data[:,3])

    ttest_results = np.column_stack((ttest_base_wash,ttest_wash_UV,ttest_UV_drug))
    
    stats_baseline = gf.get_stats(data[:,0])
    stats_washOn = gf.get_stats(data[:,1])
    stats_UV = gf.get_stats(data[:,2])
    stats_drug = gf.get_stats(data[:,3])

    
    
    df_stats = pd.DataFrame(data=ttest_results)
    df_stats.insert(0, "What", ["test stat","p-value"], True)
    df_stats.to_csv(finalAnalysisFolder + r'\\stats_threshold_trans_OCT-CAP_notIncubated_5uM_230428.csv')
    
    
    
    return


def vSteps_finalFigure():
    finalAnalysisFolder=r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\finalAnalysis'
    df_noIncubation = pd.read_csv(finalAnalysisFolder + '\\vSteps_fullData_noIncubation_averages.csv')

    steps=np.array(df_noIncubation['vSteps'])
    baseline_ni = np.array(df_noIncubation['Baseline'])
    wash_ni=np.array(df_noIncubation['washOn'])
    UV_ni = np.array(df_noIncubation['Irradiation'])
    CAP_ni = np.array(df_noIncubation['CAP'])

    df_incubation = pd.read_csv(finalAnalysisFolder + '\\vSteps_fullData_incubation_averages.csv')

    baseline_i = np.array(df_incubation['Baseline'])
    UV_i = np.array(df_incubation['Irradiation'])
    CAP_i = np.array(df_incubation['CAP'])
    
    
    #figure
    fig = plt.gcf()      
    ax = plt.subplot(111)  
    # plt.plot(steps,np.transpose(baseline_ni),'k',linewidth=2,label='Non-Tethered Baseline')
    plt.plot(steps,np.transpose(baseline_i),'--',color='k',linewidth=2,label='Tethered Baseline')   
    
    # plt.plot(steps,np.transpose(wash_ni),'#a53a6e',linewidth=2,label='Non-Tethered Wash On')
    
    # plt.plot(steps,np.transpose(UV_ni),'#f653a6',linewidth=2,label='Non-Tethered Irradiation')
    plt.plot(steps,np.transpose(UV_i),'--',color='#f653a6',linewidth=2,label='Tethered Irradiation')

    # plt.plot(steps,np.transpose(CAP_ni),'#3a6ea5',linewidth=2,label='Non-Tethered Capsaicin')
    plt.plot(steps,np.transpose(CAP_i),'--',color='#3a6ea5',linewidth=2,label='Tethered Capsaicin')


    
    ax.spines['left'].set_position(('data', 0))
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.set_size_inches(4.5,4.5)
    plt.legend(frameon=False)
    pf.saveFigure(fig,finalAnalysisFolder + r'\\figures\\','IV_tetheredd_averages_currentIn-pA')
    
    
    return

def numSpikes_figure(folder):
    df = pd.read_csv(folder + r'\\analysedData\\' + 'numSpikesThreshold_py.csv',index_col=0)
    
    y = np.array(df)
    stepValues = np.array(df.columns[0:len(y[1])-1])
    num_spikes= y[:,0:len(y[1])-1]
    
    
    colors = plt.cm.winter(np.linspace(0,1,len(num_spikes)))

    fig = plt.gcf()      
    ax = plt.subplot(111)  
    # set y axis to be integer values only
    for axis in [ax.yaxis]:
        axis.set_major_locator(ticker.MaxNLocator(integer=True))
    # for sweepNumber in range(len(num_spikes)):
    #     plt.plot(stepValues,np.transpose(num_spikes[sweepNumber]), color=colors[sweepNumber], linewidth=2, alpha=.8)
    plt.plot(stepValues,num_spikes[0], linewidth=2,color='k',label='Baseline')
    plt.plot(stepValues,num_spikes[len(num_spikes)-1], linewidth=2,color='red',label='Senktide')
   
    pf.lrBorders(ax)
    fig.set_size_inches(4.5,4.5)
    plt.xlabel('Current (pA)')
    plt.ylabel('Num. Spikes')
    ax.set_yticks([0,5,10,15,20])
    ax.set_xticks([0,5,10,15])
    ax.legend(frameon=False)
    pf.saveFigure(fig,folder + r'\\figures\\','currentInjection_vs_numSpikes_firstLast')
    
    return




def finalNumSpikes():
    folder = r'Z:\Labs\Frank Lab\Carmel\Ephys\Hippo'
    loc_dataSummary = folder + '\\ephysDataSummary-ayo.csv'
    
    df = pd.read_csv(loc_dataSummary)
 
    drug = 'CAP'

    df_drug = df.loc[df['Cond 1'] == drug]   
    
    baselineSpikes = []
    baselineThreshold = []
    drugSpikes=[]
    drugThreshold =[]
    for ind in df_drug.index:
        trialFolder = df_drug['Folder'][ind]
        print(trialFolder)
        
        df_ns = pd.read_csv(trialFolder + '\\analysedData\\numSpikesThreshold_py.csv')
        
        baselineIdx =  str(df_drug['fileNam_0'][ind])
        drugIdx =  str(df_drug['fileName_1'][ind])
        
        baseline_row = df_ns[df_ns['Unnamed: 0'].str.contains(baselineIdx)]
        baseline_row = np.array(baseline_row) 
        baselineSpikes.append(baseline_row[:,0:len(baseline_row[0])-1])
        
        baselineThreshold.append(baseline_row[:, [0, len(baseline_row[0])-1]]) 
        
        drug_row = df_ns[df_ns['Unnamed: 0'].str.contains(drugIdx)]
        drug_row = np.array(drug_row) 
        drugSpikes.append(baseline_row[:,0:len(drug_row[0])-1])
 
        drugThreshold.append(drug_row[:, [0, len(drug_row[0])-1]]) 
     
    
    # TRESHOLD figure & stats
    baselineThreshold = np.squeeze(np.array(baselineThreshold))
    drugThreshold = np.squeeze(np.array(drugThreshold))
    
    data = np.column_stack((baselineThreshold[:,1],drugThreshold[:,1]))
    
    boxprops = dict(linewidth=2)
    medianprops = dict(linewidth=2.5,color='DarkCyan')
    flierprops = dict(markersize=10,linestyle='none')
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for repeat in range(len(data)):
        plt.plot([1,2],[data[repeat,0],data[repeat,1]],linewidth=2,color='#bdccd4')
    
    bp=plt.boxplot(data,boxprops=boxprops,medianprops=medianprops,flierprops=flierprops, widths = 0.4)
    my_xticks = ['Baseline', 'Capcasisin']
    x=np.array([1,2])
    plt.xticks(x, my_xticks)
    
    for whisker in bp['whiskers']:
        whisker.set(linewidth=3)
    for cap in bp['caps']:
        cap.set(linewidth=3)  
        
    for i in range(len(data[1])):
         y = data[:,i]
         x = np.random.normal(1+i, 0.04, size=len(y))
         P.plot(x, y, 'r.', markersize=14, alpha=0.6)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    plt.ylabel('Threshold (pA)')
   # plt.xlim((0.5,3.1))
    plt.tight_layout()
    fig.set_size_inches(4.5,4.5) 
    pf.saveFigure(fig,folder + r'\\resultsCapcasisin\\','threshold_1uM-cap_poired_230412')
    
    ttest_results = stats.ttest_rel(data[:,0],data[:,1])
    
    stats_baseline = gf.get_stats(data[:,0])
    stats_drug = gf.get_stats(data[:,1])
    
    
    
    df_stats = pd.DataFrame(data=ttest_results)
    df_stats.insert(0, "What", ["test stat","p-value"], True)
    df_stats.to_csv(folder + r'\\resultsCapcasisin\\' + 'stats_threshold_1uM-cap_230412.csv')
    
    

    return

def voltage_hold_finalAnal_diffConcs(folder):
    finalAnalysisFolder=r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\finalAnalysis'
    
    
    # for controls. baseline + control e.g. UV
    df = pd.read_csv(finalAnalysisFolder + '\\holding_current_values_trans_OCT-CAP_incubated_diffConcs.csv',index_col=0)
    
    concArray = np.array(df['deltaCurrent'])
    conc125 = concArray[0:1]
    conc25 = concArray[1:4]
    conc5 = concArray[5:11]
    
    conc125_pad= np.pad(conc125,(0,len(conc5)-len(conc125)), 'constant', constant_values=np.nan)
    conc25_pad =  np.pad(conc25,(0,len(conc5)-len(conc25)), 'constant', constant_values=np.nan)

    data = np.column_stack((conc125_pad,conc25_pad,conc5))
    
    fig = plt.figure()
    sns.boxplot(data=data, width=0.75, showfliers = False,
                flierprops={"markersize": 8},
                boxprops={"facecolor": (.0, .0, .0, .0)},
                medianprops={"color": "DarkCyan"},)
    for i in range(len(data[1])):
           y = data[:,i]
           x = np.random.normal(0+i, 0.04, size=len(y))
           P.plot(x, y, 'r.', markersize=12, alpha=0.6)
    my_xticks = ['1.25 uM', '2.5 uM','5 uM']
    x=np.array([0,1,2])
    plt.xticks(x, my_xticks,rotation=30)       
    fig.set_size_inches(5,5) 
    pf.saveFigure(fig,finalAnalysisFolder + r'\\figures\\','trans_OCT-CAP_incubated_diffConcs')
    
    return


def ca2ImagingFinalSumarryFigure():
    loc_dataSummary =r'Z:\Labs\Frank Lab\Carmel\Imaging\Confocal_olympusFV1200\Live\OCTCAP_HEK_jRGECO\HEK_jRGECOimaging_finalAnalysis\avergaeEachTrial_btw70and90seconds_231130.csv'
    df = pd.read_csv(loc_dataSummary)

    dmso = np.array(df['DMSO'])
    norm = np.array(df['Norm'])
    Trpv1_neg = np.array(df['Trpv1_neg'])
    SNAPmut = np.array(df['SNAPmut'])
    
    data = np.column_stack((dmso,norm,Trpv1_neg,SNAPmut))
    
    
    fig = plt.figure()
    sns.boxplot(data=data, width=0.75, showfliers = False,
                flierprops={"markersize": 8},
                boxprops={"facecolor": (.0, .0, .0, .0)},
                medianprops={"color": "DarkCyan"},)
    for i in range(len(data[1])):
           y = data[:,i]
           x = np.random.normal(0+i, 0.04, size=len(y))
           P.plot(x, y, 'r.', markersize=12, alpha=0.6)
    my_xticks = ['DMSO', 'Teth','no TRPV1','Mutant']
    x=np.array([0,1,2,3])
    plt.xticks(x, my_xticks,rotation=30)       
    fig.set_size_inches(5,5) 
    pf.saveFigure(fig,folder + r'\\figures\\','deltaCa2_all_231110')
    
    
    
    tethNorm=stats.shapiro(data[:,0])
    
    
    ttest_baselineTeth = stats.ttest_ind(data[:,0],data[:,1],nan_policy='omit')
    ttest_baselineNoTRPV1 = stats.ttest_ind(data[:,0],data[:,2],nan_policy='omit')
    ttest_baselineSNAPmut = stats.ttest_ind(data[:,0],data[:,3],nan_policy='omit')


    stats_baseline = gf.get_stats(data[0:3,0])
    stats_Cap = gf.get_stats(data[0:6,1])
    stats_unteth = gf.get_stats(data[0:4,2])
    stats_teth = gf.get_stats(data[:,3])
    stats_cpz = gf.get_stats(data[0:4,4])
    stats_noTRPV1 = gf.get_stats(data[0:3,5])
    stats_mutant = gf.get_stats(data[0:4,6])

    ttest_results = np.column_stack((ttest_baselineCap,ttest_baselineUnteth))
    
    df_stats = pd.DataFrame(data=ttest_results)
    df_stats.insert(0, "What", ["test stat","p-value"], True)
    df_stats.to_csv(folder + r'\\resultsCapcasisin\\' + 'stats_threshold_1uM-cap_230412.csv')
    
    return



# messy needs to be fixed
def finalCurrentHoldsHEK_oneFigure():
    folder = r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\finalAnalysis'
    loc_dataSummary = folder + '\\deltaCurrent_vehicle_noInc-Inc-Cap_CPZ_UV_muant_noTRPV1.csv'
    
    df = pd.read_csv(loc_dataSummary,index_col=0)
      
    
    import seaborn as sns
    
    fig = plt.figure()
    sns.boxplot(x='what',y='diff',data=df, width=0.75, showfliers = False,
                flierprops={"markersize": 8},
                boxprops={"facecolor": (.0, .0, .0, .0)},
                medianprops={"color": "DarkCyan"},)
    for i in range(len(data[1])):
           y = data[:,i]
           x = np.random.normal(0+i, 0.04, size=len(y))
           P.plot(x, y, 'r.', markersize=12, alpha=0.6)
    my_xticks = ['DMSO', 'Cap', 'Unteth', 'Teth', 'CPZ','no TRPV1','Mutant']
    x=np.array([0,1,2,3,4,5,6])
    plt.xticks(x, my_xticks,rotation=30)       
    fig.set_size_inches(5,5) 
    
    
    
    arrayValues = np.array(df['diff'])
    baseline = arrayValues[0:3]    
    capWashOn = arrayValues[3:9] 
    washOnOCTCAP = arrayValues[9:16] 
    transIncOCTCAP = arrayValues[16:23] 
    CPZ_preInc = arrayValues[23:27] 
    snapMutant = arrayValues[27:31] 
    noTRPV1 = arrayValues[31:34] 

    baseline_pad= np.pad(baseline,(0,len(transIncOCTCAP)-len(baseline)), 'constant', constant_values=np.nan)
    capWashOn_pad =  np.pad(capWashOn,(0,len(transIncOCTCAP)-len(capWashOn)), 'constant', constant_values=np.nan)
    washOnOCTCAP_pad =  np.pad(washOnOCTCAP,(0,len(transIncOCTCAP)-len(washOnOCTCAP)), 'constant', constant_values=np.nan)
    CPZ_preInc_pad =  np.pad(CPZ_preInc,(0,len(transIncOCTCAP)-len(CPZ_preInc)), 'constant', constant_values=np.nan)
    snapMutant_pad =  np.pad(snapMutant,(0,len(transIncOCTCAP)-len(snapMutant)), 'constant', constant_values=np.nan)
    noTRPV1_pad =  np.pad(noTRPV1,(0,len(transIncOCTCAP)-len(noTRPV1)), 'constant', constant_values=np.nan)

    data = np.column_stack((baseline_pad,capWashOn_pad,washOnOCTCAP_pad,transIncOCTCAP,CPZ_preInc_pad,noTRPV1_pad,snapMutant_pad))
    data = data/1e3
    
    import seaborn as sns
    
    fig = plt.figure()
    sns.boxplot(data=data, width=0.75, showfliers = False,
                flierprops={"markersize": 8},
                boxprops={"facecolor": (.0, .0, .0, .0)},
                medianprops={"color": "DarkCyan"},)
    for i in range(len(data[1])):
           y = data[:,i]
           x = np.random.normal(0+i, 0.04, size=len(y))
           P.plot(x, y, 'r.', markersize=12, alpha=0.6)
    my_xticks = ['DMSO', 'Cap', 'Unteth', 'Teth', 'CPZ','no TRPV1','Mutant']
    x=np.array([0,1,2,3,4,5,6])
    plt.xticks(x, my_xticks,rotation=30)       
    fig.set_size_inches(5,5) 
    pf.saveFigure(fig,folder + r'\\figures\\','deltaCurrent_all_230920')
   

    data = data*1e3
    
    ttest_baselineCap = stats.ttest_ind(data[0:3,0],data[0:6,1])
    ttest_baselineUnteth = stats.ttest_ind(data[0:3,0],data[0:4,2])
    ttest_baselineTeth = stats.ttest_ind(data[0:3,0],data[:,3])
    ttest_baselineCPZ = stats.ttest_ind(data[0:3,0],data[0:4,4])
    ttest_baselineMut = stats.ttest_ind(data[0:3,0],data[0:4,6])
    ttest_baseline_noTRPV1 = stats.ttest_ind(data[0:3,0],data[0:3,5])
    ttest_teth_noTeth = stats.ttest_ind(data[:,3],data[0:4,2])


    stats_baseline = gf.get_stats(data[0:3,0])
    stats_Cap = gf.get_stats(data[0:6,1])
    stats_unteth = gf.get_stats(data[0:4,2])
    stats_teth = gf.get_stats(data[:,3])
    stats_cpz = gf.get_stats(data[0:4,4])
    stats_noTRPV1 = gf.get_stats(data[0:3,5])
    stats_mutant = gf.get_stats(data[0:4,6])

    ttest_results = np.column_stack((ttest_baselineCap,ttest_baselineUnteth))
    
    df_stats = pd.DataFrame(data=ttest_results)
    df_stats.insert(0, "What", ["test stat","p-value"], True)
    df_stats.to_csv(folder + r'\\resultsCapcasisin\\' + 'stats_threshold_1uM-cap_230412.csv')
    
    return



folder=r'Z:\Labs\Frank Lab\Carmel\Ephys\Hippo\230426_cult8_OCT-Cap_wash-on_UV-irrad_cap_cClamp_x4\cs4_cell2_OCT-Cap-1uM_wash-on_UV-irrad_cap-1uM\clampfit_anal'
def currentClamp_timeseriesFigure():
    hold_folder = folder + '\\' + r'currentHolds_drug'
    filenames=np.array(os.listdir(hold_folder))
   
    abf = pyabf.ABF(hold_folder + '\\' + filenames[0])
    
    time = abf.sweepX    
    current = abf.sweepY
    voltage = abf.sweepC
    UV_led= abf.data[2]
    
    
    tPlot, axes = plt.subplots(
        nrows=2, ncols=1, sharex=True, sharey=False, 
        gridspec_kw={'height_ratios':[0.25,2]}
        )
    
    axes[0].plot(time,np.transpose(UV_led),linewidth=2,color='#a066a5')
    # axes[0].plot([0,0],[0,2],linewidth=4.0,color='k')
    
    axes[1].plot(time,np.transpose(current),linewidth=2,color='k')
    axes[1].plot([0,0],[-65,-55],linewidth=4.0,color='k')
    axes[1].plot([0,10],[-65,-65],linewidth=4.0,color='k')
    pf.noBorders(axes[0])
    pf.noBorders(axes[1])
    tPlot.set_size_inches(8,4)
    pf.saveFigure(tPlot,folder+'\\'+'figures','timeSeries_senk_100nM_SB-10seconds_10mV')
    
    
    return
    

    
def finalInputRes():
    folder = r'Z:\Labs\Frank Lab\Carmel\Ephys\Hippo'
    loc_dataSummary = folder + '\\ephysDataSummary-ayo.csv'
    
    df = pd.read_csv(loc_dataSummary)
     
    drug = 'CAP'
    
    df_drug = df.loc[df['Cond 1'] == drug]   
    
    baseline_inputRes= []
    drug_inputRes =[]
    for ind in df_drug.index:
        trialFolder = df_drug['Folder'][ind]
        print(trialFolder)
   
        try:
            df_ns = pd.read_csv(trialFolder + '\\analysedData\\negativeStepResults_py.csv')
            baseline = df_ns['slope'][0]
            baseline_inputRes.append(baseline)
            drug_0 = df_ns['slope'][1]
            drug_inputRes.append(drug_0)
        
        except:
            print("No negative steps")
 
        
    data = np.column_stack((baseline_inputRes,drug_inputRes))
    data=data/1e3
    
    boxprops = dict(linewidth=2)
    medianprops = dict(linewidth=2.5,color='DarkCyan')
    flierprops = dict(markersize=10,linestyle='none')
 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for repeat in range(len(data)):
         plt.plot([1,2],[data[repeat,0],data[repeat,1]],linewidth=2,color='#bdccd4')
     
    bp=plt.boxplot(data,boxprops=boxprops,medianprops=medianprops,flierprops=flierprops, widths = 0.4)
    my_xticks = ['Baseline', 'Capcasisin']
    x=np.array([1,2])
    plt.xticks(x, my_xticks)
    
    for whisker in bp['whiskers']:
         whisker.set(linewidth=3)
    for cap in bp['caps']:
         cap.set(linewidth=3)  
         
    for i in range(len(data[1])):
          y = data[:,i]
          x = np.random.normal(1+i, 0.04, size=len(y))
          P.plot(x, y, 'r.', markersize=14, alpha=0.6)
     
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    plt.ylabel('Input Resistance (MOhms)')
    # plt.xlim((0.5,3.1))
    plt.tight_layout()
    fig.set_size_inches(4.5,4.5) 
    pf.saveFigure(fig,folder + r'\\resultsCapcasisin\\','inputResistance_1uM-cap_poired_230412')
      
    ttest_base_wash = stats.ttest_rel(data[:,0],data[:,1])
    ttest_wash_UV = stats.ttest_rel(data[:,1],data[:,2])
    ttest_UV_drug = stats.ttest_rel(data[:,2],data[:,3])

    ttest_results = np.column_stack((ttest_base_wash,ttest_wash_UV,ttest_UV_drug))
    
    stats_baseline = gf.get_stats(data[:,0])
    stats_washOn = gf.get_stats(data[:,1])
    stats_UV = gf.get_stats(data[:,2])
    stats_drug = gf.get_stats(data[:,3])

    
    
    df_stats = pd.DataFrame(data=ttest_results)
    df_stats.insert(0, "What", ["test stat","p-value"], True)
    df_stats.to_csv(folder + r'\\stats_threshold_trans_OCT-CAP_notIncubated_5uM_230519.csv')
    
   
    
    
    return
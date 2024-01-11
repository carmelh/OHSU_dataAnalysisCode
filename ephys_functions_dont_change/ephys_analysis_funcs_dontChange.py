# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 13:31:10 2023

@author: howeca
"""

import os
import numpy as np
import pyabf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
from scipy import stats, signal
import ephys_analysis_figure_funcs_dontChange as efig
import time


timestr = time.strftime("%y%m%d") # for saving figures and filenames with current date

# to ignore some annoying warnings
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


####### Universal Functions #######

def makeFoldersAnalysis(folder,expParam):
    ''' Generates either voltage clamp or current clamp folders 
    in the specified experiment folder '''
    
    subDirsNames=np.array(os.listdir(folder))
    
    for trial in range(len(subDirsNames)):
        cellStruct = folder + '\\' + subDirsNames[trial]
        print(cellStruct)
        if expParam == "cClamp": # if the expeirment is current clamp then generate these folders
            try:
                os.mkdir(os.path.join(cellStruct, "figures"))
                os.mkdir(os.path.join(cellStruct, "analysedData"))
                os.mkdir(os.path.join(cellStruct, "currentSteps"))
                os.mkdir(os.path.join(cellStruct, "currentSteps_negative"))
                os.mkdir(os.path.join(cellStruct, "currentHolds"))
                os.mkdir(os.path.join(cellStruct, "currentHolds_UV"))
                os.mkdir(os.path.join(cellStruct, "currentHolds_drug"))
                os.mkdir(os.path.join(cellStruct, "no_good"))
            except OSError as error:
                print(error) 
        elif expParam == "vClamp": # if the expeirment is voltage clamp then generate these folders
            try:
                os.mkdir(os.path.join(cellStruct, "figures"))
                os.mkdir(os.path.join(cellStruct, "analysedData"))
                os.mkdir(os.path.join(cellStruct, "vSteps"))
                os.mkdir(os.path.join(cellStruct, "holds"))
                os.mkdir(os.path.join(cellStruct, "holds_UV"))
                os.mkdir(os.path.join(cellStruct, "holds_drug"))
                os.mkdir(os.path.join(cellStruct, "no_good"))
            except OSError as error:
                print(error) 
                


def steps_Import(folder,filenames):  
    ''' # import steps files from .abf files '''
    # output variables or tuple: 
    # time vector
    # steps - 1D vector - applied steps(sweep x time)
    # voltage - 2D array recorded output voltage or current (sweep x time)
    
    abf = pyabf.ABF(folder + '\\' + filenames)
      
    time = abf.sweepX
    numSteps = len(abf.sweepList)
      
    voltage=np.zeros((numSteps,len(time)))
    steps=np.zeros((numSteps,len(time)))
    for sweep in range(numSteps):
        abf.setSweep(sweep)
        voltage[sweep] = abf.sweepY
        steps[sweep] = abf.sweepC
      
    return time,steps,voltage              



def absoluteStepValues(data_all):
    ''' get the steps values of the applied steps '''
    # Outputs:
    # stepsDelta = how much the current changed by each step in nA
    # stepsAbs = absolute value of current injected in pA
    
    steps = data_all[2]
    stepsDelta = np.around(np.max(steps[1]-steps[0]),decimals=2)
    # for negative steps
    if stepsDelta == 0:
        stepsDelta = np.around(np.min(steps[1]-steps[0]),decimals=2)
    stepsAbs = np.arange(0,stepsDelta*len(steps),stepsDelta)*1000 # steps in pA
    stepsAbs = stepsAbs[0:len(steps)]
    
    return stepsDelta, stepsAbs



def currentSteps_allTrials(folder,pos_neg):  
    '''# import current steps for every file (trial) within defined folder 
    Just repeats the stepsImport function for all files '''
    # Outputs:
    # filenames - of the .abf files
    # time vector
    # steps - 2D array applied steps(sweep x time)
    # voltaga_all - 3D array recorded output voltage(trial x sweep x time)
    
    if pos_neg == 'pos':
        stepsFolder =folder + '\currentSteps'
    elif pos_neg == 'neg':
        stepsFolder=folder + '\currentSteps_negative'
    
    filenames=np.array(os.listdir(stepsFolder))
    filenames.sort()

    # because unknown size at this point assign empty arrays
    time_all= []
    steps_all=[]
    voltage_all=[]
    for trial in range(len(filenames)): # runs through every file within the folder
        timeHolder,stepsHolder,voltageHolder = steps_Import(stepsFolder,filenames[trial])
        time_all.append(timeHolder)
        steps_all.append(stepsHolder)
        voltage_all.append(voltageHolder)
        
        
    # because time and steps should be the same for all repeats. compare and reduce array size
    # stack arrays, get the differentiation along the axis of stacking 
    # check if all of those differentiations are equal to zeros
    time_same = (np.diff(np.vstack(time_all).reshape(len(time_all),-1),axis=0)==0).all()
    if time_same == True:
        time=time_all[0]
    else:
        print("your time arrays are different... you should figure out why")
    
    
    steps_same = (np.diff(np.vstack(steps_all).reshape(len(steps_all),-1),axis=0)==0).any()
    if steps_same == True:
        steps=steps_all[0]
    else:
        print("your steps arrays are different... you should figure out why")
        
    return filenames, time, steps, np.asanyarray(voltage_all)



# this is what you use for data_all in the following 
def runAllSteps(folder): 
    ''' runs all steps through the folder
    Checks the files are there and if so runs processing and saves the ouput'''
    # Outputs:
    # output_posSteps - 3D array of positive current steps (time x steps x recorded signal)
    # output_negSteps - 3D array of negative current steps(time x steps x recorded signal)

    pos_neg = 'pos'
    try:
        output_posSteps = currentSteps_allTrials(folder,pos_neg)
        np.save(folder + r'\\analysedData\\' + 'processedABFfiles_{}Steps.npy'.format(pos_neg),output_posSteps)
    except:
        output_posSteps = []
        print("No positive steps")      
        
    pos_neg = 'neg'
    try:
        output_negSteps = currentSteps_allTrials(folder,pos_neg)
        np.save(folder + r'\\analysedData\\' + 'processedABFfiles_{}Steps.npy'.format(pos_neg),output_negSteps)
    except:
        output_negSteps = []
        print("No negative steps")

    return output_posSteps, output_negSteps 



def startEnd_steps(stimulated,pos_neg):
    ''' take the average hold current then find where current step starts
    average is taken from first 5000 time points. can shorten if access resistance
    step comes earlier. but probably shouldnt haha '''
    # Outputs:
    # step_start - index where the step starts (pos or neg or voltage clamp step)
    # step_end - index where the step ends 
    # step_value - the amplitude of the step
    
    if pos_neg == 'neg':
        step_start = next(x for x, val in enumerate(stimulated)
                                      if val < np.mean(stimulated[0:5000]))
        
        # negative current step amplitude 
        # in nA
        step_value = stimulated[step_start] - stimulated[0]

        # get where the negative step ends
        step_end = (next(x for x, val in enumerate(stimulated[step_start:len(stimulated)])
                                      if val > stimulated[step_start])) + step_start
    elif pos_neg == 'pos':
        step_start = next(x for x, val in enumerate(stimulated)
                                      if val > np.mean(stimulated[0:5000]))
        
        # postive current step amplitude 
        # in nA
        step_value = stimulated[step_start] - stimulated[0]
    
        # get where the  step ends
        step_end = (next(x for x, val in enumerate(stimulated[step_start:len(stimulated)])
                                      if val < stimulated[step_start])) + step_start
    
    elif pos_neg=='voltage':
        step_start = next(x for x, val in enumerate(stimulated[0])
                                      if val < -80)
        
        step_value = (stimulated[0,step_start] - stimulated[1,step_start])*-1

        # get where the negative step ends
        step_end = (next(x for x, val in enumerate(stimulated[0,step_start:len(stimulated[0])])
                                      if val > stimulated[0,step_start])) + step_start
        
    return step_start, step_end, step_value



def spikeIdx(trace,fs):
    ''' Current clamp function to find the index of spikes (action potentials) '''
    # Outputs:
    # nspikes - how many spikes there were in the recorded output for one current step
    # peaks - the indices where they occured 
    
    set_crossgi= [i for i,v in enumerate(trace) if v > -10]
    if len(set_crossgi) > 1:  #This to make sure there is a spike otherwise the code below gives problems. There is an empty else statement below.
        # height is thresold voltage in mV
        # distance is to remove noisy peaks. is about 0.5 ms between APs
        peaks, _ = signal.find_peaks(trace, height=-10, distance=np.round(0.5e-3/(1/fs))) 
        nspikes= len(peaks) # Number of spikes
       # print("you have {} spikes".format(nspikes))
    else:
        peaks=[]
        nspikes= len(peaks)

    return nspikes, peaks



def threshold(data):   
    ''' find index where num spikes first exceeds zero '''
    
    threshold_idx = np.zeros((len(data)))
    for ii in range(len(data)):
        threshold_idx[ii] = next(x for x, val in enumerate(data[ii])
                                  if val > 0)
        
    return threshold_idx



def numSpikes(output_posSteps,folder,expDrug):
    ''' count number of spikes for each current step '''
    # Outputs:
    # stepValues - the step values in pA
    # num_spikes - how many spikes for each current step
    # threshold_idx - value of the current step where the first AP occurs
    # Saves .csv in analysedData folder with num_spikes and threshold_idx for every file in the input data
    
    voltage = output_posSteps[3]
    current = output_posSteps[2]
    
    stepInfo = startEnd_steps(current[1],'pos')
    stepValues = np.arange(0,len(voltage[0]),stepInfo[2]*100)
    
    # run through every file in data_all and extracts spike information
    num = np.zeros((len(voltage[0]))) 
    num_spikes = np.zeros((len(voltage),len(voltage[0])))
    pks =[]
    for trial in range(len(voltage)):
        volt_trial = voltage[trial,...]
        for sweep in range(len(voltage[0])):
            num[sweep], peaks = spikeIdx(volt_trial[sweep,stepInfo[0]:stepInfo[1]],fs=1/output_posSteps[1][1]) 
            pks.append(peaks)
        num_spikes[trial] = num
        
    stepValues=absoluteStepValues(output_posSteps) # find the step value for each current injection
    threshold_idx = threshold(num_spikes)
    threshold_idx = (threshold_idx * stepValues[0])*1000 # in pA
    
    # generates and saves pandas dataframe as .csv
    df_ns = pd.DataFrame(data=num_spikes,columns=stepValues[1],index=output_posSteps[0])
    df_ns['Threshold'] = threshold_idx
    df_ns.to_csv(folder + r'\\analysedData\\' + 'numSpikesThreshold_py.csv')
    
    colors = plt.cm.winter(np.linspace(0,1,len(num_spikes)))

    fig = plt.gcf()      
    ax = plt.subplot(111)  
    # set y axis to be integer values only
    for axis in [ax.yaxis]:
        axis.set_major_locator(ticker.MaxNLocator(integer=True))
    for sweepNumber in range(len(num_spikes)):
        plt.plot(stepValues[1],np.transpose(num_spikes[sweepNumber]), color=colors[sweepNumber], linewidth=2, alpha=.8)
    plt.plot(stepValues[1],num_spikes[0], linewidth=2,color='k',label='Baseline')
    plt.plot(stepValues[1],num_spikes[len(num_spikes)-1], linewidth=2,color='red',label=expDrug)
   
    efig.lrBorders(ax)
    fig.set_size_inches(4.5,4.5)
    plt.xlabel('Current (pA)')
    plt.ylabel('Num. Spikes')
    ax.legend(frameon=False)
    efig.saveFigure(fig,folder + r'\\figures\\','currentInjection_vs_numSpikes_all_drug_{}_{}'.format(expDrug,timestr))
    
    print("Analysed number of spikes...")
    return stepValues, num_spikes, threshold_idx



def access_resistance(output_posSteps, folder, expDrug):
    ''' calculates the access resistance from the negative current step
    that occurs at the begining of each recording '''
    
    voltage = output_posSteps[3]
    current = output_posSteps[2]

    stepInfo = startEnd_steps(current[0],'neg')
    
    # voltage value during negative current step
    recorded_amplitude = np.mean(voltage[:,:,stepInfo[0]:stepInfo[1]],2)
    
    # R = V (mV) / I (mA) 
    access_resistance = (recorded_amplitude) / (stepInfo[2]/1e6)
    access_resistance = access_resistance/ 1e8 # in MOhms
    
    # set up pandas dataframe to save
    stepValues=absoluteStepValues(output_posSteps)
    df_ar = pd.DataFrame(data=access_resistance,columns=stepValues[1],index=output_posSteps[0])
    
    fileAverage = np.mean(access_resistance,1)
    fileDeviation = np.std(access_resistance,1) # shouldnt vary too much
    
    df_ar['average'] = fileAverage
    df_ar['deviation'] = fileDeviation
    
    df_ar.to_csv(folder + r'\\analysedData\\' + 'access_resistance_inMOhms_py.csv')
    
    # plot figure of access resistance for each current step
    fig= plt.gcf()      
    ax = plt.subplot(111)  
    plt.plot(stepValues[1],np.transpose(access_resistance))
    efig.lrBorders(ax)
    fig.set_size_inches(6,4.5)
    plt.xlabel('Current (pA)')
    plt.ylabel('Access Resistance MOhms')
    ax.legend(labels=list(output_posSteps[0]) ,frameon=False, fontsize='xx-small',loc='upper left', bbox_to_anchor=(1, 1.05))
    efig.saveFigurePNG(fig,folder + r'\\figures\\','currentInjection_vs_accessResistance_all_drug_{}_{}'.format(expDrug,timestr))
    
    print("Analysed access resistance...")
    return


def rmp(output_posSteps, folder, expDrug):
    ''' calculates resting membrane potential from first 12000 frames '''
    # Outputs:
    # .csv of resting membrane potential
    
    voltage = output_posSteps[3]
    
    stepValues=absoluteStepValues(output_posSteps)
    
    rmp=np.mean(voltage[:,:,0:12000],2) # takes baseline mV
    
    # set up pandas dataframe to save
    df_rmp = pd.DataFrame(data=rmp,columns=stepValues[1],index=output_posSteps[0])
    
    fileAverage = np.mean(rmp,1)
    fileDeviation = np.std(rmp,1) # shouldnt vary too much
    
    df_rmp['average'] = fileAverage
    df_rmp['deviation'] = fileDeviation
    
    df_rmp.to_csv(folder + r'\\analysedData\\' + 'rmp_py.csv')
    
    
    # plot figure of access resistance for each current step
    fig= plt.gcf()      
    ax = plt.subplot(111)  
    plt.plot(stepValues[1],np.transpose(rmp))
    efig.lrBorders(ax)
    fig.set_size_inches(6,4.5)
    plt.xlabel('Current (pA)')
    plt.ylabel('RMP mV')
    ax.legend(labels=list(output_posSteps[0]) ,frameon=False, fontsize='xx-small',loc='upper left', bbox_to_anchor=(1, 1.05))
    efig.saveFigurePNG(fig,folder + r'\\figures\\','currentInjection_vs_rmp_all_drug_{}_{}'.format(expDrug,timestr))
    
    
    print("Analysed resting membrane potential...")
    return


def inputResistance(output_negSteps, folder, expDrug):
    ''' calculates the input resistance from the negative current steps'''
    voltage = output_negSteps[3]
    current = output_negSteps[2]
    
    stepInfo = startEnd_steps(current[1],'neg')
    
    # take the last half of the negative current step "steady state"
    ss_Start = int((stepInfo[1]-stepInfo[0])/2 + stepInfo[0])
    
    recorded_amplitude = np.mean(voltage[:,:,ss_Start:stepInfo[1]],2)
    
    # set up dataframe to save
    stepValues=absoluteStepValues(output_negSteps)
    df = pd.DataFrame(data=recorded_amplitude,columns=stepValues[1],index=output_negSteps[0])
    
    df.to_csv(folder + r'\\analysedData\\' + 'negativeStepValues_in-mV_py.csv')
    

    linRegress_results=[] # slope in ohms, intercept in mV
    for trial in range(len(recorded_amplitude)): 
         hdr = stats.linregress(stepValues[1]*1e-6,recorded_amplitude[trial]) 
         linRegress_results.append(hdr)

    
    df_ns = pd.DataFrame(data=linRegress_results,index=output_negSteps[0])
    df_ns.to_csv(folder + r'\\analysedData\\' + 'negativeStepResults_py.csv')
    
    
    ### figure ####
    colors = plt.cm.winter(np.linspace(0,1,len(recorded_amplitude)))

    fig = plt.gcf()      
    ax = plt.subplot(111)  
    # set y axis to be integer values only
    for axis in [ax.yaxis]:
        axis.set_major_locator(ticker.MaxNLocator(integer=True))
    for sweepNumber in range(len(recorded_amplitude)):
        plt.plot(stepValues[1], linRegress_results[sweepNumber].intercept + (linRegress_results[sweepNumber].slope/1e6)*stepValues[1], '--', color=colors[sweepNumber], linewidth=2, alpha=.8)
    for sweepNumber in range(len(recorded_amplitude)):
        plt.plot(stepValues[1],recorded_amplitude[sweepNumber],'o',color=colors[sweepNumber])
    
    ax.legend(labels=list(output_negSteps[0]) ,frameon=False, fontsize='xx-small',loc='upper left', bbox_to_anchor=(1, 1.05))
    efig.lrBorders(ax)
    fig.set_size_inches(4.5,4.5)
    plt.xlabel('Current (pA)')
    plt.ylabel('Voltage (mV)')
    efig.saveFigurePNG(fig,folder + r'\\figures\\','inputResistance_drug_{}_{}'.format(expDrug,timestr))
    
    
    print("Analysed negative current steps...")
    return


def voltageSteps_allTrials(stepsFolder,filenames):
    ''' import voltage steps from all files in folder '''
    # filenames - of the .abf files
    # time vector
    # steps - 2D array applied steps(sweep x time)
    # current_all - 3D array recorded output (trial x sweep x time)
    
    # because unknown size at this point assign empty arrays
    time_all= []
    steps_all=[]
    current_all=[]
    for trial in range(len(filenames)): # runs through every file within the folder
        timeHolder,stepsHolder,currentHolder = steps_Import(stepsFolder,filenames[trial])
        time_all.append(timeHolder)
        steps_all.append(stepsHolder)
        current_all.append(currentHolder)
        
        
    # because time and steps should be the same for all repeats. compare and reduce array size
    # stack arrays, get the differentiation along the axis of stacking 
    # check if all of those differentiations are equal to zeros
    time_same = (np.diff(np.vstack(time_all).reshape(len(time_all),-1),axis=0)==0).all()
    if time_same == True:
        time=time_all[0]
       # print("Time reduced - good")
    else:
        print("your time arrays are different... you should figure out why")
    
    
    steps_same = (np.diff(np.vstack(steps_all).reshape(len(steps_all),-1),axis=0)==0).any()
    #steps_same = (np.diff(np.vstack(steps_all).reshape(len(steps_all),-1),axis=0)==0).all()
    if steps_same == True:
        steps=steps_all[0]
    #    print("Steps reduced - good")
    else:
        print("your steps arrays are different... you should figure out why")
        
    return filenames, time, steps, np.asanyarray(current_all)


def runAllSteps_voltage(folder):
    ''' runs all steps through the folder
    Checks the files are there and if so runs processing and saves the ouput'''
    # Outputs:
    # output_vSteps - 3D array of voltage steps (time x steps x recorded signal)
    
    stepsFolder =folder + r'\vSteps'
    
    filenames=np.array(os.listdir(stepsFolder))
    filenames.sort()
    
    output_vSteps = voltageSteps_allTrials(stepsFolder,filenames)
    np.save(folder + r'\\analysedData\\' + 'processedABFfiles_voltageSteps.npy',output_vSteps)
 
    return output_vSteps


def voltage_IV(output_vSteps, folder):
    ''' calculates steady state current for applied voltage steps '''
    
    current = output_vSteps[3]
    voltage = output_vSteps[2]
    
    stepInfo = startEnd_steps(voltage,'voltage')
    stepValues = np.arange(-100,len(voltage[:,0])*10,stepInfo[2])
    
    ss_Start = int((stepInfo[1]-stepInfo[0])/2 + stepInfo[0]) #takes the last half of the voltage step as "steady state"
    steadyStateCurrent = np.mean(current[:,:,ss_Start:stepInfo[1]],2)    
    
    df = pd.DataFrame(data=steadyStateCurrent,columns=stepValues,index=output_vSteps[0])
    df.to_csv(folder + r'\\analysedData\\' + 'steadyStateCurrent_perVoltageStep_in-pA_py.csv')    
    return
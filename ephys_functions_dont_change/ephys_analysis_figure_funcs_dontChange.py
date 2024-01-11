# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 12:12:23 2023

@author: howeca
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import matplotlib.colors as colors
from datetime import date
import matplotlib.cm as cm
import scipy.ndimage as ndimage
import time


timestr = time.strftime("%y%m%d") # for saving figures and filenames with current date


# sets font text
font = {'family': 'sans',
        'weight': 'normal',
        'size': 18,}

plt.rc('font',**font)    
    

def saveFigure(fig,path,keyword):
    plt.savefig(path + r'\\{}.png'.format(keyword), format='png', dpi=600, bbox_inches='tight')
    plt.savefig(path + r'\\{}.eps'.format(keyword), format='eps', dpi=1900, bbox_inches='tight')
    plt.close(fig)   
    return


def saveFigurePNG(fig,path,keyword):
    plt.savefig(path + r'\\{}.png'.format(keyword), format='png', dpi=600, bbox_inches='tight')
    plt.close(fig)   
    return


def noBorders(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(False)
    ax.spines['left'].set_linewidth(False)
    ax.axes.get_yaxis().set_ticks([])
    ax.axes.get_xaxis().set_ticks([]) 
    return


def lrBorders(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    return




def numSpikes_figure(folder, expDrug):
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
    for sweepNumber in range(len(num_spikes)):
        plt.plot(stepValues,np.transpose(num_spikes[sweepNumber]), color=colors[sweepNumber], linewidth=2, alpha=.8)
    plt.plot(stepValues,num_spikes[0], linewidth=2,color='k',label='Baseline')
    plt.plot(stepValues,num_spikes[len(num_spikes)-1], linewidth=2,color='red',label=expDrug)
   
    lrBorders(ax)
    fig.set_size_inches(4.5,4.5)
    plt.xlabel('Current (pA)')
    plt.ylabel('Num. Spikes')
    ax.legend(frameon=False)
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 2)) # change tick size from 1 to 2
    plt.xticks(rotation=30)
    saveFigure(fig,folder + r'\\figures\\','currentInjection_vs_numSpikes_all_{}'.format(timestr))
    
    print('saved num. spikes figure')
    
    return



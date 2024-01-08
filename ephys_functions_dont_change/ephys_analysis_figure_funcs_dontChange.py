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



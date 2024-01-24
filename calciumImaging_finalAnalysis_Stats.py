# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 13:33:05 2024

Final statistical analysis of the calcium imaging HEK cell OCT-CAP data

@author: Carmel Howe
"""


import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import seaborn as sns
import time


timestr = time.strftime("%y%m%d") # for saving figures and filenames with current date


folder=r'Z:\Labs\Frank Lab\Carmel\Imaging\Confocal_olympusFV1200\Live\OCTCAP_HEK_jRGECO\HEK_jRGECOimaging_finalAnalysis'
loc_dataSummary =folder + '\avergaeEachTrial_btw70and90seconds_231130.csv'
df = pd.read_csv(loc_dataSummary)


# boxplot for different conditions and the difference between baseline and
# max condition
fig = plt.figure()
sns.boxplot(x='What',y='average',data=df, width=0.75, showfliers = False,
        flierprops={"markersize": 8},
        boxprops={"facecolor": (.0, .0, .0, .0)},
        medianprops={"color": "DarkCyan"},)
sns.stripplot(x='What', y='average', data=df, alpha=0.5, color='r',size=6) # add individual points
plt.xticks(rotation=30)       
fig.set_size_inches(5,5) 

plt.savefig(folder + r'\\figures\\','finalAnalysis_boxplot_{}'.format(timestr), format='png', dpi=600, bbox_inches='tight')
plt.savefig(folder + r'\\figures\\','finalAnalysis_boxplot_{}'.format(timestr), format='eps', dpi=1900, bbox_inches='tight')
plt.close(fig)   


# summary stats
summaryStats = df.groupby(['What']).agg(['mean','median','std']) 
summaryStats.to_csv(folder + '\\summaryStats_{}.csv'.format(timestr))




# ttest between dmso control and different conditions
# cannot find a more sophisticated way to do this
dmso = df[df['What'] == 'DMSO']
teth = df[df['What'] == 'Teth']
trpv1_neg = df[df['What'] == 'Trpv1_neg']
snap_mutant = df[df['What'] == 'SNAPmut']


ttest_dmsoTeth = stats.ttest_ind(dmso['average'],teth['average'])
ttest_dmsoNoTRPV1 = stats.ttest_ind(dmso['average'],trpv1_neg['average'])
ttest_dmsoMut = stats.ttest_ind(dmso['average'],snap_mutant['average'])


# send ttests to dictionry then convert to df and save
ttestDict = {'dmsoTeth': ttest_dmsoTeth, 'dmsoNoTRPV1': ttest_dmsoNoTRPV1, 'dmsoSNAPmut': ttest_dmsoMut}

ttest_df = pd.DataFrame(ttestDict, index=['statistic', 'p-value'])
ttest_df.to_csv(folder + '\\indTtestResults_baselineTethMutantNoTRPV1_{}.csv'.format(timestr))


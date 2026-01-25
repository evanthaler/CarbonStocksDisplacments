import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


disp_df = pd.read_csv('/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/disp2use/TL27TL47Displacments_PFTopo.csv')


tl27_pfdisp = disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "1"),'dhorizontal'].dropna()
tl27_nonpfdisp = disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "0"),'dhorizontal'].dropna()
tl47_pfdisp = disp_df.loc[(disp_df['Site'] == 'TL47') & (disp_df['pf'] == "1"),'dhorizontal'].dropna()
tl47_nonpfdisp = disp_df.loc[(disp_df['Site'] == 'TL47') & (disp_df['pf'] == "0"),'dhorizontal'].dropna()

tl27_pfslope = disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "1"),'slope'].dropna()
tl27_nonpfslope = disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "0"),'slope'].dropna()
tl47_pfslope = disp_df.loc[(disp_df['Site'] == 'TL47') & (disp_df['pf'] == "1"),'slope'].dropna()
tl47_nonpfslope = disp_df.loc[(disp_df['Site'] == 'TL47') & (disp_df['pf'] == "0"),'slope'].dropna()



colors = ['blue','lightblue','black','darkgray']
fig,ax = plt.subplots(figsize=(5,5))
bplot = ax.boxplot(
    [tl27_pfslope,tl27_nonpfslope,tl47_pfslope,tl47_nonpfslope],
    tick_labels =['TL27 PF','TL27 NonPF','TL47 PF','TL47 NonPF'],
    patch_artist=True,
    showfliers=False )
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
ax.tick_params(axis='x', labelrotation=45)
ax.set_ylabel('Slope (degrees)')
fig.show()

colors = ['blue','lightblue','black','darkgray']
fig,ax = plt.subplots(figsize=(5,5))
bplot = ax.boxplot(
    [tl27_pfdisp,tl27_nonpfdisp,tl47_pfdisp,tl47_nonpfdisp],
    tick_labels =['TL27 PF','TL27 NonPF','TL47 PF','TL47 NonPF'],
    patch_artist=True,
    showfliers=False )
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
ax.tick_params(axis='x', labelrotation=45)
ax.set_ylabel('Horizontal displacement (m yr$^{-1}$)')
fig.show()





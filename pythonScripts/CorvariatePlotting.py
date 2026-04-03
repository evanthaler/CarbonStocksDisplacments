import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
plt.rcParams.update({'font.size': 14})
df_sites = pd.read_csv('/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/wDisplacement/TL47TL27StocksDisplacement_covariates_CNratio.csv')
figoutdir = '/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/figs'



x,y = df_sites.mean_disp_rate*0.5,df_sites.C_N_ratio
cn_linstats = linregress(x,y)
fig,ax = plt.subplots(figsize=(6,6))
ax.plot(x,y,'.k',markersize=12)
ax.plot(x,x*cn_linstats.slope+cn_linstats.intercept,'-k')
ax.set_ylabel('C:N')
ax.set_xlabel('Soil flux (m$^{2}$ yr$^{-1}$)')
plt.text(0.05*0.5, 13,f'$r^{2}$={cn_linstats.rvalue**2:.2f}')
plt.text(0.05*0.5, 12.5,f'$p$={cn_linstats.pvalue:.3f}')
plt.savefig(f'{figoutdir}/CNvsFluxes.jpg',dpi=300)
plt.show()
print(cn_linstats.rvalue**2,cn_linstats.pvalue)





# ###################################################################

#########
#Plotting Displacment and covars for each site
#########
cols_to_plot = ['slope','curv','drainagearea']
cols = ['Slope (degrees)',
            'Curvature (m$^{-1}$)',
            'Drainage area (m$^{2}$)']

fig, axes = plt.subplots(1,3, figsize=(10, 5), sharey=True)
axes = axes.flatten()
for ax, col,colname in zip(axes, cols_to_plot,cols):
    x,y = df_sites[col],df_sites.mean_disp_rate*0.5
    linstats = linregress(x,y)
    
    yline = linstats.slope * x + linstats.intercept

    ax.plot(x,y,  'ok')
    #ax.text(0.25*max(x),0.025,f'$r^{2}$={linstats.rvalue**2:.2f}')
    ax.text( 0.05, 0.9,  # position (top-left)
        f'$R^2$ = {linstats.rvalue**2:.2f}\np = {linstats.pvalue:.3f}',
        transform=ax.transAxes,
        verticalalignment='top')
    if linstats.pvalue < 0.05:
        ax.plot(x, yline, '-k', lw=3)
    else:
        ax.plot(x, yline, '-',color='gray', lw=1, ms=10)

    ax.set_xlabel(colname)
for ax in axes[[0]]:
    ax.set_ylabel('Soil flux (m$^{2}$ yr$^{-1}$)')


plt.tight_layout()
plt.savefig(f'{figoutdir}/DisplacementCovariates.jpg', dpi=300)
plt.show()

######################################


fig, axes = plt.subplots(1,3, figsize=(10, 5), sharey=True)
axes = axes.flatten()
for ax, col,colname in zip(axes, cols_to_plot,cols):
    x,y = df_sites[col],df_sites.Total_C_stock_kg_m2_0to50cm
    linstats = linregress(x,y)
    
    yline = linstats.slope * x + linstats.intercept

    ax.plot(x,y,  'ok')
    #ax.text(0.25*max(x),0.025,f'$r^{2}$={linstats.rvalue**2:.2f}')
    ax.text( 0.05, 0.9,  # position (top-left)
        f'$R^2$ = {linstats.rvalue**2:.2f}\np = {linstats.pvalue:.3f}',
        transform=ax.transAxes,
        verticalalignment='top')
    if linstats.pvalue < 0.05:
        ax.plot(x, yline, '-k', lw=3)
    else:
        ax.plot(x, yline, '-',color='gray', lw=1, ms=10)

    ax.set_xlabel(colname)
for ax in axes[[0]]:
    ax.set_ylabel('Soil organic carbon stock (kg m$^{-2}$)')


plt.tight_layout()
plt.savefig(f'{figoutdir}/CStocksCovariates.jpg', dpi=300)
plt.show()
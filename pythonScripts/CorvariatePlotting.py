import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
plt.rcParams.update({'font.size': 14})
df_sites = pd.read_csv('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/FinalCleanedFiles/wDisplacement/TL47TL27StocksDisplacement_covariates.csv')

#################################
#Displacement fit all data ######
#################################
def linear_exp(x, x0, y0, m, k):
    return np.where(
        x <= x0,
        m * (x - x0) + y0,
        y0 * np.exp(-k * (x - x0))
    )
def fitDispAllData(df,outfig,xcol='mean_disp_rate',ycol='Total_C_stock_kg_m2_0to50cm'):
    x = df[xcol]
    y = df[ycol]
    x0_init = x[np.argmax(y)]# pick the y peak location for x
    y0_init = y.max()

    p0 = [
        x0_init,      # x0 (peak location)
        y0_init,      # y0 (peak value)
        200,          # m (slope; adjust if needed)
        50            # k (decay rate; adjust if needed)
    ]

    popt, pcov = curve_fit(
        linear_exp, x, y, p0=p0, maxfev=10000
    )
    x0, y0, m, k = popt
    #Get peak standard error
    x0_se = np.sqrt(pcov[0, 0])          # standard error
    x0_ci = 1.96 * x0_se   

    
    x_fit = np.linspace(x.min(), x.max(), 400)
    y_fit = linear_exp(x_fit, *popt)
    
    plt.figure()
    plt.plot(df.mean_disp_rate,df.Total_C_stock_kg_m2_0to50cm,'ok')
    #plt.scatter(df.mean_disp_rate,df.Total_C_stock_kg_m2_0to50cm,c=df.curv,s=500)
    plt.colorbar(label='Slope (degrees)')
    plt.plot(x_fit, y_fit, color='k', lw=2)
    plt.axvspan(
        x0 - x0_ci, x0 + x0_ci,
        color="gray", alpha=0.2,
        label="95% CI (peak)"
    )
    plt.ylabel('Soil organic carbon stock (kg m$^{-2}$)')
    plt.xlabel('Horizontal displacement rate (m yr$^{-1}$)')
    #plt.savefig('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/figs/SOCStock_DisplacementCombinedSites.jpg',dpi=300)
    plt.show()


def normalizeDisplacement(df,site,dispcol = 'mean_disp_rate',sitecol='Site'):
    df_site = df[df[sitecol]==site]
    disp_norm = df_site[dispcol]/np.max(df_site[dispcol])
    return disp_norm
    
tl27dispnorm = normalizeDisplacement(df_sites,'TL27')
tl47dispnorm = normalizeDisplacement(df_sites,'TL47')

##################################################
#Just some renaming to make labeling prettier later
##################################################
col_dict = {'slope':'Slope (degrees)',
            'curv':'Curvature (m$^{-1}$)',
            'drainagearea':'Drainage area (m$^{2}$)',
            'ndvi':'NDVI'}
df_sites = df_sites.rename(columns = col_dict)

##########################
df_sites['normalized_disp'] = np.zeros(len(df_sites['mean_disp_rate']))
df_sites.loc[df_sites['Site'] == 'TL27', 'normalized_disp'] = tl27dispnorm
df_sites.loc[df_sites['Site'] == 'TL47', 'normalized_disp'] = tl47dispnorm
#######################
tl27_df = df_sites[df_sites.Site=='TL27']
tl47_df = df_sites[df_sites.Site=='TL47']

#########
#Plotting
#########
cols = ['Slope (degrees)',
            'Curvature (m$^{-1}$)',
            'Drainage area (m$^{2}$)',
            'NDVI']
fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
axes = axes.flatten()
for ax, col in zip(axes, cols):
    linstats_47 = linregress(tl47_df.normalized_disp, tl47_df[col])
    linstats_27 = linregress(tl27_df.normalized_disp, tl27_df[col])
    yline_47 = linstats_47.slope * tl47_df.normalized_disp + linstats_47.intercept
    yline_27 = linstats_27.slope * tl27_df.normalized_disp + linstats_27.intercept
    ax.plot(tl47_df.normalized_disp, tl47_df[col], 'ok', label='TL47')
    ax.plot(tl27_df.normalized_disp, tl27_df[col], 'ob', label='TL27')
    if linstats_47.pvalue < 0.05:
        ax.plot(tl47_df.normalized_disp, yline_47, '-k', lw=2)
    else:
        ax.plot(tl47_df.normalized_disp, yline_47, 'x-k', lw=2, ms=10)
    if linstats_27.pvalue < 0.05:
        ax.plot(tl27_df.normalized_disp, yline_27, '-b', lw=2)
    else:
        ax.plot(tl27_df.normalized_disp, yline_27, 'x-b', lw=2, ms=10)
    ax.set_ylabel(col)
    print(f"{col} | TL47 p={linstats_47.pvalue:.3g}, TL27 p={linstats_27.pvalue:.3g}")
for ax in axes[2:]:
    ax.set_xlabel('Normalized horizontal displacement')

axes[1].legend(loc='upper left')
plt.tight_layout()
plt.savefig('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/figs/DisplacementCovariates.jpg', dpi=300)
plt.show()

cols = ['Slope (degrees)',
            'Curvature (m$^{-1}$)',
            'Drainage area (m$^{2}$)',
            'NDVI']
fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
axes = axes.flatten()
for ax, col in zip(axes, cols):
    linstats_47 = linregress(tl47_df.Total_C_stock_kg_m2_0to50cm, tl47_df[col])
    linstats_27 = linregress(tl27_df.Total_C_stock_kg_m2_0to50cm, tl27_df[col])
    yline_47 = linstats_47.slope * tl47_df.Total_C_stock_kg_m2_0to50cm + linstats_47.intercept
    yline_27 = linstats_27.slope * tl27_df.Total_C_stock_kg_m2_0to50cm + linstats_27.intercept
    ax.plot(tl47_df.Total_C_stock_kg_m2_0to50cm, tl47_df[col], 'ok', label='TL47')
    ax.plot(tl27_df.Total_C_stock_kg_m2_0to50cm, tl27_df[col], 'ob', label='TL27')
    if linstats_47.pvalue < 0.05:
        ax.plot(tl47_df.Total_C_stock_kg_m2_0to50cm, yline_47, '-k', lw=2)
    else:
        ax.plot(tl47_df.Total_C_stock_kg_m2_0to50cm, yline_47, 'x-k', lw=2, ms=10)
    if linstats_27.pvalue < 0.05:
        ax.plot(tl27_df.Total_C_stock_kg_m2_0to50cm, yline_27, '-b', lw=2)
    else:
        ax.plot(tl27_df.Total_C_stock_kg_m2_0to50cm, yline_27, 'x-b', lw=2, ms=10)
    ax.set_ylabel(col)
    print(f"{col} | TL47 p={linstats_47.pvalue:.3g}, TL27 p={linstats_27.pvalue:.3g}")
for ax in axes[2:]:
    ax.set_xlabel('Soil organic carbon stocks (kg m$^{-2}$)')

axes[1].legend(loc='upper right')
plt.tight_layout()
plt.savefig('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/figs/CstocksCovariates.jpg', dpi=300)
plt.show()
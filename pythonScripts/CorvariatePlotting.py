import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
plt.rcParams.update({'font.size': 14})
df_sites = pd.read_csv('/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/wDisplacement/TL47TL27StocksDisplacement.csv')
figoutdir = '/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/figs'
#################################
#Displacement fit all data ######
#################################
def linear_exp(x, x0, y0, m, k):
    return np.where(
        x <= x0,
        m * (x - x0) + y0,
        y0 * np.exp(-k * (x - x0))
    )

def soil_prod_hump(x, a, b):
    return a * x * np.exp(-b * x)


def fitDispAllData(df,outfig,xcol='mean_disp_rate',ycol='Total_C_stock_kg_m2_0to50cm',ycollabel='Soil organic carbon stock (kg m$^{-2}$)',plotline=True):
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
    #plt.plot(df.mean_disp_rate,df[ycol],'ok')

    # plt.plot(
    # df_sites.loc[(df_sites['Site'] == 'TL47') & (df_sites['PF'] == 1), 'mean_disp_rate'],
    # df_sites.loc[(df_sites['Site'] == 'TL47') & (df_sites['PF'] == 1), ycol],
    # '*k')
    # plt.plot(
    # df_sites.loc[(df_sites['Site'] == 'TL47') & (df_sites['PF'] == 0), 'mean_disp_rate'],
    # df_sites.loc[(df_sites['Site'] == 'TL47') & (df_sites['PF'] == 0), ycol],
    # 'ok',
    # label='TL47')

    # plt.plot(
    # df_sites.loc[(df_sites['Site'] == 'TL27') & (df_sites['PF'] == 1), 'mean_disp_rate'],
    # df_sites.loc[(df_sites['Site'] == 'TL27') & (df_sites['PF'] == 1), ycol],
    # '*b')
    # plt.plot(
    # df_sites.loc[(df_sites['Site'] == 'TL27') & (df_sites['PF'] == 0), 'mean_disp_rate'],
    # df_sites.loc[(df_sites['Site'] == 'TL27') & (df_sites['PF'] == 0), ycol],
    # 'ob',
    # label='TL27')


    plt.plot(df_sites.loc[df_sites['Site'] == 'TL47','mean_disp_rate'],
    df_sites.loc[df_sites['Site'] == 'TL47',ycol],'ok',label='Teller 47')

    plt.plot(df_sites.loc[df_sites['Site'] == 'TL27','mean_disp_rate'],
    df_sites.loc[df_sites['Site'] == 'TL27',ycol],'ob',label='Teller 27')



    if plotline:
        plt.plot(x_fit, y_fit, color='k', lw=2)
        plt.axvspan(
            x0 - x0_ci, x0 + x0_ci,
            color="gray", alpha=0.2
        )
    plt.ylabel(ycollabel)
    plt.xlabel('Horizontal displacement rate (m yr$^{-1}$)')
    plt.tight_layout()
    plt.legend()
    plt.savefig(outfig,dpi=300)
    plt.show()


fitDispAllData(df_sites,f'{figoutdir}/SOCStockDisplacementCombinedSites_nocurvefit.jpg',plotline=False)

###################################################################
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
#############################
tl27_df = df_sites[df_sites.Site=='TL27']
tl47_df = df_sites[df_sites.Site=='TL47']

#########
#Plotting Displacment and covars for each site
#########
cols = ['Slope (degrees)',
            'Curvature (m$^{-1}$)',
            'Drainage area (m$^{2}$)']

fig, axes = plt.subplots(1,3, figsize=(10, 5), sharey=True)
axes = axes.flatten()
for ax, col in zip(axes, cols):
    linstats_47 = linregress(tl47_df[col],tl47_df.mean_disp_rate)
    linstats_27 = linregress(tl27_df[col],tl27_df.mean_disp_rate)
    yline_47 = linstats_47.slope * tl47_df[col] + linstats_47.intercept
    yline_27 = linstats_27.slope * tl27_df[col] + linstats_27.intercept
    ax.plot(tl47_df[col],tl47_df.mean_disp_rate,  'ok', label='TL47')
    ax.plot(tl27_df[col],tl27_df.mean_disp_rate,  'ob', label='TL27')
    if linstats_47.pvalue < 0.05:
        ax.plot(tl47_df[col], yline_47, '-k', lw=3)
    else:
        ax.plot(tl47_df[col], yline_47, '-',color='gray', lw=1, ms=10)
    if linstats_27.pvalue < 0.05:
        ax.plot(tl27_df[col], yline_27, '-b', lw=3)
    else:
        ax.plot(tl27_df[col], yline_27, '-',color='gray', lw=1, ms=10)
    ax.set_xlabel(col)
    print(f"{col} | TL47 p={linstats_47.pvalue:.3g}, TL27 p={linstats_27.pvalue:.3g}")
for ax in axes[[0]]:
    ax.set_ylabel('Horizontal displacement (m yr$^{-1}$)')

axes[0].legend(loc='best')
plt.tight_layout()
plt.savefig(f'{figoutdir}/DisplacementCovariates.jpg', dpi=300)
plt.show()

#########
#Plotting carbon and covars for each site
#########
fig, axes = plt.subplots(1,3, figsize=(10, 5), sharey=True)
axes = axes.flatten()
for ax, col in zip(axes, cols):
    linstats_47 = linregress(tl47_df[col],tl47_df.Total_C_stock_kg_m2_0to50cm )
    linstats_27 = linregress(tl27_df[col],tl27_df.Total_C_stock_kg_m2_0to50cm )
    yline_47 = linstats_47.slope * tl47_df[col] + linstats_47.intercept
    yline_27 = linstats_27.slope * tl27_df[col] + linstats_27.intercept
    ax.plot(tl47_df[col],tl47_df.Total_C_stock_kg_m2_0to50cm , 'ok', label='TL47')
    ax.plot(tl27_df[col],tl27_df.Total_C_stock_kg_m2_0to50cm, 'ob', label='TL27')
    if linstats_47.pvalue < 0.05:
        ax.plot(tl47_df[col], yline_47, '-k', lw=3)
    else:
        ax.plot(tl47_df[col], yline_47, '-',color='gray', lw=1, ms=10)
    if linstats_27.pvalue < 0.05:
        ax.plot(tl27_df[col], yline_27, '-b', lw=3)
    else:
        ax.plot(tl27_df[col], yline_27, '-',color='gray', lw=1, ms=10)
    ax.set_xlabel(col)
    print(f"{col} | TL47 p={linstats_47.pvalue:.3g}, TL27 p={linstats_27.pvalue:.3g}")
for ax in axes[[0]]:
    ax.set_ylabel('Soil organic carbon stocks (kg m$^{-2}$)')

axes[0].legend(loc='best')
plt.tight_layout()
plt.savefig(f'{figoutdir}/CstocksCovariates.jpg', dpi=300)
plt.show()


#########
#Plotting carbon and covars for combined sites
#########

fig, axes = plt.subplots(1,3, figsize=(10, 5), sharey=True)
axes = axes.flatten()
for ax, col in zip(axes, cols):
    x = df_sites[col]
    y = df_sites.Total_C_stock_kg_m2_0to50cm
    linstats_carb = linregress(x,y)
    yline_disp = linstats_carb.slope * x + linstats_carb.intercept
    ax.plot(x,y, 'ok')
    print(f"{col} | p={linstats_carb.pvalue:.3g}")
    if linstats_carb.pvalue < 0.05:
        ax.plot(x, yline_disp, '-k', lw=3)
    else:
        ax.plot(x, yline_disp, '-',color='gray', lw=1, ms=10)

    ax.set_xlabel(col)

for ax in axes[[0]]:
    ax.set_ylabel('Soil organic carbon stocks (kg m$^{-2}$)')

plt.tight_layout()
plt.savefig(f'{figoutdir}/CstocksCovariates_BothSites.jpg', dpi=300)
plt.show()


#########
#Plotting disp and covars for combined sites
#########

fig, axes = plt.subplots(1,3, figsize=(10, 5), sharey=True)
axes = axes.flatten()
for ax, col in zip(axes, cols):
    x = df_sites[col]
    y = df_sites.mean_disp_rate
    linstats_carb = linregress(x,y)
    yline_disp = linstats_carb.slope * x + linstats_carb.intercept
    ax.plot(x,y, 'ok')
    print(f"{col} | p={linstats_carb.pvalue:.3g}")
    if linstats_carb.pvalue < 0.05:
        ax.plot(x, yline_disp, '-k', lw=3)
    else:
        ax.plot(x, yline_disp, '-',color='gray', lw=1, ms=10)

    ax.set_xlabel(col)

for ax in axes[[0]]:
    ax.set_ylabel('Horizontal displacement (m yr${-1}$)')

plt.tight_layout()
#plt.savefig(f'{figoutdir}/DisplacementCovariates_BothSites.jpg', dpi=300)
plt.show()





########################
#########
#########################
ycol='Total_C_stock_kg_m2_0to50cm'
ycollabel='Soil organic carbon stock (kg m$^{-2}$)'
linreg_47 = linregress(df_sites.loc[df_sites['Site'] == 'TL47','mean_disp_rate'],
df_sites.loc[df_sites['Site'] == 'TL47',ycol])
fig,ax = plt.subplots(figsize=(5,5))
ax.plot(
    df_sites.loc[df_sites['Site'] == 'TL47','mean_disp_rate'],
    df_sites.loc[df_sites['Site'] == 'TL47',ycol],'ok',label='Teller 47')

ax.plot(
    df_sites.loc[df_sites['Site'] == 'TL47','mean_disp_rate'],
    linreg_47.slope*df_sites.loc[df_sites['Site'] == 'TL47','mean_disp_rate']+linreg_47.intercept,'-k'

)
ax.text(0.02,13,f'r$^{2}$={linreg_47.rvalue**2:.2f}')
ax.set_xlabel('Horizontal displacement rate (m yr$^{-1}$)')
ax.set_ylabel(ycollabel)
plt.tight_layout()
plt.savefig(f'{figoutdir}/TL47StocksDisp.jpg',dpi=300)
plt.show()


linreg_27 = linregress(df_sites.loc[df_sites['Site'] == 'TL27','mean_disp_rate'],
df_sites.loc[df_sites['Site'] == 'TL27',ycol])

fig,ax = plt.subplots(figsize=(5,5))
ax.plot(
    df_sites.loc[df_sites['Site'] == 'TL27','mean_disp_rate'],
    df_sites.loc[df_sites['Site'] == 'TL27',ycol],'ob',label='Teller 27')

ax.plot(
    df_sites.loc[df_sites['Site'] == 'TL27','mean_disp_rate'],
    linreg_27.slope*df_sites.loc[df_sites['Site'] == 'TL27','mean_disp_rate']+linreg_27.intercept,'-b'

)
ax.set_xlabel('Horizontal displacement rate (m yr$^{-1}$)')
ax.set_ylabel(ycollabel)
ax.text(0.016,25,f'r$^{2}$={linreg_27.rvalue**2:.2f}')
plt.tight_layout()
plt.savefig(f'{figoutdir}/TL27StocksDisp.jpg',dpi=300)
plt.show()

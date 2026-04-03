import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
plt.rcParams.update({'font.size': 14})
df_sites = pd.read_csv('/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/wDisplacement/TL47TL27StocksDisplacement_covariates_CNratio.csv')
figoutdir = '/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/figs'
def lognormal_hump(x, a, mu, sigma):
    return a * np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
def linear_exp(x, x0, y0, m, k):
    return np.where(
        x <= x0,
        m * (x - x0) + y0,
        y0 * np.exp(-k * (x - x0))
    )

def soil_prod_hump(x, a, b):
    return a * x * np.exp(-b * x)
def fitDispAllData(df,outfig,xcol='mean_disp_rate',ycol='Total_C_stock_kg_m2_0to50cm',ycollabel='Soil organic carbon stock (kg m$^{-2}$)',scat_col = 'C_N_ratio',scat_col_label = 'C:N',plotline=True):
    x = df[xcol]*0.5
    y = df[ycol]
    x0_init = x[np.argmax(y)]# pick the y peak location for x
    y0_init = y.max()
    fig,ax = plt.subplots(figsize=(7,6))
    if plotline:
        

        
        mask = x > 0  # log requires positive x

        p0 = [np.max(y), np.log(np.median(x)), 0.5]

        params, cov = curve_fit(
            lognormal_hump, x[mask], y[mask], p0=p0
        )

        xfit = np.linspace(x[mask].min(), x.max(), 200)
        yfit = lognormal_hump(xfit, *params)
        ax.plot(xfit, yfit, color='k', lw=2)


    scat1 = ax.scatter(df_sites[xcol]*0.5,df_sites[ycol],c= df_sites[scat_col],s=df_sites['slope']*20,cmap='viridis',edgecolor='k')
    #scat2 = ax.scatter(tl27_disp,tl27_ycol,c = tl27_scatcol, s=100,cmap='viridis')
    cbar = plt.colorbar(scat1,label=scat_col_label)



    handles, labels = scat1.legend_elements(prop="sizes", alpha=0.6, num=4)

    # Convert labels back to slope (divide by 10)
    # Strip LaTeX formatting and rescale
    clean_labels = [
        float(l.replace('$\\mathdefault{', '').replace('}$', '')) / 20
        for l in labels
    ]

    ax.legend(handles, [f"{round(l)}" for l in clean_labels], title="Slope")

    ax.set_ylabel(ycollabel)
    ax.set_xlabel('Soil flux (m$^{2}$ yr$^{-1}$)')
    plt.tight_layout()
    #ax.legend()
    plt.savefig(outfig,dpi=300)
    plt.show()


fitDispAllData(df_sites,f'{figoutdir}/TL27CstockDisplacementCombinedSites_flux_.jpg',plotline=False)


import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
import geopandas as gpd
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams.update({'font.size': 14})
df_sites = gpd.read_file('FinalCleanedFiles/wDisplacement/Tl47Tl27CarbonDisplacementTempsTopo.gpkg')
figoutdir = 'figs'
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
# def fitDispAllData(df,outfig,xcol='mean_disp_rate',ycol='Total_C_stock_kg_m2_0to50cm',ycollabel='Soil organic carbon stock (kg m$^{-2}$)',scat_col = 'C_N_ratio',scat_col_label = 'C:N',plotline=True):
#     x = df[xcol]*0.5
#     y = df[ycol]
#     # x0_init = x[np.argmax(y)]# pick the y peak location for x
#     # y0_init = y.max()
#     fig,ax = plt.subplots(1,2,figsize=(7,6))
#     # if plotline:
        

        
#     #     mask = x > 0  # log requires positive x

#     #     p0 = [np.max(y), np.log(np.median(x)), 0.5]

#     #     params, cov = curve_fit(
#     #         lognormal_hump, x[mask], y[mask], p0=p0
#     #     )

#     #     xfit = np.linspace(x[mask].min(), x.max(), 200)
#     #     yfit = lognormal_hump(xfit, *params)
#     #     ax.plot(xfit, yfit, color='k', lw=2)


#     scat1 = ax.scatter(df_sites[xcol]*0.5,df_sites[ycol],c= df_sites[scat_col],s=df_sites['slope']*20,cmap='viridis',edgecolor='k')
#     #scat2 = ax.scatter(tl27_disp,tl27_ycol,c = tl27_scatcol, s=100,cmap='viridis')
#     cbar = plt.colorbar(scat1,label=scat_col_label)



#     handles, labels = scat1.legend_elements(prop="sizes", alpha=0.6, num=4)

#     # Convert labels back to slope (divide by 10)
#     # Strip LaTeX formatting and rescale
#     clean_labels = [
#         float(l.replace('$\\mathdefault{', '').replace('}$', '')) / 20
#         for l in labels
#     ]

#     ax.legend(handles, [f"{round(l)}" for l in clean_labels], title="Slope")

#     ax.set_ylabel(ycollabel)
#     ax.set_xlabel('Soil flux (m$^{2}$ yr$^{-1}$)')
#     plt.tight_layout()
#     #ax.legend()
#     #plt.savefig(outfig,dpi=300)
#     plt.show()


# fitDispAllData(df_sites,f'{figoutdir}/CstockDisplacementCombinedSites_flux_.jpg',plotline=False)


def fitDispAllData(
    df,
    outfig,
    xcol='mean_disp_rate',
    ycol='Total_C_stock_kg_m2_0to50cm',
    ycollabel='Soil organic carbon stock (kg m$^{-2}$)',
    scat_col1='C_N_ratio',
    scat_col_label1='C:N mean of interval',
    scat_col2='mean_temp_50cm',
    scat_col_label2='Mean annual temperature ($^{\\circ}$C) 50 cm depth',
    plotline=True
):

    # ==========================================================
    # DATA
    # ==========================================================

    x_flux = df[xcol].astype(float) * 0.5
    y_stock = df[ycol].astype(float)
    y_cn = df[scat_col1].astype(float)
    # ==========================================================
    # FIGURE
    # ==========================================================

    fig, axes = plt.subplots(
        1, 2,
        figsize=(12, 6),
        sharex=True
    )

    ax1, ax2 = axes

    # ==========================================================
    # LEFT PANEL — SOC STOCK VS FLUX
    # ==========================================================

    scat1 = ax1.scatter(
        x_flux,
        y_stock,
        c=df[scat_col1],
        s=df['slope'] * 20,
        cmap='viridis',
        edgecolor='k'
    )


    ax1.set_ylabel(ycollabel)
    ax1.set_xlabel('Soil flux (m$^{2}$ yr$^{-1}$)')

    # panel label
    ax1.text(
        0.02, 0.98,
        'a',
        transform=ax1.transAxes,
        fontsize=14,
        fontweight='bold',
        va='top'
    )

    # ==========================================================
    # RIGHT PANEL — C:N VS FLUX
    # ==========================================================

    scat2 = ax2.scatter(
        x_flux,
        y_cn,
        c=df[scat_col2],
        s=100,
        cmap ='Blues_r',
        edgecolor='k'
    )

    # optional regression
    if plotline:
        linstats2 = linregress(x_flux, y_cn)
        yline2 = linstats2.slope * x_flux + linstats2.intercept

        ax2.plot(x_flux, yline2, '-k', lw=2)

        ax2.text(
            0.05, 0.95,
            f'$r^2$ = {linstats2.rvalue**2:.2f}\np = {linstats2.pvalue:.3f}',
            transform=ax2.transAxes,
            va='top'
        )

    ax2.set_ylabel('C:N mean of interval')
    ax2.set_xlabel('Soil flux (m$^{2}$ yr$^{-1}$)')

    # panel label
    ax2.text(
        0.02, 0.98,
        'b',
        transform=ax2.transAxes,
        fontsize=14,
        fontweight='bold',
        va='top'
    )

    # ==========================================================
    # COLORBAR
    # ==========================================================
    divider1 = make_axes_locatable(ax1)
    divider2 = make_axes_locatable(ax2)
    # cbar1 = fig.colorbar(
    #     scat1,
    #     ax=ax1,
    #     label=scat_col_label1,
    #     shrink=0.9
    # )
    cax1 = divider1.append_axes(
    "top",
    size="5%",
    pad=0.1
)

    cbar1 = fig.colorbar(
        scat1,
        cax=cax1,
        orientation="horizontal"
    )

    cbar1.set_label(scat_col_label1)
    cax1.xaxis.set_ticks_position('top')
    cax1.xaxis.set_label_position('top')

    # cbar2 = fig.colorbar(
    #     scat2,
    #     ax=ax2,
    #     label=scat_col_label2,
    #     shrink=0.9
    # )
    cax2 = divider2.append_axes(
        "top",
        size="5%",
        pad=0.1
    )
    cbar2 = fig.colorbar(
        scat2,
        cax=cax2,
        orientation='horizontal'
    )
    cbar2.set_label(scat_col_label2)
    cax2.xaxis.set_ticks_position('top')
    cax2.xaxis.set_label_position('top')
    # ==========================================================
    # SIZE LEGEND
    # ==========================================================

    handles, labels = scat1.legend_elements(
        prop="sizes",
        alpha=0.6,
        num=4
    )

    clean_labels = [
        float(l.replace('$\\mathdefault{', '').replace('}$', '')) / 20
        for l in labels
    ]

    ax1.legend(
        handles,
        [f"{round(l)}" for l in clean_labels],
        title="Slope"
    )

    # ==========================================================
    # FINAL FORMATTING
    # ==========================================================

    plt.tight_layout()

    plt.savefig(
        outfig,
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()


fitDispAllData(
    df_sites,   
    scat_col2='days_unfrozen_below_50cm',
    scat_col_label2='Unfrozen days below 50 cm depth',
    outfig=f'{figoutdir}/Cstock_CN_DisplacementCombinedSites_flux.jpg',
    plotline=True
)
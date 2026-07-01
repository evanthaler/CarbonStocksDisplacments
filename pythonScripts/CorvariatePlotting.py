import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
import pandas as pd
import geopandas as gpd
import numpy as np
plt.rcParams.update({'font.size': 14})

df_sites = gpd.read_file('FinalCleanedFiles/wDisplacement/Tl47Tl27CarbonDisplacementTempsTopo.gpkg')
df_curvgrids= pd.read_csv('FinalCleanedFiles/wDisplacement/TL47TL27StocksDisplacement_multicurv.csv')
df_curvgrids = df_curvgrids[['dem_curv_gauss_sigma9m_laplacian',
       'dem_curv_gauss_sigma15m_laplacian',
        'dem_curv_gauss_sigma30m_laplacian','SampleLocationName']]

df_sites = df_sites.merge(df_curvgrids,on='SampleLocationName')

df_sites['mean_disp_rate']=df_sites.mean_disp_rate.astype(float)
df_sites['Total_C_stock_kg_m2_0to50cm']=df_sites.Total_C_stock_kg_m2_0to50cm.astype(float)

cols_to_plot = ['slope','curv','drainagearea']
cols = ['Slope (degrees)',
            'Curvature (m$^{-1}$)',
            'Drainage area (m$^{2}$)']

labels = ['a', 'b', 'c', 'd', 'e', 'f']

save_fig = True



# ###################################################################
# Plotting displacement + C stocks together (2 rows x 3 columns)
# ###################################################################

# create figure
fig, axes = plt.subplots(
    2, 3,
    figsize=(8,8),
    sharex='col'
)

axes = axes.flatten()

# ==========================================================
# TOP ROW — Displacement
# ==========================================================

for i, (ax, col, colname, label) in enumerate(
    zip(axes[:3], cols_to_plot, cols, labels[:3])
):

    x = df_sites[col]
    y = df_sites.mean_disp_rate*0.5

    linstats = linregress(x, y)
    yline = linstats.slope * x + linstats.intercept

    # points
    ax.plot(x, y, 'ok')

    # regression line
    if linstats.pvalue < 0.05:
        ax.plot(x, yline, '-k', lw=3)
    else:
        ax.plot(x, yline, '-', color='gray', lw=1)

    # stats text
    ax.text(
        0.05, 0.90,
        f'$r^2$ = {linstats.rvalue**2:.2f}\np = {linstats.pvalue:.3f}',
        transform=ax.transAxes,
        va='top'
    )

    # subplot label
    ax.text(
        0.02, 0.98,
        label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight='bold',
        va='top'
    )

    # x labels only on bottom row
    

# y-axis label
axes[0].set_ylabel('Soil flux (m$^{2}$ yr$^{-1}$)')


# ==========================================================
# BOTTOM ROW — Carbon stocks
# ==========================================================

for i, (ax, col, colname, label) in enumerate(
    zip(axes[3:], cols_to_plot, cols, labels[3:])
):

    x = df_sites[col]
    y = df_sites.Total_C_stock_kg_m2_0to50cm

    linstats = linregress(x, y)
    yline = linstats.slope * x + linstats.intercept

    # points
    ax.plot(x, y, 'ok')

    # regression line
    if linstats.pvalue < 0.05:
        ax.plot(x, yline, '-k', lw=3)
    else:
        ax.plot(x, yline, '-', color='gray', lw=1)

    # stats text
    ax.text(
        0.05, 0.90,
        f'$r^2$ = {linstats.rvalue**2:.2f}\np = {linstats.pvalue:.3f}',
        transform=ax.transAxes,
        va='top'
    )

    # subplot label
    ax.text(
        0.02, 0.98,
        label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight='bold',
        va='top'
    )
    
    # x-axis labels
    ax.set_xlabel(colname)

# y-axis label
axes[3].set_ylabel('SOC stock (kg m$^{-2}$)')


# ==========================================================
# Final formatting
# ==========================================================
for ax in axes:
    ax.set_box_aspect(1)
#plt.tight_layout()
if save_fig:
    plt.savefig(
        'figs/Displacement_CStock_Covariates.jpg',
        dpi=300,
        bbox_inches='tight'
    )

plt.show()
plt.close()









# fig, axes = plt.subplots(1,3, figsize=(10, 5), sharey=True)
# axes = axes.flatten()

# for ax, col, colname, label in zip(axes, cols_to_plot, cols, labels):

#     x, y = df_sites[col], df_sites.C_N_ratio
#     linstats = linregress(x, y)

#     yline = linstats.slope * x + linstats.intercept

#     ax.plot(x, y, 'ok')

#     # R2 and p-value
#     ax.text(
#         0.05, 0.9,
#         f'$r^2$ = {linstats.rvalue**2:.2f}\np = {linstats.pvalue:.3f}',
#         transform=ax.transAxes,
#         verticalalignment='top'
#     )

#     # subplot label
#     ax.text(
#         0.02, 0.98,
#         f'{label}',
#         transform=ax.transAxes,
#         fontsize=14,
#         fontweight='bold',
#         va='top'
#     )

#     if linstats.pvalue < 0.05:
#         ax.plot(x, yline, '-k', lw=3)
#     else:
#         ax.plot(x, yline, '-', color='gray', lw=1, ms=10)

#     ax.set_xlabel(colname)

# axes[0].set_ylabel('C:N')

# plt.tight_layout()
# #plt.savefig(f'{figoutdir}/C_NCovariates.jpg', dpi=300)
# plt.show()
# plt.close()



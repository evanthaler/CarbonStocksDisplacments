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
labels = ['a', 'b', 'c']

cols_to_plot = ['slope','curv','drainagearea']
cols = ['Slope (degrees)',
            'Curvature (m$^{-1}$)',
            'Drainage area (m$^{2}$)']

fig, axes = plt.subplots(1,3, figsize=(10, 5), sharey=True)
axes = axes.flatten()

labels = ['a', 'b', 'c']

for ax, col, colname, label in zip(axes, cols_to_plot, cols, labels):

    x, y = df_sites[col], df_sites.mean_disp_rate
    linstats = linregress(x, y)

    yline = linstats.slope * x + linstats.intercept

    ax.plot(x, y, 'ok')

    # R2 and p-value
    ax.text(
        0.05, 0.9,
        f'$R^2$ = {linstats.rvalue**2:.2f}\np = {linstats.pvalue:.3f}',
        transform=ax.transAxes,
        verticalalignment='top'
    )

    # subplot label
    ax.text(
        0.02, 0.98,
        f'{label}',
        transform=ax.transAxes,
        fontsize=14,
        fontweight='bold',
        va='top'
    )

    if linstats.pvalue < 0.05:
        ax.plot(x, yline, '-k', lw=3)
    else:
        ax.plot(x, yline, '-', color='gray', lw=1, ms=10)

    ax.set_xlabel(colname)
axes[0].set_ylabel('Soil flux (m$^{2}$ yr$^{-1}$)')


plt.tight_layout()
plt.savefig(f'{figoutdir}/DisplacementCovariates.jpg', dpi=300)
plt.show()

######################################


fig, axes = plt.subplots(1,3, figsize=(10, 5), sharey=True)
axes = axes.flatten()

for ax, col, colname, label in zip(axes, cols_to_plot, cols, labels):

    x, y = df_sites[col], df_sites.Total_C_stock_kg_m2_0to50cm
    linstats = linregress(x, y)

    yline = linstats.slope * x + linstats.intercept

    ax.plot(x, y, 'ok')
    # R2 and p-value
    ax.text(
        0.05, 0.9,
        f'$r^2$ = {linstats.rvalue**2:.2f}\np = {linstats.pvalue:.3f}',
        transform=ax.transAxes,
        verticalalignment='top'
    )

    # subplot label
    ax.text(
        0.02, 0.98,
        f'{label}',
        transform=ax.transAxes,
        fontsize=14,
        fontweight='bold',
        va='top'
    )

    if linstats.pvalue < 0.05:
        ax.plot(x, yline, '-k', lw=3)
    else:
        ax.plot(x, yline, '-', color='gray', lw=1, ms=10)

    ax.set_xlabel(colname)

axes[0].set_ylabel('Soil organic carbon stock (kg m$^{-2}$)')

plt.tight_layout()
plt.savefig(f'{figoutdir}/CStocksCovariates.jpg', dpi=300)
plt.show()



fig, axes = plt.subplots(1,3, figsize=(10, 5), sharey=True)
axes = axes.flatten()

for ax, col, colname, label in zip(axes, cols_to_plot, cols, labels):

    x, y = df_sites[col], df_sites.C_N_ratio
    linstats = linregress(x, y)

    yline = linstats.slope * x + linstats.intercept

    ax.plot(x, y, 'ok')

    # R2 and p-value
    ax.text(
        0.05, 0.9,
        f'$r^2$ = {linstats.rvalue**2:.2f}\np = {linstats.pvalue:.3f}',
        transform=ax.transAxes,
        verticalalignment='top'
    )

    # subplot label
    ax.text(
        0.02, 0.98,
        f'{label}',
        transform=ax.transAxes,
        fontsize=14,
        fontweight='bold',
        va='top'
    )

    if linstats.pvalue < 0.05:
        ax.plot(x, yline, '-k', lw=3)
    else:
        ax.plot(x, yline, '-', color='gray', lw=1, ms=10)

    ax.set_xlabel(colname)

axes[0].set_ylabel('C:N')

plt.tight_layout()
plt.savefig(f'{figoutdir}/C_NCovariates.jpg', dpi=300)
plt.show()







# ###################################################################
# Plotting displacement + C stocks together (2 rows x 3 columns)
# ###################################################################



cols_to_plot = ['slope', 'curv', 'drainagearea']

cols = [
    'Slope (degrees)',
    'Curvature (m$^{-1}$)',
    'Drainage area (m$^{2}$)'
]

# panel labels
labels = ['a', 'b', 'c', 'd', 'e', 'f']

# create figure
fig, axes = plt.subplots(
    2, 3,
    figsize=(12, 12),
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
    y = df_sites.mean_disp_rate

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

plt.tight_layout()

plt.savefig(
    f'{figoutdir}/Displacement_CStock_Covariates.jpg',
    dpi=300,
    bbox_inches='tight'
)

plt.show()
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal
import rasterio
from rasterio.plot import show
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


# --------------------------------------------------
# PATHS (EDIT THESE)
# --------------------------------------------------
rasterpath = '/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/Rasters/TL27/hs.tif'
contourpath = '/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/Rasters/TL27/elev.tif'
joelptspath = '/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/wDisplacement/TL27StocksDisplacement.csv'


# --------------------------------------------------
# READ & PREP VECTOR DATA
# --------------------------------------------------
def vector_df():
    """
    Reads a single CSV containing displacement, direction, and SOC
    """

    df = pd.read_csv(joelptspath)

    # Convert azimuth (deg clockwise from north) → radians
    df['angle_rads'] = df['aspect'] * (math.pi / 180)

    return df


# --------------------------------------------------
# MAIN PLOTTING FUNCTION
# --------------------------------------------------
def displacement_slope_fig(ax):

    # ---- Read slope raster
    ds = gdal.Open(rasterpath)
    gt = ds.GetGeoTransform()

    xmin, ymax = gt[0], gt[3]
    xres, yres = gt[1], gt[5]
    cols, rows = ds.RasterXSize, ds.RasterYSize

    xmax = xmin + xres * cols
    ymin = ymax + yres * rows
    extent = (xmin, xmax, ymin, ymax)

    elev = ds.GetRasterBand(1).ReadAsArray().astype(float)

    # ---- Elevation contours
    elevds = gdal.Open(contourpath)
    elev = elevds.GetRasterBand(1).ReadAsArray().astype(float)

    levels = np.arange(
        np.floor(np.nanmin(elev) / 50) * 50,
        np.ceil(np.nanmax(elev) / 50) * 50 + 50,
        50
    )

    xcoords = np.arange(cols) * xres + xmin
    ycoords = np.arange(rows) * yres + ymax

    # ---- Vector dataframe
    vectors = vector_df()

    # Vector components
    mag = vectors['mean_disp_rate'] * 5000  # scale for visibility
    ang = vectors['angle_rads']

    vectors['U'] = mag * np.sin(ang)  # east-west
    vectors['V'] = mag * np.cos(ang)  # north-south

    # --------------------------------------------------
    # PLOTTING
    # --------------------------------------------------


    with rasterio.open(rasterpath) as ds:
        show(ds, ax=ax, cmap="gray", alpha=1.0, zorder=-1)
        bounds = ds.bounds 

    # cbar_slope = plt.colorbar(slope_im, ax=ax, shrink=0.8)
    # cbar_slope.set_label('Slope (degrees)')

    # ---- Elevation contours
    # ctr = ax.contour(
    #     xcoords,
    #     ycoords,
    #     elev,
    #     levels=levels,
    #     colors='k',
    #     linewidths=0.5,
    #     alpha=0.5,
    #     zorder=2
    # )
    # ax.clabel(ctr, fontsize=9)

    # ---- SOC-colored points
    sc = ax.scatter(
        vectors['X'],
        vectors['Y'],
        c=vectors['Total_C_stock_kg_m2_0to50cm'],
        cmap='viridis',
        s=60,
        edgecolor='k',
        zorder=4
    )

    cbar_soc = plt.colorbar(sc, ax=ax, shrink=0.8)
    cbar_soc.set_label('SOC stock (kg C m$^{-2}$)')

    # ---- Displacement vectors
    ax.quiver(
        vectors['X'],
        vectors['Y'],
        vectors['U'],
        vectors['V'],
        angles='xy',
        scale_units='xy',
        scale=1,
        width=0.005,
        headwidth=3,
        headlength=4,
        color='k',
        alpha=0.85,
        zorder=5
    )

    # ---- Vector scale arrow
    pad_x = 0.5 * (bounds.right - bounds.left)
    pad_y = 0.05 * (bounds.top - bounds.bottom)

    scale_len = 0.1 * 5000
    sx, sy = bounds.left+pad_x, bounds.bottom+pad_y

    bbox = FancyBboxPatch(
        (sx - 20, sy - 40),
        scale_len + 260,
        70,
        boxstyle='round,pad=0.5',
        facecolor='none',
        edgecolor='none',
        alpha=0.85,
        zorder=6
    )
    ax.add_patch(bbox)

    arrow = FancyArrowPatch(
        (sx, sy),
        (sx + scale_len, sy),
        lw=0.7,
        color='k',
        zorder=7
    )
    ax.add_patch(arrow)

    ax.text(
        sx + scale_len / 2,
        sy + 20,
        '10 cm yr$^{-1}$',
        ha='center',
        va='bottom',
        fontsize=9,
        zorder=8
    )

    # ---- Axes formatting
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')

    ax.get_xaxis().set_major_formatter(ScalarFormatter(useOffset=False))
    ax.get_yaxis().set_major_formatter(ScalarFormatter(useOffset=False))
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)


    for lab in ax.get_yticklabels():
        lab.set_rotation(45)
    ax.tick_params(axis='both', labelrotation=45)


    ax.set_xlim(bounds.left, bounds.right)
    ax.set_ylim(bounds.bottom, bounds.top)

    # ax.set_xlim(442250, 443350)
    # ax.set_ylim(7206000, 7207200)


# --------------------------------------------------
# RUN
# --------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6), tight_layout=True)
displacement_slope_fig(ax)

#plt.savefig('/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/figs/TL27displacementvectors_carbonstocks.jpg', dpi=300, bbox_inches='tight')


plt.show()

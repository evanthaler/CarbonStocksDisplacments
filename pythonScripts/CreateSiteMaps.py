import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import contextily as cx
import xyzservices.providers as xyz
import numpy as np
from osgeo import gdal
import rasterio
from rasterio.plot import show

def load_tiff(inputfile):
    tif = gdal.Open(inputfile)
    gt = tif.GetGeoTransform()
    dat = tif.GetRasterBand(1)
    nodat = dat.GetNoDataValue()
    dat = np.array(dat.ReadAsArray());dat=dat.astype('float32')
    dat[dat==nodat]=np.nan
    tif=None
    return dat, gt

def plotstudymap(hillshadepath,dempath=None,soilpath=None,legendlabel=None,radiocarbonpath=None,legendlabel2=None,statepath=None,state=None,plotpoints=True,outpath=None,dispath=None,
xlim_min=None,xlim_max=None,ylim_min = None, ylim_max=None,
contour_interval=None,contour_levels=None,contour_color='saddlebrown',contour_linewidth=0.5,contour_labels=False,
inset=True):
    if statepath is not None:
        us_shapefile = gpd.read_file(statepath)
        state_shp = us_shapefile[us_shapefile.STUSPS==state].to_crs('epsg:32603')

    if soilpath:
        soilpoints = gpd.read_file(soilpath).to_crs('epsg:32603')

    #xmin, ymin, xmax, ymax = soilpoints.total_bounds

    fig, ax = plt.subplots(figsize=(6, 6))
    if plotpoints:
        soilpoints.plot(ax=ax, color='black',marker='o', edgecolor='black', markersize=50, zorder=0,label='Soil carbon sample')
    if radiocarbonpath is not None:
        radiopoints = gpd.read_file(radiocarbonpath).to_crs('epsg:32603')
        radiopoints.plot(ax=ax, color='darkred', marker='s', edgecolor='none',markersize=50, zorder=0,label='Radiocarbon sample')
    if dispath is not None:
        dispoints = gpd.read_file(dispath).to_crs('epsg:32603')
        dispoints.plot(ax=ax,color='blue',marker='o',edgecolor='darkblue',markersize=50,zorder=1,label='Soil carbon sample\n and flux measurement')
    if inset:
        ax_inset = inset_axes(ax, width="20%", height="20%", loc='lower left', borderpad=1)
        state_shp.plot(ax=ax_inset, color='white', edgecolor='black', linewidth=0.8)
        soilpoints.iloc[[0]].plot(
            ax=ax_inset,
            color='red',
            edgecolor='red',
            marker='*',
            markersize=100,
            zorder=5
        )
    #ax.legend(loc='best')
    # ax_inset.set_xticks([]); ax_inset.set_yticks([]) #ax_inset.set_title("Oregon", fontsize=8)

    # --- Style main map ---
    #ax.axis('off')
    # padding = 10 
    # ax.set_xlim(xmin - padding, xmax + padding)
    # ax.set_ylim(ymin - padding, ymax + padding)
    ax.set_xlabel('Easting')
    ax.set_ylabel('Northing')
    ax.ticklabel_format(style='plain')
    ax.tick_params(axis='both', labelrotation=45)
    #if legendlabel is not None:
    ax.legend(fontsize=14, loc='upper right')
    with rasterio.open(hillshadepath) as ds:
        show(ds, ax=ax, cmap="gray", alpha=1.0, zorder=-1)
        raster_left, raster_bottom, raster_right, raster_top = ds.bounds
    ax.set_xlim(xlim_min if xlim_min is not None else raster_left, xlim_max if xlim_max is not None else raster_right)
    ax.set_ylim(ylim_min if ylim_min is not None else raster_bottom, ylim_max if ylim_max is not None else raster_top)
    if dempath is not None:
        dem, dem_gt = load_tiff(dempath)
        nrows, ncols = dem.shape
        x_coords = dem_gt[0] + dem_gt[1] * (np.arange(ncols) + 0.5)
        y_coords = dem_gt[3] + dem_gt[5] * (np.arange(nrows) + 0.5)
        X, Y = np.meshgrid(x_coords, y_coords)
        if contour_levels is not None:
            levels = contour_levels
        elif contour_interval is not None:
            vmin, vmax = np.nanmin(dem), np.nanmax(dem)
            levels = np.arange(np.floor(vmin / contour_interval) * contour_interval, vmax + contour_interval, contour_interval)
        else:
            levels = 10
        cs = ax.contour(X, Y, dem, levels=levels, colors=contour_color, linewidths=contour_linewidth, zorder=0)
        if contour_labels:
            ax.clabel(cs, inline=True, fontsize=14, fmt='%d')

    plt.tight_layout()
    if outpath is not None:
        plt.savefig(outpath,dpi=300,bbox_inches='tight')
    
    plt.show()
    plt.close()




kr64_elev='/Users/evanthaler/Documents/Projects/permafrost/SewardLidar/KG64/kr64_elev.tif'
kr64_hs = "/Users/evanthaler/Documents/Projects/permafrost/SewardLidar/KG64/kr64_hillshade.tif"
kr86_elev="/Users/evanthaler/Documents/Projects/permafrost/SewardLidar/KGRiver/KR86_elev.tif"
kr86_hs="/Users/evanthaler/Documents/Projects/permafrost/SewardLidar/KGRiver/KR86_hillshade.tif"
tl27_hs="/Users/evanthaler/Documents/Projects/permafrost/SewardLidar/TL27/TL27_site_hillshade.tif"
tl27_elev= "/Users/evanthaler/Documents/Projects/permafrost/SewardLidar/TL27/tl27_elev.tif"
tl47_hs="/Users/evanthaler/Documents/Projects/permafrost/SewardLidar/TL47/hs.tif"
tl47_elev="/Users/evanthaler/Documents/Projects/permafrost/SewardLidar/TL47/tl47_elev.tif"

plotstudymap(
hillshadepath = tl47_hs,
dempath = tl47_elev,
soilpath = '/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/UltimateSheet/TL47TL27KR_UltimateSheet.gpkg',
radiocarbonpath='/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/UltimateSheet/sewardpen14c_all_samples_ever_utm.gpkg',
dispath = '/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/wDisplacement/TL47TL27StocksDisplacement_multicurv.gpkg',
#legendlabel= 'Soil carbon stocks\nand displacement\nmeasurements',
plotpoints = False,
contour_interval = 10,
contour_labels = False,

outpath = '/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/figs/TL47StudyMap.jpg',
xlim_min = 442000,xlim_max=443500,#TL47

# xlim_min = 508500,xlim_max=508900,#KG64,
# ylim_min =7226750 , ylim_max =7227200,#KG64

# xlim_min = 513000,xlim_max=514500,#KG86,
# ylim_min =7256000 , ylim_max =7257000,#KG86

# xlim_min = 453500,xlim_max=456000,#TL27,
# ylim_min =7178500 , ylim_max =7180300,#TL27
inset=False)








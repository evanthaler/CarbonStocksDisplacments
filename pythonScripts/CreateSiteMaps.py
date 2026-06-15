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

def plotstudymap(dempath,soilpath=None,legendlabel=None,radiocarbonpath=None,legendlabel2=None,statepath=None,state=None,plotpoints=True,pfraster = None,outpath=None,dispath=None,
xlim_min=None,xlim_max=None,ylim_min = None, ylim_max=None,
inset=True): 
    if statepath is not None:
        us_shapefile = gpd.read_file(statepath)
        state_shp = us_shapefile[us_shapefile.STUSPS==state].to_crs('epsg:32603')

    soilpoints = gpd.read_file(soilpath).to_crs('epsg:32603')

    #xmin, ymin, xmax, ymax = soilpoints.total_bounds

    fig, ax = plt.subplots(figsize=(6, 6))
    if plotpoints:
        soilpoints.plot(ax=ax, color='black',marker='o', edgecolor='black', markersize=50, zorder=0,label='Soil carbon sample')
    if radiocarbonpath is not None:
        radiopoints = gpd.read_file(radiocarbonpath).to_crs('epsg:32603')
        radiopoints.plot(ax=ax, color='darkred', marker='*', edgecolor='none',markersize=30, zorder=0,label='Radiocarbon sample')
    if dispath is not None:
        dispoints = gpd.read_file(dispath).to_crs(state_shp.crs)
        dispoints.plot(ax=ax,color='red',edgecolor='darkred',markersize=12,zorder=1,label='Displacement measurement')
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
    if xlim_min is not None:
        ax.set_xlim(xlim_min,xlim_max)
    if ylim_min is not None:
        ax.set_ylim(ylim_min,ylim_max)
    #if legendlabel is not None:
    ax.legend(fontsize=8, loc='upper left')
    with rasterio.open(dempath) as ds:
        show(ds, ax=ax, cmap="gray", alpha=1.0, zorder=-1)
    if pfraster is not None:
            img,img_gt = load_tiff(pfraster)
            img_mask= np.where(img == 1, img, np.nan)
            extent = (img_gt[0], img_gt[0] + img.shape[1] * img_gt[1], img_gt[3] + img.shape[0] * img_gt[5], img_gt[3])
            plt.imshow(img_mask,cmap=plt.cm.Greens_r,alpha=0.5,extent=extent,interpolation='none',zorder=-1)
    plt.tight_layout()
    if outpath is not None:
        plt.savefig(outpath,dpi=300)
    
    plt.show()
    plt.close()



plotstudymap(
dempath = "/Users/evanthaler/Documents/Projects/permafrost/SewardLidar/KGRiver/KR86_hillshade.tif",
soilpath = '/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/UltimateSheet/TL47TL27KR_UltimateSheet.gpkg',
radiocarbonpath='/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/UltimateSheet/sewardpen14c_all_samples_ever_utm.gpkg',
#soilpath2 = '/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/UltimateSheet/TL47TL27KR_UltimateSheet.gpkg',
#legendlabel= 'Soil carbon stocks\nand displacement\nmeasurements',
dispath=None,
plotpoints = True,
pfraster = None,

outpath = '/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/figs/KG86_samplelocations.jpg',
#xlim_min = 442000,xlim_max=443500,#TL47
# xlim_min = 508500,xlim_max=508900,#KG64,
# ylim_min =7226750 , ylim_max =7227200,#KG64

xlim_min = 513000,xlim_max=514500,#KG86,
ylim_min =7256000 , ylim_max =7257000,#KG86
inset=False)


#plotstudymap(dempath,statepath,state,soilpath,plotpoints=True,pfraster = None,outpath=None,dispath=None,xlim_min=None,xlim_max=None,inset=True)
#/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/UltimateSheet/sewardpen14c_all_samples_ever_utm.gpkg
#/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/FinalCleanedFiles/wDisplacement/TL47TL27StocksDisplacement_covariates.gpkg
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

def plotstudymap(dempath,statepath,state,soilpath,plotpoints=True,pfraster = None,outpath=None,dispath=None,xlim_min=None,xlim_max=None,inset=True): 
    
    us_shapefile = gpd.read_file(statepath)
    state_shp = us_shapefile[us_shapefile.STUSPS==state].to_crs('epsg:32603')

    soilpoints = gpd.read_file(soilpath).to_crs(state_shp.crs)
    

    fig, ax = plt.subplots(figsize=(6, 6))
    if plotpoints:
        soilpoints.plot(ax=ax, color='blue', edgecolor='black', markersize=30, zorder=0,label='Soil carbon stocks\nand displacement\nmeasurements ')
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
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')
    ax.ticklabel_format(style='plain')
    ax.tick_params(axis='both', labelrotation=45)
    if xlim_min is not None:
        ax.set_xlim(xlim_min,xlim_max)
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



plotstudymap(
dempath = "/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/Rasters/TL27/hs.tif",
statepath='/Users/evanthaler/Documents/Projects/OSU/StateShapefiles/tl_2023_us_state/tl_2023_us_state.shp',
state='AK',
soilpath = '/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/FinalCleanedFiles/wDisplacement/TL47TL27StocksDisplacement_covariates.gpkg',
dispath='/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/disp2use/TL27TL47Displacments_PFTopo.gpkg',
plotpoints = False,
pfraster = "/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/Rasters/TL27/pf.tif",

outpath = '/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/figs/TL27DisplacementsMap_wPF.jpg',
#xlim_min = 442000,xlim_max=443500,
inset=False)


#plotstudymap(dempath,statepath,state,soilpath,plotpoints=True,pfraster = None,outpath=None,dispath=None,xlim_min=None,xlim_max=None,inset=True)
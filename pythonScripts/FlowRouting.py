import richdem as rd


# Load DEM (e.g., GeoTIFF or array)
dem = rd.LoadGDAL('/Users/evanthaler/Documents/Projects/permafrost/SewardLidar/KG64/KG64_15m.tif')

#slope = rd.TerrainAttribute(dem,attrib='slope_riserun')
#rd.SaveGDAL('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/Rasters/TL47/slope_ArcticDEM.tif',slope)
methods =['Quinn']#'Dinf','D8','Rho8',
for m in methods:
# Compute D8 flow directions
    flowacc = rd.FlowAccumulation(dem,method=m)

    # Save result
    rd.SaveGDAL(f'/Users/evanthaler/Documents/Projects/permafrost/SewardLidar/KG64/Drainage_area{m}_15m.tif', flowacc)

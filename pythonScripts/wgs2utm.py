from pyproj import Proj, Transformer
import pandas as pd
import numpy as np

def wgs84_to_utm(lat, lon, utm_zone, northern_hemisphere=True):
    """
    Convert WGS84 coordinates (latitude, longitude) to UTM coordinates in a specified zone.
    
    Parameters:
    lat (float): Latitude in decimal degrees
    lon (float): Longitude in decimal degrees
    utm_zone (int): UTM zone number (1 to 60)
    northern_hemisphere (bool): True for Northern Hemisphere, False for Southern Hemisphere
    
    Returns:
    tuple: (easting, northing)
    """
    hemisphere = "north" if northern_hemisphere else "south"
    proj_utm = Proj(proj="utm", zone=utm_zone, datum="WGS84", hemisphere=hemisphere)
    transformer = Transformer.from_proj("epsg:4326", proj_utm, always_xy=True)
    easting, northing = transformer.transform(lon, lat)
    return easting, northing


f = pd.read_csv('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/gpsData/KG64_C14Samples.csv')
x=[];y=[]
for i in np.arange(0,len(f.lat)):     

    utm_coords = wgs84_to_utm(f.lat[i], f.lon[i], 3)
    x.append(utm_coords[1]);y.append(utm_coords[0])
#Looks weird to have x-->y, but the order gets reversed by from_proj
f['X'] = np.array(y)
f['Y'] = np.array(x)
f.to_csv('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/gpsData/KG64_C14Samples_utmCoords.csv')
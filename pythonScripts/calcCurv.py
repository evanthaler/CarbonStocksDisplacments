import numpy as np
import rasterio
from scipy.signal import convolve2d
import glob
from rasterio.errors import RasterioIOError

def calc_curv_faster(tif, dx):
    try:
        out_tif = tif.replace('.tif', f'_curv{str(dx)}.tif')

        with rasterio.open(tif) as src:
            data = src.read(1).astype(np.float32)
            nodata = src.nodata
            profile = src.profile

        data[data == nodata] = np.nan

        # Simple 3x3 Laplacian kernel approximation
        kernel = np.array([[1, 1, 1],
                        [1, -8, 1],
                        [1, 1, 1]]) / (6 * dx**2)

        # Fill NaNs with nearest value or 0 (cheap trick, acceptable for non-critical edge cells)
        filled_data = np.nan_to_num(data, nan=0)
        curvature = convolve2d(filled_data, kernel, mode='same', boundary='symm')

        # Set curvature to nodata where input was nodata
        curvature[np.isnan(data)] = -9999
        curvature = np.negative(curvature)

        profile.update(dtype='float32', nodata=-9999)

        with rasterio.open(out_tif, 'w', **profile) as dst:
            dst.write(curvature.astype(np.float32), 1)
    except Exception as e:
        print (e)

for d in [9,15,30,45,90]:
    calc_curv_faster('/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/Rasters/TL47/elev.tif', dx=9)

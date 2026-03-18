import numpy as np
import rasterio as rio
import os
from scipy.ndimage import convolve

def calc_all_curvatures(tif):

    out_base = os.path.splitext(tif)[0]

    with rio.open(tif) as src:
        dem = src.read(1).astype(np.float64)
        nodata = src.nodata
        profile = src.profile
        dx = src.transform[0]
        dy = -src.transform[4]

    dem[dem == nodata] = np.nan

    ############################################################
    # First derivatives (Horn 1981 method)
    ############################################################

    kernel_dx = np.array([
        [-1, 0, 1],
        [-2, 0, 2],
        [-1, 0, 1]
    ]) / (8 * dx)

    kernel_dy = np.array([
        [-1, -2, -1],
        [0, 0, 0],
        [1, 2, 1]
    ]) / (8 * dy)

    p = convolve(dem, kernel_dx, mode="nearest")
    q = convolve(dem, kernel_dy, mode="nearest")

    ############################################################
    # Second derivatives
    ############################################################

    kernel_dxx = np.array([
        [0, 0, 0],
        [1, -2, 1],
        [0, 0, 0]
    ]) / (dx**2)

    kernel_dyy = np.array([
        [0, 1, 0],
        [0, -2, 0],
        [0, 1, 0]
    ]) / (dy**2)

    kernel_dxy = np.array([
        [1, 0, -1],
        [0, 0, 0],
        [-1, 0, 1]
    ]) / (4 * dx * dy)

    r = convolve(dem, kernel_dxx, mode="nearest")
    t = convolve(dem, kernel_dyy, mode="nearest")
    s = convolve(dem, kernel_dxy, mode="nearest")

    ############################################################
    # Terrain metrics
    ############################################################

    slope = np.sqrt(p**2 + q**2)

    denom = (p**2 + q**2)

    # Profile curvature
    profile_curv = -(r*p**2 + 2*s*p*q + t*q**2) / (denom * (1 + denom)**1.5)

    # Planform curvature
    plan_curv = (r*q**2 - 2*s*p*q + t*p**2) / (denom**1.5)

    # Mean curvature
    mean_curv = ((1 + q**2)*r - 2*p*q*s + (1 + p**2)*t) / (2*(1 + p**2 + q**2)**1.5)

    # Gaussian curvature
    gauss_curv = (r*t - s**2) / (1 + p**2 + q**2)**2

    # Total curvature (Laplacian)
    total_curv = r + t

    ############################################################
    # Replace NaNs
    ############################################################

    outputs = {
        "slope": slope,
        "profile_curv": profile_curv,
        "plan_curv": plan_curv,
        "mean_curv": mean_curv,
        "gauss_curv": gauss_curv,
        "total_curv": total_curv
    }

    profile.update(dtype="float64", nodata=-9999, count=1)

    for name, arr in outputs.items():

        arr[np.isnan(arr)] = -9999

        out_tif = f"{out_base}_{name}.tif"

        with rio.open(out_tif, "w", **profile) as dst:
            dst.write(arr, 1)


calc_all_curvatures(
    "/Users/evanthaler/Documents/Projects/permafrost/SewardLidar/KGRiver/KG_RiverSites3m.tif"
)

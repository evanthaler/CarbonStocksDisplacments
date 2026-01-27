import numpy as np
import rasterio
from scipy.ndimage import gaussian_laplace

def gaussian_laplacian_curvature(
    dem,
    pixel_size,
    sigma_m,
    nodata=None,
    scale_normalized=True
):
    """
    Gaussian Laplacian curvature.
    Convex = negative, Concave = positive.
    """

    sigma_px = sigma_m / pixel_size

    dem_filled = dem.copy()
    if nodata is not None:
        dem_filled[dem == nodata] = np.nan

    # Gaussian Laplacian
    log = gaussian_laplace(dem_filled, sigma=sigma_px, mode="reflect")
    curvature = -log / (pixel_size ** 2)

    if nodata is not None:
        curvature[np.isnan(dem_filled)] = nodata

    return curvature


def process_gaussian_curvature(
    dem_path,
    output_prefix,
    outpath,
    sigmas_m
):
    with rasterio.open(dem_path) as src:
        dem = src.read(1).astype(float)
        profile = src.profile
        pixel_size = src.res[0]
        nodata = src.nodata

    for sigma in sigmas_m:
        print(f"Gaussian curvature σ = {sigma} m")

        curv = gaussian_laplacian_curvature(
            dem,
            pixel_size,
            sigma_m=sigma,
            nodata=nodata,
            scale_normalized=True
        )

        profile.update(dtype=rasterio.float32, count=1)

        outname = f"{outpath}/{output_prefix}_gauss_sigma{int(sigma)}m_laplacian.tif"
        with rasterio.open(outname, "w", **profile) as dst:
            dst.write(curv.astype(np.float32), 1)

        print(f"Saved {outname}")


if __name__ == "__main__":
    dem_path = '/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/Rasters/combinedSiteRasters/elev.tif'
    output_prefix = "dem_curv"
    outpath = '/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/Rasters/combinedSiteRasters/'
    sigmas_m = [9, 15, 30, 60, 120]

    process_gaussian_curvature(dem_path, output_prefix, outpath,sigmas_m)




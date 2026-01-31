import numpy as np
import rasterio
from scipy.ndimage import generic_filter


# -----------------------------------------------------------
# Core quadratic fit (single output selector)
# -----------------------------------------------------------
def quadratic_curvature(window, dx, nodata, mode):
    """
    window arrives as a FLATTENED vector from generic_filter
    """
    if nodata is not None and np.any(window == nodata):
        return nodata

    z = window.astype(float)
    n = int(np.sqrt(z.size))

    if n * n != z.size:
        return nodata

    half = n // 2

    x, y = np.meshgrid(
        np.arange(-half, half + 1) * dx,
        np.arange(-half, half + 1) * dx
    )

    X = np.column_stack([
        x.ravel()**2,
        y.ravel()**2,
        x.ravel() * y.ravel(),
        x.ravel(),
        y.ravel(),
        np.ones(z.size)
    ])

    # Least squares quadratic fit
    coeffs, _, _, _ = np.linalg.lstsq(X, z, rcond=None)

    # Enforce sign convention: convex = negative
    a, b, c, d, e, _ = -coeffs

    if mode == "laplacian":
        return 2 * a + 2 * b

    elif mode == "mean":
        return a + b

    elif mode == "gaussian":
        return 4 * a * b - c**2

    elif mode == "profile":
        return a * d**2 + b * e**2 + c * d * e

    elif mode == "plan":
        return a * e**2 + b * d**2 - c * d * e

    else:
        raise ValueError("Unknown curvature mode")


# -----------------------------------------------------------
# Polynomial curvature maps
# -----------------------------------------------------------
def polynomial_curvature_map(
    dem,
    pixel_size,
    window_m,
    nodata,
    mode
):
    window_px = int(np.round(window_m / pixel_size))
    if window_px % 2 == 0:
        window_px += 1

    print(f"  Window: {window_px} × {window_px} pixels")

    return generic_filter(
        dem,
        quadratic_curvature,
        size=window_px,
        mode="reflect",
        extra_arguments=(pixel_size, nodata, mode)
    )


# -----------------------------------------------------------
# Driver
# -----------------------------------------------------------
def process_polynomial_curvature(
    dem_path,
    output_prefix,
    outpath,
    windows_m
):
    with rasterio.open(dem_path) as src:
        dem = src.read(1).astype(float)
        profile = src.profile
        pixel_size = src.res[0]
        nodata = src.nodata

    curvature_modes = [
        "laplacian",
        "mean",
        "gaussian",
        "profile",
        "plan"
    ]

    profile.update(dtype=rasterio.float32, count=1)

    for w in windows_m:
        print(f"\nPolynomial curvature — window {w} m")

        for mode in curvature_modes:
            print(f"  Computing {mode}")

            curv = polynomial_curvature_map(
                dem,
                pixel_size,
                window_m=w,
                nodata=nodata,
                mode=mode
            )

            outname = (
                f"{outpath}/"
                f"{output_prefix}_poly_w{int(w)}m_{mode}.tif"
            )

            with rasterio.open(outname, "w", **profile) as dst:
                dst.write(curv.astype(np.float32), 1)

            print(f"  Saved {outname}")


# -----------------------------------------------------------
# Run
# -----------------------------------------------------------
if __name__ == "__main__":
    dem_path = (
'/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/Rasters/combinedSiteRasters/elev.tif'
    )

    outpath = (
'/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/Rasters/combinedSiteRasters/'
    )

    output_prefix = "dem_curv"
    windows_m = [9, 15, 30, 60, 120]

    process_polynomial_curvature(
        dem_path,
        output_prefix,
        outpath,
        windows_m
    )

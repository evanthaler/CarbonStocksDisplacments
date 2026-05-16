import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from shapely.geometry import box
import numpy as np
import rasterio
from rasterio.plot import show as rioshow
import matplotlib.colors as colors
import pandas as pd
from shapely.ops import transform as shp_transform
from pyproj import Transformer
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
plt.rcParams['font.size'] = 11 

def plot_points_in_inset(ax, gdf_points, *, face="deepskyblue", edge="black", s=28, lw=0.6, z=50):
    """Plot GeoDataFrame points on an axes using marker size in display units."""
    if gdf_points.empty:
        return
    xy = np.column_stack((gdf_points.geometry.x, gdf_points.geometry.y))
    ax.scatter(xy[:, 0], xy[:, 1], s=s, c=face, edgecolors=edge, linewidths=lw, zorder=z)

# ============================================================
# Load data
# ============================================================
or_shp = gpd.read_file('/Users/evanthaler/Documents/Projects/OSU/StateShapefiles/tl_2023_us_state/tl_2023_us_state.shp').to_crs('epsg:26910')
or_shp = or_shp[or_shp.STUSPS.isin(['OR', 'WA', 'CA','ID','NV'])]
#/Users/evanthaler/Documents/Projects/OSU/StateShapefiles/tl_2023_us_state/tl_2023_us_state.shp
#/Users/evanthaler/Documents/Projects/OSU/ORgpkg.gpkg
burned = gpd.read_file('/Users/evanthaler/Documents/Projects/OSU/PimontEdits/Data/PimontGeospatialData/BurnedWatershed.shp').to_crs(or_shp.crs)
burned_points = gpd.read_file('/Users/evanthaler/Documents/Projects/OSU/PimontEdits/Data/PimontGeospatialData/Infestation Creek Sampling Locations.shp').to_crs(or_shp.crs)

unburned = gpd.read_file('/Users/evanthaler/Documents/Projects/OSU/PimontEdits/Data/PimontGeospatialData/Unburned_Catchment_Export.shp').to_crs(or_shp.crs)
unburned_points = gpd.read_file('/Users/evanthaler/Documents/Projects/OSU/PimontEdits/Data/PimontGeospatialData/False_Berry_Merged_2024.shp').to_crs(or_shp.crs)

hillshade_path = "/Users/evanthaler/Documents/Projects/OSU/PimontEdits/Data/PimontGeospatialData/USGS_Hillshade_epsg26910.tif"
burnraster = "/Users/evanthaler/Documents/Projects/OSU/PimontEdits/Data/PimontGeospatialData/MTBSBurnSeverityClasses_epsg26910.tif"

watersheds = {'Infestation\nCreek': burned, 'False Berry Creek': unburned}
points = {'Infestation Creek': burned_points, 'False Berry Creek': unburned_points}

# Burn severity legend
classes = [1, 2, 3, 4]
class_colors = ["None","#fee090", "#fc8d59", "#d73027"]
#class_colors = ["#4575b4","#fee090", "#fc8d59", "#d73027"]
class_labels = ["Unburned", "Low", "Moderate", "High"]
cmap = colors.ListedColormap(class_colors)
norm = colors.BoundaryNorm(classes + [5], cmap.N)

# Load rasters
with rasterio.open(hillshade_path) as src:
    hillshade = src.read(1)
    hs_extent = (src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top)

with rasterio.open(burnraster) as src:
    burn_data = src.read(1)
    burn_data = np.ma.masked_where(burn_data == src.nodata, burn_data)
    burn_data = np.ma.masked_where(burn_data == 0, burn_data)
    burn_extent = (src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top)


# ============================================================
# Figure layout with gridspec
# ============================================================
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(
    2, 2,
    width_ratios=[1.25, 1],
    height_ratios=[1, 1],
    wspace=0.18,
    hspace=0.12
)

axA = fig.add_subplot(gs[:, 0])   # big left panel
axB = fig.add_subplot(gs[0, 1])   # top right
axC = fig.add_subplot(gs[1, 1])   # bottom right


# ============================================================
# ------------------ PANEL A (Study Area Map) -----------------
# ============================================================



# Determine zoom
all_bounds = gpd.GeoSeries([gdf.union_all() for gdf in watersheds.values()]).total_bounds
axA.set_xlim(all_bounds[0] - 5000, all_bounds[2] + 5000)
axA.set_ylim(all_bounds[1] - 5000, all_bounds[3] + 5000)
fontprops = fm.FontProperties(size=9)

scalebar = AnchoredSizeBar(
    axA.transData,
    5000,                  # length in meters (1 km)
    '5 km',                # label
    'lower right',         # location
    pad=0.4,
    color='black',
    frameon=True,
    size_vertical=40,      # thickness of bar (in data units)
    fontproperties=fontprops
)
axA.add_artist(scalebar)

axA.text(
    0.02, 0.98, 'a',
    transform=axA.transAxes,
    fontsize=20,
    fontweight='bold',
    va='top',
    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2)
)
# Hillshade first (zorder=0)
with rasterio.open(hillshade_path) as src:
    # Plot raster
    rioshow(src, ax=axA, cmap="gray", zorder=0)

# Burn severity overlay (zorder=1)
axA.imshow(burn_data, extent=burn_extent, cmap=cmap, norm=norm,
           alpha=0.45, zorder=1)

# Watersheds and labels (fixed placement)
for name, gdf in watersheds.items():
    gdf.boundary.plot(ax=axA, color="black", linewidth=4, zorder=5)

    # Place label INSIDE map rather than outside
    minx, miny, maxx, maxy = gdf.geometry.total_bounds
    axA.text(
        minx + (maxx - minx)*0.4,      # 5% from left
        miny + (maxy - miny)*2,      # 2x up
        name,
        fontsize=13, fontweight='bold', color='black',
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'),
        zorder=10
    )
axA.set_xlabel("Easting (m)", labelpad=8)
axA.set_ylabel("Northing (m)", labelpad=8)
#################
# Inset MAP
###############
# 1) Create inset axes (your placement/size)
ax_inset = inset_axes(
    axA,
    width="100%", height="100%",
    bbox_to_anchor=(0.01, 0.15, 0.65, 0.65),  # x, y, w, h in AXES fractions
    bbox_transform=axA.transAxes,
    borderpad=0
)

# 2) Reproject the state layer to 5070 for the inset
states_5070 = or_shp.to_crs("EPSG:5070")

# --- pick which states to draw in inset ---
# Include OR by default so your star/box are clearly in context.
inset_states = states_5070[states_5070.STUSPS.isin(["CA", "NV", "ID", "WA","OR"])].copy()
# If you truly do NOT want OR drawn, use:
# inset_states = states_5070[states_5070.STUSPS.isin(["CA", "NV", "ID"])].copy()

inset_states.plot(ax=ax_inset, color="white", edgecolor="black", linewidth=0.6, zorder=1)

# 3) Set inset extent to Northern CA + western NV + western ID (defined in lon/lat, then projected)
# Tune these numbers if you want tighter/looser framing.
lon_min, lon_max = -124.5, -117.0
lat_min, lat_max = 39.5, 50.5

to5070 = Transformer.from_crs("EPSG:4326", "EPSG:5070", always_xy=True)
x0, y0 = to5070.transform(lon_min, lat_min)
x1, y1 = to5070.transform(lon_max, lat_max)

xmin, xmax = sorted([x0, x1])
ymin, ymax = sorted([y0, y1])

ax_inset.set_xlim(xmin, xmax)
ax_inset.set_ylim(ymin, ymax)

# 4) Reproject sampling points to 5070 and compute a single representative location (star)
all_pts = gpd.GeoDataFrame(
    pd.concat([burned_points, unburned_points], ignore_index=True),
    crs=or_shp.crs
).to_crs("EPSG:5070")

star_x = all_pts.geometry.x.mean()
star_y = all_pts.geometry.y.mean()

ax_inset.scatter(
    star_x, star_y, marker="*", s=100, c="red", edgecolors="gold",
    linewidths=0.8, zorder=10
)

# 5) Transform your zoom box (in 26910) into 5070 and draw it
zoom_box_26910 = box(*all_bounds)
transformer = Transformer.from_crs("EPSG:26910", "EPSG:5070", always_xy=True)
zoom_box_5070 = shp_transform(transformer.transform, zoom_box_26910)

gpd.GeoSeries([zoom_box_5070], crs="EPSG:5070").plot(
    ax=ax_inset, facecolor="none", edgecolor="red", linewidth=2.0, zorder=8
)

# 6) Add state name labels (use representative_point() so text is inside polygon)
name_map = {"CA": "California","WA":"Washington", "OR": "Oregon"}
# Optional nudges if labels overlap; values are in EPSG:5070 meters
offsets = {
    "CA": (-80000, 30000),
    "OR": (-20000,  50000),
    "WA": (-20000,  70000),
}

for _, row in inset_states.iterrows():
    abbr = row["STUSPS"]
    if abbr not in name_map:
        continue
    rp = row.geometry.representative_point()
    dx, dy = offsets.get(abbr, (0, 0))
    ax_inset.text(
        rp.x + dx, rp.y + dy, name_map[abbr],
        fontsize=10, ha="center", va="center",
        color="black",
        bbox=dict(facecolor="white", alpha=0.6, edgecolor="none", pad=1.2),
        zorder=6
    )

# 7) Label Pacific Ocean (place near coast inside the inset extent)
ax_inset.text(
    xmin + 0.10*(xmax - xmin),
    ymin + 0.65*(ymax - ymin),
    "Pacific Ocean",
    fontsize=10, fontstyle="italic",
    color="steelblue",
    rotation=90,
    ha="center", va="center",
    alpha=0.9,
    zorder=4
)

# 8) Clean inset axes
ax_inset.set_xticks([])
ax_inset.set_yticks([])
ax_inset.set_in_layout(False)
axA.tick_params(axis='x', rotation=30) 
axA.tick_params(axis='y', rotation=90) 
axA.ticklabel_format(axis='y', style='plain')

# ============================================================
# ------------- PANELS B and C (Zoomed Watersheds) -----------
# ============================================================

for ax, (ws_name, ws_gdf), (_, pts_gdf), letter in zip(
        [axB, axC],
        watersheds.items(),
        points.items(),
        ["b", "c"]
):

    # Zoom to watershed
    minx, miny, maxx, maxy = ws_gdf.geometry.total_bounds
    ax.set_xlim(minx - 200, maxx + 200)
    ax.set_ylim(miny - 200, maxy + 200)
    
    # minx, miny, maxx, maxy = ws_gdf.geometry.total_bounds
    # ax.set_xlim(minx , maxx )
    # ax.set_ylim(miny , maxy )


    fontprops = fm.FontProperties(size=9)

    scalebar = AnchoredSizeBar(
        ax.transData,
        1000,                  # length in meters (1 km)
        '1 km',                # label
        'lower right',         # location
        pad=0.4,
        color='black',
        frameon=True,
        size_vertical=40,      # thickness of bar (in data units)
        fontproperties=fontprops
    )

    ax.add_artist(scalebar)
    with rasterio.open(hillshade_path) as src:
    # Plot raster
        rioshow(src, ax=ax, cmap="gray", zorder=0)

    # Burn severity
    ax.imshow(burn_data, extent=burn_extent, cmap=cmap, norm=norm,
              alpha=0.45, zorder=1)

    # Watershed boundary
    ws_gdf.boundary.plot(ax=ax, color="black", linewidth=4, zorder=5)

    # Sampling points
    pts_gdf.plot(ax=ax, marker="o", color="None", edgecolor="black",
                 markersize=16,alpha=0.9, zorder=10)

    # Labels and ticks
    # ax.set_xlabel("Easting")
    # ax.set_ylabel("Northing")
    ax.tick_params(axis='x', rotation=50)
    ax.tick_params(axis='y', rotation=50)
    ax.ticklabel_format(axis='y', style='plain')

    # Subplot label
    ax.text(
        0.02, 0.98, letter,
        transform=ax.transAxes,
        fontsize=20,
        fontweight='bold',
        va='top',
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2)
    )


# ============================================================
# Shared colorbar
# ============================================================
cbar = fig.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap=cmap),
    ax=[axA, axB, axC],
    shrink=0.7,
    pad=0.02,
    location='right'
)
cbar.set_ticks(classes)
cbar.set_ticklabels(class_labels)
cbar.set_label("Burn Severity", fontsize=11)
#plt.tight_layout()
plt.savefig('/Users/evanthaler/Documents/Projects/OSU/PimontEdits/Figures/StudyAreaMap.jpg',dpi=500)
plt.show()

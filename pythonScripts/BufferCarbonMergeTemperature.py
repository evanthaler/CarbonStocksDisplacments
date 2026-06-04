import geopandas as gpd
import matplotlib.pyplot as plt

buffer_distance = 50

carbongdf = gpd.read_file('FinalCleanedFiles/wDisplacement/TL47TL27StocksDisplacement_wcovars.gpkg')
probegdf = gpd.read_file('FinalCleanedFiles/TemperatureProbes/tl27tl47annualtemperaturestats.gpkg')

carb_buffered = carbongdf.copy()
carb_buffered["geometry"] = carb_buffered.buffer(buffer_distance)

candidates = gpd.sjoin(
    probegdf,
    carb_buffered,
    #carb_buffered[["carbon_id", "slope", "geometry"]],
    how="inner",
    predicate="within"
)

disp_summary = (
    candidates
    .groupby("SampleLocationName", as_index=False)
    .agg(unfozendays=('days_unfrozen_below_50cm', "median"))
)
print(disp_summary)

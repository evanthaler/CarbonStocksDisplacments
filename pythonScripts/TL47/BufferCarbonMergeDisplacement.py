import geopandas as gpd
import matplotlib.pyplot as plt
from scipy.stats import linregress
plt.rcParams.update({'font.size': 14})

buffer_distance = 50
carbon_depth = 50
carbongdf = gpd.read_file(f'/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/FinalCleanedFiles/gpkgs/TL47_TotalSoilStocksSlopes_PerCore_{carbon_depth}_slope.gpkg')
dispgdf = gpd.read_file('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/displacement/gpkgs/TL47AllDisplacements_wSlope.gpkg')

col_to_plot =f'Total_C_stock_kg_m2_0to{carbon_depth}cm'
carb_buffered = carbongdf.copy()
carb_buffered["geometry"] = carb_buffered.buffer(buffer_distance)

candidates = gpd.sjoin(
    dispgdf,
    carb_buffered,
    #carb_buffered[["carbon_id", "slope", "geometry"]],
    how="inner",
    predicate="within"
)
candidates = candidates[
    (candidates["slope_left"] - candidates["slope_right"]).abs() <= 2
]
disp_summary = (
    candidates
    .groupby("SampleLocationName_right", as_index=False)
    .agg(mean_disp_rate=("dhorizontal", "mean"),
    mean_temp=("pfTemp","mean"))
)
carbon_unique = (
    carbongdf
    .drop_duplicates(subset="SampleLocationName")
)
carbon_out = carbon_unique.merge(
    disp_summary,
    left_on="SampleLocationName",
    right_on="SampleLocationName_right",
    how="left"
)
carbon_out = carbon_out.drop(columns=['c_frac', 'cstock',
       'geometry', 'SampleLocationName_right',],axis=1)
carbon_out = carbon_out.dropna()

#carbon_out = carbon_out[carbon_out.SampleLocationName!='TL47 GullyPit 1']
#carbon_out.to_csv('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/FinalCleanedFiles/wDisplacement/TL47StocksDisplacement.csv')

######################################
#Carbon vs displacement, topography###
linstats=linregress(carbon_out['mean_disp_rate'],carbon_out[col_to_plot])
print('carb, disp',linstats.rvalue**2,linstats.pvalue)
plt.figure(figsize=(5,5))
plt.ylabel('Soil organic carbon stock (kg m$^{-2}$)')
plt.xlabel('Horizontal displacement rate (m yr$^{-1}$)')
plt.plot(carbon_out['mean_disp_rate'],carbon_out[col_to_plot],'ok')
#plt.scatter(carbon_out['mean_disp_rate'],carbon_out[col_to_plot],c=carbon_out.mean_temp,s=100)
plt.plot(carbon_out['mean_disp_rate'],linstats.slope*(carbon_out['mean_disp_rate'])+linstats.intercept,'-k')
plt.savefig('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/figs/TL47_SOCStock_Displacement.jpg',dpi=300)
plt.show()

linstats_curv=linregress(carbon_out['curv'],carbon_out[col_to_plot])
print('carb, curv',linstats_curv.rvalue**2,linstats_curv.pvalue)
plt.figure(figsize=(5,5))
plt.ylabel('Soil organic carbon stock (kg m$^{-2}$)')
plt.xlabel('Topographic curvature')
plt.plot(carbon_out['curv'],carbon_out[col_to_plot],'ok')
plt.plot(carbon_out['curv'],linstats_curv.slope*(carbon_out['curv'])+linstats_curv.intercept,'-k')
plt.show()


linstats_slope=linregress(carbon_out['slope'],carbon_out[col_to_plot])
print('carb, slope',linstats_slope.rvalue**2,linstats_slope.pvalue)
plt.figure(figsize=(5,5))
plt.ylabel('Soil organic carbon stock (kg m$^{-2}$)')
plt.xlabel('Topographic slope')
plt.plot(carbon_out['slope'],carbon_out[col_to_plot],'ok')
plt.show()


# #############################
# #Displacement vs topography##
# #############################
# displinstats=linregress(carbon_out['curv'],carbon_out['mean_disp_rate'])
# print('disp, curv',displinstats.rvalue**2,displinstats.pvalue)
# plt.figure(figsize=(5,5))
# plt.ylabel('Horizontal displacement rate (m yr$^{-1}$)')
# plt.xlabel('Topographic curvature')
# plt.plot(carbon_out['curv'],carbon_out['mean_disp_rate'],'ok')
# plt.plot(carbon_out['curv'],displinstats.slope*(carbon_out['curv'])+displinstats.intercept,'-k')
# plt.show()

# displinstats_slope=linregress(carbon_out['slope'],carbon_out['mean_disp_rate'])
# print('disp, slope',displinstats_slope.rvalue**2)
# plt.figure(figsize=(5,5))
# plt.ylabel('Horizontal displacement rate (m yr$^{-1}$)')
# plt.xlabel('Topographic slope')
# plt.plot(carbon_out['slope'],carbon_out['mean_disp_rate'],'ok')
# plt.plot(carbon_out['slope'],displinstats_slope.slope*(carbon_out['slope'])+displinstats_slope.intercept,'-k')
# plt.show()
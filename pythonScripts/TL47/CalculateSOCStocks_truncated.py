import pandas as pd
import numpy as np
import geopandas as gpd

def interpolate_bulk_density(df,
                             sample_col='SampleLocationName',
                             top_col='Top_cm',
                             bottom_col='Bottom_cm',
                             bd_col='BulkDensity_gcm3'):
    """
    Interpolates bulk density values down each soil profile using depth-based linear interpolation.
    Only interpolates within profiles, no reindexing required (avoids duplicate depth issues).
    """
    out = df.copy()
    
    # depth midpoint for ordering
    out['depth_mid_cm'] = (out[top_col] + out[bottom_col]) / 2

    # sort before interpolation to ensure correct ordering
    out = out.sort_values([sample_col, 'depth_mid_cm'])

    # interpolate within each core
    out[bd_col] = (
        out.groupby(sample_col)[bd_col]
           .transform(lambda s: s.interpolate(method='linear').bfill().ffill())
    )

    return out

def truncate_to_depth(df,
                      cutoff_depth_cm=50,
                      top_col='Top_cm',
                      bottom_col='Bottom_cm',
                      c_stock_col='C_stock_kg_m2',
                      n_stock_col='N_stock_kg_m2'):
    """
    Truncates soil layer stocks to a maximum depth by proportionally scaling layers
    that cross the cutoff depth.
    """
    out = df.copy()

    # Compute layer thickness actually within cutoff depth
    out['effective_thickness_cm'] = (
        out[[top_col, bottom_col]]
        .apply(lambda r: max(0, min(r[1], cutoff_depth_cm) - r[0]), axis=1)
    )

    # Avoid division errors
    out['thickness_ratio'] = out['effective_thickness_cm'] / out['layerthickness_cm']

    # Scale stocks
    out['C_stock_trunc_kg_m2'] = out[c_stock_col] * out['thickness_ratio']
    out['N_stock_trunc_kg_m2'] = out[n_stock_col] * out['thickness_ratio']

    # Keep only layers contributing to the profile
    out = out[out['effective_thickness_cm'] > 0].copy()

    return out

def add_carbon_nitrogen_stocks(df,
                               bd_col='BulkDensity_gcm3',
                               c_col='%C',
                               n_col='%N',
                               thickness_col='layerthickness_cm'):
    """
    Adds carbon and nitrogen stock (kg/m^2) columns to the dataframe.
    Stock (kg/m^2) = BulkDensity(g/cm3) * (percentage/100) * thickness(cm) * 10
    """
    out = df.copy()
    out['C_frac'] = out[c_col] / 100.0
    out['N_frac'] = out[n_col] / 100.0

    out['C_stock_kg_m2'] = out[bd_col] * out['C_frac'] * out[thickness_col] * 10
    out['N_stock_kg_m2'] = out[bd_col] * out['N_frac'] * out[thickness_col] * 10
    
    return out



if __name__ == "__main__":
    df = pd.read_csv('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/FinalCleanedFiles/UltimateSheet/TL47_UltimateSheet.csv')

    # 1) Interpolate BD
    df_interp = interpolate_bulk_density(df)

    # 2) Compute layer-level C/N stocks
    df_stocks = add_carbon_nitrogen_stocks(df_interp)

    for c in [40,50,60,100]:
        cutoff_depth_cm = c
        # 3) Truncate to given depth
        df_trunc = truncate_to_depth(df_stocks, cutoff_depth_cm=cutoff_depth_cm)

        C_total_col = f'Total_C_stock_kg_m2_0to{cutoff_depth_cm}cm'
        N_total_col = f'Total_N_stock_kg_m2_0to{cutoff_depth_cm}cm'
        # 4) Sum truncated stocks per core
        df_profile_totals = (
            df_trunc
            .groupby('SampleLocationName', as_index=False)
            .agg(
                **{
                    C_total_col: ('C_stock_trunc_kg_m2', 'sum'),
                    N_total_col: ('N_stock_trunc_kg_m2', 'sum')
                }
            )
        )

        df_core_meta = (
            df
            .sort_values('Top_cm')        # optional: ensure surface row kept
            .drop_duplicates(subset='SampleLocationName')
        )

        merge_df = df_profile_totals.merge(
            df_core_meta,
            on='SampleLocationName',
            how='left'
        )
        merge_df = merge_df.drop(columns=['SampleName', 'Top_cm', 'Bottom_cm',
        'layerthickness_cm', 'BD_measurement_cm', 'BulkDensity_gcm3', '%C', '%N',
        'C/N', 'stdev%C', 'stdev%N', 'stdevC/N'],axis=1)
        
        merge_df=merge_df.dropna()
        merge_df=merge_df.drop_duplicates()
        # Save
        #out_gdf = gpd.GeoDataFrame(merge_df,geometry=gpd.points_from_xy(x=merge_df.X,y=merge_df.Y),crs='epsg:32603')
        #out_gdf.to_file(f'/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/FinalCleanedFiles/gpkgs/TL47_TotalSoilStocksSlopes_PerCore_{cutoff_depth_cm}.gpkg',driver='GPKG')
        #df_trunc.to_csv('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/FinalCleanedCSVs/TL47TL27_SoilCarbonNitrogen_WithLayerStocksTEST.csv', index=False)
        merge_df.to_csv(f'/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/FinalCleanedFiles/TL47_TotalSoilStocks_PerCore_{cutoff_depth_cm}cm.csv', index=False)
        



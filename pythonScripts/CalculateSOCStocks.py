import pandas as pd
import numpy as np

# === User-defined settings ===
MAX_DEPTH_CM = 50  # only calculate stocks from 0 to this depth

# === Load data ===
data = pd.read_csv("/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/KR/KR_CoreSamples_utmCoords_Corvars.csv")

# === Output containers ===
site_totals = []
horizon_records = []

# === Loop through each profile ===
for site in data['SampleLocationName'].dropna().unique():
    df_site = data[data['SampleLocationName'] == site].copy()

    # 1. Extract bulk density measurements
    bd_df = df_site[df_site['BD_measurement_cm'].notna() & df_site['BulkDensity_gcm3'].notna()]
    bd_depths = bd_df['BD_measurement_cm'].values
    bd_values = bd_df['BulkDensity_gcm3'].values

    # 2. Filter layers with C data
    c_df = df_site[df_site['Top_cm'].notna() & df_site['Bottom_cm'].notna()].copy()

    # 3. Exclude or truncate layers beyond MAX_DEPTH
    c_df = c_df[c_df['Top_cm'] < MAX_DEPTH_CM].copy()
    c_df['Bottom_cm'] = c_df['Bottom_cm'].clip(upper=MAX_DEPTH_CM)
    c_df['thickness_cm'] = c_df['Bottom_cm'] - c_df['Top_cm']

    # 4. Calculate carbon fraction and midpoint
    c_df['carbon_frac'] = c_df['%C'] / 100
    c_df['midpoint_cm'] = (c_df['Top_cm'] + c_df['Bottom_cm']) / 2

    # 5. Assign BD via interpolation or fallback
    if len(bd_values) >= 2:
        c_df['BulkDensity_gcm3'] = np.interp(c_df['midpoint_cm'], bd_depths, bd_values)
    elif len(bd_values) == 1:
        c_df['BulkDensity_gcm3'] = bd_values[0]
    else:
        c_df['BulkDensity_gcm3'] = df_site['BulkDensity_gcm3'].mean(skipna=True)

    # 6. Compute horizon-level C stock
    c_df['C_stock_kg_m2'] = (
        c_df['BulkDensity_gcm3'] *
        c_df['thickness_cm'] *
        c_df['carbon_frac']
    )

    # 7. Record each horizon
    for _, row in c_df.iterrows():
        horizon_records.append({
            'SampleLocationName': site,
            'Top_cm': row['Top_cm'],
            'Bottom_cm': row['Bottom_cm'],
            'Thickness_cm': row['thickness_cm'],
            '%C': row['%C'],
            'BulkDensity_gcm3': row['BulkDensity_gcm3'],
            'C_stock_kg_m2': row['C_stock_kg_m2']
        })

    # 8. Total carbon for profile
    total_carbon = c_df['C_stock_kg_m2'].sum()
    site_totals.append({'SampleLocationName': site, 'C_stock_kg_m2': total_carbon})

# === Save outputs ===
df_total = pd.DataFrame(site_totals)
df_horizons = pd.DataFrame(horizon_records)

df_total.to_csv("/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/KR/KR_CoreSamples_stocks.csv", index=False)
df_horizons.to_csv("/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/FinalCleanedFiles/KR/KR_CoreSamples_horizonStocks.csv", index=False)

print(f"Finished! Profiles truncated to {MAX_DEPTH_CM} cm and outputs saved.")


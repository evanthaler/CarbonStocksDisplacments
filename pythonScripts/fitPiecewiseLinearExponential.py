import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
plt.rcParams.update({'font.size': 14})


def linear_exp(x, x0, y0, m, k):
    return np.where(
        x <= x0,
        m * (x - x0) + y0,
        y0 * np.exp(-k * (x - x0))
    )

df = pd.read_csv('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/FinalCleanedFiles/wDisplacement/TL47TL27StocksDisplacement.csv')
x = df.mean_disp_rate
y = df.Total_C_stock_kg_m2_0to50cm


x0_init = x[np.argmax(y)]# pick the y peak location for x
y0_init = y.max()

p0 = [
    x0_init,      # x0 (peak location)
    y0_init,      # y0 (peak value)
    200,          # m (slope; adjust if needed)
    50            # k (decay rate; adjust if needed)
]

popt, pcov = curve_fit(
    linear_exp, x, y, p0=p0, maxfev=10000
)
x0, y0, m, k = popt
#Get peak standard error
x0_se = np.sqrt(pcov[0, 0])          # standard error
x0_ci = 1.96 * x0_se   

x_fit = np.linspace(x.min(), x.max(), 400)
y_fit = linear_exp(x_fit, *popt)



plt.figure()
plt.plot(df.mean_disp_rate,df.Total_C_stock_kg_m2_0to50cm,'ok')
#plt.scatter(df.mean_disp_rate,df.Total_C_stock_kg_m2_0to50cm,c=df.curv,s=500)
plt.colorbar(label='Slope (degrees)')
plt.plot(x_fit, y_fit, color='k', lw=2)
plt.axvspan(
    x0 - x0_ci, x0 + x0_ci,
    color="gray", alpha=0.2,
    label="95% CI (peak)"
)
plt.ylabel('Soil organic carbon stock (kg m$^{-2}$)')
plt.xlabel('Horizontal displacement rate (m yr$^{-1}$)')
#plt.savefig('/Users/evanthaler/Documents/Projects/permafrost/permafrostCarbon/figs/SOCStock_DisplacementCombinedSites.jpg',dpi=300)
plt.show()

print("Fitted model:")
print(f"For x ≤ {x0:.4f}:")
print(f"  y = {m:.2f} (x - {x0:.4f}) + {y0:.2f}")
print(f"For x > {x0:.4f}:")
print(f"  y = {y0:.2f} · exp[-{k:.2f} (x - {x0:.4f})]")
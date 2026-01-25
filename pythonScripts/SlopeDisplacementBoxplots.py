import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


disp_df = pd.read_csv('/Users/evanthaler/Documents/GitHub/CarbonStocksDisplacments/disp2use/TL27TL47Displacments_PFTopo.csv')


tl27_pfdisp = disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "1"),'dhorizontal'].dropna()
tl27_nonpfdisp = disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "0"),'dhorizontal'].dropna()
tl47_pfdisp = disp_df.loc[(disp_df['Site'] == 'TL47') & (disp_df['pf'] == "1"),'dhorizontal'].dropna()
tl47_nonpfdisp = disp_df.loc[(disp_df['Site'] == 'TL47') & (disp_df['pf'] == "0"),'dhorizontal'].dropna()

tl27_pfslope = disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "1"),'slope'].dropna()
tl27_nonpfslope = disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "0"),'slope'].dropna()
tl47_pfslope = disp_df.loc[(disp_df['Site'] == 'TL47') & (disp_df['pf'] == "1"),'slope'].dropna()
tl47_nonpfslope = disp_df.loc[(disp_df['Site'] == 'TL47') & (disp_df['pf'] == "0"),'slope'].dropna()

# pf27kde = gaussian_kde(disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "1"),'dhorizontal'].dropna())
# nonpf27kde=gaussian_kde(disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "0"),'dhorizontal'].dropna())
# pf27slopekde = gaussian_kde(disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "1"),'slope'].dropna())
# nonpf27slopekde = gaussian_kde(disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "0"),'slope'].dropna())

# pf47kde = gaussian_kde(disp_df.loc[(disp_df['Site'] == 'TL47') & (disp_df['pf'] == "0"),'dhorizontal'].dropna())
# nonpf47kde=gaussian_kde(disp_df.loc[(disp_df['Site'] == 'TL47') & (disp_df['pf'] == "1"),'dhorizontal'].dropna())
# pf47slopekde = gaussian_kde(disp_df.loc[(disp_df['Site'] == 'TL47') & (disp_df['pf'] == "1"),'slope'].dropna())
# nonpf47slopekde = gaussian_kde(disp_df.loc[(disp_df['Site'] == 'TL47') & (disp_df['pf'] == "0"),'slope'].dropna())

disp27len = np.linspace(0,0.3,100)

slope27len = np.linspace(
    min(disp_df.loc[(disp_df['Site'] == 'TL27') ,'slope']),
    max(disp_df.loc[(disp_df['Site'] == 'TL27') ,'slope']),
    100)


disp47len = np.linspace(0,0.3,100)

slope47len = np.linspace(
    min(disp_df.loc[(disp_df['Site'] == 'TL47') ,'slope']),
    max(disp_df.loc[(disp_df['Site'] == 'TL47') ,'slope']),
    100)

fig,ax = plt.subplots(figsize=(5,5))
#ax.boxplot([tl27_pfdisp,tl27_nonpfdisp,tl47_pfdisp,tl47_nonpfdisp],labels =['TL27 PF','TL27 NonPF','TL47 PF','TL47 NonPF'],showfliers=False )
ax.boxplot([tl27_pfslope,tl27_nonpfslope,tl47_pfslope,tl47_nonpfslope],labels =['TL27 PF','TL27 NonPF','TL47 PF','TL47 NonPF'],showfliers=False )

plt.show()
#ax.boxplot(disp_df.loc[(disp_df['Site'] == 'TL27') & (disp_df['pf'] == "1"),'dhorizontal'].dropna())


plt.show()

# ax.plot(disp27len,pf27kde(disp27len),'--b',label ='TL27 PF')
# ax.plot(disp27len,nonpf27kde(disp27len),'-b',label ='TL27 Non PF')

# ax.plot(disp47len,pf47kde(disp47len),'--k',label ='TL47 PF')
# ax.plot(disp47len,nonpf47kde(disp47len),'-k',label ='TL47 Non PF')
# ax.set_xlabel('Horizontal displacement (m yr$^{-1}$)')
# ax.set_ylabel('Probability distribution')
# ax.set_ylim(0); ax.set_xlim(0)
# ax.legend()
# plt.show()

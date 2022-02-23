#%% Plot pop exposure comparing hazards. Currently my version of Fig C in the paper.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['axes.linewidth'] = .4
matplotlib.rcParams['xtick.major.width'] = .4
matplotlib.rcParams['ytick.major.width'] = .4
# matplotlib.rcParams['text.color'] = 0
# matplotlib.rcParams['axes.labelcolor'] = 0

import seaborn as sns
# sns.set_style("ticks", {"axes.facecolor": ".9"})
sns.set_style("ticks")
sns.set_style({"xtick.direction": "in"})
sns.set_context("paper")
# sns.set_theme(style="whitegrid")
from matplotlib import gridspec

#%% Load
#region preproces
print('Read master exposure')
# data = pd.read_csv('csv/MASTER_exposure.csv', index_col=None).drop('Unnamed: 0', axis=1)
Data = pd.read_csv('csv/MASTER_exposure_v3.csv', index_col=None).drop('Unnamed: 0', axis=1)
# Data = Data[Data.month==0].drop('month', axis=1)

# Prepare roads
Roads = pd.read_csv('csv/MASTER_roads_v3.csv', index_col=None).drop('Unnamed: 0', axis=1)
LC = pd.read_csv('csv/MASTER_roads_v3_LC.csv', index_col=None).drop('Unnamed: 0', axis=1)
Roads = Roads.loc[Roads['hazard']!='LC'].append(LC)
Roads = Roads[Roads.month==0].drop('month', axis=1)
Roads['Motorway'] = Roads.filter(regex='length_motorway*').sum(axis=1)
Roads['Arterial'] = Roads.filter(regex='length_trunk*|length_primary*').sum(axis=1)
Roads['Collector'] = Roads.filter(regex='length_secondary*|length_tertiary*').sum(axis=1)
Roads['Local'] = Roads.filter(regex='length_unclassified*|length_residential*|length_service*|length_living_street*|length_road*|length_unknown*').sum(axis=1)
Roads = Roads.drop(Roads.filter(regex='length*').columns, axis=1)

# Read volcano index with country names
volcDB = pd.read_csv('csv/volcCoordinates.csv')
volcDB = volcDB.rename({'name': 'volcano'}, axis=1)

db = volcDB.copy()
data = Data.copy()
roads = Roads.copy()
data = data.merge(roads, how='outer').merge(db, how='outer')

data = data.replace('Philippines and SE Asia', 'Philippines')
data = data.sort_values(by=['country', 'region', 'volcano'])
data = data.set_index(['hazard', 'VEI', 'prob'])

# Format string
data['volcano'] = data['volcano'].str.replace('_',' ')
data['volcano'] = data['volcano'].str.replace('Leroboleng','Lereboleng')
data['region'] = data['region'].str.replace('Sumatera','Sumatra')
data['region'] = data['region'].str.replace('The Philippines','Philippines')
data['region'] = data['region'].str.replace('Maluku/Halmahera','Halmahera/Banda Sea')
data = data.sort_values(by=['country', 'region', 'volcano'])

# Sum roads and buildings numbers for exposure analysis and drop RNDS and cost
data['length_roads'] = data[['Motorway','Arterial','Collector','Local']].sum(axis=1)
data['n_buildings'] = data[['buildingsW', 'buildingsS']].sum(axis=1)
data = data.drop(['Motorway','Arterial','Collector','Local','buildingsW', 'buildingsS','RNDS', 'buildingsLoss'], axis=1)
data.loc[data.pop_count==0, 'pop_count'] = np.nan
data.pop_count = np.log10(data.pop_count)
#%%
f, ax = plt.subplots(ncols=3, nrows=4, figsize=(9,8),sharey='row')
# sns.despine(bottom=True, left=True)
# sns.despine()


def pltMe(DATA,ax):
    sns.swarmplot(ax=ax, x="pop_count", y="region", hue="region", data=DATA, dodge=True, alpha=.9, zorder=1, size=7, palette="pastel")
    sns.pointplot(ax=ax, x="pop_count", y="region", hue="region", data=DATA, estimator=np.median,dodge=.532, join=False, palette="deep", markers="x", scale=1, ci=None,errwidth=.57)
    ax.get_legend().remove()


data = data.set_index('mass',append=True)
pltMe(data.loc['Tephra',3,50,1].reset_index(),ax[0,0])
pltMe(data.loc['Tephra',4,50,1].reset_index(),ax[0,1])
pltMe(data.loc['Tephra',5,50,1].reset_index(),ax[0,2])
data = data.reset_index().set_index(['hazard', 'VEI', 'prob'])
pltMe(data.loc['PDC',3,50,slice(None)].reset_index(),ax[1,0])
pltMe(data.loc['PDC',4,50,slice(None)].reset_index(),ax[1,1])
pltMe(data.loc['PDC',5,50,slice(None)].reset_index(),ax[1,2])
pltMe(data.loc['LC',3,50,slice(None)].reset_index(),ax[2,0])
pltMe(data.loc['LC',4,50,slice(None)].reset_index(),ax[2,1])
pltMe(data.loc['LC',5,50,slice(None)].reset_index(),ax[2,2])
data = data.reset_index().set_index(['hazard', 'volume','buffer', 'prob'])
pltMe(data.loc['BAF',4.5e5,990,slice(None)].reset_index(),ax[3,0])
pltMe(data.loc['BAF',9.8e6,990,slice(None)].reset_index(),ax[3,1])
# pltMe(data.loc['BAF',9.8e6,990,slice(None)].reset_index(),ax[3,2])
# sns.stripplot(ax=ax[0,0], x="pop_count", y="region", hue="VEI", data=tephra3, dodge=True, alpha=.25, zorder=1, size=5, palette="pastel")
# sns.pointplot(ax=ax[0,0], x="pop_count", y="region", hue="VEI", data=tephra3, dodge=.532, join=False, palette="deep", markers="x", scale=1, ci='sd',errwidth=.57)
# tephra4 = data.loc['Tephra',4,50].reset_index()
# sns.stripplot(ax=ax[0,1], x="pop_count", y="region", hue="VEI", data=tephra4, dodge=True, alpha=.25, zorder=1, size=5, palette="pastel")
# sns.pointplot(ax=ax[0,1], x="pop_count", y="region", hue="VEI", data=tephra4, dodge=.532, join=False, palette="deep", markers="x", scale=1, ci='sd',errwidth=.57)
# tephra5 = data.loc['Tephra',5,50].reset_index()
# sns.stripplot(ax=ax[0,2], x="pop_count", y="region", hue="VEI", data=tephra5, dodge=True, alpha=.25, zorder=1, size=5, palette="pastel")
# sns.pointplot(ax=ax[0,2], x="pop_count", y="region", hue="VEI", data=tephra5, dodge=.532, join=False, palette="deep", markers="x", scale=1, ci='sd',errwidth=.57)

# PDC = data.loc['PDC',slice(None),50].reset_index()
# sns.stripplot(ax=ax[1,0], x="pop_count", y="region", hue="VEI", data=PDC, dodge=True, alpha=.25, zorder=1, size=5, palette="pastel")
# sns.pointplot(ax=ax[1,0], x="pop_count", y="region", hue="VEI", data=PDC, dodge=.532, join=False, palette="deep", markers="x", scale=1, ci='sd',errwidth=.57)

# VEI3 = data.loc[slice(None),3,50].reset_index()
# sns.stripplot(ax=ax[0,0], x="pop_count", y="region", hue="hazard", data=tephra, dodge=True, alpha=.25, zorder=1, size=5, palette="pastel")
# sns.pointplot(ax=ax[0,0], x="pop_count", y="region", hue="hazard", data=tephra, dodge=.532, join=False, palette="deep", markers="x", scale=1, ci='sd',errwidth=.57)

# PDC = data.loc['PDC',slice(None),50].reset_index()
# sns.stripplot(ax=ax[1,0], x="pop_count", y="region", hue="VEI", data=PDC, dodge=True, alpha=.25, zorder=1, size=5, palette="pastel")
# sns.pointplot(ax=ax[1,0], x="pop_count", y="region", hue="VEI", data=PDC, dodge=.532, join=False, palette="deep", markers="x", scale=1, ci='sd',errwidth=.57)

YLAB = ['Tephra ≥ 1kg/m$^2$', 'Column collapse', 'Large clasts', 'Block-and-ash flow']
for iy in range(0,4):
    if iy<3:
        ttl = ['VEI 3', 'VEI 4', 'VEI 5']
    else:
        ttl = ['4.5e$^5$ m$^3$', '9.8e$^6$ m$^3$', '']
    
    for ix in range(0,3):
        ax[iy,ix].set(title=ttl[ix])
        if ix==0:
            ax[iy,ix].set(ylabel=YLAB[iy])
        else:
            ax[iy,ix].set(ylabel=' ')
                          
        if (ix==0 and iy==3) or (ix==1 and iy==3) or (ix==2 and iy==2):
            ax[iy,ix].set(xlabel='Log$_{10}$ population')
        else:
            ax[iy,ix].set(xlabel=f'')
        

f.delaxes(ax[3,2])
f.subplots_adjust(wspace=.1,hspace=.3)

plt.savefig(f"Paper/Figures/populationHazardRegion.pdf",bbox_inches='tight', dpi=300)
plt.savefig(f"Paper/Figures/populationHazardRegion.png",bbox_inches='tight', dpi=300)

# %%

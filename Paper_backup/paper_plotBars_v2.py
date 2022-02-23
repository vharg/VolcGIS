#%%
# Plots exposure as vertical bars (look at _v2 for the first paper)
#v2 updated to plot buildings number and road length instead of RNDS and cost
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib
import matplotlib.patches as mpatches
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['axes.linewidth'] = .3
matplotlib.rcParams['xtick.major.width'] = .3
matplotlib.rcParams['ytick.major.width'] = .3
# matplotlib.rcParams['text.color'] = 0
# matplotlib.rcParams['axes.labelcolor'] = 0

import seaborn as sns
sns.set_style("ticks", {"axes.facecolor": ".9"})
sns.set_context("paper")
from matplotlib import gridspec
import numpy as np

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
Roads = Roads[Roads.month==0]#.drop('month', axis=1)
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

# Load monthly
DataM = pd.read_csv('csv/MASTER_exposure.csv', index_col=None).drop('Unnamed: 0', axis=1)
dataM = DataM.merge(db)
dataM = dataM.set_index(['hazard', 'VEI', 'prob'])
dataM['volcano'] = dataM['volcano'].str.replace('_',' ')
dataM['volcano'] = dataM['volcano'].str.replace('Leroboleng','Lereboleng')
dataM['region'] = dataM['region'].str.replace('Sumatera','Sumatra')
dataM['region'] = dataM['region'].str.replace('The Philippines','Philippines')
dataM['region'] = dataM['region'].str.replace('Maluku/Halmahera','Halmahera/Banda Sea')
dataM = dataM.sort_values(by=['country', 'region', 'volcano'])

# Regions
regions = data.region.unique()
regionsRatio = [data[data.region==reg].volcano.unique().shape[0] for reg in regions]

# Fig size
figW = 9
figType = 'png'
outPth = 'Paper/Figures/'
# Colormaps
reverseCM = True
CMAP = [['#bdd7e7','#6baed6','#2171b5'], ['#cbc9e2','#9e9ac8','#6a51a3']]
CMAPb = [['#fdbe85','#fd8d3c','#e6550d'], ['#78c679','#31a354','#006837']]
CMAPm = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
CMAPi = [['#1d91c0','#225ea8','#081d58'], ['#f768a1','#ae017e','#7a0177']]
if reverseCM:
    CMAP = [i[::-1] for i in CMAP]
    CMAPi = [i[::-1] for i in CMAPi]

# Options
plotLog = True # Plot y-axis as log
plotCI = False # Plot confidence interval
plotOutline = False # Plot bar outline
plotPCT = False # Plot seasonality as percent

#endregion


# %% TEPHRA ____________________________________________________
# region Tephra
print('Tephra exposure')

# Exposure
# expCol = ['pop_count', 'area_crops', 'area_urban', 'RNDS', 'buildingsLoss']
# expName = ['Population count\n1 kg/m$^2$', 'Crop area (km$^2$)\n5 kg/m$^2$', 'Urban area (km$^2$)\n5 kg/m$^2$', 'Road network\ndisruption score', 'Buildings loss\n($x10^6$ USD)']
expCol = ['pop_count', 'n_buildings', 'length_roads', 'area_crops', 'area_urban']
expName = ['Population (n)\n1 kg/m$^2$', 'Buildings (n)\n100 kg/m$^2$', 'Road length (km)\n1 kg/m$^2$', 'Crop area (km$^2$)\n5 kg/m$^2$', 'Urban area (km$^2$)\n1 kg/m$^2$']
massT = [1,100,1,5,1]

# tephra = data.loc['Tephra', :, 50]
tephra = data.loc['Tephra']
tephra = tephra.loc[tephra.month == 0]
# Replace regions where all rows are 0, which causes trouble for log plotting
tephra[['length_roads', 'n_buildings']] = tephra[['length_roads', 'n_buildings']].replace({0:np.nan})

fig, axes = plt.subplots(nrows=len(expName), ncols=len(regionsRatio), sharey='row',
                               gridspec_kw={'width_ratios': regionsRatio},
                               figsize=(10,6.5))

for iRow in range(0, len(expName)):
    
    # If exposure, take the 1 kg/m2
    # if (iRow == 0) or (iRow == 3):
    #     dataRow = tephra.loc[tephra.mass == 1]
    # elif (iRow == 1) or (iRow == 2):
    #     dataRow = tephra.loc[tephra.mass == 5]
    # else:
    #     dataRow = tephra.loc[tephra.mass == 100]
    dataRow = tephra.loc[tephra.mass == massT[iRow]]
    
    for iCol in range(0, len(regionsRatio)):
        # Colormap
        if iCol == len(regionsRatio)-1:
            cmap = CMAP[1] # Purples
            cmapb = CMAPb[1] # Purples
            cmapi = CMAPi[1]
        else:
            cmap = CMAP[0]
            cmapb = CMAPb[0]
            cmapi = CMAPi[0]
        
        # Datacol
        dataCol = dataRow[dataRow.region==regions[iCol]]
        dataCol = dataCol.sort_index().sort_values(['volcano'])
        
        # Plot
        vCount = 0
        for iV in [5,4,3]:
            sns.barplot(ax=axes[iRow][iCol], x='volcano', y=expCol[iRow], data = dataCol.loc[iV,50], label=f'VEI {iV}', color=cmap[vCount], ci=None)
            axes[iRow][iCol].axhline(y=dataCol.loc[iV,50][expCol[iRow]].median(), color=cmapi[vCount], linestyle=':')
            if plotCI:
                yerr = [ dataCol.loc[(iV,90), expCol[iRow]], dataCol.loc[(iV,10), expCol[iRow]]]
                axes[iRow][iCol].errorbar(x='volcano', y=expCol[iRow], data = dataCol.loc[iV,50], yerr=yerr, fmt='none', c=cmapb[vCount])
            vCount += 1


        # Control y label
        if iCol == 0:
            axes[iRow][iCol].set(ylabel=expName[iRow])
        else:
            axes[iRow][iCol].set(ylabel=None)  
        
        # Control x label
        if iRow == len(expName)-1:
            axes[iRow][iCol].set_xticklabels(dataCol.volcano.unique(), rotation=90)#, ha='right')
        else:
            axes[iRow][iCol].set_xticklabels([])
        
        # Control column title
        if iRow == 0:
            axes[iRow][iCol].set_title(regions[iCol])

sns.despine()
plt.tight_layout()
fig.subplots_adjust(wspace=.25,hspace=.25)
if plotCI:
    fig.suptitle('Exposure to tephra accumulation', y=1.02, fontweight='bold')
else:
    fig.suptitle('Exposure to tephra accumulation ($P=50\%$)', y=1.02, fontweight='bold')
    
for ax in fig.axes:
    # ax.tick_params(labelrotation=45)
    if plotLog:
        ax.set_yscale("log") 
    ax.set_xlabel(None)
    plt.setp(ax.patches, linewidth=.3)
    if plotOutline:
        plt.setp(ax.patches, linewidth=0)

if plotLog:
    plt.savefig(f"{outPth}tephra.{figType}",bbox_inches='tight', dpi=300)
else:
    plt.savefig(f"{outPth}tephra_linear.{figType}",bbox_inches='tight', dpi=300)

#endregion

# %% BUFFER ____________________________________________________
# region Buffer
print('Buffer exposure')
# Exposure
expCol = ['pop_count']
expName = ['Population (n)']
# expCol = ['pop_count', 'n_buildings', 'length_roads', 'area_crops', 'area_urban']
# expName = ['Population (n)', 'Buildings (n)', 'Road length (km)', 'Crop area (km$^2$)', 'Urban area (km$^2$)']
buffer = data.loc['BUF', :, :].reset_index().set_index(['radius'])

# tephra = data.loc['Tephra', :, 50]

fig, axes = plt.subplots(nrows=len(expName), ncols=len(regionsRatio), sharey='row',
                               gridspec_kw={'width_ratios': regionsRatio},
                               figsize=(10,3.2))
r10 = mpatches.Patch(color=cmap[0], label='Radius: 10 km')
r30 = mpatches.Patch(color=cmap[0], label='Radius: 30 km')
r100 = mpatches.Patch(color=cmap[0], label='Radius: 100 km')

storLegLine = {}
for iRow in range(0, len(expName)):

    for iCol in range(0, len(regionsRatio)):
        # Colormap
        if iCol == len(regionsRatio)-1:
            cmap = CMAP[1] # Purples
            cmapb = CMAPb[1] # Purples
            cmapi = CMAPi[1] # Purples
        else:
            cmap = CMAP[0]
            cmapb = CMAPb[0]
            cmapi = CMAPi[0]
            
        # Datacol
        dataCol = buffer[buffer.region==regions[iCol]]
        dataCol = dataCol.sort_index().sort_values(['volcano'])
            
        # Plot
        vCount = 0
        for iV in [100,30,10]:
            sns.barplot(ax=axes[iCol], x='volcano', y=expCol[iRow], data = dataCol.loc[iV], label=f'Buffer: {iV} km', color=cmap[vCount], ci=None)
            storLegLine[iV] = axes[iCol].axhline(y=dataCol.loc[iV][expCol[iRow]].median(), color=cmapi[vCount], label=f'Median: {iV} km', linestyle=':')
            vCount += 1

        # Control y label
        if iCol == 0:
            axes[iCol].set(ylabel=expName[iRow])
        else:
            axes[iCol].set(ylabel=None)  
        
        # Control x label
        if iRow == len(expName)-1:
            axes[iCol].set_xticklabels(dataCol.volcano.unique(), rotation=90)#, ha='right')
        else:
            axes[iCol].set_xticklabels([])
        
        # Control column title
        if iRow == 0:
            axes[iCol].set_title(regions[iCol])
            
# plt.legend(handles=[r10,r30,r100]+[storLegLine[i] for i in storLegLine.keys()], title='title', bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2)

sns.despine()
plt.tight_layout()
fig.subplots_adjust(wspace=.25,hspace=.25)
fig.suptitle('Exposure as concentric radii', y=1.02, fontweight='bold')
for ax in fig.axes:
    # ax.tick_params(labelrotation=45)
    if plotLog:
        ax.set_yscale("log") 
    ax.set_xlabel(None)
    plt.setp(ax.patches, linewidth=.3)
    if plotOutline:
        plt.setp(ax.patches, linewidth=0)

if plotLog:
    plt.savefig(f"{outPth}buffer.{figType}",bbox_inches='tight', dpi=300)
else:
    plt.savefig(f"{outPth}buffer_linear.{figType}",bbox_inches='tight', dpi=300)

#endregion

#%% Tephra seasonal
#region tephra seasonal
print('Tephra seasonality')

tephra = data.loc['Tephra']
tephra = tephra.loc[tephra.month == 0]
tephraM = dataM.loc['Tephra']
tephraM = tephraM.set_index('month',append=True)

# Exposure
# expCol = ['pop_count', 'n_buildings', 'length_roads', 'area_crops', 'area_urban']
# expName = ['Population (n)\n1 kg/m$^2$', 'Buildings (n)\n100 kg/m$^2$', 'Road length (km)\n1 kg/m$^2$', 'Crop area (km$^2$)\n5 kg/m$^2$', 'Urban area (km$^2$)\n5 kg/m$^2$']
# Exposure
expCol = ['pop_count']
expName = ['Population (n)\n1 kg/m$^2$']

fig, axes = plt.subplots(nrows=len(expName), ncols=len(regionsRatio), sharey='row',
                               gridspec_kw={'width_ratios': regionsRatio},
                               figsize=(10,3.2))

for iRow in range(0, len(expName)):
    
    # If exposure, take the 1 kg/m2
    dataRowM = tephraM.loc[tephraM.mass == 1]
    dataRow = tephra.loc[tephra.mass == 1]

        
    for iCol in range(0, len(regionsRatio)):
        # Colormap
        if iCol == len(regionsRatio)-1:
            cmap = CMAP[1] # Purples
            cmapb = CMAPb[1] # Purples
            cmapi = CMAPi[1] # Purples
        else:
            cmap = CMAP[0]
            cmapb = CMAPb[0]
            cmapi = CMAPi[0] # Purples
            
        # Datacol
        dataColM = dataRowM[dataRowM.region==regions[iCol]]
        dataColM = dataColM.sort_index().sort_values(['volcano'])
        # Datacol
        dataCol = dataRow[dataRow.region==regions[iCol]]
        dataCol = dataCol.sort_index().sort_values(['volcano'])
            
        # Plot
        vCount = 0
        for iV in range(1,13):
            valTmp = dataCol.loc[(4,50), expCol[iRow]].values - dataColM.loc[(4,50,iV), expCol[iRow]].values
            if plotPCT:
                valTmp = valTmp/dataCol.loc[(4,50), expCol[iRow]].values
            
            sns.lineplot(ax=axes[iCol], x=dataColM.volcano.unique(), y=valTmp, label=f'VEI {iV}', color=CMAPm[vCount], ci=None, legend=False)
            vCount += 1

        # Control y label
        if iCol == 0:
            axes[iCol].set(ylabel=expName[iRow])
        else:
            axes[iCol].set(ylabel=None)  
        
        # Control x label
        if iRow == len(expName)-1:
            axes[iCol].set_xticklabels(dataColM.volcano.unique(), rotation=90)#, ha='right')
        else:
            axes[iCol].set_xticklabels([])
        
        # Control column title
        if iRow == 0:
            axes[iCol].set_title(regions[iCol])

sns.despine()
plt.tight_layout()
fig.subplots_adjust(wspace=.25,hspace=.25)
if plotPCT:
    fig.suptitle('Relative (%) variability of exposure to tephra accumulation with month ($VEI=4$, $P=50\%$)', y=1.02, fontweight='bold')
else:
    fig.suptitle('Absolute variability of exposure to tephra accumulation with month ($VEI=4$, $P=50\%$)', y=1.02, fontweight='bold')
for ax in fig.axes:
    # ax.tick_params(labelrotation=45)
    if plotLog:
        ax.set_yscale("log") 
    ax.set_xlabel(None)
    if plotOutline:
        plt.setp(ax.patches, linewidth=0)

if plotLog:
    plt.savefig(f"{outPth}seasonal_variability.{figType}",bbox_inches='tight', dpi=300)
else:
    plt.savefig(f"{outPth}seasonal_variability_linear.{figType}",bbox_inches='tight', dpi=300)
#endregion

#%% Column collapse
# region CC

# Exposure
# expCol = ['pop_count', 'area_crops', 'area_urban', 'RNDS']
# expName = ['Population count\n$', 'Crop area (km$^2$)', 'Urban area (km$^2$)', 'Road network\ndisruption score']
expCol = ['pop_count', 'n_buildings', 'length_roads', 'area_crops', 'area_urban']
expName = ['Population (n)', 'Buildings (n)', 'Road length (km)', 'Crop area (km$^2$)', 'Urban area (km$^2$)']

pdc = data.loc['PDC', :, :]
pdc[['length_roads', 'n_buildings']] = pdc[['length_roads', 'n_buildings']].replace({0:np.nan})

fig, axes = plt.subplots(nrows=len(expName), ncols=len(regionsRatio), sharey='row',
                               gridspec_kw={'width_ratios': regionsRatio},
                               figsize=(10,6.5))

for iRow in range(0, len(expName)):
    
    dataRow = pdc.copy()
        
    for iCol in range(0, len(regionsRatio)):
        
        # Colormap
        if iCol == len(regionsRatio)-1:
            cmap = CMAP[1] # Purples
            cmapi = CMAPi[1] # Purples
        else:
            cmap = CMAP[0]
            cmapi = CMAPi[0] # Purples
            
        # Datacol
        dataCol = dataRow[dataRow.region==regions[iCol]]
            
        # Plot
        vCount = 0
        for iV in [5,4,3]:
            sns.barplot(ax=axes[iRow][iCol], x='volcano', y=expCol[iRow], data = dataCol.loc[iV,50], label=f'VEI {iV}', color=cmap[vCount], ci=None)
            axes[iRow][iCol].axhline(y=dataCol.loc[iV,50][expCol[iRow]].median(), color=cmapi[vCount], linestyle=':')
            if plotCI:
                yerr = [ dataCol.loc[(iV,90), expCol[iRow]], dataCol.loc[(iV,10), expCol[iRow]]]
                axes[iRow][iCol].errorbar(x='volcano', y=expCol[iRow], data = dataCol.loc[iV,50], yerr=yerr, fmt='none', c=cmapb[vCount])
            vCount += 1

        # Control y label
        if iCol == 0:
            axes[iRow][iCol].set(ylabel=expName[iRow])
        else:
            axes[iRow][iCol].set(ylabel=None)  
        
        # Control x label
        if iRow == len(expName)-1:
            axes[iRow][iCol].set_xticklabels(dataCol.volcano.unique(), rotation=90)#, ha='right')
        else:
            axes[iRow][iCol].set_xticklabels([])
        
        # Control column title
        if iRow == 0:
            axes[iRow][iCol].set_title(regions[iCol])

sns.despine()
plt.tight_layout()
fig.subplots_adjust(wspace=.25,hspace=.25)
fig.suptitle('Exposure to PDC from column collapse ($P=50\%$)', y=1.02, fontweight='bold')
for ax in fig.axes:
    # ax.tick_params(labelrotation=45)
    if plotLog:
        ax.set_yscale("log") 
    ax.set_xlabel(None)
    plt.setp(ax.patches, linewidth=.3)
    if plotOutline:
        plt.setp(ax.patches, linewidth=0)

if plotLog:
    plt.savefig(f"{outPth}column_collapse.{figType}",bbox_inches='tight', dpi=300)
else:
    plt.savefig(f"{outPth}column_collapse_linear.{figType}",bbox_inches='tight', dpi=300)
    
#endregion

# %% Large clasts
#region LC

# Exposure
# expCol = ['pop_count', 'area_crops', 'area_urban']
# expName = ['Population count', 'Crop area (km$^2$)', 'Urban area (km$^2$)']
expCol = ['pop_count', 'n_buildings', 'length_roads', 'area_crops', 'area_urban']
expName = ['Population (n)', 'Buildings (n)', 'Road length (km)', 'Crop area (km$^2$)', 'Urban area (km$^2$)']

lc = data.loc['LC', :, :]
lc[['length_roads', 'n_buildings']] = lc[['length_roads', 'n_buildings']].replace({0:np.nan})

fig, axes = plt.subplots(nrows=len(expName), ncols=len(regionsRatio), sharey='row',
                               gridspec_kw={'width_ratios': regionsRatio},
                               figsize=(10,6.5))

for iRow in range(0, len(expName)):
    
    dataRow = lc.copy()
        
    for iCol in range(0, len(regionsRatio)):
        
        # Colormap
        if iCol == len(regionsRatio)-1:
            cmap = CMAP[1] # Purples
            cmapi = CMAPi[1] # Purples
        else:
            cmap = CMAP[0]
            cmapi = CMAPi[0] # Purples
            
        # Datacol
        dataCol = dataRow[dataRow.region==regions[iCol]]
            
        # Plot
        vCount = 0
        for iV in [5,4,3]:
            sns.barplot(ax=axes[iRow][iCol], x='volcano', y=expCol[iRow], data = dataCol.loc[iV,50], label=f'VEI {iV}', color=cmap[vCount], ci=None)
            axes[iRow][iCol].axhline(y=dataCol.loc[iV,50][expCol[iRow]].median(), color=cmapi[vCount], linestyle=':')
            if plotCI:
                yerr = [ dataCol.loc[(iV,90), expCol[iRow]], dataCol.loc[(iV,10), expCol[iRow]]]
                axes[iRow][iCol].errorbar(x='volcano', y=expCol[iRow], data = dataCol.loc[iV,50], yerr=yerr, fmt='none', c=cmapb[vCount])
            vCount += 1
            
        # Control y label
        if iCol == 0:
            axes[iRow][iCol].set(ylabel=expName[iRow])
        else:
            axes[iRow][iCol].set(ylabel=None)  
        
        # Control x label
        if iRow == len(expName)-1:
            axes[iRow][iCol].set_xticklabels(dataCol.volcano.unique(), rotation=90)#, ha='right')
        else:
            axes[iRow][iCol].set_xticklabels([])
        
        # Control column title
        if iRow == 0:
            axes[iRow][iCol].set_title(regions[iCol])

sns.despine()
plt.tight_layout()
fig.subplots_adjust(wspace=.25,hspace=.25)
fig.suptitle('Exposure to lapilli impacts ($P=50\%$)', y=1.02, fontweight='bold')

for ax in fig.axes:
    # ax.tick_params(labelrotation=45)
    if plotLog:
        ax.set_yscale("log") 
    ax.set_xlabel(None)
    plt.setp(ax.patches, linewidth=.3)
    if plotOutline:
        plt.setp(ax.patches, linewidth=0)

if plotLog:
    plt.savefig(f"{outPth}large_clast.{figType}",bbox_inches='tight', dpi=300)
else:
    plt.savefig(f"{outPth}large_clast_linear.{figType}",bbox_inches='tight', dpi=300)
#endregion

#%% BAF
#region BAF


# Exposure
# expCol = ['pop_count', 'area_crops', 'area_urban', 'RNDS']
# expName = ['Population count', 'Crop area (km$^2$)', 'Urban area (km$^2$)', 'Road network\ndisruption score']
expCol = ['pop_count', 'n_buildings', 'length_roads', 'area_crops', 'area_urban']
expName = ['Population (n)', 'Buildings (n)', 'Road length (km)', 'Crop area (km$^2$)', 'Urban area (km$^2$)']

baf = data.loc['BAF', :, :]
baf = baf.reset_index().set_index(['volume', 'prob'])
baf = baf[baf['buffer']==990]
baf[['length_roads', 'n_buildings']] = baf[['length_roads', 'n_buildings']].replace({0:np.nan})

fig, axes = plt.subplots(nrows=len(expName), ncols=len(regionsRatio), sharey='row',
                               gridspec_kw={'width_ratios': regionsRatio},
                               figsize=(10,6.5))

for iRow in range(0, len(expName)):
    
    dataRow = baf.copy()
        
    for iCol in range(0, len(regionsRatio)):
        # Colormap
        if iCol == len(regionsRatio)-1:
            cmap = CMAP[1] # Purples
            cmapi = CMAPi[1] # Purples
        else:
            cmap = CMAP[0]
            cmapi = CMAPi[0] # Purples
            
        # Datacol
        dataCol = dataRow[dataRow.region==regions[iCol]]
        dataCol = dataCol.sort_index().sort_values(['volcano'])
        
        
        # Plot
        vCount = 0
        for iV in [9800000,450000]:
            sns.barplot(ax=axes[iRow][iCol], x='volcano', y=expCol[iRow], data = dataCol.loc[iV,50], label=f'VEI {iV}', color=cmap[vCount], ci=None)
            axes[iRow][iCol].axhline(y=dataCol.loc[iV,50][expCol[iRow]].median(), color=cmapi[vCount], linestyle=':')
            if plotCI:
                yerr = [ dataCol.loc[(iV,90), expCol[iRow]], dataCol.loc[(iV,10), expCol[iRow]]]
                axes[iRow][iCol].errorbar(x='volcano', y=expCol[iRow], data = dataCol.loc[iV,50], yerr=yerr, fmt='none', c=cmapb[vCount])
            vCount += 1
            
        # Control y label
        if iCol == 0:
            axes[iRow][iCol].set(ylabel=expName[iRow])
        else:
            axes[iRow][iCol].set(ylabel=None)  
        
        # Control x label
        if iRow == len(expName)-1:
            axes[iRow][iCol].set_xticklabels(dataCol.volcano.unique(), rotation=90)#, ha='right')
        else:
            axes[iRow][iCol].set_xticklabels([])
        
        # Control column title
        if iRow == 0:
            axes[iRow][iCol].set_title(regions[iCol])
        # if (iRow == 2) or (iRow==1):
        #     axes[iRow][iCol].set_ylim([0.9, 5])

sns.despine()
plt.tight_layout()
fig.subplots_adjust(wspace=.25,hspace=.25)
fig.suptitle('Exposure to block and ash flows ($P=50\%$, Buffer$=990 m$)', y=1.02, fontweight='bold')

for ax in fig.axes:
    # ax.tick_params(labelrotation=45)
    if plotLog:
        ax.set_yscale("log") 
    ax.set_xlabel(None)
    plt.setp(ax.patches, linewidth=.3)
    if plotOutline:
        plt.setp(ax.patches, linewidth=0)

if plotLog:
    plt.savefig(f"{outPth}baf.{figType}",bbox_inches='tight', dpi=300)
else:
    plt.savefig(f"{outPth}baf_linear.{figType}",bbox_inches='tight', dpi=300)
#endregion


# %%

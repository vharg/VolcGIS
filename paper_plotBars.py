#%%
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['axes.linewidth'] = .4
matplotlib.rcParams['xtick.major.width'] = .4
matplotlib.rcParams['ytick.major.width'] = .4
# matplotlib.rcParams['text.color'] = 0
# matplotlib.rcParams['axes.labelcolor'] = 0

import seaborn as sns
sns.set_style("ticks", {"axes.facecolor": ".9"})
sns.set_context("paper")
from matplotlib import gridspec

#%% Load
#region preproces
print('Read master exposure')
data = pd.read_csv('MASTER_exposure.csv', index_col=None).drop('Unnamed: 0', axis=1)
buffer = pd.read_csv('MASTER_buffer.csv', index_col=None)

# Read volcano index with country names
volcDB = pd.read_csv('csv/volcCoordinates.csv')
volcDB = volcDB.rename({'name': 'volcano'}, axis=1)

def mergeData(data, db):
    db = db.set_index('volcano')[['country', 'region']]
    data = data.set_index('volcano').join(db).reset_index()
    data = data.replace('Philippines and SE Asia', 'Philippines')
    data = data.sort_values(by=['country', 'region', 'volcano'])
    return data

# Merge data and set index
data = mergeData(data, volcDB)
data = data.set_index(['hazard', 'VEI', 'prob'])

# Format string
data['volcano'] = data['volcano'].str.replace('_',' ')
data['volcano'] = data['volcano'].str.replace('Leroboleng','Lereboleng')
data['region'] = data['region'].str.replace('Sumatera','Sumatra')
data['region'] = data['region'].str.replace('The Philippines','Philippines')
data['region'] = data['region'].str.replace('Maluku/Halmahera','Halmahera/Banda Sea')
data = data.sort_values(by=['country', 'region', 'volcano'])
# Region
buffer = mergeData(buffer, volcDB)
buffer = buffer.set_index('buffer')
buffer['region'] = buffer['region'].str.replace('Sumatera','Sumatra')
buffer['region'] = buffer['region'].str.replace('The Philippines','Philippines')
buffer['region'] = buffer['region'].str.replace('Maluku/Halmahera','Halmahera/Banda Sea')
buffer = buffer.sort_values(by=['country', 'region', 'volcano'])

# Regions
regions = data.region.unique()
regionsRatio = [data[data.region==reg].volcano.unique().shape[0] for reg in regions]

# Fig size
figW = 9
figType = 'png'

# Colormaps
reverseCM = True
CMAP = [['#bdd7e7','#6baed6','#2171b5'], ['#cbc9e2','#9e9ac8','#6a51a3']]
CMAPb = [['#fdbe85','#fd8d3c','#e6550d'], ['#78c679','#31a354','#006837']]
CMAPm = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']

if reverseCM:
    CMAP = [i[::-1] for i in CMAP]

# Options
plotLog = True # Plot y-axis as log
plotCI = False # Plot confidence interval
plotOutline = False # Plot bar outline
plotPCT = True # Plot seasonality as percent

#endregion


# %% TEPHRA ____________________________________________________
# region Tephra
print('Tephra exposure')

# Exposure
expCol = ['pop_count', 'area_crops', 'area_urban', 'RNDS', 'buildingsLoss']
expName = ['Population count\n1 kg/m$^2$', 'Crop area (km$^2$)\n5 kg/m$^2$', 'Urban area (km$^2$)\n5 kg/m$^2$', 'Road network\ndisruption score', 'Buildings loss\n($x10^6$ USD)']


# tephra = data.loc['Tephra', :, 50]
tephra = data.loc['Tephra']
tephra = tephra.loc[tephra.month == 0]

fig, axes = plt.subplots(nrows=len(expName), ncols=len(regionsRatio), sharey='row',
                               gridspec_kw={'width_ratios': regionsRatio},
                               figsize=(figW,7))

for iRow in range(0, len(expName)):
    
    # If exposure, take the 1 kg/m2
    if iRow == 0:
        dataRow = tephra.loc[tephra.mass == 1]
    elif (iRow == 1) or (iRow == 2):
        dataRow = tephra.loc[tephra.mass == 5]
    else:
        dataRow = tephra.loc[tephra.mass.isnull()]
        
    for iCol in range(0, len(regionsRatio)):
        # Colormap
        if iCol == len(regionsRatio)-1:
            cmap = CMAP[1] # Purples
            cmapb = CMAPb[1] # Purples
        else:
            cmap = CMAP[0]
            cmapb = CMAPb[0]
            
        # Datacol
        dataCol = dataRow[dataRow.region==regions[iCol]]
        dataCol = dataCol.sort_index().sort_values(['volcano'])
            
        # Plot
        vCount = 0
        for iV in [5,4,3]:
            sns.barplot(ax=axes[iRow][iCol], x='volcano', y=expCol[iRow], data = dataCol.loc[iV,50], label=f'VEI {iV}', color=cmap[vCount], ci=None)
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
    if plotOutline:
        plt.setp(ax.patches, linewidth=0)

if plotLog:
    plt.savefig(f"Figures/tephra.{figType}",bbox_inches='tight', dpi=300)
else:
    plt.savefig(f"Figures/tephra_linear.{figType}",bbox_inches='tight', dpi=300)

#endregion

# %% BUFFER ____________________________________________________
# region Buffer
print('Buffer exposure')
# Exposure
expCol = ['pop_count']
expName = ['Population count']


# tephra = data.loc['Tephra', :, 50]

fig, axes = plt.subplots(nrows=len(expName), ncols=len(regionsRatio), sharey='row',
                               gridspec_kw={'width_ratios': regionsRatio},
                               figsize=(10,4))

for iRow in range(0, len(expName)):

    for iCol in range(0, len(regionsRatio)):
        # Colormap
        if iCol == len(regionsRatio)-1:
            cmap = CMAP[1] # Purples
            cmapb = CMAPb[1] # Purples
        else:
            cmap = CMAP[0]
            cmapb = CMAPb[0]
            
        # Datacol
        dataCol = buffer[buffer.region==regions[iCol]]
        dataCol = dataCol.sort_index().sort_values(['volcano'])
            
        # Plot
        vCount = 0
        for iV in [100,30,10]:
            sns.barplot(ax=axes[iCol], x='volcano', y=expCol[iRow], data = dataCol.loc[iV], label=f'Buffer: {iV}', color=cmap[vCount], ci=None)
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

sns.despine()
plt.tight_layout()
fig.subplots_adjust(wspace=.25,hspace=.25)
fig.suptitle('Exposure as circular buffers', y=1.02, fontweight='bold')
for ax in fig.axes:
    # ax.tick_params(labelrotation=45)
    if plotLog:
        ax.set_yscale("log") 
    ax.set_xlabel(None)
    if plotOutline:
        plt.setp(ax.patches, linewidth=0)

if plotLog:
    plt.savefig(f"Figures/buffer.{figType}",bbox_inches='tight', dpi=300)
else:
    plt.savefig(f"Figures/buffer_linear.{figType}",bbox_inches='tight', dpi=300)

#endregion

#%% Tephra seasonal
#region tephra seasonal
print('Tephra seasonality')

tephra = data.loc['Tephra']
tephra = tephra.set_index('month',append=True)
# Exposure
expCol = ['pop_count', 'area_crops', 'area_urban', 'RNDS', 'buildingsLoss']
expName = ['Population count\nexp. to 1 kg/m$^2$', 'Crop area (km$^2$)\nexp. to 1 kg/m$^2$', 'Urban area (km$^2$)\nexp. to 1 kg/m$^2$', 'Road network\ndisruption score', 'Buildings loss\n($x10^6$ USD)']


fig, axes = plt.subplots(nrows=len(expName), ncols=len(regionsRatio), sharey='row',
                               gridspec_kw={'width_ratios': regionsRatio},
                               figsize=(figW,7))

for iRow in range(0, len(expName)):
    
    # If exposure, take the 1 kg/m2
    if iRow <= 2:
        dataRow = tephra.loc[tephra.mass == 1]
    else:
        dataRow = tephra.loc[tephra.mass.isnull()]
        
    for iCol in range(0, len(regionsRatio)):
        # Colormap
        if iCol == len(regionsRatio)-1:
            cmap = CMAP[1] # Purples
            cmapb = CMAPb[1] # Purples
        else:
            cmap = CMAP[0]
            cmapb = CMAPb[0]
            
        # Datacol
        dataCol = dataRow[dataRow.region==regions[iCol]]
        dataCol = dataCol.sort_index().sort_values(['volcano'])
            
        # Plot
        vCount = 0
        for iV in range(1,13):
            valTmp = dataCol.loc[(4,50,0), expCol[iRow]].values - dataCol.loc[(4,50,iV), expCol[iRow]].values
            if plotPCT:
                valTmp = valTmp/dataCol.loc[(4,50,0), expCol[iRow]].values
            
            sns.lineplot(ax=axes[iRow][iCol], x=dataCol.volcano.unique(), y=valTmp, label=f'VEI {iV}', color=CMAPm[vCount], ci=None, legend=False)
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
    plt.savefig(f"Figures/seasonal_variability.{figType}",bbox_inches='tight', dpi=300)
else:
    plt.savefig(f"Figures/seasonal_variability_linear.{figType}",bbox_inches='tight', dpi=300)
#endregion

#%% Column collapse
# region CC

# Exposure
expCol = ['pop_count', 'area_crops', 'area_urban', 'RNDS']
expName = ['Population count\nexp. to 1 kg/m$^2$', 'Crop area (km$^2$)\nexp. to 1 kg/m$^2$', 'Urban area (km$^2$)\nexp. to 1 kg/m$^2$', 'Road network\ndisruption score']

pdc = data.loc['PDC', :, :]

fig, axes = plt.subplots(nrows=len(expName), ncols=len(regionsRatio), sharey='row',
                               gridspec_kw={'width_ratios': regionsRatio},
                               figsize=(figW,7.5))

for iRow in range(0, len(expName)):
    
    dataRow = pdc.copy()
        
    for iCol in range(0, len(regionsRatio)):
        
        # Colormap
        if iCol == len(regionsRatio)-1:
            cmap = CMAP[1] # Purples
        else:
            cmap = CMAP[0]
            
        # Datacol
        dataCol = dataRow[dataRow.region==regions[iCol]]
            
        # Plot
        vCount = 0
        for iV in [5,4,3]:
            sns.barplot(ax=axes[iRow][iCol], x='volcano', y=expCol[iRow], data = dataCol.loc[iV,50], label=f'VEI {iV}', color=cmap[vCount], ci=None)
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
    if plotOutline:
        plt.setp(ax.patches, linewidth=0)

if plotLog:
    plt.savefig(f"Figures/column_collapse.{figType}",bbox_inches='tight', dpi=300)
else:
    plt.savefig(f"Figures/column_collapse_linear.{figType}",bbox_inches='tight', dpi=300)
    
#endregion


# %% Large clasts
#region LC

# Exposure
expCol = ['pop_count', 'area_crops', 'area_urban']
expName = ['Population count\nexp. to 1 kg/m$^2$', 'Crop area (km$^2$)\nexp. to 1 kg/m$^2$', 'Urban area (km$^2$)\nexp. to 1 kg/m$^2$']

lc = data.loc['LC', :, :]

fig, axes = plt.subplots(nrows=len(expName), ncols=len(regionsRatio), sharey='row',
                               gridspec_kw={'width_ratios': regionsRatio},
                               figsize=(10,6.5))

for iRow in range(0, len(expName)):
    
    dataRow = lc.copy()
        
    for iCol in range(0, len(regionsRatio)):
        
        # Colormap
        if iCol == len(regionsRatio)-1:
            cmap = CMAP[1] # Purples
        else:
            cmap = CMAP[0]
            
        # Datacol
        dataCol = dataRow[dataRow.region==regions[iCol]]
            
        # Plot
        vCount = 0
        for iV in [5,4,3]:
            sns.barplot(ax=axes[iRow][iCol], x='volcano', y=expCol[iRow], data = dataCol.loc[iV,50], label=f'VEI {iV}', color=cmap[vCount], ci=None)
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
    if plotOutline:
        plt.setp(ax.patches, linewidth=0)

if plotLog:
    plt.savefig(f"Figures/large_clast.{figType}",bbox_inches='tight', dpi=300)
else:
    plt.savefig(f"Figures/large_clast_linear.{figType}",bbox_inches='tight', dpi=300)
#endregion

#%% BAF
#region BAF


# Exposure
expCol = ['pop_count', 'area_crops', 'area_urban', 'RNDS']
expName = ['Population count', 'Crop area (km$^2$)', 'Urban area (km$^2$)', 'Road network\ndisruption score']

baf = data.loc['BAF', :, :]
baf = baf.reset_index().set_index(['volume', 'prob'])
baf = baf[baf['buffer']==990]

fig, axes = plt.subplots(nrows=len(expName), ncols=len(regionsRatio), sharey='row',
                               gridspec_kw={'width_ratios': regionsRatio},
                               figsize=(10,6.5))

for iRow in range(0, len(expName)):
    
    dataRow = baf.copy()
        
    for iCol in range(0, len(regionsRatio)):
        print(dataCol.area_urban.unique())
        # Colormap
        if iCol == len(regionsRatio)-1:
            cmap = CMAP[1] # Purples
        else:
            cmap = CMAP[0]
            
        # Datacol
        dataCol = dataRow[dataRow.region==regions[iCol]]
        dataCol = dataCol.sort_index().sort_values(['volcano'])
        
        
        # Plot
        vCount = 0
        for iV in [9800000,450000]:
            sns.barplot(ax=axes[iRow][iCol], x='volcano', y=expCol[iRow], data = dataCol.loc[iV,50], label=f'VEI {iV}', color=cmap[vCount], ci=None)
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
        if (iRow == 2) or (iRow==1):
            axes[iRow][iCol].set_ylim([0.9, 5])

sns.despine()
plt.tight_layout()
fig.subplots_adjust(wspace=.25,hspace=.25)
fig.suptitle('Exposure to block and ash flows ($P=50\%$, Buffer$=990 m$)', y=1.02, fontweight='bold')

for ax in fig.axes:
    # ax.tick_params(labelrotation=45)
    if plotLog:
        ax.set_yscale("log") 
    ax.set_xlabel(None)
    if plotOutline:
        plt.setp(ax.patches, linewidth=0)

if plotLog:
    plt.savefig(f"Figures/baf.{figType}",bbox_inches='tight', dpi=300)
else:
    plt.savefig(f"Figures/baf_linear.{figType}",bbox_inches='tight', dpi=300)
#endregion


# %%

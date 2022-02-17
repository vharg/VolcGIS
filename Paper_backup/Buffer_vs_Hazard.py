#%%
# Makes the comparison between population exposure from hazard footprints and circles, i.e. our Fig. radii in the paper

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("ticks", {"axes.facecolor": ".9"})
sns.set_context("paper")
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
# Roads = Roads[Roads.month==0].drop('month', axis=1)
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
# data = data.set_index(['hazard', 'VEI', 'prob'])
data['n_buildings'] = data[['buildingsW', 'buildingsS']].sum(axis=1)
data = data.drop(['buildingsW', 'buildingsS'],axis=1)

# Format string
data['volcano'] = data['volcano'].str.replace('_',' ')
data['volcano'] = data['volcano'].str.replace('Leroboleng','Lereboleng')
data['region'] = data['region'].str.replace('Sumatera','Sumatra')
data['region'] = data['region'].str.replace('The Philippines','Philippines')
data['region'] = data['region'].str.replace('Maluku/Halmahera','Halmahera/Banda Sea')
# data['hazard'] = data['hazard'].str.replace('LC','Large clast')
# data['hazard'] = data['hazard'].str.replace('BUF','Circular buffers')
data = data.sort_values(by=['country', 'region', 'volcano'])
data = data[~((data['hazard']=='Tephra') & (data['mass'].isnull()))] # Remove the row

data = data.set_index(['volcano'])





# %%

fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(6,6), sharex=True, sharey=True)


VEI=[3,4,5]
PRB=50
MASS=[1,100]
RAD=[10,30,100]
EXP='pop_count'

for iRow in range(0,3):
    for iCol in range(0,3):
        axes[iRow][iCol].plot([1e3,1e8], [1e3,1e8], color='k', linestyle=':')
        for iRad in range(0,3):
            
            if iRow<2:
                HAZ = data[(data['hazard']=='Tephra') & (data['VEI']==VEI[iCol]) & (data['prob']==PRB) & (data['mass']==MASS[iRow])][[EXP]].join(
                    data[(data['hazard']=='BUF')&(data['radius']==RAD[iRad])][[EXP]], rsuffix='B')

                sns.regplot(data=HAZ, x=EXP+'B', y=EXP,ax=axes[iRow][iCol], label = f'r={RAD[iRad]} km',ci=None, fit_reg=False)
                axes[iRow][iCol].set(xlim=(1e3, 1e8), ylim=(1e3, 1e8))
                
                YLAB = f'Tephra {MASS[iRow]} kgm2'
                
            else:
                HAZ = data[(data['hazard']=='PDC') & (data['VEI']==VEI[iCol]) & (data['prob']==PRB)][[EXP]].join(
                    data[(data['hazard']=='BUF')&(data['radius']==RAD[iRad])][[EXP]], rsuffix='B'
                )

                sns.regplot(data=HAZ, x=EXP+'B', y=EXP,ax=axes[iRow][iCol], label = f'r={RAD[iRad]} km',ci=None, fit_reg=False)
                axes[iRow][iCol].set(xlim=(1e3, 1e8), ylim=(1e3, 1e8))
                
                YLAB = f'PDC - Column collapse'
                
                # Control y label
                
        
        
        if iCol == 0 :
            axes[iRow][iCol].set(ylabel=YLAB)
        else:
            axes[iRow][iCol].set(ylabel=None)  
            
        if iRow == 2:
            axes[iRow][iCol].set(xlabel='Concentric radii')
        else:
            axes[iRow][iCol].set(xlabel=None)
              
        if iRow == 0:
            axes[iRow][iCol].set_title(f'VEI {VEI[iCol]}')
        if iCol==2 and iRow==2:
            plt.legend(borderpad=-0.3, loc='upper left', frameon=False, framealpha=0)    
            # plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
            
        axes[iRow][iCol].set_xscale("log") 
        axes[iRow][iCol].set_yscale("log") 
        

    
plt.tight_layout()
fig.subplots_adjust(wspace=.1,hspace=.1)
fig.suptitle(f'Population count - P={PRB}%', y=1.02, fontweight='bold')

plt.savefig(f"Paper/Figures/radiusVShazard.pdf",bbox_inches='tight', dpi=300)
plt.savefig(f"Paper/Figures/radiusVShazard.png",bbox_inches='tight', dpi=300)
# %%

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(2,2), sharex=True, sharey=True)
for iRad in range(0,3):
    
    HAZ = data[(data['hazard']=='BAF') & (data['volume']==9800000) & (data['prob']==50) & (data['buffer']==990)][[EXP]].join(
        data[(data['hazard']=='BUF')&(data['radius']==RAD[iRad])][[EXP]], rsuffix='B')

    sns.regplot(data=HAZ, x=EXP+'B', y=EXP,ax=axes, label = f'r={RAD[iRad]} km',ci=None, fit_reg=False)
    axes.set(xlim=(1e3, 1e8), ylim=(1e3, 1e8))

axes.set_xscale("log") 
axes.set_yscale("log") 
    
# %%

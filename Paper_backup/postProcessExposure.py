#%% [markdown]
# ## Processing of the exposure analysis
# Reads the separate files `exposure.csv`, `damage_ratio.csv`, `damage_state.csv` and `road_exposure.csv`
# both for all synoptic monthly wind profiles and per month (tephra)
#
# Seb 2021-03-16

#%%
import pandas as pd
import os

#%%
volcDB = pd.read_csv('csv/volcCoordinates.csv')

EXPOSURE = pd.DataFrame()
DS = pd.DataFrame()
DS = pd.DataFrame()
ROAD = pd.DataFrame()

processMonthly = False

for i in range(0, volcDB.shape[0]):
    name = volcDB.loc[i, 'name']
    
    pth = os.path.join('volcanoes', name, '_exposure')
    
    # General exposure - all months
    exposure = pd.read_csv(os.path.join(pth, 'exposure.csv'), index_col=False ).drop('Unnamed: 0', axis=1)
    exposure['month'] = 0
    
    ds = pd.read_csv(os.path.join(pth, 'damage_state.csv'), index_col=False )
    ds['month'] = 0
    
    dr = pd.read_csv(os.path.join(pth, 'damage_ratio.csv'), index_col=False )
    dr['month'] = 0
    
    road = pd.read_csv(os.path.join(pth, 'road_exposure.csv'), index_col=False ).drop('Unnamed: 0', axis=1)
    road['month'] = 0
    
    if i == 0:
        EXPOSURE = exposure
        DS = ds
        DR = dr
        ROAD = road
    else:
        EXPOSURE = pd.concat([EXPOSURE, exposure])
        DS = pd.concat([DS, ds])
        DR = pd.concat([DR, dr])
        ROAD = pd.concat([ROAD, road])
        
    
    # Monthly exposure
    if processMonthly:
        for iM in range(1,13):
            
            tmp = pd.read_csv(os.path.join(pth, 'exposure_month{}.csv'.format(iM)), index_col=False ).drop('Unnamed: 0', axis=1)
            tmp['month'] = iM
            EXPOSURE = pd.concat([EXPOSURE, tmp])
            
            tmp = pd.read_csv(os.path.join(pth, 'damage_state_month{}.csv'.format(iM)), index_col=False )
            tmp['month'] = iM
            DS = pd.concat([DS, tmp])
            
            tmp = pd.read_csv(os.path.join(pth, 'damage_ratio_month{}.csv'.format(iM)), index_col=False )
            tmp['month'] = iM
            DR = pd.concat([DR, tmp])
            
            tmp = pd.read_csv(os.path.join(pth, 'road_exposure_month{}.csv'.format(iM)), index_col=False ).drop('Unnamed: 0', axis=1)
            tmp['month'] = iM
            ROAD = pd.concat([ROAD, tmp])

    
# %%
EXPOSURE.to_csv('MASTER_exposure.csv')
DS.to_csv('MASTER_damage_state.csv')
DR.to_csv('MASTER_damage_ratio.csv')
ROAD.to_csv('MASTER_roads.csv')
# %%

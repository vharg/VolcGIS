#%%
# This script computes the mean exposure for each hazard per region and exports the results as excel files.
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
sns.set_style("ticks", {"axes.facecolor": ".9"})
sns.set_context("paper")
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

#%%
def getData(VEI=4,VOL=450000,PROB=50,BUFFER=990,agg='mean',RADIUS=100):
    print(f'PROB={PROB}, VEI={VEI}, VOL={VOL}, BUFFER={BUFFER}')

    # Format data for BAD
    dataBaf = data.loc['BAF',:,PROB]
    dataBaf = dataBaf[(dataBaf.buffer==BUFFER) & (dataBaf.volume==VOL)]
    dataBaf.index = ['BAF']*dataBaf.shape[0]
    dataBaf.index.rename('hazard', inplace=True)

    # Format data for Buffers
    dataBuf = data.loc['BUF',:,:]
    dataBuf = dataBuf[(dataBuf.radius==RADIUS)]
    dataBuf.index = ['BUF']*dataBuf.shape[0]
    dataBuf.index.rename('hazard', inplace=True)
    
    # Format data for other hazards
    dataOther = data.loc[:,VEI,PROB]
    # dataOther.loc[((dataOther.mass!=50) & (dataOther.mass!=100))] 
    # dataOther[(dataOther.mass!=50) & (dataOther.mass!=100)]

    # Tephra
    tephra = dataOther.loc['Tephra'].set_index('mass')
    tephra_crops = tephra.loc[5][['area_crops','region']].groupby('region').agg(agg)
    tephra_crops['hazard'] = 'Tephra'
    tephra_urban = tephra.loc[5][['area_urban','region']].groupby('region').agg(agg)
    tephra_urban['hazard'] = 'Tephra'
    tephra_pop = tephra.loc[1][['pop_count','region']].groupby('region').agg(agg)
    tephra_pop['hazard'] = 'Tephra'
    tephra_roads = tephra.loc[1][['length_roads','region']].groupby('region').agg(agg)
    tephra_roads['hazard'] = 'Tephra'
    tephra_buildings = tephra.loc[100][['n_buildings','region']].groupby('region').agg(agg)
    tephra_buildings['hazard'] = 'Tephra'
    # tephra_RNDS = tephra[np.isnan(tephra.index.values)][['RNDS','region']].groupby('region').agg(agg)
    # tephra_RNDS['hazard'] = 'Tephra'
    # tephra_buildings = tephra[np.isnan(tephra.index.values)][['buildingsLoss','region']].groupby('region').agg(agg)
    # tephra_buildings['hazard'] = 'Tephra'

    tephra_crops.set_index('hazard', append=True, inplace=True)
    tephra_urban.set_index('hazard', append=True, inplace=True)
    tephra_pop.set_index('hazard', append=True, inplace=True)
    tephra_roads.set_index('hazard', append=True, inplace=True)
    tephra_buildings.set_index('hazard', append=True, inplace=True)
    # tephra_RNDS.set_index('hazard', append=True, inplace=True)
    # tephra_buildings.set_index('hazard', append=True, inplace=True)

    TEPHRA = tephra_pop.join(tephra_urban).join(tephra_crops).join(tephra_roads).join(tephra_buildings)
    # TEPHRA = tephra_pop.join(tephra_urban).join(tephra_crops).join(tephra_RNDS).join(tephra_buildings)
    
    # BAF
    # BAF = dataBaf[['pop_count','area_crops', 'area_urban', 'RNDS','region']].groupby('region').agg(agg)
    BAF = dataBaf[['pop_count','area_crops', 'area_urban', 'length_roads', 'n_buildings','region']].groupby('region').agg(agg)
    BAF['hazard'] = 'BAF'
    BAF.set_index('hazard', append=True, inplace=True)

    # PDC
    pdc = dataOther.loc['PDC']
    pdc = pdc[['pop_count','area_crops', 'area_urban', 'length_roads', 'n_buildings','region']].groupby('region').agg(agg)
    pdc['hazard'] = 'PDC'
    pdc.set_index('hazard', append=True, inplace=True)

    # LC
    LC = dataOther.loc['LC']
    LC = LC[['pop_count','area_crops', 'area_urban','length_roads', 'n_buildings','region']].groupby('region').agg(agg)
    LC['hazard'] = 'LC'
    LC.set_index('hazard', append=True, inplace=True)

    # buffer
    BUF = dataBuf[['pop_count','area_crops', 'area_urban', 'length_roads', 'n_buildings','region']].groupby('region').agg(agg)
    BUF['hazard'] = 'BUF'
    BUF.set_index('hazard', append=True, inplace=True)

    DATA = TEPHRA.append(pdc).append(LC).append(BAF).append(BUF).sort_index()
    return DATA

# DATA.style.format({
# 'pop_count': "{:.2E}",
# 'area_urban': "{:.2f}",
# 'area_crops': "{:.2f}",
# 'RNDS': "{:.2E}",
# 'buildingsLoss': "{:.2f}",
# }).loc[:,'Tephra']

#%% 
# agg = 'std'
agg = 'median'
DATA = getData(VEI=4, agg=agg)

hazard = 'Tephra'
# hazard = 'BUF'
# hazard = 'PDC'
# hazard = 'BAF'
# hazard = 'LC'

tmp=DATA.loc[pd.IndexSlice[:, hazard],:]
DATA.loc[pd.IndexSlice[:, hazard],:].rank(ascending=False)

# with pd.ExcelWriter('csv/PerRegion.xlsx') as writer: 
#     for i in [3,4,5]:
#         DATA = getData(VEI=i, agg=agg)
#         tmp = DATA.loc[pd.IndexSlice[:, 'Tephra'],:]
#         tmp.to_excel(writer, sheet_name=f'Tephra VEI {i} 50%')

#     for i in [3,4,5]:
#         DATA = getData(VEI=i, agg=agg)
#         tmp = DATA.loc[pd.IndexSlice[:, 'PDC'],:]
#         tmp.to_excel(writer, sheet_name=f'PDC VEI {i} 50%')
    
#     for i in [3,4,5]:
#         DATA = getData(VEI=i, agg=agg)
#         tmp = DATA.loc[pd.IndexSlice[:, 'LC'],:]
#         tmp.to_excel(writer, sheet_name=f'Large Clast {i} 50%')
        
#     for i in [450000, 9800000]:
#         DATA = getData(VOL=i, agg=agg)
#         tmp = DATA.loc[pd.IndexSlice[:, 'BAF'],:]
#         tmp.to_excel(writer, sheet_name=f'BAF {i} m3 990 m50%')
        
#     for i in [10,30,100]:
#         DATA = getData(RADIUS=i, agg=agg)
#         tmp = DATA.loc[pd.IndexSlice[:, 'BUF'],:]
#         tmp.to_excel(writer, sheet_name=f'Radius {i} km')
#%%
# Same data all on one sheet
with pd.ExcelWriter('csv/PerRegion_v2.xlsx') as writer: 
    startcol = 0
    startrow = 0
    for i in [3,4,5]:
        DATA = getData(VEI=i, agg=agg)
        tmp = DATA.loc[pd.IndexSlice[:, 'Tephra'],:]
        tmp.index.rename(f'VEI {i} 50%',level=0,inplace=True)
        tmp.to_excel(writer, startrow=startrow, startcol=startcol, sheet_name='values')     
        tmp.rank(ascending=False, axis=0).to_excel(writer, startrow=startrow, startcol=startcol, sheet_name='rank')     
        startrow+=8

    startrow = 0
    startcol+=8
    for i in [3,4,5]:
        DATA = getData(VEI=i, agg=agg)
        tmp = DATA.loc[pd.IndexSlice[:, 'PDC'],:]
        tmp.index.rename(f'VEI {i} 50%',level=0,inplace=True)
        tmp.to_excel(writer, startrow=startrow, startcol=startcol, sheet_name='values')   
        tmp.rank(ascending=False, axis=0).to_excel(writer, startrow=startrow, startcol=startcol, sheet_name='rank')   
        startrow+=8

    startrow = 0
    startcol+=8

    for i in [3,4,5]:
        DATA = getData(VEI=i, agg=agg)
        tmp = DATA.loc[pd.IndexSlice[:, 'LC'],:]
        tmp.index.rename(f'VEI {i} 50%',level=0,inplace=True)
        tmp.to_excel(writer, startrow=startrow, startcol=startcol, sheet_name='values')   
        tmp.rank(ascending=False, axis=0).to_excel(writer, startrow=startrow, startcol=startcol, sheet_name='rank')   
        startrow+=8

    startrow = 0        
    startcol+=8
    for i in [450000, 9800000]:
        DATA = getData(VOL=i, agg=agg)
        tmp = DATA.loc[pd.IndexSlice[:, 'BAF'],:]
        tmp.index.rename(f'{i} m3 990 m50%',level=0,inplace=True)
        tmp.to_excel(writer, startrow=startrow, startcol=startcol, sheet_name='values')   
        tmp.rank(ascending=False, axis=0).to_excel(writer, startrow=startrow, startcol=startcol, sheet_name='rank')   
        startrow+=8

    startrow = 0       
    startcol+=8  
    for i in [10,30,100]:
        DATA = getData(RADIUS=i, agg=agg)
        tmp = DATA.loc[pd.IndexSlice[:, 'BUF'],:]
        tmp.index.rename(f'Radius {i} km',level=0,inplace=True)
        tmp.to_excel(writer, startrow=startrow, startcol=startcol, sheet_name='values')   
        tmp.rank(ascending=False, axis=0).to_excel(writer, startrow=startrow, startcol=startcol, sheet_name='rank')   
        
        startrow+=8
        
# VEI = pd.concat([VEI3, VEI4, VEI5], keys=['VEI3', 'VEI4', 'VEI5'],axis=1)
# VOL = pd.concat([VOL1, VOL2], keys=['Vol. 1', 'Vol. 2'],axis=1)
#%%



# DATA.style.format({
#     'pop_count': "{:.2E}",
#     'area_urban': "{:.2f}",
#     'area_crops': "{:.2f}",
#     'length_roads': "{:.2E}",
#     'n_buildings': "{:.2f}",
#     })

# DATA.to_excel('csv/MeanExposurePerRegion.xlsx')

#%%
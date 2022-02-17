#%% [markdown]
# This file combines MASTER_exposire and FM_SEA into tables for the paper. One
# is conditional probabilities (with 1 col per VEI), one is absolute with yearly
# probabilities (1 col for the 10, 50 and 90th percentiles of the FM relationship)
# This is the script used to make the ranking tables.

#%%
import pandas as pd

# data = pd.read_csv('csv/MASTER_exposure_v3.csv')
# # data = data[data.hazard!='BAF']
# data = data[data.month==0]
# # data = data.drop(['buffer', 'volume','month', 'Unnamed: 0'], axis=1)

# data['volcano'] = data.volcano.str.replace('_', ' ')



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
# data = data.set_index(['hazard', 'VEI', 'prob'])

# Format string
data['volcano'] = data['volcano'].str.replace('_',' ')
# data['volcano'] = data['volcano'].str.replace('Leroboleng','Lereboleng')
data['region'] = data['region'].str.replace('Sumatera','Sumatra')
data['region'] = data['region'].str.replace('The Philippines','Philippines')
data['region'] = data['region'].str.replace('Maluku/Halmahera','Halmahera/Banda Sea')
data = data.sort_values(by=['country', 'region', 'volcano'])

# Sum roads and buildings numbers for exposure analysis and drop RNDS and cost
data['length_roads'] = data[['Motorway','Arterial','Collector','Local']].sum(axis=1)
data['n_buildings'] = data[['buildingsW', 'buildingsS']].sum(axis=1)
data = data.drop(['Motorway','Arterial','Collector','Local','buildingsW', 'buildingsS','RNDS', 'buildingsLoss'], axis=1)






FM = pd.read_csv('csv/Master_FM_SEA_Exposure.csv')
FM.columns = map(lambda x: str.replace(x, "Volcano Name", "volcano"), FM.columns)

data = data.set_index('volcano')
FM = FM.set_index('volcano')
FM = FM.drop(['Volcano Number', 'GVP DB Year', 'Included eruptions',
       'Estimate method', 'Change point', 'Average method', 'Power law',
       'Probability of eruption 10th percentile',
       'Probability of eruption 50th percentile',
       'Probability of eruption 90th percentile', 'VEI <=2 10th percentile',
       'VEI <=2 50th percentile', 'VEI <=2 90th percentile', 'VEI 6 10th percentile',
       'VEI 6 50th percentile', 'VEI 6 90th percentile',
       'VEI 7 10th percentile', 'VEI 7 50th percentile',
       'VEI 7 90th percentile'],axis=1)
Data = data.join(FM)
Data = Data.sort_values(by=['country', 'region', 'volcano'])
#%%
DATA = Data.set_index(['hazard', 'prob', 'VEI', 'mass'], append=True)
# DATA = DATA.loc[:,:,50,:,:][['pop_count','area_crops','area_urban', 'RNDS','buildingsLoss']].stack()
DATA = DATA.loc[:,:,50,:,:][['pop_count','area_crops','area_urban', 'length_roads','n_buildings']].stack()

BAF = Data.set_index(['hazard', 'prob', 'buffer', 'volume'], append=True).stack()
BUF = Data.set_index(['hazard', 'radius'], append=True).stack()
# LC = LC.loc[:,'BAF',50, 990, 9800000]
# LC = LC.loc[:,'BAF',50, 990, 9800000]

# Used for RNDS and buildings costs
# con = pd.DataFrame(DATA.loc[:,'PDC', 3, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to PDC'], [3]], names=('', 'VEI'),codes=[[0],[0]]))\
#        .join(pd.DataFrame(DATA.loc[:,'PDC', 4, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to PDC'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'PDC', 5, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to PDC'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(LC.loc[:,'BAF',50, 990, 450000,'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to BAF'], [0]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'LC', 3, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to large clasts'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'LC', 4, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to large clasts'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'LC', 5, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to large clasts'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 3, 1, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to 1 kg/m2'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, 1, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to 1 kg/m2'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, 1, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to 1 kg/m2'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 3, 5, 'area_crops'], columns=pd.MultiIndex(levels=[['Crops exposed to 5 kg/m2'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, 5, 'area_crops'], columns=pd.MultiIndex(levels=[['Crops exposed to 5 kg/m2'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, 5, 'area_crops'], columns=pd.MultiIndex(levels=[['Crops exposed to 5 kg/m2'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 3, None, 'buildingsLoss'], columns=pd.MultiIndex(levels=[['Tephra fall building repair costs'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, None, 'buildingsLoss'], columns=pd.MultiIndex(levels=[['Tephra fall building repair costs'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, None, 'buildingsLoss'], columns=pd.MultiIndex(levels=[['Tephra fall building repair costs'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 3, None, 'RNDS'], columns=pd.MultiIndex(levels=[['Road Network Disruption Score'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, None, 'RNDS'], columns=pd.MultiIndex(levels=[['Road Network Disruption Score'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, None, 'RNDS'], columns=pd.MultiIndex(levels=[['Road Network Disruption Score'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\

# con = pd.DataFrame(DATA.loc[:,'PDC', 3, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to PDC'], [3]], names=('', 'VEI'),codes=[[0],[0]]))\
#        .join(pd.DataFrame(DATA.loc[:,'PDC', 4, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to PDC'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'PDC', 5, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to PDC'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(LC.loc[:,'BAF',50, 990, 450000,'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to BAF'], [0]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'LC', 3, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to large clasts'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'LC', 4, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to large clasts'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'LC', 5, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to large clasts'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 3, 1, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to 1 kg/m2'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, 1, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to 1 kg/m2'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, 1, 'pop_count'], columns=pd.MultiIndex(levels=[['Population exposed to 1 kg/m2'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 3, 5, 'area_crops'], columns=pd.MultiIndex(levels=[['Crops exposed to 5 kg/m2'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, 5, 'area_crops'], columns=pd.MultiIndex(levels=[['Crops exposed to 5 kg/m2'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
#        .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, 5, 'area_crops'], columns=pd.MultiIndex(levels=[['Crops exposed to 5 kg/m2'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       # .join(pd.DataFrame(DATA.loc[:,'Tephra', 3, None, 'buildingsLoss'], columns=pd.MultiIndex(levels=[['Tephra fall building repair costs'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       # .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, None, 'buildingsLoss'], columns=pd.MultiIndex(levels=[['Tephra fall building repair costs'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       # .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, None, 'buildingsLoss'], columns=pd.MultiIndex(levels=[['Tephra fall building repair costs'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       # .join(pd.DataFrame(DATA.loc[:,'Tephra', 3, None, 'RNDS'], columns=pd.MultiIndex(levels=[['Road Network Disruption Score'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       # .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, None, 'RNDS'], columns=pd.MultiIndex(levels=[['Road Network Disruption Score'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       # .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, None, 'RNDS'], columns=pd.MultiIndex(levels=[['Road Network Disruption Score'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\

TephraCon = pd.DataFrame(DATA.loc[:,'Tephra', 3, 1, 'pop_count'], columns=pd.MultiIndex(levels=[['Population (n) exposed to 1 kg/m2'], [3]], names=('', 'VEI'),codes=[[0],[0]]))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, 1, 'pop_count'], columns=pd.MultiIndex(levels=[['Population (n) exposed to 1 kg/m2'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, 1, 'pop_count'], columns=pd.MultiIndex(levels=[['Population (n) exposed to 1 kg/m2'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 3, 100, 'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n) exposed to 100 kg/m2'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, 100, 'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n) exposed to 100 kg/m2'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, 100, 'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n) exposed to 100 kg/m2'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 3, 1, 'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km) exposed to 1 kg/m2'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, 1, 'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km) exposed to 1 kg/m2'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, 1, 'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km) exposed to 1 kg/m2'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 3, 5, 'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2) exposed to 5 kg/m2'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, 5, 'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2) exposed to 5 kg/m2'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, 5, 'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2) exposed to 5 kg/m2'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 3, 1, 'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2) exposed to 1 kg/m2'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 4, 1, 'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2) exposed to 1 kg/m2'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'Tephra', 5, 1, 'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2) exposed to 1 kg/m2'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .round(2)
PDCCon =  pd.DataFrame(DATA.loc[:,'PDC', 3,  None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population (n)'], [3]], names=('', 'VEI'),codes=[[0],[0]]))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 4, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population (n)'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 5, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population (n)'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 3, None, 'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km)'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 4, None, 'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km)'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 5, None, 'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km)'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 3, None, 'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n)'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 4, None, 'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n)'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 5, None, 'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n)'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 3, None, 'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2)'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 4, None, 'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2)'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 5, None, 'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2)'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 3, None, 'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2)'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 4, None, 'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2)'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'PDC', 5, None, 'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2)'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .round(2)
LCCon =  pd.DataFrame(DATA.loc[:,'LC', 3,  None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population (n)'], [3]], names=('', 'VEI'),codes=[[0],[0]]))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 4, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population (n)'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 5, None, 'pop_count'], columns=pd.MultiIndex(levels=[['Population (n)'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 3, None, 'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km)'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 4, None, 'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km)'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 5, None, 'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km)'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 3, None, 'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n)'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 4, None, 'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n)'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 5, None, 'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n)'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 3, None, 'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2)'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 4, None, 'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2)'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 5, None, 'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2)'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 3, None, 'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2)'], [3]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 4, None, 'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2)'], [4]], names=('', 'VEI'),codes=[[0],[0]])))\
       .join(pd.DataFrame(DATA.loc[:,'LC', 5, None, 'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2)'], [5]], names=('', 'VEI'),codes=[[0],[0]])))\
       .round(2)
BAFCon =  pd.DataFrame(    BAF.loc[:,'BAF',50, 990,450000, 'pop_count'], columns=pd.MultiIndex(levels=[['Population (n)'], [450000]], names=('', 'Flow volume (m3)'),codes=[[0],[0]]))\
       .join(pd.DataFrame(BAF.loc[:,'BAF',50, 990, 9800000, 'pop_count'], columns=pd.MultiIndex(levels=[['Population (n)'], [9800000]], names=('', 'Flow volume (m3)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BAF.loc[:,'BAF',50, 990, 450000, 'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km)'], [450000]], names=('', 'Flow volume (m3)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BAF.loc[:,'BAF',50, 990, 9800000, 'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km)'], [9800000]], names=('', 'Flow volume (m3)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BAF.loc[:,'BAF',50, 990, 450000, 'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n)'], [450000]], names=('', 'Flow volume (m3)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BAF.loc[:,'BAF',50, 990, 9800000, 'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n)'], [9800000]], names=('', 'Flow volume (m3)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BAF.loc[:,'BAF',50, 990, 450000, 'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2)'], [450000]], names=('', 'Flow volume (m3)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BAF.loc[:,'BAF',50, 990, 9800000, 'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2)'], [9800000]], names=('', 'Flow volume (m3)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BAF.loc[:,'BAF',50, 990, 450000, 'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2)'], [450000]], names=('', 'Flow volume (m3)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BAF.loc[:,'BAF',50, 990, 9800000, 'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2)'], [9800000]], names=('', 'Flow volume (m3)'),codes=[[0],[0]])))\
       .round(2)
BUFCon =  pd.DataFrame(   BUF.loc[:,'BUF',10,  'pop_count'], columns=pd.MultiIndex(levels=[['Population (n)'], [10]], names=('', 'Radius (km)'),codes=[[0],[0]]))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',30,  'pop_count'], columns=pd.MultiIndex(levels=[['Population (n)'], [30]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',100, 'pop_count'], columns=pd.MultiIndex(levels=[['Population (n)'], [100]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',10,  'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km)'], [10]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',30,  'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km)'], [30]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',100, 'length_roads'], columns=pd.MultiIndex(levels=[['Road length (km)'], [100]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',10,  'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n)'], [10]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',30,  'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n)'], [30]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',100, 'n_buildings'], columns=pd.MultiIndex(levels=[['Buildings (n)'], [100]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',10,  'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2)'], [10]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',30,  'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2)'], [30]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',100, 'area_crops'], columns=pd.MultiIndex(levels=[['Crop areas (km2)'], [100]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',10,  'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2)'], [10]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',30,  'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2)'], [30]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .join(pd.DataFrame(BUF.loc[:,'BUF',100, 'area_urban'], columns=pd.MultiIndex(levels=[['Urban areas (km2)'], [100]], names=('', 'Radius (km)'),codes=[[0],[0]])))\
       .round(2)

def rankMe(df):
       dfR = df.rank(ascending=False)
       dfR['sum'] = dfR.sum(axis=1)
       dfR['sum'] = dfR['sum'].rank()
       dfR = dfR.sort_values('sum')
       dfR = dfR.drop('sum',axis=1,level=0)
       return df.loc[dfR.index], dfR

TephraCon, TephraConR = rankMe(TephraCon)
PDCCon, PDCConR = rankMe(PDCCon)      
LCCon, LCConR = rankMe(LCCon)      
BAFCon, BAFConR = rankMe(BAFCon)      
BUFCon, BUFConR = rankMe(BUFCon)      


#%% 
# columns = ['Population exposed to PDC','Population exposed to large clasts','Population exposed to 1 kg/m2','Crops exposed to 5 kg/m2','Tephra fall building repair costs', 'Road Network Disruption Score']
# columns = ['Population exposed to PDC','Population exposed to large clasts','Population exposed to 1 kg/m2','Crops exposed to 5 kg/m2','Roads exposed to 1 kg/m2', 'Buildings exposed to 1 kg/m2']
def annualProb(con):
       # con = TephraCon
       columns = list(con.columns.levels[0])
       prob = [10,50,90]
       VEI = [3,4,5]
       # Data['VEI'] = Data['VEI'].astype(int)
       # Data = Data.set_index(['VEI'], append=True)
       cnt = 0
       for iC in columns:
              for iP in prob:
                     tmp = con[iC][3].multiply(Data[f'VEI {3} {iP}th percentile'].drop_duplicates()) + \
                     con[iC][4].multiply(Data[f'VEI {4} {iP}th percentile'].drop_duplicates()) + \
                     con[iC][5].multiply(Data[f'VEI {5} {iP}th percentile'].drop_duplicates())
                     
                     if cnt == 0:
                            abs = pd.DataFrame(tmp, columns=pd.MultiIndex(levels=[[iC], [iP]], names=('', 'F-M percentile'),codes=[[0],[0]]))
                     else:
                            abs = abs.join(pd.DataFrame(tmp, columns=pd.MultiIndex(levels=[[iC], [iP]], names=('', 'F-M percentile'),codes=[[0],[0]])))
                     cnt += 1

       abs = abs.round(2)
       abs, absR = rankMe(abs)
       
       return abs, absR

TephraAbs, TephraAbsR = annualProb(TephraCon)
PDCAbs, PDCAbsR = annualProb(PDCCon)
LCAbs, LCAbsR = annualProb(LCCon)



#%%
with pd.ExcelWriter('csv/ExposureFormatted.xlsx') as writer: 
       TephraCon.to_excel(writer, sheet_name='Tephra conditional')
       TephraConR.to_excel(writer, sheet_name='Tephra conditional rank')
       PDCCon.to_excel(writer, sheet_name='PDC conditional')
       PDCConR.to_excel(writer, sheet_name='PDC conditional rank')
       LCCon.to_excel(writer, sheet_name='Large clasts conditional')
       LCConR.to_excel(writer, sheet_name='Large clasts conditional rank')
       BAFCon.to_excel(writer, sheet_name='BAF conditional')
       BAFConR.to_excel(writer, sheet_name='BAF conditional rank')
       BUFCon.to_excel(writer, sheet_name='Radius conditional')
       BUFConR.to_excel(writer, sheet_name='Radius conditional rank')
       TephraAbs.to_excel(writer, sheet_name='Tephra absolute')
       TephraAbsR.to_excel(writer, sheet_name='Tephra absolute rank')
       PDCAbs.to_excel(writer, sheet_name='PDC absolute')
       PDCAbsR.to_excel(writer, sheet_name='PDC absolute rank')
       LCAbs.to_excel(writer, sheet_name='Large clasts absolute')
       LCAbsR.to_excel(writer, sheet_name='Large clasts absolute rank')
# %%

TephraAbsTmp = TephraAbsR.drop([10,90],axis=1,level=1).rename({50:'A'},axis=1)
TephraF = TephraAbsTmp.join(TephraConR).sort_index(axis=1)
TephraF.columns.rename(['',''],inplace=True)

PDCAbsTmp = PDCAbsR.drop([10,90],axis=1,level=1).rename({50:'A'},axis=1)
PDCF = PDCAbsTmp.join(PDCConR).sort_index(axis=1)
PDCF.columns.rename(['',''],inplace=True)

LCAbsTmp = LCAbsR.drop([10,90],axis=1,level=1).rename({50:'A'},axis=1)
LCF = LCAbsTmp.join(LCConR).sort_index(axis=1)
LCF.columns.rename(['',''],inplace=True)
# %%

#%%
# Formats the table for appendix 3 of the paper
import pandas as pd


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
# data['volcano'] = data['volcano'].str.replace('Leroboleng','Lereboleng')
data['region'] = data['region'].str.replace('Sumatera','Sumatra')
data['region'] = data['region'].str.replace('The Philippines','Philippines')
data['region'] = data['region'].str.replace('Maluku/Halmahera','Halmahera/Banda Sea')
data['hazard'] = data['hazard'].str.replace('LC','Large clast')
data['hazard'] = data['hazard'].str.replace('BUF','Circular buffers')
data = data.sort_values(by=['country', 'region', 'volcano'])
data = data[~((data['hazard']=='Tephra') & (data['mass'].isnull()))] # Remove the row

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
Data = Data.drop(['lat','lon','xmin','xmax','ymin','ymax','RNDS','buildingsLoss'],axis=1)
Data = Data.reset_index().set_index(['country', 'region','volcano','hazard'])


# %%

col = {
    'VEI': 'VEI',
    'prob': 'Hazard probability (%)',
    'mass': 'Tephra load (kg/m2)',
    'buffer': 'BAF buffer (m)',
    'volume': 'BAF volume (m3)',
    'radius': 'Radius size (km)',
    'pop_count': 'Population count',
    'n_buildings': 'Buildings count',
    'area_crops': 'Crop area (km2)',
    'area_urban': 'Urban area (km2)',
    'Motorway': 'Motorway roads length (km)',
    'Arterial': 'Arterial roads length (km)',
    'Collector': 'Collector roads length (km)',
    'Local': 'Local roads length (km)',
    'VEI 3 10th percentile': 'VEI 3 10th pct',
    'VEI 3 50th percentile': 'VEI 3 50th pct',
    'VEI 3 90th percentile': 'VEI 3 90th pct',
    'VEI 4 10th percentile': 'VEI 4 10th pct',
    'VEI 4 50th percentile': 'VEI 4 50th pct',
    'VEI 4 90th percentile': 'VEI 4 90th pct',
    'VEI 5 10th percentile': 'VEI 5 10th pct',
    'VEI 5 50th percentile': 'VEI 5 50th pct',
    'VEI 5 90th percentile': 'VEI 5 90th pct',
}
colSup = {
    'VEI': 'Hazard properties',
    'prob': 'Hazard properties',
    'mass': 'Hazard properties',
    'buffer': 'Hazard properties',
    'volume': 'Hazard properties',
    'radius': 'Hazard properties',
    'pop_count': 'Exposure properties',
    'n_buildings': 'Exposure properties',
    'area_crops': 'Exposure properties',
    'area_urban': 'Exposure properties',
    'Motorway': 'Exposure properties',
    'Arterial': 'Exposure properties',
    'Collector': 'Exposure properties',
    'Local': 'Exposure properties',
    'VEI 3 10th percentile': 'Annual probability',
    'VEI 3 50th percentile': 'Annual probability',
    'VEI 3 90th percentile': 'Annual probability',
    'VEI 4 10th percentile': 'Annual probability',
    'VEI 4 50th percentile': 'Annual probability',
    'VEI 4 90th percentile': 'Annual probability',
    'VEI 5 10th percentile': 'Annual probability',
    'VEI 5 50th percentile': 'Annual probability',
    'VEI 5 90th percentile': 'Annual probability',
}
arrays = [list(colSup.values()), list(col.values())]
tuples = list(zip(*arrays))
idx = pd.MultiIndex.from_tuples(tuples,  names=['',''])

Data = Data[col].rename(col, axis=1)
Data.columns=idx

Data.to_excel('csv/Appendix3.xlsx')

# %%

# %%

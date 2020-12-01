#%%
import os
os.chdir('/Users/seb/Documents/Codes/VolcGIS')
import volcgis 
from rasterio.plot import show
import rasterio as rio


#%% SETUP ERUPTION
# Setup eruption dictionary. This will have to be updated from main DB
eruption = {
    'name':     'Gede-Pangrango',
    'vent':     [-6.787948, 106.981215, 2765], # Vent lat, lon and elevation
    'area':     [-7, -6.55, 106.7, 107.2], # S N W E
    'extent':   [1.e5, 1.e5, 1.e5, 1.e5], # extent around vent in meters, [minx maxx miny maxy]
    'epsg':     32748
}

# Setup eruption
erup = volcgis.eruption(eruption, 90)


#%% SETUP TEPHRA HAZARD

nameConstructor = {
    'volcano':    ['Gede-Pangrango'],
    'VEI':        ['VEI3', 'VEI4', 'VEI5'],
    'perc':       ['P1', 'P5', 'P9'],
    'format':     ['.tif']
}

tephra = {
    'hazard':   'Tephra',
    'epsg':     32748,
    'kind':     'him',
    'unit':     'kg/m2',
    'label':    [''],
    'rootDir':   'hazards/Tephra/',
    'nameConstructor': nameConstructor,
    }

erup.prepareHazard(tephra)

#%% SETUP BLOCK AND ASH FLOW
nameConstructor = {
    'volcano':    ['Gede-Pangrango'],
    'buffer':     ['0.0027deg', '0.009deg'],
    'volume':     ['9800000m3', '450000m3'],
    'format':     ['.tif']
}

BAF = {
    'hazard':   'BAF',
    'epsg':     32748,
    'kind':     'prob',
    'unit':     '%',
    'label':    [''],
    'rootDir':    'hazards/BAF/',
    'nameConstructor': nameConstructor,
    }  

erup.prepareHazard(BAF)



#%% SETUP COLUMN COLLAPSE PDC
nameConstructor = {
    'volcano':    ['Gede-Pangrango'],
    'VEI':        ['3', '4', '5'],
    'suffix':     ['output_map'],
    'format':     ['.asc']
}

PDC = {
    'hazard':   'PDC',
    'epsg':     32748,
    'kind':     'prob',
    'unit':     '%',
    'rootDir':    'hazards/PDC/',
    'nameConstructor': nameConstructor,
    }  

erup.prepareHazard(PDC)

#%% SETUP Large clast
nameConstructor = {
    'volcano':    ['Gede-Pangrango'],
    'VEI':        ['VEI3', 'VEI4', 'VEI5'],
    'format':     ['.tif']
}

LC = {
    'hazard':   'LC',
    'epsg':     4326,
    'kind':     'prob',
    'unit':     '%',
    'rootDir':  'hazards/LC/',
    'nameConstructor': nameConstructor,
    }  

erup.prepareHazard(LC)

#%% PREPARE EXPOSURE DATA
# erup.getLandcover()
erup = volcgis.eruption(eruption, 90)
erup.getLandscan()


# %%

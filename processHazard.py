#%%
import os
os.chdir('/Users/seb/Documents/Codes/VolcGIS')
import volcgis 
from rasterio.plot import show
import rasterio as rio
import pandas as pd
from pyproj import CRS
from pyproj import Transformer
import numpy as np
import utm

# Tephra thresholds
# Josh: 1, 100
# Ag/Hort: 1, 5, 50
# Forest: 5, 200, 1000

#%%
volcDB = pd.read_csv('volcCoordinates.csv')
res = 90

# Steps selection
processTephra = True
processBAF = True
processPDC = True
processLC = True

analyzeTephra = True
analyzeBAF = True
analyzePDC = True
analyzeLC = True

# Main loops for processing
VEI = [3,4,5]
probT = [10,50,90]
intT = [1,5,50,100]
buffT = ['0.009', '0.0027']
volT = ['9800000', '450000']

# Exposure storage
exposure = {}
exposure['volcano'] = []
exposure['hazard'] = []
exposure['VEI'] = []
exposure['prob'] = []
exposure['mass'] = []
exposure['buffer'] = []
exposure['volume'] = []
exposure['pop_count'] = []
exposure['area_crops'] = []
exposure['area_urban'] = []

EXPOSURE = pd.DataFrame(exposure)
 
 
def updateExposure(EXPOSURE, volcano, hazard, VEI, prob, intensity, buffer, volume, pop, crops, urban):
    
    expTmp = {
        'volcano':      [volcano],
        'hazard':       [hazard],
        'VEI':          [VEI],
        'prob':         [prob],
        'mass':         [intensity],
        'pop_count':    [pop],
        'area_crops':   [crops],
        'area_urban':   [urban],
        'buffer':       [buffer],
        'volume':       [volume],
    }
    
    expTmp = pd.DataFrame(expTmp)
    return EXPOSURE.append(expTmp)
    
    
# def updateExposure(exposure, volcano, hazard, VEI, prob, intensity, buffer, volume, pop, crops, urban):
#     exposure['volcano'].append(volcano)
#     exposure['hazard'].append(hazard)
#     exposure['VEI'].append(VEI)
#     exposure['prob'].append(prob)
#     exposure['mass'].append(intensity)
#     exposure['pop_count'].append(pop)
#     exposure['area_crops'].append(crops)
#     exposure['area_urban'].append(urban)
#     exposure['buffer'].append(buffer)
#     exposure['volume'].append(volume)
#     return exposure

def getExposure(haz_data, pop_data, LC_data, val):
    # Get hazard index
    idx = haz_data >= val
    # Population
    popTmp = np.sum(pop_data[idx])
    
    # Landcover / crops km2
    idx = (haz_data >= val) & (LC_data==40)
    cropsTmp = np.sum(idx)*res**2/1e6
    # Landcover / urban km2
    idx = (haz_data >= val) & (LC_data==50)
    urbanTmp = np.sum(idx)*res**2/1e6
    
    return popTmp, cropsTmp, urbanTmp
                        
                        
for i in range(0, volcDB.shape[0]):
    print(volcDB.loc[i, 'name'].upper())
    lat = volcDB.loc[i, 'lat']
    lon = volcDB.loc[i, 'lon']
    name = volcDB.loc[i, 'name']
    
    # Find the espg
    zone = utm.from_latlon(lat,lon)  
    if lat<0:
        crs = CRS.from_dict({'proj': 'utm', 'zone': zone[2], 'south': True})
    else:
        crs = CRS.from_dict({'proj': 'utm', 'zone': zone[2], 'south': False})
    epsg = crs.to_authority()
    epsg = int(epsg[1])
    
    # Get the utm coordinates of the vent
    transformer = Transformer.from_crs(4326, epsg)
    [xtmp, ytmp] = transformer.transform(lat, lon)
    
    # Find region extent
    if volcDB.loc[i, 'xmin'] < 0 :
        dxMin = round(xtmp + abs(volcDB.loc[i, 'xmin']))
    else:
        dxMin = round(xtmp - volcDB.loc[i, 'xmin'])
    dxMax = round(volcDB.loc[i, 'xmax'] - xtmp)
    
    if volcDB.loc[i, 'ymin'] < 0 :
        dyMin = round(ytmp + abs(volcDB.loc[i, 'ymin']))
    else:
        dyMin = round(ytmp - volcDB.loc[i, 'ymin'])
    dyMax = round(volcDB.loc[i, 'ymax'] - ytmp)
        
    # Setup eruption dictionary. 
    eruption = {
        'name':     name,
        'vent':     [lat, lon, 2765], # Vent lat, lon and elevation
        'extent':   [dxMin, dxMax, dyMin, dyMax], # extent around vent in meters, [minx maxx miny maxy]
        'epsg':     epsg
    }

    erup = volcgis.eruption(eruption, res)


    #%% SETUP TEPHRA HAZARD     
    if processTephra:
        nameConstructor = {
            'volcano':    [name],
            'VEI':        ['VEI3', 'VEI4', 'VEI5'],
            'perc':       ['P1', 'P5', 'P9'],
            'format':     ['.tif']
        }

        tephra = {
            'hazard':   'Tephra',
            'epsg':     epsg,
            'rootDir':   'hazards/Tephra/',
            'nameConstructor': nameConstructor,
            }

        erup.prepareHazard(tephra)

    # SETUP BLOCK AND ASH FLOW
    if processBAF:
        nameConstructor = {
            'volcano':    [name],
            'buffer':     ['0.0027deg', '0.009deg'],
            'volume':     ['9800000m3', '450000m3'],
            'format':     ['.tif']
        }

        BAF = {
            'hazard':   'BAF',
            'epsg':     epsg,
            'kind':     'prob',
            'unit':     '%',
            'label':    [''],
            'rootDir':    'hazards/BAF/',
            'probT':    [.1,.5,.9],
            'intensityT': None,
            'nameConstructor': nameConstructor,
            }  

        erup.prepareHazard(BAF)

    # SETUP COLUMN COLLAPSE PDC
    if processPDC:
        nameConstructor = {
            'volcano':    [name],
            'VEI':        ['3', '4', '5'],
            'suffix':     ['output_map'],
            'format':     ['.asc']
        }

        PDC = {
            'hazard':   'PDC',
            'epsg':     epsg,
            'kind':     'prob',
            'unit':     '%',
            'rootDir':    'hazards/PDC/',
            'probT':    [.1,.5,.9],
            'intensityT': None,
            'nameConstructor': nameConstructor,
            }  

        erup.prepareHazard(PDC)

    # SETUP Large clast
    if processLC:
        nameConstructor = {
            'volcano':    [name],
            'VEI':        ['VEI3', 'VEI4', 'VEI5'],
            'format':     ['.tif']
        }

        LC = {
            'hazard':   'LC',
            'epsg':     4326,
            'kind':     'prob',
            'unit':     '%',
            'rootDir':  'hazards/LC/',
            'probT':    [.1,.5,.9],
            'intensityT': None,
            'nameConstructor': nameConstructor,
            }  

        erup.prepareHazard(LC)

    # PREPARE EXPOSURE DATA
    erup.getLandcover()
    erup.getLandscan()
    
    
    # Exposure analysis
    LCf = 'volcanoes/{}/_data/Landcover.tif'.format(erup.name)
    popf = 'volcanoes/{}/_data/Landscan.tif'.format(erup.name)
    
    LC = rio.open(LCf)
    pop = rio.open(popf)
    
    LC_data = LC.read(1)
    pop_data = pop.read(1)/(1000/res)**2
    pop_data[pop_data<1] = 0
    
    if analyzeTephra:
        for iVEI in VEI:
            for iP in probT:
                fl = 'volcanoes/{}/_hazard/Tephra/{}_VEI{}_P{}.tif'.format(erup.name, erup.name, iVEI, int(10-iP/10))
                with rio.open(fl) as haz:
                    haz_data = haz.read(1)
                    for iM in intT:
                        popTmp, cropsTmp, urbanTmp = getExposure(haz_data, pop_data, LC_data, iM)
                        # exposure = updateExposure(exposure, erup.name, 'Tephra', iVEI, iP, iM, None, None, popTmp, cropsTmp, urbanTmp)
                        EXPOSURE = updateExposure(EXPOSURE, erup.name, 'Tephra', iVEI, iP, iM, None, None, popTmp, cropsTmp, urbanTmp)
    
    if analyzeBAF: 
        for iV in volT:
            for iB in buffT:
                fl = 'volcanoes/{}/_hazard/BAF/{}_{}deg_{}m3.tif'.format(erup.name, erup.name, iB, iV)
                with rio.open(fl) as haz:
                    haz_data = haz.read(1)
                    for iP in probT:
                        popTmp, cropsTmp, urbanTmp = getExposure(haz_data, pop_data, LC_data, iP/100)
                        # exposure = updateExposure(exposure, erup.name, 'BAF', None, iP, None, iB, iV, popTmp, cropsTmp, urbanTmp)
                        EXPOSURE = updateExposure(EXPOSURE, erup.name, 'BAF', None, iP, None, iB, iV, popTmp, cropsTmp, urbanTmp)

    if analyzePDC: 
        for iVEI in VEI:
            fl = 'volcanoes/{}/_hazard/PDC/{}_{}_output_map.tif'.format(erup.name, erup.name, iVEI)
            with rio.open(fl) as haz:
                haz_data = haz.read(1)
                for iP in probT:
                    popTmp, cropsTmp, urbanTmp = getExposure(haz_data, pop_data, LC_data, iP/100)
                    # # exposure = updateExposure(exposure, erup.name, 'PDC', iVEI, iP, None, None, None, popTmp, cropsTmp, urbanTmp)
                    EXPOSURE = updateExposure(EXPOSURE, erup.name, 'PDC', iVEI, iP, None, None, None, popTmp, cropsTmp, urbanTmp)
    
    if analyzeLC: 
        for iVEI in VEI:
            fl = 'volcanoes/{}/_hazard/LC/{}_VEI{}.tif'.format(erup.name, erup.name, iVEI)
            with rio.open(fl) as haz:
                haz_data = haz.read(1)
                for iP in probT:
                    popTmp, cropsTmp, urbanTmp = getExposure(haz_data, pop_data, LC_data, iP/100)
                    # exposure = updateExposure(exposure, erup.name, 'LC', iVEI, iP, None, None, None, popTmp, cropsTmp, urbanTmp)
                    EXPOSURE = updateExposure(EXPOSURE, erup.name, 'LC', iVEI, iP, None, None, None, popTmp, cropsTmp, urbanTmp)
 
    LC.close()
    pop.close()

EXPOSURE.to_csv('results.csv')
# %%
# EXPOSURE = pd.DataFrame(exposure)




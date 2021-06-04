#%%
import os
import sys
# os.chdir('/Users/seb/Documents/Codes/VolcGIS')
import volcgis 
from volcgis.exposureAnalysis import *
from rasterio.plot import show
import rasterio as rio
import pandas as pd
from pyproj import CRS
from pyproj import Transformer
import numpy as np
from printy import printy
import geopandas as gpd
from shapely.geometry import shape
# from alive_progress import alive_bar
import time

# Tephra thresholds
# Josh: 1, 100
# Ag/Hort: 1, 5, 50
# Forest: 5, 200, 1000

#%% 

# Read the master csv file containing all the volcanoes
# The file contains these columns: name,lat,lon,xmin,xmax,ymin,ymax,country
volcDB = pd.read_csv('csv/volcCoordinates.csv')

# Read the vulnerability data
vuln = pd.read_csv('csv/function_params.csv')
vuln.columns = ['Country', 'Building_type', 'load_mean', 'load_disp', 'tot_rep_cost']
vuln['load_mean'] = vuln['load_mean']*101.97 # kPa to kgm-2 
vuln = vuln.set_index(['Country', 'Building_type'])

# The main resolution used across all rasters
res = 90

# Steps selection
processTephra = True
analyzeTephra = True
removeTephra = False

# Main loops for processing
VEI = [3,4,5]
probT = [10,50,90]
intT = [1,5,50,100]

# Dictionaries for road disruption
RNDS_intensity_map = {1: 10,
                    100: 100}

RNDS_probability_map = {0.1: 100,
                    0.5: 100,
                    0.9: 100}



#%% Main processing step            

# Handling if the function is called from the command line                          
if (len(sys.argv)==1) or (len(sys.argv) > 2):
    rng = range(0, volcDB.shape[0])
else:
    rng = range(int(sys.argv[1]), int(sys.argv[1])+1)

# Debug
# rng = range(0,1)

# Loop through volcanoes          
for i in rng:
    tic = time.perf_counter() # Start counter
    printy(volcDB.loc[i, 'name'].upper(), 'cB')
    lat = volcDB.loc[i, 'lat']
    lon = volcDB.loc[i, 'lon']
    name = volcDB.loc[i, 'name']
    
    # Find the espg
    zone = int(np.ceil((lon+180)/6))
    if lat<0:
        crs = CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': True})
    else:
        crs = CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': False})
    epsg = crs.to_authority()
    epsg = int(epsg[1])
    
    # Get the utm coordinates of the vent
    transformer = Transformer.from_crs(4326, epsg)
    [xtmp, ytmp] = transformer.transform(lat, lon)
    
    # Find region extent        
    dxMin, dxMax, dyMin, dyMax = volcgis.makeZoneBoundaries(xtmp, ytmp, 
                                                        volcDB.loc[i, 'xmin'], volcDB.loc[i, 'xmax'],
                                                        volcDB.loc[i, 'ymin'], volcDB.loc[i, 'ymax'])
    # Setup eruption dictionary. 
    eruption = {
        'name':     name,
        'vent':     [lat, lon, 2765], # Vent lat, lon and elevation
        'extent':   [volcDB.loc[i, 'xmin'], volcDB.loc[i, 'xmax'], volcDB.loc[i, 'ymin'], volcDB.loc[i, 'ymax']], # extent around vent in meters, [minx maxx miny maxy]
        'epsg':     epsg
    }

    # Define the eruption
    erup = volcgis.eruption(eruption, res)

    #%% SETUP TEPHRA HAZARD     
    if processTephra:
        nameConstructor = {
            'volcano':    [name],
            'VEI':        ['VEI3', 'VEI4', 'VEI5'],
            'perc':       ['P10', 'P50', 'P90'],
            'month':      ['month{}'.format(m) for m in range(1,13)],
            'format':     ['.tif']
        }

        tephra = {
            'hazard':   'Tephra',
            'epsg':     epsg,
            'rootDir':   'hazards/Tephra/',
            'nameConstructor': nameConstructor,
            }

        erup.prepareHazard(tephra)


    # PREPARE EXPOSURE DATA
    # if processLandCover:
    #     erup.getLandcover()
    # if processLandScan:
    #     erup.getLandscan()
    # if processRoadNetwork:
    #     erup.getRoadNetwork()
    # if processBuildings:
    #     erup.getBuildingExposure(inPath = 'DATA/buildings_strong_{}.tif'.format(volcDB.iloc[i].country),
    #                              outPath = os.path.join(erup.outPath, erup.name, '_data/buildings_strong.tif'))
    #     erup.getBuildingExposure(inPath = 'DATA/buildings_weak_{}.tif'.format(volcDB.iloc[i].country),
    #                              outPath = os.path.join(erup.outPath, erup.name, '_data/buildings_weak.tif'))
        
    
    # Exposure analysis
    LCf = os.path.join(erup.outPath, erup.name, '_data/Landcover.tif')
    popf = os.path.join(erup.outPath, erup.name, '_data/Landscan.tif')
    roadf = os.path.join(erup.outPath, erup.name, '_data/roads.feather')
    buildWf = os.path.join(erup.outPath, erup.name, '_data/buildings_weak.tif')
    buildSf = os.path.join(erup.outPath, erup.name, '_data/buildings_strong.tif')
    
    # If any exposure analysis is required, load the data
    if analyzeTephra:

        # Read datasets
        LC = rio.open(LCf)
        pop = rio.open(popf)
        road = gpd.read_feather(roadf)
        buildW = rio.open(buildWf)
        buildS = rio.open(buildSf)
        
        # Read and adjust
        LC_data = LC.read(1)
        pop_data = pop.read(1)/(1000/res)**2
        buildW_data = buildW.read(1)
        buildS_data = buildS.read(1)
        
        # Adjust residuals from interpolation
        pop_data[pop_data<1] = 0
        buildW_data[buildW_data<0] = 0
        buildS_data[buildS_data<0] = 0
    
    if analyzeTephra:
        printy(' - Processing tephra...')
        for iMth in range(1,13):
            # Exposure storage
            damageRatio = pd.DataFrame({})
            damageState = pd.DataFrame({})
            exposure = pd.DataFrame({})
            roadExposure = pd.DataFrame({})
            RSDS = pd.DataFrame()
            
            for iVEI in VEI:
                for iP in probT:
                    printy('   - VEI: {} - Prob: {}...'.format(iVEI, iP))
                    
                    # Define hazard file
                    fl = os.path.join(erup.outPath, erup.name, '_hazard/Tephra/{}_VEI{}_P{}_month{}.tif'.format(erup.name, iVEI, int(100-iP), iMth))
                    # Get the road disruption
                    rnds, roadLength, rsds = getRNDS(fl, RNDS_intensity_map, road, epsg, intensity=True)

                    # In case of tephra, need to concatenate the rsds for various columns/accumulations into one. Cat, rename and drop
                    
                    # Test to handle problems on Gekko
                    if 'RSDS_100' not in rsds.columns:
                        rsds['RSDS_100'] = 0
                    if 'RSDS_1' not in rsds.columns:
                        rsds['RSDS_1'] = 0
                    
                    rsds['tephra_VEI{}_P{}_month{}'.format(iVEI, int(100-iP), iMth)] = rsds.RSDS_100.fillna(0)+rsds.RSDS_1.fillna(0)
                    rsds = rsds[rsds.columns.drop(list(rsds.filter(regex='RSDS')))]
                    RSDS = updateRSDS(RSDS, rsds)
                    
                    # Compute exposure from Landscan and landcover
                    with rio.open(fl) as haz:
                        haz_data = haz.read(1)
                        
                        # Update exposure
                        # Loss from weak buildings
                        DRw, DSw = getBuildingImpact(haz_data, buildW_data, vuln.loc[(volcDB.loc[i,'country'], 'Weak')], erup.outPath, None, erup.name, haz.profile)
                        damageRatio, damageState = updateBuildingExposure(damageRatio,damageState,erup.name,iVEI, iP, 'weak', DRw, DSw)
                        
                        # Loss from strong buildings
                        DRs, DSs = getBuildingImpact(haz_data, buildW_data, vuln.loc[(volcDB.loc[i,'country'], 'Strong')], erup.outPath, None, erup.name, haz.profile)
                        damageRatio, damageState = updateBuildingExposure(damageRatio,damageState,erup.name,iVEI, iP, 'strong', DRs, DSs)
                        
                        # General exposure
                        exposure = updateExposure(exposure, erup.name, 'Tephra', iVEI, iP, None, None, None, None, None, None, rnds, round((DRs.loss.sum()+DRw.loss.sum())/1e6,2))
                        
                        for iM in intT:
                            popTmp, cropsTmp, urbanTmp = getExposure(haz_data, pop_data, LC_data, res, iM)
                            
                            # If the value of iM is found in roadLength, then add the value of roadLength to exposure
                            if iM in roadLength.columns:
                                roadExposure = updateRoadExposure(roadExposure, erup.name, 'Tephra', iVEI, iP, iM, None, None, roadLength[[iM]])
                            
                            exposure = updateExposure(exposure, erup.name, 'Tephra', iVEI, iP, iM, None, None, popTmp, cropsTmp, urbanTmp, None, None)
                
                    # Remove hazard file to save space on Gekko
                    if removeTephra:
                        os.remove(fl)
                    
                exposure['month'] = iMth
   
                    # Write the exposure for each volcano
            
            RSDS.to_csv(os.path.join(erup.outPath, erup.name, '_exposure/RSDS_month{}.csv'.format(iMth)))
            exposure.to_csv(os.path.join(erup.outPath, erup.name, '_exposure/exposure_month{}.csv'.format(iMth)))
            roadExposure.to_csv(os.path.join(erup.outPath, erup.name, '_exposure/road_exposure_month{}.csv'.format(iMth)))
            damageRatio.to_csv(os.path.join(erup.outPath, erup.name, '_exposure/damage_ratio_month{}.csv'.format(iMth)))
            damageState.to_csv(os.path.join(erup.outPath, erup.name, '_exposure/damage_state_month{}.csv'.format(iMth)))
        
    if analyzeTephra:
        # Close connection to rasters
        LC.close()
        pop.close()
        
        # Write the exposure for each volcano
        # RSDS.to_csv(os.path.join(erup.outPath, erup.name, '_exposure/RSDS_month.csv'))
        # exposure.to_csv(os.path.join(erup.outPath, erup.name, '_exposure/exposure_month.csv'))
        # roadExposure.to_csv(os.path.join(erup.outPath, erup.name, '_exposure/road_exposure_month.csv'))
        # damageRatio.to_csv(os.path.join(erup.outPath, erup.name, '_exposure/damage_ratio_month.csv'))
        # damageState.to_csv(os.path.join(erup.outPath, erup.name, '_exposure/damage_state_month.csv'))
    
    toc = time.perf_counter() # Stop counter
    print(f" Time: {toc - tic:0.0f} seconds")
# EXPOSURE.to_csv('results.csv')




# %% Debug
roadf = 'volcanoes/Agung/_data/roads.feather'
road = gpd.read_feather(roadf)
hazardPath = 'volcanoes/Agung/_hazard/Tephra/Agung_VEI4_P50.tif'
dictMap = {1: 10,
            100: 100}
intensity = True
epsg = 32750

ax = gdf.to_crs('EPSG:3857').plot(alpha=0.5)
ctx.add_basemap(ax)

# %%

#%%
import os
import sys
# os.chdir('/Users/seb/Documents/Codes/VolcGIS')
import volcgis
from volcgis.exposureAnalysis import *
from volcgis import eruption
from rasterio.plot import show
import rasterio as rio
from rasterio import features
import rioxarray as xrio
import pandas as pd
from pyproj import CRS
from pyproj import Transformer
import numpy as npp
# import utm
from printy import printy
# import fiona
import geopandas as gpd
from shapely.geometry import shape
# from alive_progress import alive_bar
import contextily as ctx
import time

# Tephra thresholds
# Josh: 1, 100
# Ag/Hort: 1, 5, 50
# Forest: 5, 200, 1000

# Updated 2021-06-14 for compatibility with v0.2

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
processTephra = False
processBAF = False
processPDC = False
processLC = False

processLandCover = False
processLandScan = False 
processRoadNetwork = False
processBuildings = False

analyzeTephra = False
analyzeBAF = False
analyzePDC = True
analyzeLC = False
analyzeBuffer = True

# Main loops for processing
VEI = [3,4,5]
probT = [10,50,90]
intT = [1,5,50,100]
buffT = ['300', '990']
volT = ['9800000', '450000']
radii = [10,30,100]

# Dictionaries for road disruption
RNDS_intensity_map = {1: 10,
                    100: 100}

RNDS_probability_map = {0.1: 100,
                    0.5: 100,
                    0.9: 100}

# Exposure storage
damageRatio = pd.DataFrame({})
damageState = pd.DataFrame({})
exposure = pd.DataFrame({})
roadExposure = pd.DataFrame({})

#%% Main processing step            

# Handling if the function is called from the command line                          
if (len(sys.argv)==1) or (len(sys.argv) > 2):
    rng = range(0, volcDB.shape[0])
else:
    rng = range(int(sys.argv[1]), int(sys.argv[1])+1)

# Debug
# rng = range(39,40)

# Loop through volcanoes          
for i in rng:
    tic = time.perf_counter() # Start counter
    print(volcDB.loc[i, 'name'].upper())
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

    # Setup eruption dictionary. 
    E = {
        'name':     name,
        'vent':     [lat, lon, 2765], # Vent lat, lon and elevation
        'extent':   [volcDB.loc[i, 'xmin'], volcDB.loc[i, 'xmax'], volcDB.loc[i, 'ymin'], volcDB.loc[i, 'ymax']], # extent around vent in meters, [minx maxx miny maxy]
        'epsg':     epsg
    }
    path = {
        'outPath': 'volcanoes',
        'populationPath':   'DATA/Landscan.tif',
        'landcoverPath':    'DATA/LC100_2018_croped.tif',
        'roadsPath':        'DATA/SEA_roads_criticality.gpkg'
    }
 
    # Define the eruption
    erup = volcgis.eruption.eruption(E, res,path)

    # Define the storage for RSDS
    RSDS = pd.DataFrame()

    #%% SETUP TEPHRA HAZARD     
    if processTephra:
        nameConstructor = {
            'volcano':    [name],
            'VEI':        ['VEI3', 'VEI4', 'VEI5'],
            'perc':       ['P10', 'P50', 'P90'],
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
            'buffer':     ['300m_buff', '990m_buff'],
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
    if processLandCover:
        erup.getLandcover()
    if processLandScan:
        erup.getLandscan()
    if processRoadNetwork:
        erup.getRoadNetwork()
    if processBuildings:
        erup.getBuildingExposure(inPath = 'DATA/buildings_strong_{}.tif'.format(volcDB.iloc[i].country),
                                 outPath = os.path.join(erup.path['outPath'], erup.name, '_data/buildings_strong.tif'))
        erup.getBuildingExposure(inPath = 'DATA/buildings_weak_{}.tif'.format(volcDB.iloc[i].country),
                                 outPath = os.path.join(erup.path['outPath'], erup.name, '_data/buildings_weak.tif'))
        
    
    # Exposure analysis
    LCf = os.path.join(erup.path['outPath'], erup.name, '_data/Landcover.tif')
    popf = os.path.join(erup.path['outPath'], erup.name, '_data/Landscan.tif')
    roadf = os.path.join(erup.path['outPath'], erup.name, '_data/roads.feather')
    buildWf = os.path.join(erup.path['outPath'], erup.name, '_data/buildings_weak.tif')
    buildSf = os.path.join(erup.path['outPath'], erup.name, '_data/buildings_strong.tif')
    
    # If any exposure analysis is required, load the data
    if any([analyzeBAF, analyzeLC, analyzePDC, analyzeTephra, analyzeBuffer]):

        # Read datasets
        # LC = rio.open(LCf)
        LC = xrio.open_rasterio(LCf)
        # pop = rio.open(popf)
        pop = xrio.open_rasterio(popf)
        road = gpd.read_feather(roadf)
        buildW = xrio.open_rasterio(buildWf)
        buildS = xrio.open_rasterio(buildSf)
        
        # Read and adjust
        # v0.2: Changed from rio to rioxarray behaviour
        # LC_data = LC.read(1)
        LC_data = np.squeeze(LC.data)
        # pop_data = pop.read(1)/(1000/res)**2
        pop.data = pop.data/(1000/res)**2
        # buildW_data = buildW(1)
        # buildS_data = buildS.read(1)
        buildW_data = np.squeeze(buildW.data)
        buildS_data = np.squeeze(buildS.data)
        
        # Adjust residuals from interpolation
        # pop.where(pop.data<1) = 0
        buildW_data[buildW_data<0] = 0
        buildS_data[buildS_data<0] = 0
    
    if analyzeTephra:
        printy(' - Processing tephra...')
        for iVEI in VEI:
            for iP in probT:
                printy('   - VEI: {} - Prob: {}...'.format(iVEI, iP))
                
                # Define hazard file
                fl = os.path.join(erup.path['outPath'], erup.name, '_hazard/Tephra/{}_VEI{}_P{}.tif'.format(erup.name, iVEI, int(100-iP)))
                # Get the road disruption
                rnds, roadLength, rsds = getRNDS(fl, RNDS_intensity_map, road, epsg, intensity=True)
                
                # In case of tephra, need to concatenate the rsds for various columns/accumulations into one. Cat, rename and drop
                
                # Test to handle problems on Gekko
                # if 'RSDS_100' not in rsds.columns:
                #     rsds['RSDS_100'] = 0
                # if 'RSDS_1' not in rsds.columns:
                #     rsds['RSDS_1'] = 0
                
                if 'RSDS' not in rsds.columns:
                    rsds['tephra_VEI{}_P{}'.format(iVEI, int(100-iP))] = 0
                else:
                    rsds['tephra_VEI{}_P{}'.format(iVEI, int(100-iP))] = rsds['RSDS']
                # rsds['tephra_VEI{}_P{}'.format(iVEI, int(100-iP))] = rsds.RSDS_100.fillna(0)+rsds.RSDS_1.fillna(0)
                rsds = rsds[rsds.columns.drop(list(rsds.filter(regex='RSDS')))]
                # print(rsds['tephra_VEI{}_P{}'.format(iVEI, int(100-iP))].unique())
                RSDS = updateRSDS(RSDS, rsds)

                # Compute exposure from Landscan and landcover
                with rio.open(fl) as haz:
                    haz_data = haz.read(1)
                    
                    # Update exposure
                    # Loss from weak buildings
                    DRw, DSw = getBuildingImpact(haz_data, buildW_data, vuln.loc[(volcDB.loc[i,'country'], 'Weak')], erup.path['outPath'], 'weak_VEI{}_P{}.tif'.format(iVEI, int(100-iP)), erup.name, haz.profile)
                    damageRatio, damageState = updateBuildingExposure(damageRatio,damageState,erup.name,iVEI, iP, 'weak', DRw, DSw)
                    
                    # Loss from strong buildings
                    DRs, DSs = getBuildingImpact(haz_data, buildW_data, vuln.loc[(volcDB.loc[i,'country'], 'Strong')], erup.path['outPath'], 'strong_VEI{}_P{}.tif'.format(iVEI, int(100-iP)), erup.name, haz.profile)
                    damageRatio, damageState = updateBuildingExposure(damageRatio,damageState,erup.name,iVEI, iP, 'strong', DRs, DSs)
                    
                    # General exposure
                    exposure = updateExposure(exposure, erup.name, 'Tephra', iVEI, iP, None, None, None, None, None, None, rnds, round((DRs.loss.sum()+DRw.loss.sum())/1e6,2), None, None,None)
                    
                    for iM in intT:
                        exposureTmp = getExposure(haz_data, pop, LC, res, iM, buildW=buildW_data, buildS=buildS_data)
                        
                        # If the value of iM is found in roadLength, then add the value of roadLength to exposure
                        if iM in roadLength.columns:
                            roadExposure = updateRoadExposure(roadExposure, erup.name, 'Tephra', iVEI, iP, iM, None, None, roadLength[[iM]])
                        exposure = updateExposure(exposure, erup.name, 'Tephra', iVEI, iP, iM, None, None, exposureTmp['pop_count'], exposureTmp['crops'], exposureTmp['urban'], None, None, exposureTmp['buildingsW'], exposureTmp['buildingsS'],None)
    
    if analyzeBAF:
        printy(' - Processing BAF...')
        for iV in volT:
            for iB in buffT:
                printy('   - Volume: {} - Buffer: {}...'.format(iV, iB))
                # Define hazard file
                fl = os.path.join(erup.path['outPath'], erup.name, '_hazard/BAF/{}_{}m_buff_{}m3.tif'.format(erup.name, iB, iV))
                # Get the road disruption
                rnds, roadLength, rsds = getRNDS(fl, RNDS_probability_map, road, epsg, intensity=False)
                
                # Rename columns in rsds
                colMapName = {}
                for key, value in rnds.items():
                    colMapName['RSDS_{}'.format(key)] = 'BAF_{}m_buff_{}m3_P{}'.format(iB, iV, int(key*100))
                rsds = rsds.rename(colMapName,axis=1)
                # Update RSDS
                RSDS = updateRSDS(RSDS, rsds)
                
                # Compute exposure from Landscan and landcover
                with rio.open(fl) as haz:
                    haz_data = haz.read(1)
                    for iP in probT:
                        exposureTmp = getExposure(haz_data, pop, LC, res, iP/100, buildW=buildW_data, buildS=buildS_data)
                        exposure = updateExposure(exposure, erup.name, 'BAF', None, iP, None, iB, iV, exposureTmp['pop_count'], exposureTmp['crops'], exposureTmp['urban'], rnds[iP/100], None, exposureTmp['buildingsW'], exposureTmp['buildingsS'],None)
                        roadExposure = updateRoadExposure(roadExposure, erup.name, 'BAF', None, iP, None, iB, iV, roadLength[[iP/100]])

    if analyzePDC:
        printy(' - Processing PDC...')
        for iVEI in VEI:
            printy('   - VEI: {}...'.format(iVEI))
            # Define hazard file
            fl = os.path.join(erup.path['outPath'], erup.name, '_hazard/PDC/{}_{}_output_map.tif'.format(erup.name, iVEI))
            # Get the road disruption
            rnds, roadLength, rsds  = getRNDS(fl, RNDS_probability_map, road, epsg, intensity=False)
            
            # Rename columns in rsds
            colMapName = {}
            for key, value in rnds.items():
                colMapName['RSDS_{}'.format(key)] = 'PDC_VEI{}_P{}'.format(iVEI, int(key*100))
            rsds = rsds.rename(colMapName,axis=1)
            # Update RSDS
            RSDS = updateRSDS(RSDS, rsds)
            
            # Compute exposure from Landscan and landcover
            with rio.open(fl) as haz:
                haz_data = haz.read(1)
                for iP in probT:
                    exposureTmp = getExposure(haz_data, pop, LC, res, iP/100, buildW=buildW_data, buildS=buildS_data)
                    exposure = updateExposure(exposure, erup.name, 'PDC', iVEI, iP, None, None, None, exposureTmp['pop_count'], exposureTmp['crops'], exposureTmp['urban'], rnds[iP/100], None, exposureTmp['buildingsW'], exposureTmp['buildingsS'],None)
                    roadExposure = updateRoadExposure(roadExposure, erup.name, 'PDC', iVEI, iP, None, None, None, roadLength[[iP/100]])
    
    if analyzeLC: 
        printy(' - Processing lapilli...')
        for iVEI in VEI:
            printy('   - VEI: {}...'.format(iVEI))
            # Define hazard file
            fl = os.path.join(erup.path['outPath'], erup.name, '_hazard/LC/{}_VEI{}.tif'.format(erup.name, iVEI))
            _, roadLength, _  = getRNDS(fl, RNDS_probability_map, road, epsg, intensity=False)
            with rio.open(fl) as haz:
                haz_data = haz.read(1)
                for iP in probT:
                    exposureTmp = getExposure(haz_data, pop, LC, res, iP/100, buildW=buildW_data, buildS=buildS_data)
                    exposure = updateExposure(exposure, erup.name, 'LC', iVEI, iP, None, None, None, exposureTmp['pop_count'], exposureTmp['crops'], exposureTmp['urban'], None, None, exposureTmp['buildingsW'], exposureTmp['buildingsS'],None)
                    roadExposure = updateRoadExposure(roadExposure, erup.name, 'LC', iVEI, iP, None, None, None, roadLength[[iP/100]])
                    
    if analyzeBuffer: 
        printy(' - Processing buffer...')
        exposureTmp, roadsTot = getBufferExposure(erup.buffer, pop, LC, res, buildW=buildW, buildS=buildS, roadL=road)

        for iB in exposureTmp.keys():
            exposure = updateExposure(exposure, erup.name, 'BUF', None, None, None, None, None, exposureTmp[iB]['pop_count'], exposureTmp[iB]['crops'], exposureTmp[iB]['urban'], None, None, exposureTmp[iB]['buildingsW'], exposureTmp[iB]['buildingsS'],iB)
            roadExposure = updateRoadExposure(roadExposure, erup.name, 'BUF', None, None, None, None, None, roadsTot[[iB]], iB)

        # for iR in exposureTmp
        #     # mask = features.rasterize(erup.buffer.loc[iR],out=pop.data, transform=pop.rio.transform())
        #     mask = features.rasterize(erup.buffer.loc[iR],out_shape=np.squeeze(pop.data).shape)
        #     # mask[mask<0] = 0
        #     # mask[mask>0] = 1
            
    if any([analyzeBAF, analyzeLC, analyzePDC, analyzeTephra]):
        # Close connection to rasters
        LC.close()
        pop.close()
        
        # Write the exposure for each volcano
        RSDS.to_csv(os.path.join(erup.path['outPath'], erup.name, '_exposure/RSDS.csv'))
        exposure.to_csv(os.path.join(erup.path['outPath'], erup.name, '_exposure/exposure.csv'))
        roadExposure.to_csv(os.path.join(erup.path['outPath'], erup.name, '_exposure/road_exposure.csv'))
        damageRatio.to_csv(os.path.join(erup.path['outPath'], erup.name, '_exposure/damage_ratio.csv'))
        damageState.to_csv(os.path.join(erup.path['outPath'], erup.name, '_exposure/damage_state.csv'))
    
    toc = time.perf_counter() # Stop counter
    print(f" Time: {toc - tic:0.0f} seconds")
# EXPOSURE.to_csv('results.csv')


# %% Debug

# # Prepare data
# vuln = pd.read_csv('csv/function_params.csv')
# vuln.columns = ['Country', 'Building_type', 'load_mean', 'load_disp', 'tot_rep_cost']
# vuln['load_mean'] = vuln['load_mean']*101.97 # kPa to kgm-2 

# # Damage state
# dr2ds = pd.read_csv('csv/ratio_to_damage.csv')
# dr2ds.columns = ['DS', 'DS_level ', 'cdv ', 'dr_range', 'dr_lwr', 'dr_upr']

# buildWf = 'volcanoes/Cereme/_data/buildings_weak.tif'
# fl = 'volcanoes/Cereme/_hazard/Tephra/Cereme_VEI5_P5.tif'

# haz = rio.open(fl)
# haz_data = haz.read(1)

# build = rio.open(buildWf)
# build_data = build.read(1)

# iV = 0

# getBuildingImpact(haz_data, build_data, vuln, outPath, flName, erup, profile)




# # Get the damage ratio
# f_damageRatio = lambda mn, dsp, load: scipy.stats.norm(loc = math.log(mn), scale = dsp).cdf(np.log(load))
# haz_data[haz_data<1] = 1e-6
# damageRatio = f_damageRatio(vuln.iloc[iV]['load_mean'], vuln.iloc[iV]['load_disp'], haz_data)

# # Multiply the total by the damage ratio and by the number of buildings that were exposed to a specific tephra fall load
# loss = damageRatio * vuln.iloc[iV]['tot_rep_cost'] * build_data
# loss[loss<0] = 0 # Just make sure there are no negative values

# # Convert damage ratio to damage state
# damageState = np.zeros(damageRatio.shape)
# for i in range(1,6):
#     damageState[damageRatio>=dr2ds.iloc[i]['dr_lwr']] = i



# with rio.open('test.tif', 'w', **haz.profile) as dst:
#     dst.write(DS,1)
    
# # Damage ratio table
# damageRatioR = np.round(damageRatio,1)

# DR = np.round(np.linspace(0,1,11),2)

# storDR = {
#     'damageRatio': DR,
#     'nBuildings': np.zeros(DR.shape),
#     'loss': np.zeros(DR.shape),
# }

# storDR = pd.DataFrame(storDR)
# storDR = storDR.set_index('damageRatio')

# for iR in DR:
#     idx = damageRatioR==iR
#     storDR.loc[iR, 'nBuildings'] = np.sum(build_data[idx])
#     storDR.loc[iR, 'loss'] = np.sum(loss[idx])
    

# # Damage state table

# DS = np.linspace(0,5,6)
# # DS = [int(d) for d in DS]

# storDS = {
#     'damageState': DS,
#     'nBuildings': np.zeros(DS.shape),
#     'loss': np.zeros(DS.shape),
# }

# storDS = pd.DataFrame(storDS)
# storDS = storDS.set_index('damageState')

# for iR in DS:
#     idx = damageState==iR
#     storDS.loc[iR, 'nBuildings'] = np.sum(build_data[idx])
#     storDS.loc[iR, 'loss'] = np.sum(loss[idx])
    

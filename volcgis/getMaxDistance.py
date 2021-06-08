#%% [markdown]
# This script finds the maximum extent of a hazard footprint from the vent using `xarray`. Ref:
# https://xarray.pydata.org/en/v0.10.4/auto_gallery/plot_rasterio.html

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
import time
import matplotlib.pyplot as plt
import xarray as xr

#%%
# Read the master csv file containing all the volcanoes
# The file contains these columns: name,lat,lon,xmin,xmax,ymin,ymax,country
volcDB = pd.read_csv('csv/volcCoordinates.csv')

# Read the vulnerability data
vuln = pd.read_csv('csv/function_params.csv')
vuln.columns = ['Country', 'Building_type', 'load_mean', 'load_disp', 'tot_rep_cost']
vuln['load_mean'] = vuln['load_mean']*101.97 # kPa to kgm-2 
vuln = vuln.set_index(['Country', 'Building_type'])
res = 90

# Main loops for processing
VEI = [3,4,5]
probT = [10,50,90]
intT = [1,5,50,100]
buffT = ['300', '990']
volT = ['9800000', '450000']

# Storage
maxD = pd.DataFrame()

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


    # 1. Tephra
    # if analyzeTephra:
    #     for iVEI in VEI:
    #         for iP in probT:
    #             # Define hazard file
    #             fl = os.path.join(erup.outPath, erup.name, '_hazard/Tephra/{}_VEI{}_P{}.tif'.format(erup.name, iVEI, int(100-iP)))

    #             data, _, _, dist = parseDistanceMatrix(fl, erup.vent['easting'], erup.vent['northing'])
                
    #             for iM in intT:
    #                 idx = data>=iM
    #                 distTmp = {
    #                     'volcano':      erup.name,
    #                     'hazard':       'Tephra',
    #                     'VEI':          iVEI,
    #                     'prob':         iP,
    #                     'mass':         iM,
    #                     'buffer':       None,
    #                     'volume':       None,
    #                     'buildingsLoss':np.max(dist[idx]),
    #                 }
    #             maxD = maxD.append(pd.DataFrame(distTmp))


    if analyzeLC: 
        # printy(' - Processing lapilli...')
        # for iVEI in VEI:
            printy('   - VEI: {}...'.format(iVEI))
            # Define hazard file
            fl = os.path.join(erup.outPath, erup.name, '_hazard/LC/{}_VEI{}.tif'.format(erup.name, iVEI))
            data, _, _, dist = parseDistanceMatrix(fl, erup.vent['easting'], erup.vent['northing'])

            for iP in probT:
                popTmp, cropsTmp, urbanTmp = getExposure(haz_data, pop_data, LC_data, res, iP/100)
                exposure = updateExposure(exposure, erup.name, 'LC', iVEI, iP, None, None, None, popTmp, cropsTmp, urbanTmp, None, None)

    # if analyzeBAF:
    #     printy(' - Processing BAF...')
    #     for iV in volT:
    #         for iB in buffT:
    #             printy('   - Volume: {} - Buffer: {}...'.format(iV, iB))
    #             # Define hazard file
    #             fl = os.path.join(erup.outPath, erup.name, '_hazard/BAF/{}_{}m_buff_{}m3.tif'.format(erup.name, iB, iV))
    #             # Get the road disruption
    #             rnds, roadLength, rsds = getRNDS(fl, RNDS_probability_map, road, epsg, intensity=False)
                
    #             # Rename columns in rsds
    #             colMapName = {}
    #             for key, value in rnds.items():
    #                 colMapName['RSDS_{}'.format(key)] = 'BAF_{}m_buff_{}m3_P{}'.format(iB, iV, int(key*100))
    #             rsds = rsds.rename(colMapName,axis=1)
    #             # Update RSDS
    #             RSDS = updateRSDS(RSDS, rsds)
                
    #             # Compute exposure from Landscan and landcover
    #             with rio.open(fl) as haz:
    #                 haz_data = haz.read(1)
    #                 for iP in probT:
    #                     popTmp, cropsTmp, urbanTmp = getExposure(haz_data, pop_data, LC_data, res, iP/100)
    #                     exposure = updateExposure(exposure, erup.name, 'BAF', None, iP, None, iB, iV, popTmp, cropsTmp, urbanTmp, rnds[iP/100], None)
    #                     roadExposure = updateRoadExposure(roadExposure, erup.name, 'BAF', None, iP, None, iB, iV, roadLength[[iP/100]])

    # if analyzePDC:
    #     printy(' - Processing PDC...')
    #     for iVEI in VEI:
    #         printy('   - VEI: {}...'.format(iVEI))
    #         # Define hazard file
    #         fl = os.path.join(erup.outPath, erup.name, '_hazard/PDC/{}_{}_output_map.tif'.format(erup.name, iVEI))
    #         # Get the road disruption
    #         rnds, roadLength, rsds  = getRNDS(fl, RNDS_probability_map, road, epsg, intensity=False)
            
    #         # Rename columns in rsds
    #         colMapName = {}
    #         for key, value in rnds.items():
    #             colMapName['RSDS_{}'.format(key)] = 'PDC_VEI{}_P{}'.format(iVEI, int(key*100))
    #         rsds = rsds.rename(colMapName,axis=1)
    #         # Update RSDS
    #         RSDS = updateRSDS(RSDS, rsds)
            
    #         # Compute exposure from Landscan and landcover
    #         with rio.open(fl) as haz:
    #             haz_data = haz.read(1)
    #             for iP in probT:
    #                 popTmp, cropsTmp, urbanTmp = getExposure(haz_data, pop_data, LC_data, res, iP/100)
    #                 exposure = updateExposure(exposure, erup.name, 'PDC', iVEI, iP, None, None, None, popTmp, cropsTmp, urbanTmp, rnds[iP/100], None)
    #                 roadExposure = updateRoadExposure(roadExposure, erup.name, 'PDC', iVEI, iP, None, None, None, roadLength[[iP/100]])
    

    
    
    
    
#%%
# Find the maximum distance of hazard to vent
# haz = '/Users/seb/Documents/Codes/VolcGIS/volcanoes/Agung/_hazard/Tephra/Agung_VEI4_P50.tif'
# vent = [335190, 9077710]

# data = rio.open(haz, 'r')
# # %%
# Data = data.read(1)

# #%%
# da = xr.open_rasterio(haz)
# ny, nx = len(da['y']), len(da['x'])
# x, y = np.meshgrid(da['x'], da['y'])

# data = da.data[0,:,:]

# dist = np.sqrt(np.power((x-vent[0]),2)+np.power((y-vent[1]),2))

# idx1 = data>=1

# np.max(dist[idx1])

# plt.pcolor(x,y,data)
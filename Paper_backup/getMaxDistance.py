#%% [markdown]
# This script finds the maximum extent of a hazard footprint from the vent using `xarray`. Ref:
# https://xarray.pydata.org/en/v0.10.4/auto_gallery/plot_rasterio.html

#%%
import os
import glob
import sys
# os.chdir('/Users/seb/Documents/Codes/VolcGIS')
from volcgis.exposureAnalysis import *
import volcgis.eruption as volcgis
from rasterio.plot import show
import rasterio as rio
import pandas as pd
from pyproj import CRS
from pyproj import Transformer
import numpy as np
import time
import matplotlib.pyplot as plt
import xarray as xr
import rioxarray

#%%
# Read the master csv file containing all the volcanoes
# The file contains these columns: name,lat,lon,xmin,xmax,ymin,ymax,country
volcDB = pd.read_csv('csv/volcCoordinates.csv')

res = 90

# Main loops for processing
VEI = [3,4,5]
probT = [10,50,90]
intT = [1,5,50,100]
buffT = ['300', '990']
volT = ['9800000', '450000']

# Storage
analyzeLC = True
maxD_LC = pd.DataFrame()
idxLC = 0

analyzeBAF = True
maxD_BAF = pd.DataFrame()
idxBAF = 0

analyzePDC = True
maxD_PDC = pd.DataFrame()
idxPDC = 0

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
    #                   maxD = maxD.append(pd.DataFrame(distTmp))

    # xds = xr.open_dataset(fl)

    # data = da.data[0,:,:]
    # x, y = np.meshgrid(da['x'], da['y'])
    
    if analyzeLC:
        for fl in glob.glob(f'hazards/LC_individual/runsTiff/*{erup.name}_VEI*'):
            # Define hazard file
            print(fl)
            data, _, _, dist = volcgis.parseDistanceMatrix(fl, erup.vent['easting'], erup.vent['northing'], epsg=erup.EPSG_proj)
            if 'VEI3' in fl:
                iVEI = 3
            elif 'VEI4' in fl:
                iVEI = 4
            elif 'VEI5' in fl:
                iVEI = 5
                

            idx = data>0.9 # Masked data, get everything larger than 0
            if np.any(idx):
                distTmp = {
                    'volcano':      erup.name,
                    'hazard':       'LC',
                    'VEI':          iVEI,
                    'mass':         None,
                    'prob':         None,
                    'buffer':       None,
                    'volume':       None,
                    'maxD':         round(np.max(dist[idx])),
                }
                maxD_LC = maxD_LC.append(pd.DataFrame(distTmp, index=[idxLC]))
            
            idxLC += 1
                
    if analyzeBAF:
        for fl in glob.glob(f'hazards/BAF/*{erup.name}*'):
            # Define hazard file
            print(fl)
            data, _, _, dist = volcgis.parseDistanceMatrix(fl, erup.vent['easting'], erup.vent['northing'])
            if '9800000' in fl:
                iVOL = 9800000
            elif '450000' in fl:
                iVOL = 450000
                
            if '300' in fl:
                iBUF = 300
            elif '990' in fl:
                iBUF = 990

            idx = data>0.01 # Take anything non-zero
            if np.any(idx):
                distTmp = {
                    'volcano':      erup.name,
                    'hazard':       'BAF',
                    'VEI':          None,
                    'mass':         None,
                    'prob':         None,
                    'buffer':       iBUF,
                    'volume':       iVOL,
                    'maxD':         round(np.max(dist[idx])),
                }
                maxD_BAF = maxD_BAF.append(pd.DataFrame(distTmp, index=[idxBAF]))
            
            idxBAF += 1

    if analyzePDC:
        for fl in glob.glob(f'hazards/PDC/*{erup.name}*.asc'):
            # Define hazard file
            print(fl)
            data, _, _, dist = volcgis.parseDistanceMatrix(fl, erup.vent['easting'], erup.vent['northing'])
        
            if '_3_' in fl:
                iVEI = 3
            elif '_4_' in fl:
                iVEI = 4
            elif '_5_' in fl:
                iVEI = 5
                
            for iP in [.1,.5,.9]:
                idx = data>=iP # Take anything non-zero
                if np.any(idx):
                    distTmp = {
                        'volcano':      erup.name,
                        'hazard':       'PDC',
                        'VEI':          iVEI,
                        'mass':         None,
                        'prob':         iP,
                        'buffer':       None,
                        'volume':       None,
                        'maxD':         round(np.max(dist[idx])),
                    }
                    maxD_PDC = maxD_PDC.append(pd.DataFrame(distTmp, index=[idxPDC]))
                
                idxPDC += 1

    
    
    
    
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



if analyzeLC:
    maxD_LC.to_csv('csv/maxD_LC.csv')
if analyzePDC:
    maxD_PDC.to_csv('csv/maxD_PDC.csv')
if analyzeBAF:
    maxD_BAF.to_csv('csv/maxD_BAF.csv')
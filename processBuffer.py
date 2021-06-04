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
# import utm
# from printy import printy
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

# Main loops for processing
VEI = [3,4,5]
probT = [10,50,90]
intT = [1,5,50,100]
buffT = ['300', '990']
volT = ['9800000', '450000']

# Dictionaries for road disruption
RNDS_intensity_map = {1: 10,
                    100: 100}

RNDS_probability_map = {0.1: 100,
                    0.5: 100,
                    0.9: 100}

# Exposure storage
BUFFER = pd.DataFrame()


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

    # Exposure analysis
    popf = os.path.join(erup.outPath, erup.name, '_data/Landscan.tif')
    
    pop = rio.open(popf)
    
    for iB in [10,30,100]:
        geoms = erup.buffer.loc[[iB]].geometry.values # list of shapely geometries
        geometry = geoms[0] # shapely geometry
        # transform to GeJSON format
        from shapely.geometry import mapping
        geoms = [mapping(geoms[0])]    
        
        inter, _ = rio.mask.mask(pop, geoms)
        inter = inter/(1000/res)**2
        inter[inter<1] = 0
        sum(inter[inter>0])

        # full, _ = rio.mask.mask(POP, geoms)
        
        # print(f'Buffer {iB}: original = {sum(full[full>0])}, interpolated = {sum(inter[inter>0])}')
        print(f'Buffer {iB}: {sum(inter[inter>0])}')
        
        bufData = {
            'volcano': erup.name,
            'buffer': iB,
            'pop_count': sum(inter[inter>0])
        }
        BUFFER = BUFFER.append(pd.DataFrame(bufData, index=[0]).set_index(['volcano', 'buffer']))
# %%
BUFFER.to_csv('MASTER_buffer.csv')
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
from printy import printy
import fiona
import geopandas as gpd
from shapely.geometry import shape
from alive_progress import alive_bar
import time

def updateRSDS(RSDS, rsds):
    """
    Update the rsds for a given eruption. It assumes that all the variables (i.e. hazard type, VEI, probs etc) are contained
    in the column name. Just thinking along mapping them on a GIS, it will make more sense to have all the data in one single
    file so only one join is required.
    
    Arguments:
        RSDS (pd.DataFrame): Main storage, row is `Road_ID`, column is a given hazard occurrence
        rsds (pd.DataFrame): RSDS for one hazard occurrence, row is `Road_ID`,
    """
    
    rsds_tmp = rsds.copy()
    
    # If first call of updateRSDS
    if RSDS.shape[0] == 0:
        RSDS = rsds_tmp
    else:
        RSDS = RSDS.join(rsds_tmp)
        
    return RSDS
        
def updateExposure(EXPOSURE, volcano, hazard, VEI, prob, intensity, buffer, volume, pop, crops, urban, RNDS, roadLength):

    expTmp = {
        'volcano':      [volcano],
        'hazard':       [hazard],
        'VEI':          [VEI],
        'prob':         [prob],
        'mass':         [intensity],
        'buffer':       [buffer],
        'volume':       [volume],
        'pop_count':    [pop],
        'area_crops':   [crops],
        'area_urban':   [urban],
        'RNDS':         [RNDS],
    }

    expTmp = pd.DataFrame(expTmp)
    
    # If roadLength is defined, then add the columns with length_ appended to the column name
    if roadLength is not None:
        for iRoad in range(roadLength.shape[0]):
            name = roadLength.iloc[iRoad].name
            expTmp['length_' + name] = roadLength.loc[name].values[0]
    
    return EXPOSURE.append(expTmp)

def getExposure(haz_data, pop_data, LC_data, val):
    """ Get exposure 
    
    """
    
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

def getRNDS(hazardPath, dictMap, road, epsg, intensity):
    """
    Arguments:
        hazardPath (str): Path to hazard file
        dictMap (dict): 
        road (feather): path location for the pre-processed study are road file
        intensity (bool): Defines if hazard file contains probabilities (False) or hazard intensity metrics (True)

    Returns:
        rnds (float, dict): RNDN value, either as a float (if intensity==True) or a dictionary (if intensity==False)
        roadLength (pd.DataFrame): Length of each road type defined in the `highway` column of the road variable. 
                                    Each row is a road type, each column is a single value defined in dictMap
        rsds (pd.DataFrame): RSDS value of each road segment defined by `Road_ID`
    """

    # Make sure dicMap is ordered in decreasing keys
    dictMap = dict(sorted(dictMap.items(), key=lambda item: item[0], reverse=True))
    # Convert it to a tuple so we can iterate back and forth
    inpVal = [(k, v) for k, v in dictMap.items()]
    # Output
    rnds = {}
    
    # Set id
    road = road.set_index('Road_ID')
    
    # Output dataframe that will contain the length
    roadLength = pd.DataFrame(columns=dictMap.keys(), index=road.highway.unique())

    with rio.open(hazardPath) as src:
        #Read image
        image = src.read()

        # Loop through dictMap, which returns threshold and impact score
        # for threshold, score in dictMap.items():
        for i in range(0, len(inpVal)):
            # Get threshold and score
            threshold = inpVal[i][0]
            #threshold_str = str(threshold)
            score = inpVal[i][1]

            # Create a mask
            mask = image >= threshold
            mask = mask.astype('uint8')

            # In case the hazard type is tephra and the loop is not pointing to the innermost zone,
            # then we substract the previous mask
            if i > 0 and intensity:
                maskP = image >= inpVal[i - 1][0]
                maskP = maskP.astype('uint8')
                mask = mask - maskP

            shapes = rio.features.shapes(mask, transform=src.transform)
            # Extract geometry from shapes
            geometry = []
            for shapedict, value in shapes:
                if value == 1:
                    geometry.append(shape(shapedict))

            # Create gdf for clipping
            gdf = gpd.GeoDataFrame(
                {'geometry': geometry},
                crs="EPSG:{}".format(epsg))
            # In case the mask for the given threshold is empty
            if gdf.size == 0:
                rnds[threshold] = 0

            else:
                # Create GeoDataFrame
                clipped_road = gpd.clip(road, gdf)
                clipped_road['impact_score'] = score
                clipped_road['RSDS_{}'.format(threshold)] = clipped_road['Criticality score'] * clipped_road['LoR_score'] * clipped_road['impact_score']
                rnds[threshold] = clipped_road['RSDS_{}'.format(threshold)].sum()
                
                # Calculate total road length per `highway` type and append to the storage df
                roadLengthTmp = clipped_road.groupby('highway').sum()[['Length_m']]
                roadLength.loc[roadLengthTmp.index, threshold] = roadLengthTmp.loc[roadLengthTmp.index, 'Length_m']
    
                # Append the RSDS to the full road network
                road = road.join(clipped_road[['RSDS_{}'.format(threshold)]])
    # Sum all values
    if intensity==True:
        rnds = sum(rnds.values())
    else:
        pass
    
    # Remove the highway types that don't have any data
    roadLength = roadLength.dropna(how='all')
    
    # Filter the RSDS columns to save them with the `Road_ID`
    rsds = road.filter(regex=("RSDS*"))
    rsds = rsds.dropna(how='all')

    return rnds, roadLength, rsds

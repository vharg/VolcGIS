#%%
import os
os.chdir('/Users/seb/Documents/Codes/VolcGIS')
import volcgis 
from rasterio.plot import show
import rasterio as rio
import pandas as pd
from pyproj import CRS, Transformer
import numpy as np
import scipy.stats
import math
# import utm
# from printy import printy
import fiona
import geopandas as gpd
from shapely.geometry import shape, mapping
# from alive_progress import alive_bar
import time
import json

def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    return [json.loads(gdf.to_json())['features'][0]['geometry']]

def updateRSDS(RSDS, rsds):
    """
    Update the rsds for a given eruption. It assumes that all the variables (i.e. hazard type, VEI, probs etc) are contained
    in the column name. Just thinking along mapping them on a GIS, it will make more sense to have all the data in one single
    file so only one join is required.
    
    Args:
        RSDS (pd.DataFrame): Main storage, row is `Road_ID`, column is a given hazard occurrence
        rsds (pd.DataFrame): RSDS for one hazard occurrence, row is `Road_ID`

    Returns:
        pd.DataFrame: Updated RSDS dataframe
    """
    
    rsds_tmp = rsds.copy()
    
    # If first call of updateRSDS
    if RSDS.shape[0] == 0:
        RSDS = rsds_tmp
    else:
        RSDS = RSDS.join(rsds_tmp)
        
    return RSDS
        
def updateExposure(EXPOSURE, volcano, hazard, VEI, prob, intensity, buffer, volume, pop, crops, urban, RNDS, buildingsLoss):
    """
        Updates main exposure dataframe. Buildings loss is in millions.
        
        Args:
            EXPOSURE (pd.DataFrame): Main dataframe to update
            volcano (str): Volcano name
            hazard (str): Hazard name
            VEI (float): VEI
            prob (float): Probability
            intensity (float): Hazard intensity
            buffer (float): 
            volume (float):
            pop (float):
            crops (float):
            urban (float):
            RNDS (float):
            buildingsLoss (float):
            
        Returns:
            pd.DataFrame: Updated data frame
    """
    
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
        'buildingsLoss':[buildingsLoss],
    }

    expTmp = pd.DataFrame(expTmp)
    
    # If roadLength is defined, then add the columns with length_ appended to the column name
    # if roadLength is not None:
    #     for iRoad in range(roadLength.shape[0]):
    #         name = roadLength.iloc[iRoad].name
    #         expTmp['length_' + name] = roadLength.loc[name].values[0]
    
    return EXPOSURE.append(expTmp)

def updateRoadExposure(roadExposure, volcano, hazard, VEI, prob, intensity, buffer, volume, roadLength):
    """
        Updates road exposure. Length is in km
    """
    
    expTmp = {
        'volcano':      [volcano],
        'hazard':       [hazard],
        'VEI':          [VEI],
        'prob':         [prob],
        'mass':         [intensity],
        'buffer':       [buffer],
        'volume':       [volume],
    }

    expTmp = pd.DataFrame(expTmp)

    for iRoad in range(roadLength.shape[0]):
        name = roadLength.iloc[iRoad].name
        expTmp['length_' + name] = round(roadLength.loc[name].values[0] / 1e3, 2)
    
    return roadExposure.append(expTmp)

def updateBuildingExposure(damageRatio, damageState, volcano, VEI, prob, typ, DR, DS):
    """ Update building exposure and impact to tephra accumulation
    
        Arguments:
            damageRatio (pd.DataFrame): Master damage ratio dataframe
            damageState (pd.DataFrame): Master damage state dataframe
            volcano (str): Volcano name
            VEI (str): VEI
            prob (str): Probability 
            typ (str): Building type
            DR (pd.DataFrame): Damage ratio dataframe for a given volcano/vei/prob
            DS (pd.DataFrame): Damage state dataframe for a given volcano/vei/prob

            
    """
    DR['volcano'] = volcano
    DR['VEI'] = VEI
    DR['prob'] = prob
    DR['type'] = typ
    
    DS['volcano'] = volcano
    DS['VEI'] = VEI
    DS['prob'] = prob
    DS['type'] = typ

    return damageRatio.append(DR), damageState.append(DS)

def getExposure(haz_data, pop_data, LC_data, res, val, LC={'crops':40, 'urban':50}):
    """ This is the most basic exposure function. It compute the exposure of i) population from Landsat, 
            ii) urban area and iii) crop area from CGLS from a simple mask with hazard data. For the landcover 
            data, using `CGLS` so `urban==50` and `crops==40`.   
            
        - 2021-06-09: Made the function flexible to account for variable LC classes
    
        Args:
            haz_data (np.array): Hazard data
            pop_data (xarray): Population density data
            LC_data (xarray): Landcover data
            res (int): Raster resolution 
            val (float): Hazard value to mask
            LC (dict): Dictionary containing `'class_name': class_val`

        Returns:
            dict: A dictionary containing `pop_count` (float) and one exposure value for each dict entry in `LC`: 
    """
    
    exposure = {}

    # Population
    idx = (haz_data >= val) & (pop_data.data[0] > 0)
    popTmp = np.nansum(np.nansum(pop_data.where(idx).data[0]))
    exposure['pop_count'] = round(popTmp,0)
       
    # Landcover
    for LCi in LC.keys():
        # Landcover / crops km2
        idx = (haz_data >= val) & (LC_data.data[0]==LC[LCi])
        exposure[LCi] = np.sum(idx)*res**2/1e6

    return exposure 

def getBufferExposure(buffer, pop_data, LC_data, res, LC={'crops':40, 'urban':50}):
    """ Get exposure as concentric circles around the vent
    
    Args:
        buffer (pd.DataFrame): List of buffers as defined in `erup.buffer`
        pop_data (xarray): Population density data
        LC_data (xarray): Landcover data
        res (int): Raster resolution 
        val (float): Hazard value to mask
        LC (dict): Dictionary containing `'class_name': class_val`

    Returns:
        dict: A dictionary containing `pop_count` (float) and one exposure value for each dict entry in `LC`: 

    """

    bufferTot = {}

    # Loop through buffers
    for iB in buffer.index.values:
        bufferTmp = {}
        
        # Create GDF for clipping
        geoms = gpd.GeoDataFrame(buffer.loc[[iB]])

        # Population
        popB = pop_data.rio.clip(geoms.geometry)
        idx = popB.data[0] > 0
        popTmp = np.nansum(np.nansum(popB.where(idx).data[0]))
        bufferTmp['pop_count'] = round(popTmp,0)
        
        # Landcover
        LCB = LC_data.rio.clip(geoms.geometry)
        for LCi in LC.keys():
            idx = LCB.data[0]==LC[LCi]
            bufferTmp[LCi] = round(np.sum(idx)*res**2/1e6,0)

        bufferTot[iB] = bufferTmp

    return bufferTot

def getRNDS(hazardPath, dictMap, road, epsg, intensity):
    """
        Arguments:
            hazardPath (str): Path to hazard file
            dictMap (dict): 
            road (gpd): 
            epsg (int): Digits of the epsg code
            intensity (bool): Defines if hazard file contains probabilities (False) or hazard intensity metrics (True)

        Returns:
            tuple: A tuple containing: 
            
            - rnds (float, dict): RNDN value, either as a float (if intensity==True) or a dictionary (if intensity==False)  
            - roadLength (pd.DataFrame): Length of each road type defined in the `highway` column of the road variable. 
                                Each row is a road type, each column is a single value defined in dictMap  
            - rsds (pd.DataFrame): RSDS value of each road segment defined by `Road_ID`  
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

def getBuildingImpact(haz_data, build_data, vuln, outPath, flName, erup, profile):
    """
    Physical impact on buildings from tephra fallout
    
    Args:
        vuln (pd.Series): CSV file containing the parameters of the fragility curves. Must contain the following fields: 
            load_mean, load_disp, tot_rep_cost
        outPth (str): Main output directory
        flName (str): String appended to the raster names. If ``None``, no raster is written
        erup (str): Volcano name
        profile (dict): Rasterio reference profile for writing output
    """

    # Import damage state
    dr2ds = pd.read_csv('csv/ratio_to_damage.csv')
    dr2ds.columns = ['DS', 'DS_level ', 'cdv ', 'dr_range', 'dr_lwr', 'dr_upr']
    
    # Get the damage ratio
    f_damageRatio = lambda mn, dsp, load: scipy.stats.norm(loc = math.log(mn), scale = dsp).cdf(np.log(load))
    haz_data[haz_data<.1] = np.nan
    damageRatio = f_damageRatio(vuln['load_mean'], vuln['load_disp'], haz_data)

    # Multiply the total by the damage ratio and by the number of buildings that were exposed to a specific tephra fall load
    loss = damageRatio * vuln['tot_rep_cost'] * build_data
    loss[loss<0] = 0 # Just make sure there are no negative values

    # Convert damage ratio to damage state
    damageState = np.zeros(damageRatio.shape)
    for i in range(1,6):
        damageState[damageRatio>=dr2ds.iloc[i]['dr_upr']] = i
    
    ## Write the rasters 
    if flName is not None:
    # Write damage ratio as raster
        with rio.open(os.path.join(outPath, erup, '_exposure/damageRatio_{}'.format(flName)), 'w', **profile) as dst:
            dst.write(damageRatio,1)
        
        # Write loss as raster
        with rio.open(os.path.join(outPath, erup, '_exposure/loss_{}'.format(flName)), 'w', **profile) as dst:
            dst.write(loss,1)
            
    ## Write the tiffs
    # Damage ratio table
    damageRatioR = np.round(damageRatio,1)

    DR = np.round(np.linspace(0,1,11),2)

    storDR = {
        'damageRatio': DR,
        'nBuildings': np.zeros(DR.shape),
        'loss': np.zeros(DR.shape),
    }

    storDR = pd.DataFrame(storDR)
    storDR = storDR.set_index('damageRatio')

    for iR in DR:
        idx = damageRatioR==iR
        storDR.loc[iR, 'nBuildings'] = np.sum(build_data[idx])
        storDR.loc[iR, 'loss'] = np.sum(loss[idx])
        

    # Damage state table
    DS = np.linspace(0,5,6)
    # DS = [int(d) for d in DS]

    storDS = {
        'damageState': DS,
        'nBuildings': np.zeros(DS.shape),
        'loss': np.zeros(DS.shape),
    }

    storDS = pd.DataFrame(storDS)
    storDS = storDS.set_index('damageState')

    for iR in DS:
        idx = damageState==iR
        storDS.loc[iR, 'nBuildings'] = np.sum(build_data[idx])
        storDS.loc[iR, 'loss'] = np.sum(loss[idx])
        
    return storDR, storDS
# Integration of the road disruption index
import rasterio as rio
from shapely.geometry import shape
import matplotlib.pyplot as plt

# References
# [Contour value with rasterio](https://stackoverflow.com/questions/41487642/how-to-export-contours-created-in-scikit-image-find-contours-to-shapefile-or-geo)
# [rasterio features to gdf](https://gis.stackexchange.com/questions/379412/how-to-create-a-geopandas-geodataframe-from-rasterio-features)

#%%
def joshifyMe(flPath, dictMap, epsg, intensity=False):
    """
    
        Arguments:
            flPath (str): Path to hazard file
            dictMap (dict):
            epsg (int): EPSG code of hazard file
            intensity (bool): Defines if hazard file contains probabilities (False) or hazard intensity metrics (True)
    
        Returns:
            RNDN (dict): Dictionary of RNDN
    """
    
    
    # Make sure dicMap is ordered in decreasing keys
    dictMap = dict(sorted(dictMap.items(), key=lambda item: item[0], reverse=True))
    # Convert it to a tuple so we can iterate back and forth
    inpVal = [(k,v) for k,v in dictMap.items()]
    # Output
    RNDN = dict()
    
    with rio.open(flPath) as src:
        # Read image
        image = src.read()
        
        # Loop through dictMap, which returns threshold and impact score
        # for threshold, score in dictMap.items():
        for i in range(0, len(inpVal)): 
            # Get threshold and score
            threshold = inpVal[i][0]
            score = inpVal[i][1]
            
            # Create a mask
            mask = image >= threshold
            mask = mask.astype('uint8')
            
            # In case the hazard type is tephra and the loop is not pointing to the innermost zone, then we substract the previous mask 
            if i>0 and intensity:
                maskP = image >= inpVal[i-1][0]
                maskP = maskP.astype('uint8')
                mask = mask-maskP
                
            # Contour values from raster
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
                RNDN[threshold] = 0
            
            else:
                print('yey!')
                
                # Apply your function using `score` as the score and `clipMask` as the hazard extent
                
                # RNDN[threshold] = __output_of_your_function
                
        
        # if intensity:
        #     test
    return(RNDN)

#%%

# If any other hazard than tephra, define a dictionary where the impact score is mapped to a probability (i.e. 0-1)
dictMap = {0.1: 100,
           0.5: 100,
           0.9: 100}
epsg = 32748

flPath = 'volcanoes/Gede-Pangrango/_hazard/BAF/Gede-Pangrango_990m_buff_9800000m3.tif'
flPath = 'volcanoes/Gede-Pangrango/_hazard/LC/Gede-Pangrango_VEI4.tif'
flPath = 'volcanoes/Gede-Pangrango/_hazard/PDC/Gede-Pangrango_5_output_map.tif'

joshifyMe(flPath,dictMap)

#%%

# If tephra, define a dictionary where the impact score is mapped to a hazard intensity (i.e. 1, 100)
dictMap = {1: 10,
           100: 100}
epsg = 32748
flPath = 'volcanoes/Gede-Pangrango/_hazard/Tephra/Gede-Pangrango_VEI5_P5.tif'

joshifyMe(flPath, dictMap, epsg, intensity=True)
# %%

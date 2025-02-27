#%%
# import pandas as pd
# import geopandas as gpd
import numpy as np
# import matplotlib.pyplot as plt
# import itertools
# from shapely.geometry import box
# from fiona.crs import from_epsg
# import osmnx as ox
# from osgeo import gdal
# import rasterio as rio
# from rasterio.plot import show
# from rasterio.warp import calculate_default_transform, reproject, Resampling
# from rasterio.windows import from_bounds
from rasterio import shutil as rio_shutil
# from rasterio.vrt import WarpedVRT
# import affine
# import contextily as ctx
import xarray as xr
# import rioxarray as xio
from pyproj import Transformer, CRS
# import os

# #%%
# class eruption:
#     def __init__(self, eruption, masterRes, path, extent='absolute', buffer=None):
#         ''' Sets main eruption object.
        
#         Args:
#             eruption (dict): Eruption dictionary
#             masterRes (int): Cell size (m) used throughout analyses
#             path (dict): Dictionary containing paths to various file. The keys accept `'outPath'`, `'LandscanPath'`, `'LandcoverPath'`
#             extent (str): Defines information stored in `eruption.extent`. If `absolute`, `eruption.extent` contains absolute UTM 
#                 coordinates. If `relative` it contains relative distances (m) from the vent. In both cases, the `eruption.extent`
#                 should contain [`xmin`, `xmax`, `ymin`, `ymax`].
#             buffer (list[int]): List of buffers (km) to compute exposure as concentric circles around the vent.
#                 If None, set to [10, 30, 100].
        
#         Returns:
#             object: eruption object
#         '''

#         self.path = path
#         # self.path['outPath'] = 'volcanoes' if outPath == None else outPath
#         # Set default variables
#         self.name = eruption['name']
#         self.res = masterRes

#         # Set path and folders
#         if not os.path.exists(os.path.join(self.path['outPath'], self.name)):
#             os.mkdir(os.path.join(self.path['outPath'], self.name))
#             os.mkdir(os.path.join(self.path['outPath'], self.name, '_data'))
#             os.mkdir(os.path.join(self.path['outPath'], self.name, '_tmp'))
#             os.mkdir(os.path.join(self.path['outPath'], self.name, '_hazard'))
#             os.mkdir(os.path.join(self.path['outPath'], self.name, '_exposure'))
        
#         # Define projections
#         self.EPSG = eruption['epsg']
#         self.EPSG_geo = 'epsg:{}'.format(4326)
#         self.EPSG_proj = 'epsg:{}'.format(self.EPSG)

#         # Convert vent into a geometry
#         transformer = Transformer.from_crs(self.EPSG_geo, self.EPSG_proj)
#         [xtmp, ytmp] = transformer.transform(eruption['vent'][0], eruption['vent'][1])
#         self.vent = {'lat': eruption['vent'][0], 'lon': eruption['vent'][1], 'easting': round(xtmp), 'northing': round(ytmp), 'alt': eruption['vent'][2]}

#         # Create area mask
#         if extent == 'relative':
#             eruption['extent'][0] = self.vent['easting'] - eruption['extent'][0] 
#             eruption['extent'][1] = self.vent['easting'] + eruption['extent'][1] 
#             eruption['extent'][2] = self.vent['northing'] - eruption['extent'][2] 
#             eruption['extent'][3] = self.vent['northing'] + eruption['extent'][3] 
            
#         areaPl = box(eruption['extent'][0], 
#                      eruption['extent'][2],
#                      eruption['extent'][1],
#                      eruption['extent'][3])
#         self.area = gpd.GeoDataFrame({'geometry': areaPl}, index=[0], crs=self.EPSG_proj) # Area is now projected
#         self.areaG = self.area.to_crs(from_epsg(4326)) # Project it

#         # Define buffers
#         tmp = pd.DataFrame(self.vent, index=[0])
#         gdf = gpd.GeoDataFrame(tmp, geometry=gpd.points_from_xy(tmp.easting, tmp.northing), crs=self.EPSG_proj)
        
#         if buffer is None:
#             bufList = [10,30,100]
#         else:
#             bufList = buffer
            
#         buffer = pd.DataFrame()
#         for iB in bufList:
#             buffer = buffer.append(gdf.geometry.buffer(iB*1e3).rename(iB))
#         buffer.columns = ['geometry']
#         self.buffer = buffer
        
#         # Create an empty dictionary to store hazard information
#         self.hazards = {}
        
#         # Create an empty dictionary to store exposure information
#         self.exposure = {}
        
#         # Define reference for virtual wrap on which all rasters will be aligned
#         self.ref = {}
#         self.ref['bounds'] = self.area.geometry[0].bounds
#         self.ref['height'] = len(np.arange(int(self.ref['bounds'][1]), int(self.ref['bounds'][3]), self.res))
#         self.ref['width'] = len(np.arange(int(self.ref['bounds'][0]), int(self.ref['bounds'][2]), self.res))
#         self.ref['transform'] = affine.Affine(self.res, 0.0, self.ref['bounds'][0], 0.0, -self.res, self.ref['bounds'][3])
#         self.ref['EPSG'] = rio.crs.CRS.from_epsg(self.EPSG)

#     def prepareExposure(self, population=True, landcover=True, populationInt='nearest', landcoverInt='nearest',
#                         populationRes=1000, landcoverRes=100, populationScaling=True):
#         """ Prepares datasets for exposure analysis.
#             Args:
#                 population (bool): Process population data
#                 landcover (bool): Process landcover data
#                 populationInt (str): Interpolation method for population data, accepts `'nearest'` and `'cubic'`
#                 landcoverInt (str): Interpolation method for landcover data, accepts `'nearest'` and `'cubic'`
#                 populationRes (int): Original resolution of population data
#                 landcoverRes (int): Original resolution of landcover data
#                 populationScaling (bool): Corrects for the change in pixel resolution of population count data
#         """
        
#         expType = {}
#         if population:
#             self.path['pop_count'] = os.path.join(self.path['outPath'], self.name, '_data', 'population.tif')
#             profile = self.alignRaster(self.path['populationPath'], self.path['pop_count'], resampling=populationInt)
            
#             # Post processing
#             # First open the raster in read mode, retrieve data and profile and close it
#             with rio.open(self.path['pop_count'], 'r') as src:
#                 profile = src.profile
#                 data = src.read(1)
#             #Then open the raster in write mode
#             with rio.open(self.path['pop_count'], 'w', **profile) as dst:
#                 # Scale data
#                 if populationScaling:
#                     data = data/(populationRes/self.res)**2
#                 # Remove negative values
#                 data[data<1] = 0
#                 dst.write(np.round(data).astype(rio.int32), indexes=1)
#             expType['pop_count'] = 1
            
#         if landcover:
#             self.path['LC_class'] = os.path.join(self.path['outPath'], self.name, '_data', 'landcover.tif')
#             self.alignRaster(self.path['landcoverPath'], self.path['LC_class'], resampling=landcoverInt)
#             expType['LC_class'] = 1
        
#         self.exposure['expType'] = expType
    
#     def getHazardExposure(self, hazard, hazardProps, LC={'crops':40, 'urban':50}):
        
#         columns = hazardProps['columns'] + [hazardProps['varName']]
        
#         if 'pop_count' in self.exposure['expType'].keys():
#             pop_data = xio.open_rasterio(self.path['pop_count'])
#             columns = columns + ['pop_count']
#             pop_count = True
            
#         if 'LC_class' in self.exposure['expType'].keys():
#             LC_data = xio.open_rasterio(self.path['LC_class'])
#             columns = columns + list(LC.keys())
#             LC_class = True
            
#         exposure = pd.DataFrame(columns=columns)
        
#         hazPath = os.path.join(self.path['outPath'], self.name, '_hazard', hazard)
        
#         # Loop through hazard files
#         for i in range(0,self.hazards[hazard]['data'].shape[0]):
#             haz_data = xio.open_rasterio(os.path.join(hazPath, self.hazards[hazard]['data'].iloc[i]['filePth']))

#             # Loop through the thresholds
#             for iT in hazardProps['varVal']:

#                 # Create dictionary to update the exposure dataframe
#                 exposureTmp = {}
#                 exposureTmp[hazardProps['varName']] = iT
#                 for iC in hazardProps['columns']:
#                     exposureTmp[iC] = self.hazards[hazard]['data'].iloc[i][iC]
                
#                 # Population
#                 if pop_count:
#                     idx = (haz_data >= iT) & (pop_data.data[0]>0)
#                     popTmp = np.nansum(np.nansum(pop_data.where(idx).data[0]))
#                     exposureTmp['pop_count'] = np.round(popTmp,0)

#                 # Landcover
#                 if LC_class:
#                     for LCi in LC.keys():
#                         idx = (haz_data >= iT) & (LC_data.data[0]==LC[LCi])
#                         exposureTmp[LCi] = np.round(np.sum(idx.data[0])*self.res**2/1e6,0)
                
#                 # Update exposure
#                 exposure = exposure.append(
#                     pd.DataFrame(exposureTmp, columns=columns, index=[0])
#                 )
                
#         self.exposure[hazard] = exposure
                
#     def getBufferExposure(self, LC={'crops':40, 'urban':50}):
#         """ Get exposure as concentric circles around the vent for buffers defined in self.buffer for exposure
#                 datasets defined in self.exposure['expType']
        
#         Args:
#             LC (dict): Dictionary containing `'class_name': class_val`

#         Returns:
#             pd.DataFrame: A dataFrame saved in `self.exposure['buffer']`. Landcover-based exposure indices are in km2. 

#         """

#         index = []
#         if 'pop_count' in self.exposure['expType'].keys():
#             pop_data = xio.open_rasterio(self.path['pop_count'])
#             index = index + ['pop_count']
#             pop_count = True
            
#         if 'LC_class' in self.exposure['expType'].keys():
#             LC_data = xio.open_rasterio(self.path['LC_class'])
#             index = index + list(LC.keys())
#             LC_class = True
            
#         buffer = pd.DataFrame(columns=self.buffer.index.values, index=index)
            
#         # Loop through buffers
#         for iB in self.buffer.index.values:            
#             # Create GDF for clipping
#             geoms = gpd.GeoDataFrame(self.buffer.loc[[iB]])

#             # Population
#             if pop_count:
#                 popB = pop_data.rio.clip(geoms.geometry)
#                 idx = popB.data[0] > 0
#                 popTmp = np.nansum(np.nansum(popB.where(idx).data[0]))
#                 buffer.loc['pop_count', iB] = round(popTmp,0)
            
#             # Landcover
#             if LC_class:
#                 LCB = LC_data.rio.clip(geoms.geometry)
#                 for LCi in LC.keys():
#                     idx = LCB.data[0]==LC[LCi]
#                     buffer.loc[LCi, iB] = round(np.sum(idx)*self.res**2/1e6,0)

#         self.exposure['buffer'] = buffer

#     def getLandscan(self, inPath=None):
#         """ Retrieves Landscan data for the area defined by self.area.

#         Args:
#             inPath (str): Path to corresponding raster. If `None`, set to `'DATA/Landscan.tif'`

#         Returns:
#             image: A geotif cropped, aligned and projected to `self.ref['bounds']`

#         """
#         if inPath is None:
#             inPath = 'DATA/Landscan.tif'

#         outPath = os.path.join(self.path['outPath'], self.name, '_data', 'Landscan.tif')

#         originalRes = 1000 # Landscan resolution
#         scaling = originalRes / self.res # Scaling factor to correct population

#         self.alignRaster(inPath, outPath, resampling='nearest', scalingFactor=scaling)

#     def getLandcover(self, inPath=None):
#         """ Retrieves Landcover data for the area defined by self.area. Resampling is set to 'nearest' for discrete data

#        Args:
#             inPath (str): Path to corresponding raster. If `None`, set to `'DATA/LC100_2018_croped.tif'`

#         Returns:
#             image: A geotif cropped, aligned and projected to `self.ref['bounds']`

#         """
#         if inPath is None:
#             inPath = 'DATA/LC100_2018_croped.tif'

#         outPath = os.path.join(self.path['outPath'], self.name, '_data', 'Landcover.tif')
#         self.alignRaster(inPath, outPath, resampling='nearest')

#     def getBuildingExposure(self, inPath=None, outPath=None):
#         """ Retrieves building exposure from George's analysis for area defined by self.area.

#        Args:
#             inPath (str): Path to corresponding raster

#         Returns:
#             image: A geotif cropped, aligned and projected to `self.ref['bounds']`

#         """
#         self.alignRaster(inPath, outPath, resampling='nearest')
#         # if inPath is None:
            
#     def getRoadNetwork(self, inPath=None):
#         """

#         Args:
#             inPath (str): Path to the main roads dataset. If `None`, set to `'DATA/SEA_roads_criticality.gpkg'`

#         Returns:
#             roads (feather): A feather file of the roads within BBox

#         """
        
#         # Re-project self.areaG to pseudo mercator 3857
#         bbox = self.area.to_crs(from_epsg(3857))

#         # Get bbox for the geometry
#         if inPath is None:
#             inPath = 'DATA/SEA_roads_criticality.gpkg'

#         roads = gpd.read_file(inPath, bbox=bbox['geometry'])

#         # Reproject to self.EPSG
#         roads = roads.to_crs(crs="EPSG:{}".format(self.EPSG))

#         # Save as a feather to outPath
#         outPath = os.path.join(self.path['outPath'], self.name, '_data', 'roads.feather')
#         roads.to_feather(outPath)
         
#     def prepareHazard(self, hazard, noAlign=False):
#         """ Prepare the hazard layers

#             Loads the files for a given hazard defined as hazard['hazard'] and from hazard['nameConstructor'] from the hazards/ folder.

#             Args:
#                 hazard (dict): Main hazard dictionary
#                 noAlign (bool): Defines if hazard rasters are aligned

#             Returns:
#                 Geotif: A geotif image
#         """
#         print('- Clipping hazard files...')
#         def makeInputFileName(rootDir, nameConstructor):
#             ''' Create list of filenames based on nameConstructor '''
#             flName = []
#             flData = pd.DataFrame(columns=list(nameConstructor.keys())+['filePth'])
#             for p in itertools.product(*list(nameConstructor.values())):
#                 fl = '_'.join(p)
#                 flOut = fl[:-5]+fl[-4:]
#                 fl = rootDir+fl[:-5]+fl[-4:]
#                 flName.append(fl)
#                 flData = flData.append(pd.DataFrame([list(p)+[flOut]], columns=list(nameConstructor.keys())+['filePth']))
#             return flName, flData
    
#         print('Preparing hazard: {}'.format(hazard['hazard']))
#         # Set path and folders
#         targetDir = os.path.join(self.path['outPath'], self.name, '_hazard', hazard['hazard'])
#         if not os.path.exists(targetDir):
#             os.mkdir(targetDir)

#         # Create a list of file names
#         hazFl, hazard['data'] = makeInputFileName(hazard['rootDir'], hazard['nameConstructor'])

#         # Align files
#         for inPath in hazFl:
#             outPath = inPath.replace(hazard['rootDir'], '')
#             print(f'   - {outPath}...')
#             outPath = outPath.replace('.asc', '.tif')
#             outPath = '{}/{}/_hazard/{}/{}'.format(self.path['outPath'], self.name, hazard['hazard'], outPath)
            
#             if not noAlign:
#                 print('  - Processing: {}'.format(outPath))
#                 self.alignRaster(inPath, outPath, epsg=hazard['epsg'])

#         # Save hazard data
#         self.hazards[hazard['hazard']] = hazard
        
#     def alignRaster(self, inPath, outPath, epsg=None, resampling='cubic', scalingFactor=None):
#         """
#             Aligns raster to a reference grid using a virtual wrap. Use reference contained in self.ref to read a window of original file and wrap it

#             Args:
#                 inPath (str): Path to input raster
#                 outPath (str): Path to output raster
#                 resampling (str): Resampling method
#                 scalingFactor (float):

#             Returns:
#                 image: A geotif cropped, aligned and projected to `self.ref['bounds']`

#         """
#         # Virtual wrapper options
#         vrt_options = {
#             'resampling': Resampling.cubic,
#             'crs': self.ref['EPSG'],
#             'transform': self.ref['transform'],
#             'height': self.ref['height'],
#             'width': self.ref['width'],
#             'driver': 'GTiff'
#         }

#         if resampling == 'nearest':
#             vrt_options['resampling'] = Resampling.nearest

#         with rio.open(inPath, 'r') as src:
#             # In case no EPSG is specified (e.g. asc)
#             if src.crs == None:
#                  src.crs = rio.crs.CRS.from_epsg(epsg)
#             with WarpedVRT(src, **vrt_options) as vrt:
#                 rst = vrt.read(1, window=from_bounds(self.ref['bounds'][0], self.ref['bounds'][1], self.ref['bounds'][2], self.ref['bounds'][3], self.ref['transform']))
#                 rio_shutil.copy(vrt, outPath, driver='GTiff', compress='lzw')

#         # if scalingFactor is not None:
#             # with rio.open(outPath, 'w', **profile) as src:
#             #     data = np.round(src.read(1)/(scalingFactor**2))
#             #     src.write(data.astype(rio.int32))
#             # First open the raster in read mode, retrieve data and profile and close it
#             # with rio.open(outPath, 'r') as src:
#             #     profile = src.profile
#             #     data = src.read(1)
#             # Then open the raster in write mode
#             # with rio.open(outPath, 'w', **profile) as src:
#             #     data = np.round(src.read(1)/(scalingFactor**2))
#             #     src.write(data.astype(rio.int32))

#         return vrt_options
    
#     def plot(self, plotExposure=None, plotBuffer=None, plotHazard=None, hazLevels=None, hazProps=None, LC={'crops':40, 'urban':50}, figsize = [8,5]):
#         """ Plot exposure on a map
        
#             Args:
#                 plotExposure (str): Type of exposure to plot, accepts `LC`, `pop`
#                 plotBuffer (bool): Controls if buffers are plotted. They need to be defined in `self.buffer`
#                 plotHazard (str): Type of hazard to plot. Must be an entry in `self.hazards`
#                 hazLevels (list[int]): Hazard values to contour
#                 hazProps (dict): Reference dict to get the hazard file to plot, must correspond to the columns
#                     defined in `self.hazards[plotHazard]` and return a unique value
#                 LC (dict): Dictionary containing `'class_name': class_val`
#                 figsize (list[float]): Figure size

#         """
#         fig = plt.figure(figsize=figsize)
#         ax = fig.add_subplot(1, 1, 1)
#         self.areaG.plot(facecolor="none", edgecolor="black",alpha=0.5,ax=ax)

#         # Setup title
#         ttl = self.name
        
#         # Plot population
#         if plotExposure in ['population', 'pop']:
#             pop_data = xio.open_rasterio(self.path['pop_count'])
#             pop_data = pop_data.rio.reproject(self.areaG.crs.to_string(), resampling=0)
#             vmax = np.percentile(pop_data.data,90)
#             pop_data.plot(cmap='viridis', alpha=.5, ax=ax, vmin=0, vmax=vmax, cbar_kwargs={'label': 'Population count'}, levels=[0,10,25,50,75,100,150,200])

#         # Plot landcover
#         elif plotExposure in ['landcover', 'lc', 'LC']:
            
#             LC_data = xio.open_rasterio(self.path['LC_class'])
#             LC_data = LC_data.rio.reproject(self.areaG.crs.to_string(), resampling=0)

#             # Update mask
#             cLC = 0
#             for iLC in LC.keys():
#                 if cLC == 0:
#                     idx = LC_data.data[0]==LC[iLC]
#                 else:
#                     idx = np.logical_or(idx, LC_data.data[0]==LC[iLC])
#                 cLC += 1
#             LCplt = LC_data.where(idx).plot(cmap='Set1', alpha=.5, ax=ax, cbar_kwargs={'label': 'Landcover'}, levels=len(list(LC.values()))+1 )#, levels=list(LC.values()))#, vmin=min(list(LC.values())), vmax=max(list(LC.values())))
#             LCplt.colorbar.set_ticks(list(LC.values()))
#             LCplt.colorbar.set_ticklabels(list(LC.keys()))
            
#         # Plot buffer
#         if plotBuffer:
#             buf_data = gpd.GeoDataFrame(self.buffer).set_crs(epsg=self.EPSG)
#             buf_data.to_crs('EPSG:4326').plot(facecolor="none", edgecolor="yellow", alpha=0.7, linewidths=2, ax=ax)

#         # Plot hazard
#         if plotHazard is not None:
#             hazList = self.hazards[plotHazard]['data']
            
#             # Retrieve the hazard file
#             ttl = f'{ttl} - {plotHazard} ('
#             for i in hazProps.keys():
#                 hazList = hazList[hazList[i] == hazProps[i]]
#                 ttl = f'{ttl}{i}: {hazProps[i]}, '
#             ttl = f'{ttl[:-2]})'
#             hazPath = os.path.join(self.path['outPath'], self.name, '_hazard',plotHazard, hazList.iloc[0]['filePth'])
#             haz_data = xio.open_rasterio(hazPath)
#             haz_data = haz_data.rio.reproject(self.areaG.crs.to_string(), resampling=0)
            
#             if plotExposure is None:
#                 haz_data.where(haz_data.data>0).plot(cmap='cividis', alpha=.5, ax=ax)
#             CS = haz_data.squeeze().plot.contour(levels=hazLevels, ax=ax, colors='black', linewidths=2)
#             plt.clabel(CS, inline=1, fontsize=10)

#         # Plot vent
#         ax.plot(self.vent['lon'], self.vent['lat'], '^r')

#         # add basemap
#         ctx.add_basemap(ax,crs=self.areaG.crs.to_string())

#         # Setup labels
#         plt.title(ttl)
#         plt.ylabel("Latitude")
#         plt.xlabel("Longitude")
#         plt.tight_layout()
        
#     def plotArea(self):
#         """
#             Plots the extent of self.area on a map with a basemap
#         """
        
#         transformer = Transformer.from_crs(self.EPSG_geo, 'EPSG:3857')
#         [xtmp, ytmp] = transformer.transform(self.vent['lat'], self.vent['lon'])

#         fig = plt.figure(figsize=[8,10])
#         ax = fig.add_subplot(1, 1, 1)
#         self.areaG.to_crs('EPSG:3857').plot(alpha=0.5,ax=ax)
#         ax.plot(xtmp, ytmp, '^r')
#         ctx.add_basemap(ax)

# def makeZoneBoundaries(x, y, xmin, xmax, ymin, ymax):
#         """ Calculates the boundaries in m from a central point. Useful when
#             a boundary falls outside the limit of the zone defined by the central
#             point.

#         Args:
#             x (float): Easting of the central point
#             y (float): Northing of the central point
#             xmin (float): Minimum Easting of the bounding box
#             xmax (float): Maximum Easting of the bounding box
#             ymin (float): Minimum Northing of the bounding box
#             ymax (float): Maximum Northing of the bounding box

#         Returns:
#             tuple: A tuple containing `xmin`, `xmax`, `ymin`, `ymax`

#         Todo:
        
#         - Maybe obsolete?
            
#         """

#         if xmin<0:
#             xmin = -1*(round(x))+xmin
#         else:
#             xmin = round(x)-xmin

#         xmax = xmax - x

#         if ymin<0:
#             ymin = -1*(round(x))+ymin
#         else:
#             ymin = round(x)-ymin

#         ymax = ymax - y

#         return xmin, xmax, ymin, ymax

def parseDistanceMatrix(pth, cX, cY, band=0, epsg=None):
    """
    Loads a geotif with xarray and returns the matrices containing data, 
        x coordinates, y coordinates and distance from a central point.
        Adapted from https://xarray.pydata.org/en/v0.10.4/auto_gallery/plot_rasterio.html
    
        Args:
            pth (str): Path to the geotif file
            cX (float): x coordinate of the centra point from which distance is calculated. Must be in the same reference coordinate as `pth`
            cY (float): Same as `cX` for y coordinate
            band (int): The band number of the data dimension to read
            epsg (int): EPSG string (e.g.'epsg:32749') to which the raster will be reprojected before computing distance
            
        Returns:
            tuple: A tuple containing the following matrices, all of which have the same dimension: 
            
            - `data (np.array)`: The data matrix
            - `x (np.array)`: The x-coordinate matrix
            - `y (np.array)`: The y-coordinate matrix
            - `dist (np.array)`: The matrix containing distances from `[cX, cY]` in the same unit as `pth`
    
    """
    # Open raster with xarray and make meshgrid
    da = xr.open_dataset(pth)
    
    if epsg is not None:
        da = da.rio.reproject(epsg)
        
    data = da['band_data'].data[band]
    x, y = np.meshgrid(da['x'], da['y'])

    # Compute distance to vent
    dist = np.sqrt(np.power((x-cX),2)+np.power((y-cY),2))

    return data, x, y, dist

def getEPSG(lat, lon):
    """
    Finds the EPSG code associated with a pair of lat, lon coordinates
    
        Args:
            lat (float): Latitude (decimal degrees). Negative in S hemisphere
            lon (float): Longitude (decimal degrees). Negative in W hemisphere

        Returns:
            int: The integer EPSG code
    
    """
    
    # Find the espg
    zone = int(np.ceil((lon+180)/6))
    if lat<0:
        crs = CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': True})
    else:
        crs = CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': False})
    epsg = crs.to_authority()
    
    return int(epsg[1])
# %%

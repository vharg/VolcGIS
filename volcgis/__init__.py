#%%
import pandas as pd
import geopandas as gpd
import numpy as np
import itertools
from shapely.geometry import Point, LineString, Polygon, box
import fiona
from fiona.crs import from_epsg
import osmnx as ox
from osgeo import gdal
import rasterio as rio
from rasterio.plot import show
from rasterio.mask import mask
from rasterio import features
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.transform import Affine
from rasterio.windows import from_bounds
import affine
from rasterio import shutil as rio_shutil
from rasterio.vrt import WarpedVRT

# import pycrs
# from rasterstats import zonal_stats
import pyproj
from pyproj import Transformer
import shapely

import os
import json
# os.chdir('/Users/seb/Documents/Codes/VolcGIS')

# Note:
# `Backup__init__.py`: Backup of the library. There might still be some good ideas in there though a lot were simplified

#%%
class eruption:
    def __init__(self, eruption, masterRes, outPath=None):
        ''' Sets main eruption object.
        
        Params:
            eruption (dict): Eruption dictionary
            masterRes (int): Cell size (m)
            outPath (str): Master folder for output path. If None, default is './volcanoes'
        '''

        self.outPath = 'volcanoes' if outPath == None else outPath
        
        # Set default variables
        self.name = eruption['name']
        self.res = masterRes

        # Set path and folders
        if not os.path.exists(os.path.join(self.outPath, self.name)):
            os.mkdir(os.path.join(self.outPath, self.name))
            os.mkdir(os.path.join(self.outPath, self.name, '_data'))
            os.mkdir(os.path.join(self.outPath, self.name, '_tmp'))
            os.mkdir(os.path.join(self.outPath, self.name, '_hazard'))
            os.mkdir(os.path.join(self.outPath, self.name, '_exposure'))
        
        # Define projections
        self.EPSG = eruption['epsg']
        # self.EPSG_geo = pyproj.CRS('EPSG:{}'.format(4326))
        # self.EPSG_proj = pyproj.CRS('EPSG:{}'.format(self.EPSG))
        self.EPSG_geo = 'epsg:{}'.format(4326)
        self.EPSG_proj = 'epsg:{}'.format(self.EPSG)

        # Convert vent into a geometry
        # xtmp, ytmp = pyproj.Transform.transform(self.EPSG_geo, self.EPSG_proj, eruption['vent'][0], eruption['vent'][1])
        # self.vent = {'lat': eruption['vent'][0], 'lon': eruption['vent'][1], 'easting': round(xtmp), 'northing': round(ytmp), 'alt': eruption['vent'][2]}
        # Update Pyproj2
        transformer = Transformer.from_crs(self.EPSG_geo, self.EPSG_proj)
        [xtmp, ytmp] = transformer.transform(eruption['vent'][0], eruption['vent'][1])
        self.vent = {'lat': eruption['vent'][0], 'lon': eruption['vent'][1], 'easting': round(xtmp), 'northing': round(ytmp), 'alt': eruption['vent'][2]}

        # Create area mask
        areaPl = box(eruption['extent'][0], 
                     eruption['extent'][2],
                     eruption['extent'][1],
                     eruption['extent'][3])
        
        # dxMin, dxMax, dyMin, dyMax = volcgis.makeZoneBoundaries(xtmp, ytmp, 

        self.area = gpd.GeoDataFrame({'geometry': areaPl}, index=[0], crs=self.EPSG_proj) # Area is now projected
        self.areaG = self.area.to_crs(from_epsg(4326)) # Project it

        # Define buffers
        tmp = pd.DataFrame(self.vent, index=[0])
        gdf = gpd.GeoDataFrame(tmp, geometry=gpd.points_from_xy(tmp.easting, tmp.northing), crs=self.EPSG_proj)
        
        buffer = pd.DataFrame()
        for iB in [10,30,100]:
            buffer = buffer.append(gdf.geometry.buffer(iB*1e3).rename(iB))
        buffer.columns = ['geometry']
        self.buffer = buffer
        
    
        # Define reference for virtual wrap on which all rasters will be aligned
        self.ref = {}
        self.ref['bounds'] = self.area.geometry[0].bounds
        self.ref['height'] = len(np.arange(int(self.ref['bounds'][1]), int(self.ref['bounds'][3]), self.res))
        self.ref['width'] = len(np.arange(int(self.ref['bounds'][0]), int(self.ref['bounds'][2]), self.res))
        self.ref['transform'] = affine.Affine(self.res, 0.0, self.ref['bounds'][0], 0.0, -self.res, self.ref['bounds'][3])
        self.ref['EPSG'] = rio.crs.CRS.from_epsg(self.EPSG)

    def getLandscan(self, inPath=None):
        """ Retrieves Landscan data for the area defined by self.area.

        Arguments:
            inPath (str): Path to corresponding raster

        Output:

        """
        if inPath is None:
            inPath = 'DATA/Landscan.tif'

        outPath = os.path.join(self.outPath, self.name, '_data', 'Landscan.tif')

        originalRes = 1000 # Landscan resolution
        scaling = originalRes / self.res # Scaling factor to correct population

        self.alignRaster(inPath, outPath, resampling='nearest', scalingFactor=scaling)

    def getLandcover(self, inPath=None):
        """ Retrieves Landcover data for the area defined by self.area. Resampling is set to 'nearest' for discrete data

        Arguments:
            inPath (str): Path to corresponding raster

        Output:

        """
        if inPath is None:
            inPath = 'DATA/LC100_2018_croped.tif'

        outPath = os.path.join(self.outPath, self.name, '_data', 'Landcover.tif')
        self.alignRaster(inPath, outPath, resampling='nearest')

    def getBuildingExposure(self, inPath=None, outPath=None):
        """ Retrieves building exposure from George's analysis for area defined by self.area.

        Arguments:
            inPath (str): Path to corresponding raster

        Output:

        """
        self.alignRaster(inPath, outPath, resampling='nearest')
        # if inPath is None:
            
    def getRoadNetwork(self, inPath=None):
        """

        Arguments:
            inPath (str): Path to the main roads dataset

        Returns:
            roads (feather): A feather file of the roads within BBox

        """
        
        # Re-project self.areaG to pseudo mercator 3857
        bbox = self.area.to_crs(from_epsg(3857))

        # Get bbox for the geometry
        if inPath is None:
            inPath = 'DATA/SEA_roads_criticality.gpkg'

        roads = gpd.read_file(inPath, bbox=bbox['geometry'])

        # Reproject to self.EPSG
        roads = roads.to_crs(crs="EPSG:{}".format(self.EPSG))

        # Save as a feather to outPath
        outPath = os.path.join(self.outPath, self.name, '_data', 'roads.feather')
        roads.to_feather(outPath)
         
    def prepareHazard(self, hazard):
        """ Prepare the hazard layers

            Loads the files for a given hazard defined as hazard['hazard'] and from hazard['nameConstructor'] from the hazards/ folder.

            Args:
                hazard (dict): Main hazard dictionary

            Returns:
                Geotif
        """

        print('Preparing hazard: {}'.format(hazard['hazard']))
        # Set path and folders
        targetDir = os.path.join(self.outPath, self.name, '_hazard', hazard['hazard'])
        if not os.path.exists(targetDir):
            os.mkdir(targetDir)

        # Create a list of file names
        hazard['files'] = self.makeInputFileName(hazard['rootDir'], hazard['nameConstructor'])

        # Align files
        for inPath in hazard['files']:
            outPath = inPath.replace(hazard['rootDir'], '')
            outPath = outPath.replace('.asc', '.tif')
            outPath = 'volcanoes/{}/_hazard/{}/{}'.format(self.name, hazard['hazard'], outPath)
            print('  - Processing: {}'.format(outPath))

            # If asc file, need to define projection
            # if inPath[-3:] == 'asc':
            #     inPath = asc2raster(inPath, hazard
            #                        )
            self.alignRaster(inPath, outPath, epsg=hazard['epsg'])

    def alignRaster(self, inPath, outPath, epsg=None, resampling='cubic', scalingFactor=None):
        """
            Aligns raster to a reference grid using a virtual wrap. Use reference contained in self.ref to read a window of original file and wrap it

            Args:
                inPath (str): Path to input raster
                outPath (str): Path to output raster
                resampling (str): Resampling method
                scalingFactor (float):

        """
        # Virtual wrapper options
        vrt_options = {
            'resampling': Resampling.cubic,
            'crs': self.ref['EPSG'],
            'transform': self.ref['transform'],
            'height': self.ref['height'],
            'width': self.ref['width'],
            'driver': 'GTiff'
        }

        if resampling == 'nearest':
            vrt_options['resampling'] = Resampling.nearest

        with rio.open(inPath, 'r') as src:
            # In case no EPSG is specified (e.g. asc)
            if src.crs == None:
                #  src.crs = 'EPSG:{}'.format(epsg)
                 src.crs = rio.crs.CRS.from_epsg(epsg)
            with WarpedVRT(src, **vrt_options) as vrt:
                rst = vrt.read(1, window=from_bounds(self.ref['bounds'][0], self.ref['bounds'][1], self.ref['bounds'][2], self.ref['bounds'][3], self.ref['transform']))
                # rst = vrt.read(1, window=from_bounds(self.ref['bounds'][0], self.ref['bounds'][1], self.ref['bounds'][2], self.ref['bounds'][3], src.transform))

                rio_shutil.copy(vrt, outPath, driver='GTiff', compress='lzw')

        # if scalingFactor is not None:
        #     with rio.open(outPath, 'w', **profile) as src:
        #         data = np.round(src.read(1)/(scalingFactor**2))
        #         src.write(data.astype(rio.int32))

    def makeInputFileName(self, rootDir, nameConstructor):
        ''' Create list of filenames based on nameConstructor '''
        flName = []
        for p in itertools.product(*list(nameConstructor.values())):
            fl = '_'.join(p)
            fl = rootDir+fl[:-5]+fl[-4:]
            flName.append(fl)
        return flName

def makeZoneBoundaries(x, y, xmin, xmax, ymin, ymax):
        """ Calculates the boundaries in m from a central point. Useful when
            a boundary falls outside the limit of the zone defined by the central
            point.

        Args:
            x (float): Easting of the central point
            y (float): Northing of the central point
            minx (float): Minimum Easting of the bounding box
            maxx (float): Maximum Easting of the bounding box
            miny (float): Minimum Northing of the bounding box
            maxy (float): Maximum Northing of the bounding box

        """

        if xmin<0:
            xmin = -1*(round(x))+xmin
        else:
            xmin = round(x)-xmin

        xmax = xmax - x

        if ymin<0:
            ymin = -1*(round(x))+ymin
        else:
            ymin = round(x)-ymin

        ymax = ymax - y

        return xmin, xmax, ymin, ymax
             

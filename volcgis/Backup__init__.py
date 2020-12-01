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

import pycrs
from rasterstats import zonal_stats
import pyproj
import shapely

import os
import json
os.chdir('/Users/seb/Documents/Codes/VolcGIS')

#%%
class eruption:
    # def __init__(self,
    #     name = None,
    #     VEI = None,
    #     vent = None, # List containing lat, lon, alt
    #     area = None, # List containing south, north, west, east
    #     espg = None
    #     ):
    def __init__(self, eruption, masterRes):
        ''' Sets main eruption object.
        
        Params:
            eruption (dict): Eruption dictionary
            masterRes (int): Cell size (m)
            
        '''
        
        # Set default variables
        self.name = eruption['name']
        self.res = masterRes
        # self.VEI = VEI

        # Set path and folders
        if not os.path.exists(os.path.join('volcanoes', self.name)):
            os.mkdir(os.path.join('volcanoes', self.name))
            os.mkdir(os.path.join('volcanoes', self.name, '_data'))
            os.mkdir(os.path.join('volcanoes', self.name, '_tmp'))
            os.mkdir(os.path.join('volcanoes', self.name, '_hazard'))

        # Define projections
        self.EPSG = eruption['epsg']
        self.EPSG_geo = pyproj.CRS('EPSG:{}'.format(4326))
        self.EPSG_proj = pyproj.CRS('EPSG:{}'.format(self.EPSG))

        # Convert vent into a geometry
        xtmp, ytmp = pyproj.transform(self.EPSG_geo, self.EPSG_proj, eruption['vent'][0], eruption['vent'][1])
        self.vent = {'lat': eruption['vent'][0], 'lon': eruption['vent'][1], 'easting': round(xtmp), 'northing': round(ytmp), 'alt': eruption['vent'][2]}

        # Create area mask
        # areaPl = box(eruption['area'][2], eruption['area'][0], eruption['area'][3], eruption['area'][1])
        # self.area = gpd.GeoDataFrame({'geometry': areaPl}, index=[0], crs=from_epsg(4326))
        # self.areaP = self.area.to_crs(self.EPSG_proj) # Project it
        areaPl = box(round(xtmp)-eruption['extent'][0], 
                     round(ytmp)-eruption['extent'][2],
                     round(xtmp)+eruption['extent'][1],
                     round(ytmp)+eruption['extent'][3])
                    #  eruption['area'][0], eruption['area'][3], eruption['area'][1])
        self.area = gpd.GeoDataFrame({'geometry': areaPl}, index=[0], crs=self.EPSG_proj) # Area is now projected
        self.areaG = self.area.to_crs(from_epsg(4326)) # Project it

        # Define reference for virtual wrap
        self.ref = {}
        self.ref['bounds'] = self.area.geometry[0].bounds
        self.ref['height'] = len(np.arange(int(self.ref['bounds'][1]), int(self.ref['bounds'][3]), self.res))
        self.ref['width'] = len(np.arange(int(self.ref['bounds'][0]), int(self.ref['bounds'][2]), self.res))
        self.ref['transform'] = affine.Affine(self.res, 0.0, self.ref['bounds'][0], 0.0, -self.res, self.ref['bounds'][3])
        self.ref['EPSG'] = rio.crs.CRS.from_epsg(self.EPSG)

        # # Datasets
        # self.Landscan = []
        # self.OSM_buildings = []
        # self.LC = []

        # # Hazards
        # self.tephra = []
        # self.PDC = []
        # self.BAF = []
        # self.lahars = []
        # self.clasts = []



    def getLandscan(self):
        """ Retrieves Landscan data for the area defined by self.area
            force:  Clip original even if _data/Landscan.tif exists already
        """
        
        inPth = 'DATA/Landscan.tif'
        outPth = os.path.join('volcanoes', self.name, '_data', 'Landscan.tif')
        self.alignRaster(inPth, outPth)

    def getLandcover(self):
        """ Retrieves Landcover data for the area defined by self.area
        """
        
        inPth = 'DATA/LC100_2018.tif'
        outPth = os.path.join('volcanoes', self.name, '_data', 'Landcover.tif')
        self.alignRaster(inPth, outPth)
        
    def getOSM(self):
        """ Retrieves OSM data for self.area using OSMnx

        """

        # https://osmnx.readthedocs.io/en/stable/osmnx.html#osmnx.footprints.footprints_from_polygon
        self.OSM_buildings = ox.footprints_from_polygon(self.area.loc[0]['geometry'])
        self.OSM_buildings = self.OSM_buildings[['geometry', 'name', 'nodes']] # 2020-04-29 13:22:19: For now keeping only these 3 columns to attempt saving time during saving, maybe we can see if there is anything more useful later on

        # 2020-04-29 10:17:06 I haven't found a way to save a gpd in a SpatiaLite... Saving into shapefiles takes a lot of time, so I think I will stick with loading it every time for now

        # Extract only one type of road:
        # https://stackoverflow.com/questions/52231122/python-osmnx-extract-only-big-freeways-of-a-country/52412274#52412274
    
        # Convert graph to gdf
        # https://automating-gis-processes.github.io/CSC18/lessons/L3/retrieve-osm-data.html
        # nodes, edges = ox.graph_to_gdfs(graph)
        # need to convert object types

    
    @staticmethod
    def getFeatures(gdf):
        """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
        return [json.loads(gdf.to_json())['features'][0]['geometry']]

    # Clip raster
    def clipRaster(self, inPth, outPth, epsg_code=None):
        """ Clip rasters according to self.area
            inPth:      Path to input raster (str)
            outPth:     Path to output raster (str)
            epsg_code:  EPSG code (e.g. 'EPSG:4326') (str). If none, assumes EPSG:4326
        """

        with rio.open(inPth) as data:
            out_img, out_transform = mask(data, self.getFeatures(self.area), crop=True)

            out_meta = data.meta.copy()
            # if epsg_code is None:
            #     epsg_code = int(data.crs.data['init'][5:])
            #     epsg_code = pycrs.parse.from_epsg_code(epsg_code).to_proj4()
            # else:
            #     epsg_code = data.crs
            

            out_meta.update({"driver": "GTiff",
                            "height": out_img.shape[1],
                            "width": out_img.shape[2],
                            "transform": out_transform})
            # ,
            #                 "crs": pycrs.parse.from_epsg_code(epsg_code).to_proj4()}
            #                         )

            with rio.open(outPth, "w", **out_meta) as dest:
                dest.write(out_img)

    def prepareHazard(self, hazard, resampling='cubic'):
        
        # Set path and folders
        targetDir = os.path.join('volcanoes', self.name, '_hazard', hazard['hazard'])
        if not os.path.exists(targetDir):
            os.mkdir(targetDir)

        # Create a list of file names
        hazard['files'] = self.makeInputFileName(hazard['rootDir'], hazard['nameConstructor'])
        
        # Align files
        for inPth in hazard['files']:
            outPth = inPth.replace(hazard['rootDir'], '')
            outPth = outPth.replace('.asc', '.tif')
            outPth = 'volcanoes/{}/_hazard/{}/{}'.format(self.name, hazard['hazard'], outPth)
            print(inPth)
            self.alignRaster(inPth, outPth)

      
    def alignRaster(self, inPth, outPth, resampling='cubic'):
        # Virtual wrapper options
        vrt_options = {
            'resampling': Resampling.cubic,
            'crs': self.ref['EPSG'],
            'transform': self.ref['transform'],
            'height': self.ref['height'],
            'width': self.ref['width'],
        }
        print(vrt_options)
        if resampling == 'nearest':
            vrt_options['resampling'] = Resampling.nearest
            
        with rio.open(inPth) as src:  
            with WarpedVRT(src, **vrt_options) as vrt:
                rst = vrt.read(1, window=from_bounds(self.ref['bounds'][0], self.ref['bounds'][1], self.ref['bounds'][2], self.ref['bounds'][3], self.ref['transform']))
                rio_shutil.copy(vrt, outPth, driver='GTiff', compress='lzw')    
       
        # with rio.open(inPth) as src:
        #     with WarpedVRT(src, **vrt_options) as vrt:
        #                 rst = vrt.read(1, window=from_bounds(erup.ref['bounds'][0], erup.ref['bounds'][1], erup.ref['bounds'][2], erup.ref['bounds'][3], src.transform))

        #         rio_shutil.copy(vrt, outPth, driver='GTiff')
                    
    def importHazard(self, hazard):
        """ Pre-process hazard data. 
            hazard:     Dictionary
        """
        
        # Initialize
        hazard['raster'] = []

        # Set path and folders
        targetDir = os.path.join(self.name, '_hazard', hazard['hazard'])
        if not os.path.exists(targetDir):
            os.mkdir(targetDir)
        # Remove temporary files
        for root,_,files in os.walk(self.name + '/_tmp/'):
            for fl in files:
                os.remove(os.path.join(root,fl))
        
        # Create a list of file names
        hazard['files'] = self.makeInputFileName(hazard['rootDir'], hazard['nameConstructor'])
        
        for i in range(0, len(hazard['files'])):
        # for i in range(1,2):

            # Read asc file and set CRS with Rasterio
            # this is the only way I found to properly set a crs to a .asc file
            print('Processing file:' + inPth)
            inPth = hazard['files'][i]
            outPth = '{}/_tmp/crs.tif'.format(self.name)

            # Setting up projection
            print('  - Setting up projection EPSG:{}'.format(hazard['epsg']))
            with rio.open(inPth, mode='r+') as tmp:
                tmp.crs = rio.crs.CRS.from_epsg(hazard['epsg'])

                out_meta = tmp.meta.copy()
                out_meta.update({"driver": "GTiff",
                                "height": tmp.height,
                                "width": tmp.width})

                # Good to know: use  .read() to access underlying data
                with rio.open(outPth, 'w', **out_meta) as dest:
                    dest.write(tmp.read())
 
            # Reproject to EPSG:4326 -> not anymore
            # Project to utm
            old_crs = 'EPSG:{}'.format(hazard['epsg'])

            if old_crs == 'EPSG:4326':
                print('  - Projecting from {} to {}'.format(old_crs, self.EPSG))
                inPth = outPth
                outPth = '{}/_tmp/projected.tif'.format(self.name)
                # new_crs = 'EPSG:4326'
                new_crs = 'EPSG:{}'.format(self.EPSG)
                self.reProject(inPth, outPth, old_crs, new_crs)
            else:
                new_crs = old_crs

            # Clip to region
            print('  - Clipping region extent')
            inPth = outPth
            outPth = '{}/_tmp/clip_region.tif'.format(self.name)
            self.clipRaster(inPth, outPth, epsg_code=new_crs)
            
            # # Mask and clip extent
            # inPth = outPth
            # outPth = '{}/_tmp/clip_extent.tif'.format(self.name)
            # self.clipExtent(inPth, outPth, val=0)
            
            # Resample
            tmp = rio.open(outPth, 'r')
            scale = tmp.res[0]/self.res
            print('  Resampling by factor ' + round(scale,3))
            outPth = hazard['files'][i]
            outPth = outPth.replace(hazard['rootDir'], '')
            outPth = outPth.replace('.asc', '.tif')
            outPth = '{}/_hazard/{}/{}'.format(self.name, hazard['hazard'], outPth)

            self.resampleRaster(tmp, outPth, scale)
            # hazard['raster'].append(rio.open(outPth))

        setattr(self, hazard['hazard'], hazard)

    def clipExtent(self, inPth, outPth, val=0):
        # Create mask from rasterio raster value:
        # https://stackoverflow.com/questions/37898113/python-boolean-array-to-polygon
        # Combine polygons all together:
        # https://stackoverflow.com/questions/40385782/make-a-union-of-polygons-in-geopandas-or-shapely-into-a-single-geometry
        # https://gis.stackexchange.com/questions/247118/crop-a-raster-using-rasterio-and-geopandas

        with rio.open(inPth) as src:
            # Metadata
            src_meta = src.meta.copy()
            data = src.read().astype(np.float32)
            msk = data>0
            # Extract features -> don't understand everything
            shapes = rio.features.shapes(data, mask=msk, transform=src.transform)
            # Convert to shape
            polygons = [shapely.geometry.Polygon(shape[0]["coordinates"][0]) for shape in shapes if shape[1] > val]
            # Convert to GPD and combone them
            boundary = gpd.GeoSeries(shapely.ops.cascaded_union(polygons))
            # Convert to json
            json_str = boundary.to_json()
            json_dict = json.loads(json_str)
            # Create mask
            msk_poly = [feature["geometry"] for feature in json_dict["features"]]
            # Mask and crop
            masked_band, masked_transform = mask(src, msk_poly, crop=True)
            # Update metadata
            src_meta.update(dtype=rio.float32, height=int(masked_band.shape[1]), width=int(masked_band.shape[2]), nodata=0, transform=masked_transform)
                # Write
            with rio.open(outPth, 'w', **src_meta) as dst:
                dst.write(masked_band.astype(rio.float32))  
                        
                
    def resampleRaster(self, raster, outPth, scale=2):
        # https://gis.stackexchange.com/questions/329945/should-resampling-downsampling-a-raster-using-rasterio-cause-the-coordinates-t
        t = raster.transform

        # rescale the metadata
        transform = Affine(t.a / scale, t.b, t.c, t.d, t.e / scale, t.f) # <== division
        height = raster.height * scale                                   # <== multiplication
        width = raster.width * scale

        profile = raster.profile
        profile.update(transform=transform, driver='GTiff', height=height, width=width, crs=raster.crs)

        data = raster.read( # Note changed order of indexes, arrays are band, row, col order not row, col, band
                out_shape=(int(raster.count), int(height), int(width)),
                resampling=Resampling.cubic)
        data[data<0] = 0
        data[data>1.1] = 0

        with rio.open(outPth,'w', **profile) as dst:
            dst.write(data)
            # yield data

    def reProject(self, inpath, outpath, old_crs, new_crs):
        """Reproject raster 
           https://www.earthdatascience.org/courses/use-data-open-source-python/intro-raster-data-python/raster-data-processing/reproject-raster/
           
           """
        # dst_crs = rio.crs.CRS.from_epsg(new_crs) # CRS for web meractor 
        # old_crs = rio.crs.CRS.from_epsg(old_crs)
        dst_crs = new_crs
        # print(new_crs)
        # print(old_crs)

        with rio.open(inpath) as src:
            # transform, width, height = calculate_default_transform(
            #     src.crs, dst_crs, src.width, src.height, *src.bounds)
            print(src.bounds)
            transform, width, height = calculate_default_transform(
                old_crs, dst_crs, src.width, src.height, *src.bounds)
            kwargs = src.meta.copy()
            kwargs.update({
                'crs': dst_crs,
                'transform': transform,
                'width': width,
                'height': height
            })

            with rio.open(outpath, 'w', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rio.band(src, i),
                        destination=rio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=old_crs,
                        dst_transform=transform,
                        dst_crs=dst_crs,
                        resampling=Resampling.nearest)


    def makeInputFileName(self, rootDir, nameConstructor):
        ''' Create list of filenames based on nameCons '''
        flName = []
        for p in itertools.product(*list(nameConstructor.values())):
            fl = '_'.join(p)
            fl = rootDir+fl[:-5]+fl[-4:]
            flName.append(fl)
        return flName
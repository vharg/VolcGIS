#%%
import os
os.chdir('/Users/seb/Documents/Codes/VolcGIS')
import volcgis 
from rasterio.plot import show
import rasterio as rio
# import pycrs

# name = 'sakurajima'
# VEI = 5
# vent = [31.593, 130.657, 1117]
# area = [31.5, 31.8, 130.5, 130.725]
# # Setup eruption
# erup = volcgis.eruption(VEI=VEI, vent=vent, area=area, name='SK')

# All coordinates must be EPSG:4326

# Create an eruption
# name = 'gede' # Name
# VEI = 5 # VEI (unused for now)
# vent = [-6.787948, 106.981215, 2765] # Vent lat, lon and elevation
# area = [-7, -6.55, 106.7, 107.2] # S N W E
# epsg = 

eruption = {
    'name':     'Gede-Pangrango',
    'vent':     [-6.787948, 106.981215, 2765], # Vent lat, lon and elevation
    'area':     [-7, -6.55, 106.7, 107.2], # S N W E
    'extent':   [1.e5, 1.e5, 1.e5, 1.e5], # extent around vent in meters, [minx maxx miny maxy]
    'epsg':     32748
}



# Setup eruption
erup = volcgis.eruption(eruption, 90)


#%%
# Get Landscan data
# erup.getLandscan(force=True) # This is not accessible for now as the data (3 Gb) is stored locally

# Retrieve OSM data
# erup.getOSM()
# For now, this populates erup.OSM_buildings which is a gdf containing the buildings. It is really big


# YANG: Maybe you don't need that for now
# Set tephra hazard
nameConstructor = {
    'volcano':    ['Gede-Pangrango'],
    'VEI':        ['VEI3', 'VEI4', 'VEI5'],
    'perc':       ['P1', 'P5', 'P9'],
    'format':     ['.tif']
}
tephra = {
    'hazard':   'Tephra',
    'epsg':     32748,
    'kind':     'him',
    'unit':     'kg/m2',
    'label':    [''],
    'rootDir':   'hazards/Tephra/',
    'nameConstructor': nameConstructor,
    }

# erup.importHazard(tephra)
erup.prepareHazard(tephra)

#%% BAF
nameConstructor = {
    'volcano':    ['Gede-Pangrango'],
    'buffer':     ['0.0027deg', '0.009deg'],
    'volume':     ['9800000m3', '450000m3'],
    'format':     ['.tif']
}

BAF = {
    'hazard':   'BAF',
    'epsg':     32748,
    'kind':     'prob',
    'unit':     '%',
    'label':    [''],
    'rootDir':    'hazards/BAF/',
    'nameConstructor': nameConstructor,
    }  

# erup.importHazard(BAF)
erup.prepareHazard(BAF)



#%%
# PDC
nameConstructor = {
    'volcano':    ['Gede-Pangrango'],
    'VEI':        ['3', '4', '5'],
    'suffix':     ['output_map'],
    'format':     ['.asc']
}

PDC = {
    'hazard':   'PDC',
    'epsg':     32748,
    'kind':     'prob',
    'unit':     '%',
    'rootDir':    'hazards/PDC/',
    'nameConstructor': nameConstructor,
    }  
# erup.importHazard(PDC)
erup.prepareHazard(PDC)

#%%
# Large clast
nameConstructor = {
    'volcano':    ['Gede-Pangrango'],
    'VEI':        ['VEI3', 'VEI4', 'VEI5'],
    'format':     ['.tif']
}

LC = {
    'hazard':   'LC',
    'epsg':     4326,
    'kind':     'prob',
    'unit':     '%',
    'rootDir':  'hazards/LC/',
    'nameConstructor': nameConstructor,
    }  
# erup.importHazard(LC)
erup.prepareHazard(LC)

#%%
erup.getLandcover()
erup.getLandscan()



#%%
# How to do zonal stats
zone = zonal_stats(erup.OSM_buildings, erup.Landscan.read(1), affine=erup.Landscan.transform, stats=['mean'], all_touched=True)

zone = zonal_stats(erup.OSM_buildings, erup.tephra['raster'][0].read(1), affine=erup.tephra['raster'][0].transform, stats=['mean'], all_touched=True, geojson_out=True)
df = pd.DataFrame(zone)

testB = erup.OSM_buildings.iloc[1:100][['geometry', 'name', 'nodes']]

# Zonal stat
zone = zonal_stats(testB, erup.tephra['raster'][0].read(1), affine=erup.tephra['raster'][0].transform, stats=['mean'], all_touched=True, geojson_out=True)
# Get GeoJSON output into dataframe
df = pd.DataFrame(zone)
# Extract the mean value
df['mean'] = df['properties'].apply(lambda x: x.get('mean'))
# Set columns, column types and id
df = df[['id','mean']].astype({'id': 'int64', 'mean': 'float64'}).rename({'mean':'hazard_mean'}, axis=1).set_index('id')
# Join datasets
testB = testB.join(df)



## WAYS

edgesBackup = edges.copy()
edges = edgesBackup.copy()

edges = edges[['osmid', 'length', 'highway', 'geometry']]

# YES!!!!!
# Create a polygon from a value in a rasterio raster
# https://stackoverflow.com/questions/37898113/python-boolean-array-to-polygon


tmpRaster = erup.tephra['raster'][0] # Raster
tmpData = tmpRaster.read(1)   # Array

# Manually mask the data
tmpData[tmpData>=1] =1
tmpData[tmpData<1] = 0

# Mask the data and export the polygon 
shapes = rio.features.shapes(tmpData,transform = tmpRaster.transform)
polygons = [shapely.geometry.Polygon(shape[0]["coordinates"][0]) for shape in shapes if shape[1] == 1]

# This is assuming there is only one
buildClip = gpd.clip(erup.OSM_buildings, polygons[0])







## Coordinate changes
pyproj.transform(wgs84, isn2004, 63.983, -19.700)     





EPSG = eruption['epsg']
EPSG_geo = pyproj.CRS('EPSG:{}'.format(4326))
EPSG_proj = pyproj.CRS('EPSG:{}'.format(EPSG))

xtmp, ytmp = pyproj.transform(EPSG_geo, EPSG_proj, 17, 20)


#%%

refPth = 'gede/_hazard/PDC_prob_.tif'
ref = rio.open(refPth)

targetPth = 'gede_data/IM_dry_25_kgm2.txt'
target = rio.open(targetPth)


# %% Test clipping non-zero extent

img = rio.open('Gede-Pangrango/_hazard/BAF/Gede-Pangrango_0.009deg_450000m3.tif')

# %%

# Read raster
img = rio.open('test.tif')
# Get values
data = img.read().astype(np.float32)
# Create a mask
# data[data>0]=255
# https://gis.stackexchange.com/questions/268331/how-can-i-extract-contours-from-a-raster-with-python



msk = data>0
# shapes = rio.features.shapes(data, mask=mask)
shapes = rio.features.shapes(data, mask=msk, transform=img.transform)
pprint.pprint(next(shapes))

geom_list = []
for geom, value in shapes:
    # Print GeoJSON shapes to stdout
    geom = shapely.geometry.shape(geom)
    geom_list.append(geom)


out_img, out_transform = mask(img, shapes, crop=True)


#%%
##### WINNER
# Create mask from rasterio raster value:
# https://stackoverflow.com/questions/37898113/python-boolean-array-to-polygon
# Combine polygons all together:
# https://stackoverflow.com/questions/40385782/make-a-union-of-polygons-in-geopandas-or-shapely-into-a-single-geometry
# https://gis.stackexchange.com/questions/247118/crop-a-raster-using-rasterio-and-geopandas

with rio.open('test.tif') as src:
    # Metadata
    src_meta = src.meta.copy()
    data = src.read().astype(np.float32)
    # Extract features -> don't understand everything
    shapes = rio.features.shapes(data, transform=src.transform)
    # Convert to shape
    polygons = [shapely.geometry.Polygon(shape[0]["coordinates"][0]) for shape in shapes if shape[1] > 0]
    # Convert to GPD and combone them
    boundary = gpd.GeoSeries(shapely.ops.cascaded_union(polygons))
    # Convert to json
    json_str = boundary.to_json()
    json_dict = json.loads(json_str)
    # Create mask
    msk_poly = [feature["geometry"] for feature in json_dict["features"]]
    # Mask and crop
    masked_band, masked_transform = mask.mask(src, msk_poly, crop=True)
    # Update metadata
    src_meta.update(dtype=rio.float32, height=int(masked_band.shape[1]), width=int(masked_band.shape[2]), nodata=0, transform=masked_transform, compress='lzw')
        # Write
        with rio.open('test2.tif', 'w', **src_meta) as dst:
                dst.write(masked_band.astype(rio.float32))  

####
 
 
# %%
exte
# %%
contour= measure.find_contours(img.read(1),.1)
ct = rio.features.shapes(contour, transform=img.transform) #, mask=None, connectivity=4, transform=Affine(1.0, 0.0, 0.0, 0.0, 1.0, 0.0))

geom_list = []
for geom, value in ct:
    # Print GeoJSON shapes to stdout
    geom = shapely.geometry.shape(geom)
    geom_list.append(geom)

# rasterio.features.shapes(source, mask=None, connectivity=4, transform=Affine(1.0, 0.0, 0.0, 0.0, 1.0, 0.0))
# %%



#%%
import rasterio
from rasterio import features
from rasterio import mask
import json
import geopandas as gpd
import os

#%%
results = []
final_results = []

with rio.open("test.tif") as src:
    src_meta = src.meta.copy()
    src_affine = src_meta.get("transform")

    band = src.read(1).astype(np.float32)

    for geometry, raster_value in features.shapes(band, transform=src_affine):
        if raster_value > 0:
            result = {'properties': {'raster_value': raster_value}, 'geometry': geometry}
            results.append(result)

        gpd_results_filtered = gpd.GeoDataFrame.from_features(results)

        # gpd_results["area"] = gpd_results["geometry"].area

        # gpd_results_filtered = gpd_results[gpd_results["area"] > 5000]

        gpd_filtered_json_str = gpd_results_filtered.to_json()

        gpd_filtered_json_dict = json.loads(gpd_filtered_json_str)

        final_results = [feature["geometry"] for feature in gpd_filtered_json_dict["features"]]

        masked_band, masked_transform = mask.mask(src, final_results, invert=True)

        masked_band[masked_band > 255] = 0

        src_meta.update(dtype=rio.float32, height=int(masked_band.shape[1]), width=int(masked_band.shape[2]), nodata=0, transform=masked_transform, compress='lzw')

        with rio.open('test2.tid', 'w', **src_meta) as dst:
                dst.write(masked_band.astype(rio.float32))   
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
# %% Test virtual wrap

ref = {}
ref['bounds'] = erup.area.geometry[0].bounds
ref['height'] = len(np.arange(int(ref['bounds'][1]), int(ref['bounds'][3]), erup.res))
ref['width'] = len(np.arange(int(ref['bounds'][0]), int(ref['bounds'][2]), erup.res))
ref['transform'] = affine.Affine(erup.res, 0.0, ref['bounds'][0], 0.0, -erup.res, ref['bounds'][3])
ref['EPSG'] = CRS.from_epsg(erup.EPSG)

# %% Testing windowing and virtual wrap
# inPth = 'DATA/Landscan.tif'
inPth = 'hazards/BAF/Gede-Pangrango_0.009deg_450000m3.tif'
outPth = 'testLS.tif'
# bounds = erup.areaG.geometry[0].bounds

vrt_options = {
    'resampling': Resampling.cubic,
    'crs': erup.ref['EPSG'],
    'transform': erup.ref['transform'],
    'height': erup.ref['height'],
    'width': erup.ref['width'],
}
        
with rio.open(inPth) as src:  
    with WarpedVRT(src, **vrt_options) as vrt:
        rst = vrt.read(1, window=from_bounds(erup.ref['bounds'][0], erup.ref['bounds'][1], erup.ref['bounds'][2], erup.ref['bounds'][3], erup.ref['transform']))
        rio_shutil.copy(vrt, outPth, driver='GTiff', compress='lzw')
        
    # rst = src.read(1, window=from_bounds(bounds[0], bounds[1], bounds[2], bounds[3], src.transform))
# %%

# %%
EPSG = eruption['epsg']
EPSG_geo = pyproj.CRS('EPSG:{}'.format(4326))
EPSG_proj = pyproj.CRS('EPSG:{}'.format(EPSG))

# Convert vent into a geometry
xtmp, ytmp = pyproj.Transform.transform(EPSG_geo, EPSG_proj, eruption['vent'][0], eruption['vent'][1])
vent = {'lat': eruption['vent'][0], 'lon': eruption['vent'][1], 'easting': round(xtmp), 'northing': round(ytmp), 'alt': eruption['vent'][2]}

from pyproj import Transformer
transformer = Transformer.from_crs("epsg:4326", "epsg:{}".format(EPSG))
[x,y] = transformer.transform(12, 12)

transformer = Transformer.from_crs("epsg:4326", "epsg:{}".format(EPSG))
[xtmp, ytmp] = pyproj.Transform.transform(self.EPSG_geo, self.EPSG_proj, eruption['vent'][0], eruption['vent'][1])


#%%
[xtmp, ytmp] = transformer.transform(eruption['vent'][0], eruption['vent'][1])
self.vent = {'lat': eruption['vent'][0], 'lon': eruption['vent'][1], 'easting': round(xtmp), 'northing': round(ytmp), 'alt': eruption['vent'][2]}

# Create area mask
areaPl = box(round(xtmp)-eruption['extent'][0], 
                round(ytmp)-eruption['extent'][2],
                round(xtmp)+eruption['extent'][1],
                round(ytmp)+eruption['extent'][3])






#%%

with rio.open('hazards/LC/Agung_VEI3.tif') as src:
    print(src.crs)
# %%
with rio.open('hazards/PDC/Agung_3_output_map.asc') as src:
    print(src.crs)

# %%

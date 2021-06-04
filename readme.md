
## Install

Create environment and set the channel to `conda-forge`:

```

conda config --env --add channels conda-forge
conda config --env --set channel_priority strict
```

Then:
```
conda install -c conda-forge pyarrow rioxarray rasterio geopandas bokeh contextily osmnx
conda install -c conda-forge holoviews datashader panel param geoviews
```

Install these packages with `pip`:
```
pip install utm
pip install alive-progress
```

### Content
- `csv/`: Root of csv files
- `plottingApp`: Bokeh function used on WOVODAT

## Basic commands 

### Import
```python
import volcgis
```

### Defining and initializing an eruption

An eruption is defined with a dictionary, e.g.:

```python
eruption = {
    'name':     'Gede-Pangrango',
    'vent':     [-6.787948, 106.981215, 2765], # Vent lat, lon and elevation
    'extent':   [1.e5, 1.e5, 1.e5, 1.e5], # extent around vent in meters, [minx maxx miny maxy]
    'epsg':     32748
}
```

Where:
- `name`: eruption/volcano name (str)
- `vent`: geographic coordinates of the vent defined as [`lat`, `lon`, `alt`], where:
  - `lat`: Latitude in decimal degrees (WGS84, negative in S hemisphere)
  - `lon`: Longitude in decimal degrees (WGS84, negative in W hemisphere)
  - `alt`: Vent elevation (m asl)
- `extent`: Extent of the reference region, defined as [`dxMin`, `dxMax`, `dyMin`, `dyMax`]:
  - `dxMin`: Distance (m) from the vent in the W direction
  - `dxMax`: Distance (m) from the vent in the E direction 
  - `dyMin`: Distance (m) from the vent in the S direction 
  - `dyMax`: Distance (m) from the vent in the N direction 
- `epsg`: Digits of the projection to use (e.g. )

An eruption is then defined as:

```python
erup = volcgis.eruption(eruption, res)
```

Where:
- `res` is the reference resolution (m)

### Defining and pre-processing a hazard type

A hazard first requires the definition of a `nameConstructor` to retrieve hazard model outputs, e.g.:

```python
nameConstructor = {
    'volcano':    ['Gede-Pangrango'],
    'VEI':        ['VEI3', 'VEI4', 'VEI5'],
    'perc':       ['P1', 'P5', 'P9'],
    'format':     ['.tif']
}
```

`nameConstructor` is used to reconstruct the file names in a modular way, where each subcomponent is separated by an underscore `_`. In this particular case, `nameConstructor` is designed to retrieve all these files:

```
Gede-Pangrango_VEI3_P1.tif
Gede-Pangrango_VEI3_P5.tif
Gede-Pangrango_VEI3_P9.tif
Gede-Pangrango_VEI4_P1.tif
Gede-Pangrango_VEI4_P5.tif
Gede-Pangrango_VEI4_P9.tif
Gede-Pangrango_VEI5_P1.tif
Gede-Pangrango_VEI5_P5.tif
Gede-Pangrango_VEI5_P9.tif
```

The hazard is the defined with another dictionary, e.g.:

```python
tephra = {
    'hazard':    'Tephra',
    'epsg':      32748,
    'rootDir':   'hazards/Tephra/',
    'nameConstructor': nameConstructor,
    }
```

Where:
- `hazard`: Type of hazard (arbitrary but needs to be consistent)
- `epsg`: Digits of the coordinate system of the hazard model output
- `rootDir`: Location of the hazard model output files
- `nameConstructor`: See above

#### Preprocessing hazard files

The preprocessing of hazard files is the started with:

```python
erup.prepareHazard(tephra)
```

Essentially, this step:
1. Read each file found by `nameConstructor` in the `rootDir` folder
2. Use a virtual wrapper to **reproject** the hazard file to `erup.epsg`, **clip** the extent and **aling** pixels to the extent defined by `erup.area`
3. Saves the file in `volcanoes/___volcanoName/_hazard/___hazardName/`, where:
    - `___volcanoName` is `erup.name`
    - `___hazardName` is `hazard.hazard`

### Retrieving exposure data
Reproject, crop and align extent based on `erup.area`/`erup.epsg` and saves subset in `volcanoes/___volcanoName/_data/`. In all cases, a custom path to the source dataset can be specified with the named keyword argument `inPth`.

#### Landscan

```python
erup.getLandscan()
```

**IMPORTANT:** The strategy used to resample Landscan data is to use a *nearest neighbour* algorithm so as to not change pixel values and then to *divide* the pixel values by the ratio of original/new resolutions. However I haven't found a way to do the second step in the whole process, so this has to be done manually before analysis. For instance:

```python
# Path to Landscan file
popf = 'volcanoes/{}/_data/Landscan.tif'.format(erup.name)
# Open file
pop = rio.open(popf)
# Read data and normalised by square of the ratio of original to resample resolutions
pop_data = pop.read(1)/(1000/erup.res)**2
# For population count, remove noData values (i.e. -9999)
pop_data[pop_data<1] = 0
# Close original raster
pop.close()
```

#### Landcover

```python
erup.getLandcover()
```


## Tips

### Plotting geopandas basemaps with contextily

It is possible to plot basemaps to geopandas with [contextily](https://contextily.readthedocs.io). Data needs to be in Web Mercator `EPSG:3875`:

```python
import contextily as ctx
fig = plt.figure(figsize=[8,10])
ax = fig.add_subplot(1, 1, 1)
ax = erup.areaG.to_crs('EPSG:3857').plot(alpha=0.5,ax=ax)
ctx.add_basemap(ax)
```

Alternatively, it is possible to re-project the basemap to `EPSG:4326`:

```python
ax = erup.areaG.plot(alpha=0.5)
ctx.add_basemap(ax, crs=erup.areaG.crs.to_string())
```

### Join the RSDS data to the original .feather file

```python
import geopandas as gpd

name = 'Taal'

roadf = os.path.join('volcanoes', name, '_data/roads.feather')
road = gpd.read_feather(roadf)
road['Road_ID'] = road['Road_ID'].astype(int) # Make sure the ID is in the same type
road = road.set_index('Road_ID')

rsdsf = os.path.join('volcanoes', name, '_exposure/RSDS.csv')
rsds = pd.read_csv(rsdsf)
rsds = rsds.set_index('Road_ID')

ROAD = road.join(rsds)
ROAD = ROAD.fillna(0)
```

### Starting the plotting app on Wovodat

Start the bokeh plotting app in `exposure.py`, making it available at [gee.wovodat.org/exposure](gee.wovodat.org/exposure)

```bash
conda activate ee
nohup bokeh serve --allow-websocket-origin='*' exposure.py --log-level=debug
```

## Change log

### 2021-03-04

- Moved all exposure analysis functions to `volcgis.exposureAnalysis`

#### Road Exposure
Updated road exposure from Josh's commit. Namely:
- Removed `ROAD_EXPOSURE` variable and added content to `EXPOSURE`
- Remove `updateRoads` and added to `updateExposure`

#### `getRNDS`
- `getRNDS` now returns three arguments, including `roadLength` and `RSDS`
- `roadLength` is added to `EXPOSURE`
- `RSDS` is a pd.DataFrame that contains `RSDS` value for each road segment defined by `Road_ID`. Each column is a different hazard occurrence (i.e. hazard type, VEI, prob etc)

#### `Gekko`
- Prepared `processHazard.py` for parallel processing on Gekko using Job Arrays
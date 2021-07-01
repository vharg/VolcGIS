# %% [markdown]
# # VolcGIS Demo
# 
# **VolcGIS** is a Python library designed to facilitate exposure analyses to volcanic hazards. This notebook illustrates its use taking Merapi volcano as a case study, considering the probabilistic hazard assessments for tephra using [TephraProb](https://github.com/e5k/TephraProb) and for PDC using https://github.com/AlvaroAravena/ECMapProb.
#
# > Biass, S., Bonadonna, C., Connor, L., Connor, C., 2016. TephraProb: a Matlab package for probabilistic hazard assessments of tephra fallout. J. Appl. Volcanol. 5, 1–16. doi:10.1186/s13617-016-0050-5
#
# > Aravena, A., Cioni, R., Bevilacqua, A., de' Michieli Vitturi, M., Esposti Ongaro, T., and Neri, A. (2020). Tree‐branching based enhancement of kinetic energy models for reproducing channelization processes of pyroclastic density currents. Journal of Geophysical Research: Solid Earth, 125, e2019JB019271.
# %%
import volcgis
from volcgis import eruption
import seaborn as sns

# %% [markdown]
# ## Setting up eruption/volcano
# We start by defining the volcano, get the EPSG code of the target UTM zone centered on the vent and et the resolution that will be used to resample all rasters:

# %%
# Define volcano
name = 'Merapi'
lat = -7.541
lon = 110.446

#Get the EPSG of the target UTM zone centered on the vent
epsg = volcgis.getEPSG(lat, lon)

# Set the resolution that will be used to resample all rasters
res = 250

# %% [markdown]
# We now define an eruption/volcano object and plot the study area. This will create a folder in `volcanoes/Merapi/` that contains the following folders:
# - `_data`: Will contain the processed datasets (i.e. re-projected, clipped and aligned)
# - `_exposure`: Will contain the exposure data
# - `_hazard`: Will contain the processed hazard datasets (i.e. geotifs, re-projected, clipped and aligned)
# - `_tmp`: Temporary files

# %%
# Setup eruption dictionary. The extent defines the coverage of the study area. It can be specified either in *absolute* UTM coordinates or as a *relative* distance from the vent. The behaviour is controlled below.
eruption = {
    'name':     name,
    'vent':     [lat, lon, 2765], 
    'extent':   [5e4, 5e4, 4e4, 4e4],
    'epsg':     epsg
}

# Set a dictionary with the paths to i) the output folder (i.e. `outPath/eruption['name']`) and ii) all the relevant exposure datasets
path = {
    'outPath': 'volcanoes',
    'populationPath': 'DATA/Landscan.tif',
    'landcoverPath': 'DATA/LC100_2018_croped.tif',
    'OSMroadsPath': 'DATA/indonesia-latest.osm.pbf'
}
# Define the main eruption/volcano object. Here, the extent is set as relative, meaning that our study covers 50 km in E-W and 40 km in N-S.
# `buffer` computes a circular polygon at the various distances (km), which is useful for assessing the exposure at concentric circles.
E = volcgis.eruption.eruption(eruption, res, path, extent='relative', buffer=[10,25,40])

# We can now plot the study area and, as an option, add the concentric circles with `plotBuffer`. In this case, it is obvious that the 100 km buffer is outside the study area - that is something one should pay attention to.
# E.plot(plotBuffer=True)

#%% [markdown]
# ### Note on the use of OSM data
# In order to work, the OSM data requires an additional dependency called OSMOSIS and used to clip pre-processed OSM data in `.pbf` format. Before using OSM data, make sure to:
# 1. Download the `osm.pbf` file for your region/country from [Geofabrik](http://download.geofabrik.de) (i.e. the file to which `OSMroadsPath` must point)
# 2. Install [Osmosis](https://wiki.openstreetmap.org/wiki/Osmosis/Installation)

# %% [markdown]
# ## Setting up hazards
# This step assumes that all hazard models have been run and are stored into readable raster formats (e.g. Geotif, ascii rasters).
# 
# We first need a `nameConstructor` dictionary to retrieve hazard model outputs. This is used to reconstruct the file names in a modular way, where each subcomponent is separated by an underscore `_`. In the case of the tephra example below, it will gather all these files:
# 
# ```
#  Merapi_VEI3_P50.tif
#  Merapi_VEI3_P90.tif
#  Merapi_VEI4_P50.tif
#  Merapi_VEI4_P90.tif
# ```

# %%
# Tephra hazard
nameConstructorTephra = {
    'volcano':    ['Merapi'],
    'VEI':        ['VEI3', 'VEI4'],
    'prob':       ['P50', 'P90'],
    'format':     ['.tif']
}
# PDC hazard
nameConstructorPDC = {
    'volcano':    ['Merapi'],
    'VEI':        ['3', '4'],
    'suffix':     ['output_map'],
    'format':     ['.asc']
}

# %% [markdown]
# We now define a hazard dictionary, specifying:
# - `hazard` as the name
# - `epsg` as the int code of the original coordinate system of the hazard files and 
# - `roodDir` as the directory where the files are located

#%%
tephra = {
    'hazard':   'Tephra',
    'epsg':     epsg,
    'rootDir':  'hazards/Tephra/',
    'nameConstructor': nameConstructorTephra,
    }

PDC = {
    'hazard':   'PDC',
    'epsg':     epsg,
    'rootDir':  'hazards/PDC/',
    'nameConstructor': nameConstructorPDC,
    }  


# %% [markdown]
# ### Pre-processing hazard files
# Essentially, this step:
# 
# 1. Reads each file found by `nameConstructor` in the `rootDir` folder
# 2. Uses a virtual wrapper to reproject the hazard file to `E.epsg`, clip the extent and align pixels to the extent defined by `E.area`
# 3. Saves the file in `E.path['outPath']`/`E.name`/_hazard/`tephra['hazard']`
# 
# Use `noAlign=True` to skip the align part - which should only be used in debug mode. 

# %%
E.prepareHazard(tephra, noAlign=True)
E.prepareHazard(PDC, noAlign=True)

# %% [markdown]
# ### Preparing hazard data for exposure analyses
# In addition, this step populates `E.hazards['hazardName']` with different variables, one useful one being `E.hazards['hazardName']['data']`, which will become handy during exposure analyses. It is stored as a `pd.DataFrame()`, which we need to clean a bit. For instance, for the `Tephra` example above:
# 
# 1. Only keep the digit of the `VEI` column
# 2. In this case the `perc` column represents a percentile rank, which we want to convert in a survivor function
# 3. Dropping the column `format`

# %%
tmp = E.hazards['Tephra']['data']
tmp['VEI'] = tmp['VEI'].str.extract('(\d+)').astype(int)
tmp['prob'] = 100-tmp['prob'].str.extract('(\d+)').astype(int)
tmp = tmp.drop(['format'], axis=1)
E.hazards['Tephra']['data'] = tmp

#%% [markdown]
# We do the same for the PDC hazard, where we just need to get rid of the `suffix` and `format` columns. Although the `VEI` column only contains a number, it is stored as a string, which we need to convert to an `int`.
#%%
tmp = E.hazards['PDC']['data']
print(tmp.dtypes) # See that VEI is stored as an object
tmp['VEI'] = tmp['VEI'].astype(int)
E.hazards['PDC']['data'] = E.hazards['PDC']['data'].drop(['suffix', 'format'], axis=1)

# %% [markdown]
# ### Plotting hazard data
# We can now plot the hazard on a map. `plotHazard` controls type of hazard to plot - which must correspond to an entry in `E.hazards`. `plotHazard` requires two additional arguments: `hazLevels` - which is a list containing the raster values to contour and `hazProps` - which is a dictionary that will control what hazard file to plot. Each pair of key/value must correspond to the columns defined in `E.hazards[plotHazard]` and return a unique value.

# %%
# E.plot(plotHazard='Tephra', hazLevels=[1,10,100,300], hazProps={'VEI': 5,'prob': 90})
# E.plot(plotHazard='PDC', hazLevels=[0.1, 0.5, 0.9], hazProps={'VEI': 4})

# %% [markdown]
# #### Knowing the data
# Contours defined in `hazLevels` must reflect the value of the underlying raster. In case of doubt, the argument `plotContours=False` can become useful to plot the raster only.

# %% [markdown]
# ## Prepare exposure data
# This step reprojects global exposure datasets to `E.epsg`, clip the extent and align pixels to the extent defined by `E.area`. In the case of population data, the population count is scaled by a ratio of `originalRasterResolution`/`E.res`. Processed files are saved in `E.path['outPath']`/`E.name`/_data/. 
# Setting `noOSM=True` is for the purpose of this notebook only, please remove.

# %%
E.prepareExposure(noOSM=True)

# %% [markdown]
# Let's plot the population dataset and overlay buffers:

#%%
E.plot(plotExposure='pop', plotBuffer=True)

# %% [markdown]
# Let's plot the road network. By default, `Local` roads are not displayed to avoid memory issues. It is possible to customize the appearance of the road network, so check the documentation.

#%%
E.plot(plotExposure='roads')

# %% [markdown]
# ## Get exposure data
# ### Exposure from concentric circles
# By default, `E.prepareExposure()` processes population and landcover files, which updates the property `self.exposure['expType']`. `E.getBufferExposure()` reads what exposure datasets have been defined and analyse the exposure for the buffers defined in `E.buffer`.
# 
# Here, we use the COPERNICUS CGLS-100 Landcover dataset to assess the exposure to crops (`LC==40`) and urban (`LC==50`) areas. To specify these, we first create a dictionary:

# %%
LC = {
    'crops':40,
    'urban':50
    }

# %%
# We then feed it to `E.getBufferExposure()`:
E.getBufferExposure(LC=LC)

# %% [markdown]
# `E.exposure['buffer']` now contains the exposure information as a function of concentric circles. We can plot the population exposure as a distance from the vent:

# %%
sns.barplot(x=E.exposure['buffer'].columns, y=E.exposure['buffer'].loc['pop_count'], color='coral')

# %% [markdown]
# ## Exposure from hazard footprints
# Similarly, we can get the exposure associated with the hazard footprints - here those in `E.hazards['Tephra']['data']`. In this particular example, each file contains *tephra accumulations* extracted at a given percentile. So whilst looping through the files will vary `VEI` and `prob`, we also want - in each file - to assess the exposure associated with various thresholds of tephra accumulations. To achieve that, we setup a `hazardProps` dictionary, where `columns` are the columns in `E.hazards['Tephra']['data']` that contain the relevant information. `varName` is just a variable name for the various variable values defined in `varVal`.

# %%
TephraProps= {
    'columns':  ['VEI', 'prob'],
    'varName':  'massT',
    'varVal':   [1,10,100]
}

PDCProps= {
    'columns':  ['VEI'],
    'varName':  'prob',
    'varVal':   [.25, .5, .25]
}

# %% [markdown]
# Now we can process the hazard, specfiying the hazard type as the first argument (note that it must be a valid key contained in `E.hazards`):

# %%
E.getHazardExposure('Tephra', TephraProps)
E.getHazardExposure('PDC', PDCProps)
           
# %% [markdown]
# The associated exposure data is now recorded in `E.exposure['Tephra']`. Let's plot the population exposure to 1 kg/m2, comparing VEIs and probabilities of occurrence:

# %%
data = E.exposure['Tephra']
data = data[data.massT == 1]
g = sns.FacetGrid(data, col="VEI")
g.map(sns.barplot, "prob", 'pop_count', order=[10,50,90])

print(data.head(10).to_markdown())

# %% [markdown]
# ## Plotting maps
# We can now combine exposure and hazards on one map

# %%
E.plot(plotExposure='pop', plotBuffer=True, plotHazard='Tephra', hazLevels=[1,10,25,50,100], hazProps={'VEI': 4,'prob': 50})
E.plot(plotExposure='LC', plotBuffer=True, plotHazard='PDC', hazLevels=[.1, .5, .9], hazProps={'VEI': 4})

# %% [markdown]
# ## Save and load projects
# We can now save the result to a file. By default, the project is saved to `E.name` in the root folder.

#%%
E.save(outPth='volcanoes/Merapi/Merapi.volc')

# %% [markdown]
# This project can now be loaded as such:

#%%
E = volcgis.eruption.load('volcanoes/Merapi/Merapi.volc')
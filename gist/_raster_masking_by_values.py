#%%
import os
os.chdir('/Users/seb/Documents/Codes/VolcGIS')
import volcgis 
from rasterio.plot import show
import rasterio as rio
import numpy as np

#%%
haz = 'volcanoes/Gede-Pangrango/_hazard/Tephra/Gede-Pangrango_VEI5_P5.tif'
exp = 'volcanoes/Gede-Pangrango/_data/Landscan.tif'

# with rio.open(haz) as HAZ:
HAZ = rio.open(haz)
EXP = rio.open(exp)

dataHaz = HAZ.read(1)
dataHazM = data>100

dataExp = EXP.read(1)
dataExp[dataExp<1] = 0
np.sum(dataExp[dataHazM])/11.1111111**2

# %%
dataExp2 = dataExp
dataExp2[dataHazM] = 1e6
show(dataExp2)
# %%

test =  rio.open('test_clip.tif')
test2 = test.read(1)
test2[test2<1] = 0
show(test)
np.min(test2)
np.max(test2)
np.sum(test2)
    print(np.sum(test.read(1)))
# %%

lc = rio.open('volcanoes/Cereme/_data/Landcover.tif')
haz = rio.open('volcanoes/Cereme/_hazard/Tephra/Cereme_VEI5_P5.tif')

LC = lc.read(1)
HAZ = haz.read(1)

idx = HAZ>=10

# %%

#%%
import os
import glob
import pandas as pd

volcDB = pd.read_csv('csv/volcCoordinates.csv')

for i in range(0, volcDB.shape[0]):
    name = volcDB.loc[i, 'name']
    
    tephra = os.path.join('volcanoes', name, '_hazard', 'Tephra', '*')
    for f in glob.glob(tephra):
        os.remove(f)
        
    exposure = os.path.join('volcanoes', name, '_exposure', '*')
    for f in glob.glob(exposure):
        os.remove(f)
        
# %%

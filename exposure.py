#%%
import os
# os.chdir('/Users/seb/Documents/Codes/VolcGIS/plottingApp/')

from bokeh.io import curdoc
from bokeh.layouts import column, layout
from bokeh.models import ColumnDataSource, Div, Select, Slider, RadioGroup, TextInput, HoverTool
from bokeh.plotting import figure
from bokeh.models.widgets import DataTable
from bokeh.tile_providers import get_provider

import numpy as np
import pandas as pd
pd.set_option('mode.chained_assignment', None)
# from bokeh.io import show


# Volcano reference
volcDB = pd.read_csv('csv/volcCoordinates.csv')
volcDB = volcDB.rename({'name': 'volcano'}, axis=1)
# volcDB = volcDB.drop(['Unnamed: 8'], axis=1)

def mergeData(data, db):
    db = db.set_index('volcano')[['country']]
    data = data.set_index('volcano').join(db).reset_index()
    data = data.replace('Philippines and SE Asia', 'Philippines')
    data = data.sort_values(by=['country', 'volcano'])
    return data

def wgs84_to_web_mercator(df, lon="lon", lat="lat"):
      k = 6378137
      df["x"] = df[lon] * (k * np.pi/180.0)
      df["y"] = np.log(np.tan((90 + df[lat]) * np.pi/360.0)) * k
      return df

#%%

# Prepare exposure data 
ALL = pd.read_csv('csv/MASTER_exposure_v2.csv')
ALL = ALL.drop(['Unnamed: 0'], axis=1)
ALL = ALL.rename({'pop_count': 'Population', 'area_crops': 'Crops', 'area_urban': 'Urban', 'buildingsLoss': 'Buildings'}, axis=1)
ALL = ALL[ALL.month==0] # No seasonal analysis
ALL = ALL.set_index('volcano').join(volcDB.set_index('volcano')[['lat', 'lon']]).reset_index() # 2021-05-13: add lat and lon
ALL = wgs84_to_web_mercator(ALL)

# Prepare x_range
test = ALL[(ALL.VEI==5) & (ALL.hazard=='Tephra') & (ALL.prob==50) & (ALL.mass==1)]
xRange = mergeData(test, volcDB)

# yaxis reference:
ylabel = {
    'Population': 'Population (n)',
    'Crops': 'Crops area (km2)',
    'Urban': 'Urban area (km2)',
    'RNDS': 'Road network disruption score',
    'Buildings': 'Buildings losses (million USD)'
}

# Color references:
clr = {
    'ind': ['#a1dab4','#41b6c4','#225ea8'], 
    'phi': ['#fbb4b9','#f768a1','#ae017e'],
}

## TEPHRA 
# Create input controls
T_hazard = Select(title='Hazard', options=['Tephra', 'Large Clasts', 'Dome collapse', 'Column collapse'], value='Tephra')
T_exp = Select(title='Exposure', options=['Population', 'Crops', 'Urban', 'RNDS', 'Buildings'], value='Population')
T_VEI = Select(title='VEI', options=['3', '4', '5'], value='4')
T_mass = Select(title='Isomass (kg/m2)', options=['1', '5', '50', '100'], value='1')
T_prob = Select(title='Probability', options=['10', '50', '90'], value='50')
T_bar = Select(title='Variability', options=['VEI', 'Probability', 'Mass'], value='VEI')
T_volume = Select(title='BAF volume (x10^6 m3)', options=['0.45', '9.8'], value='0.45')
T_buffer = Select(title='BAF buffer (m)', options=['300', '990'], value='300')

# Create Column Data Source that will be used by the plot
T_low = ColumnDataSource(data=dict(x=[], y=[], lat=[], lon=[], rIn=[], rOut=[], color=[], rng=[]))
T_mid = ColumnDataSource(data=dict(x=[], y=[], lat=[], lon=[], rIn=[], rOut=[], color=[], rng=[]))
T_top = ColumnDataSource(data=dict(x=[], y=[], lat=[], lon=[], rIn=[], rOut=[], color=[], rng=[]))

# Define behaviour for hazard TEPHRA
def select_Tephra():
    # If VEI
    if T_exp.value not in ['RNDS', 'Buildings']:
        if T_bar.value == 'VEI':
            selected = ALL[(ALL.hazard=='Tephra') & (ALL.prob==int(T_prob.value)) & (ALL.mass==int(T_mass.value))]
        elif T_bar.value == 'Probability':
            selected = ALL[(ALL.hazard=='Tephra') & (ALL.VEI==int(T_VEI.value)) & (ALL.mass==int(T_mass.value))]
        elif T_bar.value == 'Mass':
            selected = ALL[(ALL.hazard=='Tephra') & (ALL.VEI==int(T_VEI.value)) & (ALL.prob==int(T_prob.value))]
    else:
        if T_bar.value == 'VEI':
            selected = ALL[(ALL.hazard=='Tephra') & (ALL.prob==int(T_prob.value))]
        elif T_bar.value == 'Probability':
            selected = ALL[(ALL.hazard=='Tephra') & (ALL.VEI==int(T_VEI.value))]
        elif T_bar.value == 'Mass':
            selected = ALL[(ALL.hazard=='Tephra') & (ALL.VEI==int(T_VEI.value)) & (ALL.prob==int(T_prob.value))]
            
    selected = mergeData(selected, volcDB)
    return selected

# Define behaviour for hazard LARGE CLASTS
def select_LC():
    # If VEI
    if T_bar.value == 'VEI':
        selected = ALL[(ALL.hazard=='LC') & (ALL.prob==int(T_prob.value))]
    elif T_bar.value == 'Probability':
        selected = ALL[(ALL.hazard=='LC') & (ALL.VEI==int(T_VEI.value))]
    selected = mergeData(selected, volcDB)
    return selected

# Define behaviour for hazard BAF
def select_BAF():
    if T_bar.value == 'Probability':
        selected = ALL[(ALL.hazard=='BAF') & (ALL.buffer==int(T_buffer.value)) & (ALL.volume==float(T_volume.value)*1e6)]
    elif T_bar.value == 'Buffer':
        selected = ALL[(ALL.hazard=='BAF') & (ALL.prob==int(T_prob.value)) & (ALL.volume==float(T_volume.value)*1e6)]
    elif T_bar.value == 'Volume':
        selected = ALL[(ALL.hazard=='BAF') & (ALL.buffer==int(T_buffer.value)) & (ALL.prob==int(T_prob.value))]
    
    selected = mergeData(selected, volcDB)
    return selected

# Define behaviour for hazard PDC
def select_PDC():
    if T_bar.value == 'VEI':
        selected = ALL[(ALL.hazard=='PDC') & (ALL.prob==int(T_prob.value))]
    elif T_bar.value == 'Probability':
        selected = ALL[(ALL.hazard=='PDC') & (ALL.VEI==int(T_VEI.value))]

    selected = mergeData(selected, volcDB)
    return selected

def update():
    
    # Retrieve data per hazard
    controls = [T_hazard, T_exp, T_bar, T_VEI, T_prob, T_mass, T_buffer, T_volume]
    for control in controls:
        control.disabled = False
        
    T_buffer.disabled = True
    T_volume.disabled = True
    
    if T_hazard.value == 'Tephra':
        if T_bar.options != ['VEI', 'Probability', 'Mass']:
            T_bar.options = ['VEI', 'Probability', 'Mass']
            T_bar.value = 'Probability'
            
        df = select_Tephra()

    elif T_hazard.value == 'Large Clasts':
        T_mass.disabled = True
        
        if T_bar.options != ['VEI', 'Probability']:
            T_bar.options = ['VEI', 'Probability']
            T_bar.value = 'Probability'
        
        df = select_LC()
        
    elif T_hazard.value == 'Column collapse':
        T_mass.disabled = True
        
        if T_bar.options != ['VEI', 'Probability']:
            T_bar.options = ['VEI', 'Probability']
            T_bar.value = 'Probability'
        
        df = select_PDC()
        
    elif T_hazard.value == 'Dome collapse':
        T_buffer.disabled = False
        T_volume.disabled = False
        T_mass.disabled = True
        T_VEI.disabled = True
        
        if T_bar.options != ['Probability', 'Buffer', 'Volume']:
            T_bar.options = ['Probability', 'Buffer', 'Volume']
            T_bar.value = 'Probability'
            
        df = select_BAF()
        
    if T_bar.value == 'Probability':
        T_prob.disabled = True
    elif T_bar.value == 'VEI':
        T_VEI.disabled = True
    elif T_bar.value == 'Mass':
        T_mass.disabled = True
    elif T_bar.value == 'Volume':
        T_volume.disabled = True
    elif T_bar.value == 'Buffer':
        T_buffer.disabled = True
        
    if T_exp.value in ['RNDS', 'Buildings']:
        T_mass.disabled = True
    else:
        T_mass.disabled = False
        
    # Define source data
    if T_bar.value == 'VEI':
        low = df[df.VEI == 3]
        mid = df[df.VEI == 4]
        top = df[df.VEI == 5]
        ttl = 'VEI'
        rngCol = 'VEI'
    elif T_bar.value == 'Probability':
        low = df[df.prob == 90]
        mid = df[df.prob == 50]
        top = df[df.prob == 10]
        ttl = 'Probability'
        rngCol = 'prob'
    elif T_bar.value == 'Mass':
        low = df[df.mass == 100]
        mid = df[df.mass == 10]
        top = df[df.mass == 1]
        if T_exp.value in ['RNDS', 'Buildings']:
            low = df[df.mass.isna()]
            mid = df[df.mass.isna()]
            top = df[df.mass.isna()]
            
            
        ttl = 'Mass threshold'
        rngCol = 'mass'
    elif T_bar.value == 'Volume':
        low = df[df.volume == 450000]
        mid = df[df.volume == 10]
        top = df[df.volume == 9800000]
        ttl = 'BAF volume'
        rngCol = 'volume'
    elif T_bar.value == 'Buffer':
        low = df[df.buffer == 300]
        mid = df[df.buffer == 10]
        top = df[df.buffer == 990]
        ttl = 'BAF buffer'
        rngCol = 'buffer'
    
    # Update figure
    TOOLTIPS=["(ttl, @y)"]

    low['color'] = np.where(low['country']=='Indonesia', clr['ind'][0], clr['phi'][0])
    mid['color'] = np.where(mid['country']=='Indonesia', clr['ind'][1], clr['phi'][1])
    top['color'] = np.where(top['country']=='Indonesia', clr['ind'][2], clr['phi'][2])

    # Scale radii for map plotting
    sMax = 1e5 # Maximum size to plot
    rMax = top[T_exp.value].max()
    
    rtop = top[T_exp.value] * sMax / rMax
    rmid = mid[T_exp.value] * sMax / rMax
    rlow = low[T_exp.value] * sMax / rMax
    
    T_low.data = dict(x=low['volcano'], y=low[T_exp.value], lat=low['y'], lon=low['x'], rIn=[0]*low.shape[0], rOut=rlow, color=low['color'], rng=low[rngCol])
    T_mid.data = dict(x=mid['volcano'], y=mid[T_exp.value], lat=mid['y'], lon=mid['x'], rIn=rlow, rOut=rmid, color=mid['color'], rng=mid[rngCol])
    T_top.data = dict(x=top['volcano'], y=top[T_exp.value], lat=top['y'], lon=top['x'], rIn=rmid, rOut=rtop, color=top['color'], rng=top[rngCol])
    
    p.yaxis.axis_label = ylabel[T_exp.value]
    p.title.text = ttl
    
    hover = HoverTool(tooltips=[(ttl, "@rng")])
    p.add_tools(hover)
    
    hoverMap = HoverTool(tooltips=[('Volcano', "@x")])
    map.add_tools(hoverMap)
    
 
# Bar plot    
p = figure(x_range=xRange['volcano'], plot_height=500,plot_width=1000, toolbar_location=None, sizing_mode="scale_both")
p.vbar(x='x', top='y', color='color', width=0.9, source=T_top)    
p.vbar(x='x', top='y', color='color', width=0.9, source=T_mid)
p.vbar(x='x', top='y', color='color', width=0.9, source=T_low)

p.xgrid.grid_line_color = None
p.y_range.start = 0
p.xaxis.major_label_orientation = 1
p.x_range.range_padding = 0.01

# Map
tile_provider = get_provider('CARTODBPOSITRON_RETINA')
# range bounds supplied in web mercator coordinates
map = figure(x_axis_type="mercator", y_axis_type="mercator", plot_height=500 ,plot_width=1000)
map.add_tile(tile_provider)
map.annulus(x='lon', y='lat', inner_radius='rIn', outer_radius='rOut', color='color', source=T_top, alpha=.5)
map.annulus(x='lon', y='lat', inner_radius='rIn', outer_radius='rOut', color='color', source=T_mid, alpha=.5)
map.annulus(x='lon', y='lat', inner_radius='rIn', outer_radius='rOut', color='color', source=T_low, alpha=.5)
map.triangle(x='lon', y='lat', size=5, fill_color='red',line_color='black', line_width=.5, source=T_top)

# Controls
controls = [T_hazard, T_exp, T_bar, T_VEI, T_prob, T_mass, T_buffer, T_volume]
for control in controls:
    control.on_change('value', lambda attr, old, new: update())

inputs = column(*controls, width=220, height=1000)

# Layout
inputs.sizing_mode = "fixed"
l = layout([
    [inputs, [p, map]],
], sizing_mode="scale_both")

# print(T_hazard.value)
# print(T_exp.value)
# print(T_VEI.value)
# print(T_mass.value)
# print(T_prob.value)
# print(T_bar.value)

update()  # initial load of the data

curdoc().add_root(l)
curdoc().title = "Exposure"

# %%

#%%
from bokeh.io import output_file, show, output_notebook
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, FactorRange, Panel, Tabs
from bokeh.themes import built_in_themes
from bokeh.io import curdoc
from bokeh.layouts import column

output_file("exposure.html")
curdoc().theme = 'dark_minimal'

import pandas as pd

# Volcano reference
volcDB = pd.read_csv('volcCoordinates.csv')
volcDB = volcDB.rename({'name': 'volcano'}, axis=1)
volcDB = volcDB.drop(['Unnamed: 8'], axis=1)

# Prepare exposure data 
ALL = pd.read_csv('results.csv')
ALL = ALL.drop(['Unnamed: 0'], axis=1)


#%%
def mergeData(data, db):
    db = db.set_index('volcano')[['country']]
    data = data.set_index('volcano').join(db).reset_index()
    data = data.replace('Philippines and SE Asia', 'Philippines')
    data = data.sort_values(by=['country', 'volcano'])
    return data






#%% Large clasts
# test = ALL[(ALL.VEI==5) & (ALL.hazard=='Tephra') & (ALL.prob==50) & (ALL.mass==1)]
# prb50 = ALL[(ALL.VEI==5) & (ALL.hazard=='LC') & (ALL.prob==50)]
# prb10 = ALL[(ALL.VEI==5) & (ALL.hazard=='LC') & (ALL.prob==10)]
# prb90 = ALL[(ALL.VEI==5) & (ALL.hazard=='LC') & (ALL.prob==90)]
prb = ALL[(ALL.VEI==5) & (ALL.hazard=='LC')]

colors = ["#c9d9d3", "#718dbf", "#e84d60"]

target = 'pop_count'



prb50 = mergeData(prb50, volcDB)
prb10 = mergeData(prb10, volcDB)
prb10 = mergeData(prb10, volcDB)

x = list(zip(prb.country, prb.volcano))

#%%
p = figure(x_range=FactorRange(*x), plot_height=350,plot_width=1000, title="Tephra VEI5",
           toolbar_location=None)
# p = figure(x_range=test.volcano, plot_height=350,plot_width=900, title="Tephra VEI5",
#            toolbar_location=None)

# p.vbar(x=test.volcano, top=test.pop_count, width=0.9)
# p.vbar(x=prb10.volcano, top=prb10.pop_count, width=0.9, color=colors[2])
# p.vbar(x=prb50.volcano, top=prb50.pop_count, width=0.9, color=colors[1])
p.vbar(x=x, top=prb10.pop_count, width=0.9, color=colors[2])
p.vbar(x=x, top=prb50.pop_count, width=0.9, color=colors[0])
p.vbar(x=x, top=prb90.pop_count, width=0.9, color=colors[1])

p.xgrid.grid_line_color = None
p.y_range.start = 0
p.xaxis.major_label_orientation = 1
p.x_range.range_padding = 0.01
p.group_padding = 1
show(p)
# %%

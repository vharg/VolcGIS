# Tips

## Plotting geopandas basemaps with contextily

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

## Join the RSDS data to the original .feather file

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

## Starting the plotting app on Wovodat

Start the bokeh plotting app in `exposure.py`, making it available at [gee.wovodat.org/exposure](gee.wovodat.org/exposure)

```bash
conda activate ee
nohup bokeh serve --allow-websocket-origin='*' exposure.py --log-level=debug
```

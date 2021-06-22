# Change log

## 2021-06-22

- Replaced `gpd.clip()` by `gpd.sjoin`
- 
## 2021-06-21

- Added the example for `PDC` in `MASTER_example.py` and fixed bugs
- Added the option for roads in `getExposure`
- Added the exposure analysis to roads in `getBufferExposure`

## 2021-06-15

- Fixed a problem in `volcgis.exposureAnalysis.getRNDS()` that returned wrong results without raising an error. This was due to polygons with complex shapes (e.g. holes in the polygon), which caused problem when clipping with `gpd`. The more I dug, the more I found that small hazard footprints (e.g. BAF) or those with irregular outlines (e.g. PDC) could produce "noisy" polygons at the periphery when contoured with `rio`. I added some filters, which are hopefully flexible enough to handle most case studies.

## 2021-06-10

- Added `v0.2.0`
- Documented module is now `volcgis.eruption.py`
- Added documentation with `MASTER_example.py` notebook
- `volcgis.exposureAnalysis` is still here for the paper, but it will be replaced by the new functions
- `processHazard_v2.py` is an update for compatibility with `v0.2.0`

## 2021-03-04

- Moved all exposure analysis functions to `volcgis.exposureAnalysis`

### Road Exposure
Updated road exposure from Josh's commit. Namely:
- Removed `ROAD_EXPOSURE` variable and added content to `EXPOSURE`
- Remove `updateRoads` and added to `updateExposure`

### `getRNDS`
- `getRNDS` now returns three arguments, including `roadLength` and `RSDS`
- `roadLength` is added to `EXPOSURE`
- `RSDS` is a pd.DataFrame that contains `RSDS` value for each road segment defined by `Road_ID`. Each column is a different hazard occurrence (i.e. hazard type, VEI, prob etc)

### `Gekko`
- Prepared `processHazard.py` for parallel processing on Gekko using Job Arrays
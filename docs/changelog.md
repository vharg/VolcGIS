# Change log

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
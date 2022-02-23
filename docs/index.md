# VolcGIS

Welcome to the documentation page for **VolcGIS**!

**VolcGIS** is a Python package that implements several functions for exposure analyses to volcanic hazards. It is developed by the [Volcanic Hazards and Risk Group](https://www.earthobservatory.sg/research-group/volcanic-hazards-and-risk-susanna-jenkins) group at the [Earth Observatory of Singapore](https://www.earthobservatory.sg), Nanyang Technological University. In a nutshell, **VolcGIS** provides high-level functions to streamline the process of:

- Setting up a GIS to fully exploit the spatial relationship between hazard and exposure datasets
- Pre-process hazard and exposure datasets
- Extract global exposure data and figures for various hazard types and scenarios
- Implement methodologies to estimate *some* aspects of vulnerability, impact and risk

**VolcGIS** is available on [Github](https://github.com/vharg/VolcGIS), documentation is available at [vharg.github.io/VolcGIS/](https://vharg.github.io/VolcGIS/). Check out the [demo notebook](https://vharg.github.io/VolcGIS/MASTER_example/) to get started.

## Install

**VolcGIS** was written in Python 3.9 and was tested using `Rasterio 1.1.8`, `rioxarray 0.4`, `GDAL 3.1.4` and `Numpy 1.19.4` for raster analyses and `GeoPandas 0.8.1` and `Shapely 1.7.1` for vector analyses.

### Setup conda environment

Create environment and set the channel to `conda-forge`:

```
conda config --env --add channels conda-forge
conda config --env --set channel_priority strict
```

Then:

```
conda install -c conda-forge pyarrow rioxarray rasterio geopandas bokeh openpyxl contextily osmnx pyrosm tabulate
```

### Additional pip packages

Install these packages with `pip`:

```
pip install utm
pip install alive-progress
```

### Documentation

To install the documentation:

```
pip install mkdocs-material mkdocstrings livereload mkdocs-jupyter
```

## Notes on the use of VolcGIS

Anyone who has ever had a critical glance at exposure data derived from global datasets in volcanology or any other natural hazards will have noticed *how variable* exposure figures from different studies can be. This variability can be explained by various aspects. The first one is *the variability of global datasets*. For population, different datasets (e.g. GHSL, Landscan) rely on different assumptions (see [Freire et al. 2019](https://doi.org/10.3390/ijgi8080341) for more detail). Similarly global land-cover maps, use different resolutions, classification algorithms and training/validation datasets, which induces a discrepancy in the results. The second one is *the manipulation of geospatial data* which, through such tasks as re-projection between different coordinate systems and interpolation, add a layer of *noise* to the data. Whereas it is important to keep these operations under as much control as possible, the same population count achieved using two different Python libraries shows a 3-5% variability in the result.

As a caveat, a critical interpretation of the results is required. Although not systematically quantified, a ±10% variability of the results is a sensible baseline.

> Freire, S., Florczyk, A., Pesaresi, M., Sliuzas, R., 2019. An Improved Global Analysis of Population Distribution in Proximity to Active Volcanoes, 1975–2015. ISPRS International Journal of Geo-Information 8, 341. https://doi.org/10.3390/ijgi8080341

## License

**VolcGIS** is published under a GNU GPL v3 License

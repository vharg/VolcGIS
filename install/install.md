# Install

## Setup conda environment

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

## Additional pip packages

Install these packages with `pip`:

```
pip install utm
pip install alive-progress
```

## Documentation 

To install the documentation:

```
pip install mkdocs-material mkdocstrings livereload mkdocs-jupyter
```
# Raster-Basics
Basic geoprocessing tools to handle geotiff and shapefile data

## Installation
Activate the desired python or conda environment, then enter `pip install raster_basics`

To import the modules, type `import raster_basics` or import specific functions with `from raster_basics.RasterBasics import function`

All of the functions can be imported with:

 - `from raster_basics.RasterBasics import show_fig, rasterLike`
 - `from raster_basics.RasterBasics import shpReprojection, shpClip, tifReprojectionResample,fillHole, mosaic_files`
 - `from raster_basics.RasterBasics import extract_along_line,points_along_lines, end_points`

## Functions
 - `show_fig`: show geotiff as a matplotlib figure
 - `rasterLike`: save array as a raster
 - `shpReprojection`: reproject shapefile
 - `shpClip`: clip geotiff by a shapefile
 - `tifReprojectionResample`: reproject, resample, and/or clip geotiff
 - `fillHole`: fill missing data in geotiff
 - `mosaic_files`: mosaic two or more geotiffs together
 - `extract_along_line`: sample for points on a line (this uses geocoordinate points as line)
 - `points_along_lines`: sample for points on a line (this uses a shapefile line)
 - `end_points`: obtain raster values at point locations

## Tutorial
Follow the Jupyter Notebook `rasterio_basics.ipynb` to see how functions can be used. The notebook `rasterio_basics-package.ipynb` is identical to `rasterio_basics.ipynb` but uses the installed Python Raster_Basics module. Access to the [sample data used in the tutorial can be found here](https://drive.google.com/file/d/1lNiQBo-rNe2_VC6vUCM2gfp-Z-F2Q49c/view?usp=share_link). Run the notebook code in the same folder as the data.

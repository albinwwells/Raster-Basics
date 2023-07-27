# Raster-Basics
Basic geoprocessing tools to handle geotiff and shapefile data

Contains some additional tools for basic glacier calculations and plotting

## Installation
Activate the desired python or conda environment, then enter `pip install raster_basics`

To import the modules, type `import raster_basics` or import specific functions with `from raster_basics.RasterBasics import function`

If you've already downloaded the module, enter `sudo pip install raster_basics --upgrade` while in the conda environment to upgrade to the newest version

All of the functions can be imported with:

 - `from raster_basics.RasterBasics import show_fig, rasterLike`
 - `from raster_basics.RasterBasics import shpReprojection, shpClip, tifReprojectionResample,fillHole, mosaic_files`
 - `from raster_basics.RasterBasics import extract_along_line,points_along_lines, end_points`
 
### Installing Dependencies
Dependencies to run `raster_basics` can be installed manually or via the `requirements.txt` file. Download the file, navigate to it's directory in the terminal, and install dependencies by typing `pip install -r requirements.txt` or `conda install --file requirements.txt`. If creating a new conda environment for `raster_basics`, you may also need to `conda install pip` initially. This should install all the dependencies needed to run the `raster_basics` module.

## Functions
 - `rOpen`: open a geotiff file as an array, and/or get the geotiff resolution and coordinate system
 - `show_fig`: show geotiff as a matplotlib figure
 - `rasterLike`: save array as a raster
 - `shpReprojection`: reproject shapefile
 - `shpBuffer`: add a buffer to shapefile
 - `shpClip`: clip geotiff by a shapefile
 - `tifReprojectionResample`: reproject, resample, and/or clip geotiff
 - `rasterMath`: perform basic arithmetic on two raster files (even if they have difference crs, res, and/or extents)
 - `fillHole`: fill missing data in geotiff
 - `fillArrayHole`: helper function for `fillHole`, but works directly with arrays
 - `mosaic_files`: mosaic two or more geotiffs together
 - `extract_along_line`: sample for points on a line (this uses geocoordinate points as line)
 - `points_along_lines`: sample for points on a line (this uses a shapefile line)
 - `end_points`: obtain raster values at point locations
 - `distance_to_shp`: gets the distance from every pixel to a shapefile

## Tutorial
Follow the Jupyter Notebook `rasterio_basics.ipynb` to see how functions can be used. The notebook `rasterio_basics-package.ipynb` is identical to `rasterio_basics.ipynb` but uses the installed Python Raster_Basics module. Access to the [sample data used in the tutorial can be found here](https://drive.google.com/file/d/1lNiQBo-rNe2_VC6vUCM2gfp-Z-F2Q49c/view?usp=share_link). Run the notebook code in the same folder as the data.

`wolverine_index_site_method.ipynb` is a notebook that shows an example comparing glaciologic and geodetic mass balance calculations using the `raster-basics` module. This includes an uncertainty estimate of the geodetic mass balance following a Monte Carlo approach.


# Other functions in the package
There are a few more functions to help produce specific plots and glacier calculations. These can be imported with:

 - `from raster_basics.GlacierFunctions import FUNCTION_NAME`
 - `from raster_basics.BaseFunctions import FUNCTION_NAME`
 - `from raster_basics.DataPlots import FUNCTION_NAME`
 - `from raster_basics.SmoothingFunctions import FUNCTION_NAME`
 - `from raster_basics.NoiseFunctions import FUNCTION_NAME`
 
 where `FUNCTION_NAME` is the name of the desired function.
 
## Example Glacier Functions
  - `glacierArea`: get the area of a glacier
  - `totalMassBalance`: get the total mass change on a glacier (pixel by pixel)
  - `totalMassBalanceValue`: get the total mass change on a glacier (area-sum of pixels)
  - `divQ`: get the flux divergence
  - `glacierSlope`: alternate function to obtain the riserun slope from an array
  - `demHillshade`: create a hillshade from a DEM
  - `velocityAspect`: returns the aspect of velocity based on vx and vy
  - `velAspectAngle`: get the angle between two arrays of velocity aspect
  - `velFlowlineAspect`: return the aspect of every pixel based on a shapely geometry flowline
  -  `velFlowlineOutlineAspect`: return the aspect of every pixel based on a shapely geometry flowline and outline
  - `distance_from_line`: get an array representing the distance of every pixel from any True value in an array
 
## Example Base Functions
  - `glacierOutline`: uses a shapefile and raster array of ones to create a binary array of glacier terrain (1 is glacier, 0 is off-glacier)
  - `altitudeAggregation`: returns statistics for a desired statistic in an elevation bin, bin boundaries, count, number, standar deviation, min, and max
  - `binPercentile`: returns binned statistics value at a desired percentile of a bin (e.g. 25th percentile MB value per elevation bin)
  - `latlonTiffIndex`: obtain raster array index values from a lat/lon coordinate pair

## Example Plotting Functions
 Note that all of these plotting functions only save the image if `savefig=True` (by default, it is set to `False`).
  - `plotData`: another method to plot an array. Can handle quiver inputs for arrows on velocity plots (quiver input from `velPlot` function)
  - `plotData3`: same as `plotData` but places 3 plots side-by-side
  - `plotMany`: an arbitrary number of plots side-by-side in a grid
  - `elevationBinPlot`: plots two sets of data: line plot for MB and horizontal bar graph for elevation bin size
  - `elevationBinPlot3Subfigs`: plots `elevationBinPlot` along with 3 more subplots from desired input arrays
  - `plotDataPoints`: plots labeled data points on top of an array basemap. This is used to plot stake locations, for example
  - `velPlot`: this does no actually plot anything. Rather, it returns the quiver input needed to plot velocity arrows in other functions such as `plotData`
  - `plotClassify`: plot array with discrete values
  - `plotContinuous`: plot array with continuous values
  - `plot_binned_data`: returns a data array with values corresponding to alitudinally-aggregated elevation bins
  
## Example Smoothing Functions
  - `sgolay2d`: Savitzky-Golay smoothing filter to eliminate high frequency noise via moving average; a low-pass filter
  - `gaussianFilter`: Gaussian (normal) smoothing filter
  - `dynamicSmoothing`: Gaussian filter, dynamic smoothing with changing pixel window size
  - `dynamicSmoothingExponential`: Exponential filter, dynamic smoothing with changing pixel window size
  - `add_array_border`: Helper function to add a border around an array in `dynamicSmoothing` and `dynamicSmoothingExponential`
  - `smoothingCorrection`: Apply a correction to smoothed data products, such that smoothing does not decrease the mean values Correction is based on glacier centerline values
  - `distance_scaling_correction`: corrects a smoothed data products based on a scaling from relative distance from centerline

## Example Noise Functions
  - `add_grf_noise`: Introduce spatially-correlated noise via Gaussian random fields
  - `add_ice_thickness_bias`: Add ice thickness bias based on parabolic vs V-shape bed cross-section assumption



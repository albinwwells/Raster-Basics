import rasterio
import rasterio.plot
from rasterio.plot import show
import rasterio.mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.merge import merge
from rasterio.windows import from_bounds
from rasterio.fill import fillnodata

import rioxarray
import xarray

import glob, os
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib import cm
import fiona
import shapely.geometry
import geopandas as gpd
import pandas as pd



""" Simple way to plot raster file """

def show_fig(image, title=None, color='Spectral', ctitle='colorbar title', bounds=None, res=None, vmin=None, vmax=None, savefig=False):
    """
    Plotting a raster file
	image: input array (array-like) (e.g. rasterio.open('raster.tif').read(1))
	title: figure title (str) (e.g. 'velocity plot')
	color: figure colormap (str) (e.g. 'Spectral', 'RdBu', etc)
	ctitle: colorbar title (str) (e.g. 'm/yr')
	bounds: raster plot bounds, input: (left, right, bottom, top) (floats) 
	res: raster pixel resolution (int or float)
	vmin, vmax: plot minimum and maximum (int or float)
	savefig: If true, figure is saved as a .jpg in the current directory with the filename being the title input (boolean)
    """
    fig, ax = plt.subplots(figsize=(12,6))
    c = ax.imshow(image, cmap=color, extent=bounds, vmin=vmin, vmax=vmax)
    if res != None:
        ax.add_artist(ScaleBar(dx=res, units='m')) # add scalebar
    fig.colorbar(c, label=ctitle)
    fig.suptitle(title)
    plt.show()
    if savefig == True:
        fig.savefig(title + '.jpg', dpi=1000) # to save the plot as a jpg image


""" Geoprocessing tools:

	Reproject shapefile
	Clip shapefile
	Reproject and clip geotiff
	Fill missing geotiff data
	Mosaic geotiff
	Save array as raster with identical properties
		
"""

def shpReprojection(shapefile, crs, dst='reprojected_shape.shp'):
    """
    Reproject a shapefile
	shapefile: input shapefile path (str) (e.g. 'input_shapefile.shp')
	crs: output coordinate system (str) (e.g. 'EPSG:4326')
	dst: output filename (str) (e.g. 'output_shapefile.shp')
    """
    src = gpd.read_file(shapefile)
    src = src.to_crs(crs)
    src.to_file(dst)
    

def shpClip(geotiff, shapefile, destination, nan_val=0, fill=True, pad_size=0):
    """
    Clip a geotiff with a shapefile
	geotiff: input raster to clip (str) (e.g. 'input.tif')
	shapefile: shapefile used for clipping (str) (e.g. 'input_shape.shp')
	destination: output clipped raster filename (str) (e.g. 'output.tif')
	nan_val: NaN value (int or float)
	fill: If True, the pixels outside the features will be set to nan_val. If False, the output array will contain the original pixel data, and only the mask will be based on shapes (boolean) 
	pad_size: width of padding in pixels (int or float)    
    """
    
    with rasterio.open(geotiff) as src:
        # we need to make a temporary shapefile with the same crs as the geotiff 
        shpReprojection(shapefile, src.crs, dst='temp.shp')
        with fiona.open('temp.shp', 'r') as openshape:
            shapes = [feature['geometry'] for feature in openshape]

        out_image, out_transform = rasterio.mask.mask(
            src, 
            shapes, 
            crop=True, 
            filled=fill, 
            nodata=nan_val, 
            pad=True, 
            pad_width=pad_size
        )
        
        out_image[np.isnan(out_image)] = 0
        kwargs = src.meta
        kwargs.update({'driver': 'GTiff',
                       'height': out_image.shape[1],
                       'width': out_image.shape[2],
                       'transform': out_transform})
        with rasterio.open(destination, 'w', **kwargs) as dst:
            dst.write(out_image)
            
    for f in glob.glob('temp.*'): # remove the 'temp.*' shapefiles
        os.remove(f)


def tifReprojectionResample(file, reprojected_tif, crs, res, interp, extent_file=None):
    """
    Reproject, resample, and/or clip raster extent
	file: input raster to resample, reproject, and/or clip (str) (e.g. 'input.tif')
	reprojected_tif: output filename (str) (e.g. 'output_raster.tif')
	crs: target raster coordinate system (str) (e.g. 'EPGS:4326')
	res: target raster resolution (int or float)
	interp: resampling interpolation method (Resampling method from rasterio.warp module)  (e.g. rasterio.warp.Resampling.cubic_spline)
	extent_file: raster file to match extent (str)    
    """
    with rasterio.open(file) as src:
        if extent_file is None:
            # keep initial file bounds
            transform, width, height = calculate_default_transform(
                src.crs,
                crs,
                src.width,
                src.height,
                *src.bounds,
                resolution=res
            )
        else:
            # this readjusts the bounds to be that of extent_file
            dst = rasterio.open(extent_file)
            transform, width, height = calculate_default_transform(
                src.crs,
                crs,
                dst.width,
                dst.height,
                *dst.bounds,
                resolution=res
            )

        kwargs = src.meta.copy()
        kwargs.update({
            'crs': crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rasterio.open(reprojected_tif, 'w', **kwargs) as dst:
            reproject(
                source=rasterio.band(src, 1),
                destination=rasterio.band(dst, 1),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=crs,
                resampling=interp   # Resampling.cubic_spline
            )


def fillHole(file, dest='output_filled.tif', dist=10, iters=1):
    """
    Fill hole in raster
	file: input raster to fill (str) (e.g. 'input.tif')
	dest: output filename (str) (e.g. 'output_raster.tif')
	dist: the maximum number of pixels to search in all directions to find values to interpolate from (int)
	iters: the number of 3x3 smoothing filter passes to run (int)
    """
    with rasterio.open(file) as src:
        profile = src.profile
        inputs = src.read(1)
        
        fillmask = inputs.copy() # fillnodata is applied where the mask=0
        fillmask[inputs>=0] = 1
        fillmask[fillmask!=1] = 0
        
        inputFilled = fillnodata(inputs, mask=fillmask, max_search_distance=dist, smoothing_iterations=iters)
        inputFilled[pd.isnull(inputFilled) == True] = 0

    with rasterio.open(dest, 'w', **profile) as dst:
        dst.write_band(1, inputFilled)


def mosaic_files(files, mosaic_output):
    """
   Create a mosiac of the obtained raster files. reproject files to the same crs if needed
	files: input rasters to mosaic (list of str) (e.g. ['input1.tif', 'input2.tif', 'input3.tif', 'input4.tif'])
	mosaic_output: the number of 3x3 smoothing filter passes to run (int)
    """
    src_files_to_mosaic = []
    dst_crs = rasterio.open(files[0]).crs
    for file in files:
        src = rasterio.open(file)   
        if src.crs != dst_crs: # NOTE: files must have the same crs
            raise 'Both datasets do not share a common coordinate system. One must be reprojected.'
        src_files_to_mosaic.append(src)
        
    mosaic, out_trans = merge(src_files_to_mosaic)
    out_meta = src.meta.copy()
    out_meta.update({'driver': 'GTiff',
                             'height': mosaic.shape[1],
                             'width': mosaic.shape[2],
                             'transform': out_trans})
    with rasterio.open(mosaic_output, 'w', **out_meta) as dest:
        dest.write(mosaic)


def rasterLike(array, destination, geotiff):
    """
    Save array as raster, assuming it is 'like' the raster (same res, crs, size, etc)
	array: input array to save as a raster (array-like)
	destination: output filename (str) (e.g. 'output_raster.tif')
	geotiff: raster (str) (e.g. 'raster.tif)
    """
    with rasterio.open(geotiff) as src:
        kwargs = src.meta.copy()
    with rasterio.open(destination, 'w', **kwargs) as dst:
        dst.write(array, 1)


""" Extra functions for sampling along lines """


def extract_along_line(xarr, line, n_samples=512):
    # samples line with n_samples number of points
    profile = []
    dist = []
    for i in range(n_samples):
        point = line.interpolate(i / n_samples - 1., normalized=True) # get next point on the line
        value = xarr.sel(x=point.x, y=point.y, method="nearest").data # access the nearest pixel in the xarray
        profile.append(value)
        dist.append([point.x, point.y])
    return profile, dist


def points_along_lines(geotiff, line, ID='SEGMENT_ID'):
    # extract raster values along the line
    xarr = rioxarray.open_rasterio(geotiff).squeeze()
    profile, coordinate = [], []
    ids = line[ID]
    for g in line['geometry']:
        prof, coord = extract_along_line(xarr, g)
        profile.append(prof)
        coordinate.append(coord)

    # get the raster values along the line and the distances along the transect
    vals, line_dists, dist = [], [], []
    for n in range(len(profile)):
        v = [float(x) for x in profile[n]]
        line_d = [[coordinate[n][i], coordinate[n][i+1]] for i in range(len(coordinate[n])-1)]
        d = [math.dist(pair[0], pair[1]) for pair in line_d]
        d = [sum(d[:i]) for i in range(len(d)+1)]
        vals.append(v)
        line_dists.append(line_d)
        dist.append(d)
    return vals, line_dists, dist, ids.to_numpy()


def end_points(xarr, points):
    # obtains raster values at our point selections themselves (not interpolated lines)
    point_spot = []
    for p in points:
        point = xarr.sel(x=p[0], y=p[1], method="nearest").data # access the nearest pixel in the xarray
        point_spot.append(point)
    return point_spot




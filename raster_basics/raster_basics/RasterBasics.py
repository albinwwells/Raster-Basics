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

def show_fig(image, title=None, color='spectral', ctitle='colorbar title', bounds=None, res=None, vmin=None, vmax=None, savefig=False):
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
    src = gpd.read_file(shapefile)
    src = src.to_crs(crs)
    src.to_file(dst)
    

def shpClip(geotiff, shapefile, destination, nan_val=0, fill=True, pad_size=0):
    # clip a geotiff with a shapefile
    with rasterio.open(geotiff) as src:
        # we need to make a temporary shapefile with the same crs as the geotiff 
        shpReprojection(shapefile, src.crs, dst='temp.shp')
        with fiona.open('temp.shp', 'r') as shapefile:
            shapes = [feature['geometry'] for feature in shapefile]

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


def fillHole(file, dist=10, iters=1):
    with rasterio.open(file) as src:
        profile = src.profile
        inputs = src.read(1)
        
        fillmask = inputs.copy() # fillnodata is applied where the mask=0
        fillmask[inputs>=0] = 1
        fillmask[fillmask!=1] = 0
        
        inputFilled = fillnodata(inputs, mask=fillmask, max_search_distance=dist, smoothing_iterations=iters)
        inputFilled[pd.isnull(inputFilled) == True] = 0

    destination = 'gulkana_filled.tif'
    with rasterio.open(destination, 'w', **profile) as dst:
        dst.write_band(1, inputFilled)


def mosaic_files(files, mosaic_output):
    # create a mosiac of the obtained raster files. reproject files to the same crs if needed
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




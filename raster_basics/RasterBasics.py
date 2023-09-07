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

import glob, os, warnings
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib import cm
import fiona
import shapely.geometry
import geopandas as gpd
import pandas as pd


""" Simple open a raster file """

def rOpen(geotiff, band=1, returnArray=True, returnRes=False, returnCrs=False):
    """
    Open a raster file as an array. Options to output raster resolution and coordinate system
        geotiff: input raster filename (e.g. 'raster.tif')
        band: the band to read. Defaults to 1 (int)
        returnArray: return the raster as a numpy array. Default is True (bool)
        returnRes: return the raster resolution; this returns a tuple. Default is False (bool)
        returnCrs: return the raster coordinate system. Default is False (bool)
    """
    vals = []
    if returnArray == True:
        array = rasterio.open(geotiff).read(band)
        vals.append(array)
    if returnRes == True:
        res = rasterio.open(geotiff).res
        vals.append(res)
    if returnCrs == True:
        crs = rasterio.open(geotiff).crs
        vals.append(crs)
    if len(vals) == 1:
        vals = vals[0] 
    return vals

""" Simple way to plot raster file """

def show_fig(image, title=None, color='Spectral', ctitle='', bounds=None, res=None, vmin=None, vcenter=None, vmax=None, savefig=False):
    """
    Plotting a raster file
	image: input array (array-like) (e.g. rasterio.open('raster.tif').read(1))
	title: figure title (str) (e.g. 'velocity plot')
	color: figure colormap (str) (e.g. 'Spectral', 'RdBu', etc)
	ctitle: colorbar title (str) (e.g. 'm/yr')
	bounds: raster plot bounds, input: (left, right, bottom, top) (floats) 
	res: raster pixel resolution (int or float)
	vmin, vmax: plot minimum and maximum (int or float)
	vcenter: center of colormap
	savefig: If true, figure is saved as a .jpg in the current directory with the filename being the title input (boolean)
    """
    fig, ax = plt.subplots(figsize=(12,6))
    if vcenter != None:
    	divnorm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    	c = ax.imshow(image, cmap=color, extent=bounds, norm=divnorm)
    else:
    	c = ax.imshow(image, cmap=color, extent=bounds, vmin=vmin, vmax=vmax)
    if res != None:
        ax.add_artist(ScaleBar(dx=res, units='m')) # add scalebar
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    fig.colorbar(c, label=ctitle)
    fig.suptitle(title)
    plt.show()
    if savefig == True:
    	fig.savefig(title + '.jpg', dpi=500) # to save the plot as a jpg image


def show_fig_subplot(images, titles=None, colors=None, ctitles=None, bounds=None, res=None, vmins=None, vcenters=None, vmaxs=None,
                     suptitle=None, ncols=2, savefig=False):
    """
    Plotting multiple raster files
	images: list of input arrays (list of array-like) (e.g. [rasterio.open('raster.tif').read(1), rasterio.open('raster2.tif').read(1)])
	titles: list of figure titles (list of str) (e.g. ['velocity plot', 'thickness plot'])
	colors: list of figure colormaps (list of str) (e.g. ['Spectral', 'RdBu'])
	ctitles: list of colorbar titles (list of str) (e.g. ['m/yr', 'm'])
	bounds: list of raster plot bounds, input: (list of [left, right, bottom, top]) (list of floats) 
 	res: list of raster pixel resolutions (list of int or float)
	vmins, vcenters, vmaxs: list of plot minimums, centers, and maximums (list of int or float)
	suptitle: title of the plot (str)
	ncols: number of columns for subplots (int)
	savefig: If true, figure is saved as a .jpg in the current directory with the filename being the suptitle input (boolean)
    """
    num_images = len(images)
    
    # Calculate the number of rows and columns for the subplot grid
    num_rows = 1 if num_images <= ncols else (num_images + ncols - 1) // ncols
    num_cols = min(num_images, ncols)
    
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 6 * num_rows))
    if num_images == 1:
        axes = [[axes]]  # Handle the case with only one image
    
    for i in range(num_rows):
        for j in range(num_cols):
            if num_rows > 1:
                ax = axes[i][j]
            else:
                ax = axes[j]  # For a single row, access axes directly

            index = i * num_cols + j
            if index >= num_images:
                ax.axis('off')  # Hide unused subplots if necessary
                continue
            
            image = images[index]
            title = titles[index] if titles is not None else None
            color = colors[index] if colors is not None else 'Spectral'
            ctitle = ctitles[index] if ctitles is not None else ''
            vmin = vmins[index] if vmins is not None else None
            vcenter = vcenters[index] if vcenters is not None else None
            vmax = vmaxs[index] if vmaxs is not None else None
            bound = bounds[index] if bounds is not None else None
            
            if vcenter != None:
                divnorm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
                c = ax.imshow(image, cmap=color, extent=bound, norm=divnorm)
            else:
                c = ax.imshow(image, cmap=color, extent=bound, vmin=vmin, vmax=vmax)
            
            if res is not None:
                if res[index] is not None:
                    ax.add_artist(ScaleBar(dx=res[index], units='m'))
                    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

            ax.set_title(title)
            fig.colorbar(c, ax=ax, label=ctitle)
    
    fig.suptitle(suptitle)
    plt.tight_layout()
    plt.show()
    if savefig == True:
        fig.savefig(suptitle + '.jpg', dpi=1000) # to save the plot as a jpg image


""" Geoprocessing tools:

	Reproject shapefile
 	Buffer shapefile
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
    

def shpBuffer(shapefile, dist=100, dst='buffered_shape.shp'):
    """
    Add a buffer around a shapefile. Note that the shapefile must be in a projected coordinate system (e.g. UTM)
        shapefile: input shapefile path (str) (e.g. 'input_shapefile.shp')
        dist: buffer distance, m (int, float) (500)
        dst: output filename (str) (e.g. 'output_shapefile.shp')
    """
    shape = gpd.read_file(shapefile)
    
    if shape.crs.is_geographic == True: 
        warnings.warn("The shapefile is in a geographic coordinate system. Reproject the shapefile and try again")
        return None

    # Create a buffer
    buffered_data = shape.buffer(dist)

    # Create a new GeoDataFrame with the buffered data
    buffered_gdf = gpd.GeoDataFrame(geometry=buffered_data, crs=shape.crs)
    buffered_gdf.to_file(dst) # save the buffered shapefile


def shpClip(geotiff, shapefile, destination, nan_val=0, fill=True, pad_size=0, crop=True):
    """
    Clip a geotiff with a shapefile
	geotiff: input raster to clip (str) (e.g. 'input.tif')
	shapefile: shapefile used for clipping (str) (e.g. 'input_shape.shp')
	destination: output clipped raster filename (str) (e.g. 'output.tif')
	nan_val: NaN value (int or float)
	fill: If True, the pixels outside the features will be set to nan_val. If False, the output array will contain the original pixel data, and only the mask will be based on shapes (boolean) 
	pad_size: width of padding in pixels (int or float)    
	crop: whether to crop the raster to the extent of the shapes
    """
    
    with rasterio.open(geotiff) as src:
        # we need to make a temporary shapefile with the same crs as the geotiff 
        shpReprojection(shapefile, src.crs, dst='temp.shp')
        with fiona.open('temp.shp', 'r') as openshape:
            shapes = [feature['geometry'] for feature in openshape]

        out_image, out_transform = rasterio.mask.mask(
            src, 
            shapes, 
            crop=crop, 
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


def tifReprojectionResample(file, reprojected_tif, crs, res, interp, extent_file=None, fill_val=None):
    """
    Reproject, resample, and/or clip raster extent
	file: input raster to resample, reproject, and/or clip (str) (e.g. 'input.tif')
	reprojected_tif: output filename (str) (e.g. 'output_raster.tif')
	crs: target raster coordinate system (str) (e.g. 'EPGS:4326')
	res: target raster resolution (int or float)
	interp: resampling interpolation method (Resampling method from rasterio.warp module)  (e.g. rasterio.warp.Resampling.cubic_spline)
	extent_file: raster file to match extent (str)    
 	fill_val: nodata value to fill all areas not covered by the reprojected source (numeric)
    """

    if extent_file != None:
        tifReprojectionResample(file, 'temp_output.tif', crs, res, interp, extent_file=None)
        file = 'temp_output.tif'
    
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
                dst_nodata=fill_val,
                resampling=interp   # Resampling.cubic_spline
            )
            
        if extent_file != None:
            os.remove('temp_output.tif')


def fillArrayHoles(arr, fillmask, dist=10, iters=1):
    '''
    Fill hole in array
    	arr: array to fill holes (np array)
     	fillmask: mask where to fill holes. hole-filling applied where mask=0 (np array)
      	dist: the maximum number of pixels to search in all directions to find values to interpolate from (int)
        iters: the number of 3x3 smoothing filter passes to run (int)
    '''
    # use rasterio fill fillnodata to fill holes in our array
    inputFilled = fillnodata(arr, mask=fillmask, max_search_distance=dist, smoothing_iterations=iters)
    inputFilled[pd.isnull(inputFilled) == True] = 0
    return inputFilled


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
        filledRaster = fillArrayHoles(inputs, mask=fillmask, dist=dist, iters=iters)

    with rasterio.open(dest, 'w', **profile) as dst:
        dst.write_band(1, filledRaster)


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


def rasterMath(geotiff1, geotiff2, eqn='-', interp=Resampling.cubic_spline, outfile=None):
    '''
    Perform basic math with two raster files, regardless of coordinate system, resolution, and extent
        geotiff1: input raster filename. Note that the second raster will be reprojected and clipped based on this raster (e.g. 'raster1.tif')
        geotiff2: input raster filename (e.g. 'raster2.tif')
        eqn: raster math to perform; addition, subtraction, multiplication, or division. Add 'r' to reverse the order of input files for subtraction or division (to do geotiff2 - geotiff1 or geotiff2 / geotiff1) (e.g. '+', '-', '*', '/', '-r', '/r' OR 'add', 'sub', 'mult', 'div', 'subr', or 'divr')
        interp: resampling interpolation method (Resampling method from rasterio.warp module)  (e.g. rasterio.warp.Resampling.cubic_spline)
        outfile: output file to save raster. Output filename to save. If none, only the array is return (default: None) (str, e.g. 'output.tif')
    '''
    res = rOpen(geotiff1, returnArray=False, returnRes=True, returnCrs=False) # geotiff resolution
    crs = rOpen(geotiff1, returnArray=False, returnRes=False, returnCrs=True) # geotiff coordinate system
    ext = rasterio.open(geotiff1).bounds # geotiff bounds
    arr = rOpen(geotiff1) # geotiff array
    
    res2 = rOpen(geotiff2, returnArray=False, returnRes=True, returnCrs=False) # geotiff resolution
    crs2 = rOpen(geotiff2, returnArray=False, returnRes=False, returnCrs=True) # geotiff coordinate system
    ext2 = rasterio.open(geotiff2).bounds # geotiff bounds
    
    if res == res2 and crs == crs2 and ext == ext2:
        arr2 = rOpen(geotiff2) # if our geotiffs are aligned, we just take the second array
    else: 
        # otherwise we need to reproject and resample our array temporarily
        tifReprojectionResample(geotiff2, 'temp_math_reproj.tif', crs=crs, res=res, 
                                interp=interp, extent_file=geotiff1)
        arr2 = rOpen('temp_math_reproj.tif')
        os.remove('temp_math_reproj.tif')

    # perform the arithmetic
    if eqn == '-' or eqn == 'sub':
        math_arr = np.subtract(arr, arr2)
    elif eqn == '+' or eqn == 'add':
        math_arr = np.add(arr, arr2)
    elif eqn == '*' or eqn == 'mult':
        math_arr = np.multiply(arr, arr2)
    elif eqn == '/' or eqn == 'div':
        math_arr = np.divide(arr, arr2)
    elif eqn == '-r' or eqn == 'subr':
        math_arr = np.subtract(arr2, arr)
    elif eqn == '/r' or eqn == 'divr':
        math_arr = np.divide(arr2, arr)    
        
    if outfile != None: # save as geotiff
        rasterLike(math_arr, outfile, geotiff1)      
    return math_arr


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

def distance_to_shp(geotiff, gpd_shapefile, common_crs='EPSG:32606'):
    '''
    Gets the distance from every pixel to a shapefile
    :param geotiff: input raster file (geotiff as string)
    :param gpd_shapefile: shapefile object to calculate distance from, e.g. centerlines (gpd read shapefile)
    :param common_crs: common crs for geotiff and shapefile data, e.g. EPSG:32606 (coordinate as string)
    :return: array of values with distance from shapefile
    '''
    with rasterio.open(geotiff) as dataset:
        val = dataset.read(1)
        data = [(dataset.xy(x,y)[0],dataset.xy(x,y)[1],val[x,y]) if (np.isnan(val[x,y]) == False) 
                else (np.nan,np.nan,np.nan) for x,y in np.ndindex(val.shape)]
        lon = [i[0] for i in data]
        lat = [i[1] for i in data]
        d = [i[2] for i in data]
        res = pd.DataFrame({"lon":lon,'lat':lat,"data":d})
        res_geodf = gpd.GeoDataFrame(res, geometry=gpd.points_from_xy(res['lon'], res['lat'], crs=dataset.crs))

    # Convert both line and points to the same projected CRS.
    gpd_shapefile_com = gpd_shapefile.to_crs(common_crs)
    res_geodf_com = res_geodf.to_crs(common_crs)

    # The track line shapefile may contain several lines. To measure the shortest distance do a union of all tracks
    single_gpd_shapefile = gpd_shapefile_com.unary_union
    res_geodf_com['distance_m'] = res_geodf_com.distance(single_gpd_shapefile) # Measure the distance
    distances = res_geodf_com['distance_m'].to_numpy()
    distances_array = np.reshape(distances, val.shape)
    return distances_array




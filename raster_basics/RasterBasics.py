import rasterio
import rasterio.plot
from rasterio.plot import show
import rasterio.mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.merge import merge
from rasterio.windows import from_bounds
from rasterio.fill import fillnodata
from rasterio.transform import from_bounds
from rasterio.features import rasterize

import rioxarray
import xarray

import glob, os, warnings
import numpy as np
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import TwoSlopeNorm

import fiona
import shapely.geometry
from shapely.geometry import Point
import geopandas as gpd
import pandas as pd
from pyproj import Transformer, CRS


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
    imshape_ratio = images[0].shape[0]/images[0].shape[1]
    
    # Calculate the number of rows and columns for the subplot grid
    num_rows = 1 if num_images <= ncols else (num_images + ncols - 1) // ncols
    num_cols = min(num_images, ncols)
    
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(10*num_cols, 10*imshape_ratio*num_rows))
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
            fig.colorbar(c, ax=ax, label=ctitle, aspect=4*ax.get_aspect())
    
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
    buffered_data = shape['geometry'].buffer(dist)

    # Create a new GeoDataFrame with the buffered data
    buffered_gdf = gpd.GeoDataFrame(geometry=buffered_data, crs=shape.crs)
    buffered_data.to_file(dst) # save the buffered shapefile


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


def tifReprojectionResample(file, reprojected_tif, crs, res, interp=Resampling.cubic_spline, extent_file=None, fill_val=None):
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


def fillHole(file, dest='output_filled.tif', dist=10, iters=1, mask=None):
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

        if mask == None:
            fillmask = inputs.copy() # fillnodata is applied where the mask=0
            fillmask[inputs>=0] = 1
            fillmask[fillmask!=1] = 0
        else:
            fillmask = mask.copy()

        filledRaster = fillArrayHoles(inputs, fillmask=fillmask, dist=dist, iters=iters)

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


def extract_along_line(xarr, line, n_samples=512, method='nearest'):
    # samples line with n_samples number of points
    profile = []
    dist = []
    for i in range(n_samples):
        point = line.interpolate(i / n_samples - 1., normalized=True) # get next point on the line
        value = xarr.sel(x=point.x, y=point.y, method=method).data # access the nearest pixel in the xarray
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


def end_points(xarr, points, method='nearest'):
    # obtains raster values at our point selections themselves (not interpolated lines)
    point_spot = []
    for p in points:
        point = xarr.sel(x=p[0], y=p[1], method=method).data # access the nearest pixel in the xarray
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


def reproject_velocity(vx_fp, vy_fp, from_epsg, to_epsg, output_vx_fp, output_vy_fp, nan_val=np.nan, apply_correction=False, plot_changes=False, **kwargs):
    '''
    Reproject velocity:
        1. Extract velocity data from geotiff
        2. Convert geotiff to points: get start and end point of velocity vector
        3. Reproject start and end points to new coordinate system
        4. Rasterize reprojected points in new coordinate system
        5. OPTIONAL: apply_correction: Adjust velocity vectors based on magnitude of a direct reprojection
        6. OPTIONAL: plot_changes: Plot velocity vector reproject (contains optional arguments for plotting)

    Arguments:
        vx_fp, vy_fp: Input velocity vectors
        from_epsg, to_epsg: Initial and target coordinate systems in EPSG number (e.g., 32606)
        output_vx_fp, output_vy_fp: Output filepaths
        nan_val: data NaN value to remove (default: np.nan. NOTE: absolute values greater than 1e30 are automatically removed from data)
        apply_correction: whether to apply a correction that matches velocity magnitude with the direct reprojection (default: False)
        plot_changes: show plot of reprojection (default: False)
        **kwargs: Optional arguments passed along for plotting (quiv_scale, quiv_width, title, ctitle, base_arr_color, vmin, vcenter, vmax) 
    '''
        
    # ----------------  1. extract velocity data from existing geotiff ---------------- 
    proj_from = CRS.from_user_input(from_epsg)
    proj_to = CRS.from_user_input(to_epsg)
    with rasterio.open(vx_fp) as vx_src, rasterio.open(vy_fp) as vy_src:
        vx_data = vx_src.read(1)  # Read the vx component
        vy_data = vy_src.read(1)  # Read the vy component
        v_res = vx_src.res
        transform_from = vx_src.transform  # Get the affine transform of the original raster
        
    # ---------------- 2. convert geotiff values to points ---------------- 
    rows, cols = np.indices(vx_data.shape)
    lon, lat = rasterio.transform.xy(transform_from, rows, cols)
    points_from = np.column_stack((np.array(lon).flatten(), np.array(lat).flatten()))

    # calculate the velocity endpoints
    vx_data[np.abs(vx_data) > 1e30] = 0 # remove unrealistic values
    vy_data[np.abs(vy_data) > 1e30] = 0
    if np.isnan(nan_val): # remove NaN values
        vx_data[np.isnan(vx_data)] = 0
        vy_data[np.isnan(vy_data)] = 0
    else:
        vx_data[vx_data == nan_val] = 0
        vy_data[vy_data == nan_val] = 0
    vv_data = (vx_data**2 + vy_data**2)**0.5
    endpoints_from = np.column_stack((points_from[:,0] + vx_data.flatten(), points_from[:,1] + vy_data.flatten()))

    # ---------------- 3. reproject start and end points ---------------- 
    transformer = Transformer.from_crs(proj_from, proj_to, always_xy=True)
    points_to = np.array([transformer.transform(x, y) for x, y in points_from])
    endpoints_to = np.array([transformer.transform(x, y) for x, y in endpoints_from])

    # recalculate the velocity magnitude and components
    dx = endpoints_to[:,0] - points_to[:,0]  # reprojected vx
    dy = endpoints_to[:,1] - points_to[:,1]  # reprojected vy
    angles = np.arctan2(dy, dx)
    mask = ~((np.isnan(dx) | (dx == 0)) & (np.isnan(dy) | (dy == 0))) # mask off-glacier
    
    # ---------------- 4. rasterize points ----------------     
    transform_to, width, height = calculate_default_transform(proj_from, proj_to, vx_data.shape[1], vx_data.shape[0],
                                                              *rasterio.open(vx_fp).bounds, resolution=rasterio.open(vx_fp).res[0])

    pixel_size_x, pixel_size_y = v_res[0], v_res[1]
    x_min, x_max = min(points_to[:,0])-(pixel_size_x/2), max(points_to[:,0])+(pixel_size_x/2)
    y_min, y_max = min(points_to[:,1])-(pixel_size_y/2), max(points_to[:,1])+(pixel_size_y/2)
    
    # raster dimensions
    width = int((x_max - x_min) / pixel_size_x) + 1
    height = int((y_max - y_min) / pixel_size_y) + 1

    # define affine transform
    transform = from_bounds(x_min, y_min, x_max, y_max, width, height)

    # prepare the data for rasterization
    points_vx = [(Point(coord), value) for coord, value in zip(zip(points_to[:,0][mask], points_to[:,1][mask]), dx[mask])]
    points_vy = [(Point(coord), value) for coord, value in zip(zip(points_to[:,0][mask], points_to[:,1][mask]), dy[mask])]

    # rasterize the data
    vx_arr = rasterize(points_vx, out_shape=(height, width), transform=transform, dtype='float32')
    vy_arr = rasterize(points_vy, out_shape=(height, width), transform=transform, dtype='float32')
    vv_arr = np.sqrt(vx_arr**2 + vy_arr**2)

    with rasterio.open(output_vx_fp, 'w', driver='GTiff', height=height, width=width, count=1, dtype='float32', 
                       crs=proj_to, transform=transform) as dst:
        dst.write(vx_arr, 1)
    with rasterio.open(output_vy_fp, 'w', driver='GTiff', height=height, width=width, count=1, dtype='float32', 
                       crs=proj_to, transform=transform) as dst:
        dst.write(vy_arr, 1)

    # ------------------- plotting ----------------------
    if plot_changes:
        quiv_scale = kwargs.get('quiv_scale', 200)
        quiv_width = kwargs.get('quiv_width', 0.005)
        quiv_color= kwargs.get('quiv_color', 'autumn')
    
        title = kwargs.get('title', None)
        ctitle = kwargs.get('ctitle', None)
        base_arr_color = kwargs.get('base_arr_color', 'BrBG')
        base_arr_vmin = kwargs.get('base_arr_vmin', 0)
        base_arr_vcenter = kwargs.get('base_arr_vcenter', 10)
        base_arr_vmax = kwargs.get('base_arr_vmax', 50)
    
        # original plot
        rasterLike(vv_data, 'tmp_vv_plt.tif', vx_fp)
        quiv_plot(points_from[:,0], points_from[:,1], vx_data, vy_data, 1, quiv_scale, quiv_width=quiv_width, quiv_color=quiv_color,
                  base_tiff='tmp_vv_plt.tif', title=title, ctitle=ctitle, base_arr_color=base_arr_color,  
                  base_arr_vmin=base_arr_vmin, base_arr_vcenter=base_arr_vcenter, base_arr_vmax=base_arr_vmax)

        # reprojected plot
        rasterLike(vv_arr, 'tmp_vv_plt.tif', output_vx_fp)
        quiv_plot(points_to[:,0][mask], points_to[:,1][mask], dx[mask], dy[mask], 1, quiv_scale, quiv_width=quiv_width, 
                  quiv_color=quiv_color, base_tiff='tmp_vv_plt.tif', title=title, ctitle=ctitle, base_arr_color=base_arr_color,  
                  base_arr_vmin=base_arr_vmin, base_arr_vcenter=base_arr_vcenter, base_arr_vmax=base_arr_vmax)

    # ---------------- 5. apply correction based on straight reprojection magnitude ----------------
    if apply_correction:
        assert np.abs(pixel_size_x-pixel_size_y) < 0.1, f'Pixels are not square: ({pixel_size_x:.2f}, {pixel_size_y:.2f})'
        # straight reprojection of velocity data
        tifReprojectionResample(vx_fp, 'temp_vx.tif', proj_to, pixel_size_x, Resampling.cubic_spline) # reproject
        tifReprojectionResample(vy_fp, 'temp_vy.tif', proj_to, pixel_size_x, Resampling.cubic_spline)
        vx_arr_tmp = rasterio.open('temp_vx.tif').read(1)
        vy_arr_tmp = rasterio.open('temp_vy.tif').read(1)
        vx_arr_tmp[(np.abs(vx_arr_tmp) >= 1e10) | np.isnan(vx_arr_tmp)] = 0
        vy_arr_tmp[(np.abs(vy_arr_tmp) >= 1e10) | np.isnan(vy_arr_tmp)] = 0
        vv_dir_arr = (vx_arr_tmp**2 + vy_arr_tmp**2)**0.5 # derive velocity magnitude 

        # adjust point-to-point velocity based on magnitude
        mask_mag = (vv_dir_arr != 0) & (vv_arr != 0)
        vx_arr_corr = np.copy(vx_arr)
        vy_arr_corr = np.copy(vx_arr)
        vx_arr_corr[mask_mag] = vx_arr[mask_mag] * vv_dir_arr[mask_mag] / vv_arr[mask_mag]
        vy_arr_corr[mask_mag] = vy_arr[mask_mag] * vv_dir_arr[mask_mag] / vv_arr[mask_mag]
        rasterLike(vx_arr_corr, output_vx_fp, 'temp_vx.tif')
        rasterLike(vy_arr_corr, output_vy_fp, 'temp_vy.tif')
        os.remove('temp_vx.tif')
        os.remove('temp_vy.tif')
        
    # ------------------- plotting ----------------------
    if plot_changes:
        # reprojected plot, with correction
        if apply_correction:
            vv_arr_corr = np.sqrt(vx_arr_corr**2 + vy_arr_corr**2)
            rasterLike(vv_arr_corr, 'tmp_vv_plt.tif', output_vx_fp)
            quiv_plot(points_to[:,0][mask], points_to[:,1][mask], dx[mask], dy[mask], 1, quiv_scale, quiv_width=quiv_width, 
                      quiv_color=quiv_color, base_tiff='tmp_vv_plt.tif', title=title, ctitle=ctitle, base_arr_color=base_arr_color,
                      base_arr_vmin=base_arr_vmin, base_arr_vcenter=base_arr_vcenter, base_arr_vmax=base_arr_vmax)
        
        os.remove('tmp_vv_plt.tif')




def quiv_plot(q1, q2, q3, q4, q5, q6, base_tiff=None, **kwargs):
    """
    Plots a velocity field with arrows and an optional background array.
    
    Parameters:
        q1, q2, q3, q4: Coordinates and components for quiver plot.
        q5: Magnitude of arrows for quiver.
        q6: Scale factor for quiver arrows.
        base_tiff: Geotiff for the background plot (optional).
        **kwargs: Additional parameters for base array customization:
            - `quiv_width`: Quiver arrow width
            - `quiv_color`: Quiver arrow color
            - `title`: Plot title
            - `ctitle`: Colorbar title
            - `base_arr_color`: Colormap for the base array (default is 'RdBu').
            - `base_arr_ext`: Extent of the base array (default is None).
            - `base_arr_vmin`, `base_arr_vcenter`, `base_arr_vmax`: Min, center, and max for colorbar (default is -50, 0, 50)
    """    
    fig = plt.figure()
    ax = fig.add_subplot(111)  
    
    quiv_width = kwargs.get('quiv_width', 0.005)
    quiv_color= kwargs.get('quiv_color', 'autumn')
    ax.quiver(q1, q2, q3, q4, q5, cmap=quiv_color, scale=q6, width=quiv_width)

    title = kwargs.get('title', None)
    ctitle = kwargs.get('ctitle', None)
    base_arr_color = kwargs.get('base_arr_color', 'RdBu')
    base_arr_vmin = kwargs.get('base_arr_vmin', 0)
    base_arr_vcenter = kwargs.get('base_arr_vcenter', 10)
    base_arr_vmax = kwargs.get('base_arr_vmax', 50)
    if base_tiff is not None:
        divnorm = TwoSlopeNorm(vmin=base_arr_vmin, vcenter=base_arr_vcenter, vmax=base_arr_vmax)
        base_ext = [rasterio.open(base_tiff).bounds[0], rasterio.open(base_tiff).bounds[2], 
                    rasterio.open(base_tiff).bounds[1], rasterio.open(base_tiff).bounds[3]]
        im = ax.imshow(rasterio.open(base_tiff).read(1), cmap=mpl.colormaps[base_arr_color], norm=divnorm, extent=base_ext)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, label=ctitle)
    
    ax.set_title(title, weight='bold', pad=10)
    fig.tight_layout(pad=3, w_pad=-3.0, h_pad=2.0)
    plt.show()


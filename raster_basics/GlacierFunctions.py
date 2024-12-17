# Glacier functions. Mass balance, area, slope, aspect, hillshade calculations.
# Also some functions correcting velocities
import os
import rasterio
import rasterio.plot
import rasterio.mask
import earthpy.spatial as es
import numpy as np
import math
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy import stats
from scipy.ndimage import distance_transform_edt
from shapely.geometry import Point
from raster_basics.SmoothingFunctions import dynamicSmoothing, dynamicSmoothingExponential


def glacierArea(outline, res):
    """
    Area of a glacier
	outline: glacier outline, binary array with 1 denoting on-glacier terrain and 0 denoting off-glacier terrain (array-like)
	res: raster resolution (int or float)
    """
    totalArea = np.sum(outline) * res * res  # area in m2
    return totalArea

def totalMassBalance(dhdt, outline, density):
    """
    Total (dhdt) mass balance (kg/m2-yr) = dhdht (m/yr) * density (kg/m3)
	dhdt: elevation change, m/yr (array-like)
	outline: glacier outline, binary array with 1 denoting on-glacier terrain and 0 denoting off-glacier terrain (array-like)
	density: density assumption, kg/cubic.m (int or float)
    """
    massBalanceTotal = dhdt * density * outline
    return massBalanceTotal

def totalMassBalanceValue(totalMB, area, res):
    """
    Total (dhdt) mass balance value
	totalMB: pixel-wise total mass balance, kg/m2-yr (mm w.e.) (array-like)
	area: glacier area, m2 (int or float)
	res: pixel resolution (int or float)
    """
    # totalMB[totalMB == 0] = np.nan
    totalMassBalanceVal = np.nansum(totalMB * res * res) / area  # total mass balance in kg/m2-yr
    return totalMassBalanceVal

def divQ(vx, vy, h, res, vCol=0.8):
    """
    Flux divergence
	vx, vy: x- and y-direction velocity, m/yr (array-like)
	h: glacier thickness, m (array-like)
	res: pixel resolution (int or float)
	vCol: vertical-average velocity scaling factor, typically 0.8 (int or float)
    """
    divQarray = ((np.gradient(vx, res)[1] * h) + (np.gradient(h, res)[1] * vx) +
                  (np.gradient(vy * -1, res)[0] * h) + (np.gradient(h * -1, res)[0] * vy)) * vCol
    return divQarray

def flux_div_comps(vx, vy, h, res, vCol=0.8, filt=True, filter_type='Gauss', filter_factor=4):
    """
    Flux divergence - return each component as well
	vx, vy: x- and y-direction velocity, m/yr (array-like)
	h: glacier thickness, m (array-like)
	res: pixel resolution (int or float)
	vCol: vertical-average velocity scaling factor, typically 0.8 (int or float)
	filt: whether or not to smooth derivatives
	filter_type: type of filter, either 'Gauss' or 'Exp'
	filter_factor: smoothing factor for filtering
    """
    dvx = (np.gradient(vx, res)[1] * h) * vCol
    dvy = (np.gradient(vy * -1, res)[0] * h) * vCol
    dhx = (np.gradient(h, res)[1] * vx) * vCol
    dhy = (np.gradient(h * -1, res)[0] * vy) * vCol
    if filt == True:
        if filter_type == 'Gauss':
            dvx = dynamicSmoothing(dvx, h, res, filter_factor)
            dvy = dynamicSmoothing(dvy, h, res, filter_factor)
            dhx = dynamicSmoothing(dhx, h, res, filter_factor)
            dhy = dynamicSmoothing(dhy, h, res, filter_factor)
        elif filter_type == 'Exp':
            dvx = dynamicSmoothingExponential(dvx, h, res, filter_factor)
            dvy = dynamicSmoothingExponential(dvy, h, res, filter_factor)
            dhx = dynamicSmoothingExponential(dhx, h, res, filter_factor)
            dhy = dynamicSmoothingExponential(dhy, h, res, filter_factor)
    divQ = dvx + dvy + dhx + dhy
    return divQ, dvx, dvy, dhx, dhy


# import richdem as rd
# def glacierAttributes(dem_rast, attrFormat):
#     # return desired attribute (attrFormat is 'slope_degree' or 'slope_percentage', 'slope_riserun', 'aspect')
#     # https://richdem.readthedocs.io/en/latest/python_api.html#richdem.rdarray
#     # gdal.DEMProcessing('slope.tif', dem_rast, 'slope', slopeFormat=slopeFormat)
#     # with rasterio.open('slope.tif') as dataset:
#     #     slope = dataset.read(1)
#     #     slope[slope == -9999.0] = 0
#     no_data_val = rasterio.open(dem_rast).nodatavals[0]
#     if no_data_val == None:
#         dem_array = rd.LoadGDAL(dem_rast, no_data=-9999.0)  # assign a No Data value if none is prescribed
#     else:
#         dem_array = rd.LoadGDAL(dem_rast)
#     attr_array = rd.TerrainAttribute(dem_array, attrib=attrFormat)
#     return attr_array

def glacierSlope(array, res, attrFormat='riserun'):
    """
    Rise/run or degree slope of a raster array
	array: input array (array-like)
	res: pixel resolution (int or float)
 	attrFormat: slope format (str, either 'riserun' or 'degree')
    """
    # alternate way to obtain the riserun slope, with an array input instead of raster
    px, py = np.gradient(array, res)
    if attrFormat == 'riserun':
        slope = np.sqrt(px ** 2 + py ** 2)
    elif attrFormat == 'degree':
        slope = np.degrees(np.arctan2(np.sqrt(px ** 2 + py ** 2), np.ones_like(array)))
    else:
        print("Invalid attrFormat input. Please choose either 'riserun' (default) or 'degree'")
    return slope

def demHillshade(demArray, az=135, alt=30):
    # https://earthpy.readthedocs.io/en/latest/gallery_vignettes/plot_dem_hillshade.html
    # demArray is an input raster DEM opened and read as an array
    # azimuth 90 is east, 0 is north. default value 135 is SE for Alaska
    # altitude is from 0 to 90, 90 being directly overhead
    hillshade = es.hillshade(demArray, azimuth=az, altitude=alt)     # Create and plot the hillshade with earthpy
    return hillshade

def demAspect(dem_array, res):
    """
    Aspect of a DEM array – takes an array input instead of raster
	dem_array: input array (array-like)
	res: pixel resolution (int or float)
    """
    aspect_array = np.zeros_like(dem_array)
    py, px = np.gradient(dem_array, res)

    for i in range(len(aspect_array)):
        for j in range(len(aspect_array[0])):
            pixel_deg = math.atan2(py[i][j], px[i][j]) * 180 / math.pi
            if pixel_deg >= 0 and pixel_deg <= 90:
                pixel_aspect = 90 - pixel_deg
            elif pixel_deg < 0:
                pixel_aspect = abs(pixel_deg) + 90
            elif pixel_deg > 90:
                pixel_aspect = 450 - pixel_deg
            aspect_array[i][j] = pixel_aspect
    return aspect_array


def velFlowlineAspect(linestring, array, arr_extent, mask=None):
    '''
    Use flowline as velocity direction. Returns array of aspects, based on the flowline
    # linestring: center flowline json object. Type: shapely.geometry.linestring.LineString
    # array: 2D numpy array, output aspect will be the same size
    # arr_extent: Array extent in shapely geometry coordinates (left, bottom, right, top) 
    # mask: Glacier outline mask, 2D array (default None considers all values in arry)
    # returns: 2D array of velocity aspect based on flowline
    '''
    aspect_array = np.zeros_like(array)
    x, y = np.meshgrid(np.arange(arr_extent[0], arr_extent[2]+1), np.arange(arr_extent[1], arr_extent[3]+1))
    if mask is None:
        mask = np.ones_like(array)
    
    # iterate through velocity array
    for i in range(len(array)): 
        for j in range(len(array[0])):
            if mask[i][j] == 1:
                point = Point(x[i,j], y[i,j]) # find our point in json coordinates
                closest_point = linestring.interpolate(linestring.project(point)) # find the nearest point on linestring
                coords = np.array(linestring.coords)
                dists = np.linalg.norm(coords - closest_point.coords, axis=1)
                closest_idx = np.argmin(dists) # find the index of the closest point on the line

                if closest_idx == 0: # check if nearest point is endpoint
                    prev_point = coords[closest_idx]
                    next_point = coords[closest_idx + 2]
                elif closest_idx == len(np.array(linestring.coords)) - 1:
                    prev_point = coords[closest_idx - 2]
                    next_point = coords[closest_idx]
                else:
                    prev_point = coords[closest_idx - 1]
                    next_point = coords[closest_idx + 1]

                # angle of perpendicular line -- perpendicular line is -1/m, so we input (-x, y) instead of (y, x)
                    # perpendicular line is important because we shift aspect from normal coordinates (east is 0)
                    # to DEM aspect (north is 0)
                aspect_angle = np.rad2deg(np.arctan2(prev_point[0] - next_point[0], next_point[1] - prev_point[1]))
                aspect_angle %= 360 # adjust range from 0 to 360
                aspect_array[i][j] = aspect_angle # add to array of aspects
    return aspect_array


def weighted_average(value1, value2, weight1, weight2):
    # Convert directions to radians
    value1_rad = np.deg2rad(value1)
    value2_rad = np.deg2rad(value2)

    # Convert directions to complex numbers
    value1_complex = np.exp(1j * value1_rad)
    value2_complex = np.exp(1j * value2_rad)

    # Calculate the weighted average of complex numbers
    average_complex = (weight1 * value1_complex + weight2 * value2_complex) / (weight1 + weight2)

    # Convert the average complex number back to direction in radians
    average_rad = np.angle(average_complex)

    # Convert the average direction to degrees
    average_deg = np.rad2deg(average_rad)

    # Wrap the average direction back to the range of -180 to 180
    if average_deg > 180:
        average_deg -= 360

    return average_deg


def velFlowlineOutlineAspect(linestring, outlinestring, relative_distance_array, arr_extent, mask):
    '''
    Use flowline/outline as velocity direction. Returns array of aspects, based on the flowline and outline.
    The direction is a combination of the flowline/outline direction, weighted by the distance from each
    # linestring: center flowline json object. Type: shapely.geometry.linestring.LineString
    # outlinestring: glacier outline json object. Type: shapely.geometry.linestring.LineString
    # relative_distance_array: 2D numpy array, array of relative distance from centerline to outline (0 to 1)
    # arr_extent: Array extent in shapely geometry coordinates (left, bottom, right, top) 
    # mask: Glacier outline mask, 2D array (default None considers all values in arry)
    # returns: 2D array of velocity aspect based on flowline and outline
    '''
    aspect_array = np.zeros_like(relative_distance_array)
    
    x, y = np.meshgrid(np.arange(arr_extent[0], arr_extent[2]+1), np.arange(arr_extent[1], arr_extent[3]+1))
    if mask is None:
        mask = np.ones_like(relative_distance_array)
    
    # iterate through velocity array
    for i in range(len(relative_distance_array)): 
        for j in range(len(relative_distance_array[0])):
            if mask[i][j] == 1:
                point = Point(x[i,j], y[i,j]) # find our point in json coordinates

                # # aspect of the centerline
                closest_point_cl = linestring.interpolate(linestring.project(point)) # find the nearest point on line
                coords_cl = np.array(linestring.coords)
                dists_cl = np.linalg.norm(coords_cl - closest_point_cl.coords, axis=1)
                closest_idx_cl = np.argmin(dists_cl) # find the index of the closest point on the line

                if closest_idx_cl == 0: # check if nearest point is endpoint
                    prev_point_cl = coords_cl[closest_idx_cl]
                    next_point_cl = coords_cl[closest_idx_cl + 2]
                elif closest_idx_cl == len(np.array(linestring.coords)) - 1:
                    prev_point_cl = coords_cl[closest_idx_cl - 2]
                    next_point_cl = coords_cl[closest_idx_cl]
                else:
                    prev_point_cl = coords_cl[closest_idx_cl - 1]
                    next_point_cl = coords_cl[closest_idx_cl + 1]

                # angle of perpendicular line -- perpendicular line is -1/m, so we input (-x, y) instead of (y, x)
                aspect_angle_cl = np.rad2deg(np.arctan2(prev_point_cl[0] - next_point_cl[0], 
                                                        next_point_cl[1] - prev_point_cl[1]))
                aspect_angle_cl %= 360 # adjust range from 0 to 360

                # # aspect of the outline
                closest_point_ol = outlinestring.interpolate(outlinestring.project(point))
                coords_ol = np.array(outlinestring.coords)
                dists_ol = np.linalg.norm(coords_ol - closest_point_ol.coords, axis=1)
                closest_idx_ol = np.argmin(dists_ol)
                
                if closest_idx_ol == len(np.array(outlinestring.coords)) - 1:
                    prev_point_ol = coords_ol[closest_idx_ol - 1]
                    next_point_ol = coords_ol[0]
                else:
                    prev_point_ol = coords_ol[closest_idx_ol - 1]
                    next_point_ol = coords_ol[closest_idx_ol + 1]
                # angle of perpendicular line -- perpendicular line is -1/m, so we input (-x, y) instead of (y, x)
                aspect_angle_olr = np.arctan2(prev_point_ol[0] - next_point_ol[0],
                                              next_point_ol[1] - prev_point_ol[1])

                # get aspect of the angle and adjust range from 0 to 360
                aspect_angle_ol = np.rad2deg(aspect_angle_olr)
                aspect_angle_ol %= 360
                
                # make sure aspect is in the right direction (the angle should be closer to the centerline angle)
                if abs(aspect_angle_ol - aspect_angle_cl + 180) < abs(aspect_angle_ol - aspect_angle_cl):
                    aspect_angle_ol += 180
                elif abs(aspect_angle_ol - aspect_angle_cl - 180) < abs(aspect_angle_ol - aspect_angle_cl):
                    aspect_angle_ol -= 180
                aspect_angle_ol %= 360

                # weighted average based on relative distance  
                aspect_array[i][j] = weighted_average(aspect_angle_ol, aspect_angle_cl, relative_distance_array[i][j], 
                                                      (1 - relative_distance_array[i][j]))
    return aspect_array


def velocityAspect(vel_x, vel_y):
    # gets the aspect of velocity vectors, given input x- and y-direction velocity arrays
    vel_aspect = np.zeros_like(vel_x)
    for i in range(len(vel_aspect)):
        for j in range(len(vel_aspect[0])):
            pixel_deg = math.atan2(vel_y[i][j], vel_x[i][j]) * 180 / math.pi
            if pixel_deg >= 0 and pixel_deg <= 90:
                pixel_aspect = 90 - pixel_deg
            elif pixel_deg < 0:
                pixel_aspect = abs(pixel_deg) + 90
            elif pixel_deg > 90:
                pixel_aspect = 450 - pixel_deg
            elif math.isnan(pixel_deg) == True:
                pixel_aspect = math.nan
            vel_aspect[i][j] = pixel_aspect
    return vel_aspect

def velocityAspectAngle(vel_aspect1, vel_aspect2):
    # gets the angle between an aspect array and an ARRAY of aspect arrays. returns the maximum
    for i in range(len(vel_aspect2)):
        if i == 1:
            b = a.copy()
        a = np.arccos(np.cos((vel_aspect1 - vel_aspect2[i]) * math.pi / 180)) * 180 / math.pi
        if i > 0:
            b = np.maximum(a, b)
    return b


def velAspectCorrection(dem_aspect, vel_x, vel_y, threshold=90):
    # replace velocities that go against the aspect ('uphill') with 0
    # first, find the aspect based on the velocity (save this in vel_aspect)
    vel_aspect = velocityAspect(vel_x, vel_y)

    # now, compare and replace values where aspect differs beyond the threshold
    vel_x_cor = np.zeros_like(vel_x)
    vel_y_cor = np.zeros_like(vel_y)
    for i in range(len(dem_aspect)):
        for j in range(len(dem_aspect[0])):
            dem_aspect_high = dem_aspect[i][j] + 360
            dem_aspect_low = dem_aspect[i][j] - 360
            # find the locations where we are within the threshold
            if np.array([abs(dem_aspect[i][j] - vel_aspect[i][j]) <= threshold,
                         abs(dem_aspect_high - vel_aspect[i][j]) <= threshold,
                         abs(dem_aspect_low - vel_aspect[i][j]) <= threshold]).any():
                vel_x_cor[i][j] = vel_x[i][j]
                vel_y_cor[i][j] = vel_y[i][j]
            else:       # if we are beyond the threshold, we return 0
                # vel_x_cor[i][j] = vel_x[i][j] * 0.000001
                # vel_y_cor[i][j] = vel_y[i][j] * 0.000001
                vel_x_cor[i][j] = vel_x[i][j] * 0.0000001
                vel_y_cor[i][j] = vel_y[i][j] * 0.0000001
    return vel_x_cor, vel_y_cor

def velAspectDirection(dem_aspect, vel):
    # use velocity magnitudes, but use aspect for the direction
    vel_x_cor = np.zeros_like(vel)
    vel_y_cor = np.zeros_like(vel)
    for i in range(len(dem_aspect)):
        for j in range(len(dem_aspect[0])):
            # convert from aspect and vel magnitude to vx and vy
            vel_x_cor[i][j] = np.sin(dem_aspect[i][j] * math.pi / 180) * vel[i][j]  # this is our new vel_x
            vel_y_cor[i][j] = np.cos(dem_aspect[i][j] * math.pi / 180) * vel[i][j]  # this is our new vel_y
    return vel_x_cor, vel_y_cor

def velAspectSlopeThreshold(dem_aspect, vel_x, vel_y, dem_slope, slope_threshold):
    # use velocity magnitudes, but use aspect for the direction IF SLOPE EXCEEDS A THRESHOLD
    # compare and replace values where aspect differs beyond the threshold
    vel_x_cor = np.zeros_like(vel_x)
    vel_y_cor = np.zeros_like(vel_y)
    vel = np.power((np.square(vel_x) + np.square(vel_y)), 0.5)
    for i in range(len(dem_aspect)):
        for j in range(len(dem_aspect[0])):
            # convert from aspect and vel magnitude to vx and vy
            if dem_slope[i][j] > slope_threshold:
                vel_x_cor[i][j] = np.sin(dem_aspect[i][j] * math.pi / 180) * vel[i][j]  # this is our new vel_x
                vel_y_cor[i][j] = np.cos(dem_aspect[i][j] * math.pi / 180) * vel[i][j]  # this is our new vel_y
            else:
                vel_x_cor[i][j] = vel_x[i][j]
                vel_y_cor[i][j] = vel_y[i][j]
    return vel_x_cor, vel_y_cor

def velAspectSlopeAverage(dem_aspect, vel_x, vel_y, dem_slope, slope_weight):
    # use velocity magnitudes, but calculate direction as an average of aspect and raw data weighted based on slope
    # slope_weight is the slope value where raw data and DEM-based aspect have even weight. higher values attribute
    #   more weight to the raw data based aspect
    # compare and replace values where aspect differs beyond the threshold
    vel_x_cor = np.zeros_like(vel_x)
    vel_y_cor = np.zeros_like(vel_y)
    vel = np.power((np.square(vel_x) + np.square(vel_y)), 0.5)
    for i in range(len(dem_aspect)):
        for j in range(len(dem_aspect[0])):
            # first, find the aspect based on the velocity (save this in vel_aspect)
            pixel_deg = math.atan2(vel_y[i][j], vel_x[i][j]) * 180 / math.pi
            if pixel_deg >= 0 and pixel_deg <= 90:
                pixel_aspect = 90 - pixel_deg
            elif pixel_deg < 0:
                pixel_aspect = abs(pixel_deg) + 90
            elif pixel_deg > 90:
                pixel_aspect = 450 - pixel_deg
            # wieghted aspect of velocity using DEM and raw velocity product direction
            weighted_aspect = np.average([pixel_aspect, dem_aspect[i][j]], weights=[1, dem_slope[i][j]/slope_weight])
            vel_x_cor[i][j] = np.sin(weighted_aspect * math.pi / 180) * vel[i][j]  # this is our new vel_x
            vel_y_cor[i][j] = np.cos(weighted_aspect * math.pi / 180) * vel[i][j]  # this is our new vel_y
    return vel_x_cor, vel_y_cor

def slope_vel_plot(dem_slope, vel, title, slope_threshold=60, showPlot=False):
    # plot the difference in velocity magnitudes vs slope
    # filter out slope values larger than 60, where ice can't stick to, or a input threshold
    dem_slope_masked = np.ma.masked_where(dem_slope > min(slope_threshold, 60), dem_slope)
    vel_masked = np.ma.masked_where(np.ma.getmask(dem_slope_masked), vel)  # apply the mask to the velocity data too
    new_vel = vel_masked / np.cos(dem_slope_masked * math.pi / 180)  # correction for map-velocity vs in-plane velocity

    fig, ax = plt.subplots()
    ax.scatter(dem_slope_masked.flatten(), new_vel.flatten(), s=1, c='r', alpha=0.25)

    m, b = np.polyfit(dem_slope_masked.flatten(), new_vel.flatten(), 1)
    plt.plot(dem_slope_masked.flatten(), m * dem_slope_masked.flatten() + b, color='k', lw=0.5)

    mean, boundary, number = stats.binned_statistic(dem_slope_masked.flatten(), new_vel.flatten(),
                                                    statistic='mean', bins=np.linspace(0, 60, 24))
    ax.hlines(mean, boundary[1:], boundary[:-1], colors='b', alpha=1)

    ax.set_xlabel('Slope (degrees)')
    ax.set_ylabel('Velocity Magnitude (In-Plane) (m/a)')
    ax.set_title(title, weight='bold')
    ax.set_xlim(left=0, right=dem_slope_masked.max())
    ax.set_ylim(bottom=0, top=new_vel.max())
    plt.grid()
    plt.show(block=showPlot)

def h_vel_plot(dem_slope, thickness, vel, title, slope_threshold=60, showPlot=False):
    # plot the difference in velocity magnitudes vs slope
    # filter out slope values larger than 60, where ice can't stick to, or a input threshold
    dem_slope_masked = np.ma.masked_where(dem_slope > min(slope_threshold, 60), dem_slope)
    vel_masked = np.ma.masked_where(np.ma.getmask(dem_slope_masked), vel)  # apply the mask to the velocity data too
    thickness_masked = np.ma.masked_where(np.ma.getmask(dem_slope_masked), thickness)
    new_vel = vel_masked / np.cos(dem_slope_masked * math.pi / 180)  # correction for map-velocity vs in-plane velocity

    fig, ax = plt.subplots()
    ax.scatter(thickness_masked.flatten(), new_vel.flatten(), s=1, c='r', alpha=0.25)

    m, b = np.polyfit(thickness_masked.flatten(), new_vel.flatten(), 1)
    plt.plot(thickness_masked.flatten(), m * thickness_masked.flatten() + b, color='k', lw=0.5)

    mean, boundary, number = stats.binned_statistic(thickness_masked.flatten(), new_vel.flatten(),
                                                    statistic='mean', bins=np.linspace(0, 220, 21))
    ax.hlines(mean, boundary[1:], boundary[:-1], colors='b', alpha=1)

    ax.set_xlabel('Thickness (m)')
    ax.set_ylabel('Velocity Magnitude (In-Plane) (m/a)')
    ax.set_title(title, weight='bold')
    ax.set_xlim(left=0, right=thickness_masked.max())
    ax.set_ylim(bottom=0, top=new_vel.max())
    plt.grid()
    plt.show(block=showPlot)

def stress_vel_plot(dem_slope, thickness, vel, title, slope_threshold=60, showPlot=False):
    # plot the difference in velocity magnitudes vs driving stress (slope * thickness)
    # filter out slope values larger than 60, where ice can't stick to, or a input threshold
    dem_slope_masked = np.ma.masked_where(dem_slope > min(slope_threshold, 60), dem_slope)
    vel_masked = np.ma.masked_where(np.ma.getmask(dem_slope_masked), vel)  # apply the mask to the velocity data too
    thickness_masked = np.ma.masked_where(np.ma.getmask(dem_slope_masked), thickness)
    driving_stress = dem_slope_masked * thickness_masked
    new_vel = vel_masked / np.cos(dem_slope_masked * math.pi / 180)  # correction for map-velocity vs in-plane velocity

    fig, ax = plt.subplots()
    ax.scatter(driving_stress.flatten(), new_vel.flatten(), s=1, c='r', alpha=0.25)

    mean, boundary, number = stats.binned_statistic(driving_stress.flatten(), new_vel.flatten(),
                                                    statistic='mean', bins=np.linspace(0, 4000, 20))
    ax.hlines(mean, boundary[1:], boundary[:-1], colors='b', alpha=1)

    m, b = np.polyfit(driving_stress.flatten(), new_vel.flatten(), 1)
    plt.plot(driving_stress.flatten(), m * driving_stress.flatten() + b, color='k', lw=0.5)

    ax.set_xlabel('Driving Stress (m-deg)')
    ax.set_ylabel('Velocity Magnitude (In-Plane) (m/a)')
    ax.set_title(title, weight='bold')
    ax.set_xlim(left=0, right=driving_stress.max())
    ax.set_ylim(bottom=0, top=new_vel.max())
    plt.grid()
    plt.show(block=showPlot)


def particle_flow(vx_fp, vy_fp, t0_loc, time_steps=range(0,11), num_substep=10000, dist_xy=0.1, xfact=1, yfact=1, verbose=True):
    '''
    Takes velocity x- and y- files and a point location. Returns the location of the point after each time step
        vx_fp, vy_fp: velocity x and y filepaths (str, make sure files are in UTM coordinates)
        t0_loc: starting point coordiantes in UTM coordinates (tuple)
        time_steps: time steps for analysis (list or range() of values, same time unit as velocity)
        num_substep: number of substeps (int, default to 1000)
        dist_xy: distance to move each substep (int, default 1, meters)
        xfact, yfact: multiplication factor for vx and vy (for unit conversion or reversal of its_live vx result)
        verbose: whether to print informational statements (boolean, default to True)
    '''
    if verbose:
        print('Calculating flowpath for', os.path.basename(vx_fp), 'and', os.path.basename(vy_fp), '...')
    current_x, current_y = t0_loc # stake starting point location
    err_tot = 0

    point_x_coords = []
    point_y_coords = []
    res = rasterio.open(vx_fp).res[0] # raster resolution
    break_tf = False

    for time_step in tqdm(time_steps):
        # append the UTM coordinates to lists
        point_x_coords.append(current_x)
        point_y_coords.append(current_y)
        t_step = 0 # start with 0 time each step

        # open the velocity files and extract the velocity values at the current UTM coordinates
        for i in range(num_substep):
            try:
                with rasterio.open(vx_fp) as vx_dataset, rasterio.open(vy_fp) as vy_dataset:
                    vx_value = list(vx_dataset.sample([(current_x, current_y)], indexes=1))[0][0] * xfact
                    vy_value = list(vy_dataset.sample([(current_x, current_y)], indexes=1))[0][0] * yfact

                # move dist_xy meters, then updates velocity vals
                v_value = (vx_value**2 + vy_value**2)**0.5 # velocity magnitude
                dist_xy = max(v_value * 0.01, 0.05) # set cutoff minimum value of 0.05
                current_x += (dist_xy * vx_value / v_value) # normalize vx, vy values to desired distance
                current_y += (dist_xy * vy_value / v_value)
                t_step += dist_xy/v_value # time taken to move that distance

                assert i != num_substep-1, 'WARNING: NOT ENOUGH "NUM_SUBSTEP" OR TOO SMALL "DIST_XY"'
                if t_step >= 0.995: # break when time is 1 step (year)
                    err_step = (t_step - 1) * v_value
                    err_tot += err_step
                    break
            except AssertionError as err:
                print(err)  # print assertation errors
                    
            except:
                print('\tParticle moved off glacier during time_step', time_step)
                break_tf = True
                break

        if break_tf:
            break
    
    if verbose:
        print('\tTotal uncertainty from cumulating time_step residuals (m):', round(err_tot, 4))
        print('\tDONE\n')
    return point_x_coords, point_y_coords


def distance_from_line(line_array, res, mask=None):
    '''
    Calculate the Euclidean distance from centerline array
    :param line_array: Boolean array input to calculate distance from. Distance is calculate from true values
    :param res: Array Euclidean spacing (resoluion) -- numerical
    :param mask: Boolean array mask for on/off glacier terrain (default None takes entire array as on-glacier)
    :return: Array of distance values from line array  
    '''
    
    if mask is None: # if not mask is input, we create a mask that includes all values in array
        mask = np.ones_like(line_array) 
        
    distance_array = distance_transform_edt(np.logical_not(line_array), sampling=[res, res]) # calculate distance
    distance_array[mask == 0] = 0 # remove off-glacier terrain
    return distance_array


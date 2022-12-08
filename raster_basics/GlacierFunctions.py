# Glacier functions. Mass balance, area, slope, aspect, hillshade calculations.
# Also some functions correcting velocities
import rasterio
import rasterio.plot
import rasterio.mask
import earthpy.spatial as es
import richdem as rd
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import stats


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
    totalMassBalanceVal = np.sum(totalMB * res * res) / area  # total mass balance in kg/m2-yr
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


def glacierAttributes(dem_rast, attrFormat):
    # return desired attribute (attrFormat is 'slope_degree' or 'slope_percentage', 'slope_riserun', 'aspect')
    # https://richdem.readthedocs.io/en/latest/python_api.html#richdem.rdarray
    # gdal.DEMProcessing('slope.tif', dem_rast, 'slope', slopeFormat=slopeFormat)
    # with rasterio.open('slope.tif') as dataset:
    #     slope = dataset.read(1)
    #     slope[slope == -9999.0] = 0
    no_data_val = rasterio.open(dem_rast).nodatavals[0]
    if no_data_val == None:
        dem_array = rd.LoadGDAL(dem_rast, no_data=-9999.0)  # assign a No Data value if none is prescribed
    else:
        dem_array = rd.LoadGDAL(dem_rast)
    attr_array = rd.TerrainAttribute(dem_array, attrib=attrFormat)
    return attr_array

def glacierSlope(array, res):
    """
    Rise/run slope of a raster array
	array: input array (array-like)
	res: pixel resolution (int or float)
    """
    # alternate way to obtain the riserun slope, with an array input instead of raster
    px, py = np.gradient(array, res)
    slope = np.sqrt(px ** 2 + py ** 2)
    return slope

def demHillshade(demArray, az=135, alt=30):
    # https://earthpy.readthedocs.io/en/latest/gallery_vignettes/plot_dem_hillshade.html
    # demArray is an input raster DEM opened and read as an array
    # azimuth 90 is east, 0 is north. default value 135 is SE for Alaska
    # altitude is from 0 to 90, 90 being directly overhead
    hillshade = es.hillshade(demArray, azimuth=az, altitude=alt)     # Create and plot the hillshade with earthpy
    return hillshade


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


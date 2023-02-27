# Glacier basic functions. Outline, altitude aggregations, binning calculations
import rasterio
import rasterio.plot
import rasterio.mask
import pyproj
import pandas
import numpy as np
from scipy import stats
from .RasterBasics import shpClip

def glacierOutline(ones_raster, shape, dest):
    """
	Glacier outline
		ones_raster: input array of ones to clip (array-like)
		shape: glacier outline shapefline filepath (str) (e.g. 'outline.shp')
		dest: output filename (str) (e.g. 'output.tif')
    """
    
    shpClip(ones_raster, shape, dest, nan_val=0, fill=True, pad_size=0, crop=False)
    return rasterio.open(dest).read(1)

# def glacierOutline(outlineFile):
#     outlineArray = rasterio.open(outlineFile).read(1)
#     currentGlacierOutline = pandas.DataFrame(outlineArray)  # ice thickness file only has values on the glacier
#     currentGlacierOutline[currentGlacierOutline != 0] = 1  # make mask binary (0 and 1)
#     return currentGlacierOutline.to_numpy()

def altitudeAggregation(calcFile, dem, outline, stat, bin_z=50):
    '''
    Returns elevation-binned mass balance means
    :param calcFile: file on which to compute the statistic (array-like)
    :param dem: DEM for elevation binning (array-like)
    :param outline: glacier outline (array-like)
    :param stat: statistic to be calculated ('mean', 'count', 'sum')
    :param bin_z: interval for binning, which is a DEM height interval or list of elevation boundaries
    :return: array of the altitudinally-binned statistic and elevation bin boundaries
    '''
    demGlacier = np.multiply(dem, outline)
    demGlacierArray = demGlacier[~np.isnan(demGlacier)]
    demGlacier_findMin = demGlacier[demGlacier != 0]
    if str(bin_z).isnumeric() == True:
        z_bin_range = range(int(demGlacier_findMin.min()), int(demGlacier_findMin.max() + bin_z), bin_z)
    else:
        z_bin_range = bin_z

    calcFileArray = np.multiply(calcFile, outline)
    calcFileArray = calcFileArray[~np.isnan(calcFileArray)]

    # bin from dem raster, but operate on massBalance raster
    # returns: stats array for each bin, the bin edges, and the indices of each value in the bin edges
    demBinCount = stats.binned_statistic(demGlacierArray, calcFileArray, statistic='count', bins=z_bin_range)
    demBin_std = stats.binned_statistic(demGlacierArray, calcFileArray, statistic='std', bins=z_bin_range)
    demBin_min = stats.binned_statistic(demGlacierArray, calcFileArray, statistic='min', bins=z_bin_range)
    demBin_max = stats.binned_statistic(demGlacierArray, calcFileArray, statistic='max', bins=z_bin_range)
    demBinStat = stats.binned_statistic(demGlacierArray, calcFileArray, statistic=stat, bins=z_bin_range)
    binStat = demBinStat[0]
    binBoundaries = demBinStat[1]
    binCount = demBinCount[0]
    binNumber = demBinStat[2]
    return binStat, binBoundaries, binCount, binNumber, demBin_std[0], demBin_min[0], demBin_max[0]

def binPercentile(calcFile, dem, outline, lower_percentile, bin_z=None):
    '''
    Returns elevation-binned mass balance percentiles
    :param calcFile: mass balance file (array-like)
    :param dem: DEM (array-like)
    :param outline: glacier outline (array-like)
    :param lower_percentile: lower percentile to be calculated
    :return:
    '''

    demGlacier = dem * outline
    demGlacierArray = demGlacier[~np.isnan(demGlacier)]
    demGlacier_findMin = demGlacier[demGlacier != 0]

    calcFileArray = calcFile * outline
    calcFileArray = calcFileArray[~np.isnan(calcFileArray)]

    # bin from dem raster, but operate on massBalance raster
    # returns: stats array for each bin, the bin edges, and the indices of each value in the bin edges
    demBinPercentile_low = stats.binned_statistic(demGlacierArray, calcFileArray,
                                                  statistic=lambda calcFileArray: np.percentile(calcFileArray,
                                                                                                lower_percentile),
                                                  bins=range(int(demGlacier_findMin.min()),
                                                             int(demGlacier_findMin.max() + bin_z), bin_z))
    demBinPercentile_high = stats.binned_statistic(demGlacierArray, calcFileArray,
                                                   statistic=lambda calcFileArray: np.percentile(calcFileArray,
                                                                                                 100 - lower_percentile),
                                                   bins=range(int(demGlacier_findMin.min()),
                                                              int(demGlacier_findMin.max() + bin_z), bin_z))
    binStat_low = demBinPercentile_low[0]
    binStat_high = demBinPercentile_high[0]
    # binBoundaries = demBinPercentile_low[1]
    # binNumber = demBinPercentile_low[2]
    return binStat_low, binStat_high


def latlonTiffIndex(geotiff, coordinates, crs):
    """
	obtain raster array index values at a lat / lon coordinate location
	takes in a raster geotiff and tuple of lat lon coordinates
		geotiff: raster filepath (str) (e.g. 'raster.tif')
	 	coordinates: coordinates at which to sample the raster (tuple)
		crs: coordiante system of the coordinates (str) (e.g. 'EPSG:4326')
    """
    geo_array = rasterio.open(geotiff)
    transformer = pyproj.Transformer.from_crs(pyproj.CRS('EPSG:4326'), pyproj.CRS(crs))
    east, north = transformer.transform(coordinates[0], coordinates[1])         # find utm coordinates from lat/lon
    row, col = geo_array.index(east, north)                 # obtain raster row/col closest to utm coordinates
    return row, col



## -------- THESE ARE A WORK IN PROGRESS / ABANDONED ---------
def velAdjusted(self, *args):
    # adjust velocity products so they have the same mean values
    # input of 2 or 3 velocity arrays
    num_zeros = []
    for ar in args:
        # get the number of zeros in each array
        num_zeros = np.append(num_zeros, ar.size - np.count_nonzero(ar))

    if len(args) == 2:
        # get velocity arrays
        vel1 = args[0]
        vel2 = args[1]
        # find average values, excluding zero values (off-glacier terrain)
        avg_vel1 = np.sum(vel1) / (vel1.size - min(num_zeros))
        avg_vel2 = np.sum(vel2) / (vel2.size - min(num_zeros))
        # adjust velocity by applying a factor based on other products
        vel1_adj = vel1 * max([avg_vel1, avg_vel2]) / avg_vel1
        vel2_adj = vel2 * max([avg_vel1, avg_vel2]) / avg_vel2
        adjusted_vels = [vel1_adj, vel2_adj]
        return adjusted_vels

    elif len(args) == 3:
        # get velocity arrays
        vel1 = args[0]
        vel2 = args[1]
        vel3 = args[2]
        # find average values, excluding zero values (off-glacier terrain)
        avg_vel1 = np.sum(vel1) / (vel1.size - min(num_zeros))
        avg_vel2 = np.sum(vel2) / (vel2.size - min(num_zeros))
        avg_vel3 = np.sum(vel3) / (vel3.size - min(num_zeros))
        # adjust velocity by applying a factor based on other products
        vel1_adj = vel1 * max([avg_vel1, avg_vel2, avg_vel3]) / avg_vel1
        vel2_adj = vel2 * max([avg_vel1, avg_vel2, avg_vel3]) / avg_vel2
        vel3_adj = vel3 * max([avg_vel1, avg_vel2, avg_vel3]) / avg_vel3
        adjusted_vels = [vel1_adj, vel2_adj, vel3_adj]            # factor by which to multiply velocity
        return adjusted_vels

def velDiscrepancyRegions(self, adjusted_vels, discrepancy_threshold):
    # takes the adjusted velocities and returns a binary mask identifying regions where velocity products differ
    # 'differ' depends on the prescribed discrepancy threshold
    if len(adjusted_vels) == 2:
        # returns boolean array where two arrays are within a tolerance: true are locations of product deviance
        masked_array = np.invert(np.isclose(adjusted_vels[0], adjusted_vels[1], rtol=discrepancy_threshold))
        return masked_array

    elif len(adjusted_vels) == 3:
        # returns locations where arrays are within a tolerance: true are locations of product deviance
        masked_array1 = np.invert(np.isclose(adjusted_vels[0], adjusted_vels[1], rtol=discrepancy_threshold))
        masked_array2 = np.invert(np.isclose(adjusted_vels[0], adjusted_vels[2], rtol=discrepancy_threshold))
        masked_array3 = np.invert(np.isclose(adjusted_vels[1], adjusted_vels[2], rtol=discrepancy_threshold))
        # returns array with locations where ANY two velocity products differ
        masked_array = ((masked_array1 == masked_array2) & (masked_array1 == masked_array3) &
                        (masked_array2 == masked_array3))
        return masked_array

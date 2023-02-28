import numpy as np
import scipy.signal
import scipy.ndimage
from scipy.stats import binned_statistic
import itertools
import rasterio
from .RasterBasics import points_along_lines


def sgolay2d(z, window_size, order, derivative=None):
    '''
    Savitzky-Golay smoothing filter for eliminating high frequency noise from data via moving average technique
    Type of low-pass filter
    :param z: input file array
    :param window_size: odd number: window size of smoothing
    :param order: order of filtering polynomial (cannot be greater than window_size squared)
    :param derivative: calculate the derivative as well: gradient. entries can be 'col', 'row', or 'both'
    :return: smoothed signal
    '''
    # number of terms in the polynomial expression
    n_terms = (order + 1) * (order + 2) / 2.0

    if window_size % 2 == 0:
        raise ValueError('window_size must be odd')
    if window_size ** 2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2

    # exponents of the polynomial.
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ...
    # this line gives a list of two item tuple. Each tuple contains
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [(k - n, n) for k in range(order + 1) for n in range(k + 1)]

    # coordinates of points
    ind = np.arange(-half_size, half_size + 1, dtype=np.float64)
    dx = np.repeat(ind, window_size)
    dy = np.tile(ind, [window_size, 1]).reshape(window_size ** 2, )

    # build matrix of system of equation
    A = np.empty((window_size ** 2, len(exps)))

    for i, exp in enumerate(exps):
        A[:, i] = (dx ** exp[0]) * (dy ** exp[1])

    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2 * half_size, z.shape[1] + 2 * half_size
    Z = np.zeros((new_shape))
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] = band - np.abs(np.flipud(z[1:half_size + 1, :]) - band)
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band + np.abs(np.flipud(z[-half_size - 1:-1, :]) - band)
    # left band
    band = np.tile(z[:, 0].reshape(-1, 1), [1, half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs(np.fliplr(z[:, 1:half_size + 1]) - band)
    # right band
    band = np.tile(z[:, -1].reshape(-1, 1), [1, half_size])
    Z[half_size:-half_size, -half_size:] = band + np.abs(np.fliplr(z[:, -half_size - 1:-1]) - band)
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z
    # top left corner
    band = z[0, 0]
    Z[:half_size, :half_size] = band - np.abs(np.flipud(np.fliplr(z[1:half_size + 1, 1:half_size + 1])) - band)
    # bottom right corner
    band = z[-1, -1]
    Z[-half_size:, -half_size:] = band + np.abs(np.flipud(np.fliplr(z[-half_size - 1:-1, -half_size - 1:-1])) - band)
    # top right corner
    band = Z[half_size, -half_size:]
    Z[:half_size, -half_size:] = band - np.abs(np.flipud(Z[half_size + 1:2 * half_size + 1, -half_size:]) - band)
    # bottom left corner
    band = Z[-half_size:, half_size].reshape(-1, 1)
    Z[-half_size:, :half_size] = band - np.abs(np.fliplr(Z[-half_size:, half_size + 1:2 * half_size + 1]) - band)

    # solve system and convolve
    if derivative == None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -c, mode='valid')
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid')
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid'), scipy.signal.fftconvolve(Z, -c, mode='valid')

def gaussianFilter(z, st_dev, truncate=4):
    '''
    Gaussian smoothing filter
    :param z: input file (open, read raster file; array-like)
    :param st_dev: standard deviation for filter. higher value indicates greater blur
    :return: new, smoothed array
    '''
    pixel_w = 2 * int(truncate * st_dev + 0.5) + 1      # compute pixel window width
    #print('Gaussian Filter Window Width is', pixel_w)
    gaussian = scipy.ndimage.filters.gaussian_filter(z, st_dev, truncate=truncate)
    return gaussian

def dynamicSmoothing(z, window_weight, pixel_size, f):
    '''
    Gaussian filter, dynamic smoothing with changing pixel window size
    :param z: input file (array)
    :param window_weight: weight of window size (thickness file)
    :param pixel_size: size of pixel (resolution)
    :param f: factor of window size from window weight (4 for vx, vy, h. 1 for divQ)
    :return: array of smoothed values
    '''

    # Apply a border of zeros based on the half size
    half_size_max = int(np.max(window_weight) / pixel_size) * f
    pad = np.max(half_size_max) + 2
    r = z.shape[0]
    c = z.shape[1]
    lr_border = np.zeros((r, pad), dtype=int)                    # left/right border
    tb_border = np.zeros((pad, c + (2*pad)), dtype=int)          # top/bottom border
    x = z.copy()
    x = np.concatenate((lr_border, x), axis=1)
    x = np.concatenate((x, lr_border), axis=1)
    x = np.concatenate((tb_border, x), axis=0)
    x = np.concatenate((x, tb_border), axis=0)

    smooth_file = np.zeros(z.shape)         # initiate array to store each smooth pixel
    for r in range(z.shape[0]):
        for c in range(z.shape[1]):
            if window_weight[r, c] != 0:
                # half of window size (must be at least f * pixels)
                half_size = int(f * max(window_weight[r, c], pixel_size) / pixel_size)
                # add the border size to start at actual values
                pixel_window = x[r+pad-half_size:r+pad+half_size+1, c+pad-half_size:c+pad+half_size+1]

                # Gaussian Filter Params
                truncate = 3
                st_dev = (half_size - 0.5) / truncate          # filter st.dev depends on desired window size at a pixel
                smooth_array = scipy.ndimage.filters.gaussian_filter(pixel_window, st_dev, truncate=truncate)
                smooth_file[r, c] = smooth_array[half_size, half_size]      # select the center value as the smoothed value
            else:
                smooth_file[r, c] = 0           # we don't want to smooth 0 regions so we leave these as 0
    return smooth_file
    
def smoothingCorrection(orig_data, smooth_data, centerlines):
    '''
    Apply a correction to smoothed data products, such that smoothing does not decrease the mean values
    Correction is based on glacier centerline values
    :param orig_data: originial data, e.g. unsmoothed velocity or ice thickness (raster)
    :param smooth_data: smoothed data, e.g. output after applying Gaussian filter (raster)
    :param centerlines: centerline(s) for the data, e.g. from OGGM (shapefile)
    :return: smooth array that has been corrected (array-like) and scaling factor (float)
    '''
    orig_out = points_along_lines(orig_data, centerlines) # get values along all centerlines
    smooth_out = points_along_lines(smooth_data, centerlines)
    orig_vals = np.array(list(itertools.chain.from_iterable(orig_out[0]))) # combine centerline values to one array
    smooth_vals = np.array(list(itertools.chain.from_iterable(smooth_out[0])))
    
    orig_vals_mean = np.mean(orig_vals) # find the mean of all centerline values
    smooth_vals_mean = np.mean(smooth_vals)
    scaling = orig_vals_mean/smooth_vals_mean # scaling factor between smooth and raw data based on centerline mean
    smooth_data_scaled = rasterio.open(smooth_data).read(1)*scaling
    return smooth_data_scaled, scaling

def distance_scaling_correction(data_vals1, data_vals2, distance, polyfit_order=2, filter_std=None):
    '''
    Correct smoothed data products based on a scaling from relative distance from centerline
    :param data_vals1: input array/data product 1, e.g. raw data (array-like)
    :param data_vals2: input data product 2--product to be corrected, e.g. smoothed data (array-like)
    :param distance: relative distance of every point from feature(s), e.g. centerlines (array-like)
    :param polyfit_order: order of the polynomial used to fit data, default 2 (int)
    :param filter_std: filtering outliers based on st.dev, default None (numeric)
    :return:
    '''
    x = distance.flatten()
    y_ratio = np.divide(data_vals1, data_vals2)
    y = y_ratio.flatten()

    # Filter outliers, if desired
    if filter_std != None:
        mean = np.nanmean(y)
        std = np.nanstd(y)
        y[y > mean + (3*std)] = np.nan
        
    # Find binned average values
    stat, edges, __ = binned_statistic(x, y, statistic=np.nanmean, bins=30, range=(x.min(), x.max()))
    count, __, __ = binned_statistic(x, y, statistic='count', bins=30, range=(x.min(), x.max()))
    bin_centers = edges[1:] - (abs(edges[0] - edges[1]))/2
    
    # Fit polynomial to values and use polynomial to derive values at each point
    polyfit = np.polyfit(bin_centers[1:], stat[1:], polyfit_order)
    scale_array = 0
    for n in range(polyfit_order+1):
        scale_array += (polyfit[::-1][n]*(distance**n))
    
    # Correct input data array
    data_vals2_corr = data_vals2*scale_array
    return data_vals2_corr


import numpy as np


def add_grf_noise(array, array_res, correlation_length, noise_strength, mask_val='nan', thresh=None):
    '''
    Returns the array with spatially-correlated noise added from a Gaussian Random Field
    :param array: np array that we want to add noise to (2D np.array)
    :param array_res: array resolution (numeric, represents meters)
    :param correlation_length: length of spatial correlation in noise (numeric, represents meters)
    :param noise_strength: magnitude of noise. This is multiplied by a GRF with zero mean and unit variance (numeric)
    :param mask_val: off-glacier value to mask. By default, function masks np.nan. Otherwise, this should be numeric
    :param thresh: threshold for maximum pixel-based noise relative to original value (default None, otherwise value in range [0, 1])
    '''
    if mask_val == 'nan':
        mask = np.isnan(array) # off-glacier mask (nan values)
    else:
        mask = np.equal(array, mask_val) # off-glacier mask *(number; e.g. 0)
    correlation_length_eff = correlation_length / array_res # effective correlation length (pixel length)
    
    # --- Generate the Gaussian Random Field ---
    distances = np.meshgrid(*[np.arange(dim) for dim in array.shape], indexing='ij') # Generate a grid of distances
    squared_distances = sum(dist**2 for dist in distances)

    covariance = np.exp(-0.5 * squared_distances / correlation_length_eff**2) # Compute the covariance function
    random_field = np.random.normal(size=array.shape) # Generate a random field with zero mean and unit variance
    fft_field = np.fft.fftn(random_field) # Perform a Fast Fourier Transform
    filtered_fft_field = np.fft.fftn(covariance**0.5) * fft_field # sq.rt of the covariance function in Fourier space
    grf_noise = np.fft.ifftn(filtered_fft_field).real # Perform an inverse Fast Fourier Transform for GRF
    # --- --- --- --- --- --- --- --- --- --- ---
    
    grf_noise_norm = (grf_noise - np.mean(grf_noise)) / np.std(grf_noise) # normalize GRF (zero mean, unit variance)

    scaled_grf = noise_strength * grf_noise_norm # Scale the noise field by the desired noise strength
    if thresh != None:
        # calculate the pixel-wise noise threshold from the data, then clip the noise where it exceeds the threshold
        max_threshold = thresh * array
        min_threshold = -thresh * array
        scaled_grf = np.clip(scaled_grf, min_threshold, max_threshold)
        
    grf_array_noise = np.where(mask, array, array + scaled_grf) # Add the noise to vel array
    array_noisy = np.where(grf_array_noise < 0, 0, grf_array_noise) # Ensure that velocity values remain positive
    return array_noisy


def add_ice_thickness_bias(ice_thickness, relative_distance, comp_factor=2, direction='UV'):
    """
    Add bias to the ice thickness data based on the relative distance from the glacier center. This assumes bias in ice thickness 
    relative to distance from centerline resulting from inaccurate bed profile assumption during ice thickness inversion

    :param ice_thickness: array of ice thickness values. (np array)
    :param relative_distance: array of relative distances from the glacier centerline to edge (np array)
    :param comp_factor: comparison factor for the V to U-shaped cross-sections. Default of 2 gives max of ~30% thickness change. Lower values add greater bias (numeric)
    :param direction: direction of bias, either 'UV' or 'VU'. 'UV' assumes our modeled ice thickness is a parabola and reality is a V. 'VU' assumes the opposite case (str)
    
    Returns: array of ice thickness values with added bias, and the bias factor
    """

    # apply bias based on the relative distance 
    # we use a 2nd order parabolic function, following from the shallow-ice approximation (consistent with OGGM)
    # we use a V-shape function with the same mean as the parabola
    # bias_factor = ((relative_distance ** 4) + (3/5) + comp_factor) / ((8/5) * relative_distance + comp_factor)
    bias_factor = ((relative_distance ** 2) + (1/3) + comp_factor) / ((4/3) * relative_distance + comp_factor)

    # Add the bias to the ice thickness data
    if direction == 'VU':
        bias_factor = 1 / bias_factor
    
    ice_thickness_bias = ice_thickness * bias_factor

    return ice_thickness_bias, bias_factor


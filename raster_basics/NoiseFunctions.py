import numpy as np


def add_grf_noise(array, array_res, correlation_length, noise_strength):
    '''
    Returns the array with spatially-correlated noise added from a Gaussian Random Field
    :param array: np array that we want to add noise to (2D np.array)
    :param array_res: array resolution (numeric, represents meters)
    :param correlation_length: length of spatial correlation in noise (numeric, represents meters)
    :param noise_strength: magnitude of noise. This is multiplied by a GRF with zero mean and unit variance (numeric)
    '''
    mask = np.isnan(array) # off-glacier mask (nan values)
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
    grf_array_noise = np.where(mask, array, array + scaled_grf) # Add the noise to vel array
    array_noisy = np.where(grf_array_noise < 0, 0, grf_array_noise) # Ensure that velocity values remain positive
    return array_noisy


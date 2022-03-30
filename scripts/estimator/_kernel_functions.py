import numpy as np

VALID_KERNELS = ['quartic','linear','exponential','tophat', 'gaussian']


def get_kernel_function(kernel_name):
    """
    Parameters
    ----------
    kernel_name : str
    Name of kernel function to return
    Valid kernels are ['quartic','linear','exponential', 'tophat', 'gaussian']
    """
    if kernel_name not in VALID_KERNELS:
        raise ValueError('Invalid kernel function:', kernel_name)
    return _kernels_dict[kernel_name]

def tophat(distances, bandwidth):
    return np.where(distances <= bandwidth, 1/2, 0)

def gaussian(distances, bandwidth):
    return (1./np.sqrt(2*np.pi)) * (np.e ** (-(1/2  * ((distances / bandwidth) ** 2))))

def quartic_gaussian_distance(distances, bandwidth):
    return np.where(distances <= bandwidth, 15. / 16 * (1 - (distances / bandwidth) ** 2) ** 2, 0)

def linear_kernel(distances, bandwidth):
    return np.where(distances <= bandwidth , (bandwidth - distances) / (bandwidth ** 2), 0)

def exponential_decay(distances, bandwidth):
    return (1./bandwidth) * np.exp(-(distances/bandwidth))

_kernels_dict = {
                'quartic': quartic_gaussian_distance,
                'linear': linear_kernel,
                'exponential': exponential_decay,
                'tophat': tophat, 
                'gaussian': gaussian,
                }

import numpy as np

from py_vbc.config import Cosmology


def hubble(z, cosmology: Cosmology):
    """Calculate the dimensionless Hubble parameter at redshift *z*.

    :param z: redshift at which to evaluate the Hubble parameter
    :param cosmology: validated cosmological parameters
    """
    return np.sqrt(
        (1 - (cosmology.omega_m + cosmology.omega_r))
        + cosmology.omega_m * (1 + z) ** 3
        + cosmology.omega_r * (1 + z) ** 4
    )

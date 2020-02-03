import numpy as np

from py_vbc.constants import *

def hubble(z):
    """Calculates the value of the Hubble parameter at redshift z in terms
    of H_0, given omega=1. 

    :param z: redshift to evaluate H at
    """
    return np.sqrt((1-(omega_m + omega_r)) + (omega_m*(1+z)**3) +
                   (omega_r*(1+z)**4))

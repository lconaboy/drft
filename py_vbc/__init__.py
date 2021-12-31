import numpy as np

from py_vbc.spectra import *
from py_vbc.constants import *
from py_vbc.derivatives import calc_derivs

"""
A Python version of the vbc_transfer module written by Matt
McQuinn & Ryan O'Leary (astro-ph/1204.1344). Almost identical save for
the removal of some obsolete features and Python optimisation.

TODO

- should costh really be a constant input at runtime? or more
physically motivated, perhaps even random 
"""

def run_pyvbc(vbc, zstart, zend, dz, kmin=1.0, kmax=1.0e3, n=64, delta=False, verbose=False):
    """
    Runs py_vbc and returns either the power spectrum or dimensionless power
    spectrum.

    :param vbc:
        (float)
        Magnitude of v_bc at redshift zstart
    :param zstart:
        (float)
        Redshift to start integrating the evoltuion equations form, e.g. z~1000
        N.B. you will need the appropriate transfer functions at zstart, zstart+dz
        and zstart-dz
    :param zend:
        (float)
        Redshift to evolve the equations to, e.g. the start of the simulation
    :param dz:
        (float)
        Interval in redshift to calculate the derivative of the transfer functions
        over
    :param kmin:
        (float)
        Minimum k to calculate the evolution for (units of Mpc^-1)
    :param kmax:
        (float)
        Maximum k to calculate the evolution for (units of Mpc^-1)
    :param n:
        (int)
        Number of k-values to calculate for, will be equally distributed in log_10
        space between k_min and k_max
    :param delta:
        (bool)
        If True, will return the dimensionless power spectrum, if False will return
        the usual power spectrum
    """
    k = np.logspace(np.log10(kmin), np.log10(kmax), num=n)

    g = calc_derivs(k, vbc, zstart, zend, dz, verbose=verbose)

    p_c, p_b, p_vc, p_vb = calc_power_spec(k, g, zstart)

    if delta is False:
        return k, (p_c, p_b, p_vc, p_vb)

    elif delta is True:
        d_c = calc_delta(k, p_c)
        d_b = calc_delta(k, p_b)
        d_vc = calc_delta(k, p_vc)
        d_vb = calc_delta(k, p_vb)

        return k, (d_c, d_b, d_vc, d_vb)


def run_tests():
    # from py_vbc.tests.itf_test import itf_test
    # from py_vbc.tests.irf_test import irf_test
    # from py_vbc.tests.g_test import g_test
    # from py_vbc.tests.p_test import p_test
    # from py_vbc.tests.pv_test import pv_test
    # from py_vbc.tests.g_test import ratio_test
    from py_vbc.tests.bias_test import bias_test

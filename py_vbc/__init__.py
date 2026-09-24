import os
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

def run_pyvbc(vbc, zstart, zend, dz, k=None, kmin=1.0, kmax=1.0e3, n=64, delta=False, verbose=False, transfer=False, isothermal=False):
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
    :param k:
        (array-like, optional)
        Array of k values (units of Mpc^-1). If provided, kmin, kmax, and n are ignored.
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
    if k is None:
        # k = np.logspace(np.log10(kmin), np.log10(kmax), num=n)
        lkmi = np.log10(kmin)
        lkma = np.log10(kmax)
        dlk = (lkma - lkmi) / float(n - 1.0)
        k = 10.0 ** (np.arange(n, dtype=float) * dlk + lkmi)
    else:
        k = np.asarray(k, dtype=float)

    g = calc_derivs(k, vbc, zstart, zend, dz, verbose=verbose, isothermal=isothermal)


    if transfer:
        t_c, t_b, t_vc, t_vb = calc_tf(k, g, zstart)
        return k, (t_c, t_b, t_vc, t_vb)
        
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
    """Run the py_vbc test suite using pytest."""
    import pytest
    test_file = os.path.join(os.path.dirname(__file__), 'tests', 'test_vbc.py')
    return pytest.main(['-v', test_file])

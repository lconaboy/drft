import os
from pathlib import Path

import numpy as np

from py_vbc.config import ConfigError, RuntimeConfig, load_config
from py_vbc.constants import (
    ASB,
    BOLTZK,
    CRITDENSITY,
    DELB,
    DELC,
    DELT,
    DELTAB,
    DELTAC,
    DELTAT,
    IMAG,
    LIGHTSPEED,
    MEL_EV,
    MELECTRON,
    MPCTOCM,
    MPCTOKM,
    MPROTON,
    MUB,
    REAL,
    SIGMAT,
    TCMB,
    T_GAMMA,
    U_CMB,
    VEL,
    VELIMAG,
    VELREAL,
    YHE,
    verbose,
)
from py_vbc.derivatives import calc_derivs
from py_vbc.spectra import calc_delta, calc_norm, calc_power_spec, calc_tf, sigma


def run_pyvbc(
    vbc,
    zstart,
    zend,
    dz,
    config_file: str | Path | RuntimeConfig,
    *,
    k=None,
    kmin=1.0,
    kmax=1.0e3,
    n=64,
    delta=False,
    verbose=False,
    transfer=False,
    isothermal=False,
):
    """Run py_vbc and return power or transfer functions.

    :param vbc: magnitude of the baryon--dark-matter relative velocity at
        ``zstart`` in km/s
    :param zstart: redshift at which transfer-function initial conditions are set
    :param zend: redshift at which to evaluate the evolved fields
    :param dz: half-width of the redshift interval used to differentiate TFs
    :param config_file: YAML path or preloaded RuntimeConfig selected at runtime
    :param k: optional wavenumber array in Mpc^-1
    :param kmin: minimum wavenumber when ``k`` is not provided
    :param kmax: maximum wavenumber when ``k`` is not provided
    :param n: number of logarithmically spaced wavenumbers
    :param delta: return dimensionless power spectra when true
    :param verbose: print integration progress when true
    :param transfer: return transfer functions instead of power spectra
    :param isothermal: omit baryon sound-speed coupling when true
    """
    if k is None:
        log_kmin = np.log10(kmin)
        log_kmax = np.log10(kmax)
        delta_log_k = (log_kmax - log_kmin) / float(n - 1.0)
        k = 10.0 ** (np.arange(n, dtype=float) * delta_log_k + log_kmin)
    else:
        k = np.asarray(k, dtype=float)

    if isinstance(config_file, RuntimeConfig):
        config = config_file
    else:
        config = load_config(config_file)
    growth = calc_derivs(
        k,
        vbc,
        zstart,
        zend,
        dz,
        config,
        verbose=verbose,
        isothermal=isothermal,
        require_z0=not transfer,
    )

    if transfer:
        values = calc_tf(k, growth, zstart, config)
        return k, values

    power = calc_power_spec(k, growth, zstart, config)
    if delta:
        return k, tuple(calc_delta(k, component) for component in power)

    return k, power


def run_tests():
    """Run the py_vbc test suite using pytest."""
    import pytest

    tests_dir = os.path.join(os.path.dirname(__file__), "tests")
    return pytest.main(["-v", tests_dir])


__all__ = [
    "ConfigError",
    "RuntimeConfig",
    "load_config",
    "run_pyvbc",
    "run_tests",
    "verbose",
    "MPCTOKM",
    "MPCTOCM",
    "DELTAC",
    "DELTAB",
    "DELTAT",
    "REAL",
    "IMAG",
    "VELREAL",
    "VELIMAG",
    "DELC",
    "DELB",
    "VEL",
    "DELT",
    "TCMB",
    "MPROTON",
    "MELECTRON",
    "MEL_EV",
    "LIGHTSPEED",
    "CRITDENSITY",
    "SIGMAT",
    "YHE",
    "BOLTZK",
    "MUB",
    "ASB",
    "U_CMB",
    "T_GAMMA",
    "sigma",
    "calc_norm",
    "calc_power_spec",
    "calc_tf",
    "calc_delta",
]

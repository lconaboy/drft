import os
from pathlib import Path

import numpy as np

from py_vbc.config import ConfigError, EvolutionConfig, RuntimeConfig, load_config
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
    zstart=None,
    zend=None,
    dz=None,
    config_file: str | Path | RuntimeConfig | None = None,
    *,
    k=None,
    kmin=1.0,
    kmax=1.0e3,
    n=64,
    delta=False,
    verbose=False,
    transfer=False,
    isothermal=False,
    mu=None,
):
    """Run py_vbc and return power or transfer functions.

    :param vbc: magnitude of the baryon--dark-matter relative velocity at
        ``zstart`` in km/s
    :param zstart: optional override for the redshift at which transfer-function
        initial conditions are set; defaults to ``evolution.zstart``
    :param zend: optional override for the output redshift; defaults to
        ``evolution.zend`` and is required if that configuration value is null
    :param dz: optional override for the transfer-function derivative offset;
        defaults to ``evolution.dz``
    :param config_file: YAML path or preloaded RuntimeConfig selected at runtime
    :param k: optional wavenumber array in Mpc^-1
    :param kmin: minimum wavenumber when ``k`` is not provided
    :param kmax: maximum wavenumber when ``k`` is not provided
    :param n: number of logarithmically spaced wavenumbers
    :param delta: return dimensionless power spectra when true
    :param verbose: print integration progress when true
    :param transfer: return transfer functions instead of power spectra
    :param isothermal: omit baryon sound-speed coupling when true
    :param mu: directional cosine(s) between the relative velocity and each
        wavenumber; defaults to ``cosmology.mu_isotropic``
    """
    if config_file is None:
        raise ValueError("run_pyvbc requires a py_vbc runtime configuration")
    if isinstance(config_file, RuntimeConfig):
        config = config_file
    else:
        config = load_config(config_file)

    evolution = config.evolution
    zstart = evolution.zstart if zstart is None else float(zstart)
    zend = evolution.zend if zend is None else float(zend)
    if zend is None:
        raise ConfigError(
            "an evolution endpoint is required for run_pyvbc; set evolution.zend "
            "or pass zend explicitly"
        )
    dz = evolution.dz if dz is None else float(dz)
    if dz <= 0:
        raise ValueError("dz must be greater than zero")

    if k is None:
        log_kmin = np.log10(kmin)
        log_kmax = np.log10(kmax)
        delta_log_k = (log_kmax - log_kmin) / float(n - 1.0)
        k = 10.0 ** (np.arange(n, dtype=float) * delta_log_k + log_kmin)
    else:
        k = np.asarray(k, dtype=float)

    if mu is None:
        mu = config.cosmology.mu_isotropic
    mu = np.asarray(mu, dtype=float)
    if mu.ndim == 0:
        mu = np.full(k.shape, float(mu))
    if mu.shape != k.shape:
        raise ValueError("mu must be scalar or have the same shape as k")
    if not np.all(np.isfinite(mu)) or np.any(np.abs(mu) > 1.0 + 1e-12):
        raise ValueError("mu values must be finite and lie in [-1, 1]")

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
        mu=mu,
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
    "EvolutionConfig",
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

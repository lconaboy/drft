from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d

from py_vbc.config import RuntimeConfig


_TF_INDEX = {"c": 1, "b": 2, "g": 3, "t": 6, "vc": 10, "vb": 11}


def _load_transfer_function(path: Path) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(f"configured transfer function does not exist: {path}")
    return np.loadtxt(path)


def interpolate_tf(flag: str, z: float, config: RuntimeConfig, dz=None):
    """Interpolate CAMB transfer functions in wavenumber at redshift *z*.

    When *dz* is supplied, also calculate a centred finite-difference transfer
    function in redshift.

    :param flag: one of ``c``, ``b``, ``g``, ``t``, ``vc``, or ``vb``
    :param z: redshift at which to interpolate
    :param config: runtime configuration containing explicit TF paths
    :param dz: half-width of the redshift interval for the derivative
    :returns: a spline, or a pair of splines when ``dz`` is supplied
    """
    if flag not in _TF_INDEX:
        raise ValueError("flag must be one of 'c', 'vc', 'b', 'vb', 'g', or 't'")

    i_flag = _TF_INDEX[flag]
    z = int(z)
    tf_path = config.transfer_function_path(z)
    tfz = _load_transfer_function(tf_path)

    # CAMB outputs k in h/Mpc, while py_vbc's public API uses Mpc^-1.
    kh = tfz[:, 0] * config.cosmology.h
    tfz_spline = interp1d(kh, tfz[:, i_flag], kind="cubic")

    if dz is None:
        return tfz_spline

    if dz <= 0:
        raise ValueError("dz must be greater than zero")

    tfzm = _load_transfer_function(config.transfer_function_path(int(z - dz)))
    tfzp = _load_transfer_function(config.transfer_function_path(int(z + dz)))

    if not (tfz.shape == tfzm.shape == tfzp.shape):
        raise ValueError("input transfer functions have different lengths")

    dtf = (tfzm - tfzp) / (2.0 * dz)
    dtf_spline = interp1d(kh, dtf[:, i_flag], kind="cubic")
    return tfz_spline, dtf_spline


def interpolate_tf2d(flag: str, zs, config: RuntimeConfig):
    """Interpolate CAMB transfer functions in wavenumber and redshift.

    :param flag: one of ``c``, ``b``, ``g``, ``t``, ``vc``, or ``vb``
    :param zs: redshifts to interpolate
    :param config: runtime configuration containing explicit TF paths
    :returns: a two-dimensional spline in ``(k, z)``
    """
    if flag not in _TF_INDEX:
        raise ValueError("flag must be one of 'c', 'vc', 'b', 'vb', 'g', or 't'")

    redshifts = np.asarray(zs)
    if redshifts.ndim != 1 or len(redshifts) < 4:
        raise ValueError(
            "zs must be a one-dimensional array with at least four entries"
        )

    i_flag = _TF_INDEX[flag]
    if np.any(np.diff(redshifts) <= 0):
        raise ValueError("zs must be strictly increasing")
    first_tf = _load_transfer_function(config.transfer_function_path(int(redshifts[0])))
    kh = first_tf[:, 0] * config.cosmology.h
    if np.any(np.diff(kh) <= 0):
        raise ValueError("transfer-function wavenumbers must be strictly increasing")
    tfkz = np.zeros((len(redshifts), first_tf.shape[0]))

    for i, redshift in enumerate(redshifts):
        tfz = _load_transfer_function(config.transfer_function_path(int(redshift)))
        if tfz.shape != first_tf.shape:
            raise ValueError("input transfer functions have different lengths")
        tfkz[i, :] = tfz[:, i_flag]

    return RectBivariateSpline(redshifts, kh, tfkz, kx=3, ky=3, s=0)


def interpolate_recfast(config: RuntimeConfig):
    """Interpolate RECFAST temperature and electron-fraction data in redshift."""
    if not config.rf_base.is_file():
        raise FileNotFoundError(
            f"configured RECFAST file does not exist: {config.rf_base}"
        )

    values = np.loadtxt(config.rf_base)
    # Reverse order so redshift is increasing.
    values = values[::-1]

    z = values[:, 0]
    xe = values[:, 1]
    temperature = values[:, -1]  # RECFAST++ can include extra columns

    temperature_spline = interp1d(z, temperature, kind="cubic")
    xe_spline = interp1d(z, xe, kind="cubic")
    return temperature_spline, xe_spline

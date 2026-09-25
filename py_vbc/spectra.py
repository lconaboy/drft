import numpy as np

from py_vbc.config import RuntimeConfig
from py_vbc.interpolations import interpolate_tf


def sigma(k, tf_spline, cosmology, R=None):
    """Return the mass-fluctuation integrand up to a normalization constant.

    :param k: wavenumbers in Mpc^-1
    :param tf_spline: interpolated total-matter transfer function
    :param cosmology: validated cosmological parameters
    :param R: real-space top-hat radius; defaults to 8 h^-1 Mpc
    """
    if R is None:
        R = 8.0 / cosmology.h

    def window(x):
        return (3 / x**3) * (np.sin(x) - x * np.cos(x))

    x = k * R
    transfer = tf_spline(k)
    return k ** (2 + cosmology.ns) * transfer**2 * window(x) ** 2


def calc_norm(config: RuntimeConfig):
    """Calculate the normalization required to match ``sigma_8``."""
    from scipy.integrate import quad

    tf0_spline = interpolate_tf(flag="t", z=0, config=config)
    kmin = np.min(tf0_spline.x)
    kmax = np.max(tf0_spline.x)

    # Sigma is highly oscillatory above k ~ 5 Mpc^-1. Splitting the integral
    # improves convergence; 10 Mpc^-1 works for the bundled TFs.
    kmid = 10.0
    epsrel = 1.0e-6
    limit = int(1e6)

    sigma_8c1 = quad(
        lambda k: sigma(k, tf0_spline, config.cosmology),
        kmin,
        kmid,
        limit=limit,
        epsrel=epsrel,
    )[0]
    sigma_8c2 = quad(
        lambda k: sigma(k, tf0_spline, config.cosmology),
        kmid,
        kmax,
        limit=limit,
        epsrel=epsrel,
    )[0]

    sigma_8c = np.sqrt(sigma_8c1 + sigma_8c2)
    return config.cosmology.sigma_8 / sigma_8c


def calc_power_spec(k, g, zstart, config: RuntimeConfig):
    """Calculate CDM, baryon, and velocity power spectra at *zstart*."""
    tf_c_spline = interpolate_tf("c", zstart, config)
    tf_c = tf_c_spline(k)
    norm = calc_norm(config)

    p_c = (
        2 * np.pi**2 * norm**2 * g[:, 0] ** 2 * tf_c**2 * k**config.cosmology.ns
    )
    p_b = (
        2 * np.pi**2 * norm**2 * g[:, 2] ** 2 * tf_c**2 * k**config.cosmology.ns
    )

    # Velocity power spectra are not renormalised.
    p_vc = (
        2 * np.pi**2 * norm**2 * g[:, 5] ** 2 * tf_c**2 * k**config.cosmology.ns
    )
    p_vb = (
        2 * np.pi**2 * norm**2 * g[:, 6] ** 2 * tf_c**2 * k**config.cosmology.ns
    )

    return p_c, p_b, p_vc, p_vb


def calc_tf(k, g, zstart, config: RuntimeConfig):
    """Calculate CDM, baryon, and velocity transfer functions at *zstart*."""
    tf_c_spline = interpolate_tf("c", zstart, config)
    tf_c = tf_c_spline(k)

    t_c = g[:, 0] * tf_c
    t_b = g[:, 2] * tf_c
    t_vc = g[:, 5] * tf_c
    t_vb = g[:, 6] * tf_c
    return t_c, t_b, t_vc, t_vb


def calc_delta(k, power):
    """Return the dimensionless power spectrum ``Delta^2(k)``."""
    return power * k**3 / (2 * np.pi**2)

import numpy as np

from py_vbc.constants import *
from py_vbc.interpolations import interpolate_tf


def sigma(k, tf_spline, R=8.0/hconst):
    """Integrand to calculate the mass fluctuations in a sphere of radius
    R, up to some constant of proportionality C, using transfer
    functions. Uses the fact that

    sigma^2 = int dk/k Delta^2(k) w(k)
    
    and
    
    Delta^2(k) = C k^(3+ns) T^2(k)
    
    w(x) is the window function, defined as the Fourier transform of a
    real-space top hat function.

    :param k: should be in units of Mpc^-1
    :param Tk: value of total matter transfer function at k and at z=0
    :param R: radius of sphere to calculate fluctuations, for sigma_8 this is 
              8 h^-1 Mpc
    :returns: dsigma^2/C where C is a normalization constant
    :rtype: float

    """
    
    def w(x):  
        return (3/x**3)*(np.sin(x) - x*np.cos(x))

    x = k*R
    Tk = tf_spline(k)
    
    return k**(2+ns) * Tk**2 * w(x)**2


def calc_norm():
    """This calculates the value of the normalization constant, with
    respect to sigma_8. The idea is that we already calculated
    sigma_8/C from the transfer functions, so by dividing the
    (specified) value of sigma_8 (at z=0) by our calculated sigma_8/C
    we get sqrt(C), which we can use to go from transfer functions to
    power spectra.

    :returns: normalisation constant
    :rtype: float

    """
    from scipy.integrate import quad
    
    tf0_spline = interpolate_tf(flag='t', z=0)

    # Need to check limits on the spline
    kmin = np.min(tf0_spline.x)
    kmax = np.max(tf0_spline.x)
    
    # Sigma is highly oscillatory above k ~ 5 Mpc^-1, so best to split
    # the integral into two parts to improve the convergence --
    # ideally kmid would be dynamically defined but 10 Mpc^-1 seems to
    # work
    kmid = 10.0

    # Arguments for quad
    epsrel = 1.0e-6
    limit = int(1e6)

    sigma_8c1 = quad(lambda k: sigma(k, tf0_spline), kmin, kmid, limit=limit, epsrel=epsrel)[0]
    sigma_8c2 = quad(lambda k: sigma(k, tf0_spline), kmid, kmax, limit=limit, epsrel=epsrel)[0]
    
    sigma_8c = np.sqrt(sigma_8c1 + sigma_8c2)
    
    return sigma_8/sigma_8c


def calc_power_spec(k, g, zstart):
    """Calculates the power spectra at z=zinit. First evolves the z=1000
    transfer functions forward using the linear growth factors
    calculated earlier, then converts to power spectra by using the
    normalization constant and

    P(k) propto T(k)^2 k^ns

    where ns is the tilt of the power spectrum.

    :param k: (array) k-values for which the growth factors were calculated
    :param g: (array) growth factors as produced by calc_derivs(), where the 
                      first column is for CDM perturbations and the third
                      column is the baryon perturbations 
    :returns: the CDM and baryon power spectra
    :rtype: arrays

    """
    # Density transfer function
    tf_c_spline = interpolate_tf('c', zstart)
    tf_c = tf_c_spline(k)

    norm = calc_norm()

    # Density power spectra
    p_c = 2*np.pi**2 * norm**2 * g[:, 0]**2 * tf_c**2 * k**ns
    p_b = 2*np.pi**2 * norm**2 * g[:, 2]**2 * tf_c**2 * k**ns

    # Velocity power spectra, not normalising
    p_vc = 2*np.pi**2 * norm**2 * g[:, 5]**2 * tf_c**2 * k**ns
    p_vb = 2*np.pi**2 * norm**2 * g[:, 6]**2 * tf_c**2 * k**ns
    
    return p_c, p_b, p_vc, p_vb


def calc_tf(k, g, zstart):
    """Calculates the transfer functions at z=zstart. This is done by
    undoing the previous normalisation by the CDM transfer functions.
    
    :param k: (array) k-values for which the growth factors were calculated
    :param g: (array) growth factors as produced by calc_derivs(), where the 
                      first column is for CDM perturbations and the third column
                      is the baryon perturbations 
    :returns: the CDM and baryon transfer functions
    :rtype: arrays

    """
    # Density transfer function
    tf_c_spline = interpolate_tf('c', zstart)
    tf_c = tf_c_spline(k)

    # Density
    t_c = g[:, 0] * tf_c
    t_b = g[:, 2] * tf_c

    # Velocity
    t_vc = g[:, 5] * tf_c
    t_vb = g[:, 6] * tf_c

    return t_c, t_b, t_vc, t_vb


def calc_delta(k, p):
    """Calculates the dimensionless power spectrum Delta^2 given a power
spectrum P(k)

    :param k: k values of P(k) 
    :param p: power spectrum
    :returns: Delta^2(k)
    :rtype: array

    """
    return p*k**3/(2*np.pi**2)

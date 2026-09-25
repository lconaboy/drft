import numpy as np
from scipy.integrate import solve_ivp

from py_vbc.config import RuntimeConfig
from py_vbc.constants import (
    BOLTZK,
    DELB,
    DELC,
    DELT,
    DELTAB,
    DELTAC,
    DELTAT,
    IMAG,
    MPCTOCM,
    MPCTOKM,
    MPROTON,
    MUB,
    REAL,
    TCMB,
    T_GAMMA,
    VEL,
    VELIMAG,
    VELREAL,
    YHE,
)
from py_vbc.interpolations import interpolate_recfast, interpolate_tf
from py_vbc.utils import hubble

"""
TODO
- incorporate photon density fluctuations
- check whether the transfer function derivatives have the right sign
"""


def set_ics(k, zstart, dz, config: RuntimeConfig):
    cosmology = config.cosmology

    # Compute splines
    tf_c_spline, dtf_c_spline = interpolate_tf("c", zstart, config, dz=dz)
    tf_b_spline, dtf_b_spline = interpolate_tf("b", zstart, config, dz=dz)
    tf_g_spline, dtf_g_spline = interpolate_tf("g", zstart, config, dz=dz)

    temperature_spline, xe_spline = interpolate_recfast(config)

    # Go from splines to actual values
    tf_c = tf_c_spline(k)
    tf_b = tf_b_spline(k)
    tf_g = tf_g_spline(k)
    dtf_c = dtf_c_spline(k)
    dtf_b = dtf_b_spline(k)
    dtf_g = dtf_g_spline(k)
    temperature = temperature_spline(zstart)
    xe = xe_spline(zstart)

    # Calculate value of Hubble parameter at z=zstart
    H = 100.0 * cosmology.h * hubble(zstart, cosmology) / MPCTOKM

    # Calculate T_CMB(zstart)
    Tcmb_z = TCMB * (1 + zstart)

    # ICs for the baryon overdensities and their derivatives

    # Pretty sure the derivatives should include a minus sign to go
    # from dy/dz to dy/dt...

    # ICs as in CICsASS
    delta_c = 1.0
    delta_b = (tf_b / tf_c) * delta_c
    delta_c_dot = (dtf_c / tf_c) * H * (1 + zstart) * delta_c
    delta_b_dot = (dtf_b / tf_b) * H * (1 + zstart) * delta_b
    # delta_c_dot = (tf_vc / tf_c) * H * k / (1+zstart)
    # delta_b_dot = (tf_vb / tf_b) * H * k / (1+zstart)

    # delta_c = tf_c * k * k
    # delta_b = tf_b * k * k
    # delta_c_dot = tf_vc *k *k # * H / (1+zstart)
    # delta_b_dot = tf_vb *k* k# H / (1+zstart)

    # Assume that initially temperature perturbations are coupled to
    # radiation perturbations, this is not strictly a good guess as we
    # are looking at zstart < 1000
    delta_g = (tf_g / tf_c) * delta_c
    delta_g_dot = (dtf_g / tf_g) * H * (1 + zstart) * delta_g

    astart = 1 / (1 + zstart)
    Tr = temperature / Tcmb_z
    delta_t = Tr * (
        delta_g * (1.25 * Tr - 1)
        - (T_GAMMA * astart**4 / xe) * (0.25 * delta_g_dot - 2.0 * delta_b_dot / 3.0)
    )  # d(delta_T - delta_Tg) = 0

    # LC isothermal for plots
    # print('-- warning temperatures are isothermal')
    # delta_t = 0

    # Store IC values in array y
    y = np.zeros(shape=(k.shape[0], 10))

    # Assign ICs
    y[:, DELTAC + REAL] = delta_c  # dark matter overdensity (REAL)
    y[:, DELTAC + IMAG] = 0.0  # dark matter overdensity (IMAG)
    y[:, DELTAC + VELREAL] = -delta_c_dot  # dark matter velocity divergence (REAL)
    y[:, DELTAC + VELIMAG] = 0.0  # dark matter velocity divergence (IMAG)
    y[:, DELTAB + REAL] = delta_b  # baryon overdensity (REAL)
    y[:, DELTAB + IMAG] = 0.0  # baryon overdensity (IMAG)
    y[:, DELTAB + VELREAL] = -delta_b_dot  # baryon velocity divergence (REAL)
    y[:, DELTAB + VELIMAG] = 0.0  # baryon velocity divergence (IMAG)
    y[:, DELTAT + REAL] = delta_t  # temperature fluctuations (REAL)
    y[:, DELTAT + IMAG] = 0.0  # temperature fluctuations (IMAG)

    return y


def set_ics_cicsass(k, zstart, dz, config: RuntimeConfig):
    """Set initial conditions using the legacy CICsASS convention."""
    cosmology = config.cosmology
    temperature_spline, xe_spline = interpolate_recfast(config)

    tf_c_spline, dtf_c_spline = interpolate_tf("c", zstart, config, dz=dz)
    tf_b_spline, dtf_b_spline = interpolate_tf("b", zstart, config, dz=dz)
    interpolate_tf("g", zstart, config, dz=dz)

    tf_c = tf_c_spline(k)
    tf_b = tf_b_spline(k)
    dtf_c = dtf_c_spline(k)
    dtf_b = dtf_b_spline(k)

    # Keep the legacy initial-condition convention, but use the central
    # constants and the runtime cosmology.
    del temperature_spline, xe_spline
    H = 100.0 * cosmology.h * hubble(zstart, cosmology) / MPCTOKM
    delta0 = 1.0
    delta_b = tf_b / tf_c * delta0
    delta_c_dot = dtf_c / tf_c * H * (1 + zstart) * delta0
    delta_b_dot = dtf_b / tf_b * delta_b * H * (1 + zstart)

    y = np.zeros(shape=(k.shape[0], 10))
    y[:, DELTAC + REAL] = delta0
    y[:, DELTAC + IMAG] = 0.0
    y[:, DELTAC + VELREAL] = -delta_c_dot
    y[:, DELTAC + VELIMAG] = 0.0
    y[:, DELTAB + REAL] = delta_b
    y[:, DELTAB + IMAG] = 0.0
    y[:, DELTAB + VELREAL] = -delta_b_dot
    y[:, DELTAB + VELIMAG] = 0.0
    y[:, DELTAT + REAL] = 0.0
    y[:, DELTAT + IMAG] = 0.0
    return y


def set_ics_comb(k, zstart, dz, config: RuntimeConfig):
    cosmology = config.cosmology

    # Compute splines
    tf_c_spline, dtf_c_spline = interpolate_tf("c", zstart, config, dz=dz)
    tf_b_spline, dtf_b_spline = interpolate_tf("b", zstart, config, dz=dz)
    tf_g_spline, dtf_g_spline = interpolate_tf("g", zstart, config, dz=dz)

    temperature_spline, xe_spline = interpolate_recfast(config)

    # Go from splines to actual values
    tf_c = tf_c_spline(k)
    tf_b = tf_b_spline(k)
    tf_g = tf_g_spline(k)
    # tf_vc = tf_vc_spline(k)
    # tf_vb = tf_vb_spline(k)
    dtf_c = dtf_c_spline(k)
    dtf_b = dtf_b_spline(k)
    dtf_g = dtf_g_spline(k)
    temperature = temperature_spline(zstart)
    xe = xe_spline(zstart)

    # Calculate value of Hubble parameter at z=zstart
    H = 100.0 * cosmology.h * hubble(zstart, cosmology) / MPCTOKM

    # Calculate T_CMB(zstart)
    Tcmb_z = TCMB * (1 + zstart)

    # ICs for the baryon overdensities and their derivatives

    # Pretty sure the derivatives should include a minus sign to go
    # from dy/dz to dy/dt...

    # ICs as in CICsASS
    delta_c = 1.0
    delta_b = (tf_b / tf_c) * delta_c
    delta_c_dot = (dtf_c / tf_c) * H * (1 + zstart) * delta_c
    delta_b_dot = (dtf_b / tf_b) * H * (1 + zstart) * delta_b

    # Assume that initially temperature perturbations are coupled to
    # radiation perturbations, this is not strictly a good guess as we
    # are looking at zstart < 1000
    delta_g = (tf_g / tf_c) * delta_c
    delta_g_dot = (dtf_g / tf_g) * H * (1 + zstart) * delta_g

    astart = 1 / (1 + zstart)
    Tr = temperature / Tcmb_z
    delta_t = Tr * (
        delta_g * (1.25 * Tr - 1)
        - (T_GAMMA * astart**4 / xe) * (0.25 * delta_g_dot - 2.0 * delta_b_dot / 3.0)
    )  # d(delta_T - delta_Tg) = 0

    # LC isothermal for plots
    # print('-- warning temperatures are isothermal')
    # delta_t = 0

    # Store IC values in array y
    y = np.zeros(shape=(k.shape[0], 5), dtype=complex)

    # Assign ICs
    y[:, DELC] = delta_c  # dark matter overdensity (complex)
    y[:, DELC + VEL] = -delta_c_dot  # dark matter velocity divergence (complex)
    y[:, DELB] = delta_b  # baryon overdensity (complex)
    y[:, DELB + VEL] = -delta_b_dot  # baryon velocity divergence (complex)
    y[:, DELT] = delta_t  # temperature fluctuations (complex)

    return y


def derivs_isothermal(z, y, k, T_spline, xe_spline, vbc, zstart, config: RuntimeConfig):
    cosmology = config.cosmology

    # Calculate Hubble parameter
    H = 100.0 * cosmology.h * hubble(z, cosmology) / MPCTOKM

    # Put redshift dependence into density parameters
    z1 = 1 + z
    o_m = cosmology.omega_m * z1**3
    f_b = cosmology.omega_b / cosmology.omega_m
    f_c = 1.0 - f_b
    o_r = cosmology.omega_r * z1**4
    o_mz = o_m / (o_m + o_r)
    o_b = f_b * o_mz
    o_c = f_c * o_mz

    # o_c = (omega_m-omega_b)/(omega_m + omega_r*z1)
    # o_b = omega_b/(omega_m + omega_r*z1)

    # Convert vstream to Mpc to it is in the same units as k and
    # insert decaying dependence with redshift
    vbck = vbc / ((1 + zstart) / z1) / MPCTOKM

    # Get temperature and electron fraction values
    T = T_spline(z)
    xe = xe_spline(z)

    # Get photon density fluctuations over CDM fluctuations (this is
    # how the ICs are defined)
    # g = g_spline(k, z)
    # c = c_spline(k, z)
    # gc = g / c

    fhe = 0.25 * YHE / ((1.0 - YHE) + 0.25 * YHE)

    # Convenience variables
    mu = vbck * k * cosmology.costh * z1
    alpha = 1.5 * H**2 * o_c
    beta = 1.5 * H**2 * o_b
    tau = (BOLTZK * T / (MUB * MPROTON * MPCTOCM**2)) * k**2 * z1**2
    # gamma = (xe/T_GAMMA)*(z1**5)
    gamma = (xe / T_GAMMA) * (z1**4)  # LC was the a-dependence overestimated?
    Tcmb_z = TCMB * z1
    Tr = Tcmb_z / T
    eta = 1.0 + fhe + xe

    # dy contains equations for delta_c_dot, theta_c_dot, delta_b_dot,
    # theta_b_dot, delta_t_dot split into REAL and imaginary parts
    dy = np.zeros(10)

    dy[DELTAC + REAL] = -y[DELTAC + VELREAL]
    dy[DELTAC + IMAG] = -y[DELTAC + VELIMAG]
    dy[DELTAC + VELREAL] = (
        -(alpha * y[DELTAC + REAL] + beta * y[DELTAB + REAL])
        - 2.0 * H * y[DELTAC + VELREAL]
    )
    dy[DELTAC + VELIMAG] = (
        -(alpha * y[DELTAC + IMAG] + beta * y[DELTAB + IMAG])
        - 2.0 * H * y[DELTAC + VELIMAG]
    )
    dy[DELTAB + REAL] = mu * y[DELTAB + IMAG] - y[DELTAB + VELREAL]
    dy[DELTAB + IMAG] = -mu * y[DELTAB + REAL] - y[DELTAB + VELIMAG]
    dy[DELTAB + VELREAL] = (
        mu * y[DELTAB + VELIMAG]
        - 2 * H * y[DELTAB + VELREAL]
        - (alpha * y[DELTAC + REAL] + beta * y[DELTAB + REAL])
        + tau * y[DELTAB + REAL]
    )
    dy[DELTAB + VELIMAG] = (
        -mu * y[DELTAB + VELREAL]
        - 2 * H * y[DELTAB + VELIMAG]
        - (alpha * y[DELTAC + IMAG] + beta * y[DELTAB + IMAG])
        + tau * y[DELTAB + IMAG]
    )

    dy /= -H * z1

    return dy


def derivs(z, y, k, T_spline, xe_spline, vbc, zstart, config: RuntimeConfig):
    cosmology = config.cosmology

    # Calculate Hubble parameter
    H = 100.0 * cosmology.h * hubble(z, cosmology) / MPCTOKM

    # Put redshift dependence into density parameters
    z1 = 1 + z
    o_m = cosmology.omega_m * z1**3
    f_b = cosmology.omega_b / cosmology.omega_m
    f_c = 1.0 - f_b
    o_r = cosmology.omega_r * z1**4
    o_mz = o_m / (o_m + o_r)
    o_b = f_b * o_mz
    o_c = f_c * o_mz

    # o_c = (omega_m-omega_b)/(omega_m + omega_r*z1)
    # o_b = omega_b/(omega_m + omega_r*z1)

    # Convert vstream to Mpc to it is in the same units as k and
    # insert decaying dependence with redshift
    vbck = vbc / ((1 + zstart) / z1) / MPCTOKM

    # Get temperature and electron fraction values
    T = T_spline(z)
    xe = xe_spline(z)

    # Get photon density fluctuations over CDM fluctuations (this is
    # how the ICs are defined)
    # g = g_spline(k, z)
    # c = c_spline(k, z)
    # gc = g / c

    fhe = 0.25 * YHE / ((1.0 - YHE) + 0.25 * YHE)

    # Convenience variables
    mu = vbck * k * cosmology.costh * z1
    alpha = 1.5 * H**2 * o_c
    beta = 1.5 * H**2 * o_b
    tau = (BOLTZK * T / (MUB * MPROTON * MPCTOCM**2)) * k**2 * z1**2
    # gamma = (xe/T_GAMMA)*(z1**5)
    gamma = (xe / T_GAMMA) * (z1**4)  # LC was the a-dependence overestimated?
    Tcmb_z = TCMB * z1
    Tr = Tcmb_z / T
    eta = 1.0 + fhe + xe

    # dy contains equations for delta_c_dot, theta_c_dot, delta_b_dot,
    # theta_b_dot, delta_t_dot split into REAL and imaginary parts
    dy = np.zeros(10)

    dy[DELTAC + REAL] = -y[DELTAC + VELREAL]
    dy[DELTAC + IMAG] = -y[DELTAC + VELIMAG]
    dy[DELTAC + VELREAL] = (
        -(alpha * y[DELTAC + REAL] + beta * y[DELTAB + REAL])
        - 2.0 * H * y[DELTAC + VELREAL]
    )
    dy[DELTAC + VELIMAG] = (
        -(alpha * y[DELTAC + IMAG] + beta * y[DELTAB + IMAG])
        - 2.0 * H * y[DELTAC + VELIMAG]
    )
    dy[DELTAB + REAL] = mu * y[DELTAB + IMAG] - y[DELTAB + VELREAL]
    dy[DELTAB + IMAG] = -mu * y[DELTAB + REAL] - y[DELTAB + VELIMAG]
    dy[DELTAB + VELREAL] = (
        mu * y[DELTAB + VELIMAG]
        - 2 * H * y[DELTAB + VELREAL]
        - (alpha * y[DELTAC + REAL] + beta * y[DELTAB + REAL])
        + tau * (y[DELTAB + REAL] + y[DELTAT + REAL])
    )
    dy[DELTAB + VELIMAG] = (
        -mu * y[DELTAB + VELREAL]
        - 2 * H * y[DELTAB + VELIMAG]
        - (alpha * y[DELTAC + IMAG] + beta * y[DELTAB + IMAG])
        + tau * (y[DELTAB + IMAG] + y[DELTAT + IMAG])
    )

    # Writing the temperature fluctuations exactly as in CA
    # dy[DELTAT+REAL] = mu*y[DELTAT+IMAG] - (2.0/3.0)*y[DELTAB+VELREAL] - gamma*y[DELTAT+REAL]/eta
    # dy[DELTAT+REAL] = mu*y[DELTAT+IMAG] - (2.0/3.0)*(dy[DELTAB+REAL] - mu*y[DELTAB+IMAG]) - gamma*y[DELTAT+REAL]
    # dy[DELTAT+IMAG] = -mu*y[DELTAT+REAL] - (2.0/3.0)*y[DELTAB+VELIMAG] - gamma*y[DELTAT+IMAG]/eta
    # dy[DELTAT+IMAG] = -mu*y[DELTAT+REAL] - (2.0/3.0)*(dy[DELTAB+IMAG] - mu*y[DELTAB+REAL]) - gamma*y[DELTAT+IMAG]

    # Naoz+ (2005) say this is only valid for z<=200, but Ahn showed
    # that actually its fine for z < 1000
    dy[DELTAT + REAL] = (2.0 / 3.0) * dy[DELTAB + REAL] - gamma * y[DELTAT + REAL]
    dy[DELTAT + IMAG] = (2.0 / 3.0) * dy[DELTAB + IMAG] - gamma * y[DELTAT + IMAG]

    # Proper set for z > 200 (i.e. Naoz+ (2005))
    # dy[DELTAT+REAL] = (2.0/3.0)*dy[DELTAB+REAL] + gamma*(gc * (1.25*Tr - 1) -
    #                                                      Tr*y[DELTAT+REAL])
    # dy[DELTAT+IMAG] = (2.0/3.0)*dy[DELTAB+IMAG] + gamma*(gc * (1.25*Tr - 1) -
    #                                                      Tr*y[DELTAT+IMAG])

    # dy[DELTAT+REAL] = (2.0/3.0)*dy[DELTAB+REAL] + gamma*(g * (1.25*Tr - 1) -
    #                                                      Tr*y[DELTAT+REAL])
    # dy[DELTAT+IMAG] = (2.0/3.0)*dy[DELTAB+IMAG] + gamma*(g * (1.25*Tr - 1) -
    #                                                      Tr*y[DELTAT+IMAG])

    # Convert dy from time to redshift derivative
    dy /= -H * z1

    return dy


def derivs_cicsass(a, y, k, T_spline, xe_spline, vbc, zstart, config: RuntimeConfig):
    """Evaluate the evolution equations using the legacy CICsASS convention."""
    cosmology = config.cosmology
    z = 1.0 / a - 1.0

    omega_r0 = cosmology.omega_r
    omega_b0 = cosmology.omega_b
    omega_c0 = cosmology.omega_m - cosmology.omega_b
    temperature = T_spline(z)
    xe = xe_spline(z)
    sound_speed = np.sqrt(BOLTZK * temperature / MUB / MPROTON) / MPCTOCM
    helium_fraction = 0.25 * YHE / ((1.0 - YHE) + 0.25 * YHE)
    H = 100.0 * cosmology.h * hubbleZ(z, cosmology) / MPCTOKM
    v_stream_kpc = vbc / MPCTOKM / (a * (1.0 + zstart))
    omega_c = omega_c0 / (omega_c0 + omega_b0 + omega_r0 / a)
    omega_b = omega_b0 / (omega_c0 + omega_b0 + omega_r0 / a)
    theta = cosmology.costh

    deriv = np.zeros(10)
    deriv[DELTAC + REAL] = -y[DELTAC + VELREAL]
    deriv[DELTAC + IMAG] = -y[DELTAC + VELIMAG]
    deriv[DELTAC + VELREAL] = (
        -1.5 * H**2 * (omega_c * y[DELTAC + REAL] + omega_b * y[DELTAB + REAL])
        - 2.0 * H * y[DELTAC + VELREAL]
    )
    deriv[DELTAC + VELIMAG] = (
        -1.5 * H**2 * (omega_c * y[DELTAC + IMAG] + omega_b * y[DELTAB + IMAG])
        - 2.0 * H * y[DELTAC + VELIMAG]
    )

    deriv[DELTAB + REAL] = (
        v_stream_kpc * k * theta * y[DELTAB + IMAG] / a - y[DELTAB + VELREAL]
    )
    deriv[DELTAB + IMAG] = (
        -v_stream_kpc * k * theta * y[DELTAB + REAL] / a - y[DELTAB + VELIMAG]
    )
    deriv[DELTAB + VELREAL] = (
        v_stream_kpc * k * theta * y[DELTAB + VELIMAG] / a
        - 1.5 * H**2 * (omega_c * y[DELTAC + REAL] + omega_b * y[DELTAB + REAL])
        - 2.0 * H * y[DELTAB + VELREAL]
        + sound_speed**2 * k**2 / a**2 * (y[DELTAB + REAL] + y[DELTAT + REAL])
    )
    deriv[DELTAB + VELIMAG] = (
        -v_stream_kpc * k * theta * y[DELTAB + VELREAL] / a
        - 1.5 * H**2 * (omega_c * y[DELTAC + IMAG] + omega_b * y[DELTAB + IMAG])
        - 2.0 * H * y[DELTAB + VELIMAG]
        + sound_speed**2 * k**2 / a**2 * (y[DELTAB + IMAG] + y[DELTAT + IMAG])
    )

    eta = 1.0 + helium_fraction + xe
    deriv[DELTAT + REAL] = (
        v_stream_kpc * k * theta * y[DELTAT + IMAG] / a
        + (2.0 / 3.0)
        * (deriv[DELTAB + REAL] - v_stream_kpc * k * theta * y[DELTAB + IMAG] / a)
        - xe / (eta * T_GAMMA * a**4) * (TCMB / a / temperature) * y[DELTAT + REAL]
    )
    deriv[DELTAT + IMAG] = (
        -v_stream_kpc * k * theta * y[DELTAT + REAL] / a
        + (2.0 / 3.0)
        * (deriv[DELTAB + IMAG] + v_stream_kpc * k * theta * y[DELTAB + REAL] / a)
        - xe / (eta * T_GAMMA * a**4) * (TCMB / a / temperature) * y[DELTAT + IMAG]
    )

    deriv /= H * a
    return deriv


def hubbleZ(zCurrent, cosmology):
    """Legacy name for the dimensionless Hubble parameter."""
    one_plus_z = 1.0 + zCurrent
    return np.sqrt(
        (1.0 - cosmology.omega_m - cosmology.omega_r)
        + cosmology.omega_m * one_plus_z**3
        + cosmology.omega_r * one_plus_z**4
    )


def derivs_comb(z, y, k, T_spline, xe_spline, vbc, zstart, config: RuntimeConfig):
    cosmology = config.cosmology

    # Calculate Hubble parameter
    H = 100.0 * cosmology.h * hubble(z, cosmology) / MPCTOKM

    z1 = 1 + z
    o_m = cosmology.omega_m * z1**3
    f_b = cosmology.omega_b / cosmology.omega_m
    f_c = 1.0 - f_b
    o_r = cosmology.omega_r * z1**4
    o_mz = o_m / (o_m + o_r)
    o_b = f_b * o_mz
    o_c = f_c * o_mz

    Tcmb_z = TCMB * z1

    # Convert vstream to Mpc to it is in the same units as k and
    # insert decaying dependence with redshift
    vbck = (z1 / (1 + zstart)) * vbc / MPCTOKM

    # Get temperature and electron fraction values
    T = T_spline(z)
    xe = xe_spline(z)

    # Get photon density fluctuations over CDM fluctuations (this is
    # how the ICs are defined)
    # g = g_spline(k, z)
    # c = c_spline(k, z)
    # gc = g / c

    fhe = 0.25 * YHE / ((1.0 - YHE) + 0.25 * YHE)

    # Convenience variables
    mu = vbck * k * cosmology.costh * z1
    alpha = 1.5 * H**2 * o_c
    beta = 1.5 * H**2 * o_b
    tau = (BOLTZK * T / (MUB * MPROTON * MPCTOCM**2)) * k**2 * z1**2
    # gamma = (xe/T_GAMMA)*(z1**5)
    gamma = (xe / T_GAMMA) * (z1**4)  # LC was the a-dependence overestimated?
    Tr = Tcmb_z / T
    eta = 1.0 + fhe + xe

    # dy contains equations for delta_c_dot, theta_c_dot, delta_b_dot,
    # theta_b_dot, delta_t_dot split into REAL and imaginary parts
    dy = np.zeros(5, dtype=complex)

    # DELC = 0
    # DELB = 2
    # VEL = 1
    # DELT = 4

    dy[DELC] = -y[DELC + VEL]
    dy[DELC + VEL] = -(alpha * y[DELC] + beta * y[DELB]) - 2.0 * H * y[DELC + VEL]
    dy[DELB] = -1j * mu * y[DELB] - y[DELB + VEL]
    dy[DELB + VEL] = (
        -((alpha * y[DELC] + beta * y[DELB]))
        - (2 * H * y[DELB + VEL])
        - (1j * mu * y[DELB + VEL])
        + tau * (y[DELB] + y[DELT])
    )

    # Writing the temperature fluctuations exactly as in CA
    # dy[DELTAT+REAL] = mu*y[DELTAT+IMAG] - (2.0/3.0)*y[DELTAB+VELREAL] - gamma*y[DELTAT+REAL]/eta
    # dy[DELTAT+REAL] = mu*y[DELTAT+IMAG] - (2.0/3.0)*(dy[DELTAB+REAL] - mu*y[DELTAB+IMAG]) - gamma*y[DELTAT+REAL]
    # dy[DELTAT+IMAG] = -mu*y[DELTAT+REAL] - (2.0/3.0)*y[DELTAB+VELIMAG] - gamma*y[DELTAT+IMAG]/eta
    # dy[DELTAT+IMAG] = -mu*y[DELTAT+REAL] - (2.0/3.0)*(dy[DELTAB+IMAG] - mu*y[DELTAB+REAL]) - gamma*y[DELTAT+IMAG]

    # Naoz+ (2005) say this is only valid for z<=200, but Ahn showed
    # that actually its fine for z < 1000
    # dy[DELT] = (2.0/3.0)*dy[DELB] - gamma * Tr * y[DELT]
    dy[DELT] = (2.0 / 3.0) * (-1j * mu * y[DELB] - y[DELB + VEL]) - gamma * Tr * y[DELT]

    dy /= -H * z1

    return dy


def calc_derivs(
    k,
    vbc,
    zstart,
    zend,
    dz,
    config: RuntimeConfig,
    verbose=False,
    integrator="LSODA",
    isothermal=False,
    require_z0=False,
):
    if dz <= 0:
        raise ValueError("dz must be greater than zero")
    redshifts = [int(zstart), int(zstart - dz), int(zstart + dz)]
    if require_z0:
        redshifts.append(0)
    config.validate_files(redshifts)

    T_spline, xe_spline = interpolate_recfast(config)
    y0 = set_ics(k, zstart, dz, config)
    y = np.zeros(shape=y0.shape)
    nk = y0.shape[0]
    z = np.zeros(shape=(nk, 7))

    if verbose:
        print(
            "Solving evolution equations from z={0:4.2f} to z={1:4.2f} with v_bc={2:4.2f} ...".format(
                zstart, zend, vbc
            )
        )
    for i, ik in enumerate(k):
        # Progress, remove end='\r' for Python 2 compatibility
        if verbose:
            if i + 1 != nk:
                print(
                    "    k = {0:.3f} Mpc^-1 [{1:d}/{2:d}]".format(ik, i + 1, nk),
                    end="\r",
                )
            else:
                print("    k = {0:.3f} Mpc^-1 [{1:d}/{2:d}]".format(ik, i + 1, nk))

        if isothermal:
            derivs_ = derivs_isothermal
        else:
            derivs_ = derivs

        r = solve_ivp(
            fun=lambda z, y: derivs_(
                z, y, ik, T_spline, xe_spline, vbc, zstart, config
            ),
            t_span=(zstart, zend),
            t_eval=np.array([zend]),
            rtol=1.0e-6,
            y0=y0[i, :],
            method=integrator,
            jac=None,
        )

        y[i, :] = np.transpose(r.y)

    if verbose:
        print("    ... done.")

    # Convert from REAL and imaginary components to magnitude
    z[:, 0] = np.sqrt(y[:, 0] ** 2 + y[:, 1] ** 2)  # delta_c
    z[:, 1] = np.sqrt(y[:, 2] ** 2 + y[:, 3] ** 2)  # theta_c
    z[:, 2] = np.sqrt(y[:, 4] ** 2 + y[:, 5] ** 2)  # delta_b
    z[:, 3] = np.sqrt(y[:, 6] ** 2 + y[:, 7] ** 2)  # theta_b
    z[:, 4] = np.sqrt(y[:, 8] ** 2 + y[:, 9] ** 2)  # delta_t
    z[:, 5] = theta_to_vel(y[:, [2, 3]], k, zend)  # v_c
    z[:, 6] = theta_to_vel(y[:, [6, 7]], k, zend)  # v_b

    return z


def calc_derivs_cicsass(
    k,
    vbc,
    zstart,
    zend,
    dz,
    config: RuntimeConfig,
    verbose=False,
    integrator="LSODA",
):
    config.validate_files([int(zstart), int(zstart - dz), int(zstart + dz)])
    T_spline, xe_spline = interpolate_recfast(config)
    y0 = set_ics_cicsass(k, zstart, dz, config)
    y = np.zeros(shape=y0.shape)
    nk = y0.shape[0]
    z = np.zeros(shape=(nk, 7))

    if verbose:
        print(
            "Solving evolution equations from z={0:4.2f} to z={1:4.2f} with v_bc={2:4.2f} ...".format(
                zstart, zend, vbc
            )
        )
    for i, ik in enumerate(k):
        # Progress, remove end='\r' for Python 2 compatibility
        if verbose:
            if i + 1 != nk:
                print(
                    "    k = {0:.3f} Mpc^-1 [{1:d}/{2:d}]".format(ik, i + 1, nk),
                    end="\r",
                )
            else:
                print("    k = {0:.3f} Mpc^-1 [{1:d}/{2:d}]".format(ik, i + 1, nk))

        r = solve_ivp(
            fun=lambda a, y: derivs_cicsass(
                a, y, ik, T_spline, xe_spline, vbc, zstart, config
            ),
            t_span=(1.0 / (1.0 + zstart), 1.0 / (1.0 + zend)),
            t_eval=np.array([1.0 / (1.0 + zend)]),
            atol=1e-8,
            rtol=1e-18,
            first_step=1e-3,
            y0=y0[i, :],
            method=integrator,
            jac=None,
        )

        y[i, :] = np.transpose(r.y)

    if verbose:
        print("    ... done.")

    # Convert from REAL and imaginary components to magnitude
    z[:, 0] = np.sqrt(y[:, DELTAC + REAL] ** 2 + y[:, DELTAC + IMAG] ** 2)  # delta_c
    z[:, 1] = np.sqrt(
        y[:, DELTAC + VELREAL] ** 2 + y[:, DELTAC + VELIMAG] ** 2
    )  # theta_c
    z[:, 2] = np.sqrt(y[:, DELTAB + REAL] ** 2 + y[:, DELTAB + IMAG] ** 2)  # delta_b
    z[:, 3] = np.sqrt(
        y[:, DELTAB + VELREAL] ** 2 + y[:, DELTAB + VELIMAG] ** 2
    )  # theta_b
    z[:, 4] = np.sqrt(y[:, DELTAT + REAL] ** 2 + y[:, DELTAT + IMAG] ** 2)  # delta_t
    z[:, 5] = theta_to_vel(y[:, [DELTAC + VELREAL, DELTAC + VELIMAG]], k, zend)  # v_c
    z[:, 6] = theta_to_vel(y[:, [DELTAB + VELREAL, DELTAB + VELIMAG]], k, zend)  # v_b

    return z


def calc_derivs_comp(
    k,
    vbc,
    zstart,
    zend,
    dz,
    config: RuntimeConfig,
    verbose=False,
    integrator="LSODA",
):
    config.validate_files([int(zstart), int(zstart - dz), int(zstart + dz)])
    T_spline, xe_spline = interpolate_recfast(config)
    y0 = set_ics_comb(k, zstart, dz, config)
    y = np.zeros(shape=y0.shape)
    nk = y0.shape[0]
    z = np.zeros(shape=(nk, 7))

    if verbose:
        print(
            "Solving evolution equations from z={0:4.2f} to z={1:4.2f} with v_bc={2:4.2f} ...".format(
                zstart, zend, vbc
            )
        )
    for i, ik in enumerate(k):
        # Progress, remove end='\r' for Python 2 compatibility
        if verbose:
            if i + 1 != nk:
                print(
                    "    k = {0:.3f} Mpc^-1 [{1:d}/{2:d}]".format(ik, i + 1, nk),
                    end="\r",
                )
            else:
                print("    k = {0:.3f} Mpc^-1 [{1:d}/{2:d}]".format(ik, i + 1, nk))

        r = solve_ivp(
            fun=lambda z, y: derivs_comb(
                z, y, ik, T_spline, xe_spline, vbc, zstart, config
            ),
            t_span=(zstart, zend),
            t_eval=np.array([zend]),
            rtol=1.0e-6,
            y0=y0[i, :],
            method=integrator,
            jac=None,
        )

        y[i, :] = np.transpose(r.y)

    if verbose:
        print("    ... done.")

    # Convert from REAL and imaginary components to magnitude
    z[:, 0] = np.abs(y[:, 0])  # delta_c
    z[:, 1] = np.abs(y[:, 1])  # theta_c
    z[:, 2] = np.abs(y[:, 2])  # delta_b
    z[:, 3] = np.abs(y[:, 3])  # theta_b
    z[:, 4] = np.abs(y[:, 4])  # delta_t
    z[:, 5] = theta_to_vel_comp(y[:, 1], k, zend)  # v_c
    z[:, 6] = theta_to_vel_comp(y[:, 3], k, zend)  # v_b

    return z


def theta_to_vel(theta, k, z):
    """
    Converts velocity divergence (in k-space) to proper velocity by

    v = - i a \theta / k

    :param theta:
        (array) velocity divergence
    :param k:
        (array) values of k to calculate velocity at
    :param z:
        (float) redshift to calculate velocity at
    """

    def convert(theta, k, z):
        a = 1.0 / (1.0 + z)
        return (-1j * a * theta) / k

    real_theta = theta[:, 0]
    imag_theta = theta[:, 1]
    comp_theta = real_theta + 1j * imag_theta

    comp_vel = convert(comp_theta, k, z)
    real_vel = comp_vel.real
    imag_vel = comp_vel.imag

    return np.sqrt(real_vel**2 + imag_vel**2)


def theta_to_vel_comp(comp_theta, k, z):
    """
    Converts velocity divergence (in k-space) to proper velocity by

    v = - i a \theta / k

    :param theta:
        (array, complex) velocity divergence
    :param k:
        (array) values of k to calculate velocity at
    :param z:
        (float) redshift to calculate velocity at
    """

    def convert(theta, k, z):
        a = 1.0 / (1.0 + z)
        return (-1j * a * theta) / k

    comp_vel = convert(comp_theta, k, z)
    vel = np.abs(comp_vel)

    return vel

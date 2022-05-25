import numpy as np
# import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from py_vbc.constants import *
from py_vbc.utils import hubble
from py_vbc.interpolations import interpolate_tf, interpolate_tf2d, interpolate_recfast
    
"""
TODO
- incorporate photon density fluctuations
- check whether the transfer function derivatives have the right sign
"""

# TCMB = 2.726                   # // CMB temperature at z=0
# MPROTON = 1.6726e-24           # // in grams
# MELECTRON = 9.11e-28
# MEL_EV = 5.11e5
# LIGHTSPEED = 3.0e10
# CRITDENSITY = 1.8791e-29   # // Current critical density, 
#                                         # // in g/cm^3, not including h^2!
# SIGMAT = 0.665e-24   # //Thomson cross section in cm^2
# YHE = 0.25
# BOLTZK = 1.3806e-16            # // in erg/K
# MUB = 1.22           # // Mean molecular weight, for primordial
# aSB = 7.56e-15 # //u = aSB T^4 


def set_ics(k, zstart, dz):    
    # Compute splines
    tf_c_spline, dtf_c_spline = interpolate_tf('c', zstart, dz)
    tf_b_spline, dtf_b_spline = interpolate_tf('b', zstart, dz)
    tf_g_spline, dtf_g_spline = interpolate_tf('g', zstart, dz)
    tf_vc_spline = interpolate_tf('vc', zstart)
    tf_vb_spline = interpolate_tf('vb', zstart)
    
    T_spline, xe_spline = interpolate_recfast()

    # Go from splines to actual values
    tf_c = tf_c_spline(k)
    tf_b = tf_b_spline(k)
    tf_g = tf_g_spline(k)
    tf_vc = tf_vc_spline(k)
    tf_vb = tf_vb_spline(k)
    dtf_c = dtf_c_spline(k)
    dtf_b = dtf_b_spline(k)
    dtf_g = dtf_g_spline(k)
    T = T_spline(zstart)
    xe = xe_spline(zstart)

    # Calculate value of Hubble parameter at z=zstart
    H = (100.0*hconst)*hubble(zstart)/mpctokm

    # Calculate T_CMB(zstart)
    Tcmb_z = Tcmb*(1+zstart)

    # ICs for the baryon overdensities and their derivatives

    # Pretty sure the derivatives should include a minus sign to go
    # from dy/dz to dy/dt...

    # ICs as in CICsASS
    delta_c = 1.0
    delta_b = (tf_b/tf_c) * delta_c
    delta_c_dot = (dtf_c/tf_c)*H*(1+zstart)*delta_c
    delta_b_dot = (dtf_b/tf_b)*H*(1+zstart)*delta_b
    # delta_c_dot = (tf_vc / tf_c) * H * k / (1+zstart)
    # delta_b_dot = (tf_vb / tf_b) * H * k / (1+zstart)

    # delta_c = tf_c * k * k
    # delta_b = tf_b * k * k
    # delta_c_dot = tf_vc *k *k # * H / (1+zstart)
    # delta_b_dot = tf_vb *k* k# H / (1+zstart)
    
    
    # Assume that initially temperature perturbations are coupled to
    # radiation perturbations, this is not strictly a good guess as we
    # are looking at zstart < 1000
    delta_g = (tf_g/tf_c)*delta_c
    delta_g_dot = (dtf_g/tf_g)*H*(1+zstart)*delta_g

    
    astart= 1 / (1 + zstart)
    Tr = T / Tcmb_z
    delta_t = Tr * (
        delta_g * (1.25*Tr - 1) -
        (t_gamma * astart**4 / xe) * 
        (0.25*delta_g_dot - 2.*delta_b_dot/3.))  # d(delta_T - delta_Tg) = 0

    # LC isothermal for plots
    # print('-- warning temperatures are isothermal')
    # delta_t = 0
    
    # Store IC values in array y
    y = np.zeros(shape=(k.shape[0], 10))

    # Assign ICs
    y[:, deltac+real] = delta_c          # dark matter overdensity (real)
    y[:, deltac+imag] = 0.0              # dark matter overdensity (imag)
    y[:, deltac+velreal] = -delta_c_dot  # dark matter velocity divergence (real)
    y[:, deltac+velimag] = 0.0           # dark matter velocity divergence (imag)
    y[:, deltab+real] = delta_b          # baryon overdensity (real)
    y[:, deltab+imag] = 0.0              # baryon overdensity (imag)
    y[:, deltab+velreal] = -delta_b_dot  # baryon velocity divergence (real)
    y[:, deltab+velimag] = 0.0           # baryon velocity divergence (imag)
    y[:, deltat+real] = delta_t          # temperature fluctuations (real)
    y[:, deltat+imag] = 0.0              # temperature fluctuations (imag)

    return y


def set_ics_cicsass(k, zstart, dz):
    MPCTOKM = 3.08568e19
    MPCTOCM = 3.08568e24
    DELTAC = 0
    DELTAB = 4
    REAL = 0
    IMAG = 1
    VELREAL = 2
    VELIMAG = 3
    DELTAT = DELTAB+VELIMAG+1
    
    T_spline, xe_spline = interpolate_recfast()
    T = T_spline(zstart)
    xe = xe_spline(zstart)
    
    OR0 = omega_r
    OB0 = omega_b
    OC0 = omega_m - omega_b
    ZSTART = zstart
    a = 1. / (1. + zstart)
    cs = np.sqrt(BOLTZK*T/MUB/MPROTON)/MPCTOCM #; //sound speed for iosthermal gas
    YHE = 0.25
    fHe = .25*YHE/((1. - YHE) + .25*YHE);
    H = 100.*hconst*hubbleZ(1./a - 1.) /MPCTOKM;

    OMEGAC = OC0/(OC0 +OB0 + OR0/a);
    OMEGAB = OB0/(OC0 +OB0 + OR0/a);
    
    # Compute splines
    tf_c_spline, dtf_c_spline = interpolate_tf('c', zstart, dz)
    tf_b_spline, dtf_b_spline = interpolate_tf('b', zstart, dz)
    tf_g_spline, dtf_g_spline = interpolate_tf('g', zstart, dz)
    # tf_vc_spline = interpolate_tf('vc', zstart)
    # tf_vb_spline = interpolate_tf('vb', zstart) 

    # Go from splines to actual values
    tf_c = tf_c_spline(k)
    tf_b = tf_b_spline(k)
    tf_g = tf_g_spline(k)
    # tf_vc = tf_vc_spline(k)
    # tf_vb = tf_vb_spline(k)
    dtf_c = dtf_c_spline(k)
    dtf_b = dtf_b_spline(k)
    dtf_g = dtf_g_spline(k)
    T = T_spline(zstart)
    xe = xe_spline(zstart)

    TFb = tf_b
    TFc = tf_c
    TFdirc = dtf_c
    TFdirb = dtf_b

    # ICs for the baryon overdensities and their derivatives

    # Pretty sure the derivatives should include a minus sign to go
    # from dy/dz to dy/dt...

    
    # Store IC values in array y
    y = np.zeros(shape=(k.shape[0], 10))

    DELTA0 = 1.; # //sets amplitude to unity


    H = 100.*hconst*hubbleZ(ZSTART)/MPCTOKM;

    deltab = TFb/TFc*DELTA0;
    deltacdot = TFdirc/TFc*H*(1+ZSTART)*DELTA0;
    deltabdot = TFdirb/TFb*deltab*H*(1+ZSTART);

    y[:, 0] = DELTA0; # //real part of dark matter density at specified k
    y[:, 1] = 0.;  # //imaginary part (start with zero for simplicity)
    y[:, 2] = -deltacdot; # //1st derivative of delta_dm, real part
    y[:, 3] = 0.; # //1st derivative of delta_dm, imaginary part
    y[:, 4] = deltab;  # //real part of baryonic overdensity
    y[:, 5] = 0.; # //imaginary part of baryonic overdensity
    y[:, 6] = -deltabdot; # //real part of baryonic velocity
    y[:, 7] = 0.;        # //imaginary part
    y[:, 8] = 0.; # //temperature
    # // initially assumed isothermal (which should be a good assumption on simulated scales -- a better assumption would be to set
    # //equal to temperature inhomogenety as the two are tightly coupled (as in Bovi and Dvorkin 2012) 
    y[:, 9] = 0.;


    return y


def set_ics_comb(k, zstart, dz):    
    # Compute splines
    tf_c_spline, dtf_c_spline = interpolate_tf('c', zstart, dz)
    tf_b_spline, dtf_b_spline = interpolate_tf('b', zstart, dz)
    tf_g_spline, dtf_g_spline = interpolate_tf('g', zstart, dz)
    # tf_vc_spline = interpolate_tf('vc', zstart)
    # tf_vb_spline = interpolate_tf('vb', zstart)
    
    T_spline, xe_spline = interpolate_recfast()

    # Go from splines to actual values
    tf_c = tf_c_spline(k)
    tf_b = tf_b_spline(k)
    tf_g = tf_g_spline(k)
    # tf_vc = tf_vc_spline(k)
    # tf_vb = tf_vb_spline(k)
    dtf_c = dtf_c_spline(k)
    dtf_b = dtf_b_spline(k)
    dtf_g = dtf_g_spline(k)
    T = T_spline(zstart)
    xe = xe_spline(zstart)

    # Calculate value of Hubble parameter at z=zstart
    H = (100.0*hconst)*hubble(zstart)/mpctokm

    # Calculate T_CMB(zstart)
    Tcmb_z = Tcmb*(1+zstart)

    # ICs for the baryon overdensities and their derivatives

    # Pretty sure the derivatives should include a minus sign to go
    # from dy/dz to dy/dt...

    # ICs as in CICsASS
    delta_c = 1.0
    delta_b = (tf_b/tf_c) * delta_c
    delta_c_dot = (dtf_c/tf_c)*H*(1+zstart)*delta_c
    delta_b_dot = (dtf_b/tf_b)*H*(1+zstart)*delta_b

    # Assume that initially temperature perturbations are coupled to
    # radiation perturbations, this is not strictly a good guess as we
    # are looking at zstart < 1000
    delta_g = (tf_g/tf_c)*delta_c
    delta_g_dot = (dtf_g/tf_g)*H*(1+zstart)*delta_g

    
    astart= 1 / (1 + zstart)
    Tr = T / Tcmb_z
    delta_t =  Tr * (
        delta_g * (1.25*Tr - 1) -
        (t_gamma * astart**4 / xe) * 
        (0.25*delta_g_dot - 2.*delta_b_dot/3.))  # d(delta_T - delta_Tg) = 0

    # LC isothermal for plots
    # print('-- warning temperatures are isothermal')
    # delta_t = 0
    
    # Store IC values in array y
    y = np.zeros(shape=(k.shape[0], 5), dtype=complex)

    # Assign ICs
    y[:, delc] = delta_c          # dark matter overdensity (complex)
    y[:, delc+vel] = -delta_c_dot  # dark matter velocity divergence (complex)
    y[:, delb] = delta_b          # baryon overdensity (complex)
    y[:, delb+vel] = -delta_b_dot  # baryon velocity divergence (complex)
    y[:, delt] = delta_t          # temperature fluctuations (complex)

    return y


# def set_ics_comb(k, zstart, dz):    
#     # Compute splines
#     tf_c_spline, dtf_c_spline = interpolate_tf('c')
#     tf_b_spline, dtf_b_spline = interpolate_tf('b')

#     # Go from splines to actual values
#     tf_c = tf_c_spline(k)
#     tf_b = tf_b_spline(k)
#     dtf_c = dtf_c_spline(k)
#     dtf_b = dtf_b_spline(k)

#     # Calculate value of Hubble parameter at z=zstart
#     H = (100.0*hconst)*hubble(zstart)/mpctokm

#     # ICs for the baryon overdensities and their derivatives
#     delta_c = 1.0
#     # delta_c = (tf_c/tf_b)*delta_b
#     delta_b = (tf_b/tf_c)*delta_c
#     delta_c_dot = (dtf_c/tf_c)*H*(1+zstart)*delta_c
#     delta_b_dot = (dtf_b/tf_b)*H*(1+zstart)*delta_b

#     # Store IC values in array y
#     y = np.zeros(shape=(k.shape[0], 5), dtype=complex)

#     # Assign ICs
#     y[:, delc] = delta_c           # dark matter overdensity
#     y[:, delc+vel] = -delta_c_dot  # time derivative of dark matter overdensity
#     y[:, delb] = delta_b           # baryon overdensity
#     y[:, delb+vel] = -delta_b_dot  # time derivative of baryon overdensity
#     y[:, delt] = 0.0               # temperature fluctuations

#     return y


# Original complex differential equations as in O'Leary and McQuinn
# (2012)
# def derivs_comb(z, y, k, T_spline, xe_spline):
#     from utils import hubble

#     # Calculate Hubble parameter
#     H = (100.0*hconst)*hubble(z)/mpctokm

#     # Put redshift dependence into density paramters
#     o_c = (omega_m-omega_b)/(omega_m + omega_r*(1+z))
#     o_b = omega_b/(omega_m + omega_r*(1+z))

#     # Convert vstream to Mpc to it is in the same units as k and
#     # insert decaying dependence with redshift
#     vbck = vstream/(((1+zstart)/(1+z))*mpctokm)

#     # # Get temperature and electron fraction values
#     # # T_spline, xe_spline = interpolate_recfast()
#     T = T_spline(z)
#     xe = xe_spline(z)

#     # Calculate (isothermal) sound speed for that redshift
#     cs = np.sqrt(boltzk*T/(mub*mproton))/mpctocm    
    
#     # dy contains the equations for delta_c_dot, theta_c_dot,
#     # delta_b_dot, theta_b_dot, delta_t_dot
#     dy = np.zeros(5, dtype=complex)
#     dy[delc] = -y[delc+vel]  # Eq. (A1)
#     dy[delc+vel] = -1.5*H**2*(o_c*y[delc] + o_b*y[delb]) - 2*H*y[delc+vel]  # Eq. (A2)
#     dy[delb] = -vbck*k*costh*y[delb]*(1+z)*1j - y[delb+vel]  # Eq. (A3)
#     dy[delb+vel] = (- vbck*k*costh*y[delc+vel]*(1+z)*1j - 1.5*H**2*(o_c*y[delc] + o_b*y[delb])
#                     - 2*H*y[delb+vel] + cs**2*k**2*(1+z)**2*(y[delb] + y[delt]))  # Eq. (A4)
#     dy[delt] = (-(2.0/3.0)*(vbck*k*costh*y[delb]*(1+z)*1j - y[delb+vel])
#                 -(xe/t_gamma)*(Tcmb/T)*(1+z)**5*y[delt])# - vbck*k*costh*y[delt]*(1+z)*1j)  # Ahn

#     # Convert dy from time to redshift derivative
#     dy /= -H*(1+z)
 
#     return dy


def derivs(z, y, k, T_spline, xe_spline, vbc, zstart):
    from py_vbc.utils import hubble
    
    # Calculate Hubble parameter
    H = (100.0*hconst)*hubble(z)/mpctokm
    
    # Put redshift dependence into density parameters
    z1 = 1+z
    o_m = omega_m * (z1 ** 3.)
    f_b = omega_b / omega_m
    f_c = 1. - f_b
    # o_b = omega_b * (z1 ** 3.)
    # o_c = o_m - o_b
    o_r = omega_r * (z1 ** 4.)
    o_mz = o_m / (o_m + o_r)
    o_b = f_b * o_mz
    o_c = f_c * o_mz
    
    # o_c = (omega_m-omega_b)/(omega_m + omega_r*z1)
    # o_b = omega_b/(omega_m + omega_r*z1)

    # Convert vstream to Mpc to it is in the same units as k and
    # insert decaying dependence with redshift
    vbck = vbc/((1+zstart)/z1)/mpctokm

    # Get temperature and electron fraction values
    T = T_spline(z)
    xe = xe_spline(z)


    # Get photon density fluctuations over CDM fluctuations (this is
    # how the ICs are defined)
    # g = g_spline(k, z)
    # c = c_spline(k, z)
    # gc = g / c
    
    yhe = 0.25
    fhe = 0.25*yhe/((1.0-yhe) + 0.25*yhe)
    
    # Convenience variables
    mu = vbck*k*costh*z1
    alpha = 1.5*H**2*o_c
    beta = 1.5*H**2*o_b
    tau = (boltzk*T/(mub*mproton*mpctocm**2)) * k**2 * z1**2 
    # gamma = (xe/t_gamma)*(z1**5)
    gamma = (xe/t_gamma)*(z1**4)  # LC was the a-dependence overestimated?
    Tcmb_z = Tcmb * z1
    Tr = (Tcmb_z/T)
    eta = 1.0 + fhe +xe
    
    # dy contains equations for delta_c_dot, theta_c_dot, delta_b_dot,
    # theta_b_dot, delta_t_dot split into real and imaginary parts
    dy = np.zeros(10)
    
    dy[deltac+real] = -y[deltac+velreal]
    dy[deltac+imag] = -y[deltac+velimag]
    dy[deltac+velreal] = (- (alpha*y[deltac+real] + beta*y[deltab+real])
                          - 2.0*H*y[deltac+velreal])
    dy[deltac+velimag] = (- (alpha*y[deltac+imag] + beta*y[deltab+imag])
                          - 2.0*H*y[deltac+velimag])
    dy[deltab+real] = mu*y[deltab+imag] - y[deltab+velreal]
    dy[deltab+imag] = -mu*y[deltab+real] - y[deltab+velimag]
    dy[deltab+velreal] = (mu*y[deltab+velimag] - 2*H*y[deltab+velreal]
                          - (alpha*y[deltac+real] + beta*y[deltab+real])
                          + tau*(y[deltab+real] + y[deltat+real]))
    dy[deltab+velimag] = (- mu*y[deltab+velreal] - 2*H*y[deltab+velimag]
                          - (alpha*y[deltac+imag] + beta*y[deltab+imag])
                          + tau*(y[deltab+imag] + y[deltat+imag]))

    # Writing the temperature fluctuations exactly as in CA
    # dy[deltat+real] = mu*y[deltat+imag] - (2.0/3.0)*y[deltab+velreal] - gamma*y[deltat+real]/eta
    # dy[deltat+real] = mu*y[deltat+imag] - (2.0/3.0)*(dy[deltab+real] - mu*y[deltab+imag]) - gamma*y[deltat+real]
    # dy[deltat+imag] = -mu*y[deltat+real] - (2.0/3.0)*y[deltab+velimag] - gamma*y[deltat+imag]/eta
    # dy[deltat+imag] = -mu*y[deltat+real] - (2.0/3.0)*(dy[deltab+imag] - mu*y[deltab+real]) - gamma*y[deltat+imag]

    # Naoz+ (2005) say this is only valid for z<=200, but Ahn showed
    # that actually its fine for z < 1000
    dy[deltat+real] = (2.0/3.0)*dy[deltab+real] - gamma*y[deltat+real]
    dy[deltat+imag] = (2.0/3.0)*dy[deltab+imag] - gamma*y[deltat+imag]

    # Proper set for z > 200 (i.e. Naoz+ (2005))
    # dy[deltat+real] = (2.0/3.0)*dy[deltab+real] + gamma*(gc * (1.25*Tr - 1) -
    #                                                      Tr*y[deltat+real])
    # dy[deltat+imag] = (2.0/3.0)*dy[deltab+imag] + gamma*(gc * (1.25*Tr - 1) -
    #                                                      Tr*y[deltat+imag])

    # dy[deltat+real] = (2.0/3.0)*dy[deltab+real] + gamma*(g * (1.25*Tr - 1) -
    #                                                      Tr*y[deltat+real])
    # dy[deltat+imag] = (2.0/3.0)*dy[deltab+imag] + gamma*(g * (1.25*Tr - 1) -
    #                                                      Tr*y[deltat+imag])

    # Convert dy from time to redshift derivative
    dy /= (-H*z1)
    
    return dy


def derivs_cicsass(a, y, k, T_spline, xe_spline, vbc, zstart):
    MPCTOKM = 3.08568e19
    MPCTOCM = 3.08568e24
    DELTAC = 0
    DELTAB = 4
    REAL = 0
    IMAG = 1
    VELREAL = 2
    VELIMAG = 3
    DELTAT = DELTAB+VELIMAG+1

    z = 1./a - 1.
    
    OR0 = omega_r
    OB0 = omega_b
    OC0 = omega_m - omega_b
    vstream = vbc
    ZSTART = zstart
    T = T_spline(z)
    xe = xe_spline(z)
    # a = 1. / (1. + z)
    cs = np.sqrt(BOLTZK*T/MUB/MPROTON)/MPCTOCM #; //sound speed for iosthermal gas
    YHE = 0.25
    fHe = .25*YHE/((1. - YHE) + .25*YHE);
    H = 100.*hconst*hubbleZ(z)/MPCTOKM;
    vbck = vstream/MPCTOKM/(a*(1.+ZSTART)); # //VSTREAM is in km/s
    OMEGAC = OC0/(OC0 +OB0 + OR0/a);
    OMEGAB = OB0/(OC0 +OB0 + OR0/a);
    h = hconst
    nH = CRITDENSITY*OB0*h*h*(1. - YHE)/(MPROTON);  # //Hydrogen density
    ucmb0 = aSB*TCMB*TCMB*TCMB*TCMB;
    tgamma0 = 3.*MELECTRON*LIGHTSPEED/(8.*SIGMAT*ucmb0);

    
    deriv = np.zeros(10)
    
    deriv[DELTAC+REAL] =  -y[DELTAC+VELREAL];
    deriv[DELTAC+IMAG] =  -y[DELTAC+VELIMAG];
    deriv[DELTAC+VELREAL] = - 1.5*H*H*(OMEGAC*y[DELTAC+REAL]+OMEGAB*y[DELTAB+REAL]) - 2.*H*y[DELTAC+VELREAL];
    deriv[DELTAC+VELIMAG] = - 1.5*H*H*(OMEGAC*y[DELTAC+IMAG]+OMEGAB*y[DELTAB+IMAG]) - 2.*H*y[DELTAC+VELIMAG];

    #baryons
    deriv[DELTAB+REAL] = vbck*k*costh*y[DELTAB+IMAG]/a - y[DELTAB+VELREAL];
    deriv[DELTAB+IMAG] = -vbck*k*costh*y[DELTAB+REAL]/a - y[DELTAB+VELIMAG];
    deriv[DELTAB+VELREAL] =  vbck*k*costh*y[DELTAB+VELIMAG]/a - 1.5*H*H*(OMEGAC*y[DELTAC+REAL]+OMEGAB*y[DELTAB+REAL]) - 2.*H*y[DELTAB+VELREAL] + cs*cs*k*k/(a*a)*(y[DELTAB+REAL] + y[DELTAT+REAL]);
    deriv[DELTAB+VELIMAG] =  -vbck*k*costh*y[DELTAB+VELREAL]/a - 1.5*H*H*(OMEGAC*y[DELTAC+IMAG]+OMEGAB*y[DELTAB+IMAG]) - 2.*H*y[DELTAB+VELIMAG] + cs*cs*k*k/(a*a)*(y[DELTAB+IMAG] + y[DELTAT+IMAG]); 
    

    eta1 = (1.+ fHe + xe);

    #//EQNS 9 and, if RECFAST is on, 10 are not currently used
    # deriv[9] = 0. #-2.*H*y[9] + y[10]/(eta1*tgamma0*(a*a*a*a)) *(TCMB/a - y[9]) #;//mean temperature (adiabatic and compton cooling)
    # deriv[10] = 0. #-alphaB_recomb(T, 0)*y[10]*y[10]*nH/(a*a*a);

    deriv[DELTAT+REAL] =  vbck*k*costh*y[DELTAT+IMAG]/a + 2./3.*(deriv[DELTAB+REAL]- vbck*k*costh*y[DELTAB+IMAG]/a) - xe/(eta1*tgamma0*(a*a*a*a)) *(TCMB/a/T)*y[DELTAT+REAL] # //temperature fluctuations --uses either recfast temp or y[9]
    deriv[DELTAT+IMAG] = -vbck*k*costh*y[DELTAT+REAL]/a + 2./3.*(deriv[DELTAB+IMAG] + vbck*k*costh*y[DELTAB+REAL]/a) - xe/(eta1*tgamma0*(a*a*a*a))*(TCMB/a/T)*y[DELTAT+IMAG] #; //temperature fluctuations (imaginary)

    
    deriv /= (H * a)  # dividing by aH does not work

    #print(deriv)
    
    return deriv


def hubbleZ(zCurrent):
    OmegaM = omega_m
    OmegaR = omega_r
    zCurrent += 1.0;
    temp = ((1.0 - OmegaM-OmegaR) 
	    + OmegaM*(zCurrent**3.0) + OmegaR*(zCurrent**4.0));
    return np.sqrt(temp)




def derivs_comb(z, y, k, T_spline, xe_spline, vbc, zstart):
    from py_vbc.utils import hubble
    
    # Calculate Hubble parameter
    H = (100.0*hconst)*hubble(z)/mpctokm
    
    # Put redshift dependence into density parameters
    # z1 = 1+z
    # o_c = (omega_m-omega_b)/(omega_m + omega_r*z1)
    # o_b = omega_b/(omega_m + omega_r*z1)

    z1 = 1+z
    o_m = omega_m * (z1 ** 3.)
    f_b = omega_b / omega_m
    f_c = 1. - f_b
    # o_b = omega_b * (z1 ** 3.)
    # o_c = o_m - o_b
    o_r = omega_r * (z1 ** 4.)
    o_mz = o_m / (o_m + o_r)
    o_b = f_b * o_mz
    o_c = f_c * o_mz

    Tcmb_z = Tcmb * z1  # Tcmb goes as 1/a

    # Convert vstream to Mpc to it is in the same units as k and
    # insert decaying dependence with redshift
    vbck = (z1 / (1+zstart)) * vbc / mpctokm

    # Get temperature and electron fraction values
    T = T_spline(z)
    xe = xe_spline(z)

    # Get photon density fluctuations over CDM fluctuations (this is
    # how the ICs are defined)
    # g = g_spline(k, z)
    # c = c_spline(k, z)
    # gc = g / c
    
    yhe = 0.25
    fhe = 0.25*yhe/((1.0-yhe) + 0.25*yhe)
    
    # Convenience variables
    mu = vbck*k*costh*z1
    alpha = 1.5*H**2*o_c
    beta = 1.5*H**2*o_b
    tau = (boltzk*T/(mub*mproton*mpctocm**2)) * k**2 * z1**2 
    # gamma = (xe/t_gamma)*(z1**5)
    gamma = (xe/t_gamma)*(z1**4)  # LC was the a-dependence overestimated?
    Tr = (Tcmb_z/T)
    eta = 1.0 + fhe +xe
    
    # dy contains equations for delta_c_dot, theta_c_dot, delta_b_dot,
    # theta_b_dot, delta_t_dot split into real and imaginary parts
    dy = np.zeros(5, dtype=complex)

    # delc = 0
    # delb = 2
    # vel = 1
    # delt = 4
    
    dy[delc] = -y[delc+vel]
    dy[delc+vel] = (- (alpha*y[delc] + beta*y[delb])
                      - 2.0*H*y[delc+vel])
    dy[delb] = -1j * mu*y[delb] - y[delb+vel]
    dy[delb+vel] = (- ((alpha*y[delc] + beta*y[delb]))
                    - (2*H*y[delb+vel])
                    - (1j * mu*y[delb+vel])
                    + tau*(y[delb] + y[delt]))

    # Writing the temperature fluctuations exactly as in CA
    # dy[deltat+real] = mu*y[deltat+imag] - (2.0/3.0)*y[deltab+velreal] - gamma*y[deltat+real]/eta
    # dy[deltat+real] = mu*y[deltat+imag] - (2.0/3.0)*(dy[deltab+real] - mu*y[deltab+imag]) - gamma*y[deltat+real]
    # dy[deltat+imag] = -mu*y[deltat+real] - (2.0/3.0)*y[deltab+velimag] - gamma*y[deltat+imag]/eta
    # dy[deltat+imag] = -mu*y[deltat+real] - (2.0/3.0)*(dy[deltab+imag] - mu*y[deltab+real]) - gamma*y[deltat+imag]

    # Naoz+ (2005) say this is only valid for z<=200, but Ahn showed
    # that actually its fine for z < 1000
    # dy[delt] = (2.0/3.0)*dy[delb] - gamma * Tr * y[delt]
    dy[delt] = (2.0/3.0)* (-1j * mu*y[delb] - y[delb+vel]) - gamma * Tr * y[delt]

    dy /= (-H*z1)
    
    return dy


def calc_derivs(k, vbc, zstart, zend, dz, verbose, integrator='LSODA'):
    T_spline, xe_spline = interpolate_recfast()
    # g_spline = interpolate_tf2d('g', zs)
    # c_spline = interpolate_tf2d('c', zs)
    
    y0 = set_ics(k, zstart, dz)
    y = np.zeros(shape=y0.shape)
    nk = y0.shape[0]
    z = np.zeros(shape=(nk, 7))

    if verbose: print('Solving evolution equations from z={0:4.2f} to z={1:4.2f} with v_bc={2:4.2f} ...'.format(zstart, zend, vbc))
    for i, ik in enumerate(k):
        # Progress, remove end='\r' for Python 2 compatibility
        if verbose:
            if i+1 != nk:
                print('    k = {0:.3f} Mpc^-1 [{1:d}/{2:d}]'.format(ik, i+1, nk), end='\r')
            else:
                print('    k = {0:.3f} Mpc^-1 [{1:d}/{2:d}]'.format(ik, i+1, nk))

        r = solve_ivp(fun=lambda z, y: derivs(z, y, ik, T_spline,
                                         xe_spline, vbc, zstart),
                      t_span=(zstart, zend),
                      t_eval=np.array([zend]),
                      rtol=1.0e-6, y0=y0[i, :],
                      method=integrator, jac=None)

        y[i, :] = np.transpose(r.y)

    if verbose:
        print('    ... done.')

    # Convert from real and imaginary components to magnitude
    z[:, 0] = np.sqrt(y[:, 0]**2 + y[:, 1]**2)  # delta_c
    z[:, 1] = np.sqrt(y[:, 2]**2 + y[:, 3]**2)  # theta_c
    z[:, 2] = np.sqrt(y[:, 4]**2 + y[:, 5]**2)  # delta_b
    z[:, 3] = np.sqrt(y[:, 6]**2 + y[:, 7]**2)  # theta_b
    z[:, 4] = np.sqrt(y[:, 8]**2 + y[:, 9]**2)  # delta_t
    z[:, 5] = theta_to_vel(y[:, [2, 3]], k, zend)  # v_c 
    z[:, 6] = theta_to_vel(y[:, [6, 7]], k, zend)  # v_b
    
    return z


def calc_derivs_cicsass(k, vbc, zstart, zend, dz, verbose, integrator='LSODA'):
    T_spline, xe_spline = interpolate_recfast()
    # g_spline = interpolate_tf2d('g', zs)
    # c_spline = interpolate_tf2d('c', zs)
    DELTAC = 0
    DELTAB = 4
    REAL = 0
    IMAG = 1
    VELREAL = 2
    VELIMAG = 3
    DELTAT = DELTAB+VELIMAG+1

    y0 = set_ics_cicsass(k, zstart, dz)
    y = np.zeros(shape=y0.shape)
    nk = y0.shape[0]
    z = np.zeros(shape=(nk, 7))

    if verbose: print('Solving evolution equations from z={0:4.2f} to z={1:4.2f} with v_bc={2:4.2f} ...'.format(zstart, zend, vbc))
    for i, ik in enumerate(k):
        # Progress, remove end='\r' for Python 2 compatibility
        if verbose:
            if i+1 != nk:
                print('    k = {0:.3f} Mpc^-1 [{1:d}/{2:d}]'.format(ik, i+1, nk), end='\r')
            else:
                print('    k = {0:.3f} Mpc^-1 [{1:d}/{2:d}]'.format(ik, i+1, nk))

        r = solve_ivp(fun=lambda a, y: derivs_cicsass(a, y, ik, T_spline,
                                                 xe_spline, vbc, zstart),
                      t_span=(1./(1. + zstart), 1./(1. + zend)),
                      t_eval=np.array([1./(1. + zend)]),
                      atol=1e-8, rtol=1e-18, first_step=1e-3, y0=y0[i, :],
                      method=integrator, jac=None)

        y[i, :] = np.transpose(r.y)

    if verbose:
        print('    ... done.')

    # Convert from real and imaginary components to magnitude
    z[:, 0] = np.sqrt(y[:, DELTAC+REAL]**2 +
                      y[:, DELTAC+IMAG]**2)  # delta_c
    z[:, 1] = np.sqrt(y[:, DELTAC+VELREAL]**2 +
                      y[:, DELTAC+VELIMAG]**2)  # theta_c
    z[:, 2] = np.sqrt(y[:, DELTAB+REAL]**2 +
                      y[:, DELTAB+IMAG]**2)  # delta_b
    z[:, 3] = np.sqrt(y[:, DELTAB+VELREAL]**2 +
                      y[:, DELTAB+VELIMAG]**2)  # theta_b
    z[:, 4] = np.sqrt(y[:, DELTAT+REAL]**2 + y[:, DELTAT+IMAG]**2)  # delta_t
    z[:, 5] = theta_to_vel(y[:, [DELTAC+VELREAL, DELTAC+VELIMAG]], k, zend)  # v_c 
    z[:, 6] = theta_to_vel(y[:, [DELTAB+VELREAL, DELTAB+VELIMAG]], k, zend)  # v_b
    
    return z



def calc_derivs_comp(k, vbc, zstart, zend, dz, verbose, integrator='LSODA'):
    T_spline, xe_spline = interpolate_recfast()
    # g_spline = interpolate_tf2d('g', zs)
    # c_spline = interpolate_tf2d('c', zs)
    
    y0 = set_ics_comb(k, zstart, dz)
    y = np.zeros(shape=y0.shape)
    nk = y0.shape[0]
    z = np.zeros(shape=(nk, 7))

    if verbose: print('Solving evolution equations from z={0:4.2f} to z={1:4.2f} with v_bc={2:4.2f} ...'.format(zstart, zend, vbc))
    for i, ik in enumerate(k):
        # Progress, remove end='\r' for Python 2 compatibility
        if verbose:
            if i+1 != nk:
                print('    k = {0:.3f} Mpc^-1 [{1:d}/{2:d}]'.format(ik, i+1, nk), end='\r')
            else:
                print('    k = {0:.3f} Mpc^-1 [{1:d}/{2:d}]'.format(ik, i+1, nk))

        r = solve_ivp(fun=lambda z, y: derivs_comb(z, y, ik, T_spline,
                                              xe_spline, vbc, zstart),
                      t_span=(zstart, zend),
                      t_eval=np.array([zend]),
                      rtol=1.0e-6, y0=y0[i, :],
                      method=integrator, jac=None)

        y[i, :] = np.transpose(r.y)

    if verbose:
        print('    ... done.')

    # Convert from real and imaginary components to magnitude
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
        return (- 1j * a * theta) / k
    
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
        return (- 1j * a * theta) / k
    

    comp_vel = convert(comp_theta, k, z)
    vel = np.abs(comp_vel)
    

    return vel

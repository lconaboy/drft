import os
import py_vbc

# Definitions
mpctokm = 3.08568e19
mpctocm = mpctokm*1e5
deltac = 0
deltab = 4
deltat = 8
real = 0
imag = 1
velreal = 2
velimag = 3

delc = 0
delb = 2
vel = 1
delt = 4

# Constants
Tcmb = 2.726              # CMB temperature at z=0
mproton = 1.6726e-24      # proton mass in g
melectron = 9.11e-28      # electron mass in g
mel_ev = 5.11e5           # electron mass in eV/c^2
lightspeed = 3.0e10       # c in cm/s
critdensity = 1.8791e-29  # rho_c in g/cm^3
sigmat = 0.665e-24        # Thomson cross-section in cm^2
yhe = 0.25                # Helium fraction?
boltzk = 1.3806e-16       # Boltzmann's constant in erg/K
mub = 1.22                # Mean molecular weight, primordial (?)
aSB = 7.56e-15            # u = aSB T^4 => Stefan's constant in ?
u_cmb = aSB*(Tcmb**4)     # energy stored CMB
t_gamma = 3*lightspeed*melectron/(8*sigmat*u_cmb)  # t_gamma in s


# Parameters (runtime, CICsASS)
boxsize = 0.2
size = 128
vstream = 0.0
zstart = 1000.0
zinit = 50.0
dz = 3.0
base, _ = os.path.split(py_vbc.constants.__file__)
camb_base = 'initSB_transfer_out'
tf_base = os.path.join(base, 'tfs', '{0}_z{1:03d}.dat')
rf_base = os.path.join(base, 'recfast', 'xeTrecfast.out')
justcamb = 0

# # Parameters (cosmological, Planck 2018)
# hconst = 0.673                # little h
# omega_m = 0.312                # matter density parameter
# omega_b = 0.049                # baryon density parameter
# sigma_8 = 0.812                # sigma_8
# ns = 0.965                  # n_s
# omega_r = 4.15e-5/(hconst**2)  # radiation density parameter, from
#                               # Dodelson (2.87) assuming 3 species of
#                               # relativistic neutrinos
# costh = 1

# Parameters (cosmological, CICsASS)
hconst = 0.71                  # little h
omega_m = 0.27                 # matter density parameter
omega_b = 0.046                # baryon density parameter
sigma_8 = 0.8                  # sigma_8
ns = 0.95                      # tilt of the power spectrum
omega_r = 4.15e-5/(hconst**2)  # radiation density parameter, from
                               # Dodelson (2.87) assuming 3 species of
                               # relativistic neutrinos
costh = 1                      # cos(theta) of the angle between the wavevector
                               # and v_bc 

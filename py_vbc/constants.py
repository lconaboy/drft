import os
import configparser

# config_fname = 'test_params.ini'
config_fname = 'planck2018_params.ini'
config_path, _ = os.path.split(__file__)

config = configparser.ConfigParser()
config.read(config_path + '/' + config_fname)

# Want printout?
verbose = True

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
# Only a few of the commented ones are needed, and the rest are
# input dynamically
# boxsize = 0.2
# size = 128
# vstream = 0.0
# zstart = 1000.0
# zinit = 50.0
# dz = 3.0

# THIS WILL BREAK IF ARCHIVED
import os
base, _ = os.path.split(__file__)
camb_base = config.get('filenames', 'camb_base')
tf_base = os.path.join(base, 'tfs', '{0}_z{1:03d}.dat')
rf_base = os.path.join(base, 'recfast', config.get('filenames', 'rf_base'))
justcamb = 0

# Parameters (cosmological)
hconst = float(config.get('cosmology', 'hconst'))
omega_m = float(config.get('cosmology', 'omega_m'))
omega_b = float(config.get('cosmology', 'omega_b'))
sigma_8 = float(config.get('cosmology', 'sigma_8'))
ns = float(config.get('cosmology', 'ns'))
omega_r = 4.15e-5/(hconst**2)  # Dodelson (2002) Eq. 2.86
costh = float(config.get('cosmology', 'costh'))

"""Physical constants and array-layout definitions used by py_vbc."""

# Want printout?
verbose = True

# Definitions
MPCTOKM = 3.08568e19
MPCTOCM = MPCTOKM * 1e5
DELTAC = 0
DELTAB = 4
DELTAT = 8
REAL = 0
IMAG = 1
VELREAL = 2
VELIMAG = 3

DELC = 0
DELB = 2
VEL = 1
DELT = 4

# Constants
TCMB = 2.726  # CMB temperature at z=0
MPROTON = 1.6726e-24  # proton mass in g
MELECTRON = 9.11e-28  # electron mass in g
MEL_EV = 5.11e5  # electron mass in eV/c^2
LIGHTSPEED = 3.0e10  # c in cm/s
CRITDENSITY = 1.8791e-29  # rho_c in g/cm^3
SIGMAT = 0.665e-24  # Thomson cross-section in cm^2
YHE = 0.25  # Helium fraction?
BOLTZK = 1.3806e-16  # Boltzmann's constant in erg/K
MUB = 1.22  # Mean molecular weight, primordial (?)
ASB = 7.56e-15  # u = ASB T^4 => Stefan's constant in ?
U_CMB = ASB * (TCMB**4)  # energy stored in the CMB
T_GAMMA = 3 * LIGHTSPEED * MELECTRON / (8 * SIGMAT * U_CMB)

# Interpolation redshifts for photon density fluctuations, not used currently
# zs = [1000, 900, 800, 700, 600, 500, 400, 300, 200]


__all__ = [
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
]

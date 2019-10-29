"""
TODO
- _power_spectrum
  + fft_sample_spacing
  + fft_sample_spacing_components

- linear_velocity_ps
"""
import numpy as np


class Cosmology:
    def __init__(self, a, params='planck2018'):
        """
        Class for instantiating a given set of cosmological
        parameters, and using those parameters to calculate some
        useful cosmological quantities.

        :param: a (float) - scale factor to calculate values at
        :param: params (str) - which set of cosmological parameters 
                               to use
        """

        # Check that the scale factor is realistic (this is for my
        # purposes, nothing wrong with a scale factor outside this
        # range, if you are looking to the future)
        assert (a >= 0.0) and (a <= 1.0), 'Scale factor outside [0.0, 1.0]'
        self.a = a

        # Instantiate some parameters
        self.p = self.parameters(params)

        # Calculate all the useful cosmological quantities at scale
        # factor a
        self.H0 = self.p['h'] * 100.0  # km/s/Mpc
        self.E_a = self.E(self.a)
        self.H_a = self.H()
        self.D_a = self.D()
        self.f_a = self.f()


    def E(self, a):
        """
        Calculates the prefactor to the Hubble parameter which contains
        the density parameters, at a scale factor a. H(a) = H0 *
        E(a). See e.g. Peacock (1999) (7.77).

        """
        return np.sqrt(self.p['omega_m']/self.a**3 + self.p['omega_l'])

    
    def H(self):
        """
        Calculates the Hubble parameter at a scale factor a. See
        e.g. Peacock (1999) (3.38).

        """
        return self.H0 * self.E_a


    def D(self):
        """
        Unnormalised growth factor D1(a). See e.g. Dodelson (2003) (7.77).
        
        """
        from scipy.integrate import quad
        i = quad(lambda _a: 1 / (_a * self.E(_a))**3.0, self.a, 0.0)[0]

        return (5.0 * self.p['omega_m'] * self.E_a * i)/2.0 


    def f(self):
        """
        The equation for the growth factor as a function of scale factor,
        using f ~ \Omega^0.6. Used in the continuity equation. See
        e.g. Lahav et al. (1991).
        
        """
        omega = (self.H0 / self.H_a)**2 * (self.a ** -3)  * self.p['omega_m']

        return omega ** 0.6


    def a_to_z(self, a):
        return 1.0/a - 1.0


    def z_to_a(self, z):
        return 1.0/(1.0 + z)


    def parameters(self, params):
        if params == 'planck2015':
            return {'omega_m': 0.3156,
                    'omega_b': 0.04917,
                    'omega_l': 0.6844,
                    'h': 0.6727}

        elif params == 'planck2018':
            return {'omega_m': 0.317,
                    'omega_b': 0.049,
                    'omega_l': 0.683,
                    'h': 0.673,
                    'ns': 0.965,
                    's8': 0.812,
                    'As': 2.10e-9}

################################################################################
###############################  OLD BELOW HERE  ###############################
################################################################################

def _fft_sample_spacing(N, boxsize):
    """
    Return the sample spacing in Fourier space, given some symmetric 3D box in real space
    with N elements per dimension and length L.
    See https://gitorious.org/bubble/szclusters/commit/da1402ef95f4d40c28f53f88c99bf079063308c7
    """
    kk = np.zeros([N, N, N], dtype=np.float32)
    kx = np.zeros([N, N, N], dtype=np.float32)
    ky = np.zeros([N, N, N], dtype=np.float32)
    kz = np.zeros([N, N, N], dtype=np.float32)
    kx, ky, kz = _fft_sample_spacing_components(N)
    fac = (2. * np.pi / boxsize)
    kk = np.sqrt(kx ** 2. + ky ** 2. + kz ** 2.) * fac
    return kk

def _fft_sample_spacing_components(N):
    """
    Return the sample spacing in Fourier space, given some symmetric 3D box in real space
    with N elements per dimension and length L.
    See https://gitorious.org/bubble/szclusters/commit/da1402ef95f4d40c28f53f88c99bf079063308c7
    """
    
    NN = np.zeros(N, dtype=np.int32)
    kx = np.zeros([N, N, N], dtype=np.float32)
    ky = np.zeros([N, N, N], dtype=np.float32)
    kz = np.zeros([N, N, N], dtype=np.float32)

    NN = (N * np.fft.fftfreq(N, 1.)).astype(np.int32)
    for i in range(N):
        kx[NN[i], :, :] = NN[i]
        ky[:, NN[i], :] = NN[i]
        kz[:, :, NN[i]] = NN[i]
    return kx, ky, kz


def _window_function(N, p):
    ''' Calculate CIC smoothing window function for PS estimation.
    '''
    W = np.zeros([N, N, N], dtype=np.float32)
    kx = np.zeros([N, N, N], dtype=np.float32)
    ky = np.zeros([N, N, N], dtype=np.float32)
    kz = np.zeros([N, N, N], dtype=np.float32)
    fp = float(p)
    kny = N/2  # Nyquist frequency
    kx, ky, kz = _fft_sample_spacing_components(N)
    
    #W = ( np.sinc( (np.pi*kx) / (2.*kny) ) * np.sinc( (np.pi*ky) / (2.*kny) ) * np.sinc( (np.pi*kz) / (2.*kny) ) )**p
    W = ( np.sinc( np.pi*kx/2.*kny ) * np.sinc( np.pi*ky/2.*kny ) * np.sinc( np.pi*kz/2.*kny ) ) ** fp
    return W

def _cic_window_function(N):
    ''' Calculate CIC smoothing window function for PS estimation.
    '''
    W = np.zeros([N, N, N], dtype=np.float32)
    kx = np.zeros([N, N, N], dtype=np.float32)
    ky = np.zeros([N, N, N], dtype=np.float32)
    kz = np.zeros([N, N, N], dtype=np.float32)
    kny = N/2  # Nyquist frequency
    kx, ky, kz = _fft_sample_spacing_components(N)
    
    #W = ( np.sinc( (np.pi*kx) / (2.*kny) ) * np.sinc( (np.pi*ky) / (2.*kny) ) * np.sinc( (np.pi*kz) / (2.*kny) ) )**p
    W = ( 1. - (2./3.) * np.sin((np.pi*kx)/(2.*kny))**2 ) * ( 1. - (2./3.) * np.sin((np.pi*ky)/(2.*kny))**2 ) * ( 1. - (2./3.) * np.sin((np.pi*kz)/(2.*kny))**2 )
    return W


def linear_velocity_ps(k, delta_k, **cosmo):
    '''
    Compute velocity power spectrum using equations in Iliev et al. 2007
    k - wavenumber
    delta_k - dimensionless density power spectrum
    H0 - Hubble constant today
    omegam0 - Total matter density today
    omegal0 - Totel dark energy density today
    z - redshift to compute velocity power spectrum
    '''
    # First, compute the growth factor
    D_z_growth = D_z(**cosmo)
    H0 = cosmo['h'] * 100.
    # Compute E(z) = H(z)/H0 term
    E_z = hzoverh0(**cosmo)
    # Compute and return linear power spectrum
    t1 = delta_k ** 2 / k ** 2
    t2 = (9 * H0 ** 2 * cosmo['omega_M_0'] **
          2 * (1. + cosmo['z']) ** 4) / (4 * E_z ** 2)
    t3 = (1 - (5 / (3 * (1. + cosmo['z']) * D_z_growth))) ** 2

    return np.sqrt(t1 * t2 * t3)


def D_z(**cosmo):
    """
    Unnormalised linear growth factor
    D(a) = 5 omegam / 2 H(a) / H(0) * integral[0:a] [da / [a H(a) H0]**3]
    equation  from  peebles 1980 (or e.g. eq. 8 in lukic et al. 2008) """

    import scipy.integrate

    omegam0 = cosmo['omega_M_0']
    omegal0 = cosmo['omega_lambda_0']

    if (abs(omegam0 + omegal0 - 1.) > 1.e-4):
        raise RuntimeError(
            "Linear growth factors can only be calculated for flat cosmologies")

    a = 1. / (1. + cosmo['z'])

    # 1st calc. for z=z
    lingrowth = scipy.integrate.quad(_lingrowthintegrand, 0., a, (omegam0))[0]
    lingrowth *= 5. / 2. * omegam0 * hzoverh0(**cosmo)
    return lingrowth


def hzoverh0(**cosmo):
    """ returns: H(a) / H0  = [omegam/a**3 + (1-omegam)]**0.5 """
    return np.sqrt(cosmo['omega_M_0'] * np.power(cosmo['aexp'], -3) + (1. - cosmo['omega_M_0']))


def _lingrowthintegrand(a, omegam):
    """ (e.g. eq. 8 in lukic et al. 2008)   returns: da / [a*H(a)/H0]**3 """
    return np.power((a * hzoverh0(**{'aexp': a, 'omega_M_0': omegam})), -3)


def lingrowthfac(red, return_norm=False, **cosmo):
    """
    returns: linear growth factor, b(a) normalized to 1 at z=0, good for flat lambda only
    a = 1/1+z
    b(a) = Delta(a) / Delta(a=1)   [ so that b(z=0) = 1 ]
    (and b(a) [Einstein de Sitter, omegam=1] = a)

    Delta(a) = 5 omegam / 2 H(a) / H(0) * integral[0:a] [da / [a H(a) H0]**3]
    equation  from  peebles 1980 (or e.g. eq. 8 in lukic et al. 2008) """
    # need to add w ~= , nonflat, -1 functionality

    import scipy.integrate

    if (abs(cosmo['omega_M_0'] + cosmo['omega_lambda_0'] - 1.) > 1.e-4):
        raise RuntimeError(
            "Linear growth factors can only be calculated for flat cosmologies")

    # 1st calc. for z=z
    lingrowth = scipy.integrate.quad(
        _lingrowthintegrand, 0., cosmo['aexp'], (cosmo['omega_M_0']))[0]
    lingrowth *= 5. / 2. * cosmo['omega_M_0'] * \
        hzoverh0(**cosmo)

    # then calc. for z=0 (for normalization)
    a0 = 1.
    lingrowtha0 = scipy.integrate.quad(
        _lingrowthintegrand, 0., a0, (cosmo['omega_M_0']))[0]
    lingrowtha0 *= 5. / 2. * \
        cosmo['omega_M_0'] * \
        hzoverh0(**{'aexp': a0, 'omega_M_0': cosmo['omega_M_0']})

    lingrowthfactor = lingrowth / lingrowtha0
    if return_norm:
        return lingrowthfactor, lingrowtha0
    else:
        return lingrowthfactor

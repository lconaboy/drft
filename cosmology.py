import numpy as np

class Cosmology:
    def __init__(self, a, params='planck2018'):
        """
        Class for instantiating a given set of cosmological
        parameters, and using those parameters to calculate some
        useful cosmological quantities.

        :param: a (float) - scale factor to calculate values at
        :param: params (str or list) - which set of cosmological parameters 
                                       to use, either a string or list of
                                       parameters from grafic header
        """

        # Check that the scale factor is realistic (this is for my
        # purposes, nothing wrong with a scale factor outside this
        # range, if you are looking to the future)
        assert (a >= 0.0) and (a <= 1.0), 'Scale factor outside [0.0, 1.0]'
        self.a = a

        # Instantiate the parameters
        self.p = self.parameters(params)

        # Calculate all the useful cosmological quantities at scale
        # factor a
        self.H0 = self.p['h'] * 100.0  # km/s/Mpc
        self.E_a = self.E(self.a)  # H(a)/H0
        self.H_a = self.H()  # Hubble parameter
        self.D_a = self.D()  # growth factor
        self.d_a = self.d()  # normalised growth factor
        self.f_a = self.f()  # f(\Omega)


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

        return 2.5 * self.p['omega_m'] * self.E_a * i 


    def d(self):
        """Growth factor normalised to 1 at z=0"""
        from scipy.integrate import quad
        i = quad(lambda _a: 1 / (_a * self.E(_a))**3.0, 1.0, 0.0)[0]
        i *= 2.5 * self.p['omega_m'] * self.E_a
        
        return self.D_a / i

    
    def f(self):
        """
        The equation for the growth factor as a function of scale factor,
        using f ~ \Omega^0.6. Used in the continuity equation. See
        e.g. Lahav et al. (1991).
        
        """
        omega = self.p['omega_m'] / self.a**3 / self.E_a**2
        
        return omega ** 0.6


    def a_to_z(self, a):
        return 1.0/a - 1.0


    def z_to_a(self, z):
        return 1.0/(1.0 + z)


    def parameters(self, params):
        if type(params).__name__ == 'str':
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
        else:
            return {'omega_m': params[0],
                    'omega_l': params[1],
                    'h': params[2]}


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
    """Return the sample spacing in Fourier space, given some symmetric 3D box in real space
    with N elements per dimension and length L.
    See https://gitorious.org/bubble/szclusters/commit/da1402ef95f4d40c28f53f88c99bf079063308c7

    LC (26/4) -- after a quick check, these seem to return the kx, ky
    and kz components in the expected order (so no need to swap kx and
    ky as in linear.py)

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

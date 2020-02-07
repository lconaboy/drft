import numpy as np
import matplotlib.pyplot as plt
# import scipy.fftpack as fft
import numpy.fft as fft

# from grafic_tools import Snapshot
from cosmology import Cosmology

class Continuity:
    def __init__(self, ic, tf_path=None, params=None, approx=False, real=True):
        """
        Class to calculate peculiar velocities from overdensities using
        the continuity equation

        :param ic: 
            ic_deltab file, created using grafic_tools
        :type ic: 
            Snapshot object, created using grafic_tools
        :param tf_path: 
            path to CAMB transfer functions used to generate ICs,
            otherwise assume that transfer functions/fitting
            functions with the same amplitude for baryons and dark
            matter perturbations are used
        :param params: 
            either None, to get parameters from Snapshot object or 
            string to use defined parameters
        :type ic: 
            None or str
        :param approx: 
            whether to use the approximate or full version of the 
            continuity equation
        :type approx: 
            bool
        :param real: 
            whether to use a real FFT or regular FFT
        :type real: 
            bool
        :returns: 
        :rtype:

        """
        
        # Set up some initial constants
        self.z = ic.cosmo['z']     # Redshift
        self.a = ic.cosmo['aexp']  # Scale factor
        self.N = ic.N              # Number of samples
        self.dx = ic.dx            # Initial grid spacing
        self.L = ic.boxsize        # Box size in Mpc
        self.tf_path = tf_path     # Path to CAMB TFs

        if params is None:
            # Get parameters from grafic header
            omega_m = ic.cosmo['omega_m']
            omega_l = ic.cosmo['omega_l']
            h = ic.cosmo['h']
            
            # Collect parameters to pass to Cosmology
            params = [omega_m, omega_l, h]
            
            # Calculate the relevant cosmological quantities
            self.c = Cosmology(self.a, params=params)
        else:
            assert type(params).__name__ == 'str', 'To get parameters from the grafic file set params=None, otherwise specify cosmology with a string'
            self.c = Cosmology(self.a, params=params)

        # Extract overdensity field and take the Fourier transform
        # assert ic.field == 'deltab', 'Should be looking at deltab field'
        self.delta_x = ic.load_box()

        print('mean of delta_x', np.mean(self.delta_x))

        # Use real FFT? \delta(x) should be purely real, so should be
        # fine to use RFFT.
        self.real = real
        if self.real:
            self.delta_k = fft.rfftn(self.delta_x)
        else:
            self.delta_k = fft.fftn(self.delta_x)

        # Calculate the FFT sample spacing
        self.fft_sample_spacing()

        # Generate the TF ratios
        self.transfer_functions()
        
        # Compute the unscaled velocities
        self.unscaled_velocity()

        # Scale the velocities using the appropriate method
        if approx:
            self.scaled_velocity_approx()
        else:
            self.scaled_velocity_full()

        # Transform to real velocities
        self.realise_velocity()


    def fft_sample_spacing(self):
        """
        Calculates the sample spacing for a Fourier transform. If RFFT is
        used (real=True), then the last axis contains half as many
        components as the other two. Calculates self.k_min
        (fundamental mode) and self.k_max (Nyquist frequency).
        """

        # Calculate the fundamental and Nyquist frequencies
        self.k_min = 2*np.pi / self.L
        self.k_max = np.pi / self.dx

        # Sample spacing in real space
        # d = self.dx / (2*np.pi)
        d = 1.0 / (self.k_min * self.N)

        # Frequencies along one axis, for a full FFT
        self.k_lin = fft.fftfreq(self.N, d=d)
 
        print('max abs(k_lin)', np.max(np.abs(self.k_lin)))
        print('k_max', self.k_max)
        
        if self.real:
            # Frequencies along the half axis, for a real FFT
            self.k_rlin = fft.rfftfreq(self.N, d=d)
        else:
            self.k_rlin = self.k_lin

        # Generate grid of k values, for use as the vector later
        self.kx, self.ky, self.kz = np.meshgrid(self.k_lin, self.k_lin, self.k_rlin)

        # Magnitude of k values at each point
        self.k2 = self.kx**2. + self.ky**2. + self.kz**2.
        self.k = np.sqrt(self.k2)


    def transfer_functions(self, species='b'):
        """
        Interpolate CAMB transfer functions to calculate the offset factor
        for two-component fluids
        """
        from scipy.interpolate import interp1d

        if self.tf_path is not None:
            # TODO - generalise for older versions of CAMB i.e. different
            # number of columns
            kh, tc, tb, tt = np.loadtxt(self.tf_path, unpack=True, usecols=(0, 1, 2, 6))
            # Try using velocity TFs
            kh, tc, tb = np.loadtxt(self.tf_path, unpack=True, usecols=(0, 10, 11))
            tt = ((self.c.p['omega_m'] - self.c.p['omega_b'])/self.c.p['omega_m'] * tc) + (self.c.p['omega_b']/self.c.p['omega_m'] * tb)

            # Convert kh from h/Mpc to 1/Mpc
            k = kh * self.c.p['h']

            # Generate splines of TFs
            tc_spline = interp1d(k, tc)  # CDM
            tb_spline = interp1d(k, tb)  # baryons
            tt_spline = interp1d(k, tt)  # total

            k_lin = self.k_lin
            k_rlin = self.k_rlin

            # Set zero k values to very small
            k_lin[0] = 1.001*min(k)
            k_rlin[0] = 1.001*min(k)
            # Also take absolute values of k
            k_lin = np.abs(k_lin)
            k_rlin = np.abs(k_rlin)

            # Calculate the values at the sampled k-values
            tc_kl = tc_spline(k_lin)
            tc_krl = tc_spline(k_rlin)
            tb_kl = tb_spline(k_lin)
            tb_krl = tb_spline(k_rlin)
            tt_kl = tt_spline(k_lin)
            tt_krl = tt_spline(k_rlin)

            # Calculate the ratios
            if species == 'b':
                # ratio_l = tb_kl/tt_kl
                # ratio_l = tc_kl - tb_kl
                ratio_l = np.mean(tb_kl/tt_kl)
                # ratio_l = np.mean(tb/tt)
                # ratio_rl = tb_krl/tt_krl
                # ratio_rl = tc_krl - tb_krl
                ratio_rl = np.mean(tb_krl/tt_krl)
                # ratio_rl = np.mean(tb/tt)
            else:
                # ratio_l = np.mean(tc_kl/tt_kl)
                ratio_l = np.mean(tc/tt)
                # ratio_rl = np.mean(tc_krl/tt_krl)
                ratio_rl = np.mean(tc/tt)


        else:
            # If CAMB TFs aren't used then the ratio between TFs is just one
            ratio_l = 1.0 # np.ones(len(self.k_lin))
            ratio_rl = 1.0 # np.ones(len(self.k_rlin))

        # Now generate the grid of values
        # self.rx, self.ry, self.rz = np.meshgrid(ratio_l, ratio_l, ratio_rl)
        self.rx, self.ry, self.rz = (ratio_l, ratio_l, ratio_rl)
        # print(self.rx)
        # print(self.ry)
        # print(self.rz)
    
    def unscaled_velocity(self):
        """
        Calculate the common factor to the velocities i.e.

        v \propto i * \delta_k * k / k^2

        """
        # Components of the velocity field without any cosmological
        # quantities being involved
        vx = 1j * self.delta_k * self.kx / self.k2
        vy = 1j * self.delta_k * self.ky / self.k2
        vz = 1j * self.delta_k * self.kz / self.k2

        # Above, v(k=0) will be NaN, so convert any NaNs to zeros
        vx = np.nan_to_num(vx)
        vy = np.nan_to_num(vy)
        vz = np.nan_to_num(vz)

        # Return velocities
        self.velocity_k = [vx, vy, vz]


    def scaled_velocity_approx(self):
        """
        Calculate the scaled velocity using the approximate expression for
        the continuity equation as in the draft
        """
        # Prefactor for the scaling
        p1 = self.a * self.c.H_a * self.c.f_a

        # Now calculate the scaled velocity
        self.velocity_ks = [p1 * v for v in self.velocity_k]


    def scaled_velocity_full(self):
        """
        Calculate the scaled velocity using the full expression for the
        continuity equation as in Iliev at al. (2007)
        """
        # p1 and p2 are part of d/dt(D) / D = p1*p2
        p1 = - (3.0 * self.c.H0 * self.c.p['omega_m']) / (self.c.E_a * self.a**3.0)
        p2 = 1 - (5.0 * self.a)/(3.0 * self.c.D_a)

        # Now calculate the scaled velocity
        self.velocity_ks = [p1 * p2 * self.a * v for v in self.velocity_k]


    def realise_velocity(self):
        """
        Transform the velocity field from Fourier space to real space
        """
        if self.real:
            self.velocity_r = [fft.irfftn(v) for v in self.velocity_ks]
            # self.velocity_r = [np.nan_to_num(v) for v in self.velocity_r]
        else:
            self.velocity_r = [fft.ifftn(v).real for v in self.velocity_ks]
            

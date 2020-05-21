"""
NOTES
   - swapped vx/px and vy/py
"""

import sys
import numpy as np
import numpy.fft as fft
from cosmology import Cosmology

class Linear:
    def __init__(self, ic, qty, params=None, approx=False, real=True):
        """
        Class to calculate peculiar velocities from overdensities using
        the continuity equation

        :param ic: 
            ic_deltab file, created using grafic_tools
        :type ic: 
            Snapshot object, created using grafic_tools
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

        self.approx = approx
        self.real = real
          
        # Set up some initial constants
        self.z = ic.cosmo['z']     # Redshift
        self.a = ic.cosmo['aexp']  # Scale factor
        self.N = ic.N              # Number of samples
        self.dx = ic.dx            # Initial grid spacing
        self.L = ic.boxsize        # Box size in Mpc

        if params is None:
            print('Using parameters from grafic header')
            
            # Get parameters from grafic header
            omega_m = ic.cosmo['omega_m']
            omega_l = ic.cosmo['omega_l']
            h = ic.cosmo['h']

            print('---- h      :', h)
            print('---- omega_m:', omega_m)
            print('---- omega_l:', omega_l)
            
            # Collect parameters to pass to Cosmology
            params = [omega_m, omega_l, h]
            
            # Calculate the relevant cosmological quantities
            self.c = Cosmology(self.a, params=params)
        else:
            assert type(params).__name__ == 'str', 'To get parameters from the grafic file set params=None, otherwise specify cosmology with a string'

            print('Using ' + params + ' params file')
            self.c = Cosmology(self.a, params=params)

            print('---- h      :', self.c.p['h'])
            print('---- omega_m:', self.c.p['omega_m'])
            print('---- omega_l:', self.c.p['omega_l'])
            
            
        # Extract overdensity field and take the Fourier transform
        self.delta_x = ic.load_box()

        # Use real FFT? \delta(x) should be purely real, so should be
        # fine to use RFFT.
        if self.real:
            self.delta_k = fft.rfftn(self.delta_x)
        else:
            self.delta_k = fft.fftn(self.delta_x)

        # Calculate the FFT sample spacing
        self.set_sample_spacing()
        
        # Compute the unscaled velocities
        assert qty in ['vel', 'pos'], 'qty should be vel or pos'
        if qty == 'vel':
            print('Calculating velocities')

            if self.approx:
                print('---- using f(\Omega) for scaling')
            else:
                print('---- using full expression for scaling')
            
            self.set_velocity()
            self.scale_velocity()
            self.realise_velocity()

        elif qty == 'pos':
            print('Calculating positions')

            self.set_position()
            self.realise_position()


    def set_sample_spacing(self):
        """Calculates the sample spacing for a Fourier transform. If RFFT is
        used (real=True), then the last axis contains half as many
        components as the other two. Calculates self.k_min
        (fundamental mode) and self.k_max (Nyquist frequency). 

        kx and ky are swapped from the way you might expect since the
        files are read and written as Fortran binary files, which are
        column-major (as opposed to Python which is in the 'C' style,
        row-major)

        """

        # Calculate the fundamental and Nyquist frequencies
        self.k_min = 2*np.pi / self.L
        self.k_max = np.pi / self.dx

        # Sample spacing in real space
        d = 1.0 / (self.k_min * self.N)

        # Frequencies along one axis, for a full FFT
        self.k_lin = fft.fftfreq(self.N, d=d)

        # If we are using real FFTs then the thrid axis is half as big
        if self.real:
            # Frequencies along the half axis, for a real FFT
            self.k_rlin = fft.rfftfreq(self.N, d=d)
        else:
            self.k_rlin = self.k_lin

        # Generate grid of k values, for use as the vector later
        self.ky, self.kx, self.kz = np.meshgrid(self.k_lin, self.k_lin, self.k_rlin)

        # Magnitude of k values at each point
        self.k2 = self.kx**2. + self.ky**2. + self.kz**2.
        self.k = np.sqrt(self.k2)

    
    def set_velocity(self):
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
        self.vel = [vx, vy, vz]


    def scale_velocity(self):
        """Scale the velocity using the approximate f(\Omega) method or using
        the full expression as in Iliev+ (2007)
        """

        # Use f(\Omega)
        if self.approx:
            # Prefactor for the scaling
            p1 = self.a * self.c.H_a * self.c.f_a

            # Now calculate the scaled velocity
            for i in range(3):
                self.vel[i] *= p1
    
        # Use Eq. (10) in Iliev+ (2007)
        else:
            # p1 and p2 are part of d/dt(D) / D = p1*p2
            p1 = - (3.0 * self.c.H0 * self.c.p['omega_m']) / (self.c.E_a * self.a**3.0)
            p2 = 1 - (5.0 * self.a)/(3.0 * self.c.D_a)

            # Now calculate the scaled velocity
            for i in range(3):
                self.vel[i] *= p1 * p2 * self.a


    def realise_velocity(self):
        """Transform the velocity field from Fourier space to real space, in
        units of km/s

        """
        if self.real:
            for i in range(3):
                self.vel[i] = fft.irfftn(self.vel[i])
            # self.vel_r = [np.nan_to_num(v) for v in self.vel_r]
        else:
            for i in range(3):
                self.vel[i] = fft.ifftn(self.vel[i]).real


    def set_position(self):
        """
        Calculate the common factor to the positions i.e.

        p \propto i * \delta_k * k / k^2

        """
        # Components of the velocity field without any cosmological
        # quantities being involved
        px = 1j * self.delta_k * self.kx / self.k2
        py = 1j * self.delta_k * self.ky / self.k2
        pz = 1j * self.delta_k * self.kz / self.k2

        # Above, v(k=0) will be NaN, so convert any NaNs to zeros
        px = np.nan_to_num(px)
        py = np.nan_to_num(py)
        pz = np.nan_to_num(pz)

        # Return positions
        self.pos = [px, py, pz]

    
    def realise_position(self):
        """Transform the position field from Fourier space back to real space,
        and convert to units of comoving Mpc/h for output

        """
        if self.real:
            for i in range(3):
                self.pos[i] = fft.irfftn(self.pos[i])
            # self.pos_r = [np.nan_to_num(v) for v in self.pos_r]
        else:
            for i in range(3):
                self.pos[i] = fft.ifftn(self.pos[i]).real * self.c.p['h']
        

def generate(ic, qty, params=None, approx=False, real=True):
    """Convenience function for calculating the relevant
    quantity. Returns the quantity."""

    assert qty in ['vel', 'pos'], 'qty should be vel or pos'

    c = Linear(ic=ic, qty=qty, params=params,
                   approx=approx, real=real)
    
    if qty == 'vel':
        return c.vel
    
    elif qty == 'pos':
        assert ic.field == 'deltac', 'Should be using deltac field to calculate positions'
        return c.pos

    

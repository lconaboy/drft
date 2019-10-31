import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fft

from grafic_ics import Snapshot
from cosmology import Cosmology

# Load the overdensity ICs
ic = Snapshot('test_ics/', 7, 'deltab')
field = ic.load_box()

# Fourier tranform of overdensity
delta_k = fft.fftn(field)

# Compare to original fields, should be identical
vx_music = Snapshot('test_ics/', 7, 'velbx').load_box()
vy_music = Snapshot('test_ics/', 7, 'velby').load_box()
vz_music = Snapshot('test_ics/', 7, 'velbz').load_box()
v_music = np.sqrt(vx_music**2 + vy_music**2 + vz_music**2)

class Continuity:
    def __init__(self, ic, z, approx=False):
        # Set up some initial constants
        self.z = z                 # Redshift
        self.a = 1.0 / (1.0 + z)   # Scale factor
        self.N = ic.N              # Number of samples
        self.dx = ic.dx            # Initial grid spacing
        self.L = ic.boxsize        # Box size in Mpc

        # Calculate the relevant cosmological quantities
        self.c = Cosmology(self.a, params='planck2015')

        # Extract overdensity field and take the Fourier transform
        self.delta_x = ic.load_box()
        self.delta_k = fft.fftn(self.delta_x)

        # Calculate the FFT sample spacing
        self.set_fft_sample_spacing()

        # Compute the velocities
        self.unscaled_velocity()
        
        if approx:
            self.scaled_velocity_approx()
        else:
            self.scaled_velocity_full()


    def set_fft_sample_spacing(self):
        """
        The following was originally authored by Phil Bull.

        https://gitlab.com/cosmobubble/szclusters/tree/master

        Calculate the sample spacing in Fourier space, given some symmetric 3D 
        box in real space, with 1D grid point coordinates 'x'.
        """
        self.kx = np.zeros(shape=(self.N, self.N, self.N))
        self.ky = np.zeros(shape=(self.N, self.N, self.N))
        self.kz = np.zeros(shape=(self.N, self.N, self.N))
        NN = (self.N * fft.fftfreq(self.N, 1.)).astype(int)
        NNN = fft.fftfreq(self.N, self.dx / (2*np.pi))
        # NNN = fft.fftfreq(self.N, self.dx)

        # LC - here I have swapped the order kx and kz
        for i, j in zip(NN, NNN):
                self.kx[:,:,i] = j
                self.ky[:,i,:] = j
                self.kz[i,:,:] = j

        # LC - multiply all modes by scaling factor
        # self.fac = 2*np.pi / self.L

        self.k = np.sqrt(self.kx**2. + self.ky**2. + self.kz**2.)
        # self.k *= self.fac



    def unscaled_velocity(self):
        """
        Realise the (unscaled) velocity field in Fourier space. See
        Dodelson Eq. 9.18 for an expression; we factor out the
        time-dependent quantities here. They can be added at a later
        stage.
        """

        # If the FFT has an even number of samples, the most negative frequency 
        # mode must have the same value as the most positive frequency mode. 
        # However, when multiplying by 'i', allowing this mode to have a 
        # non-zero real part makes it impossible to satisfy the reality 
        # conditions. As such, we can set the whole mode to be zero, make sure 
        # that it's pure imaginary, or use an odd number of samples. Different 
        # ways of dealing with this could change the answer!
        # if self.N % 2 == 0: # Even no. samples
                # Set highest (negative) freq. to zero
                # self.kx[self.kx == np.min(self.kx)] = 0.0
                # self.ky[self.ky == np.min(self.ky)] = 0.0
                # self.kz[self.kz == np.min(self.kz)] = 0.0

        # Get squared k-vector in k-space (and factor in scaling from kx, ky, kz)
        k2 = self.k ** 2.0

        # Calculate components of A (the unscaled velocity)
        Ax = 1j * self.delta_k * self.kx / k2 # * self.fac / k2
        Ay = 1j * self.delta_k * self.ky / k2 # * self.fac / k2
        Az = 1j * self.delta_k * self.kz / k2 # * self.fac / k2
        Ax = np.nan_to_num(Ax)
        Ay = np.nan_to_num(Ay)
        Az = np.nan_to_num(Az)
        self.velocity_k = (Ax, Ay, Az)


    def scaled_velocity_approx(self):
        """
        The following was originally authored by Phil Bull.

        https://gitlab.com/cosmobubble/szclusters/tree/master

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


# Use the approximate or full continuty equation?
approx = False

# For plots
fn = 'full'
if approx: fn = 'approx'

# Calculate the velocity using the continutity equation
c = Continuity(ic, z=200, approx=approx)

# Take the real part of the IFFT, this is fine since the imaginary
# parts are due to numerical noise and are ~1e-6
v = [np.fft.ifftn(v_i).real for v_i in c.velocity_ks]
v_m = [vx_music, vy_music, vz_music]

# Fit a line to the values to check their relation
vx_r = v[0][:, :, 0].ravel()
vx_mr = vx_music[:, :, 0].ravel()

from scipy.optimize import curve_fit

def y(x, m, c):
    return m*x + c

popt, pcov = curve_fit(y, vx_r, vx_mr)


"""
Plotting
"""

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))
im0 = axes[0].imshow(v[0][:, :, 0])
plt.colorbar(ax=axes[0], mappable=im0, fraction=0.046, pad=0.04)
axes[0].set_title('vx[:, :, 0]')
im1 = axes[1].imshow(vx_music[:, :, 0])
plt.colorbar(ax=axes[1], mappable=im1, fraction=0.046, pad=0.04)
axes[1].set_title(r'vx\_music[:, :, 0]')
im2 = axes[2].imshow(np.log10(v[0][:, :, 0]/vx_music[:, :, 0]), cmap='hot')
plt.colorbar(ax=axes[2], mappable=im2, fraction=0.046, pad=0.04)
axes[2].set_title('log$_{10}$(vx/vx\_music)')
plt.suptitle('{}'.format(fn))
fig.tight_layout() 
fig.savefig('velocity_slice_comparison_{}.pdf'.format(fn))
plt.show()

labels = ['x', 'y', 'z']
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))

for i, ax in enumerate(axes.ravel()):
    v_r = v[i].ravel()
    v_mr = v_m[i].ravel()

    popt, pcov = curve_fit(y, v_mr, v_r)
    xmin = np.min(v_mr)
    xmax = np.max(v_mr)

    x1 = [xmin, xmax]
    y1 = [y(xmin, popt[0], popt[1]), y(xmax, popt[0], popt[1])]
    
    ax.plot(v_mr[::100], v_r[::100], 'kx', label='__nolabel__')
    ax.plot(x1, y1, c='seagreen', label='y = {0:3.4f}x + {1:3.4f}'.format(popt[0], popt[1]))
    ax.set_title(labels[i])
    ax.legend()
    ax.set_xlabel('music')
    ax.set_ylabel('continuity')
plt.suptitle('{}'.format(fn))
plt.savefig('continuity_{}.png'.format(fn), dpi=600)
plt.show()


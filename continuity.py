import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fft
from grafic_ics import Snapshot

omega_b = 0.049
omega_m = 0.317
omega_l = 0.683
h = 0.673
H0 = 100.0 * h

# Load the overdensity ICs
ic = Snapshot('test_ics/', 7, 'deltab')
field = ic.load_box()

# Plot each k vector
# fig, axes = plt.subplots(1, 3, figsize=(18, 6))
# im0 = axes[0].imshow(kx[:, :, 0])
# cax0 = fig.colorbar(ax=axes[0], mappable=im0, fraction=0.046, pad=0.04)
# axes[0].set_title('kx[:, :, 0]')
# im1 = axes[1].imshow(ky[:, :, 0])
# cax1 = fig.colorbar(ax=axes[1], mappable=im1, fraction=0.046, pad=0.04)
# axes[1].set_title('ky[:, :, 0]')
# im2 = axes[2].imshow(kz[0, :, :])
# cax2 = fig.colorbar(ax=axes[2], mappable=im2, fraction=0.046, pad=0.04)
# axes[2].set_title('kz[0, :, :]')
# plt.savefig('k_vals.pdf')

# Plot the magnitude of k
# plt.figure()
# plt.imshow(k[:, :, 0])
# plt.colorbar()
# plt.title('k[:, :, 0]')
# plt.savefig('k_mod.pdf')
# plt.show()

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
        self.L = self.N * self.dx  # Box size

        # Calculate the relevant cosmological quantities
        self.c = Cosmology(self.a)

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
            self.kx = np.zeros(shape=(self.N,self.N,self.N))
            self.ky = np.zeros(shape=(self.N,self.N,self.N))
            self.kz = np.zeros(shape=(self.N,self.N,self.N))
            NN = ( self.N*fft.fftfreq(self.N, 1.) ).astype("i")
            # LC - here I have swapped the order kx and kz
            for i in NN:
                    self.kz[i,:,:] = i
                    self.ky[:,i,:] = i
                    self.kx[:,:,i] = i

            # LC - multiply all modes by scaling factor
            fac = (2.0*np.pi / self.L)
            self.kx *= fac
            self.ky *= fac
            self.kz *= fac
            self.ki = [self.kx, self.ky, self.kz]
            self.k = np.sqrt(self.kx**2. + self.ky**2. + self.kz**2.) * fac

        
    def unscaled_velocity(self):
        """
        Realise the (unscaled) velocity field in Fourier space. See Dodelson
        Eq. 9.18 for an expression; we factor out the time-dependent
        quantities here. They can be added at a later stage.

        """

        # If the FFT has an even number of samples, the most negative frequency 
        # mode must have the same value as the most positive frequency mode. 
        # However, when multiplying by 'i', allowing this mode to have a 
        # non-zero real part makes it impossible to satisfy the reality 
        # conditions. As such, we can set the whole mode to be zero, make sure 
        # that it's pure imaginary, or use an odd number of samples. Different 
        # ways of dealing with this could change the answer!
        if self.N % 2 == 0: # Even no. samples
                # Set highest (negative) freq. to zero
                self.kx[self.kx == np.min(self.kx)] = 0.0
                self.ky[self.ky == np.min(self.ky)] = 0.0
                self.kz[self.kz == np.min(self.kz)] = 0.0

        # Get squared k-vector in k-space (and factor in scaling from kx, ky, kz)
        k2 = self.k ** 2.0

        # Calculate components of A (the unscaled velocity)
        Ax = 1j * self.delta_k * self.kx / k2  # * (2.0*np.pi/self.L) / k2
        Ay = 1j * self.delta_k * self.ky / k2  # * (2.0*np.pi/self.L) / k2
        Az = 1j * self.delta_k * self.kz / k2  # * (2.0*np.pi/self.L) / k2
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
        # # Calculate the unscaled velocity
        # self.velocity_k = self.unscaled_velocity()

        # Prefactor for the scaling
        p1 = self.a * self.c.H_a * omega_m**0.6 * h
        # p1 = self.c.H_a * omega_m**0.6
        
        # # Calculate the unscaled_velocity
        # A = [(1j * self.delta_k * _k)/self.k**2 for _k in self.ki]
        # A = np.nan_to_num(A)
        # print(A.shape)

        # Now calculate the scaled velocity
        self.velocity_ks = [p1 * v for v in self.velocity_k]


    def scaled_velocity_full(self):
        """
        Calculate the scaled velocity using the full expression for the
        continuity equation as in Iliev at al. (2007)
        """
        # # Calculate the unscaled velocity
        # self.velocity_k = self.unscaled_velocity()

        # p1 and p2 are part of d/dt(D) / D
        p1 = - (3.0 * H0 * omega_m) / (self.c.E_a * self.a**3.0)
        p2 = 1 - (5.0 * self.a)/(3.0 * self.c.D_a)
        
        # Calculate the unscaled_velocity
        # A = [(1j * self.delta_k * _k)/self.k**2 for _k in self.ki]

        # Now calculate the scaled velocity
        # self.velocity_ks = [p1 * p2 * np.nan_to_num(_A) for _A in A]
        self.velocity_ks = [p1 * p2 * v for v in self.velocity_k]


class Cosmology:
    def __init__(self, a):
        # Calculate all the useful cosmological quantities at scale
        # factor a
        self.a = a
        self.H_a = self.H(a)
        self.E_a = self.E(a)
        self.D_a = self.D(a)


    def H(self, a):
        return H0 * np.sqrt(omega_m/self.a**3 + omega_l)


    def E(self, a):
        return self.H(a) / H0


    def D(self, a):
        from scipy.integrate import quad
        i = quad(lambda _a: 1 / (_a * self.E(_a)**3.0), 0.0, a)[0]

        # Normalise?
        # i0 = quad(lambda _a: 1 / (_a * self.E(_a)**3.0), 0.0, 1.0)[0]
        # i /= i0

        return (5.0 * omega_m * self.E(a) * i)/2.0 

    
    def a_to_z(self, a):
        return (1.0/a - 1.0)


    def z_to_a(self, z):
        return 1.0/(1.0 + z)

# Use the approximate or full continuty equation?
approx = True

# For plots
fn = 'full'
if approx: fn = 'approx'

c = Continuity(ic, z=200, approx=approx)
v_k = c.velocity_ks
v = [np.fft.ifftn(v_i).real for v_i in c.velocity_ks]
v_m = [vx_music, vy_music, vz_music]

vx_r = v[0][:, :, 0].ravel()
vx_mr = vx_music[:, :, 0].ravel()

from scipy.optimize import curve_fit

def y(x, m, c):
    return m*x + c

popt, pcov = curve_fit(y, vx_r, vx_mr)

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))
im0 = axes[0].imshow(v[0][:, :, 0])
plt.colorbar(ax=axes[0], mappable=im0, fraction=0.046, pad=0.04)
axes[0].set_title('vx[:, :, 0]')
im1 = axes[1].imshow(vx_music[:, :, 0])
plt.colorbar(ax=axes[1], mappable=im1, fraction=0.046, pad=0.04)
axes[1].set_title(r'vx\_music[:, :, 0]')
im2 = axes[2].imshow(v[0][:, :, 0]-vx_music[:, :, 0], cmap='jet')
plt.colorbar(ax=axes[2], mappable=im2, fraction=0.046, pad=0.04)
axes[2].set_title('vx-vx\_music')
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

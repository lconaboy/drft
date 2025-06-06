"""
Utility functions to include drift velocity in grafIC ics by computing/convolving density
power spectrum k dependent bias. Contains routines to run CICsASS
"""
import sys
import numpy as np
from py_vbc import run_pyvbc

def fft_sample_spacing(N, boxsize):
    from cosmology import _fft_sample_spacing
    return _fft_sample_spacing(N, boxsize)


def fft_sample_spacing_components(N):
    from cosmology import _fft_sample_spacing_components
    return _fft_sample_spacing_components(N)

def vbc_rms(vbc_field):
    '''
    Computes the rms vbc in the box.  Do not use for interpolated v_bc
    -- gives too much credance to the unphysical values at the edges
    of the box.  Use vbc_med instead.
    '''
    rms = np.sqrt(np.mean(vbc_field ** 2))
    return rms


def vbc_med(vbc_field):
    """
    Computes the median vbc in the box.  Better for interpolated
    (particularly non-periodic) interpolated values.
    """
    med = np.median(vbc_field)

    return med


def vbc_patch_dist(vbc_field, origin):
    import matplotlib.pyplot as plt
    
    nbins = (vbc_field.shape[0] * vbc_field.shape[1] *
             vbc_field.shape[2]) // 4000  # keep ~4000 elements in each bin

    rms = np.sqrt(np.mean(vbc_field ** 2))
    rms_cut = np.sqrt(np.mean(vbc_field[vbc_field<50] ** 2))
    med = np.median(vbc_field)
    # rmm = np.sqrt(np.median(vbc_field**2))

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.hist(vbc_field.ravel(), bins=nbins, histtype='step')
    ax.axvline(rms, c='r', ls='solid', label='rms')
    ax.axvline(med, c='m', ls='dotted', label='median')
    # ax.axvline(rmm, c='g', ls='dotted', label='rmm')
    ax.axvline(rms_cut, c='saddlebrown', ls='dashed', label='rms $<$ 50 km s$^{{-1}}$ cut')
    ax.set_xlabel('v$_{{\\sf bc}}$ (km s$^{{-1}}$)')
    ax.set_ylabel('N')
    ax.legend()
    fig.savefig('vbc_dist_{0}_{1}_{2}.pdf'.format(
        origin[0], origin[1], origin[2]), bbox_inches='tight')
    
    
def msg(rank, s, verbose=True):
    if verbose:
        print('[rank {0:03d}]: {1}'.format(rank, s), flush=True)
    else:
        return

def apply_density_bias(ics, k_bias, b, N, delta_x=None):
    ''' Apply a bias to the realisations power spectrum, and recompute the 3D field.
    Parameters:
        b (array): bias to deconvolve with the delta_x field, such that:
        delta_x = ifft(delta_k/b)
    '''
    import scipy.fftpack as fft
    import scipy.interpolate as si

    if delta_x is None:
        delta_x = ics

    # The actual shape of the delta_x region
    shape = delta_x.shape
    # The shape of the symmetric box
    shape0 = (shape[0], shape[0], shape[0])

    boxsize = float(ics.boxsize) * \
        (float(N) / float(ics.N))  # this boxsize is the reduced one

    # print "boxsize = ", boxsize, delta_x.shape[0]

    # k = None
    if boxsize != ics.boxsize:
        # Resample k as we may be using a subregion
        # k = fft_sample_spacing(delta_x.shape[0], boxsize).flatten()
        pass
    else:
        # k = ics.k.flatten()
        print('ics.boxsize:', ics.boxsize)
        print('ics.N:', ics.N)
        print('N:', N)
        print('boxsize:', boxsize)
    # Resample k as we may be using a subregion
    k = fft_sample_spacing(delta_x.shape[0], boxsize).flatten()
    
    k[k == 0.] = (2. * np.pi) / boxsize

    # Interpolate/extrapolate the bias to the 3D grid
    def log_interp1d(xx, yy, kind='linear'):
        logx = np.log10(xx)
        logy = np.log10(yy)
        lin_interp = si.InterpolatedUnivariateSpline(logx, logy)
        log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
        return log_interp

    f = log_interp1d(k_bias, b)
    b = f(k)


    # print('LC testing')
    # print(np.any(np.isnan(b)))
    # print('LC testing before fft')
    # print(np.any(np.isnan(delta_x)))
    delta_k = fft.fftn(delta_x)
    # print('LC testing after fft')
    # print(np.any(np.isnan(delta_k)))

    # print('max before mult', delta_k.max())
    # print('max b.reshape(delta_k.shape)', b.reshape(delta_k.shape).max())
    # if b.reshape(delta_k.shape).max() > 1:
    #     np.savetxt('k.dat', k)
    #     np.savetxt('b.dat', b)
    
    # Apply the bias
    # delta_k *= np.sqrt(b.reshape(delta_k.shape))
    delta_k = delta_k * np.sqrt(b.reshape(delta_k.shape))
    # print('max after mult', delta_k.max())

    
    # print('LC testing after mult')
    print(np.any(np.isnan(delta_k)))

    
    # Inverse FFT to compute the realisation

    delta_x = fft.ifftn(delta_k).real.reshape(shape)
    
    return delta_x


def apply_velocity_bias(ics, k_bias, b, N, vel=None):
    """
    Calculate the velocity bias. Here, the bias has to be vectorised since we 
    have vx, vy and vz.

    :param ics: 
    :param k_bias: 
    :param b: 
    :param N: 
    :param vel: 
    :returns: 
    :rtype: 

    """
    
    return


def compute_bias(ics, vbc, zstart=1000, kmin=0.1, kmax=10000, n=100, delta=False):
    """
    Computes the bias to /both/ density and velocity fields.  Assumes
    v_bc is constant at z=zstart.

    :param ics: (Snapshot) Snapshot object containing grafic ICs

    :param vbc: (array) Array containing the v_bc (i.e. |v_b - v_c|)
        field.

    :param zstart: float, redshift of recombination
    :param kmin: float, minimum k-value to solve for in py_vbc
        (Mpc^-1)
    :param kmax: float, maximum k-value to solve for in py_vbc
        (Mpc^-1)
    :param n: int, number of k-values in total if positive, ~number per
              log10(k) if negative
    """
    # Compute size of grid and boxsize
    N = vbc.shape[0]
    boxsize = float(ics.boxsize) * (float(N) / float(ics.N))

    # v_bc redshifts away, so calculate the v_bc at z=zstart
    z = ics.z
    # zstart=1000
    
    # LC - switched to the median instead, as the rms was giving too
    # extreme values, particularly near the edge in non-periodic
    # interpolated vbc fields

    # Also LC: now that I'm actually properly interpolating with yt
    # I'll switch back to using the rms, which gives pretty much the
    # same result as the median now
    
    # rms = vbc_med(vbc)
    rms = vbc_rms(vbc)
    rms_recom = rms * (1. + zstart) / (1. + z)

    print(f'            v_bc rms {rms:.2f} km/s recom {rms_recom:.2f} km/s')
    
    # Calculate how many samples we need for the given per log10(k)
    if (n < 0):
        dlk = np.log10(kmax) - np.log10(kmin)
        n = max([np.int(np.ceil(np.abs(n) * dlk)), 100])
        # print('compute_bias')
        # print('kmin', kmin, 'kmax', kmax, 'dlk', dlk, 'n', n)
        # sys.exit(0)
    
    # Boxsize doesn't make a difference when calculating the power
    # spectra using py_vbc. The power spectrum tuple contains (p_c, p_b, p_vc,
    # p_vb) and k is in units of Mpc^-1.
    k, ps_vbc0 = run_pyvbc(vbc=0.0, zstart=zstart, zend=z, dz=3, kmin=kmin,
                           kmax=kmax, n=n, delta=delta)
    k, ps_vbcrecom = run_pyvbc(vbc=rms_recom, zstart=zstart, zend=z, dz=3, kmin=kmin,
                               kmax=kmax, n=n, delta=delta)

    # Calculate the biases
    b_c = ps_vbcrecom[0] / ps_vbc0[0]
    b_b = ps_vbcrecom[1] / ps_vbc0[1]
    b_vc = ps_vbcrecom[2] / ps_vbc0[2]
    b_vb = ps_vbcrecom[3] / ps_vbc0[3]

    return k, b_c, b_b, b_vc, b_vb


def cube_positions(ics, n, N=None):
    cubes = []
    if N is None:
        # Beware that ics.n is a tuple and ics.N is an int!
        # N = ics.N
        N = ics.n

    # if (N % n != 0):
    #     raise Exception(
    #         "Cannot fit %d cubes into grid with size %d" % (n, N))

    if ~np.all(np.mod(n, N)):
        raise Exception('Cannot fit {0} cubes in grid with size {1}'.format(n, N))
    
    dx_cells = N / n

    for i in range(int(n[0])):
        cen_i = dx_cells[0] * (i + 0.5)

        for j in range(int(n[1])):
            cen_j = dx_cells[1] * (j + 0.5)

            for k in range(int(n[2])):
                cen_k = dx_cells[2] * (k + 0.5)

                cubes.append([cen_i, cen_j, cen_k])

    return cubes, dx_cells


def divisors(number, mode='print'):
    n = 1
    while(n < number):
        if(number % n == 0):
            if mode is 'print':
                print(n)
            elif mode is 'yield':
                yield n
        else:
            pass
        n += 1


# Python version of bash which
def which(program):
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def clean(level):
    """Clean up after execution has finished. Anything that needs to be
    got rid of can be done here.
    """
    import shutil

    # Remove patches/
    shutil.rmtree('./patches/level_{0:03d}'.format(level))
    


# def vbc_ps_fname(rms, z, boxsize):
#     import os
#     cwd = os.getcwd()
#     if not os.path.isdir("%s/vbc_TFs_out" % cwd):
#         os.mkdir("%s/vbc_TFs_out" % cwd)
#     return '%s/vbc_TFs_out/vbc_%f_z%f_B%1.2f.dat' % (cwd, rms, z, boxsize)


# def run_cicsass_lc(boxsize, z, rms_vbc_z1000, N=256):
#     import subprocess, os, gc

#     exe = which('transfer.x')

#     if exe is None:
#         raise Exception("Unable to locate transfer.x executable")

#     # Example execution for RMS vbc=30km/s @ z=1000.:
#     # ./transfer.x -B0.2 -N128 -V30 -Z100 -D3 -SinitSB_transfer_out

#     CICsASS_home = os.getenv("CICSASS_HOME")
#     if CICsASS_home is None:
#         raise Exception("Env var CICSASS_HOME not set")

#     # Run with N=256
#     # CICsASS_home = "/lustre/scratch/astro/ds381/CICsASS/matt/Dropbox/CICASS/vbc_transfer/"
#     cmd = 'cd %s && %s -B%1.2f -N%d -V%f -Z%f -D3 -SinitSB_transfer_out' % (
#         CICsASS_home, exe, boxsize, N, rms_vbc_z1000, z)
#     # print 'Running:\n%s' % cmd

#     gc.collect() # Collect garbage

#     # Run CICsASS and wait for output, check_output is blocking and
#     # will return an Exception if cmd fails
#     output = subprocess.check_output(cmd, shell=True)
#     output = output.decode("ascii")
#     output = output.splitlines()
    
#     vals = np.zeros(shape=(64, 4))

#     # This is slow but perhaps unavoidable
#     for i in range(64):
#         vals[i, :] = output[i].split()

#     # Transpose to match original code
#     vals = np.transpose(vals)
#     # and unpack
#     vals = [vals[0, :], vals[1, :], vals[2, :]]

#     gc.collect() # Collect garbage
    
#     return vals 


# def run_cicsass(boxsize, z, rms_vbc_z1000, out_fname, N=256):
#     import subprocess, os
 
#     exe = which('transfer.x')

#     if exe is None:
#         raise Exception("Unable to locate transfer.x executable")

#     # Example execution for RMS vbc=30km/s @ z=1000.:
#     # ./transfer.x -B0.2 -N128 -V30 -Z100 -D3 -SinitSB_transfer_out

#     CICsASS_home = os.getenv("CICSASS_HOME")
#     if CICsASS_home is None:
#         raise Exception("Env var CICSASS_HOME not set")

#     # Run with N=256
#     # CICsASS_home = "/lustre/scratch/astro/ds381/CICsASS/matt/Dropbox/CICASS/vbc_transfer/"
#     cmd = 'cd %s && %s -B%1.2f -N%d -V%f -Z%f -D3 -Splanck2018_transfer_out > %s' % (
#         CICsASS_home, exe, boxsize, N, rms_vbc_z1000, z, out_fname)
#     # print 'Running:\n%s' % cmd
#     # Run CICsASS and wait for output
#     code = subprocess.check_call(cmd, shell=True)
#     if code != 0:
#         raise Exception("CICsASS returned non-zero exit code: %d", code)
#     return code


# def compute_velocity_bias(ics, vbc):
#     import os, time
#     # print 'AVERAGE INSTEAD OF RMS'
#     # Init fields
#     if vbc is None:
#         vbc = ics['vbc']

#     # Compute size of grid and boxsize
#     N = vbc.shape[0]
#     boxsize = float(ics.boxsize) * \
#         (float(N) / float(ics.N))

#     # Compute vbc @ z=1000
#     # vbc_norm = ics.vbc_rms_norm(vbc=vbc)
#     # vbc_rms = vbc_norm * (1001.)  # vbc_rms prop (1 + z)
#     # Compute vbc @ z=1000
#     z = ics.z
#     rms = vbc_rms(vbc)
#     rms_recom = rms * (1001./(z + 1.0))

#     # Check for PS and run CICsASS if necessary
#     fname_vbc0 = vbc_ps_fname(0., z, boxsize)
#     if os.path.isfile(fname_vbc0) is False:
#         exit_code = run_cicsass(boxsize, z, 0., fname_vbc0)

#     fname_vbcrecom = vbc_ps_fname(rms_recom, z, boxsize)
#     if os.path.isfile(fname_vbcrecom) is False:
#         exit_code = run_cicsass(boxsize, z, rms_recom, fname_vbcrecom)

#     # Load the power spectra and compute the bias
#     # LC - might be too quick for CICASS, check for empty files
#     ps_vbc0 = np.loadtxt(fname_vbc0, unpack=True)
#     ps_vbcrecom = np.loadtxt(fname_vbcrecom, unpack=True)
#     count = 0
#     while ((len(ps_vbc0) == 0) or (len(ps_vbcrecom) == 0)):
#         count += 1
#         if count > 10:
#             raise Exception("Reached sleep limit. File still empty.")
#             print("Caught exception (fname_vbc0): {0}".format(fname_vbc0))
#             print("Caught exception (fname_vbcrecom): {0}".format(fname_vbcrecom))
#         time.sleep(5)
#         ps_vbc0 = np.loadtxt(fname_vbc0, unpack=True)
#         ps_vbcrecom = np.loadtxt(fname_vbcrecom, unpack=True)
    
#     # Should have same lenghts if finished writing
#     count = 0
#     try:
#         while len(ps_vbcrecom[1]) != len(ps_vbc0[1]):
#             count += 1
#             if count > 10:
#                 raise Exception("Reached sleep limit. Filesizes still differ.")
#             time.sleep(5)
#             ps_vbc0 = np.loadtxt(fname_vbc0, unpack=True)
#             ps_vbcrecom = np.loadtxt(fname_vbcrecom, unpack=True)
#     except Exception as e:
#         print("Caught exception (fname_vbc0): {0}".format(fname_vbc0))
#         print("Caught exception (fname_vbcrecom): {0}".format(fname_vbcrecom))

#     cosmo = ics.cosmo

#     import cosmology
#     vdeltab0 = cosmology.linear_velocity_ps(
#         ps_vbc0[0], np.sqrt(ps_vbc0[2]), **cosmo)
#     vdeltab = cosmology.linear_velocity_ps(
#         ps_vbcrecom[0], np.sqrt(ps_vbcrecom[2]), **cosmo)

#     vdeltac0 = cosmology.linear_velocity_ps(
#         ps_vbc0[0], np.sqrt(ps_vbc0[1]), **cosmo)
#     vdeltac = cosmology.linear_velocity_ps(
#         ps_vbcrecom[0], np.sqrt(ps_vbcrecom[1]), **cosmo)

#     #CDM bias
#     b_cdm = vdeltac / vdeltac0
#     # Baryon bias/p/scratch/chpo22/hpo22i/bd/cicass/vbc_transfer/vbc_TFs_out/vbc_22.435140_z200.000005_B3.52.dat
#     b_b = vdeltab / vdeltab0
#     # Wavenumber
#     k_bias = ps_vbcrecom[0] / ics.cosmo["h"]

#     return k_bias, b_cdm, b_b


# def compute_velocity_bias_lc(ics, vbc):
#     import os, time
#     # print 'AVERAGE INSTEAD OF RMS'
#     # Init fields
#     if vbc is None:
#         vbc = ics

#     # Compute size of grid and boxsize
#     N = vbc.shape[0]
#     boxsize = float(ics.boxsize) * \
#         (float(N) / float(ics.N))

#     # Compute vbc @ z=1000
#     # vbc_norm = ics.vbc_rms_norm(vbc=vbc)
#     # vbc_rms = vbc_norm * (1001.)  # vbc_rms prop (1 + z)
#     # Compute vbc @ z=1000
#     z = ics.z
#     zstart=1000
#     rms = vbc_rms(vbc)
#     rms_recom = rms * (1001./(1.0 + z))

#     ps_vbc0 = run_cicsass_lc(boxsize, z, 0.)
#     ps_vbcrecom = run_cicsass_lc(boxsize, z, rms_recom)

#     # Boxsize doesn't make a difference when calculating the power spectra
#     # ps_vbc0 = run_pyvbc(vbc=0.0, zstart=zstart, zend=z, dz=3)
#     # ps_vbcrecom = run_pyvbc(vbc=rms_recom, zstart=zstart, zend=z, dz=3)

#     cosmo = ics.cosmo

#     import cosmology
#     vdeltab0 = cosmology.linear_velocity_ps(
#         ps_vbc0[0], np.sqrt(ps_vbc0[2]), **cosmo)
#     vdeltab = cosmology.linear_velocity_ps(
#         ps_vbcrecom[0], np.sqrt(ps_vbcrecom[2]), **cosmo)

#     vdeltac0 = cosmology.linear_velocity_ps(
#         ps_vbc0[0], np.sqrt(ps_vbc0[1]), **cosmo)
#     vdeltac = cosmology.linear_velocity_ps(
#         ps_vbcrecom[0], np.sqrt(ps_vbcrecom[1]), **cosmo)

#     #CDM bias
#     b_cdm = vdeltac / vdeltac0
#     # Baryon bias/p/scratch/chpo22/hpo22i/bd/cicass/vbc_transfer/vbc_TFs_out/vbc_22.435140_z200.000005_B3.52.dat
#     b_b = vdeltab / vdeltab0
#     # Wavenumber
#     k_bias = ps_vbcrecom[0] / ics.cosmo["h"]  # "h Mpc**-1"

#     return k_bias, b_cdm, b_b


# def compute_cicsass(ics, vbc):
#     """Function used to calculate all the cicass power spectra before
#     doing anything else. Not very efficient, but might be necessary."""
#     import os, time
   
#     # Compute size of grid and boxsize (for this patch)
#     N = vbc.shape[0]
#     boxsize = ics.boxsize * (float(N) / float(ics.N))  # "Mpc a h**-1" 

#     # Compute vbc @ z=1000
#     z = ics.z
#     rms = vbc_rms(vbc)
#     rms_recom = rms * (1001./(1.0 + z))

#         # Check for PS and run CICsASS if needed
#     fname_vbc0 = vbc_ps_fname(0., z, boxsize)
#     if not os.path.isfile(fname_vbc0):
#         exit_code = run_cicsass(boxsize, z, 0., fname_vbc0)

#     fname_vbcrecom = vbc_ps_fname(rms_recom, z, boxsize)
#     if not os.path.isfile(fname_vbcrecom):
#         exit_code = run_cicsass(boxsize, z, rms_recom, fname_vbcrecom)


# def compute_bias(ics, vbc):
#     """ Calculate the bias to the density power spectrum assuming
#     COHERENT vbc at z=1000. """
#     import os, time
   
#     # Compute size of grid and boxsize (for this patch)
#     N = vbc.shape[0]
#     boxsize = ics.boxsize * (float(N) / float(ics.N))  # "Mpc a h**-1"

#     # Compute vbc @ z=1000
#     z = ics.z
#     rms = vbc_rms(vbc)
#     rms_recom = rms * (1001./(1.0 + z))

#     # Check for PS and run CICsASS if needed
#     fname_vbc0 = vbc_ps_fname(0., z, boxsize)
#     if not os.path.isfile(fname_vbc0):
#         exit_code = run_cicsass(boxsize, z, 0., fname_vbc0)

#     fname_vbcrecom = vbc_ps_fname(rms_recom, z, boxsize)
#     if not os.path.isfile(fname_vbcrecom):
#         exit_code = run_cicsass(boxsize, z, rms_recom, fname_vbcrecom)

#     # Load the power spectra and compute the bias
#     # LC - might be too quick for CICASS, check for empty files
#     ps_vbc0 = np.loadtxt(fname_vbc0, unpack=True)
#     ps_vbcrecom = np.loadtxt(fname_vbcrecom, unpack=True)
#     count = 0
#     while ((len(ps_vbc0) == 0) or (len(ps_vbcrecom) == 0)):
#         count += 1
#         if count > 10:
#             raise Exception("Reached sleep limit. File still empty.")
#             print("Caught exception (fname_vbc0): {0}".format(fname_vbc0))
#             print("Caught exception (fname_vbcrecom): {0}".format(fname_vbcrecom))
#         time.sleep(5)
#         ps_vbc0 = np.loadtxt(fname_vbc0, unpack=True)
#         ps_vbcrecom = np.loadtxt(fname_vbcrecom, unpack=True)

#     # Should have same lenghts if finished writing
#     count = 0
#     try:
#         while len(ps_vbcrecom[1]) != len(ps_vbc0[1]):
#             count += 1
#             if count > 10:
#                 raise Exception("Reached sleep limit. Filesizes still differ")
#             time.sleep(5)
#             ps_vbc0 = np.loadtxt(fname_vbc0, unpack=True)
#             ps_vbcrecom = np.loadtxt(fname_vbcrecom, unpack=True)
#     except Exception as e:
#         print("Caught exception (fname_vbc0): {0}".format(fname_vbc0))
#         print("Caught exception (fname_vbcrecom): {0}".format(fname_vbcrecom))

#     #CDM bias
#     b_cdm = ps_vbcrecom[1] / ps_vbc0[1]
#     # Baryon bias
#     b_b = ps_vbcrecom[2] / ps_vbc0[2]
#     # Wavenumber
#     k_bias = ps_vbcrecom[0] / ics.cosmo["h"]

#     return k_bias, b_cdm, b_b


# def compute_bias_lc(ics, vbc):
#     """ Calculate the bias to the density power spectrum assuming
#     COHERENT vbc at z=1000. """
#     import os, time
   
#     # Compute size of grid and boxsize (for this patch)
#     N = vbc.shape[0]
#     boxsize = ics.boxsize * (float(N) / float(ics.N))  # "Mpc a h**-1"

#     # Compute vbc @ z=1000
#     z = ics.z
#     zstart = 1000
#     rms = vbc_rms(vbc)
#     rms_recom = rms * (1001./(1.0 + z))

#     ps_vbc0 = run_cicsass_lc(boxsize, z, 0.)
#     ps_vbcrecom = run_cicsass_lc(boxsize, z, rms_recom)
    
#     # Boxsize doesn't make a difference when calculating the power
#     # spectra using py_vbc
#     # ps_vbc0 = run_pyvbc(vbc=0.0, zstart=zstart, zend=z, dz=3)
#     # ps_vbcrecom = run_pyvbc(vbc=rms_recom, zstart=zstart, zend=z, dz=3)

#     #CDM bias
#     b_cdm = ps_vbcrecom[1] / ps_vbc0[1]
#     # Baryon bias
#     b_b = ps_vbcrecom[2] / ps_vbc0[2]
#     # Wavenumber
#     k_bias = ps_vbcrecom[0] / ics.cosmo["h"]# "h Mpc**-1"
    
#     return k_bias, b_cdm, b_b

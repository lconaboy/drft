"""
Utility functions to include drift velocity in grafIC ics by computing/convolving density
power spectrum k dependent bias. Contains routines to run CICsASS
"""
import sys
import warnings
from dataclasses import dataclass

import numpy as np
from py_vbc import RuntimeConfig, load_config, run_pyvbc


AHN_MU_DMU = 0.05
AHN_N_MU = int(round(1.0 / AHN_MU_DMU)) + 1


@dataclass
class BiasGrid:
    """Wavenumber bias, optionally resolved by direction cosine."""

    k: np.ndarray
    b_c: np.ndarray
    b_b: np.ndarray
    b_vc: np.ndarray
    b_vb: np.ndarray
    mu: np.ndarray | None = None
    direction: np.ndarray | None = None

    def __iter__(self):
        """Keep the historical five-value unpacking API working."""
        yield self.k
        yield self.b_c
        yield self.b_b
        yield self.b_vc
        yield self.b_vb


def _directional_bias_table(
    k_values,
    mu_values,
    ps_vbc0,
    rms_recom,
    zstart,
    zend,
    dz,
    config,
    delta,
):
    """Evaluate all four bias components on a ``(k, mu)`` grid.

    ``mu_values`` contains the non-negative half of the directional grid.
    The negative half is recovered by symmetry when the bias is applied.
    """
    mu_values = np.asarray(mu_values, dtype=float)
    if mu_values.ndim != 1 or len(mu_values) < 2:
        raise ValueError("directional mu grid must contain at least two samples")
    if not np.all(np.isfinite(mu_values)) or np.any(np.diff(mu_values) <= 0):
        raise ValueError("directional mu grid must be finite and strictly increasing")
    if mu_values[0] < -1.0e-12 or mu_values[-1] > 1.0 + 1.0e-12:
        raise ValueError("directional mu values must lie in [0, 1]")

    n_k = len(k_values)
    n_mu = len(mu_values)
    k_grid = np.repeat(k_values, n_mu)
    mu_grid = np.tile(mu_values, n_k)
    _, ps_vbcrecom = run_pyvbc(
        vbc=rms_recom,
        zstart=zstart,
        zend=zend,
        dz=dz,
        config_file=config,
        k=k_grid,
        mu=mu_grid,
        delta=delta,
    )
    shape = (n_k, n_mu)
    return tuple(
        ps_vbcrecom[index].reshape(shape) / ps_vbc0[index][:, None]
        for index in range(4)
    )


def _directional_bias_error(
    k_values,
    mu_current,
    mu_previous,
    current,
    previous,
    mu_atol,
):
    """Return the maximum relative change caused by refining the mu grid."""
    import scipy.interpolate as si

    log_k = np.log10(k_values)
    kx = min(3, len(k_values) - 1)
    ky = min(3, len(mu_previous) - 1)
    if kx < 1 or ky < 1:
        raise ValueError("directional bias needs at least two k and mu samples")

    query_k, query_mu = np.meshgrid(log_k, mu_current, indexing="ij")
    error = 0.0
    for current_values, previous_values in zip(current, previous):
        interpolator = si.RectBivariateSpline(
            log_k,
            mu_previous,
            previous_values,
            kx=kx,
            ky=ky,
            s=0,
        )
        predicted = interpolator(query_k, query_mu, grid=False).reshape(
            current_values.shape
        )
        scale = np.maximum(np.abs(current_values), mu_atol)
        error = max(
            error,
            float(np.max(np.abs(current_values - predicted) / scale)),
        )
    return error


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


def vbc_vector(velb, velc):
    """Return the three-component relative velocity field ``v_b - v_c``."""
    components = []
    for axis in range(3):
        components.append(np.asarray(velb[axis]) - np.asarray(velc[axis]))
    return np.stack(components, axis=-1)


def mean_vbc_direction(vbc_field, trim=0):
    """Return the normalized mean relative-velocity direction."""
    field = np.asarray(vbc_field, dtype=float)
    if trim:
        field = field[tuple(slice(trim, -trim) for _ in range(3))]
    direction = np.mean(field, axis=(0, 1, 2))
    magnitude = np.linalg.norm(direction)
    if not np.isfinite(magnitude) or magnitude == 0:
        raise ValueError("relative-velocity direction is zero or non-finite")
    return direction / magnitude


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

def apply_density_bias(ics, bias, component=None, N=None, delta_x=None):
    """Apply a scalar or directional bias to a real-space field.

    ``bias`` is normally a :class:`BiasGrid`. The legacy calling convention
    ``apply_density_bias(ics, k_bias, b, N, delta_x=...)`` is also accepted.
    """
    import scipy.fftpack as fft
    import scipy.interpolate as si

    if isinstance(bias, BiasGrid):
        component_names = {
            "c": "b_c",
            "b": "b_b",
            "vc": "b_vc",
            "vb": "b_vb",
        }
        if component not in component_names:
            raise ValueError("unknown bias component")
        k_bias = np.asarray(bias.k, dtype=float)
        b = np.asarray(getattr(bias, component_names[component]), dtype=float)
        mu_bias = bias.mu
        direction = bias.direction
    else:
        k_bias = np.asarray(bias, dtype=float)
        b = np.asarray(component, dtype=float)
        mu_bias = None
        direction = None

    if delta_x is None:
        delta_x = ics
    delta_x = np.asarray(delta_x)
    shape = delta_x.shape
    if len(shape) != 3 or not (shape[0] == shape[1] == shape[2]):
        raise ValueError("directional bias currently requires a cubic patch")

    if N is None:
        N = shape[0]
    boxsize = float(ics.boxsize) * (float(N) / float(ics.N))

    kx, ky, kz = fft_sample_spacing_components(shape[0])
    fac = 2.0 * np.pi / boxsize
    k_vectors = np.stack((kx * fac, ky * fac, kz * fac), axis=-1)
    k_magnitude = np.linalg.norm(k_vectors, axis=-1)
    zero_mode = k_magnitude == 0.0
    k_flat = k_magnitude.copy().reshape(-1)
    k_flat[zero_mode.reshape(-1)] = 2.0 * np.pi / boxsize

    if mu_bias is None:
        interpolator = si.interp1d(
            np.log10(k_bias),
            np.log10(b),
            bounds_error=False,
            fill_value="extrapolate",
        )
        b_grid = np.power(10.0, interpolator(np.log10(k_flat)))
    else:
        if direction is None:
            raise ValueError("directional bias requires a velocity direction")
        direction = np.asarray(direction, dtype=float)
        if direction.shape != (3,):
            raise ValueError("velocity direction must have three components")
        direction_norm = np.linalg.norm(direction)
        if not np.isfinite(direction_norm) or direction_norm == 0:
            raise ValueError("velocity direction must be finite and non-zero")
        direction = direction / direction_norm

        # The power-spectrum bias is even in mu, so only the
        # non-negative half of the directional table is needed.
        mu_flat = np.abs(np.sum(k_vectors * direction, axis=-1).reshape(-1))
        mu_flat[zero_mode.reshape(-1)] = 0.0
        # Use the highest supported spline degree (up to cubic) along each axis.
        logk_degree = min(3, len(k_bias) - 1)
        mu_degree = min(3, len(mu_bias) - 1)
        if logk_degree < 1 or mu_degree < 1:
            raise ValueError("directional bias needs at least two k and mu samples")
        interpolator = si.RectBivariateSpline(
            np.log10(k_bias),
            np.asarray(mu_bias, dtype=float),
            b.reshape(len(k_bias), len(mu_bias)),
            kx=logk_degree,
            ky=mu_degree,
            s=0,
        )
        b_grid = interpolator(
            np.clip(np.log10(k_flat), np.log10(k_bias[0]), np.log10(k_bias[-1])),
            np.clip(mu_flat, np.min(mu_bias), np.max(mu_bias)),
            grid=False,
        )

    b_grid = b_grid.reshape(shape)
    if np.any(~np.isfinite(b_grid)) or np.any(b_grid < 0):
        raise ValueError("interpolated bias is invalid")

    delta_k = fft.fftn(delta_x)
    delta_k = delta_k * np.sqrt(b_grid)
    return fft.ifftn(delta_k).real.reshape(shape)

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


def compute_bias(
    ics,
    vbc,
    zstart=None,
    dz=None,
    kmin=0.1,
    kmax=10000,
    n=100,
    delta=False,
    config=None,
    vbc_direction=None,
    n_mu=None,
    mu_rtol=1.0e-2,
    mu_atol=1.0e-6,
    max_n_mu=81,
    zend=None,
):
    """Compute density and velocity bias, optionally resolved by direction.

    If ``vbc_direction`` is supplied, bias is evaluated on a ``(k, mu)``
    grid, where ``mu = k_hat . v_bc_hat``. By default the calculation starts
    from Ahn's 21 non-negative samples (``dmu=0.05``) and refines to a finer
    grid only if the recalculated bias is not converged to ``mu_rtol``. The
    negative half is reconstructed by symmetry. ``zstart`` and ``dz`` default
    to the runtime configuration's ``evolution`` values. For snapshot-based
    biasing, ``zend`` is always taken from ``ics.z``; a supplied ``zend`` is
    ignored and emits a warning. ``dz`` is the redshift offset used to
    differentiate transfer functions at ``zstart``. Set ``n_mu`` to a
    positive integer to request one fixed non-negative grid instead.
    Otherwise the historical isotropic bias is returned.
    """
    if config is None:
        raise ValueError("compute_bias requires a py_vbc runtime configuration")
    if not isinstance(config, RuntimeConfig):
        config = load_config(config)

    z_snapshot = ics.z
    evolution = config.evolution
    zstart = evolution.zstart if zstart is None else float(zstart)
    if zend is not None:
        warnings.warn(
            "compute_bias ignores the supplied zend and uses the grafic "
            f"snapshot redshift ics.z={z_snapshot}",
            RuntimeWarning,
            stacklevel=2,
        )
    # Snapshot-based biasing must always evolve to the IC redshift.
    zend = z_snapshot
    dz = evolution.dz if dz is None else float(dz)
    if dz <= 0:
        raise ValueError("dz must be greater than zero")

    rms = vbc_rms(vbc)
    rms_recom = rms * (1.0 + zstart) / (1.0 + z_snapshot)
    print(f'            v_bc rms {rms:.2f} km/s recom {rms_recom:.2f} km/s')

    if n < 0:
        dlk = np.log10(kmax) - np.log10(kmin)
        n = max([int(np.ceil(np.abs(n) * dlk)), 100])

    log_kmin = np.log10(kmin)
    log_kmax = np.log10(kmax)
    delta_log_k = (log_kmax - log_kmin) / float(n - 1.0)
    k_values = 10.0 ** (np.arange(n, dtype=float) * delta_log_k + log_kmin)
    mu_isotropic = config.cosmology.mu_isotropic

    if vbc_direction is None:
        _, ps_vbc0 = run_pyvbc(
            vbc=0.0,
            zstart=zstart,
            zend=zend,
            dz=dz,
            config_file=config,
            k=k_values,
            mu=mu_isotropic,
            delta=delta,
        )
        _, ps_vbcrecom = run_pyvbc(
            vbc=rms_recom,
            zstart=zstart,
            zend=zend,
            dz=dz,
            config_file=config,
            k=k_values,
            mu=mu_isotropic,
            delta=delta,
        )
        return BiasGrid(
            k=k_values,
            b_c=ps_vbcrecom[0] / ps_vbc0[0],
            b_b=ps_vbcrecom[1] / ps_vbc0[1],
            b_vc=ps_vbcrecom[2] / ps_vbc0[2],
            b_vb=ps_vbcrecom[3] / ps_vbc0[3],
        )

    direction = np.asarray(vbc_direction, dtype=float)
    if direction.shape != (3,):
        raise ValueError("vbc_direction must have three components")
    direction_norm = np.linalg.norm(direction)
    if not np.isfinite(direction_norm) or direction_norm == 0:
        raise ValueError("vbc_direction must be finite and non-zero")
    direction = direction / direction_norm
    if mu_rtol <= 0 or mu_atol <= 0:
        raise ValueError("mu_rtol and mu_atol must be positive")
    if max_n_mu < AHN_N_MU:
        raise ValueError(f"max_n_mu must be at least {AHN_N_MU}")
    if len(k_values) < 2:
        raise ValueError("directional bias needs at least two k samples")

    # v_bc=0 is mu-independent, so use a fixed reference value here.
    _, ps_vbc0 = run_pyvbc(
        vbc=0.0,
        zstart=zstart,
        zend=zend,
        dz=dz,
        config_file=config,
        k=k_values,
        mu=1.0,
        delta=delta,
    )

    if n_mu is not None:
        n_mu = int(n_mu)
        if n_mu < 2:
            raise ValueError("n_mu must be at least 2 for directional interpolation")
        mu_values = np.linspace(0.0, 1.0, n_mu)
        bias_tables = _directional_bias_table(
            k_values,
            mu_values,
            ps_vbc0,
            rms_recom,
            zstart,
            zend,
            dz,
            config,
            delta,
        )
    else:
        candidate = AHN_N_MU
        previous = None
        previous_mu = None
        bias_tables = None
        while True:
            mu_values = np.linspace(0.0, 1.0, candidate)
            bias_tables = _directional_bias_table(
                k_values,
                mu_values,
                ps_vbc0,
                rms_recom,
                zstart,
                zend,
                dz,
                config,
                delta,
            )
            if previous is not None:
                error = _directional_bias_error(
                    k_values,
                    mu_values,
                    previous_mu,
                    bias_tables,
                    previous,
                    mu_atol,
                )
                if error <= mu_rtol:
                    selected_n_mu = len(previous_mu)
                    print(
                        "            v_bc mu grid converged: "
                        f"n_mu={selected_n_mu}, "
                        f"dmu={1.0 / (selected_n_mu - 1):.3g}, "
                        f"validated against n_mu={candidate}, "
                        f"max relative change={error:.3e}"
                    )
                    mu_values = previous_mu
                    bias_tables = previous
                    break
            if candidate >= max_n_mu:
                if previous is None:
                    message = "no refinement comparison was possible"
                else:
                    message = f"max relative change={error:.3e}"
                raise RuntimeError(
                    "directional bias did not converge before max_n_mu="
                    f"{max_n_mu} ({message})"
                )
            previous = bias_tables
            previous_mu = mu_values
            candidate = 2 * candidate - 1

    return BiasGrid(
        k=k_values,
        b_c=bias_tables[0],
        b_b=bias_tables[1],
        b_vc=bias_tables[2],
        b_vb=bias_tables[3],
        mu=mu_values,
        direction=direction,
    )

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
            if mode == 'print':
                print(n)
            elif mode == 'yield':
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
#     # ps_vbc0 = run_pyvbc(vbc=0.0, zstart=zstart, zend=z, dz=3,
#     #                    config_file='params.yaml')
#     # ps_vbcrecom = run_pyvbc(vbc=rms_recom, zstart=zstart, zend=z, dz=3,
#     #                        config_file='params.yaml')

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
#     # ps_vbc0 = run_pyvbc(vbc=0.0, zstart=zstart, zend=z, dz=3,
#     #                    config_file='params.yaml')
#     # ps_vbcrecom = run_pyvbc(vbc=rms_recom, zstart=zstart, zend=z, dz=3,
#     #                        config_file='params.yaml')

#     #CDM bias
#     b_cdm = ps_vbcrecom[1] / ps_vbc0[1]
#     # Baryon bias
#     b_b = ps_vbcrecom[2] / ps_vbc0[2]
#     # Wavenumber
#     k_bias = ps_vbcrecom[0] / ics.cosmo["h"]# "h Mpc**-1"
    
#     return k_bias, b_cdm, b_b

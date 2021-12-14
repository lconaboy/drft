import os
import yt
import sys
import numpy as np
from unyt import unyt_array
import matplotlib.pyplot as plt
from grafic_tools import load_snapshot


def vbc_calc(ad, u='km/s'):
    vcf = (('deposit', 'DM_cic_velocity_x'),
           ('deposit', 'DM_cic_velocity_y'),
           ('deposit', 'DM_cic_velocity_z'))
    vbf = (('gas', 'velocity_x'),
           ('gas', 'velocity_y'),
           ('gas', 'velocity_z'))

    b = ad[vbf[0]]  # get size
    # the ds has the units, unyt does not
    units = ad[vbf[0]].in_units(u).units
    vbc = unyt_array(np.zeros(b.shape), units ** 2.)

    for b, c in zip(vbf, vcf):
        vbc += (ad[b].in_units(u) - ad[c].in_units(u)) ** 2.

    vbc = np.sqrt(vbc)
    vbc = vbc.astype(np.float32)

    return vbc


def yt_vbc(unbiased_path, ic_path):
    
    ds = yt.load(os.path.join(unbiased_path, 'output_00001'))
    levelmin = ds['levelmin']
    levelmax = ds['levelmax']  # this is simulation levelmax, not ICs levelmax

    # Do levelmin first. Could be done inside the loop, but I had it
    # outside from when I was plotting vbc slices.
    s = load_snapshot(ic_path, ds['levelmin'], 'deltab')

    ad = ds.covering_grid(level=0, left_edge=[0, 0, 0],
                          dims=ds.domain_dimensions)

    u = 'km/s'
    vbc = vbc_calc(ad, u)
    s.write_field(vbc, 'vbc')
    
    # vmin = vbc.min()
    # vmax = vbc.max()
    # plot_vbc_slices(vbc, levelmin, vmin, vmax)
    
    # Now do the other levels. To find out the actual levelmin we
    # would have to load up all the particle masses and do a min/max,
    # which is quite expensive. Instead, we just go until we run out
    # of files.
    for ilevel in range(levelmin+1, levelmax+1):
        try:
            s = load_snapshot(ic_path, ilevel, 'deltab')
            dims = s.n
            left_edge = (s.xoff / s.dx) / 2**ilevel
            ad = ds.covering_grid(ilevel-levelmin, left_edge=left_edge,
                                  dims=dims)
            vbc = vbc_calc(ad, u)
            s.write_field(vbc, 'vbc')

            # plot_vbc_slices(vbc, ilevel, vmin, vmax)
            
        except FileNotFoundError:
            return


def plot_vbc_slices(vbc, ilevel, vmin=None, vmax=None):

    if vmin is None: vmin = 0
    if vmax is None: vmax = vbc.max()

    for i in range(vbc.shape[2]):
        plt.imsave(f'yt_{ilevel}_vbc{i}.png', vbc[:, :, i],
                   vmin=vmin, vmax=vmax)                   

    
# def
# # Now do grafic
# vbc = np.zeros(ds.domain_dimensions, dtype=np.float32)

# for i in 'xyz':
#     vbc += (load_snapshot('./', 7, 'velb'+i).load_box() -
#             load_snapshot('./', 7, 'velc'+i).load_box()) ** 2.

# vbc = np.sqrt(vbc)
# for i in range(vbc.shape[2]):
#     plt.imsave(f'vbc{i}.png', vbc[:, :, i], vmin=vmin, vmax=vmax)

# s.write_field(vbc, 'vbc')


if __name__ == '__main__':
    unbiased_path = sys.argv[1]
    ic_path = sys.argv[2]

    yt_vbc(unbiased_path, ic_path)

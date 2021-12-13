import os
import yt
import sys
import numpy as np
from unyt import unyt_array
import matplotlib.pyplot as plt
from grafic_tools import load_snapshot


def yt_vbc(unbiased_path, ic_path):
    
    ds = yt.load(os.path.join(unbiased_path, 'output_00001'))
    s = load_snapshot(ic_path, ds['levelmin'], 'deltab')

    ad = ds.covering_grid(level=0, left_edge=[0, 0, 0],
                          dims=ds.domain_dimensions)

    vcf = (('deposit', 'DM_cic_velocity_x'),
           ('deposit', 'DM_cic_velocity_y'),
           ('deposit', 'DM_cic_velocity_z'))
    vbf = (('gas', 'velocity_x'),
           ('gas', 'velocity_y'),
           ('gas', 'velocity_z'))

    u = 'km/s'
    b = ad[vbf[0]]  # get size
    # the ds has the units, unyt does not
    units = ad[vbf[0]].in_units(u).units
    vbc = unyt_array(np.zeros(b.shape), units ** 2.)

    for b, c in zip(vbf, vcf):
        vbc += (ad[b].in_units(u) - ad[c].in_units(u)) ** 2.

    vbc = np.sqrt(vbc)
    vbc = vbc.astype(np.float32)

    s.write_field(vbc, 'vbc')

    plot_vbc_slices(vbc)
    
def plot_vbc_slices(vbc):
    vmin = 0
    vmax = vbc.max()

    for i in range(vbc.shape[2]):
        plt.imsave(f'yt_vbc{i}.png', vbc[:, :, i], vmin=vmin, vmax=vmax)

    
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

import numpy as np
import matplotlib.pyplot as plt

import grafic_tools as grafic
from mpl_toolkits.axes_grid1 import ImageGrid

level = 14
# path1 = '/snap7/scratch/dp004/dc-cona1/bd/ics_new/halo7248_972/halo7248_972/'
# path2 = path1 + 'level_{0:03d}/ics_ramses_vbc/'.format(level)
# field = 'vbc'
path1 = '/cosma6/data/dp004/dc-cona1/bd/ics/halo7248_972/halo7248_972/'
path2 = '/cosma6/data/dp004/dc-cona1/bd/ics/halo7248_972/halo7248_972_b/'
field = 'deltab'

# Get data
dat1 = grafic.load_snapshot(path1, level, field).plot_slice(nslc=3, ret_im=True)
dat2 = grafic.load_snapshot(path2, level, field).plot_slice(nslc=3, ret_im=True)
dat = [dat1, dat2]
ttl = ['unbiased', 'biased']
vmin = np.min(dat)
vmax = np.max(dat)

# Hacky
l = 100.*972./16384.
ex = (50. - l, 50. + l, 50.-l, 50.+l)

fig = plt.figure(figsize=(12, 9))

grd = ImageGrid(fig, 111,
                nrows_ncols=(1, 2),
                axes_pad=0.15,
                share_all=True,
                cbar_location="bottom",
                cbar_mode="single",
                cbar_size="2%",
                cbar_pad=0.5)

for i, ax in enumerate(grd):
    im = ax.imshow(dat[i], vmin=vmin, vmax=vmax, cmap='gist_stern', origin='lower', extent=ex)
    ax.set_xlabel('h$^{{-1}}$ Mpc')
    ax.set_ylabel('h$^{{-1}}$ Mpc')
    # ax.set_title(ttl[i])
    ax.text(0.1, 0.9, ttl[i], fontsize=22, color='w', transform=ax.transAxes)

cb = ax.cax.colorbar(im)
# cb.set_label_text('|v_\mathsf{{bc}}|')
cb.set_label_text('$\delta_{{\\rm b}}$')
cb.solids.set_rasterized(True) 
ax.cax.toggle_label(True)


plt.tight_layout(rect=[0,0,0.9,1])
# plt.tight_layout()
plt.savefig(field+'_'+str(level)+'_slice_comp.pdf', bbox_inches='tight')
plt.show()

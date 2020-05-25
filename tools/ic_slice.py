import numpy as np
import matplotlib.pyplot as plt

import grafic_tools as grafic
from mpl_toolkits.axes_grid1 import ImageGrid

level = 14q
path1 = '/path/to/ics/'
path2 = path1 + 'level_{0:03d}/ics_ramses_vbc/'.format(level)
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

fig = plt.figure(figsize=(12.2, 6))

grd = ImageGrid(fig, 111,
                nrows_ncols=(1, 2),
                axes_pad=0.15,
                share_all=True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="7%",
                cbar_pad=0.15)

for i, ax in enumerate(grd):
    im = ax.imshow(dat[i], vmin=vmin, vmax=vmax, cmap='gist_stern', origin='lower', extent=ex)
    ax.set_xlabel('h$^{{-1}}$ Mpc')
    ax.set_ylabel('h$^{{-1}}$ Mpc')
    ax.set_title(ttl[i])

cb = ax.cax.colorbar(im)
cb.set_label_text('$\delta_b$')
cb.solids.set_rasterized(True) 
ax.cax.toggle_label(True)


plt.tight_layout(rect=[0,0,0.9,1])
# plt.tight_layout()
plt.savefig(field+'_'+str(level)+'_slice_comp.pdf')
plt.show()

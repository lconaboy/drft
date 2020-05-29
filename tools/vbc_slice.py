import numpy as np
import matplotlib.pyplot as plt
from grafic_tools import load_snapshot

def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar


paths = ['./']#, '../test_ics/ics_camb_novel/']
level = 14

# Hacky
l = 100.*972./16384.
ex = (50. - l, 50. + l, 50.-l, 50.+l)

# fields = ['deltab']
# fields = ['poscx', 'poscy', 'poscz', 'velcx', 'velcy', 'velcz', 'velbx', 'velby', 'velbz', 'deltab', 'vbc']
# fields = ['deltab', 'deltac', 'poscx_cont', 'poscy_cont', 'poscz_cont', 'velcx_cont', 'velcy_cont', 'velcz_cont', 'velbx_cont', 'velby_cont', 'velbz_cont', 'vbc_cont']
# fields = ['velcx', 'velbx', 'velcgx', 'velcy', 'velby', 'velcgy', 'velcz', 'velbz', 'velcgz', 'vbc']
fields = ['vbc']

for path in paths:
    print(path)
    for field in fields:
        ics = load_snapshot(path, level, field)
        slc = ics.plot_slice(nslc=20, ret_im=True)
        
        fig, ax = plt.subplots(figsize=(6, 6))
        im = ax.imshow(slc, cmap='viridis', origin='lower', extent=ex)
        cb = colorbar(im)
        cb.set_label(r'$\mid$v$_\mathsf{{bc}}$$\mid$ (km s$^{{-1}}$)')
        cb.solids.set_rasterized(True) 
        ax.set_xlabel('h$^{{-1}}$ Mpc')
        ax.set_ylabel('h$^{{-1}}$ Mpc')
        plt.tight_layout()
        plt.savefig(field+'_'+str(level)+'_slice.pdf')

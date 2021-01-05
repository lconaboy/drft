import numpy as np
import numpy.ma as ma
import grafic_tools as grafic
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def diff(x0, x1):
    """Calculates the fractional difference between x0 and x1.

    :param x0: (array) first array
    :param x1: (array) second array
    :returns: fractional difference
    :rtype: (array)

    """
    
    return (x1 - x0) / x0


def get_2dhist(v, m):

    h, ve, me = np.histogram2d(v, m, bins=128)
    
    return h, ve, me


def plot_2dhist(v, m, h, fn='hist.png', log=False):
    fig, ax = plt.subplots(figsize=(6, 6))
    
    if log:
        # Mask out zero values for LogNorm
        mh = ma.masked_array(h, mask=h==0)
        im = ax.pcolormesh(ve, me, mh,
                           norm=colors.LogNorm(vmin=mh.min(), vmax=mh.max()))
    else:
        im = ax.pcolormesh(ve, me, h)

    ax.set_xlabel('v$_{{\\sf bc}}$ (km s$^{{-1}}$)')
    ax.set_ylabel('$\Delta$M/M$_{{\\sf i}}$')
    cb = fig.colorbar(mappable=im)
    cb.set_label('counts')
    fig.savefig(fn, bbox_inches='tight')    

    
level = 14
b_path = '/cosma6/data/dp004/dc-cona1/bd/ics/halo7248_972/biased'
u_path = '/cosma6/data/dp004/dc-cona1/bd/ics/halo7248_972/halo7248_972'

us = grafic.load_snapshot(u_path, level, 'deltab')
bs = grafic.load_snapshot(b_path, level, 'deltab')

# \delta + 1
N = us.n  # [64, 64, 64]
u = (us.load_patch([0, 0, 0], N) + 1).ravel()
b = (bs.load_patch([0, 0, 0], N) + 1).ravel()

dm = diff(u, b)

del u, b

# Raw v_bc
vs = grafic.load_snapshot(u_path, level, 'vbc')
v = vs.load_patch([0, 0, 0], N).ravel()
h, ve, me = get_2dhist(v, dm)
plot_2dhist(ve, me, h, fn='{0}_{1}_hist.png'.format('mass_diff', level))
plot_2dhist(ve, me, h, fn='{0}_{1}_hist.png'.format('mass_diff_log', level),
            log=True)

# Patch v_bc
# vs = grafic.load_snapshot(u_path, level, 'vbc_patch')
# v = vs.load_patch([0, 0, 0], N).ravel()
# fig = plot_2dhist(v, dm)
# fig.savefig('{0}_{1}_hist.png'.format('mass_diff_patch', level),
#             bbox_inches='tight')

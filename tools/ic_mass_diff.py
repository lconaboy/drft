import grafic_tools as grafic
import matplotlib.pyplot as plt


def diff(x0, x1):
    """Calculates the fractional difference between x0 and x1.

    :param x0: (array) first array
    :param x1: (array) second array
    :returns: fractional difference
    :rtype: (array)

    """
    
    return (x1 - x0) / x0


def plot_2dhist(v, m):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.hist2d(v, m, bins=128)
    ax.set_xlabel('v$_{{\\sf bc}}$ (km s$^{{-1}}$)')
    ax.set_ylabel('$\Delta$M/M')
    cb = fig.colorbar()
    cb.set_label('counts')
    return fig


level = 14
b_path = ''
u_path = ''

us = grafic.load_snapshot(u_path, level, 'deltab')
bs = grafic.load_snapshot(b_path, level, 'deltab')

# \delta + 1
u = us.load_patch([0, 0, 0], 64) + 1
b = bs.load_patch([0, 0, 0], 64) + 1 

dm = diff(u, b)

del u, b

vs = grafic.load_snapshot(u_path, level, 'vbc')
v = vs.load_patch([0, 0, 0], 64)

fig = plot_2dhist(v, dm)
fig.savefig('{0}_{1}_hist.pdf'.format('delta', level))

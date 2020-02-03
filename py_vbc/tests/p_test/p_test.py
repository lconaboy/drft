import os
import numpy as np
import matplotlib.pyplot as plt

from py_vbc.constants import *
from py_vbc.derivatives import calc_derivs
from py_vbc.spectra import calc_power_spec, calc_delta

path, _ = os.path.split(__file__)

def plot_comparison(x0, yc0, yb0, x1, yc1, yb1, camb=None):
    """
    Function to plot the transfer functions of both my (py_vbc) and
    CICsASS outputs, as well as the ratio of the two functions.

    :param x0: py_vbc logk values
    :param yc0: py_vbc TF values for dark matter
    :param yb0: py_vbc TF values for baryons
    :param x1: CICsASS logk values
    :param yc1: CICsASS TF values for dark matter
    :param yb1: CICsASS TF values for baryons
    :param d: for labels, plot filename and setting log scale for TFs -- should be either '' (empty string) or 'd', depending on whether the differential is being plotted
    :param yl: adjust ylims, useful for plotting TFs -- should be either True or False
    :returns: 
    :rtype:

    """

    from matplotlib.lines import Line2D

    fig = plt.figure(figsize=(6, 9))
#    x1 = x0

    # Makes the legend clearer
    if camb is not None:
        tf_legend = [Line2D([0], [0], color='m', linestyle='-'),
                     Line2D([0], [0], color='r', linestyle='-'),
                     Line2D([0], [0], color='c', linestyle='-'),
                     Line2D([0], [0], color='k', linestyle='-'),
                     Line2D([0], [0], color='k', linestyle='--'),
                     Line2D([0], [0], color='k', linestyle=':'),]
        tf_labels = ['CDM', r'baryons', 'total', 'py\_vbc', 'CICsASS', 'CAMB']
    else:
        tf_legend = [Line2D([0], [0], color='m', linestyle='-'),
                 Line2D([0], [0], color='r', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='--')]
        tf_labels = ['CDM', r'baryons', 'py\_vbc', 'CICsASS']

        
    # For splitting up the subplots into neat ratios
    ax0 = plt.subplot2grid((6, 4), (0, 0), colspan=4, rowspan=4)
    ax1 = plt.subplot2grid((6, 4), (4, 0), colspan=4, rowspan=1, sharex=ax0)
    ax2 = plt.subplot2grid((6, 4), (5, 0), colspan=4, rowspan=2, sharex=ax0)

    # Plot the py_vbc output
    ax0.plot(x0, yc0, 'm')
    ax0.plot(x0, yb0, 'r')
    # ax0.plot(x0, yt0, 'c')
    # Plot the CICsASS output
    ax0.plot(x1, yc1, color='m', linestyle='--')
    ax0.plot(x1, yb1, color='r', linestyle='--')
    # ax0.plot(x1, yt1, color='c', linestyle='--')
    if camb is not None:
        # Plot the CAMB comparison
        ax0.plot(camb[0], camb[1], color='c', linestyle=':')
        ax0.set_xlim([np.min(x0), np.max(x0)])

    ax0.legend(tf_legend, tf_labels)

    ax1.plot(x0, yc0/yc0, color='m')
    ax1.plot(x0, yc0/yc1, color='m', linestyle='--')
    ax2.plot(x0, yb0/yb0, color='r')
    ax2.plot(x0, yb0/yb1, color='r', linestyle='--')

    # There was some weird formatting of the y-axis going on, so turn this off
    ax1.get_yaxis().get_major_formatter().set_useOffset(False)
    ax2.get_yaxis().get_major_formatter().set_useOffset(False)

    ax0.set_ylabel(r'k$^3$P(k)/2$\pi^2$')
    ax1.set_ylabel(r'Ratio')
    ax2.set_ylabel(r'Ratio')
    ax2.set_xlabel(r'k (Mpc$^{-1}$)')
    plt.xscale('log')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.savefig(path+'/ca_pyvbc_delta.pdf')


p_c = np.loadtxt(path+'/p.dat')
p_c_t = p_c[:, 1] + p_c[:, 2]

g_pv = calc_derivs(p_c[:, 0], vbc=0.0, zstart=1000.0, zend=50.0, dz=3.0)
p_pv_c, p_pv_b, _, _ = calc_power_spec(p_c[:, 0], g_pv, zstart=1000.0)
d_pv_c = calc_delta(p_c[:, 0], p_pv_c)
d_pv_b = calc_delta(p_c[:, 0], p_pv_b)

#d_pv_t = d_pv_c + d_pv_b

# c_ps = np.loadtxt(path+'test_matterpower_z50.dat')
# kh_c_ps = c_ps[:, 0]
# d_c_ps = calc_delta(c_ps[:, 0], c_ps[:, 1])

# camb = [kh_c_ps, d_c_ps]

plot_comparison(p_c[:, 0], d_pv_c, d_pv_b, p_c[:, 0], p_c[:, 1], p_c[:, 2])# , p_c_t, camb)

plt.show()

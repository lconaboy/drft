import numpy as np
import matplotlib.pyplot as plt

from py_vbc.constants import *
from py_vbc.derivatives import calc_derivs

path = 'py_vbc/tests/g_test/'

def plot_comparison(x0, yc0, yb0, yc1, yb1):
    """
    Function to plot the transfer functions of both my (py_vbc) and
    CICsASS outputs, as well as the ratio of the two functions.

    The two 

    :param x0: py_vbc x values
    :param yc0: py_vbc TF values for dark matter
    :param yb0: py_vbc TF values for baryons
    :param yc1: CICsASS TF values for dark matter
    :param yb1: CICsASS TF values for baryons
    :returns: the plot (which has been saved)
    :rtype: Figure

    """

    from matplotlib.lines import Line2D

    fig = plt.figure(figsize=(6, 9))

    # Makes the legend clearer
    tf_legend = [Line2D([0], [0], color='seagreen', linestyle='-'),
                 Line2D([0], [0], color='orange', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='--'),]
    tf_labels = ['g$_\mathrm{c}$', r'g$_\mathrm{b}$', 'py\_vbc', 'CICsASS']

    # For splitting up the subplots into neat ratios
    ax0 = plt.subplot2grid((6, 4), (0, 0), colspan=4, rowspan=4)
    ax1 = plt.subplot2grid((6, 4), (4, 0), colspan=4, rowspan=1, sharex=ax0)
    ax2 = plt.subplot2grid((6, 4), (5, 0), colspan=4, rowspan=2, sharex=ax0)

    # Plot the py_vbc output
    ax0.plot(x0, yc0, 'seagreen')
    ax0.plot(x0, yb0, 'orange')
    # Plot the CICsASS output
    ax0.plot(x0, yc1, color='seagreen', linestyle='--')
    ax0.plot(x0, yb1, color='orange', linestyle='--')
    ax0.legend(tf_legend, tf_labels)

    ax1.plot(x0, yc0/yc0, color='seagreen')
    ax1.plot(x0, yc0/yc1, color='seagreen', linestyle='--')
    ax2.plot(x0, yb0/yb0, color='orange')
    ax2.plot(x0, yb0/yb1, color='orange', linestyle='--')

    # There was some weird formatting of the y-axis going on, so turn this off
    ax1.get_yaxis().get_major_formatter().set_useOffset(False)
    ax2.get_yaxis().get_major_formatter().set_useOffset(False)

    ax0.set_ylabel(r'g(k)')
    ax1.set_ylabel(r'Ratio')
    ax2.set_ylabel(r'Ratio')
    ax2.set_xlabel(r'k (Mpc$^{-1}$)')
    plt.xscale('log')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.savefig(path+'/ca_pyvbc_g.pdf')

g_c = np.loadtxt(path+'g.dat')
g_pv = calc_derivs(g_c[:, 0], vbc=0.0, zstart=1000.0, zend=50.0, dz=3.0)

plot_comparison(g_c[:, 0], g_pv[:, 0], g_pv[:, 2], g_c[:, 1], g_c[:, 2])
plt.show()

import os
import numpy as np
import matplotlib.pyplot as plt

from py_vbc.constants import *
from py_vbc.interpolations import interpolate_tf

path, _ = os.path.split(__file__)

# Utility function for plotting
def plot_comparison(x0, yc0, yb0, x1, yc1, yb1, d='', yl=False):
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

    # Makes the legend clearer
    tf_legend = [Line2D([0], [0], color='m', linestyle='-'),
                 Line2D([0], [0], color='c', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='--'),]
    tf_labels = [r'{}T$_\mathrm{{k, c}}$'.format(d), r'{}T$_\mathrm{{k, b}}$'.format(d),
                 r'py\_vbc', r'CICsASS']

    # For splitting up the subplots into neat ratios
    ax0 = plt.subplot2grid((6, 4), (0, 0), colspan=4, rowspan=4)
    ax1 = plt.subplot2grid((6, 4), (4, 0), colspan=4, rowspan=1, sharex=ax0)
    ax2 = plt.subplot2grid((6, 4), (5, 0), colspan=4, rowspan=2, sharex=ax0)

    # Plot the py_vbc output
    ax0.plot(x0, yc0, 'm')
    ax0.plot(x0, yb0, 'c')
    # Plot the CICsASS output
    ax0.plot(x1, yc1, 'm--')
    ax0.plot(x1, yb1, 'c--')
    if d == '':
        ax0.set_yscale('log')
    if yl:
        ax0.set_ylim([1e-4, ax0.get_ylim()[1]])
    ax0.legend(tf_legend, tf_labels)

    ax1.plot(x0, yc0/yc0, 'm')
    ax1.plot(x0, yc0/yc1, 'm--')
    ax2.plot(x0, yb0/yb0, 'c')
    ax2.plot(x0, yb0/yb1, 'c--')

    # There was some weird formatting of the y-axis going on, so turn this off
    ax1.get_yaxis().get_major_formatter().set_useOffset(False)
    ax2.get_yaxis().get_major_formatter().set_useOffset(False)

    ax0.set_ylabel(r'{}T$_\mathrm{{k}}$'.format(d))
    ax1.set_ylabel(r'Ratio')
    ax2.set_ylabel(r'Ratio')
    ax2.set_xlabel(r'log$_{10}$(k [Mpc$^{-1}$])')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.savefig(path+'/ca_pyvbc_{}tf.pdf'.format(d))


# Interpolate the TFs using py_vbc
c_spline, dc_spline = interpolate_tf('c', z=1000.0, dz=3.0)
b_spline, db_spline = interpolate_tf('b', z=1000.0, dz=3.0)

# Load up the CICsASS data
vals_c = np.loadtxt(path+'/c.dat')
k_c0 = vals_c[:, 0]
c_c = vals_c[:, 1]
dc_c = vals_c[:, 2]

vals_b = np.loadtxt(path+'/b.dat')
k_c1 = vals_b[:, 0]
b_c = vals_b[:, 1]
db_c = vals_b[:, 2]

# lk_c = k_c0

# Compute the values of the splines using the CICsASS lk values as the
# new x_values
c = c_spline(k_c0)
b = b_spline(k_c0)
dc = dc_spline(k_c0)
db = db_spline(k_c0)

# Convert to log10
lk_c = np.log10(k_c0)

# Plot the TFs
plot_comparison(lk_c, c, b, lk_c, c_c, b_c, d='', yl=True)
# Plot the dTFs
plot_comparison(lk_c, dc, db, lk_c, dc_c, db_c, d='d', yl=False)

plt.show()

import os
import numpy as np
import matplotlib.pyplot as plt

from py_vbc.constants import *
from py_vbc.interpolations import interpolate_tf

path, _ = os.path.split(__file__)

# Utility function for plotting
def plot_comparison(x_pv, yc_pv, yb_pv, x_c, yc_c, yb_c, d='', yl=False):
    """
    Function to plot the transfer functions of both my (py_vbc) and
    CICsASS outputs, as well as the ratio of the two functions.

    :param x_pv: py_vbc logk values
    :param yc_pv: py_vbc TF values for dark matter
    :param yb_pv: py_vbc TF values for baryons
    :param x_c: CICsASS logk values
    :param yc_c: CICsASS TF values for dark matter
    :param yb_c: CICsASS TF values for baryons
    :param d: for labels, plot filename and setting log scale for TFs -- should be either '' (empty string) or 'd', depending on whether the differential is being plotted
    :param yl: adjust ylims, useful for plotting TFs -- should be either True or False
    :returns: 
    :rtype:

    """

    from matplotlib.lines import Line2D

    fig = plt.figure(figsize=(4.5, 6))

    # Makes the legend clearer
    tf_legend = [Line2D([0], [0], color='seagreen', linestyle='-'),
                 Line2D([0], [0], color='orange', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='--'),]
    if d == '':
        ylabel = r'T(k)'
        y1label = r'$\Delta$T/T (\%)'

    else:
        ylabel = r'$\dot{\sf T}=\delta$T(k)$/\delta$z'
        y1label = r'$\Delta \dot{\sf T}/\dot{\sf T}$ (\%)'


    tf_labels = [r'c', r'b',
                 r'CICsASS', r'py\_vbc']

    # For splitting up the subplots into neat ratios
    ax0 = plt.subplot2grid((6, 4), (0, 0), colspan=4, rowspan=4)
    ax1 = plt.subplot2grid((6, 4), (4, 0), colspan=4, rowspan=1, sharex=ax0)
    ax2 = plt.subplot2grid((6, 4), (5, 0), colspan=4, rowspan=2, sharex=ax0)

    # Plot the py_vbc output
    ax0.plot(x_c, yc_c, 'seagreen')
    # Plot the CICsASS output
    ax0.plot(x_pv, yc_pv, c='seagreen', ls='--')
    ax0.plot(x_c, yb_c, 'orange')
    ax0.semilogx(x_pv, yb_pv, c='orange', ls='--')
    if d == '':
        ax0.set_yscale('log')
    # if yl:
    #     ax0.set_ylim([1e-4, ax0.get_ylim()[1]])
    ax0.legend(tf_legend, tf_labels)

    # ax1.plot(x0, yc_pv/yc_pv, 'seagreen')
    ax1.axhline(0., c='seagreen')
    ax1.semilogx(x_pv, 100 * (yc_pv-yc_c)/yc_c, c='seagreen', ls='--')
    ax2.axhline(0., c='orange')
    # ax2.plot(x_pv, yb_pv/yb_pv, 'c')
    ax2.semilogx(x_pv, 100 * (yb_pv-yb_c)/yb_c, c='orange', ls='--')

    # There was some weird formatting of the y-axis going on, so turn this off
    ax1.get_yaxis().get_major_formatter().set_useOffset(False)
    ax2.get_yaxis().get_major_formatter().set_useOffset(False)

    ax0.set_ylabel(ylabel)
  
    ax_inv = fig.add_subplot(313, frameon=False)
    ax_inv.tick_params(which='both', labelbottom=False,
                       bottom=False, top=False,
                       left=False, right=False)
    ax_inv.set_yticks([])
    
    ax_inv.set_ylabel(y1label, labelpad=35)
    ax0.tick_params(labelbottom=False)
    ax1.tick_params(labelbottom=False)

    mima = [x_pv.min(), x_pv.max()]
    
    ax0.set_xlim(mima)
    ax1.set_xlim(mima)
    ax2.set_xlim(mima)
    
    ax1.set_ylim([-0.012, 0.012])
    if d == '':
        ax2.set_ylim([-0.25, 0.25])
    else:
        ax2.set_ylim([-1.25, 0.5])
        

    # ax2.set_ylabel(r'Ratio')
    ax2.set_xlabel(r'k (Mpc$^{-1}$)')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.1)
    fig.savefig(path+'/ca_pyvbc_{}tf.pdf'.format(d),
                bbox_inches='tight', pad_inches=0.02)


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
lk_c = k_c0

# Plot the TFs
plot_comparison(lk_c, c, b, lk_c, c_c, b_c, d='', yl=True)
# Plot the dTFs
plot_comparison(lk_c, dc, db, lk_c, dc_c, db_c, d='d', yl=False)

plt.show()

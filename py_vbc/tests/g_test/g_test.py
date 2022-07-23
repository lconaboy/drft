import numpy as np
import matplotlib.pyplot as plt

from py_vbc.constants import *
from py_vbc.derivatives import calc_derivs

path = 'py_vbc/tests/g_test/'

def plot_comparison(x0, yc_pv, yb_pv, yc_c, yb_c):
    """
    Function to plot the transfer functions of both my (py_vbc) and
    CICsASS outputs, as well as the ratio of the two functions.

    The two 

    :param x0: py_vbc x values
    :param yc_pv: py_vbc TF values for dark matter
    :param yb_pv: py_vbc TF values for baryons
    :param yc_c: CICsASS TF values for dark matter
    :param yb_c: CICsASS TF values for baryons
    :returns: the plot (which has been saved)
    :rtype: Figure

    """

    from matplotlib.lines import Line2D

    fig = plt.figure(figsize=(4.5, 6), constrained_layout=True)

    # Makes the legend clearer
    tf_legend = [Line2D([0], [0], color='seagreen', linestyle='-'),
                 Line2D([0], [0], color='orange', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='--'),]
    tf_labels = ['c', 'b', 'CICsASS', 'py\_vbc']

    # For splitting up the subplots into neat ratios
    ax0 = plt.subplot2grid((6, 4), (0, 0), colspan=4, rowspan=4)
    ax1 = plt.subplot2grid((6, 4), (4, 0), colspan=4, rowspan=1, sharex=ax0)
    ax2 = plt.subplot2grid((6, 4), (5, 0), colspan=4, rowspan=2, sharex=ax0)

    ax_inv = fig.add_subplot(313, frameon=False)
    ax_inv.tick_params(which='both', labelbottom=False,
                       bottom=False, top=False,
                       left=False, right=False)
    ax_inv.set_ylabel(r'$\Delta \delta/\delta$ (\%)', labelpad=35)
    ax_inv.set_yticks([])

    
    # Plot the py_vbc output
    ax0.plot(x0, yc_pv, 'seagreen', linestyle='--')
    ax0.plot(x0, yb_pv, 'orange', linestyle='--')
    # Plot the CICsASS output
    ax0.plot(x0, yc_c, color='seagreen')
    ax0.plot(x0, yb_c, color='orange')
    ax0.legend(tf_legend, tf_labels)

    ax1.axhline(0., color='seagreen')
    ax1.plot(x0, 100. * (yc_pv - yc_c)/yc_c, color='seagreen', linestyle='--')
    # ax2.plot(x0, (yb_c - yb_c)/yb_c, color='orange')
    ax2.axhline(0., color='orange')
    ax2.plot(x0, 100. * (yb_pv - yb_c)/yb_c, color='orange', linestyle='--')

    # There was some weird formatting of the y-axis going on, so turn this off
    ax1.get_yaxis().get_major_formatter().set_useOffset(False)
    ax2.get_yaxis().get_major_formatter().set_useOffset(False)

    mi = x0.min()
    ma = x0.max()
    xmima = [mi, ma]

    ax0.tick_params(labelbottom=False)
    ax1.tick_params(labelbottom=False)
    ax0.set_xlim(xmima)
    ax1.set_xlim(xmima)
    ax2.set_xlim(xmima)
    ax1.set_ylim([-0.01, 0.01])
    ax2.set_ylim([-3, 3])
    ax0.set_ylabel(r'$\delta$(k)')
    # ax1.set_ylabel(r'$\Delta/\delta_{\sf CICsASS}$ (\%)')
    # ax2.set_ylabel(r'$\Delta/\delta_{\sf CICsASS}$ (\%)')
    # fig.text(0.04, 0.1667, s=r'\Delta',
    #          va='center', ha='left', rotation='vertical')
    ax2.set_xlabel(r'k (Mpc$^{-1}$)')
    ax0.set_xscale('log')
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax_inv.set_xscale('log')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.1)
    fig.savefig(path+'/ca_pyvbc_g.pdf', bbox_inches='tight', pad_inches=0.04,
                va='center', ha='left')

g_c = np.loadtxt(path+'g.dat')
print(g_c.min(), g_c.max())
# g_pv = calc_derivs(g_c[:, 0], vbc=0.0, zstart=1000.0, zend=50.0,
#                    dz=3.0, verbose=True)

# np.savetxt('g_pv_isothermal.dat', g_pv)

g_pv = np.loadtxt('./g_pv_non_isothermal.dat')

plot_comparison(g_c[:, 0], g_pv[:, 0], g_pv[:, 2], g_c[:, 1], g_c[:, 2])
plt.show()

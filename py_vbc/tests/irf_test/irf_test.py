import os
import numpy as np
import matplotlib.pyplot as plt

from py_vbc.interpolations import interpolate_recfast

""" Script to plot the comparison between the recfast interpolation
done by py_vbc and CICsASS."""

path, _ = os.path.split(__file__)

def plot_comparison(x0, y0, x1, y1, qty):
    """
    Function to plot the recfast interpolations of both my (py_vbc) and
    CICsASS outputs, as well as the ratio of the two functions.

    :param x0: py_vbc z values
    :param y0: py_vbc recfast values
    :param x1: CICsASS z values
    :param y1: CICsASS recfast values
    :param ylabel: (str) quantity being plotted on the y-axis, either "T" or "xe"
    :returns: 
    :rtype:

    """

    from matplotlib.lines import Line2D

    fig = plt.figure(figsize=(6, 7.5))

    # Dict for ylabel
    ylabel = {'T':'T (K)', 'xe':r'x$_\mathrm{e}$'}
    
    # Makes the legend clearer
    rf_legend = [Line2D([0], [0], color='k', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='--'),]
    rf_labels = ['py\_vbc', 'CICsASS']

    # For splitting up the subplots into neat ratios
    ax0 = plt.subplot2grid((5, 4), (0, 0), colspan=4, rowspan=4)
    ax1 = plt.subplot2grid((5, 4), (4, 0), colspan=4, rowspan=1, sharex=ax0)


    # Plot the py_vbc output
    ax0.plot(x0, y0, 'k')
    # Plot the CICsASS output
    ax0.plot(x1, y1, 'k--')
    ax0.legend(rf_legend, rf_labels)

    ax1.plot(x0, y0/y0, 'k')
    ax1.plot(x0, y0/y1, 'k--')

    # There was some weird formatting of the y-axis going on, so turn this off
    ax1.get_yaxis().get_major_formatter().set_useOffset(False)

    # # Log scale for redshift axis
    ax0.set_xscale('log')
    ax1.set_xscale('log')

    # Log scale for quantity axis too
    ax0.set_yscale('log')

    ax0.set_ylabel(ylabel[qty])
    ax1.set_ylabel(r'Ratio')
    ax1.set_xlabel(r'z')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.savefig(path+'/ca_pyvbc_{}.pdf'.format(qty))


# Load up the CICsASS data, which contains a lot of extraneous values
r_c = np.loadtxt(path+'/r.dat')
z_c = r_c[:, 0]
T_c = r_c[:, 1]
xe_c = r_c[:, 2]

# Calculate splines and generate range of data
T_spline, xe_spline = interpolate_recfast()
T = T_spline(z_c)
xe = xe_spline(z_c)

# Plot
plot_comparison(z_c, T, z_c, T_c, 'T')
plot_comparison(z_c, xe, z_c, xe_c, 'xe')

plt.show()

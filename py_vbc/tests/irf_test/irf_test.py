import os
import numpy as np
import matplotlib.pyplot as plt

from py_vbc.interpolations import interpolate_recfast

""" Script to plot the comparison between the recfast interpolation
done by py_vbc and CICsASS."""

path, _ = os.path.split(__file__)

def plot_comparison(x_pv, y_pv, x_c, y_c, qty):
    """
    Function to plot the recfast interpolations of both my (py_vbc) and
    CICsASS outputs, as well as the ratio of the two functions.

    :param x_pv: py_vbc z values
    :param y_pv: py_vbc recfast values
    :param x_c: CICsASS z values
    :param y_c: CICsASS recfast values
    :param ylabel: (str) quantity being plotted on the y-axis, either "T" or "xe"
    :returns: 
    :rtype:

    """

    from matplotlib.lines import Line2D

    fig = plt.figure(figsize=(4, 5))

    # Dict for ylabel
    ylabel = {'T':'T (K)', 'xe':r'x$_{\sf e}$'}
    y1label = {'T':'T', 'xe': r'x$_{\sf e}$'}
    
    # Makes the legend clearer
    rf_legend = [Line2D([0], [0], color='k', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='--'),]
    rf_labels = ['CICsASS', 'py\_vbc']

    # For splitting up the subplots into neat ratios
    ax0 = plt.subplot2grid((5, 4), (0, 0), colspan=4, rowspan=4)
    ax1 = plt.subplot2grid((5, 4), (4, 0), colspan=4, rowspan=1, sharex=ax0)


    # Plot the py_vbc output
    ax0.plot(x_c, y_c, 'k')
    # Plot the CICsASS output
    ax0.plot(x_pv, y_pv, 'k--')
    ax0.legend(rf_legend, rf_labels)

    ax0.tick_params(labelbottom=False)

    # ax1.plot(x0, y_pv/y_pv, 'k')
    ax1.axhline(0., color='k')
    ax1.plot(x_pv, 100 * (y_pv-y_c)/y_c, 'k--')

    # There was some weird formatting of the y-axis going on, so turn this off
    ax1.get_yaxis().get_major_formatter().set_useOffset(False)

    # # Log scale for redshift axis
    ax0.set_xscale('log')
    ax1.set_xscale('log')

    # Log scale for quantity axis too
    ax0.set_yscale('log')

    ax1.set_ylim([-0.01, 0.01])
    ax0.set_ylabel(ylabel[qty])
    ax1.set_ylabel(r'$\Delta$'+y1label[qty]+'/'+y1label[qty]+r' (\%)')
    ax1.set_xlabel(r'z')
    ax0.set_xlim([x_pv.min(), x_pv.max()])
    ax1.set_xlim([x_pv.min(), x_pv.max()])
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.1)
    fig.savefig(path+'/ca_pyvbc_{}.pdf'.format(qty),
                bbox_inches='tight', pad_inches=0.02)


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

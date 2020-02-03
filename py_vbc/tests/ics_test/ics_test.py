import numpy as np
import matplotlib.pyplot as plt

from py_vbc.utils import hubble
from py_vbc.derivatives import set_ics
from py_vbc.tests.test_constants import *

""" Script to plot the comparison between the ICs set by py_vbc and
CICsASS """

path = 'tests/ics_test/'

def plot_comparison(x0, y0, y1, qty):
    """
    Function to plot the recfast interpolations of both my (py_vbc) and
    CICsASS outputs, as well as the ratio of the two functions.

    :param x0: py_vbc z values
    :param y0: tuple containing the real and imaginary parts of the function to be plotted (py_vbc)
    :param y1: tuple containing the real and imaginary parts of the function to be plotted (CICsASS)
    :param ylabel: (str) quantity being plotted on the y-axis, either "T" or "xe"
    :returns: 
    :rtype:

    """

    from matplotlib.lines import Line2D
    from scipy.interpolate import spline

    fig = plt.figure(figsize=(6, 7.5))

    # Dict for ylabel
    ylabel = {'dc':'\delta_\mathrm{c}', 'db':'\delta_\mathrm{b}',
              'dc_dot':'\dot{\delta_\mathrm{c}}',
              'db_dot':'\dot{\delta_\mathrm{b}}',
              'xe':'\mathrm{x_e}', 'dt':'\delta_\mathrm{T}'}
    
    # Makes the legend clearer
    rf_legend = [Line2D([0], [0], color='k', linestyle='-'),
                 Line2D([0], [0], color='k', linestyle='--'),
                 Line2D([0], [0], color='r', linestyle='-'),
                 Line2D([0], [0], color='g', linestyle='-')]
    rf_labels = ['py\_vbc', 'CICsASS', r'Re(${0}$)'.format(ylabel[qty]), r'Im(${0}$)'.format(ylabel[qty])]

    # For splitting up the subplots into neat ratios
    ax0 = plt.subplot2grid((5, 4), (0, 0), colspan=4, rowspan=4)
    ax1 = plt.subplot2grid((5, 4), (4, 0), colspan=4, rowspan=1, sharex=ax0)


    # Plot the py_vbc output
    ax0.plot(x0, y0[0], 'r')
    ax0.plot(x0, y0[1], 'g')
    # Plot the CICsASS output
    ax0.plot(x0, y1[0], 'r--')
    ax0.plot(x0, y1[1], 'g--')
    ax0.legend(rf_legend, rf_labels)

    ax1.plot(x0, y0[0]/y0[0], 'r')
    ax1.plot(x0, y0[0]/y1[0], 'r--')
    ax1.plot(x0, y0[1]/y0[1], 'g')
    ax1.plot(x0, y0[1]/y1[1], 'g--')

    # There was some weird formatting of the y-axis going on, so turn this off
    ax1.get_yaxis().get_major_formatter().set_useOffset(False)

    # # Log scale for redshift axis
    ax0.set_xscale('log')
    ax1.set_xscale('log')

    # Log scale for quantity axis too
    #ax0.set_yscale('log')

    ax0.set_ylabel('$'+ylabel[qty]+'$')
    ax1.set_ylabel(r'Ratio')
    ax1.set_xlabel(r'log$_{10}$(k/h [Mpc$^{-1}$])')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.savefig(path+'ca_pyvbc_{}.pdf'.format(qty))

ics_c = np.loadtxt(path+'ic.dat')
k_c = ics_c[:, 0]
lk_c = np.log10(k_c/hconst)
delc_c_r = ics_c[:, 1]
delc_c_i = ics_c[:, 2]
delc_dot_c_r = ics_c[:, 3]
delc_dot_c_i = ics_c[:, 4]
delb_c_r = ics_c[:, 5]
delb_c_i = ics_c[:, 6]
delb_dot_c_r = ics_c[:, 7]
delb_dot_c_i = ics_c[:, 8]
tcmb_c = ics_c[:, 9]
xe_c = ics_c[:, 10]
delt_c_r = ics_c[:, 11]
delt_c_i = ics_c[:, 12]

ics = set_ics(lk_c)
qty = ['dc', 'db', 'dc_dot', 'db_dot', 'dt']
for q, i in enumerate([0, 2, 4, 6, 10]):
    plot_comparison(lk_c, (ics[:, i], ics[:, i]), (ics_c[:, i+1], ics_c[:, i+1]), qty=qty[q])

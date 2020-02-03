import os
import numpy as np
import matplotlib.pyplot as plt

from py_vbc.constants import *
from py_vbc.derivatives import calc_derivs
from py_vbc.spectra import calc_power_spec, calc_delta

path, _ = os.path.split(__file__)

def plot_comparison(x0, yc0, yb0, camb=None):
    """
    Function to plot the transfer functions of both my (py_vbc) and
    CICsASS outputs, as well as the ratio of the two functions.

    :param x0: py_vbc logk values
    :param yc0: py_vbc TF values for dark matter
    :param yb0: py_vbc TF values for baryons
    :param d: for labels, plot filename and setting log scale for TFs -- should be either '' (empty string) or 'd', depending on whether the differential is being plotted
    :param yl: adjust ylims, useful for plotting TFs -- should be either True or False
    :returns: 
    :rtype:

    """

    from matplotlib.lines import Line2D

    
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
                     Line2D([0], [0], color='r', linestyle='-'),]
        tf_labels = ['CDM', r'baryons']#, 'py\_vbc', 'CICsASS']

        
    # For splitting up the subplots into neat ratios
    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot the py_vbc output
    ax.plot(x0, yc0, 'm')
    ax.plot(x0, yb0, 'r')
    if camb is not None:
        # Plot the CAMB comparison
        ax.plot(camb[0], camb[1], color='c', linestyle=':')
        ax.set_xlim([np.min(x0), np.max(x0)])

    ax.legend(tf_legend, tf_labels)

    ax.set_ylabel(r'P(v$_\mathsf{bc}$=30)/P(v$_\mathsf{bc}$=0)')
    ax.set_xlabel(r'k (Mpc$^{-1}$)')
    ax.set_title(r'velocity perturbations')
    plt.xscale('log')
    plt.yscale('log')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.savefig(path+'/ca_pyvbc_delta_v.pdf')


# Load up CICsASS power spec, we only use the k-values in the following
p_c = np.loadtxt(path+'/p.dat')
p_c_t = p_c[:, 1] + p_c[:, 2]

# Calculate v_bc = 0.0
g_pv = calc_derivs(p_c[:, 0], vbc=0.0, zstart=1000.0, zend=200.0, dz=3.0)
_, _, p_pv_c0, p_pv_b0 = calc_power_spec(p_c[:, 0], g_pv, zstart=1000.0)

# Calculate v_bc = 30.0
g_pv = calc_derivs(p_c[:, 0], vbc=30.0, zstart=1000.0, zend=200.0, dz=3.0)
_, _, p_pv_c30, p_pv_b30 = calc_power_spec(p_c[:, 0], g_pv, zstart=1000.0)

# Calculate the bias = v_bc0/v_bc30
bias_c = p_pv_c30 / p_pv_c0
bias_b = p_pv_b30 / p_pv_b0

#d_pv_t = d_pv_c + d_pv_b

# c_ps = np.loadtxt(path+'test_matterpower_z50.dat')
# kh_c_ps = c_ps[:, 0]
# d_c_ps = calc_delta(c_ps[:, 0], c_ps[:, 1])

# camb = [kh_c_ps, d_c_ps]

plot_comparison(p_c[:, 0], bias_c, bias_b)

plt.show()

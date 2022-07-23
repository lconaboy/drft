import numpy as np
import matplotlib.pyplot as plt

from py_vbc.constants import *
from py_vbc.derivatives import calc_derivs

path = 'py_vbc/tests/g_test/'

def plot_comparison(x0, yc0, yb0, z):
    """Plots the ratio of baryon to dark matter perturbations as in Naoz
    and Barkana (2005)

    :param x0: py_vbc x values
    :param yc0: py_vbc TF values for dark matter
    :param yb0: py_vbc TF values for baryons
    :returns: the plot (which has been saved)
    :rtype: Figure

    """

    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot the py_vbc output
    ax.plot(x0, yb0/yc0, 'seagreen', label='z = {0:5.3f}'.format(z))
    ax.legend()

    ax.set_ylabel(r'$\delta_b/\delta_c$')

    ax.set_xlabel(r'k (Mpc$^{-1}$)')
    plt.xscale('log')
    plt.yscale('log')
    fig.tight_layout()

    fig.savefig(path+'/ca_pyvbc_ratio.pdf')

g_c = np.loadtxt(path+'g.dat')
zend = 400.0
g_pv = calc_derivs(g_c[:, 0], vbc=0.0, zstart=1000.0, zend=zend,
                   dz=3.0, verbose=True)

plot_comparison(g_c[:, 0], g_pv[:, 0], g_pv[:, 2], z=zend)
plt.show()

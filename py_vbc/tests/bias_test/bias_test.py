import numpy as np
import matplotlib.pyplot as plt

from py_vbc.spectra import calc_power_spec
from py_vbc.constants import *
from py_vbc.derivatives import calc_derivs

path = 'py_vbc/tests/bias_test'

def plot_comparison(x0, yc0, yb0, z, title=None):
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
    ax.plot(x0, yc0, 'seagreen', label='CDM')
    ax.plot(x0, yb0, 'cyan', label='baryons')
    ax.legend()

    ax.set_ylabel(r'P(v$_\mathsf{bc}=30$)/P(v$_\mathsf{bc}=0$)')
    ax.set_xlabel(r'k (Mpc$^{-1}$)')

    if title is not None: ax.set_title('{0} z = {1:3.2f}'.format(title, z))
    
    plt.xscale('log')
    # plt.yscale('log')
    fig.tight_layout()

    fig.savefig(path+'/ca_pyvbc_bias_'+title+'.pdf')

# Load up CICsASS power spec, we only use the k-values in the following
p_c = np.loadtxt(path+'/p.dat')
p_c_t = p_c[:, 1] + p_c[:, 2]

zend = 200.0

# Calculate v_bc = 0.0
g_pv = calc_derivs(p_c[:, 0], vbc=0.0, zstart=1000.0, zend=zend, dz=3.0)
p_pv_c0, p_pv_b0, p_pv_vc0, p_pv_vb0 = calc_power_spec(p_c[:, 0], g_pv, zstart=1000.0)

# Calculate v_bc = 30.0
g_pv = calc_derivs(p_c[:, 0], vbc=30.0, zstart=1000.0, zend=zend, dz=3.0)
p_pv_c30, p_pv_b30, p_pv_vc30, p_pv_vb30 = calc_power_spec(p_c[:, 0], g_pv, zstart=1000.0)

# Calculate the bias = v_bc0/v_bc30
bias_c = p_pv_c30 / p_pv_c0
bias_vc = p_pv_vc30 / p_pv_vc0
bias_b = p_pv_b30 / p_pv_b0
bias_vb = p_pv_vb30 / p_pv_vb0

#d_pv_t = d_pv_c + d_pv_b

# c_ps = np.loadtxt(path+'test_matterpower_z50.dat')
# kh_c_ps = c_ps[:, 0]
# d_c_ps = calc_delta(c_ps[:, 0], c_ps[:, 1])

# camb = [kh_c_ps, d_c_ps]

plot_comparison(p_c[:, 0], bias_c, bias_b, z=zend, title='density')
plot_comparison(p_c[:, 0], bias_vc, bias_vb, z=zend, title='velocity')
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from py_vbc.derivatives import calc_derivs

def calc_bias(g0, gv):
    return gv[:, [0, 2, 5, 6]] / g0[:, [0, 2, 5, 6]]

kmin = 1
kmax = 4000
n = 128
zstart = 1000
zend = 200
dz = 3
k = 10. ** np.linspace(np.log10(kmin), np.log10(kmax), num=n)
vr = 42.03
verbose = True

# integrators = ['RK45', 'DOP853', 'BDF', 'LSODA']
integrators = ['LSODA']

for integrator in integrators:
    g0 = calc_derivs(k, 0, zstart, zend, dz, verbose, integrator)
    gv = calc_derivs(k, vr, zstart, zend, dz, verbose, integrator)
    b = calc_bias(g0, gv)
    out = np.zeros((n, 5))
    out[:, 0] = k
    out[:, 1:] = b
    np.savetxt(f'bias_42_03_kms_{integrator}.txt', out, header='k  delta_c  delta_v  v_c  v_b')

    k *= 0.673
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.semilogx(k, b[:, 0], ls=(0, (8, 2)), c='k',
            label=r'$\delta_{\sf c}$')
    ax.semilogx(k, b[:, 2], ls=(0, (8, 2, 1, 2)), c='k',
                label=r'v$_{\sf c}$')
    ax.semilogx(k, b[:, 1], ls=(0, (4, 2)), c='k',
                label=r'$\delta_{\sf b}$')
    ax.semilogx(k, b[:, 3], ls=(0, (4, 2, 1, 2)), c='k',
                label=r'v$_{\sf b}$')

    ax.text(0.05, 0.95, 'v$_{{\\sf bc,~rec}}$ = {0:.2f} km s$^{{-1}}$'.format(vr),
            transform=ax.transAxes, horizontalalignment='left',
            verticalalignment='top', fontsize='small')

    # ax.text(0.95, 0.95, f'{integrator}',
    #         transform=ax.transAxes, horizontalalignment='right',
    #         verticalalignment='top', fontsize='small')

    ax.set_xlim([5, 1200])
    ax.set_ylim([0.0, 1.2])
    ax.set_ylabel(r'b(k, v$_{\sf bc}$)')
    ax.set_xlabel(r'k (h Mpc$^{-1}$)')
    ax.legend()
    pdf_file = f'bias_vbc_rms_42_03.pdf'
    fig.savefig(pdf_file, bbox_inches='tight', pad_inches=0.02)

    
    

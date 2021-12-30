from py_vbc import run_pyvbc
import numpy as np
import matplotlib.pyplot as plt

n = 250
kmin = 0.1
kmax = 2000
zstart = 1000
zend = 200
dz = 3
vi = 8.44021320
vr = vi * (zstart + 1) / (zend + 1)
v0 = 0.0

# k0, p0 = run_pyvbc(v0, zstart, zend, dz, kmin, kmax, n, verbose=True)
# kv, pv = run_pyvbc(vr, zstart, zend, dz, kmin, kmax, n, verbose=True)

# np.savetxt(f'p_vbc_rms_{vr:.2f}.txt', np.array([kv, pv[0], pv[1],
#                                                 pv[2], pv[3]]).T)
# np.savetxt(f'p_vbc_rms_{v0:.2f}.txt', np.array([k0, p0[0], p0[1],
#                                                 p0[2], p0[3]]).T)

kv, pv0, pv1, pv2, pv3 = np.loadtxt(f'p_vbc_rms_{vr:.2f}.txt', unpack=True)
k0, p00, p01, p02, p03 = np.loadtxt(f'p_vbc_rms_{v0:.2f}.txt', unpack=True)

p0 = np.array([p00, p01, p02, p03])
pv = np.array([pv0, pv1, pv2, pv3])

k0 *= 0.673
fig, ax = plt.subplots(figsize=(4, 4))
ax.semilogx(k0, np.sqrt(pv[0] / p0[0]), ls=(0, (8, 2)), c='#44AA99',
            label=r'$\delta_{\sf c}$')
ax.semilogx(k0, np.sqrt(pv[1] / p0[1]), ls=(0, (4, 2)), c='#CC6677',
            label=r'$\delta_{\sf b}$')
# ax.text(0.05, 0.95, 'v$_{{\\sf bc,~rec}}$ = {0:.2f} km s$^{{-1}}$'.format(vr),
#         transform=ax.transAxes, horizontalalignment='left',
#         verticalalignment='top')
ax.set_xlim([0.5, 1200])
ax.set_ylim([0.0, 1.1])
ax.set_ylabel('b')
ax.set_xlabel(r'k (h Mpc$^{-1}$)')
ax.legend()
pdf_file = f'bias_vbc_rms_{vr:.2f}'
pdf_file = pdf_file.replace('.', '_')
fig.savefig(pdf_file + '.pdf', bbox_inches='tight', pad_inches=0.02)

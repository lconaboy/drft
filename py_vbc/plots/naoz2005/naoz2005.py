import numpy as np
import matplotlib.pyplot as plt
from py_vbc.constants import *
from py_vbc.derivatives import calc_derivs

k = 10.**np.linspace(-4, 3.5, num=400)
delta = {}
for z in (900, 400, 200):
    delta[str(z)] = calc_derivs(k, vbc=0.0, zstart=1000.0, zend=z, dz=3.0, verbose=True)
    
yl0 = [0.0015, 1.5]
yl1 = [0.01, 1.5]
xl = [k.min(), k.max()]

fig = plt.figure(figsize=(6, 8))
ax0 = plt.subplot2grid((6, 8), (0, 0), colspan=8, rowspan=3)
ax1 = plt.subplot2grid((6, 8), (3, 0), colspan=8, rowspan=3, sharex=ax0)

ax0.loglog(k, delta['900'][:, 2] / delta['900'][:, 0], c='c', ls='dashed')
ax1.loglog(k, delta['900'][:, 4] / delta['900'][:, 2], c='c', ls='dashed')

ax0.loglog(k, delta['400'][:, 2] / delta['400'][:, 0], c='m', ls='dotted')
ax1.loglog(k, delta['400'][:, 4] / delta['400'][:, 2], c='m', ls='dotted')

ax0.loglog(k, delta['200'][:, 2] / delta['200'][:, 0],
           c='g', ls='solid')
ax1.loglog(k, delta['200'][:, 4] / delta['200'][:, 2],
           c='g', ls='solid')


ax0.set_ylabel(r'$\delta_{\sf b}/\delta_{\sf c}$')
ax1.set_ylabel(r'$\delta_{\sf T}/\delta_{\sf b}$')
ax1.set_xlabel(r'k (Mpc$^{-1}$)')
# ax0.set_ylim(yl0)
# ax1.set_ylim(yl1)
ax0.set_xlim(xl)
ax1.set_xlim(xl)
ax0.tick_params(axis='x', labelbottom=False)
fig.subplots_adjust(hspace=0)
fig.savefig('naoz2005.pdf', bbox_inches='tight')

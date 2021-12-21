import os
import numpy as np
import matplotlib.pyplot as plt

from py_vbc.constants import *
from py_vbc.interpolations import interpolate_tf, interpolate_tf2d

zs = np.array([1000, 900, 800, 700, 600, 500, 400, 300, 200])

g_spline = interpolate_tf2d('g', zs=zs)
t_spline = interpolate_tf2d('t', zs=zs)

k = 10.**np.linspace(-2, 0, num=100)

# Plot all g
fig, ax = plt.subplots(figsize=(5, 5))
for z in zs:
    ax.loglog(k, g_spline(k, z), label=f'z = {z:d}')
ax.legend()
ax.set_ylabel(r'T$_{\gamma}$')
ax.set_xlabel(r'k (Mpc$^{-1})$')
fig.savefig('all_g.pdf', bbox_inches='tight')

# Plot all g ratios
fig, ax = plt.subplots(figsize=(5, 5))
for z in zs:
    gz_spline = interpolate_tf('g', z)  # k -interpolation from actual z
    y = g_spline(k, z) / gz_spline(k)
    print(f'z={z:d} ratio min/max', y.min(), y.max())
    ax.loglog(k, y, label=f'z = {z:d}')

ax.legend()
# ax.set_ylim([0.9999, 1.0001])
ax.set_ylabel(r'T$_{\sf \gamma,\ int}$')
ax.set_xlabel(r'k (Mpc$^{-1})$')
fig.savefig('ratio_g.pdf', bbox_inches='tight')

# Now compare interpolated zs
fig, ax = plt.subplots(figsize=(5, 5))
for i in range(0, 7, 2):
    z0 = zs[i]
    z1 = zs[i+2]
    zi = zs[i+1]

    print('z0, z1, zi', z0, z1, zi)

    gint = g_spline(k, zi)
    gz = interpolate_tf('g', int(zi))

    y = gint/gz(k)
    print(f'zi={zi:d} ratio min/max', y.min(), y.max())
    ax.loglog(k, y, label=f'z = {zi:d}')

ax.legend()
# ax.set_ylim([0.9999, 1.0001])
ax.set_ylabel(r'T$_{\sf \gamma,\ int}$/T$_{\sf \gamma}$')
ax.set_xlabel(r'k (Mpc$^{-1})$')
fig.savefig('int_g.pdf', bbox_inches='tight')

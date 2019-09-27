# import numpy as np

# def pn():
#     import pynbody
#     import pynbody.filt as filt

#     # Load the snapshot
#     base = '/p/scratch/chpo22/hpo22i/bd/runs/nbody_cubic_512_16384/'
#     i_snap = 51
#     path = base + 'output_{0:05d}'.format(i_snap)
#     f = pynbody.load(path)

#     # Create the cubic filter
#     l = 0.025
#     r = l/2
#     cen = (0.5, 0.5, 0.5)
#     cube = filt.Cuboid(cen[0]-r, cen[1]-r, cen[2]-r,
#                        cen[0]+r, cen[1]+r, cen[2]+r)

#     # Create the subsnap
#     sub = f[cube]
#     # and pick out the dark matter particles
#     sub = sub.d

#     # Now filter out the refined mass particles, to get a map of where
#     # the higher mass particles are
#     ref_mass = sub['mass'].min()
#     contam = sub[filt.HighPass('mass', ref_mass)]

#     # Produce an image of the contamination
#     pynbody.plot.image(contam, cmap='Greys', filename='contam_map.png')

import numpy as np
import matplotlib.pyplot as plt

import pymses
from pymses.analysis import * # Camera, ScalarOperator
from pymses.analysis.splatting import *
from pymses.utils.regions import Box
from pymses.utils import constants as C
from pymses.filters import RegionFilter, PointFunctionFilter

import seren3
from seren3.analysis.visualization import EngineMode, engines

# Load the snapshot
base = '/p/scratch/chpo22/hpo22i/bd/runs/cubic_512_16384/'
i_snap = 41
ro = pymses.RamsesOutput(base, i_snap, metals=True, order='=')
half = True
# s = seren3.load_snapshot(base, i_snap)

# Create the cubic filter
l = 0.0589
filename = 'contam'

# For looking at half of the box, centred on the middle
if half:
    l /= 2
    filename += '_half'

r = l/2
cen = (0.5, 0.5, 0.5)
x_min = [cen[i] - r for i in range(3)]
x_max = [cen[i] + r for i in range(3)]
bounds = (x_min, x_max)
box_filt = Box(bounds)

# Filter the data to the zoom region
dset = ro.particle_source(['level', 'mass', 'epoch'], select_stars=False) # try using level, might have to revert to mass
box_dset = RegionFilter(box_filt, dset)
# sub = s[box_filt]
# sub = sub.d

# dark_filt = lambda dset: (dset['epoch'] == 0.0)
# dark_dset = PointFunctionFilter(dark_filt, box_dset)
dark_dset = box_dset

# Now generate a new filter, this time to include all particles with a
# refinement level less than l_max
l_max = 14
coarse_filt = lambda dset: (dset['level'] < l_max)
coarse_dset = PointFunctionFilter(coarse_filt, dark_dset)
# coarse_dset_iter = coarse_dset.iter_dsets()
# coarse_points = np.zeros(shape=(100000000, 3), dtype=float)
# cc = 0

# # Iterate over the coarse points
# for p in coarse_dset_iter:
#     print(p.npoints)
#     coarse_points[cc:cc+p.npoints, :] = np.array(p.points)

#     cc += p.npoints

# # Remove empty part of array
# coarse_points = coarse_points[0:cc, :]

# coarse_sub = sub[ref_filt]

# # If the above fails (or for comparison), try
# mass_min = np.min(filt_dset['mass'])
# mass_filt = lambda dset: (dset['mass'] > mass_min)


def splatter_projection(dset, cen, l=1.0, los_ax='z', ver_ax='y', opr='mass'):
    r = l/2.0
    op = ScalarOperator(lambda dset: dset[opr], C.Msun)
    mp = SplatterProcessor(dset, ro.info, op)
    cam = Camera(center=cen, line_of_sight_axis=los_ax,
                 up_vector=ver_ax, region_size=(l, l), distance=r,
                 far_cut_depth=r, map_max_size=512)
    proj = mp.process(cam)

    return proj

def projection_plot(proj, proj1=None, log=True, cmap='Greys', ex=None, fig=None, ax=None, filename='plot', title=None, ax_label=None):
    # Extract the map from the projection
    m = proj.map.T

    if log:
        filename = filename + '_log'
        m = np.log10(m)
        cb_label = r'log$_10$(M/M$_\odot$)'
    else:
        cb_label = r'M/M$_\odot$'

    # Create a Figure object if one isn't supplied
    if fig is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    # else use the supplied one
    im = ax.imshow(m, cmap=cmap, extent=ex, vmin=m.min(), vmax=m.max())
    # cb = fig.colorbar(im, ax=ax)
    # cb.set_label(cb_label)
    ax.set_xlabel(ax_label[0])
    ax.set_ylabel(ax_label[1])

    # Do the overplot if one is specified (could generalise this if
    # given a list of projections)
    if proj1 is not None:
        m1 = proj1.map.T
        if log:
            m1 = np.log10(m1)
        ax.imshow(m1, cmap='Reds', extent=ex, alpha=0.5, vmin=m.min(), vmax=m.max())

    # Set title, if specified
    if title is not None:
        ax.set_title(title)

    # Save figure
    fig.tight_layout()
    fig.savefig(filename+'.png')


def get_axes(ax_conf='a'):
    axes = {'a':['z', 'y'], 'b':['x', 'z'], 'c':['x', 'y']}
    return axes[ax_conf]

# eng = engines.SplatterEngine(coarse_sub, 'mass')
# pro = mp.process(cam)

# Position of coarse particles
# coarse_flat = coarse_dset.flatten()
# coarse_points = coarse_flat.points

ex = [cen[0]-r, cen[0]+r, cen[1]-r, cen[1]+r]
ax_confs = ['a', 'b', 'c']
ax_labels = [['x', 'y'], ['y', 'z'], ['x', 'z']]

for i, log in enumerate([True, False]):
    for ax_conf, ax_label in zip(ax_confs, ax_labels):
        fig, ax = plt.subplots(figsize=(6, 6))
        los_ax, ver_ax = get_axes(ax_conf)
        proj = splatter_projection(box_dset, cen=cen, l=l, los_ax=los_ax, ver_ax=ver_ax)
        proj1 = splatter_projection(coarse_dset, cen=cen, l=l, los_ax=los_ax, ver_ax=ver_ax)
        projection_plot(proj, proj1=proj1, log=log, ex=ex, filename=filename+'_'+los_ax+ver_ax, title='z = {0:3.2f}'.format(1.0/ro.info['aexp'] - 1.0), ax_label=ax_label)


# Some AMR stuff
from pymses.filters import CellsToPoints
amr = ro.amr_source(["rho", "vel", "P", "g"])
amr_source = CellsToPoints(amr)

# Try reading individual particle datasets
import numpy as np
import matplotlib.pyplot as plt

import pymses
from pymses.analysis import * # Camera, ScalarOperator
from pymses.analysis.splatting import *
from pymses.utils.regions import Box
from pymses.utils import constants as C
from pymses.filters import RegionFilter, PointFunctionFilter

import seren3
from seren3.analysis.visualization import EngineMode, engines

# Load the snapshot
base = '/p/scratch/chpo22/hpo22i/bd/runs/cubic_512_16384/'
i_snap = 41
ro = pymses.RamsesOutput(base, i_snap, metals=True, order='=')
part = ro.particle_source(["level", "mass"])

for i in range(1, 385):
    loaded = False
    try:
        dset = part.get_domain_dset(i)
        loaded = True
    except:
        print("Cannot load out_{0:05d}".format(i))

    if loaded is True:
        print("Success!")
